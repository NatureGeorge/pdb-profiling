# @Created Date: 2020-01-12 01:27:18 pm
# @Filename: api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-11 04:22:22 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import os
import numpy as np
import pandas as pd
from pyexcel import get_sheet, Sheet
import tablib
from tablib import InvalidDimensions, UnsupportedFormat
from typing import Union, Optional, Iterator, Iterable, Set, Dict, List, Any, Generator, Callable, Tuple
from json import JSONDecodeError
import orjson as json
from pathlib import Path
import aiofiles
from collections import OrderedDict, defaultdict
from unsync import unsync, Unfuture
from pdb_profiling.utils import decompression, related_dataframe, flatten_dict, pipe_out
from pdb_profiling.log import Abclog
from pdb_profiling.fetcher.webfetch import UnsyncFetch

API_LYST: List = sorted(['summary', 'molecules', 'experiment', 'ligand_monomers',
                   'modified_AA_or_NA', 'mutated_AA_or_NA', 'status',
                   'polymer_coverage', 'secondary_structure',
                   'residue_listing', 'binding_sites', 'files', 'observed_residues_ratio',
                   'assembly', 'electron_density_statistics',
                   'cofactor', 'drugbank', 'related_experiment_data'])

BASE_URL: str = 'https://www.ebi.ac.uk/pdbe/'

FTP_URL: str = 'ftp://ftp.ebi.ac.uk/'

FTP_DEFAULT_PATH: str = 'pub/databases/msd/sifts/flatfiles/tsv/uniprot_pdb.tsv.gz'

PDB_ARCHIVE_URL_EBI: str = 'http://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/'
PDB_ARCHIVE_URL_WWPDB: str = 'https://ftp.wwpdb.org/pub/pdb/data/structures/'
# https://ftp.wwpdb.org/pub/pdb/data/structures/obsolete/mmCIF/a0/2a01.cif.gz
# http://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/obsolete/mmCIF/a0/2a01.cif.gz

FUNCS = list()

def dispatch_on_set(keys: Set):
    '''
    Decorator to add new dispatch functions
    '''
    def register(func):
        FUNCS.append((func, set(keys)))
        return func
    return register


def traversePDBeData(query: Any, *args):
    for func, keySet in FUNCS:
        if query in keySet:
            return func(*args)
    else:
        raise ValueError(f'Invalid query: {query}')


def convertJson2other(
        data: Union[List, str, None], 
        append_data: Union[Iterable, Iterator],
        converter: Optional[tablib.Dataset] = None, 
        export_format: str = 'tsv', 
        ignore_headers: Union[bool, int] = False, 
        log_func=print) -> Any:
    '''
    Convert valid json-string/dict into specified format via `tablib.Dataset.export`
    '''
    if converter is None:
        converter = tablib.Dataset()
    try:
        if isinstance(data, str):
            converter.json = data
        elif isinstance(data, List):
            converter.dict = data
        elif data is None:
            pass
        else:
            log_func(f'Invalid type for data`: {type(data)}')
            return None
        for data_to_append in append_data:
            converter.append_col(*data_to_append)
        if ignore_headers:
            converter.headers = None
        return converter.export(export_format)
    except KeyError:
        log_func('Not a valid json-string/dict to convert format via `tablib.Dataset.export`')
    except JSONDecodeError:
        log_func(f'Invalid json string')
    except InvalidDimensions:
        log_func('Invalid data or append_data')
    except UnsupportedFormat:
        log_func(f'Invalid export_format: {export_format}')


class SeqRangeReader(object):
    def __init__(self, name_group):
        self.name = name_group  # ('pdb_id', 'chain_id', 'UniProt')
        self.pdb_range = []
        self.unp_range = []

    def output(self):
        if self.pdb_range:
            pdb_range = json.dumps(self.pdb_range)
            unp_range = json.dumps(self.unp_range)
            return pdb_range, unp_range
        else:
            return self.default_pdb_range, self.default_unp_range

    def check(self, name_group_to_check, data_group):
        self.default_pdb_range = '[[%s, %s]]' % data_group[:2]
        self.default_unp_range = '[[%s, %s]]' % data_group[2:4]

        if self.name == name_group_to_check:
            self.pdb_range.append([int(data_group[0]), int(data_group[1])])
            self.unp_range.append([int(data_group[2]), int(data_group[3])])
        else:
            self.name = name_group_to_check
            self.pdb_range = [[int(data_group[0]), int(data_group[1])]]
            self.unp_range = [[int(data_group[2]), int(data_group[3])]]

        return self.output()


class ProcessPDBe(Abclog):

    converters = {
        'pdb_id': str,
        'chain_id': str,
        'struct_asym_id': str,
        'entity_id': str,
        'author_residue_number': int,
        'residue_number': int,
        'author_insertion_code': str,
        'id': int,
        'interface_id': int,
        'interface_number': int,
        'pdb_code': str,
        'assemble_code': int,
        'assembly_id': int,
        'oper_expression': str
        }
    
    use_existing: bool = False

    @staticmethod
    def yieldTasks(pdbs: Union[Iterable, Iterator], suffix: str, method: str, folder: str, chunksize: int = 25, task_id: int = 0) -> Generator:
        file_prefix = suffix.replace('/', '%')
        method = method.lower()
        if method == 'post':
            url = f'{BASE_URL}{suffix}'
            for i in range(0, len(pdbs), chunksize):
                params = {'url': url, 'data': ','.join(pdbs[i:i+chunksize])}
                yield method, params, os.path.join(folder, f'{file_prefix}+{task_id}+{i}.json')
        elif method == 'get':
            for pdb in pdbs:
                pdb = pdb.lower()
                identifier = pdb.replace('/', '%')
                yield method, {'url': f'{BASE_URL}{suffix}{pdb}'}, os.path.join(folder, f'{file_prefix}+{identifier}.json')
        else:
            raise ValueError(f'Invalid method: {method}, method should either be "get" or "post"')

    @classmethod
    def single_retrieve(cls, pdb: str, suffix: str, method: str, folder: Union[Path, str], semaphore, rate: float = 1.5, **kwargs):
        return UnsyncFetch.single_task(
            task=next(cls.yieldTasks((pdb, ), suffix, method, folder)),
            semaphore=semaphore,
            to_do_func=kwargs.get('to_do_func', cls.process),
            rate=rate)

    @classmethod
    def retrieve(cls, pdbs: Union[Iterable, Iterator], suffix: str, method: str, folder: str, chunksize: int = 20, concur_req: int = 20, rate: float = 1.5, task_id: int = 0, ret_res:bool=True, **kwargs):
        # t0 = time.perf_counter()
        res = UnsyncFetch.multi_tasks(
            cls.yieldTasks(pdbs, suffix, method, folder, chunksize, task_id), 
            cls.process, 
            concur_req=concur_req, 
            rate=rate, 
            logger=cls.logger,
            ret_res=ret_res,
            semaphore=kwargs.get('semaphore', None))
        # elapsed = time.perf_counter() - t0
        # cls.logger.info('{} ids downloaded in {:.2f}s'.format(len(res), elapsed))
        return res
    
    @classmethod
    @unsync
    async def process(cls, path: Union[str, Path, Unfuture]):
        cls.logger.debug('Start to decode')
        if not isinstance(path, (str, Path)):
            path = await path # .result()
        if path is None:
            return path
        path = Path(path)
        suffix = path.name.replace('%', '/').split('+')[0]
        new_path = str(path).replace('.json', '.tsv')
        if Path(new_path).exists() and cls.use_existing:
            return new_path
        async with aiofiles.open(path) as inFile:
            try:
                data = json.loads(await inFile.read())
            except Exception as e:
                cls.logger.error(f"Error in {path}")
                raise e
        res = PDBeDecoder.pyexcel_io(suffix=suffix, data=data)
        if res is not None:
            await pipe_out(
                df=res, path=new_path, 
                format='tsv', mode='w')
            cls.logger.debug(f'Decoded file in {new_path}')
            return new_path
        else:
            cls.logger.warning(f"Without Expected Data ({suffix}): {data}")
            return None
        

class ProcessSIFTS(ProcessPDBe):
    @classmethod
    def related_UNP_PDB(cls, filePath: Union[str, Path], related_unp: Optional[Iterable] = None, related_pdb: Optional[Iterable] = None):
        '''
        Reference
        
            * http://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
            * A summary of the UniProt to PDB mappings showing the UniProt accession
              followed by a semicolon-separated list of PDB four letter codes.
        '''
        filePath = Path(filePath)
        if filePath.is_dir():
            url = FTP_URL+FTP_DEFAULT_PATH
            task = ('ftp', {'url': url}, str(filePath))
            res = UnsyncFetch.multi_tasks([task]).result()
            filePath = decompression(res[0], remove=False, logger=cls.logger)
        elif filePath.is_file() and filePath.exists():
            filePath = str(filePath)
        else:
            raise ValueError('Invalid value for filePath')

        dfrm = pd.read_csv(filePath, sep='\t', header=1)
        pdb_list = list()
        if related_unp is not None:
            dfrm = dfrm[dfrm['SP_PRIMARY'].isin(related_unp)]
        for i in dfrm.index:
            pdb_list.extend(dfrm.loc[i, 'PDB'].split(';'))
        if related_pdb is not None:
            return set(pdb_list) & set(related_pdb), set(dfrm['SP_PRIMARY'])
        else:
            return set(pdb_list), set(dfrm['SP_PRIMARY'])

    @classmethod
    def reformat(cls, path: str) -> pd.DataFrame:
        dfrm = pd.read_csv(path, sep='\t', converters=cls.converters)
        group_info_col = ['pdb_id', 'chain_id', 'UniProt']
        range_info_col = ['pdb_start', 'pdb_end', 'unp_start', 'unp_end']
        reader = SeqRangeReader(group_info_col)
        dfrm[['sifts_pdb_range', 'sifts_unp_range']] = pd.DataFrame(dfrm.apply(
            lambda x: reader.check(tuple(x[i] for i in group_info_col), tuple(
                x[i] for i in range_info_col)),
            axis=1).values.tolist(), index=dfrm.index)
        dfrm = dfrm.drop(columns=range_info_col).drop_duplicates(
            subset=group_info_col, keep='last').reset_index(drop=True)
        dfrm["Entry"] = dfrm["UniProt"].apply(lambda x: x.split('-')[0])
        return dfrm

    @staticmethod
    def dealWithInDe(dfrm: pd.DataFrame) -> pd.DataFrame:
        def get_gap_list(li: List):
            return [li[i+1][0] - li[i][1] - 1 for i in range(len(li)-1)]

        def get_range_diff(lyst_a: List, lyst_b: List):
            array_a = np.array([ran[1] - ran[0] + 1 for ran in lyst_a])
            array_b = np.array([ran[1] - ran[0] + 1 for ran in lyst_b])
            return (array_a - array_b).tolist()

        def add_tage_to_range(df: pd.DataFrame, tage_name: str):
            # ADD TAGE FOR SIFTS
            df[tage_name] = 'Safe'
            # No Insertion But Deletion[Pure Deletion]
            df.loc[df[(df['group_info'] == 1) & (
                df['sifts_unp_pdb_var'] > 0)].index, tage_name] = 'Deletion'
            # Insertion & No Deletion
            df.loc[df[
                (df['group_info'] != 1) &
                (df['var_0_count'] == df['group_info']) &
                (df['unp_GAP_0_count'] == (df['group_info'] - 1))].index, tage_name] = 'Insertion'
            # Insertion & Deletion
            df.loc[df[
                (df['group_info'] != 1) &
                ((df['var_0_count'] != df['group_info']) |
                 (df['unp_GAP_0_count'] != (df['group_info'] - 1)))].index, tage_name] = 'Insertion & Deletion'

        dfrm['pdb_GAP_list'] = dfrm.apply(lambda x: json.dumps(
            get_gap_list(json.loads(x['sifts_pdb_range']))), axis=1)
        dfrm['unp_GAP_list'] = dfrm.apply(lambda x: json.dumps(
            get_gap_list(json.loads(x['sifts_unp_range']))), axis=1)
        dfrm['var_list'] = dfrm.apply(lambda x: json.dumps(get_range_diff(
            json.loads(x['sifts_unp_range']), json.loads(x['sifts_pdb_range']))), axis=1)
        dfrm['delete'] = dfrm.apply(
            lambda x: '-' in x['var_list'], axis=1)
        dfrm['delete'] = dfrm.apply(
            lambda x: True if '-' in x['unp_GAP_list'] else x['delete'], axis=1)
        dfrm['var_0_count'] = dfrm.apply(
            lambda x: json.loads(x['var_list']).count(0), axis=1)
        dfrm['unp_GAP_0_count'] = dfrm.apply(
            lambda x: json.loads(x['unp_GAP_list']).count(0), axis=1)
        dfrm['group_info'] = dfrm.apply(lambda x: len(
            json.loads(x['sifts_pdb_range'])), axis=1)
        dfrm['sifts_unp_pdb_var'] = dfrm.apply(
            lambda x: json.loads(x['var_list'])[0], axis=1)
        add_tage_to_range(dfrm, tage_name='sifts_range_tage')
        return dfrm
    
    '''
    @staticmethod
    def update_range(dfrm: pd.DataFrame, fasta_col: str, unp_fasta_files_folder: str, new_range_cols=('new_sifts_unp_range', 'new_sifts_pdb_range')) -> pd.DataFrame:
        def getSeq(fasta_path: str):
            unpSeq = None
            try:
                unpSeqOb = SeqIO.read(fasta_path, "fasta")
                unpSeq = unpSeqOb.seq
            except ValueError:
                unpSeqOb = SeqIO.parse(fasta_path, "fasta")
                for record in unpSeqOb:
                    if unp_id in record.id.split('|'):
                        unpSeq = record.seq
            return unpSeq

        focus = ('Deletion', 'Insertion & Deletion')
        focus_index = dfrm[dfrm['sifts_range_tage'].isin(focus)].index
        updated_pdb_range, updated_unp_range = list(), list()
        seqAligner = SeqPairwiseAlign()
        for index in focus_index:
            pdbSeq = dfrm.loc[index, fasta_col]
            unp_entry = dfrm.loc[index, "Entry"]
            unp_id = dfrm.loc[index, "UniProt"]
            try:
                fasta_path = os.path.join(
                    unp_fasta_files_folder, f'{unp_id}.fasta')
                unpSeq = getSeq(fasta_path)
            except FileNotFoundError:
                try:
                    fasta_path = os.path.join(
                        unp_fasta_files_folder, f'{unp_entry}.fasta')
                    unpSeq = getSeq(fasta_path)
                except FileNotFoundError:
                    unpSeq = None
            res = seqAligner.makeAlignment(unpSeq, pdbSeq)
            updated_unp_range.append(res[0])
            updated_pdb_range.append(res[1])

        updated_range_df = pd.DataFrame(
            {new_range_cols[0]: updated_unp_range, new_range_cols[1]: updated_pdb_range}, index=focus_index)
        dfrm = pd.merge(dfrm, updated_range_df, left_index=True,
                        right_index=True, how='left')
        dfrm[new_range_cols[0]] = dfrm.apply(lambda x: x['sifts_unp_range'] if pd.isna(
            x[new_range_cols[0]]) else x[new_range_cols[0]], axis=1)
        dfrm[new_range_cols[1]] = dfrm.apply(lambda x: x['sifts_pdb_range'] if pd.isna(
            x[new_range_cols[1]]) else x[new_range_cols[1]], axis=1)
        return dfrm

    @classmethod
    def main(cls, filePath: Union[str, Path], folder: str, related_unp: Optional[Iterable] = None, related_pdb: Optional[Iterable] = None):
        pdbs, _ = cls.related_UNP_PDB(filePath, related_unp, related_pdb)
        res = cls.retrieve(pdbs, 'mappings/all_isoforms/', 'get', folder)
        # return pd.concat((cls.dealWithInDe(cls.reformat(route)) for route in res if route is not None), sort=False, ignore_index=True)
        return res
    '''


class ProcessEntryData(ProcessPDBe):
    @staticmethod
    def related_PDB(pdb_col: str, **kwargs) -> pd.Series:
        dfrm = related_dataframe(**kwargs)
        return dfrm[pdb_col].drop_duplicates()

    @classmethod
    def main(cls, **kwargs):
        pdbs = cls.related_PDB(**kwargs)
        if len(pdbs) > 0:
            res = cls.retrieve(pdbs, **kwargs)
            try:
                return pd.concat((pd.read_csv(route, sep=kwargs.get('sep', '\t'), converters=cls.converters) for route in res if route is not None), sort=False, ignore_index=True)
            except ValueError:
                cls.logger.error('Non-value to concat')
        else:
            return None

    @classmethod
    def unit(cls, pdbs, **kwargs):
        if len(pdbs) > 0:
            res = cls.retrieve(pdbs, **kwargs)
            try:
                return pd.concat((pd.read_csv(route, sep=kwargs.get('sep', '\t'), converters=cls.converters) for route in res if route is not None), sort=False, ignore_index=True)
            except ValueError:
                cls.logger.warning('Non-value to concat')
        else:
            return None

    @staticmethod
    def yieldObserved(dfrm: pd.DataFrame) -> Generator:
        groups = dfrm.groupby(['pdb_id', 'entity_id', 'chain_id'])
        for i, j in groups:
            mod = j.dropna(subset=['chem_comp_id'])
            yield i, len(j[j.observed_ratio.gt(0)]), len(mod[mod.observed_ratio.gt(0)])

    @staticmethod
    def traverse(data: Dict, cols: Tuple, cutoff=50):
        '''
        temp
        '''
        observed_res_count, observed_modified_res_count = cols
        for pdb in data:
            count = 0
            cleaned = 0
            for entity in data[pdb].values():
                for chain in entity.values():
                    if chain[observed_res_count] - chain[observed_modified_res_count] < cutoff:
                        cleaned += 1
                    count += 1
            yield pdb, count, cleaned

    @classmethod
    def pipeline(cls, pdbs: Iterable, folder: str, chunksize: int = 1000):
        for i in range(0, len(pdbs), chunksize):
            related_pdbs = pdbs[i:i+chunksize]
            molecules_dfrm = ProcessEntryData.unit(
                related_pdbs,
                suffix='api/pdb/entry/molecules/',
                method='post',
                folder=folder,
                task_id=i)
            res_listing_dfrm = ProcessEntryData.unit(
                related_pdbs,
                suffix='api/pdb/entry/residue_listing/',
                method='get',
                folder=folder,
                task_id=i)
            modified_AA_dfrm = ProcessEntryData.unit(
                related_pdbs,
                suffix='api/pdb/entry/modified_AA_or_NA/',
                method='post',
                folder=folder,
                task_id=i)
            if modified_AA_dfrm is not None:
                res_listing_dfrm.drop(columns=['author_insertion_code'], inplace=True)
                modified_AA_dfrm.drop(columns=['author_insertion_code'], inplace=True)
                res_listing_mod_dfrm = pd.merge(res_listing_dfrm, modified_AA_dfrm, how='left')
            else:
                res_listing_mod_dfrm = res_listing_dfrm
                res_listing_mod_dfrm['chem_comp_id'] = np.nan
            pro_dfrm = molecules_dfrm[molecules_dfrm.molecule_type.isin(['polypeptide(L)', 'polypeptide(D)'])][['pdb_id', 'entity_id']].reset_index(drop=True)
            pro_res_listing_mod_dfrm = pd.merge(res_listing_mod_dfrm, pro_dfrm)
            data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
            for (pdb_id, entity_id, chain_id), observed_res_count, observed_modified_res_count in cls.yieldObserved(pro_res_listing_mod_dfrm):
                data[pdb_id][entity_id][chain_id]['ob_res'] = observed_res_count
                data[pdb_id][entity_id][chain_id]['ob_moded_res'] = observed_modified_res_count
            with Path(folder, f'clean_pdb_statistic+{i}.tsv').open(mode='w+') as outFile:
                for tup in cls.traverse(data, ('ob_res', 'ob_moded_res')):
                    outFile.write('%s\t%s\t%s\n' % tup)
            with Path(folder, f'clean_pdb_statistic+{i}.json').open(mode='w+') as outFile:
                json.dump(data, outFile)


class PDBeDecoder(object):    

    @staticmethod
    def sync_with_pyexcel(*args) -> Sheet:
        records, *remain = args
        if not len(records):
            return None
        keys = {key for record in records for key in record.keys()}
        head = records[0]
        records[0] = {key: head.get(key, None) for key in keys}
        sheet = get_sheet(records=records)
        if len(remain) > 1:
            append_header, append_value = remain
            append_value = tuple(i if i is not None else '' for i in append_value)
            append_data = [append_header] + [append_value for _ in range(len(records))]
            sheet.column += Sheet(append_data)
        sheet.name_columns_by_row(0)
        return sheet

    @classmethod
    def pyexcel_io(cls, suffix: str, data: Dict) -> Sheet:
        cur_sheet = None
        for res in traversePDBeData(suffix, data):
            try:
                cur_cols = set(cur_sheet.colnames)
                new_sheet = cls.sync_with_pyexcel(*res)
                new_cols = set(new_sheet.colnames)
                # assert cur_cols >= new_cols, f"Unexpected new columns: {new_cols-cur_cols}"
                if not cur_cols >= new_cols:
                    append_header = tuple(new_cols-cur_cols)
                    append_data = [append_header] + [
                        tuple('' for _ in range(len(append_header))) for _ in range(len(cur_sheet))]
                    # NOTE: cur_sheet and new_sheet have colnames while the new Sheet does not
                    cur_sheet.column += Sheet(append_data)
                cur_sheet = get_sheet(
                    records=list(cur_sheet.records)+list(new_sheet.records), 
                    name_columns_by_row=0)
            except AttributeError:
                cur_sheet = cls.sync_with_pyexcel(*res)
            # except TypeError:
            #     continue
        return cur_sheet

    @staticmethod
    def sync_with_tablib(*args) -> tablib.Dataset:
        records, *remain = args
        if not len(records):
            return None
        keys = {key for record in records for key in record.keys()}
        ob = tablib.Dataset()
        records = [OrderedDict(sorted((key, record.get(key, None)) for key in keys))
                   for record in records]
        ob.dict = records
        if len(remain) > 1:
            append_header, append_value = remain
            for i in range(len(append_value)):
                ob.append_col([append_value[i]]*len(records), append_header[i])
        return ob

    @classmethod
    def tablib_io(cls, suffix: str, data: Dict) -> tablib.Dataset:
        cur_ob = None
        for res in traversePDBeData(suffix, data):
            try:
                cur_cols = set(cur_ob.headers)
                new = cls.sync_with_tablib(*res)
                new_cols = set(new.headers)
                assert cur_cols == new_cols, f"Unequal columns: new_cols-cur_cols={new_cols-cur_cols}, cur_cols-new_cols={cur_cols-new_cols}"
                cur_ob.dict += new.dict
            except AttributeError:
                cur_ob = cls.sync_with_tablib(*res)
        return cur_ob

    @staticmethod
    @dispatch_on_set({'api/pdb/entry/status/', 'api/pdb/entry/summary/', 'api/pdb/entry/modified_AA_or_NA/',
                      'api/pdb/entry/mutated_AA_or_NA/', 'api/pdb/entry/cofactor/', 'api/pdb/entry/molecules/',
                      'api/pdb/entry/ligand_monomers/', 'api/pdb/entry/experiment/', 'api/pdb/entry/carbohydrate_polymer/',
                      'api/pdb/entry/electron_density_statistics/',
                      'api/pdb/entry/related_experiment_data/', 'api/pdb/entry/drugbank/'})
    def yieldCommon(data: Dict) -> Generator:
        for pdb in data:
            values = data[pdb]
            for value in values:
                for key in value:
                    if isinstance(value[key], (Dict, List)):
                        value[key] = json.dumps(value[key]).decode('utf-8')
            yield values, ('pdb_id',), (pdb,)

    @staticmethod
    @dispatch_on_set({'api/pdb/entry/polymer_coverage/'})
    def yieldPolymerCoverage(data: Dict) -> Generator:
        for pdb in data:
            molecules = data[pdb]['molecules']
            for entity in molecules:
                chains = entity['chains']
                for chain in chains:
                    observed = chain['observed']
                    for fragement in observed:
                        for key in ('start', 'end'):
                            fragement[key] = json.dumps(fragement[key]).decode('utf-8')
                    yield observed, ('chain_id', 'struct_asym_id', 'entity_id', 'pdb_id'), (chain['chain_id'], chain['struct_asym_id'], entity['entity_id'], pdb)

    @staticmethod
    @dispatch_on_set({'api/pdb/entry/observed_residues_ratio/'})
    def yieldObservedResiduesRatio(data: Dict) -> Generator:
        for pdb in data:
            for entity_id, entity in data[pdb].items():
                yield entity, ('entity_id', 'pdb_id'), (entity_id, pdb)

    @staticmethod
    @dispatch_on_set({'api/pdb/entry/residue_listing/'})
    def yieldResidues(data: Dict) -> Generator:
        for pdb in data:
            molecules = data[pdb]['molecules']
            for entity in molecules:
                chains = entity['chains']
                for chain in chains:
                    residues = chain['residues']
                    for res in residues:
                        if 'multiple_conformers' not in res:
                            res['multiple_conformers'] = None
                        else:
                            res['multiple_conformers'] = json.dumps(res['multiple_conformers']).decode('utf-8')
                    yield residues, ('chain_id', 'struct_asym_id', 'entity_id', 'pdb_id'), (chain['chain_id'], chain['struct_asym_id'], entity['entity_id'], pdb)

    @staticmethod
    @dispatch_on_set({'api/pdb/entry/secondary_structure/'})
    def yieldSecondaryStructure(data: Dict) -> Generator:
        for pdb in data:
            molecules = data[pdb]['molecules']
            for entity in molecules:
                chains = entity['chains']
                for chain in chains:
                    secondary_structure = chain['secondary_structure']
                    for name in secondary_structure:
                        fragment = secondary_structure[name]
                        for record in fragment:
                            for key in record:
                                if isinstance(record[key], (Dict, List)):
                                    record[key] = json.dumps(record[key]).decode('utf-8')
                            if 'sheet_id' not in record:
                                record['sheet_id'] = None
                        yield fragment, ('secondary_structure', 'chain_id', 'struct_asym_id', 'entity_id', 'pdb_id'), (name, chain['chain_id'], chain['struct_asym_id'], entity['entity_id'], pdb)

    @staticmethod
    @dispatch_on_set({'api/pdb/entry/binding_sites/'})
    def yieldBindingSites(data: Dict) -> Generator:
        for pdb in data:
            for site in data[pdb]:
                for tage in ('site_residues', 'ligand_residues'):
                    residues = site[tage]
                    for res in residues:
                        if 'symmetry_symbol' not in res:
                            res['symmetry_symbol'] = None
                    yield residues, ('residues_type', 'details', 'evidence_code', 'site_id', 'pdb_id'), (tage, site['details'], site['evidence_code'], site['site_id'], pdb)

    @staticmethod
    @dispatch_on_set({'api/pdb/entry/assembly/'})
    def yieldAssembly(data: Dict) -> Generator:
        for pdb in data:
            for biounit in data[pdb]:
                entities = biounit['entities']
                for entity in entities:
                    for key in entity:
                        if isinstance(entity[key], (Dict, List)):
                            entity[key] = json.dumps(entity[key]).decode('utf-8')
                keys = list(biounit)
                keys.remove('entities')
                yield entities, tuple(keys)+('pdb_id',), tuple(biounit[key] for key in keys)+(pdb, )

    @staticmethod
    @dispatch_on_set({'api/pdb/entry/files/'})
    def yieldAssociatedFiles(data: Dict) -> Generator:
        for pdb in data:
            for key in data[pdb]:
                for innerKey in data[pdb][key]:
                    record = data[pdb][key][innerKey]
                    if record:
                        yield record, ('innerKey', 'key', 'pdb_id'), (innerKey, key, pdb)
                    else:
                        continue

    @staticmethod
    @dispatch_on_set({'api/mappings/all_isoforms/'})
    def yieldSIFTSRange(data: Dict) -> Generator:
        top_root = next(iter(data))  # PDB_ID or UniProt Isoform ID
        sec_root = next(iter(data[top_root]))  # 'UniProt' or 'PDB'
        child = data[top_root][sec_root]
        thi_root = next(iter(child))
        test_value = child[thi_root]
        # from PDB to UniProt
        if isinstance(test_value, Dict) and sec_root == 'UniProt':
            for uniprot in child:
                name = child[uniprot]['name']
                identifier = child[uniprot]['identifier']
                chains = child[uniprot]['mappings']
                for chain in chains:
                    chain['start'] = json.dumps(chain['start']).decode('utf-8')
                    chain['end'] = json.dumps(chain['end']).decode('utf-8')
                    chain['pdb_id'] = top_root
                    chain[sec_root] = uniprot
                    chain['identifier'] = identifier
                    chain['name'] = name
                yield chains, None
        # from UniProt to PDB
        elif isinstance(test_value, List) and sec_root == 'PDB':
            for pdb in child:
                chains = child[pdb]
                for chain in chains:
                    chain['start'] = json.dumps(chain['start']).decode('utf-8')
                    chain['end'] = json.dumps(chain['end']).decode('utf-8')
                yield chains, ('pdb_id', 'UniProt'), (pdb, top_root)
        else:
            raise ValueError(f'Unexpected data structure for inputted data: {data}')

    @staticmethod
    @dispatch_on_set({'api/pisa/interfacelist/'})
    def yieldPISAInterfaceList(data: Dict):
        for pdb in data:
            records = data[pdb]['interfaceentries']
            for record in records:
                flatten_dict(record, 'structure_1')
                flatten_dict(record, 'structure_2')
            flatten_dict(data[pdb], 'page_title', False)
            cols = sorted(i for i in data[pdb].keys() if i != 'interfaceentries')
            yield records, cols, tuple(data[pdb][col] for col in cols)
    
    @staticmethod
    @dispatch_on_set({'api/pisa/interfacedetail/'})
    def yieldPISAInterfaceDetail(data: Dict):
        edge_cols1 = ('structure', 'interface_atoms', 'interface_residue', 'interface_area', 'solvation_energy')
        edge_cols2 = ('structure', 'interface_atoms', 'interface_residues', 'interface_area', 'solvation_energy')
        for pdb in data:
            try:
                records = data[pdb]['interface_detail']
            except KeyError:
                raise ValueError(f"Without Expected interface_detail: {data}")
            del records['bonds']
            for col in edge_cols1: flatten_dict(records['interface_structure_1'], col)
            for col in edge_cols2: flatten_dict(records['interface_structure_2'], col)
            flatten_dict(data[pdb], 'page_title', False)
            flatten_dict(records, 'interface_structure_1')
            flatten_dict(records, 'interface_structure_2')
            flatten_dict(data[pdb], 'interface_detail')
            cols = sorted(i for i in data[pdb].keys() if i != 'interface_detail.residues')
            yield data[pdb]['interface_detail.residues']['residue1']['residue']['residue_array'], cols, tuple(data[pdb][col] for col in cols)
            yield data[pdb]['interface_detail.residues']['residue2']['residue']['residue_array'], cols, tuple(data[pdb][col] for col in cols)

    @staticmethod
    @dispatch_on_set({'swissmodel/repository/uniprot/'})
    def yieldSMR(data: Dict):
        cols = ('sequence_length',
                'ac', 'id', 'isoid')
        uniprot_entries = data['result']['uniprot_entries']
        
        assert len(uniprot_entries) == 1, f"Unexpected length of uniprot_entries: {uniprot_entries}"
        
        for col in ('ac', 'id', 'isoid'):
            data['result'][col] = uniprot_entries[0].get(col, None)
        
        yield data['result']['structures'], cols, tuple(data['result'][col] for col in cols)


class PDBeModelServer(Abclog):
    '''
    Implement ModelServer API
    '''
    
    root = 'model-server/v1/'
    headers =  {'accept': 'text/plain', 'Content-Type': 'application/json'}
    api_sets = {'atoms', 'residueInteraction', 'assembly', 'full', 'ligand'
                'residueSurroundings', 'symmetryMates', 'query-many'}
    
    @classmethod
    def yieldTasks(cls, pdbs, suffix: str, method: str, folder: str, data_collection, params) -> Generator:
        if data_collection is None:
            assert method == 'get', 'Invalid method!'
            for pdb in pdbs:
                args = dict(
                    url=f'{BASE_URL}{cls.root}{pdb}/{suffix}?',
                    headers=cls.headers,
                    params=params)
                yield method, args, os.path.join(folder, f'{pdb}_subset.{params.get("encoding", "cif")}')
        else:
            assert method == 'post', 'Invalid method!'
            for pdb, data in zip(pdbs, data_collection):
                args = dict(
                    url=f'{BASE_URL}{cls.root}{pdb}/{suffix}?',
                    headers=cls.headers,
                    params=params,
                    data=data)
                yield method, args, os.path.join(folder, f'{pdb}_subset.{params.get("encoding", "cif")}')

    @classmethod
    def retrieve(cls, pdbs, suffix: str, method: str, folder: str, data_collection=None, params=None, concur_req: int = 20, rate: float = 1.5, ret_res:bool=True, **kwargs):
        if params is None:
            params = {'model_nums': 1, 'encoding': 'cif'}
        res = UnsyncFetch.multi_tasks(
            cls.yieldTasks(pdbs, suffix, method, folder,
                           data_collection, params),
            concur_req=concur_req,
            rate=rate,
            logger=cls.logger,
            ret_res=ret_res,
            semaphore=kwargs.get('semaphore', None))
        return res
    
    @classmethod
    def single_retrieve(cls, pdb: str, suffix: str, method: str, folder: Union[Path, str], semaphore, data_collection=None, params=None, rate: float = 1.5):
        if params is None:
            params = {'encoding': 'cif'}
        if data_collection is not None:
            data_collection = (data_collection, )
        return UnsyncFetch.single_task(
            task=next(cls.yieldTasks((pdb, ), suffix, method, folder,
                                     data_collection, params)),
            semaphore=semaphore,
            rate=rate)


class PDBArchive(Abclog):
    '''
    Download files from PDB Archive

    * wwPDB/RCSB: PDB_ARCHIVE_URL_WWPDB: str = 'https://ftp.wwpdb.org/pub/pdb/data/structures/'
    * EBI: PDB_ARCHIVE_URL_EBI: str = 'http://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/'
    '''
    root = PDB_ARCHIVE_URL_EBI
    api_sets = {f'{i}/{j}/' for i in ('obsolete', 'divided')
                for j in ('mmCIF', 'pdb', 'XML')}

    @classmethod
    def yieldTasks(cls, pdbs, suffix: str, file_suffix: str, folder: Path) -> Generator:
        for pdb in pdbs:
            args = dict(url=f'{cls.root}{suffix}{pdb[1:3]}/{pdb}{file_suffix}')
            yield 'get', args, folder/f'{pdb}{file_suffix}'

    @classmethod
    def retrieve(cls, pdbs, suffix: str, folder: Path, file_suffix: str = '.cif.gz', concur_req: int = 20, rate: float = 1.5, ret_res:bool=True, **kwargs):
        res = UnsyncFetch.multi_tasks(
            cls.yieldTasks(pdbs, suffix, file_suffix, folder),
            concur_req=concur_req,
            rate=rate,
            logger=cls.logger,
            ret_res=ret_res,
            semaphore=kwargs.get('semaphore', None))
        return res
    
    @classmethod
    def single_retrieve(cls, pdb: str, suffix: str, folder: Path, semaphore, file_suffix: str = '.cif.gz', rate: float = 1.5):
        return UnsyncFetch.single_task(
            task=next(cls.yieldTasks((pdb, ), suffix, file_suffix, folder)),
            semaphore=semaphore,
            rate=rate)

# TODO: Chain UniProt ID Mapping -> ProcessSIFTS -> ProcessPDBe
# TODO: Deal with oligomeric PDB

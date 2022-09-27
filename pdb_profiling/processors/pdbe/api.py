# @Created Date: 2020-01-12 01:27:18 pm
# @Filename: api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2022-09-14 09:01:12 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Union, Optional, Iterator, Iterable, Dict, List, Any, Generator, Tuple
import orjson as json
from pathlib import Path
from aiofiles import open as aiofiles_open
from unsync import unsync, Unfuture
from random import choice
from pdb_profiling.processors.recordbase import IdentifierBase
from pdb_profiling.utils import flatten_dict, pipe_out, dumpsParams
from pdb_profiling.log import Abclog
from pdb_profiling.fetcher.webfetch import UnsyncFetch
from pdb_profiling.processors.transformer import Dict2Tabular
from pdb_profiling.exceptions import WithoutExpectedKeyError, InvalidFileContentError
from pdb_profiling.ensure import EnsureBase
from tenacity import wait_random, stop_after_attempt, retry_if_exception_type, RetryError

ensure = EnsureBase()
msc_rt_kw = dict(wait=wait_random(max=1), stop=stop_after_attempt(3), retry=retry_if_exception_type(InvalidFileContentError))

PDBE_URL: str = 'https://www.ebi.ac.uk/pdbe/'
EBI_FTP_URL: str = 'ftp://ftp.ebi.ac.uk/'
HTTP_EBI_FTP_URL: str = 'https://ftp.ebi.ac.uk/'
HTTP_WWPDB_FTP_URL: str = 'https://ftp.wwpdb.org/'

SIFTS_FTP_DEFAULT_PATH: str = f'{HTTP_EBI_FTP_URL}pub/databases/msd/sifts/flatfiles/tsv/uniprot_pdb.tsv.gz'

# https://ftp.wwpdb.org/pub/pdb/data/structures/obsolete/mmCIF/a0/2a01.cif.gz
# http://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/obsolete/mmCIF/a0/2a01.cif.gz
# http://ftp-versioned.wwpdb.org/pdb_versioned/data/entries/wm/pdb_00002wmg/pdb_00002wmg_xyz_v1-2.cif.gz

FUNCS = []


def mask_ib(i, default='', raise_error=False):
    if i.source == 'PDB' and i.level == 'entry':
        return 'pdb_id'
    elif i.source == 'UniProt':
        return 'UniProt'
    elif raise_error:
        raise AssertionError('Unexpected Case!')
    else:
        return default

def str_number_converter(x):
    try:
        return int(x)
    except ValueError:
        return -100000


def dispatch_on_set(*keys):
    '''
    Decorator to add new dispatch functions
    '''
    def register(func):
        FUNCS.append((func, frozenset(keys)))
        return func
    return register


def traverseSuffixes(query: Any, *args):
    for func, keySet in FUNCS:
        if query in keySet:
            return func(*args)
    else:
        raise ValueError(f'Invalid query: {query}')


class ProcessPDBe(Abclog):

    headers = {'Connection': 'close', 'Content-Type': 'application/json'}

    converters = {
        'pdb_id': str,
        'chain_id': str,
        'struct_asym_id': str,
        'entity_id': str_number_converter,
        'author_residue_number': int,
        'residue_number': str_number_converter,
        'author_insertion_code': str,
        'id': int,
        'interface_id': int,
        'interface_number': int,
        'pdb_code': str,
        'assemble_code': int,
        'assembly_id': int,
        'oper_expression': str,
        'structure_1.range': str,
        'structure_2.range': str,
        'alt_code': str,
        'sheet_id': str_number_converter
    }

    @classmethod
    def yieldTasks(cls, pdbs: Union[Iterable, Iterator], suffix: str, method: str, folder: Union[str, Path], chunksize: int = 25, task_id: int = 0) -> Generator:
        file_prefix = suffix.replace('/', '%')
        method = method.lower()
        if method == 'post':
            url = f'{PDBE_URL}{suffix}'
            for i in range(0, len(pdbs), chunksize):
                params = {'headers': cls.headers, 'url': url, 'data': ','.join(pdbs[i:i+chunksize])}
                yield method, params, folder/f'{file_prefix}+{task_id}+{i}.json'
        elif method == 'get':
            for pdb in pdbs:
                identifier = pdb.replace('/', '%')
                yield method, {'headers': cls.headers, 'url': f'{PDBE_URL}{suffix}{pdb}'}, folder/f'{file_prefix}+{identifier}.json'
        else:
            raise ValueError(
                f'Invalid method: {method}, method should either be "get" or "post"')

    @classmethod
    def single_retrieve(cls, pdb: str, suffix: str, method: str, folder: Union[Path, str], semaphore, rate: float = 1.5, **kwargs):
        return UnsyncFetch.single_task(
            task=next(cls.yieldTasks((pdb, ), suffix, method, folder)),
            semaphore=semaphore,
            to_do_func=kwargs.get('to_do_func', cls.process),
            rate=rate)

    @classmethod
    def retrieve(cls, pdbs: Union[Iterable, Iterator], suffix: str, method: str, folder: Union[str, Path], chunksize: int = 20, concur_req: int = 20, rate: float = 1.5, task_id: int = 0, ret_res: bool = True, **kwargs):
        # t0 = time.perf_counter()
        res = UnsyncFetch.multi_tasks(
            cls.yieldTasks(pdbs, suffix, method, folder, chunksize, task_id),
            cls.process,
            concur_req=concur_req,
            rate=rate,
            ret_res=ret_res,
            semaphore=kwargs.get('semaphore', None))
        # elapsed = time.perf_counter() - t0
        # cls.logger.info('{} ids downloaded in {:.2f}s'.format(len(res), elapsed))
        return res

    @classmethod
    @unsync
    @ensure.make_sure_complete(**msc_rt_kw)
    async def json2tsv(cls, suffix:str, ori_path: Union[str, Path], path: Union[str, Path]):
        cls.logger.debug('Start to decode')
        async with aiofiles_open(ori_path) as handle:
            try:
                data = json.loads(await handle.read())
            except Exception as e:
                cls.logger.error(f"Error in '{ori_path}'")
                raise e
        res = Dict2Tabular.pyexcel_io(traverseSuffixes(suffix, data))
        if res is not None:
            if isinstance(res, Generator):
                count = 0
                for r in res:
                    if r is not None:
                        await pipe_out(df=r, path=path, format='tsv', mode='a' if count else 'w')
                        count += 1
                if not count:
                    cls.logger.debug(f"Without Expected Data ({suffix}): {data}")
                    return None
            else:
                await pipe_out(df=res, path=path, format='tsv', mode='w')
            cls.logger.debug(f"Decoded file in '{path}'")
            return path
        else:
            cls.logger.debug(f"Without Expected Data ({suffix}): {data}")
            return None

    @classmethod
    @unsync
    async def process(cls, path: Union[str, Path, Unfuture]):
        if not isinstance(path, (str, Path)):
            path = await path
        if path is None:
            return
        path = Path(path)
        suffix = path.name.replace('%', '/').split('+')[0]
        new_path = Path(str(path).replace('.json', '.tsv'))
        try:
            return await cls.json2tsv(suffix=suffix, ori_path=path, path=new_path)
        except RetryError:
            cls.logger.error(f"Retry failed for: {path.name} -> {new_path.name}")
            raise


class PDBeDecoder(object):

    @staticmethod
    @dispatch_on_set('api/pdb/entry/status/', 'api/pdb/entry/summary/', 'api/pdb/entry/modified_AA_or_NA/',
                     'api/pdb/entry/mutated_AA_or_NA/', 'api/pdb/entry/cofactor/', 'api/pdb/entry/molecules/',
                     'api/pdb/entry/entities/',
                     'api/pdb/entry/ligand_monomers/', 'api/pdb/entry/experiment/', 'api/pdb/entry/carbohydrate_polymer/',
                     'api/pdb/entry/electron_density_statistics/', 'api/pdb/entry/related_experiment_data/',
                     'api/pdb/entry/drugbank/', 'api/mappings/best_structures/',
                     'graph-api/pdb/mutated_AA_or_NA/', 'graph-api/pdb/modified_AA_or_NA/',
                     'graph-api/mappings/best_structures/', 'graph-api/compound/atoms/',
                     'graph-api/compound/bonds/', 'graph-api/compound/summary/',
                     'graph-api/compound/cofactors/', 'graph-api/pdb/funpdbe/',
                     'graph-api/pdb/bound_excluding_branched/',
                     'graph-api/pdb/bound_molecules/', 'graph-api/pdb/ligand_monomers/',
                     'api/validation/global-percentiles/entry/', 'api/validation/summary_quality_scores/entry/',
                     'api/validation/key_validation_stats/entry/', 'api/validation/xray_refine_data_stats/entry/',
                     'api/validation/vdw_clashes/entry/', 'api/validation/outliers/all/',
                     'api/validation/nmr_cyrange_cores/entry/', # TODO: 2tablar
                     'api/validation/nmr_ensemble_clustering/entry/'
                     )
    def yieldCommon(data: Dict) -> Generator:
        for pdb in data:
            values = data[pdb]
            for value in values:
                for key in value:
                    if isinstance(value[key], (Dict, List)):
                        value[key] = json.dumps(value[key]).decode('utf-8')
            yield values, (mask_ib(IdentifierBase(pdb), '_code_'),), (pdb,)

    @staticmethod
    @dispatch_on_set('api/pdb/entry/polymer_coverage/')
    def yieldPolymerCoverage(data: Dict) -> Generator:
        for pdb in data:
            molecules = data[pdb]['molecules']
            for entity in molecules:
                chains = entity['chains']
                for chain in chains:
                    observed = chain['observed']
                    for fragement in observed:
                        for key in ('start', 'end'):
                            fragement[key] = json.dumps(
                                fragement[key]).decode('utf-8')
                    yield observed, ('chain_id', 'struct_asym_id', 'entity_id', 'pdb_id'), (chain['chain_id'], chain['struct_asym_id'], entity['entity_id'], pdb)

    @staticmethod
    @dispatch_on_set('api/pdb/entry/observed_residues_ratio/')
    def yieldObservedResiduesRatio(data: Dict) -> Generator:
        for pdb in data:
            for entity_id, entity in data[pdb].items():
                yield entity, ('entity_id', 'pdb_id'), (entity_id, pdb)

    @staticmethod
    @dispatch_on_set('api/pdb/entry/residue_listing/')
    def yieldResidues(data: Dict) -> Generator:
        for pdb in data:
            molecules = data[pdb]['molecules']
            for entity in molecules:
                chains = entity['chains']
                for chain in chains:
                    residues = chain['residues']
                    for res in residues:
                        if 'multiple_conformers' not in res:
                            res['multiple_conformers'] = ''
                        else:
                            res['multiple_conformers'] = json.dumps(
                                res['multiple_conformers']).decode('utf-8')
                    yield residues, ('chain_id', 'struct_asym_id', 'entity_id', 'pdb_id'), (chain['chain_id'], chain['struct_asym_id'], entity['entity_id'], pdb)

    @staticmethod
    @dispatch_on_set('api/pdb/entry/secondary_structure/', 'graph-api/pdb/secondary_structure/')
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
                                    record[key] = json.dumps(
                                        record[key]).decode('utf-8')
                            if 'sheet_id' not in record:
                                record['sheet_id'] = None
                        yield fragment, ('secondary_structure', 'chain_id', 'struct_asym_id', 'entity_id', 'pdb_id'), (name, chain['chain_id'], chain['struct_asym_id'], entity['entity_id'], pdb)

    @staticmethod
    @dispatch_on_set('api/pdb/entry/binding_sites/')
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
    @dispatch_on_set('api/pdb/entry/assembly/')
    def yieldAssembly(data: Dict) -> Generator:
        for pdb in data:
            for biounit in data[pdb]:
                entities = biounit['entities']
                for entity in entities:
                    for key in entity:
                        if isinstance(entity[key], (Dict, List)):
                            entity[key] = json.dumps(
                                entity[key]).decode('utf-8')
                keys = list(biounit)
                keys.remove('entities')
                yield entities, tuple(keys)+('pdb_id',), tuple(biounit[key] for key in keys)+(pdb, )

    @staticmethod
    @dispatch_on_set('api/pdb/entry/files/')
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
    @dispatch_on_set('api/mappings/all_isoforms/', 'api/mappings/uniprot/',
                     'api/mappings/uniprot_segments/', 'api/mappings/isoforms/',
                     'api/mappings/uniref90/', 'api/mappings/homologene_uniref90/',
                     'api/mappings/interpro/', 'api/mappings/pfam/',
                     'api/mappings/cath/', 'api/mappings/cath_b/',
                     'api/mappings/scop/', 'api/mappings/go/',
                     'api/mappings/ec/', 'api/mappings/ensembl/',
                     'api/mappings/hmmer/', 'api/mappings/sequence_domains/',
                     'api/mappings/structural_domains/', 'api/mappings/homologene/',
                     'api/mappings/uniprot_to_pfam/', 'api/mappings/uniprot_publications/',
                     'graph-api/mappings/uniprot/', 'graph-api/mappings/uniprot_segments/',
                     'graph-api/mappings/all_isoforms/', 'graph-api/mappings/',
                     'graph-api/mappings/isoforms/', 'graph-api/mappings/ensembl/',
                     'graph-api/mappings/homologene/', 'graph-api/mappings/sequence_domains/',
                     'api/mappings/', 'api/nucleic_mappings/', 'api/nucleic_mappings/rfam/', 
                     'api/nucleic_mappings/sequence_domains/', 'graph-api/mappings/ec/'
                     # 'graph-api/uniprot/'
                     )
    def yieldSIFTSAnnotation(data: Dict) -> Generator:
        valid_annotation_set = {'UniProt', 'Ensembl', 'Pfam', 'CATH',
                                'CATH-B', 'SCOP', 'InterPro', 'GO', 'EC', 'Homologene', 'HMMER', 'Rfam'}
        for top_root in data:
            # top_root: PDB_ID or else ID
            if data[top_root].keys() <= valid_annotation_set:
                # from PDB to ('UniProt', 'Ensembl', 'Pfam', 'CATH', 'CATH-B', 'SCOP', 'InterPro', 'GO', 'EC', 'Homologene', 'HMMER')
                # from PDB_ENTITY (i.e. graph-api/mappings/homologene/)
                # OR: from Uniprot (i.e. api/mappings/uniprot_to_pfam/)
                for sec_root in data[top_root]:
                    child = data[top_root][sec_root]
                    for annotation in child:
                        chains = child[annotation]['mappings']
                        for chain in chains:
                            for key, value in chain.items():
                                chain[key] = json.dumps(value).decode(
                                    'utf-8') if isinstance(value, Dict) else value
                            for key, value in child[annotation].items():
                                if key == 'mappings':
                                    continue
                                chain[key] = json.dumps(value).decode(
                                    'utf-8') if isinstance(value, Dict) else value
                            chain[mask_ib(IdentifierBase(top_root), raise_error=True)] = top_root
                            chain[sec_root] = annotation
                        yield chains, None
            elif len(data[top_root].keys()) == 1 and 'PDB' in data[top_root].keys():
                # from UniProt to PDB
                for sec_root in data[top_root]:
                    child = data[top_root][sec_root]
                    for pdb in child:
                        chains = child[pdb]
                        for chain in chains:
                            chain['start'] = json.dumps(
                                chain['start']).decode('utf-8')
                            chain['end'] = json.dumps(
                                chain['end']).decode('utf-8')
                        yield chains, ('pdb_id', 'UniProt'), (pdb, top_root)
            else:
                raise ValueError(
                    f'Unexpected data structure for inputted data: {data}')

    @staticmethod
    @dispatch_on_set('api/pisa/asiscomponent/')
    def yield_pisa_asiscomponent(data: Dict):
        for pdb in data:
            if data[pdb]['status'] != 'Ok' or 'assembly_detail' not in data[pdb]:
                raise WithoutExpectedKeyError(f"Without Expected interfacelist info: {data}")
            try:
                records = data[pdb]['assembly_detail']['engaged_interfaces_list']['engaged_interfaces_array']
            except KeyError:
                raise WithoutExpectedKeyError(f"Without Expected interfacelist info: {data}")
            yield records, ('pdb_id',), (pdb,)

    @staticmethod
    @dispatch_on_set('api/pisa/interfacelist/')
    def yieldPISAInterfaceList(data: Dict):
        for pdb in data:
            try:
                records = data[pdb]['interfaceentries']
            except KeyError:
                raise WithoutExpectedKeyError(
                    f"Without Expected interface_detail: {data}")
            for record in records:
                flatten_dict(record, 'structure_1')
                flatten_dict(record, 'structure_2')
            yield records, ('pdb_id', 'assembly_id'), (pdb, data[pdb]['page_title']['assemble_code'])

    @staticmethod
    @dispatch_on_set('api/pisa/interfacedetail/')
    def yieldPISAInterfaceDetail(data: Dict):
        usecols = (
            'pdb_code', 'assemble_code', 'interface_number',
            'interface_detail.interface_structure_1.structure.selection',
            'interface_detail.interface_structure_2.structure.selection')
        # 'interface_atoms', 'interface_residue', 'interface_area', 'solvation_energy'
        edge_cols1 = ('structure',)
        # 'interface_atoms', 'interface_residues', 'interface_area', 'solvation_energy'
        edge_cols2 = ('structure',)
        for pdb in data:
            try:
                records = data[pdb]['interface_detail']
            except KeyError:
                raise WithoutExpectedKeyError(
                    f"Without Expected interface_detail: {data}")
            del records['bonds']
            for col in edge_cols1:
                flatten_dict(records['interface_structure_1'], col)
            for col in edge_cols2:
                flatten_dict(records['interface_structure_2'], col)
            flatten_dict(data[pdb], 'page_title', False)
            flatten_dict(records, 'interface_structure_1')
            flatten_dict(records, 'interface_structure_2')
            flatten_dict(data[pdb], 'interface_detail')
            # cols = sorted(i for i in data[pdb].keys() if i != 'interface_detail.residues')
            yield data[pdb]['interface_detail.residues']['residue1']['residue']['residue_array'], usecols, tuple(data[pdb][col] for col in usecols)
            yield data[pdb]['interface_detail.residues']['residue2']['residue']['residue_array'], usecols, tuple(data[pdb][col] for col in usecols)

    @staticmethod
    @dispatch_on_set('graph-api/residue_mapping/')
    def graph_api_residue_mapping(data: Dict):
        '''
        * <https://www.ebi.ac.uk/pdbe/graph-api/residue_mapping/:pdbId/:entityId/:residueNumber>
        * <https://www.ebi.ac.uk/pdbe/graph-api/residue_mapping/:pdbId/:entityId/:residueStart/:residueEnd>

        NOTE: only yield UniProt Residue Related Data
        '''
        cols = (
            'pdb_id', 'entity_id', 'chain_id', 'struct_asym_id',
            'residue_number', 'author_residue_number',
            'author_insertion_code', 'observed', 'UniProt')

        for pdb_id in data:
            assert len(data[pdb_id]) == 1, f"Unexpected Cases: {pdb_id}"
            molecules = data[pdb_id][0]
            for chain in molecules['chains']:
                for residue in chain['residues']:
                    yield list({**dict(zip(cols, (
                        pdb_id, molecules['entity_id'], chain['auth_asym_id'],
                        chain['struct_asym_id'], residue['residue_number'],
                        residue['author_residue_number'], residue['author_insertion_code'],
                        residue['observed'], feature_tag))), **feature} for feature_tag, feature in residue['features']['UniProt'].items()), None

    @staticmethod
    @dispatch_on_set('graph-api/pdb/sequence_conservation/')
    def sequence_conservation(data: Dict):
        for pdb in data:
            yield [{
                'pdb_id': pdb,
                'entity_id': data[pdb]['entity_id'],
                'length': data[pdb]['length'],
                'residue_number': val['start'],
                'conservation_score': val['conservation_score'],
                'letter_array': json.dumps(tuple(i['oneLetterCode'] for i in val['amino'])).decode('utf-8'),
                'proba_array': json.dumps(tuple(i['probability'] for i in val['amino'])).decode('utf-8')}
                for val in data[pdb]['data']], None
            # letter_array, proba_array = zip(*((i['letter'], i['proba']) for i in val['amino']))

    @staticmethod
    @dispatch_on_set('graph-api/pdb/funpdbe_annotation/depth/',
                     'graph-api/pdb/funpdbe_annotation/cath-funsites/',
                     'graph-api/pdb/funpdbe_annotation/3Dcomplex/',
                     'graph-api/pdb/funpdbe_annotation/akid/',
                     'graph-api/pdb/funpdbe_annotation/3dligandsite/',
                     'graph-api/pdb/funpdbe_annotation/camkinet/',
                     'graph-api/pdb/funpdbe_annotation/canSAR/',
                     'graph-api/pdb/funpdbe_annotation/ChannelsDB/',
                     'graph-api/pdb/funpdbe_annotation/dynamine/',
                     'graph-api/pdb/funpdbe_annotation/FoldX/',
                     'graph-api/pdb/funpdbe_annotation/MetalPDB/',
                     'graph-api/pdb/funpdbe_annotation/M-CSA/',
                     'graph-api/pdb/funpdbe_annotation/p2rank/',
                     'graph-api/pdb/funpdbe_annotation/Missense3D/',
                     'graph-api/pdb/funpdbe_annotation/POPScomp_PDBML/',
                     'graph-api/pdb/funpdbe_annotation/ProKinO/',
                     'graph-api/pdb/funpdbe_annotation/14-3-3-pred/',
                     'graph-api/pdb/funpdbe_annotation/'
                     )
    def funpdbe_resources(data: Dict):
        for pdb in data:
            info = data[pdb]
            for val in info:
                for annotation in val['annotations']:
                    yield annotation['site_residues'], ('pdb_id', 'origin', 'evidence_codes', 'site_id', 'label'), (pdb, val['origin'], val['evidence_codes'], annotation['site_id'], annotation['label'])

    @staticmethod
    @dispatch_on_set('graph-api/pdbe_pages/rfam/',
                     'graph-api/pdbe_pages/annotations/',
                     'graph-api/pdbe_pages/uniprot_mapping/',
                     'graph-api/pdbe_pages/binding_sites/',
                     'graph-api/pdbe_pages/interfaces/',
                     'graph-api/pdbe_pages/secondary_structure/',
                     'graph-api/pdbe_pages/domains/',
                     'graph-api/uniprot/unipdb/',
                     'graph-api/uniprot/annotations/',
                     'graph-api/uniprot/interface_residues/',
                     'graph-api/uniprot/ligand_sites/',
                     'graph-api/uniprot/secondary_structures/',
                     'graph-api/uniprot/domains/',
                     'graph-api/uniprot/sequence_conservation/')
    def graph_api_data_common(data: Dict):
        for pdb in data:
            id_type = 'pdb_id' if len(pdb) == 4 else 'UniProt'
            for info in data[pdb]['data']:
                if 'additionalData' in info:
                    flatten_dict(info, 'additionalData')
                com_keys = tuple(key for key in info.keys()
                                 if key != 'residues')
                yield info['residues'], (id_type,)+com_keys, (pdb,)+tuple(info[key] for key in com_keys)

    @staticmethod
    @dispatch_on_set('graph-api/pdb/bound_molecule_interactions/')
    def graph_api_bound(data: Dict):
        for pdb in data:
            info = data[pdb]
            for interactions in info:
                ret = [{j: json.dumps(i[j]).decode('utf-8') for j in i.keys()}
                       for i in interactions['interactions']]
                yield ret, ('pdb_id', 'bm_id'), (pdb, interactions['bm_id'])

    @staticmethod
    @dispatch_on_set('api/validation/protein-ramachandran-sidechain-outliers/entry/', 'api/validation/RNA_pucker_suite_outliers/entry/')
    def yield_protein_ramachandran_sidechain_outlier(data):
        for pdb in data:
            for tage in data[pdb]:
                residues = data[pdb][tage]
                yield residues, ('_type_', 'pdb_id'), (tage, pdb)

    @staticmethod
    @dispatch_on_set('api/validation/rama_sidechain_listing/entry/', 'api/validation/residuewise_outlier_summary/entry/',
                     'api/validation/protein-RNA-DNA-geometry-outlier-residues/entry/')
    def yield_rama_sidechain_listing(data):
        for pdb in data:
            molecules = data[pdb]['molecules']
            for entity in molecules:
                chains = entity['chains']
                for chain in chains:
                    models = chain['models']
                    for model in models:
                        residues = model['residues']
                        yield residues, ('chain_id', 'struct_asym_id', 'model_id', 'entity_id', 'pdb_id'), (chain['chain_id'], chain['struct_asym_id'], model['model_id'], entity['entity_id'], pdb)

    @staticmethod
    @dispatch_on_set('graph-api/uniprot/superposition/')
    def yield_unp_pdb_struct_cluster(data):
        for unp in data:
            for segment_id, segment in enumerate(data[unp]):
                clusters = segment['clusters']
                for sub_cluster_id, sub_cluster in enumerate(clusters):
                    yield sub_cluster, ('pdbekb_cluster', 'segment_start', 'segment_end', 'UniProt'), (f'{segment_id}_{sub_cluster_id}', segment['segment_start'], segment['segment_end'], unp)


class PDBeModelServer(object):
    '''
    Implement ModelServer API
    '''

    pdbe_root = f'{PDBE_URL}model-server/v1/'
    rcsb_root = 'https://models.rcsb.org/v1/'
    root = rcsb_root
    headers = {'Connection': 'close', 'accept': 'text/plain', 'Content-Type': 'application/json'}
    api_set = frozenset(('atoms', 'residueInteraction', 'assembly', 'full', 'ligand'
                         'residueSurroundings', 'symmetryMates', 'query-many'))

    @classmethod
    def task_unit(cls, pdb, suffix, method, folder, data_collection, params, filename='_subset'):
        if data_collection is None:
            assert method == 'get', 'Invalid method!'
            args = dict(
                url=f'{cls.root}{pdb}/{suffix}?{dumpsParams(params)}',
                headers=cls.headers)
        else:
            assert method == 'post', 'Invalid method!'
            args = dict(
                url=f'{cls.root}{pdb}/{suffix}?{dumpsParams(params)}',
                headers=cls.headers,
                data=data_collection)
        return method, args, folder/f'{pdb}{filename}.{params.get("encoding", "cif")}'

    @classmethod
    def single_retrieve(cls, pdb: str, suffix: str, method: str, folder: Union[Path, str], semaphore, params=None, data_collection=None, rate: float = 1.5, filename='_subset'):
        if params is None or len(params) == 0:
            params = {'model_nums': 1, 'encoding': 'cif'}
        return UnsyncFetch.single_task(
            task=cls.task_unit(pdb, suffix, method, folder,
                               data_collection, params, filename=filename),
            semaphore=semaphore,
            rate=rate)


class PDBeCoordinateServer(object):

    roots = (f'{PDBE_URL}coordinates/', 'https://cs.litemol.org/')
    headers = {'Connection': 'close', 'accept': 'text/plain'}
    api_set = frozenset(('ambientResidues', 'assembly', 'backbone', 'cartoon', 'chains'
                         'entities', 'full', 'het', 'ligandInteraction', 'residueRange',
                         'residues', 'sidechain', 'symmetryMates', 'trace', 'water'))

    def __init__(self, root: str = 'litemol'):
        if root == 'random':
            self.root = choice(self.roots)
        elif root == 'ebi':
            self.root = self.roots[0]
        elif root == 'litemol':
            self.root = self.roots[1]
        else:
            raise ValueError("root should be (ebi, litemol, random)")

    def __repr__(self):
        return f'<CoordinateServerAPI: {self.root}>'

    def task_unit(self, pdb_id, suffix: str, params, folder):
        args = dict(
            url=f'{self.root}{pdb_id}/{suffix}?',
            headers=self.headers,
            params=params)
        return 'get', args, Path(folder)/f'{pdb_id}_{dumpsParams(params)}.{params.get("encoding", "cif")}'

    def single_retrieve(self, pdb_id: str, suffix: str, params: Dict, folder: Union[Path, str], semaphore, rate: float = 1.5):
        return UnsyncFetch.single_task(
            task=self.task_unit(pdb_id, suffix, params, folder),
            semaphore=semaphore,
            rate=rate)


class FetchBase:

    root_dispatch: Dict[str, str] = {
        'wwPDB': HTTP_WWPDB_FTP_URL,
        'RCSB': HTTP_WWPDB_FTP_URL,
        'EBI': HTTP_EBI_FTP_URL,
        'PDBe': PDBE_URL,
        'UNP_FTP': 'https://ftp.uniprot.org/',
        'UNP': 'https://www.uniprot.org/',
    }

    def __init__(self, root: str, path: str, root_dispatch: bool = True):
        if root_dispatch:
            self.url: str = self.root_dispatch[root] + path
        else:
            self.url: str = self.root + path

    @unsync
    async def download(self, semaphore, path: Union[str, Path], rate: float = 1.5):
        return await UnsyncFetch.fetch_file(semaphore, method='get', info=dict(url=self.url), path=path, rate=rate)


class PDBArchive(object):
    '''
    Download files from PDB Archive

    * wwPDB/RCSB: PDB_ARCHIVE_URL_WWPDB: str = 'https://ftp.wwpdb.org/pub/pdb/data/structures/'
    * EBI: PDB_ARCHIVE_URL_EBI: str = 'http://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/'
    '''
    pdbe_root = f'{HTTP_EBI_FTP_URL}pub/databases/pdb/data/structures/'
    rcsb_root = f'{HTTP_WWPDB_FTP_URL}pub/pdb/data/structures/'
    root = rcsb_root
    api_set = frozenset(f'{i}/{j}/' for i in ('obsolete', 'divided')
                        for j in ('mmCIF', 'pdb', 'XML'))
    file_dict = {
        'mmCIF': '.cif.gz',
        'pdb': '.ent.gz',
        'XML': '.xml.gz'
    }

    @staticmethod
    def wrap_id(pdb_id, suffix):
        if suffix.endswith('pdb/'):
            return f"pdb{pdb_id}"
        else:
            return pdb_id

    @classmethod
    def get_file_suffix(cls, api_suffix):
        for key, value in cls.file_dict.items():
            if key in api_suffix:
                return value
        raise AssertionError(
            f"Unexpected Case for api_suffix: {api_suffix}, {cls.file_dict}")

    @classmethod
    def task_unit(cls, pdb: str, suffix: str, file_suffix: str, folder: Path):
        args = dict(
            url=f'{cls.root}{suffix}{pdb[1:3]}/{cls.wrap_id(pdb, suffix)}{cls.get_file_suffix(suffix)}')
        return 'get', args, folder/f'{pdb}{file_suffix}'

    @classmethod
    def yieldTasks(cls, pdbs, suffix: str, file_suffix: str, folder: Path) -> Generator:
        for pdb in pdbs:
            yield cls.task_unit(pdb, suffix, file_suffix, folder)

    """@classmethod
    def retrieve(cls, pdbs, suffix: str, folder: Path, file_suffix: Optional[str] = None, concur_req: int = 20, rate: float = 1.5, ret_res: bool = True, **kwargs):
        res = UnsyncFetch.multi_tasks(
            cls.yieldTasks(pdbs, suffix, file_suffix, folder),
            concur_req=concur_req,
            rate=rate,
            ret_res=ret_res,
            semaphore=kwargs.get('semaphore', None))
        return res"""

    @classmethod
    def single_retrieve(cls, pdb, suffix: str, folder: Path, semaphore, file_suffix: Optional[str] = None, rate: float = 1.5):
        if file_suffix is None:
            file_suffix = cls.get_file_suffix(suffix)
        return UnsyncFetch.single_task(
            task=cls.task_unit(pdb, suffix, file_suffix, folder),
            semaphore=semaphore,
            rate=rate)


class PDBVersioned(PDBArchive):
    '''
    Download files from PDB Versioned

    * wwPDB Versioned: PDB_ARCHIVE_VERSIONED_URL: str = 'http://ftp-versioned.wwpdb.org/pdb_versioned/data/entries/'

    >>> PDBVersioned.single_retrieve(
        ('2wmg', '_v1-2'), 'entries/', 
        init_folder_from_suffix(Base.get_folder(), 'pdb-versioned/entries'), 
        Base.get_web_semaphore()).result()
    '''
    root = 'http://ftp-versioned.wwpdb.org/pdb_versioned/data/'
    api_set = frozenset(('entries/', 'removed/'))

    @classmethod
    def task_unit(cls, pdb_with_version: Tuple, suffix: str, file_suffix: str, folder: Path):
        pdb, version_info = pdb_with_version
        file_name = f'pdb_0000{pdb}_xyz{version_info}{file_suffix}'
        args = dict(url=f'{cls.root}{suffix}{pdb[1:3]}/pdb_0000{pdb}/{file_name}')
        return 'get', args, folder/file_name


class PDBeKBAnnotations(object):
    ftp_root = f"{EBI_FTP_URL}pub/databases/pdbe-kb/annotations/"
    https_root = f"{HTTP_EBI_FTP_URL}pub/databases/pdbe-kb/annotations/"
    root = https_root
    api_set = frozenset({
        '14-3-3-pred/', '3DComplex/',
        '3DLigandSite/', 'AKID/',
        'COSPI-Depth/', 'CamKinet/',
        'ChannelsDB/', 'Covalentizer/',
        'DynaMine/', 'FireProtDB/',
        'FoldX/', 'KnotProt/',
        'M-CSA/', 'MetalPDB/',
        'Missense3D/', 'P2rank/',
        'POPScomp_PDBML/', 'ProKinO/',
        'Scop3P/', 'canSAR/',
        'cath-funsites/', 'webNMA/'})
    
    @staticmethod
    def wrap_id(pdb_id, suffix):
        if suffix == 'M-CSA/':
            return f"{pdb_id}-mcsa"
        else:
            return pdb_id
    
    @classmethod
    def task_unit(cls, pdb: str, suffix: str, folder: Path):
        pdb_ = cls.wrap_id(pdb, suffix)
        args = dict(
            url=f'{cls.root}{suffix}{pdb[1:3]}/{pdb_}.json')
        return 'ftp' if cls.root == cls.ftp_root else 'get', args, folder/f'{pdb_}.json'
    
    @classmethod
    def single_retrieve(cls, pdb, suffix: str, folder: Path, semaphore, rate: float = 1.5):
        return UnsyncFetch.single_task(
            task=cls.task_unit(pdb, suffix, folder),
            semaphore=semaphore,
            rate=rate)
    
    @staticmethod
    def yieldPDBeKBAnnotations(data):
        for chain in data['chains']:
            yield chain['residues'], ('data_resource', 'pdb_id', 'chain_id'), (data['data_resource'], data['pdb_id'], chain['chain_label'])

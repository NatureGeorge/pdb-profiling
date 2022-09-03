# @Created Date: 2020-09-17 12:02:37 am
# @Filename: api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2022-09-03 05:08:33 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Union, Optional, Iterable, Dict, List, Generator, Tuple
from pathlib import Path
from numpy import nan
from pandas import DataFrame, Series
from pdb_profiling.log import Abclog
from pdb_profiling.fetcher.webfetch import UnsyncFetch
from pdb_profiling.cython.cyrange import range_len
from pdb_profiling.utils import dumpsParams
from pdb_profiling.processors.uniprot.record import UniProts
from pdb_profiling.warnings import PossibleObsoletedUniProtIsoformWarning
from urllib.parse import quote
from slugify import slugify
from unsync import unsync
from warnings import warn


BASE_URL = 'https://www.ebi.ac.uk/proteins/api/'


class ProteinsAPI(Abclog):
    '''
    Implement The Proteins REST API

    * <https://www.ebi.ac.uk/proteins/api/doc/index.html>
    '''
    
    headers = {'Accept': 'application/json'}
    api_set = frozenset((
        'proteins', 'proteins/covid-19/entries', 'proteins/interaction/', 'proteins/',
        'coordinates', 'coordinates/', 'coordinates/location/',
        'uniparc', 'uniparc/accession/', 'uniparc/best/guess',
        'uniparc/dbreference/', 'uniparc/proteome/', 'uniparc/sequence',  # NOTE: uniparc/sequence use POST method!
        'uniparc/upi/', 'variation/dbsnp/', 'variation/hgvs/', 'variation/',
        'proteomics', 'proteomics/', 'proteomics-ptm', 'proteomics-ptm/',
        'mutagenesis', 'mutagenesis/',
        ))
    
    @classmethod
    def get_file_suffix(cls) -> str:
        res = cls.headers["Accept"].split('/')[1]
        assert res in ('json', 'xml', 'x-gff'), f"Unexpected Case: {cls.headers}"
        return res

    @classmethod
    def task_unit(cls, suffix: str, params: Dict, folder: Path, identifier:Optional[str]=None) -> Tuple:
        args = dict(
            url=f'{BASE_URL}{suffix}' if identifier is None else f'{BASE_URL}{suffix}{quote(identifier)}',
            headers=cls.headers,
            params=params)
        return 'get', args, folder/(suffix if suffix[-1] == '/' else suffix+'_')/f'{slugify(identifier)+"_"+dumpsParams(params) if identifier is not None else dumpsParams(params)}.{cls.get_file_suffix()}'

    @classmethod
    def yieldTasks(cls, suffix: str, params_collection: Iterable[Dict], folder: Path, identifiers: Optional[Iterable[str]]) -> Generator:
        # https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=-1&ensembl=ENST00000554444
        if identifiers is None:
            for params in params_collection:
                yield cls.task_unit(suffix, params, folder)
        else:
            for identifier, params in zip(identifiers, params_collection):
                yield cls.task_unit(suffix, params, folder, identifier)

    """@classmethod
    def retrieve(cls, suffix: str, params_collection: Iterable[Dict], folder: Union[Path, str], identifiers: Optional[Iterable[str]] = None, concur_req: int = 20, rate: float = 1.5, ret_res: bool = True, **kwargs):
        assert suffix in cls.api_set, f"Invalid suffix! Valid set is \n{cls.api_set}"
        folder = Path(folder)
        res = UnsyncFetch.multi_tasks(
            cls.yieldTasks(suffix, params_collection, folder, identifiers),
            concur_req=concur_req,
            rate=rate,
            ret_res=ret_res,
            semaphore=kwargs.get('semaphore', None))
        return res"""
    
    @classmethod
    def single_retrieve(cls, suffix: str, params:Dict, folder: Union[Path, str], semaphore, identifier:Optional[str]=None, rate: float = 1.5):
        assert suffix in cls.api_set, f"Invalid suffix! Valid set is \n{cls.api_set}"
        folder = Path(folder)
        return UnsyncFetch.single_task(
            task=cls.task_unit(suffix, params, folder, identifier),
            semaphore=semaphore,
            rate=rate)

    @classmethod
    def query_sequence(cls, params: Dict, data: Dict, folder: Union[Path, str], fileName:str, semaphore, rate: float = 1.5):
        '''
        Implement `uniparc/sequence`
        '''
        folder = Path(folder)
        args = dict(
            url=f'{BASE_URL}uniparc/sequence',
            headers=cls.headers,
            params=params,
            data=data)
        return UnsyncFetch.single_task(
            task=('post', args, folder/f'{fileName}.{cls.get_file_suffix()}'),
            semaphore=semaphore,
            rate=rate)
        
    @classmethod
    def flat_data(cls, data, surface=False, except_keys=[]):
        if isinstance(data, dict):
            if len(data) == 1 and 'value' in data:
                return cls.flat_data(data['value'], except_keys=except_keys)
            else:
                new = dict()
                for key, value in data.items():
                    if key in except_keys:
                        new[key] = value
                    else:
                        new[key] = cls.flat_data(value, except_keys=except_keys)
                return new
        elif isinstance(data, list):
            if len(data) == 1 and not surface:
                return cls.flat_data(data[0], except_keys=except_keys)
            else:
                return [cls.flat_data(i, except_keys=except_keys) for i in data]
        else:
            return data

    @classmethod
    @unsync
    async def pipe_summary(cls, data: List):
        #"""
        if len(data) > 1:
            flag = []
            for idx in range(len(data)):
                if '-' not in data[idx]["accession"]:
                    flag.append(idx)
            if flag:
                for idx in range(len(data)):
                    if idx not in flag:
                        del data[idx]
        #"""
        if len(data) > 1:
            cls.logger.warning(f"Unexpected Length from ProteinsAPI.pipe_summary: {len(data)}, {[data[i]['accession'] for i in range(len(data))]}")
            '''
            Special Case Like: 
                                P04745 (AMY1A_HUMAN)
                                P0DTE7 (AMY1B_HUMAN)
                                P0DTE8 (AMY1C_HUMAN)
                                with identical sequence
                                Causion from UniProt-KB: 'Three distinct genes (AMY1A, AMY1B and AMY1C), located in a gene cluster on 1p21, encode proteins sharing the same peptidic sequence.'
            Another Special Case Like: NP_001316851
                                P0DP23 (CALM1_HUMAN)
                                P0DP24 (CALM2_HUMAN)
                                P0DP25 (CALM3_HUMAN)
                                with identical sequence
            TODO: deal with these cases
            '''
            # raise AssertionError("With Unexpected length!")
        elif (len(data) == 0) or (data[0] is None):
            return
        data = data[0]
        dbReferences_lyst = []
        common_cols = ('type', 'isoform')
        if 'dbReferences' not in data.keys():
            data['dbReferences'] = []
        for i in data['dbReferences']:
            if i['type'] == 'RefSeq':
                dbReferences_lyst.append({
                    **{key: i.get(key, nan) for key in common_cols},
                    **dict(accession=data['accession'], protein=i['id'], transcript=i['properties']['nucleotide sequence ID'], gene=nan),
                })
            elif i['type'] == 'Ensembl':
                dbReferences_lyst.append({
                    **{key: i.get(key, nan) for key in common_cols},
                    **dict(accession=data['accession'], protein=i['properties']['protein sequence ID'], transcript=i['id'], gene=i['properties']['gene ID'])})
        dbReferences_df = DataFrame(dbReferences_lyst)
        other_dbReferences_df = DataFrame((i for i in data['dbReferences'] if i['type'] not in ('RefSeq', 'Ensembl')))
        other_dbReferences_df['accession'] = data['accession']
        if 'isoform' not in other_dbReferences_df.columns:
            other_dbReferences_df['isoform'] = ''
        else:
            other_dbReferences_df['isoform'].fillna('', inplace=True)
        if 'features' in data.keys():
            features_df = DataFrame(data['features']); features_df['accession'] = data['accession']
            if 'molecule' not in features_df.columns:
                features_df['molecule'] = ''
        else:
            features_df = None
        iso_df, int_df = None, None
        if 'comments' in data.keys():
            for comment in data['comments']:
                if comment['type'] == 'ALTERNATIVE_PRODUCTS':
                    iso_df = DataFrame(
                        cls.flat_data(comment['isoforms'], surface=True, except_keys={'sequence'}))
                    iso_df['accession'] = data['accession']
                    iso_df.rename(columns={'ids': 'isoform'}, inplace=True)
                    iso_df['else_iso'] = iso_df.isoform.apply(lambda x: x[1:] if isinstance(x, List) else nan)
                    focus_index = iso_df[(~iso_df.else_iso.isnull())].index
                    iso_df.loc[focus_index, 'isoform'] = iso_df.loc[focus_index, 'isoform'].apply(lambda x: x[0])
                    iso_df.rename(columns={'sequence': 'VAR_SEQ'}, inplace=True)
                    alt_seq_df = await UniProts.fetch_VAR_SEQ_from_DB((data['accession'], ), "%s+VAR_SEQ" % data['accession'], with_name_suffix=False)
                    
                    ac_len = data['sequence']['length']
                    ac_seq = data['sequence']['sequence']
                    
                    if len(alt_seq_df) > 0:
                        assert all(alt_seq_df.ftId.ne('')), f"{alt_seq_df[alt_seq_df.ftId.eq('')]}"
                        alt_seq_dict = alt_seq_df[["ftId", "before_len", "after_len", "after", "begin", "end"]].to_dict('list')
                        iso_alt_info = iso_df.VAR_SEQ.apply(lambda x: UniProts.get_alt_interval_base(x, alt_seq_dict) if not isinstance(x, float) else x)
                        outdated = iso_df[iso_alt_info.apply(lambda x: len(x) if not isinstance(x, float) else x).eq(0)]
                        if len(outdated) > 0:
                            outdated_isos = outdated.isoform.tolist()
                            warn(str(outdated_isos), PossibleObsoletedUniProtIsoformWarning)
                            iso_df = iso_df.drop(index=outdated.index).reset_index(drop=True)
                            iso_alt_info = iso_alt_info.drop(index=outdated.index).reset_index(drop=True)
                        iso_df[["iso_range", "length", "sequence"]] = iso_alt_info.apply(lambda x: UniProts.get_affected_interval(x, ac_len, ac_seq)).apply(Series)
                    else:
                        iso_df['iso_range'] = nan
                        iso_df['length'] = nan
                        iso_df['sequence'] = nan
                    iso_df['iso_range_len'] = iso_df.iso_range.apply(range_len)
                    iso_df.loc[iso_df.sequenceStatus.eq('displayed').idxmax(), ["length", "sequence"]] = (ac_len, ac_seq)

                if comment['type'] == 'INTERACTION':
                    int_df = DataFrame(comment['interactions'])
                    if 'accession1' not in int_df.columns:
                        int_df['accession1'] = data['accession']
                    else:
                        int_df['accession1'].fillna(data['accession'], inplace=True)
                    assert int_df[['accession1', 'accession2', 'interactor1', 'interactor2']].isnull(
                    ).sum().sum() == 0, int_df[int_df.accession2.isnull() | int_df.interactor1.isnull() | int_df.interactor2.isnull()]
                    int_df_cols = frozenset(int_df.columns)
                    assert int_df_cols <= {'accession1', 'accession2', 'gene', 'interactor1',
                                           'interactor2', 'organismDiffer', 'experiments', 'chain1', 'chain2'}, int_df_cols
                    int_df['accession'] = data['accession']

        info_data = {key: data.get(key, nan) for key in (
            'accession', 'id', 'proteinExistence', 'info', 'organism', 'secondaryAccession', 'protein', 'gene')}
        info_data['sequence'] = data['sequence']['sequence']
        info_data['length'] = data['sequence']['length']

        return dbReferences_df, other_dbReferences_df, iso_df, features_df, int_df, info_data

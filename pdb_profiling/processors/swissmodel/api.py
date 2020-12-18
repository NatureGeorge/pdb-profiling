# @Created Date: 2020-08-27 10:35:34 am
# @Filename: api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-08-27 10:35:41 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pdb_profiling.fetcher.webfetch import UnsyncFetch
from pdb_profiling.utils import pipe_out, init_folder_from_suffix, init_semaphore, a_read_csv, dumpsParams
from pdb_profiling.processors.transformer import Dict2Tabular
from typing import Union, Dict, Generator, Set, Any, Optional, List
from pathlib import Path
from pdb_profiling.log import Abclog
from unsync import unsync
from aiofiles import open as aiofiles_open
import orjson as json
from pandas import Series, concat, isna

BASE_URL: str = 'https://swissmodel.expasy.org/'


class SMR(Abclog):
    '''
    Implement SWISS-MODEL Repository API

        * <https://swissmodel.expasy.org/docs/smr_openapi>
        * <https://swissmodel.expasy.org/docs/repository_help#smr_api>
    
    >>> SMR.retrieve(
        ('Q6NZ36', 'P12755'), 
        init_folder_from_suffix(yourfolder, 'swissmodel/repository/uniprot/'))
    '''

    root = 'repository/uniprot/'
    headers = {'accept': 'text/plain'}
    use_existing: bool = True

    @classmethod
    def set_folder(cls, folder: Union[Path, str]):
        cls.folder = init_folder_from_suffix(Path(folder), 'swissmodel/repository/uniprot/')

    @classmethod
    @unsync
    async def set_web_semaphore(cls, web_semaphore_values):
        cls.web_semaphore = await init_semaphore(web_semaphore_values)

    @staticmethod
    def yieldSMR(data: Dict):
        cols = ('sequence_length',
                'ac', 'id', 'isoid')
        uniprot_entries = data['result']['uniprot_entries']

        assert len(
            uniprot_entries) == 1, f"Unexpected length of uniprot_entries: {uniprot_entries}"

        for col in ('ac', 'id', 'isoid'):
            data['result'][col] = uniprot_entries[0].get(col, None)
        
        for val in data['result']['structures']:
            keys = tuple(val.keys())
            for key in keys:
                if isinstance(val[key], (List, Dict)):
                    val[key] = json.dumps(val[key]).decode('utf-8')

        yield data['result']['structures'], cols, tuple(data['result'][col] for col in cols)

    @classmethod
    def task_unit(cls, unp, params, file_format, folder):
        args = dict(
            url=f'{BASE_URL}{cls.root}{unp}.{file_format}',
            headers=cls.headers,
            params=params)
        return 'get', args, Path(folder)/f'{unp}_{dumpsParams(params)}.{file_format}'

    @classmethod
    def yieldTasks(cls, unps, params: Dict, file_format: str, folder: Union[Path, str]) -> Generator:
        for unp in unps:
            yield cls.task_unit(unp, params, file_format, folder)

    @classmethod
    def retrieve(cls, unps, folder: Optional[Union[Path, str]]=None, params: Dict = dict(provider='swissmodel'), concur_req: int = 20, rate: float = 1.5, file_format: str = 'json', ret_res: bool = True, **kwargs):
        assert file_format in ('json', 'pdb'), "Invalid file format"
        res = UnsyncFetch.multi_tasks(
            cls.yieldTasks(unps, params, file_format, cls.folder if folder is None else folder),
            cls.process,
            concur_req=concur_req,
            rate=rate,
            ret_res=ret_res,
            semaphore=kwargs.get('semaphore', cls.web_semaphore))
        return res

    @classmethod
    def single_retrieve(cls, unp: str, folder: Optional[Union[Path, str]]=None, semaphore=None, params: Dict = dict(provider='swissmodel'), rate: float = 1.5, file_format: str = 'json'):
        assert file_format in ('json', 'pdb'), "Invalid file format"

        return UnsyncFetch.single_task(
            task=cls.task_unit(unp, params, file_format, cls.folder if folder is None else folder),
            semaphore=cls.web_semaphore if semaphore is None else semaphore,
            to_do_func=cls.process,
            rate=rate)

    @classmethod
    @unsync
    async def process(cls, path):
        cls.logger.debug('Start to decode SMR JSON')
        if not isinstance(path, (str, Path)):
            path = await path
        if path is None or str(path).endswith('.pdb'):
            return path
        new_path = Path(str(path).replace('.json', '.tsv'))
        if new_path.exists() and cls.use_existing and (new_path.stat().st_size > 0):
            return new_path
        async with aiofiles_open(path) as inFile:
            try:
                data = json.loads(await inFile.read())
            except Exception as e:
                cls.logger.error(f"Error in '{path}'")
                raise e
        res = Dict2Tabular.pyexcel_io(cls.yieldSMR(data))
        if res is not None:
            if isinstance(res, Generator):
                count = 0
                for r in res:
                    if r is not None:
                        await pipe_out(df=r, path=new_path, format='tsv', mode='a' if count else 'w')
                        count += 1
                if not count:
                    cls.logger.debug(f"Without Expected Data (swissmodel/repository/uniprot/): {data}")
                    return None
            else:
                await pipe_out(df=res, path=new_path, format='tsv', mode='w')
            cls.logger.debug(f"Decoded file in '{new_path}'")
            return new_path
        else:
            cls.logger.warning(f"Without Expected Data (swissmodel/repository/uniprot/): {data}")
            return None

    @classmethod
    @unsync
    async def to_dataframe(cls, path):
        path = await path
        if path is None:
            return None
        df = await a_read_csv(path, sep="\t", keep_default_na=False, na_values=['NULL', 'null', ''])
        df.rename(columns={'ac': 'UniProt', 'sequence_length': 'unp_len', 'id': 'identifier'}, inplace=True)
        df['unp_range'] = df.apply(lambda x: f"[[{x['from']},{x['to']}]]", axis=1)
        if 'identity' in df.columns:
            df.identity = df.identity/100
        if 'qmean' in df.columns:
            return concat((
                df.drop(columns=['from', 'to', 'qmean']),
                df.qmean.apply(lambda x: json.loads(x) if not isna(x) else default_qmean).apply(Series)), axis=1)
            # add_prefix('qmean.')
        else:
            return df.drop(columns=['from', 'to'])


default_qmean = {
    "acc_agreement_norm_score": None,
    "acc_agreement_z_score": None,
    "avg_local_score": None,
    "avg_local_score_error": None,
    "cbeta_norm_score": None,
    "cbeta_z_score": None,
    "interaction_norm_score": None,
    "interaction_z_score": None,
    "packing_norm_score": None,
    "packing_z_score": None,
    "qmean4_norm_score": None,
    "qmean4_z_score": None,
    "qmean6_norm_score": None,
    "qmean6_z_score": None,
    "ss_agreement_norm_score": None,
    "ss_agreement_z_score": None,
    "torsion_norm_score": None,
    "torsion_z_score": None
}

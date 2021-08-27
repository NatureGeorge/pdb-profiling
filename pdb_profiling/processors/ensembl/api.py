# @Created Date: 2020-09-26 01:54:56 pm
# @Filename: api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-09-26 01:55:04 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Dict, Tuple, Iterable, Generator, Union, Optional
from pathlib import Path
from pdb_profiling.log import Abclog
from pdb_profiling.fetcher.webfetch import UnsyncFetch


BASE_URL = 'https://rest.ensembl.org/'


class EnsemblAPI(Abclog):
    '''
    Implement The Ensembl REST API

    * <https://rest.ensembl.org/documentation/>
    * DEMO: <https://rest.ensembl.org/sequence/id/ENST00000288602?content-type=text/x-fasta;type=protein>
    '''
    headers = {"Content-Type": "text/x-fasta"}
    api_set = frozenset(('sequence/id/', 'archive/id/'))

    @classmethod
    def get_file_suffix(cls, headers: Optional[Dict]) -> str:
        res = headers["Content-Type"].split('/')[1]
        assert res in ('plain', 'x-seqxml+xml', 'x-fasta', 'json'), f"Unexpected Case: {cls.headers}"
        return res.replace('x-', '').replace('seqxml+', '')

    @classmethod
    def task_unit(cls, suffix: str, identifier: str, params: Optional[Dict], folder: Path, headers:Optional[Dict]) -> Tuple:
        headers = cls.headers if headers is None else headers
        if params is not None:
            args = dict(
                url=f'{BASE_URL}{suffix}{identifier}',
                headers=headers,
                params=params)
        else:
            args = dict(
                url=f'{BASE_URL}{suffix}{identifier}',
                headers=headers)
        return 'get', args, folder/f'{identifier}.{cls.get_file_suffix(headers)}'

    @classmethod
    def yieldTasks(cls, suffix: str, identifiers: Iterable[str], params_collection: Iterable[Optional[Dict]], folder: Path, headers: Optional[Dict]) -> Generator:
        for identifier, params in zip(identifiers, params_collection):
            yield cls.task_unit(suffix, identifier, params, folder, headers)
    
    """@classmethod
    def retrieve(cls, suffix: str, identifiers: Iterable[str], params_collection: Iterable[Optional[Dict]], folder: Union[Path, str], concur_req: int = 20, rate: float = 1.5, ret_res: bool = True, headers: Optional[Dict] = None, **kwargs):
        assert suffix in cls.api_set, f"Invalid suffix! Valid set is \n{cls.api_set}"
        folder = Path(folder)
        res = UnsyncFetch.multi_tasks(
            cls.yieldTasks(suffix, identifiers, params_collection, folder, headers),
            concur_req=concur_req,
            rate=rate,
            ret_res=ret_res,
            semaphore=kwargs.get('semaphore', None))
        return res"""

    @classmethod
    def single_retrieve(cls, suffix: str, identifier: str, params: Optional[Dict], folder: Union[Path, str], semaphore, rate: float = 1.5, headers: Optional[Dict] = None):
        assert suffix in cls.api_set, f"Invalid suffix! Valid set is \n{cls.api_set}"
        folder = Path(folder)
        return UnsyncFetch.single_task(
            task=cls.task_unit(suffix, identifier, params, folder, headers),
            semaphore=semaphore,
            rate=rate)

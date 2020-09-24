# @Created Date: 2020-09-17 12:02:37 am
# @Filename: api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-09-17 12:02:45 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Union, Optional, Iterator, Iterable, Set, Dict, List, Any, Generator, Callable, Tuple
from pathlib import Path
from pdb_profiling.log import Abclog
from pdb_profiling.fetcher.webfetch import UnsyncFetch


BASE_URL = 'https://www.ebi.ac.uk/proteins/api/'


class ProteinsAPI(Abclog):
    '''
    Implement The Proteins REST API

    * <https://www.ebi.ac.uk/proteins/api/doc/index.html>
    '''
    
    headers = {'Accept': 'application/json'}
    api_sets = frozenset((
                'coordinates', 'coordinates/', 'coordinates/location/',
                'uniparc', 'uniparc/accession/', 'uniparc/best/guess',
                'uniparc/dbreference/', 'uniparc/proteome/', 'uniparc/sequence',  # NOTE: uniparc/sequence use POST method!
                'uniparc/upi/'))
    
    @classmethod
    def get_file_suffix(cls) -> str:
        res = cls.headers["Accept"].split('/')[1]
        assert res in ('json', 'xml', 'x-gff'), f"Unexcepted Case: {cls.headers}"
        return res

    @classmethod
    def dumpsParams(cls, params: Dict) -> str:
        return '&'.join(f'{key}={value}' for key, value in params.items())

    @classmethod
    def task_unit(cls, suffix: str, params: Dict, folder: Path, identifier:Optional[str]=None) -> Tuple:
        args = dict(
            url=f'{BASE_URL}{suffix}' if identifier is None else f'{BASE_URL}{suffix}{identifier}',
            headers=cls.headers,
            params=params)
        return 'get', args, folder/f'{identifier if identifier is not None else cls.dumpsParams(params)}.{cls.get_file_suffix()}'

    @classmethod
    def yieldTasks(cls, suffix: str, params_collection: Iterable[Dict], folder: Path, identifiers: Optional[Iterable[str]]) -> Generator:
        # https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=-1&ensembl=ENST00000554444
        if identifiers is None:
            for params in params_collection:
                yield cls.task_unit(suffix, params, folder)
        else:
            for identifier, params in zip(identifiers, params_collection):
                yield cls.task_unit(suffix, params, folder, identifier)

    @classmethod
    def retrieve(cls, suffix: str, params_collection: Iterable[Dict], folder: Union[Path, str], identifiers: Optional[Iterable[str]] = None, concur_req: int = 20, rate: float = 1.5, ret_res: bool = True, **kwargs):
        assert suffix in cls.api_sets, f"Invalid suffix! Valid set is \n{cls.api_sets}"
        folder = Path(folder)
        res = UnsyncFetch.multi_tasks(
            cls.yieldTasks(suffix, params_collection, folder, identifiers),
            concur_req=concur_req,
            rate=rate,
            logger=cls.logger,
            ret_res=ret_res,
            semaphore=kwargs.get('semaphore', None))
        return res
    
    @classmethod
    def single_retrieve(cls, suffix: str, params:Dict, folder: Union[Path, str], semaphore, identifier:Optional[str]=None, rate: float = 1.5):
        assert suffix in cls.api_sets, f"Invalid suffix! Valid set is \n{cls.api_sets}"
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
        

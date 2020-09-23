# @Created Date: 2020-09-17 12:02:37 am
# @Filename: api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-09-17 12:02:45 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Union, Optional, Iterator, Iterable, Set, Dict, List, Any, Generator, Callable, Tuple
from pdb_profiling.log import Abclog

BASE_URL = 'https://www.ebi.ac.uk/proteins/api/'

class ProteinsAPI(Abclog):
    '''
    Implement The Proteins REST API

    * <https://www.ebi.ac.uk/proteins/api/doc/index.html>
    '''
    
    headers =  {'Content-Type': 'application/json'}
    
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
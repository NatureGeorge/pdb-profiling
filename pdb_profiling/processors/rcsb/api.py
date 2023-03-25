# @Created Date: 2020-12-24 01:21:57 pm
# @Filename: api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2023-03-25 12:18:49 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pdb_profiling.fetcher.webfetch import UnsyncFetch
from hashlib import sha1
from pathlib import Path
from typing import Union


class RCSBDataAPI(object):
    root = 'https://data.rcsb.org/'
    rest_api_root = f'{root}rest/v1/core/'
    graphql_root = f'{root}graphql'
    headers = {'Connection': 'close', 'Content-Type': 'application/json;charset=UTF-8'}
    api_set = frozenset(('entry/', 'assembly/', 'polymer_entity/', 'branched_entity/', 'nonpolymer_entity/',
                         'polymer_entity_instance/', 'branched_entity_instance/', 'nonpolymer_entity_instance/', 'interface/'))
    # TODO: add <https://data.rcsb.org/redoc/index.html#tag/Repository-Holdings-Service>

    @classmethod
    def task_unit(cls, identifier, suffix: str, folder):
        return 'get', dict(url=f'{cls.rest_api_root}{suffix}{identifier}', headers=cls.headers), Path(folder)/f'{identifier.replace("/", "%")}.json'

    @classmethod
    def single_retrieve(cls, identifier: str, suffix: str, folder: Union[Path, str], semaphore, to_do_func=None, rate: float = 1.5):
        return UnsyncFetch.single_task(task=cls.task_unit(identifier, suffix, folder), semaphore=semaphore, to_do_func=to_do_func, rate=rate)

    @classmethod
    def graphql_retrieve(cls, query, folder, semaphore, to_do_func=None, rate: float = 1.5):
        return UnsyncFetch.single_task(task=('get', dict(url=cls.graphql_root, params=dict(query=query), headers=cls.headers), Path(folder)/f'{sha1(bytes(query, encoding="utf-8")).hexdigest()}.json'), semaphore=semaphore, to_do_func=to_do_func, rate=rate)


class RCSBSearchAPI(object):
    root = 'https://search.rcsb.org/rcsbsearch/v1/query'
    headers = {'Connection': 'close', 'Content-Type': 'application/json;charset=UTF-8'}

    @classmethod
    def single_retrieve(cls, query, folder, semaphore, to_do_func=None, rate: float = 1.5):
        return UnsyncFetch.single_task(task=('get', dict(url=cls.root, params=dict(json=query), headers=cls.headers), Path(folder)/f'{sha1(bytes(query, encoding="utf-8")).hexdigest()}.json'), semaphore=semaphore, to_do_func=to_do_func, rate=rate)


class RCSB1DCoordinatesAPI(object):
    root = 'https://1d-coordinates.rcsb.org/'
    graphql_root = f'{root}graphql'

    headers = {'Connection': 'close', 'Content-Type': 'application/json;charset=UTF-8'}

    @classmethod
    def graphql_retrieve(cls, query, folder, semaphore, to_do_func=None, rate: float = 1.5):
        return UnsyncFetch.single_task(task=('get', dict(url=cls.graphql_root, params=dict(query=query), headers=cls.headers), Path(folder)/f'{sha1(bytes(query, encoding="utf-8")).hexdigest()}.json'), semaphore=semaphore, to_do_func=to_do_func, rate=rate)

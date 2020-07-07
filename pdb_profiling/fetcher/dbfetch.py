# @Created Date: 2020-05-25 10:21:37 pm
# @Filename: dbfetch.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-05-26 12:17:51 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from __future__ import absolute_import
import asyncio
from unsync import unsync, Unfuture
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable
from typing import Callable, List, Optional
from pandas import DataFrame, concat

"""
Currently support DB:

* Neo4j
* Sqlite
"""


class Neo4j:
    """
    Neo4j database API

    Add support for asyncio/unsync to neo4j-python-driver (Based on savv's `Neo4j Class`)

    Reference
    ---------
    * <https://github.com/neo4j/neo4j-python-driver/issues/180#issuecomment-380056816>
    """

    # How long to wait after each successive failure.
    RETRY_WAITS: List = [0, 1, 4]

    def __init__(self, config, concur_req: int = 200, semaphore = None, log_func: Callable = print):
        self._config = config
        self.log_func = log_func
        if semaphore is None:
            self._semaphore = asyncio.Semaphore(concur_req)
        else:
            self.log_func("Neo4j.init: pass asyncio.Semaphore")
            self._semaphore = semaphore
        
    
    @unsync
    async def connnect(self) -> Unfuture:
        for retry_wait in self.RETRY_WAITS:
            try:
                await self.init_driver()
                break
            except Exception as e:
                if retry_wait == self.RETRY_WAITS[-1]:
                    raise
                else:
                    self.log_func(f'!!!: retrying to Init DB; err: {e}')
                    # wait for 0, 1, 3... seconds.
                    await asyncio.sleep(retry_wait)
        return self

    @unsync
    def init_driver(self):
        self.log_func("Init Neo4j DataBase Driver")
        self.driver = GraphDatabase.driver(self._config['url'], auth=(self._config['user'], self._config['pass']))

    @unsync
    def afetch_start(self, query, **kwargs):
        with self.driver.session() as session:
            return [dict(i) for i in session.run(query, **kwargs)]
    
    @unsync
    async def afetch(self, query, **kwargs):
        for retry_wait in self.RETRY_WAITS:
            try:
                async with self._semaphore:
                    res = await self.afetch_start(query, **kwargs)
                break
            except (BrokenPipeError, ServiceUnavailable) as e:
                if retry_wait == self.RETRY_WAITS[-1]:
                    raise e
                else:
                    self.log_func('BrokenPipeError or ServiceUnavailable. Retrying...')
                    await asyncio.sleep(retry_wait)
                    self.init_driver().result()
        return res
         
    def close(self):
        """Close Neo4j DataBase Driver"""
        self.log_func("Close Neo4j DataBase Driver")
        self.driver.close()


if __name__ == "__main__":
    import sys
    sys.path.append("C:/GitWorks/pdb-profiling")
    from pdb_profiling.processers.pdbe.neo4j_api import Entry
    config = {'user': '...', 'pass': '...',
              'url': 'bolt://...'}
    queries = [
        Entry.get_residues('1a01', '1', 'A', observed_only=True),
        Entry.get_residues('2xyn', '1', 'A', observed_only=True),
        Entry.get_residues('1z8g', '1', 'A', observed_only=True)
        ]
    demo = Neo4j(config).connnect().result()
    ress = [demo.afetch(query, **kwargs).result() for query, kwargs in queries]
    result = concat((DataFrame(dict(i) for i in res) for res in ress), ignore_index=True, sort=False)
    print(result)
    demo.close()

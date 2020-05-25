import asyncio
from unsync import unsync
from neo4j import GraphDatabase, basic_auth, READ_ACCESS
from neo4j.exceptions import ServiceUnavailable
from typing import Callable, List, Optional
from pandas import DataFrame, concat

"""
Currently support DB:

* Neo4j
* ...
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

    def __init__(self, config, concur_req: int = 100, log_func: Callable = print):
        self._config = config
        self._semaphore = asyncio.Semaphore(concur_req)
        self.log_func = log_func
    
    @unsync
    async def connnect(self):
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
    
    @unsync
    def init_driver(self):
        self.log_func("Init Neo4j DataBase Driver")
        auth = basic_auth(self._config['user'], self._config['pass'])
        self.driver = GraphDatabase.driver(self._config['url'], auth=auth, encrypted=True)

    @unsync
    def afetch_start(self, query, **kwargs):
        with self.driver.session(access_mode=READ_ACCESS) as session:
            return session.run(query, **kwargs).records()
    
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
                    await self.init_driver()
        return res
         
    def close(self):
        """Close Neo4j DataBase Driver"""
        self.log_func("Close Neo4j DataBase Driver")
        self.driver.close()


if __name__ == "__main__":
    config = {}
    queries = [
        "MATCH(entry: Entry{ID: '5akd'})-[:HAS_ENTITY] -> (entity: Entity{POLYMER_TYPE: 'P'})-[:HAS_RESIDUE_CONFLICT] -> (resCon: ResidueConflict) RETURN entry.ID, entity.ID, resCon.DETAILS, resCon.ID",
        "MATCH(entry: Entry{ID: '1a01'})-[:HAS_ENTITY] -> (entity: Entity{POLYMER_TYPE: 'P'})-[:HAS_RESIDUE_CONFLICT] -> (resCon: ResidueConflict) RETURN entry.ID, entity.ID, resCon.DETAILS, resCon.ID",
        "MATCH(entry: Entry{ID: '2xyn'})-[:HAS_ENTITY] -> (entity: Entity{POLYMER_TYPE: 'P'})-[:HAS_RESIDUE_CONFLICT] -> (resCon: ResidueConflict) RETURN entry.ID, entity.ID, resCon.DETAILS, resCon.ID"
        ]
    demo = Neo4j(config)
    demo.connnect().result()
    ress = [demo.afetch(query).result() for query in queries]
    result = concat((DataFrame(dict(i) for i in res) for res in ress), ignore_index=True, sort=False)
    print(result)
    demo.close()

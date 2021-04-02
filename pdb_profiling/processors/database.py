# @Created Date: 2020-09-17 12:02:31 am
# @Filename: database.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-09-27 03:18:19 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from asyncio import sleep as asyncio_sleep
import databases
import orm
import sqlalchemy
from sqlite3 import OperationalError, InterfaceError
from unsync import unsync
from typing import Dict, Iterable
from pdb_profiling.log import Abclog
from random import uniform
# from queue import Queue


class SqliteDB(Abclog):

    def init_table_model(self):
        pass

    def __init__(self, url: str, drop_all: bool = False):
        self.metadata = sqlalchemy.MetaData()
        self.database = databases.Database(url)
        self.engine = sqlalchemy.create_engine(url)
        self.engine.execute("PRAGMA journal_mode=WAL")  # auto_vacuum=FULL
        self.init_table_model()
        if drop_all:
            self.metadata.drop_all(self.engine, checkfirst=True)
        self.metadata.create_all(self.engine, checkfirst=True)
        # self.queue = Queue()

    def close(self):
        """TODO: make it real"""
        self.engine.dispose()

    def sync_insert(self, table, values: Iterable[Dict], prefix_with: str = "OR IGNORE"):
        try:
            self.engine.execute(
                table.__table__.insert().prefix_with(prefix_with),
                values)
        except InterfaceError as e:
            self.logger.error(f"{e}\n{values}")
            raise e

    '''
    @unsync
    async def async_insert_queue(self, table, values: Iterable[Dict], prefix_with: str = "OR IGNORE", insert_sleep: float = 30.5):
        self.queue.put_nowait(self.database.execute(table.__table__.insert().values(values).prefix_with(prefix_with)))
        now_task = self.queue.get_nowait()
        await now_task
        self.queue.task_done()
    '''

    @unsync
    async def async_insert_chunk(self, table, values: Iterable[Dict], prefix_with: str = "OR IGNORE", insert_sleep_range=(10, 30), chunksize=10000):
        for i in range(0, len(values), chunksize):
            await self.async_insert(table=table, values=values[i:i+chunksize], prefix_with=prefix_with, insert_sleep_range=insert_sleep_range)

    @unsync
    async def async_insert(self, table, values: Iterable[Dict], prefix_with: str = "OR IGNORE", insert_sleep_range=(10,30)):
        while True:
            try:
                query = table.__table__.insert().values(values).prefix_with(prefix_with)
                await self.database.execute(query)
                break
            except OperationalError as e:
                insert_sleep = uniform(*insert_sleep_range)
                self.logger.error(f"OperationalError: sleep {insert_sleep}s and try again")
                await asyncio_sleep(insert_sleep)
            except InterfaceError as e:
                self.logger.error(f"{e}\n{values}")
                raise e

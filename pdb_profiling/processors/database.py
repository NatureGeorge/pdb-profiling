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
from sqlite3 import OperationalError
from unsync import unsync
from typing import Dict, Iterable
from pdb_profiling.log import Abclog


class SqliteDB(Abclog):

    def init_table_model(self):
        pass

    def __init__(self, url: str, drop_all: bool = False, insert_sleep: float = 45.5):
        self.metadata = sqlalchemy.MetaData()
        self.insert_sleep = insert_sleep
        self.database = databases.Database(url)
        self.engine = sqlalchemy.create_engine(url)
        self.engine.execute("PRAGMA journal_mode=WAL")
        self.init_table_model()
        if drop_all:
            self.metadata.drop_all(self.engine, checkfirst=True)
        self.metadata.create_all(self.engine, checkfirst=True)

    def close(self):
        """TODO: make it real"""
        self.engine.dispose()

    def sync_insert(self, table, values: Iterable[Dict], prefix_with: str = "OR IGNORE"):
        self.engine.execute(
            table.__table__.insert().prefix_with(prefix_with),
            values)

    @unsync
    async def async_insert(self, table, values: Iterable[Dict], prefix_with: str = "OR IGNORE"):
        while True:
            try:
                query = table.__table__.insert().values(values).prefix_with(prefix_with)
                await self.database.execute(query)
                break
            except OperationalError as e:
                self.logger.error(
                    f"{e}, sleep {self.insert_sleep}s and try again")
                await asyncio_sleep(self.insert_sleep)

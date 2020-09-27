# @Created Date: 2020-09-17 12:02:31 am
# @Filename: __init__.py
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
from pdb_profiling.processors.database import SqliteDB


class ProteinsDB(SqliteDB):

    def init_table_model(self):
        class ALTERNATIVE_PRODUCTS(orm.Model):
            __tablename__ = 'ALTERNATIVE_PRODUCTS'
            __metadata__ = self.metadata
            __database__ = self.database
            isoform = orm.String(max_length=50, primary_key=True)
            name = orm.String(max_length=50)
            synonyms = orm.Text(allow_null=True, allow_blank=True)
            sequenceStatus = orm.String(max_length=50)
            sequence = orm.JSON(allow_null=True)
            Entry = orm.String(max_length=50, primary_key=True)
        
        class DB_REFERENCES(orm.Model):
            __tablename__ = 'dbReferences'
            __metadata__ = self.metadata
            __database__ = self.database
            type = orm.String(max_length=100)
            isoform = orm.String(max_length=50, allow_null=True)
            Entry = orm.String(max_length=50, primary_key=True)
            protein = orm.String(max_length=50, primary_key=True)
            transcript = orm.String(max_length=50)
            gene = orm.String(max_length=50, allow_null=True)

        self.DB_REFERENCES = DB_REFERENCES
        self.ALTERNATIVE_PRODUCTS = ALTERNATIVE_PRODUCTS
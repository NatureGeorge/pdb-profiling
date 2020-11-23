# @Created Date: 2020-11-23 10:29:05 am
# @Filename: __init__.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-11-23 10:29:40 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import orm
from pdb_profiling.processors.database import SqliteDB


class CustomDB(SqliteDB):

    def init_table_model(self):
        class Mutation(orm.Model):
            __tablename__ = 'Mutation'
            __metadata__ = self.metadata
            __database__ = self.database
            ftId = orm.String(max_length=20, primary_key=True)
            Ref = orm.String(max_length=3, primary_key=True)
            Pos = orm.Integer(primary_key=True)
            Alt = orm.String(max_length=3, primary_key=True)

        class IDMapping(orm.Model):
            __tablename__ = 'IDMapping'
            __metadata__ = self.metadata
            __database__ = self.database
            ftId = orm.String(max_length=20, primary_key=True)
            Entry = orm.String(max_length=50, primary_key=True)
            isoform = orm.String(max_length=50, primary_key=True)
            is_canonical = orm.Boolean()


        self.Mutation = Mutation
        self.IDMapping = IDMapping

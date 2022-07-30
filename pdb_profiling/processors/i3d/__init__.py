# @Created Date: 2020-10-17 08:02:13 am
# @Filename: __init__.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-10-17 08:02:20 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import orm
from pdb_profiling.processors.database import SqliteDB


class I3DDB(SqliteDB):

    def init_table_model(self):

        class InteractionMeta(orm.Model):
            __tablename__ = 'InteractionMeta'
            __metadata__ = self.metadata
            __database__ = self.database
            pdb_id = orm.String(max_length=4, primary_key=True)
            Entry_1 = orm.String(max_length=100, primary_key=True)
            Entry_2 = orm.String(max_length=100, primary_key=True)
            assembly_id = orm.Integer(primary_key=True)
            chain_id_1 = orm.String(max_length=4, primary_key=True)
            model_id_1 = orm.Integer(primary_key=True)
            chain_id_2 = orm.String(max_length=4, primary_key=True)
            model_id_2 = orm.Integer(primary_key=True)
            organism = orm.String(max_length=100)
            interaction_type = orm.String(max_length=2)

        self.InteractionMeta = InteractionMeta

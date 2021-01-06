# @Created Date: 2020-12-24 01:21:44 pm
# @Filename: __init__.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-12-24 01:28:45 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pdb_profiling.processors.database import SqliteDB
import orm


class RCSBDB(SqliteDB):

    def init_table_model(self):
        class rcsb_cluster_membership(orm.Model):
            __tablename__ = 'rcsb_cluster_membership'
            __metadata__ = self.metadata
            __database__ = self.database
            rcsb_id = orm.Text()
            score = orm.Float()
            identity = orm.Integer(primary_key=True)
            cluster_id = orm.Integer(primary_key=True)
            pdb_id = orm.String(max_length=4, primary_key=True)
            entity_id = orm.Integer(primary_key=True)

        self.rcsb_cluster_membership = rcsb_cluster_membership

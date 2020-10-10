# @Created Date: 2020-08-12 11:28:01 pm
# @Filename: __init__.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-10-10 10:37:46 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import orm
from unsync import unsync
from pdb_profiling.processors.database import SqliteDB


class PDBeDB(SqliteDB):

    def init_table_model(self):
        class ResidueMapping(orm.Model):
            __tablename__ = 'ResidueMapping'
            __metadata__ = self.metadata
            __database__ = self.database
            UniProt = orm.String(max_length=50, primary_key=True)
            author_insertion_code = orm.String(max_length=50, allow_null=True, allow_blank=True, default='')
            author_residue_number = orm.Integer()
            chain_id = orm.String(max_length=10)
            entity_id = orm.Integer(primary_key=True)
            identifier = orm.Text(allow_null=True, allow_blank=True, default='')
            name = orm.Text(allow_null=True, allow_blank=True, default='')
            observed = orm.String(max_length=1)
            pdb_id = orm.String(max_length=4, primary_key=True)
            pdb_one_letter_code = orm.String(max_length=1, primary_key=True)
            residue_number = orm.Integer(primary_key=True)
            struct_asym_id = orm.String(max_length=10, primary_key=True)
            unp_one_letter_code = orm.String(max_length=1, primary_key=True)
            unp_residue_number = orm.Integer(primary_key=True)

        self.ResidueMapping = ResidueMapping

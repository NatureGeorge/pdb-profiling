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
        
        class ResidueMapping(orm.Model):
            __tablename__ = 'ResidueMapping'
            __metadata__ = self.metadata
            __database__ = self.database
            UniProt = orm.String(max_length=50, primary_key=True)
            author_insertion_code = orm.String(max_length=50, allow_null=True, allow_blank=True, default='')
            author_residue_number = orm.Integer()
            chain_id = orm.String(max_length=10)
            struct_asym_id = orm.String(max_length=10, primary_key=True)
            entity_id = orm.Integer(primary_key=True)
            pdb_id = orm.String(max_length=4, primary_key=True)
            residue_number = orm.Integer(primary_key=True)
            unp_residue_number = orm.Integer(primary_key=True)
            residue_name = orm.String(max_length=10)
            observed_ratio = orm.Float()
            multiple_conformers = orm.JSON(allow_null=True)
            conflict_code = orm.String(max_length=3, allow_null=True)


        self.Mutation = Mutation
        self.IDMapping = IDMapping
        self.ResidueMapping = ResidueMapping

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
        class AAThree2one(orm.Model):
            __tablename__ = 'AAThree2one'
            __metadata__ = self.metadata
            __database__ = self.database
            three_letter_code = orm.String(max_length=3, primary_key=True)
            one_letter_code = orm.String(max_length=1, primary_key=True)
        
        class UniProtSeq(orm.Model):
            __tablename__ = 'UniProtSeq'
            __metadata__ = self.metadata
            __database__ = self.database
            isoform = orm.String(max_length=20, primary_key=True)
            Pos = orm.Integer(primary_key=True)
            Ref = orm.String(max_length=1, primary_key=True)
        
        """class PDBSeq(orm.Model):
            __tablename__ = 'PDBSeq'
            __metadata__ = self.metadata
            __database__ = self.database
            pdb_id = orm.String(max_length=20, primary_key=True)
            entity_id = orm.Integer(primary_key=True)
            chain_id = orm.String(max_length=10)
            struct_asym_id = orm.String(max_length=10, primary_key=True)
            residue_number = orm.Integer(primary_key=True)
            author_residue_number = orm.Integer(primary_key=True)
            author_insertion_code = orm.String(primary_key=True, max_length=3, allow_blank=True, default='')
            Ref = orm.String(max_length=1, primary_key=True)"""

        class Mutation(orm.Model):
            __tablename__ = 'Mutation'
            __metadata__ = self.metadata
            __database__ = self.database
            ftId = orm.String(max_length=20, primary_key=True)
            Ref = orm.String(max_length=3, primary_key=True)
            Pos = orm.Integer(primary_key=True)
            Alt = orm.String(max_length=3, primary_key=True)

        class PDBMutation(orm.Model):
            __tablename__ = 'PDBMutation'
            __metadata__ = self.metadata
            __database__ = self.database
            pdb_id = orm.String(max_length=4, primary_key=True)
            entity_id = orm.Integer(allow_null=True)
            chain_id = orm.String(max_length=10, allow_null=True)
            struct_asym_id = orm.String(max_length=10, allow_null=True)
            Ref = orm.String(max_length=1, primary_key=True)
            Alt = orm.String(max_length=1, primary_key=True)
            residue_number = orm.Integer(primary_key=True)

        class PDBAuthMutation(orm.Model):
            __tablename__ = 'PDBAuthMutation'
            __metadata__ = self.metadata
            __database__ = self.database
            pdb_id = orm.String(max_length=4, primary_key=True)
            entity_id = orm.Integer(allow_null=True)
            chain_id = orm.String(max_length=10, allow_null=True)
            struct_asym_id = orm.String(max_length=10, allow_null=True)
            Ref = orm.String(max_length=1, primary_key=True)
            Alt = orm.String(max_length=1, primary_key=True)
            author_residue_number = orm.Integer(primary_key=True)
            author_insertion_code = orm.String(primary_key=True, max_length=3, allow_blank=True, default='')

        class IDMapping(orm.Model):
            __tablename__ = 'IDMapping'
            __metadata__ = self.metadata
            __database__ = self.database
            ftId = orm.String(max_length=20, primary_key=True)
            Entry = orm.String(max_length=50, primary_key=True)
            isoform = orm.String(max_length=50, primary_key=True)
            is_canonical = orm.Boolean()
        
        """class ResidueMapping(orm.Model):
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
            conflict_code = orm.String(max_length=3, allow_null=True)"""

        class ResidueMappingRange(orm.Model):
            __tablename__ = 'ResidueMappingRange'
            __metadata__ = self.metadata
            __database__ = self.database
            UniProt = orm.String(max_length=20, primary_key=True)
            pdb_id = orm.String(max_length=4, primary_key=True)
            entity_id = orm.Integer(primary_key=True)
            chain_id = orm.String(max_length=10)
            struct_asym_id = orm.String(max_length=10, primary_key=True)
            auth_pdb_beg = orm.Integer()
            auth_pdb_end = orm.Integer()
            pdb_beg = orm.Integer(primary_key=True)
            pdb_end = orm.Integer(primary_key=True)
            unp_beg = orm.Integer(primary_key=True)
            unp_end = orm.Integer(primary_key=True)
            # Special:
            author_insertion_code = orm.String(max_length=10, allow_null=True, allow_blank=True, default='')
            residue_name = orm.String(max_length=10, allow_null=True, allow_blank=True, default='')
            observed_ratio = orm.Float()
            multiple_conformers = orm.JSON(allow_null=True)
            conflict_code = orm.String(max_length=3, allow_null=True)

        class SelectedMappingMeta(orm.Model):
            __tablename__ = 'SelectedMappingMeta'
            __metadata__ = self.metadata
            __database__ = self.database
            UniProt = orm.String(max_length=50, primary_key=True)
            pdb_id = orm.String(max_length=4, primary_key=True)
            entity_id = orm.Integer(primary_key=True)
            struct_asym_id = orm.String(max_length=10, primary_key=True)
            chain_id = orm.String(max_length=10)
            bs_score = orm.Float()
            select_rank = orm.Integer()
            select_tag = orm.Boolean()
            after_select_rank = orm.Integer()

        class ResidueAnnotation(orm.Model):
            __tablename__ = 'ResidueAnnotation'
            __metadata__ = self.metadata
            __database__ = self.database
            pdb_id = orm.String(max_length=4, primary_key=True)
            entity_id = orm.Integer(primary_key=True)
            struct_asym_id = orm.String(max_length=10, primary_key=True)
            chain_id = orm.String(max_length=10)
            resource = orm.String(max_length=100, primary_key=True)
            resource_id = orm.String(max_length=200, primary_key=True)
            pdb_beg = orm.Integer()
            pdb_end = orm.Integer()

        class UniProtAnnotation(orm.Model):
            __tablename__ = 'UniProtAnnotation'
            __metadata__ = self.metadata
            __database__ = self.database
            UniProt = orm.String(max_length=50, primary_key=True)
            resource = orm.String(max_length=100, primary_key=True)
            resource_id = orm.String(max_length=200, primary_key=True)
            unp_beg = orm.Integer()
            unp_end = orm.Integer()
        
        class SMRModel(orm.Model):
            __tablename__ = 'SMRModel'
            __metadata__ = self.metadata
            __database__ = self.database
            UniProt = orm.String(max_length=50, primary_key=True)
            coordinates = orm.String(max_length=500, primary_key=True)
            unp_beg = orm.Integer()
            unp_end = orm.Integer()
            identity = orm.Float()
            similarity = orm.Float()
            coverage = orm.Float()
            oligo_state = orm.String(max_length=50)
            with_ligand = orm.Boolean()
            select_rank = orm.Integer()
            select_tag = orm.Boolean()
        
        class MappedMutation(orm.Model):
            __tablename__ = 'MappedMutation'
            __metadata__ = self.metadata
            __database__ = self.database
            UniProt = orm.String(max_length=50, primary_key=True)
            Ref = orm.String(max_length=3, primary_key=True)
            Pos = orm.Integer(primary_key=True)
            Alt = orm.String(max_length=3, primary_key=True)


        class PI(orm.Model):
            __tablename__ = 'PI'
            __metadata__ = self.metadata
            __database__ = self.database

            UniProt = orm.String(max_length=50, primary_key=True)
            pdb_id = orm.String(max_length=4, primary_key=True)
            entity_id = orm.Integer(primary_key=True)
            struct_asym_id = orm.String(max_length=10, primary_key=True)
            chain_id = orm.String(max_length=10)
            assembly_id = orm.Integer(primary_key=True)
            model_id = orm.Integer()
            struct_asym_id_in_assembly = orm.String(max_length=10, primary_key=True)
            interface_id = orm.Integer(primary_key=True)
            css = orm.Float()
            i_select_tag = orm.Boolean()
            i_select_rank = orm.Integer()
            pdb_beg = orm.Integer()
            pdb_end = orm.Integer()


        self.AAThree2one = AAThree2one
        self.UniProtSeq = UniProtSeq
        self.Mutation = Mutation
        self.IDMapping = IDMapping
        self.UniProtAnnotation = UniProtAnnotation
        self.ResidueMappingRange = ResidueMappingRange
        self.SelectedMappingMeta = SelectedMappingMeta
        self.ResidueAnnotation = ResidueAnnotation
        self.SMRModel = SMRModel
        self.MappedMutation = MappedMutation
        self.PI = PI
        self.PDBMutation = PDBMutation
        self.PDBAuthMutation = PDBAuthMutation
        #self.PDBSeq = PDBSeq

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

        class StatsProteinEntitySeq(orm.Model):
            __tablename__ = 'StatsProteinEntitySeq'
            __metadata__ = self.metadata
            __database__ = self.database
            pdb_id = orm.String(max_length=4, primary_key=True)
            molecule_type = orm.Text()
            entity_id = orm.Integer(primary_key=True)
            ca_p_only = orm.Boolean()
            SEQRES_COUNT = orm.Integer()
            STD_INDEX = orm.JSON()
            STD_COUNT = orm.Integer()
            NON_INDEX = orm.JSON()
            NON_COUNT = orm.Integer()
            UNK_INDEX = orm.JSON()
            UNK_COUNT = orm.Integer()
            ARTIFACT_INDEX = orm.JSON()
        
        class StatsNucleotideEntitySeq(orm.Model):
            __tablename__ = 'StatsNucleotideEntitySeq'
            __metadata__ = self.metadata
            __database__ = self.database
            pdb_id = orm.String(max_length=4, primary_key=True)
            molecule_type = orm.Text()
            entity_id = orm.Integer(primary_key=True)
            ca_p_only = orm.Boolean()
            SEQRES_COUNT = orm.Integer()
            dNTP_INDEX = orm.JSON()
            dNTP_COUNT = orm.Integer()
            NTP_INDEX = orm.JSON()
            NTP_COUNT = orm.Integer()
            NON_INDEX = orm.JSON()
            NON_COUNT = orm.Integer()
            UNK_INDEX = orm.JSON()
            UNK_COUNT = orm.Integer()

        class StatsChainSeq(orm.Model):
            __tablename__ = 'StatsChainSeq'
            __metadata__ = self.metadata
            __database__ = self.database
            pdb_id = orm.String(max_length=4, primary_key=True)
            entity_id = orm.Integer(primary_key=True)
            chain_id = orm.String(max_length=10)
            struct_asym_id = orm.String(max_length=10, primary_key=True)
            OBS_INDEX = orm.JSON()
            OBS_COUNT = orm.Integer()
            OBS_RATIO_SUM = orm.Float()
            BINDING_LIGAND_INDEX = orm.JSON()
            BINDING_LIGAND_COUNT = orm.Integer()

        self.ResidueMapping = ResidueMapping
        self.StatsProteinEntitySeq = StatsProteinEntitySeq
        self.StatsNucleotideEntitySeq = StatsNucleotideEntitySeq
        self.StatsChainSeq = StatsChainSeq

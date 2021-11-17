# @Created Date: 2020-09-17 12:02:31 am
# @Filename: __init__.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-09-27 03:18:19 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import orm
from pdb_profiling.processors.database import SqliteDB


class ProteinsDB(SqliteDB):

    def init_table_model(self):
        class ALTERNATIVE_PRODUCTS(orm.Model):
            __tablename__ = 'ALTERNATIVE_PRODUCTS'
            __metadata__ = self.metadata
            __database__ = self.database
            isoform = orm.String(max_length=50, primary_key=True)
            name = orm.JSON(allow_null=True)
            synonyms = orm.JSON(allow_null=True)
            sequenceStatus = orm.String(max_length=50)
            VAR_SEQ = orm.JSON(allow_null=True)
            accession = orm.String(max_length=50, primary_key=True)
            else_iso = orm.JSON(allow_null=True)
            iso_range = orm.JSON(allow_null=True)
            iso_range_len = orm.Integer(allow_null=True)
            sequence = orm.Text(allow_null=True)
            length = orm.Integer(allow_null=True)
        
        class INTERACTION(orm.Model):
            __tablename__ = 'INTERACTION'
            __metadata__ = self.metadata
            __database__ = self.database
            accession1 = orm.String(max_length=100, primary_key=True)
            accession2 = orm.String(max_length=100, primary_key=True)
            gene = orm.String(max_length=100, allow_null=True, allow_blank=True)
            interactor1 = orm.String(max_length=100, primary_key=True)
            interactor2 = orm.String(max_length=100, primary_key=True)
            organismDiffer = orm.Boolean()
            experiments = orm.Integer()
            chain1 = orm.Text(allow_null=True, allow_blank=True)
            chain2 = orm.Text(allow_null=True, allow_blank=True)
            accession = orm.String(max_length=50, primary_key=True)

        class DB_REFERENCES(orm.Model):
            __tablename__ = 'dbReferences'
            __metadata__ = self.metadata
            __database__ = self.database
            type = orm.String(max_length=100)
            isoform = orm.String(max_length=50, allow_null=True, allow_blank=True)
            accession = orm.String(max_length=50, primary_key=True)
            protein = orm.String(max_length=50, primary_key=True)
            transcript = orm.String(max_length=50)
            gene = orm.String(max_length=50, allow_null=True, allow_blank=True)

        class OTHER_DB_REFERENCES(orm.Model):
            __tablename__ = 'other_dbReferences'
            __metadata__ = self.metadata
            __database__ = self.database
            type = orm.String(primary_key=True, max_length=100)
       	    id = orm.String(primary_key=True, max_length=100)
            properties = orm.JSON(allow_null=True)
            isoform = orm.String(primary_key=True, max_length=100, allow_blank=True)
            evidences = orm.JSON(allow_null=True)
            accession = orm.String(primary_key=True, max_length=100)

        class FEATURES(orm.Model):
            __tablename__ = 'features'
            __metadata__ = self.metadata
            __database__ = self.database
            type = orm.String(primary_key=True, max_length=100)
            category = orm.String(primary_key=True, max_length=100)
            description = orm.Text(allow_null=True, allow_blank=True)
            begin = orm.Integer(primary_key=True)
            end = orm.Integer(primary_key=True)
            molecule = orm.String(primary_key=True, max_length=100, allow_blank=True)
            ftId = orm.String(max_length=100, allow_null=True, allow_blank=True)
            evidences = orm.JSON(allow_null=True)
            alternativeSequence = orm.Text(allow_null=True, allow_blank=True)
            accession = orm.String(primary_key=True, max_length=100)

        class INFO(orm.Model):
            __tablename__ = 'INFO'
            __metadata__ = self.metadata
            __database__ = self.database
            accession = orm.String(max_length=100, primary_key=True)
            id = orm.String(max_length=100)
            proteinExistence = orm.Text(allow_null=True, allow_blank=True)
            info = orm.JSON()
            organism = orm.JSON()
            secondaryAccession = orm.JSON(allow_null=True)
            protein = orm.JSON()
            gene = orm.JSON()
            sequence = orm.Text()
            length = orm.Integer()
        
        """class CoordinateMapping(orm.Model):
            __tablename__ = 'CoordinateMapping'
            __metadata__ = self.metadata
            __database__ = self.database
            taxid = orm.Integer(primary_key=True)
            chromosome = orm.String(max_length=20, primary_key=True)
            geneStart = orm.Integer(primary_key=True)
            geneEnd = orm.Integer(primary_key=True)
            proteinStart = orm.Integer(primary_key=True)
            proteinEnd = orm.Integer(primary_key=True)
            UniProt = orm.String(max_length=50, primary_key=True)"""

        self.DB_REFERENCES = DB_REFERENCES
        self.OTHER_DB_REFERENCES = OTHER_DB_REFERENCES
        self.ALTERNATIVE_PRODUCTS = ALTERNATIVE_PRODUCTS
        self.FEATURES = FEATURES
        self.INTERACTION = INTERACTION
        self.INFO = INFO
        #self.CoordinateMapping = CoordinateMapping

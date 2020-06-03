# @Created Date: 2020-05-26 08:02:13 pm
# @Filename: sqlite_api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-05-26 08:02:16 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import asyncio
import databases
import orm
import sqlalchemy
from time import perf_counter
import pandas as pd
from unsync import unsync
from typing import Dict, Iterable

converters = {
    "pdb_id": str,
    "entity_id": str,
    "chain_id": str,
    "PDB_REV_DATE_ORIGINAL": str,
    "FIRST_REV_DATE": str,
    "PDB_REV_DATE": str,
    "REVISION_DATE": str,
    "resolution": float,
    "METHOD_CLASS": str,
    "BOUND_LIGAND_COUNT": int,
    "BOUND_MOL_COUNT": int,
    "nucleotides_entity_type": str,
    "has_hybrid_nucleotides": bool,
    # SEQRES_Info
    "SEQRES_COUNT": int,
    "AVG_OBS_RATIO": float,
    "AVG_OBS_OBS_RATIO": float,
    "NON_INDEX": str,
    "UNK_INDEX": str,
    "MIS_INDEX": str,
    "UNK_COUNT": int,
    "PURE_SEQRES_COUNT": int,
    "OBS_RECORD_COUNT": int,
    "OBS_UNK_COUNT": int,
    "ATOM_RECORD_COUNT": int,
    # SIFTS_Info
    "UniProt": str,
    "identity": float,
    "pdb_range": str,
    "unp_range": str,
    "group_info": int,
    "pdb_gap_list": str,
    "unp_gap_list": str,
    "var_list": str,
    "repeated": bool,
    "var_0_count": int,
    "unp_gap_0_count": int,
    "unp_pdb_var": int,
    "sifts_range_tag": str,
    "new_unp_range": str,
    "new_pdb_range": str,
    # Res_Info
    "residue_number": int,
    "residue_name": str,
    "obs_ratio": float,
    "author_residue_number": int,
    "author_insertion_code": str,
    # Site_Info
    "from_id": str,
    "Ref": str,
    "Pos": int,
    "Alt": str
    }


class Sqlite_API(object):

    metadata = sqlalchemy.MetaData()

    def init_table_model(self):
        class Entry_Info(orm.Model):
            __tablename__ = 'entry_info'
            __metadata__ = self.metadata
            __database__ = self.database
            pdb_id = orm.String(max_length=4, primary_key=True)
            PDB_REV_DATE_ORIGINAL = orm.String(max_length=19)
            FIRST_REV_DATE = orm.String(max_length=19)
            PDB_REV_DATE = orm.String(max_length=19)
            REVISION_DATE = orm.String(max_length=19)
            resolution = orm.Float()
            METHOD_CLASS = orm.String(max_length=5)
            BOUND_LIGAND_COUNT = orm.Integer()
            BOUND_MOL_COUNT = orm.Integer()
            nucleotides_entity_type = orm.String(
                max_length=100, allow_null=True, allow_blank=True, default='')
            has_hybrid_nucleotides = orm.Boolean()

        class SEQRES_Info(orm.Model):
            __tablename__ = 'seqres_info'
            __metadata__ = self.metadata
            __database__ = self.database
            pdb_id = orm.String(max_length=4, primary_key=True)
            entity_id = orm.String(max_length=4, primary_key=True)
            chain_id = orm.String(max_length=4, primary_key=True)
            SEQRES_COUNT = orm.Integer()
            AVG_OBS_RATIO = orm.Float()
            AVG_OBS_OBS_RATIO = orm.Float()
            NON_INDEX = orm.Text(allow_null=True, allow_blank=True, default='')
            UNK_INDEX = orm.Text(allow_null=True, allow_blank=True, default='')
            MIS_INDEX = orm.Text(allow_null=True, allow_blank=True, default='')
            UNK_COUNT = orm.Integer()
            PURE_SEQRES_COUNT = orm.Integer()
            OBS_RECORD_COUNT = orm.Integer()
            OBS_UNK_COUNT = orm.Integer()
            ATOM_RECORD_COUNT = orm.Integer()

        class SIFTS_Info(orm.Model):
            __tablename__ = 'sifts_info'
            __metadata__ = self.metadata
            __database__ = self.database
            UniProt = orm.String(max_length=20, primary_key=True)
            pdb_id = orm.String(max_length=4, primary_key=True)
            entity_id = orm.String(max_length=4, primary_key=True)
            chain_id = orm.String(max_length=4, primary_key=True)
            identity = orm.Float()
            pdb_range = orm.Text()
            unp_range = orm.Text()
            group_info = orm.Integer()
            pdb_gap_list = orm.Text()
            unp_gap_list = orm.Text()
            var_list = orm.Text()
            repeated = orm.Boolean()
            var_0_count = orm.Integer()
            unp_gap_0_count = orm.Integer()
            unp_pdb_var = orm.Integer()
            sifts_range_tag = orm.String(max_length=25)
            new_unp_range = orm.Text()
            new_pdb_range = orm.Text()

        class PDBRes_Info(orm.Model):
            '''
            entry.ID as pdb_id, 
            entity.ID as entity_id, 
            chain.AUTH_ASYM_ID as chain_id, 
            res.CHEM_COMP_ID as residue_name, 
            toInteger(res.ID) as residue_number, 
            tofloat(inChain.OBSERVED_RATIO) as obs_ratio, 
            toInteger(inChain.AUTH_SEQ_ID) as author_residue_number, 
            inChain.PDB_INS_CODE as author_insertion_code ORDER BY residue_number
            '''
            __tablename__ = 'pdbres_info'
            __metadata__ = self.metadata
            __database__ = self.database
            pdb_id = orm.String(max_length=4, primary_key=True)
            entity_id = orm.String(max_length=4, primary_key=True)
            chain_id = orm.String(max_length=4, primary_key=True)
            residue_number = orm.Integer(primary_key=True)
            residue_name = orm.String(max_length=10)
            obs_ratio = orm.Float()
            author_residue_number = orm.Integer()
            author_insertion_code = orm.String(
                max_length=4, allow_null=True, allow_blank=True, default='')
            
        class Site_Info(orm.Model):
            __tablename__ = 'site_info'
            __metadata__ = self.metadata
            __database__ = self.database
            from_id = orm.String(max_length=20, primary_key=True)
            Ref = orm.String(max_length=3, primary_key=True)
            Pos = orm.Integer(primary_key=True)
            Alt = orm.String(max_length=3, primary_key=True)
            # code = orm.String(max_length=3, primary_key=True, allow_blank=True, allow_null=True, default='')
    
        self.Entry_Info = Entry_Info
        self.SEQRES_Info = SEQRES_Info
        self.SIFTS_Info = SIFTS_Info
        self.PDBRes_Info = PDBRes_Info
        self.Site_Info = Site_Info

    def __init__(self, url: str, drop_all: bool=False):
        self.database = databases.Database(url)
        self.engine = sqlalchemy.create_engine(url)
        self.init_table_model()
        if drop_all:
            self.metadata.drop_all(self.engine, checkfirst=True)
        self.metadata.create_all(self.engine, checkfirst=True)

    def sync_insert(self, table, values: Iterable[Dict], prefix_with: str = "OR IGNORE"):
        self.engine.execute(
            table.__table__.insert().prefix_with(prefix_with),
            values)
    
    @unsync
    async def async_insert(self, table, values: Iterable[Dict], prefix_with: str = "OR IGNORE"):
        '''
        # BAD CODE
        await self.database.execute_many(
            query=table.__table__.insert().prefix_with(prefix_with),
            values=values)
        '''
        query = table.__table__.insert().values(values).prefix_with(prefix_with)
        await self.database.execute(query)


@unsync
async def main():
    start = perf_counter()
    # sqlite_api = Sqlite_API("sqlite:///../../../test/db/orm_db_test.db")
    sqlite_api = Sqlite_API("sqlite:///C:\\Download\\20200601\\DB\\local_sqlite_demo.db")
    print('init db: {:.5f}s'.format(perf_counter()-start))
    '''
    for index, (Info, file) in enumerate(zip((sqlite_api.Entry_Info, sqlite_api.SEQRES_Info, sqlite_api.SIFTS_Info, sqlite_api.PDBRes_Info), (  #
            "C:/Download/20200525/LHY/entry_info.tsv",
            "C:/Download/20200525/LHY/seqres_info.tsv",
            "C:/Download/20200525/LHY/sifts_mapping.tsv",
            "C:/Download/20200525/LHY/pdbres_info.tsv"))):
        if index:
            start = perf_counter()
            values = pd.read_csv(file, sep="\t",
                             converters=converters).drop_duplicates().to_dict('records')
        else:
            start = perf_counter()
            values = pd.read_csv(file, sep="\t",
                                 converters=converters).drop_duplicates(subset=['pdb_id'], keep='last').to_dict('records')
        print('init read_csv: {:.5f}s'.format(perf_counter()-start))
        # query = Info.__table__.insert().prefix_with("OR IGNORE")
        start = perf_counter()
        # await database.execute_many(query=query, values=values)
        # engine.execute(query, values)
        sqlite_api.sync_insert(Info, values)
        print('init insert: {:.5f}s'.format(perf_counter()-start))
    '''
    '''
    site_df = pd.read_csv(
        r'C:\Download\20200525\LHY\site_info_demo.tsv', sep='\t', converters=converters)
    start = perf_counter()
    sqlite_api.sync_insert(sqlite_api.Site_Info, site_df.to_dict('records'))
    print('init insert: {:.5f}s'.format(perf_counter()-start))
    '''
    start = perf_counter()
    # entries = await SIFTS_Info.objects.all()
    # example = await SIFTS_Info.objects.get(UniProt='Q92793')
    # example = await Entry_Info.objects.filter(pdb_id__in=('4u7t', '6g6j')).all()
    # example = await SEQRES_Info.objects.filter(pdb_id__in=('4u7t', '6g6j')).all()
    # example = await SIFTS_Info.objects.filter(UniProt__in=('P12270', 'Q14980')).all()
    # example = await SIFTS_Info.objects.filter(UniProt='P12270').all()
    # example = await PDBRes_Info.objects.filter(pdb_id='4loe', entity_id="1", chain_id='C', obs_ratio__lt=1).all()
    # filter(pdb_id__in=['3wwq', '5v8w'])
    example = await sqlite_api.Entry_Info.objects.limit(1000).all()
    # example = await sqlite_api.PDBRes_Info.objects.filter(obs_ratio__lt=1).distinct()
    # example = await sqlite_api.Site_Info.objects.filter(from_id__in=('ENST00000379410', 'ENST00000379409')).all()
    # .filter(obs_ratio__lt=1).
    print('init select: {:.5f}s'.format(perf_counter()-start))
    res = pd.DataFrame(example)
    print(res)
    print(len(res), len(res.drop_duplicates()))
    # sqlite_api.sync_insert(sqlite_api.PDBRes_Info, [])


if __name__ == '__main__':
    main().result()

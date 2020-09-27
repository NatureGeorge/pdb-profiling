# @Created Date: 2020-08-12 11:27:46 pm
# @Filename: __init__.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-09-26 04:27:17 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pdb_profiling.processors.pdbe.record import (
    Base, 
    PDB, 
    PDBAssemble,
    PDBInterface,
    SIFTS,
    Compounds
    )
from pdb_profiling.processors.pdbe.api import PDBeModelServer, PDBArchive,PDBVersioned
from pdb_profiling.processors.uniprot.api import UniProtFASTA
from pdb_profiling.processors.proteins.api import ProteinsAPI
from pdb_profiling.processors.ensembl.api import EnsemblAPI
from pdb_profiling.processors.eutils.api import EutilsAPI
from pdb_profiling.processors.swissmodel.api import SMR
from pdb_profiling.processors.proteins import ProteinsDB
from pdb_profiling.log import Abclog
from pdb_profiling.utils import init_folder_from_suffix,a_seq_reader, a_load_json
from re import compile as re_compile
from pathlib import Path
from typing import Union, Optional
from unsync import unsync


class Identifier(Abclog):
    suffix = r'[0-9]+)[\.]*([0-9]*)'
    pats = {
        ('RefSeq', 'transcript'): re_compile(f'(NM_{suffix}'),
        ('RefSeq', 'protein'): re_compile(f'(NP_{suffix}'),
        ('Ensembl', 'gene'): re_compile(f'(ENSG{suffix}'),
        ('Ensembl', 'transcript'): re_compile(f'(ENST{suffix}'),
        ('Ensembl', 'protein'): re_compile(f'(ENSP{suffix}')}

    @classmethod
    def set_folder(cls, folder: Union[Path, str]):
        cls.folder = Path(folder)
        cls.sqlite_api = ProteinsDB(
            "sqlite:///%s" % (init_folder_from_suffix(cls.folder, 'local_db')/"proteinsAPI.db"))
        cls.proteins_api_folder = init_folder_from_suffix(
            cls.folder, 'proteins/api/proteins')
        cls.seq_folder = dict(
            RefSeq=init_folder_from_suffix(cls.folder, 'eutils/efetch'),
            Ensembl=init_folder_from_suffix(cls.folder, 'ensembl/sequence/id'))
        cls.ensembl_archive_folder = init_folder_from_suffix(
            cls.folder, 'ensembl/archive/id')

    @classmethod
    def get_type(cls, identifier: str):
        for key, pat in cls.pats.items():
            res = pat.fullmatch(identifier)
            if bool(res):
                return key, res.groups()

    def __init__(self, identifier: str, folder: Optional[Union[Path, str]] = None):
        try:
            (self.source, self.level), (self.identifier,
                                        self.version) = self.get_type(identifier)
            self.raw_identifier = identifier
            if folder is not None:
                self.set_db(folder)
            getattr(self, 'sqlite_api')
        except TypeError:
            raise ValueError(f"Unexcepted identifier type: {identifier}")
        except AttributeError:
            raise AttributeError(
                "Please specify class variable `sqlite_api` via set_folder() first or pass `folder` in this method!")
        self.ensembl_status = None
        self.refseq_status = None

    def __repr__(self):
        return f'<{self.source} {self.level} {self.identifier} {self.version}>'

    @unsync
    async def set_status(self):
        if self.source == 'Ensembl':
            headers = {'Content-Type': 'application/json'}
            self.ensembl_status = await EnsemblAPI.single_retrieve('archive/id/',
                                                                   self.identifier,
                                                                   headers,
                                                                   self.ensembl_archive_folder,
                                                                   Base.get_web_semaphore(),
                                                                   headers=headers).then(a_load_json)

    @unsync
    async def fetch_from_ProteinsAPI(self):
        dbReferences_df, iso_df = ProteinsAPI.pipe_summary(await ProteinsAPI.single_retrieve(
            'proteins/',
            dict(offset=0, size=-1, reviewed='true', isoform=0),
            self.proteins_api_folder,
            Base.get_web_semaphore(),
            identifier=f'{self.source}:{self.identifier}'
        ).then(a_load_json))
        if (dbReferences_df is not None) and (len(dbReferences_df) > 0):
            await self.sqlite_api.async_insert(
                self.sqlite_api.DB_REFERENCES, dbReferences_df.to_dict('records'))
            if iso_df is not None:
                await self.sqlite_api.async_insert(
                    self.sqlite_api.ALTERNATIVE_PRODUCTS, iso_df.rename(columns={'ids': 'isoform'}).to_dict('records'))
            else:
                self.logger.warning(
                    f"Can't find ALTERNATIVE_PRODUCTS with {self.identifier}")
        else:
            self.logger.warning(
                f"Can't find dbReference with {self.identifier}")

    @unsync
    async def map2unp(self):
        try:
            entry, isoform = await self.sqlite_api.database.fetch_one(
                query=f"""
                    SELECT Entry,isoform FROM dbReferences
                    WHERE type == '{self.source}' AND {self.level} == '{self.raw_identifier}'""")
        except TypeError:
            return

        if isoform is not None:
            sequenceStatus = await self.sqlite_api.database.fetch_val(
                query=f"""
                    SELECT sequenceStatus FROM ALTERNATIVE_PRODUCTS 
                    WHERE Entry == '{entry}' AND isoform == '{isoform}'""")
            if sequenceStatus == 'displayed':
                return entry, entry
        else:
            res = await self.sqlite_api.database.fetch_one(
                query=f"""
                    SELECT sequenceStatus FROM ALTERNATIVE_PRODUCTS 
                    WHERE Entry == '{entry}'""")
            if res is None:
                return entry, entry
        return entry, isoform

    @unsync
    async def fetch_sequence(self, newest: bool = True):
        if self.source == 'RefSeq':
            return await EutilsAPI.single_retrieve('efetch.fcgi',
                                                   dict(
                                                       db='sequences', id=self.identifier if newest else self.raw_identifier, rettype='fasta'),
                                                   self.seq_folder['RefSeq'],
                                                   Base.get_web_semaphore()).then(a_seq_reader)
        elif self.source == 'Ensembl':
            if self.ensembl_status is None:
                await self.set_status()
            if self.ensembl_status['is_current'] != '1':
                self.logger.warning(
                    f'Not exists in current archive: \n{self.ensembl_status}')
                return
            elif self.ensembl_status is False:
                self.logger.warning(
                    f'Invalid Identifier!')
            if not newest:
                self.logger.warning(
                    "Can't retrieve older version Ensembl Sequence via Ensembl REST API!")
            return await EnsemblAPI.single_retrieve('sequence/id/',
                                                    self.identifier,
                                                    dict(type='protein'),
                                                    self.seq_folder['Ensembl'],
                                                    Base.get_web_semaphore()).then(a_seq_reader)

# @Created Date: 2020-09-28 05:42:53 pm
# @Filename: record.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-09-28 05:43:34 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pdb_profiling.processors.ensembl.api import EnsemblAPI
from pdb_profiling.processors.eutils.api import EutilsAPI
from pdb_profiling.processors.proteins.api import ProteinsAPI
from pdb_profiling.processors.proteins import ProteinsDB
from pdb_profiling.utils import init_semaphore
from pdb_profiling.log import Abclog
from pdb_profiling.utils import init_folder_from_suffix,a_seq_reader, a_load_json
from re import compile as re_compile
from pathlib import Path
from typing import Union, Optional, Tuple
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
    @unsync
    async def set_web_semaphore(cls, *web_semaphore_values):
        if len(web_semaphore_values) == 0:
            proteins_api, ensembl_api, eutils_api = 20, 20, 20
        elif len(web_semaphore_values) < 3:
            num = web_semaphore_values[0]
            proteins_api, ensembl_api, eutils_api = num, num, num
        else:
            proteins_api, ensembl_api, eutils_api = web_semaphore_values[:3]
        cls.proteins_api_web_semaphore = await init_semaphore(proteins_api)
        cls.ensembl_api_web_semaphore = await init_semaphore(ensembl_api)
        cls.eutils_api_web_semaphore = await init_semaphore(eutils_api)

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
            raise ValueError(f"Unexpected identifier type: {identifier}")
        except AttributeError:
            raise AttributeError(
                "Please specify class variable `folder` via set_folder() first or pass `folder` in this method!")
        self.status = None

    def __repr__(self):
        return f'<{self.source} {self.level} {self.identifier} {self.version}>'

    @unsync
    async def set_status(self):
        if self.source == 'Ensembl':
            res = await EnsemblAPI.single_retrieve('archive/id/',
                self.identifier,
                None,
                self.ensembl_archive_folder,
                self.ensembl_api_web_semaphore,
                headers={'Content-Type': 'application/json'})
            if res is None:
                self.status = False
            else:
                self.status = await a_load_json(res)

    @unsync
    async def fetch_from_ProteinsAPI(self, reviewed='true'):
        dbReferences_df, iso_df = ProteinsAPI.pipe_summary(await ProteinsAPI.single_retrieve(
            'proteins/',
            dict(offset=0, size=-1, reviewed=reviewed, isoform=0),
            self.proteins_api_folder,
            self.proteins_api_web_semaphore,
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
                f"Can't find (reviewed) dbReference with {self.identifier}")

    @unsync
    async def get_all_level_identifiers(self):
        try:
            return dict(zip(('protein', 'transcript', 'gene'), await self.sqlite_api.database.fetch_one(
                query=f"""
                    SELECT protein,transcript,gene FROM dbReferences
                    WHERE type == '{self.source}' AND {self.level} == '{self.raw_identifier}'""")))
        except TypeError:
            return

    @unsync
    async def map2unp_from_localDB(self):
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
            res = await EutilsAPI.single_retrieve('efetch.fcgi',
                dict(db='sequences', id=self.identifier if newest else self.raw_identifier, rettype='fasta'),
                self.seq_folder['RefSeq'],
                self.eutils_api_web_semaphore)
            if res is not None:
                self.status = True
                return await a_seq_reader(res)
            else:
                self.status = False
                self.logger.warning(f'Invalid Identifier!')
        
        elif self.source == 'Ensembl':
            if self.status is None:
                await self.set_status()
            if self.status is False:
                self.logger.warning(f'Invalid Identifier!')
                return
            elif self.status['is_current'] != '1':
                self.logger.warning(
                    f'Not exists in current archive: \n{self.status}')
                return
            if not newest:
                self.logger.warning(
                    "Can't retrieve older version Ensembl Sequence via Ensembl REST API!")
            return await EnsemblAPI.single_retrieve('sequence/id/',
                self.identifier,
                dict(type='protein'),
                self.seq_folder['Ensembl'],
                self.ensembl_api_web_semaphore).then(a_seq_reader)

    @unsync
    async def map2unp(self, **kwargs):
        res = await self.map2unp_from_localDB()
        if res is None:
            await self.fetch_from_ProteinsAPI(**kwargs)
            res = await self.map2unp_from_localDB()
        return res

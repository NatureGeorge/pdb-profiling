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
from pdb_profiling.log import Abclog
from pdb_profiling.warnings import PossibleObsoletedUniProtWarning
from pdb_profiling.utils import init_folder_from_suffix, a_seq_reader, a_load_json, init_semaphore, unsync_wrap, unsync_run
from re import compile as re_compile
from pathlib import Path
from typing import Union, Optional, Tuple, Iterable, Callable, List
from unsync import unsync
from asyncio import as_completed
from pandas import DataFrame, concat
from collections import OrderedDict
from warnings import warn


class Identifier(Abclog):
    suffix = r'[0-9]+)[\.]*([0-9]*)'
    pats = OrderedDict({
        ('RefSeq', 'model'): re_compile('(X[A-Z]{1}_%s' % suffix),
        ('RefSeq', 'transcript'): re_compile(f'(NM_{suffix}'),
        ('RefSeq', 'protein'): re_compile(f'(NP_{suffix}'),
        ('Ensembl', 'gene'): re_compile(f'(ENSG{suffix}'),
        ('Ensembl', 'transcript'): re_compile(f'(ENST{suffix}'),
        ('Ensembl', 'protein'): re_compile(f'(ENSP{suffix}'),
        ('UniProt', 'isoform'): re_compile(r'^((?=.*[A-Za-z])(?=.*\d)[A-Za-z\d]{6,})[\-]*([0-9]*)$')
        })

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
                self.set_folder(folder)
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

    @classmethod
    def save_ProteinsAPI_data_to_DB(cls, res, identifier=''):
        dbReferences_df, other_dbReferences_df, iso_df, features_df, int_df, info_data = res

        cls.sqlite_api.sync_insert(
            cls.sqlite_api.INFO, [info_data])
        if features_df is not None:
            cls.sqlite_api.sync_insert(
                cls.sqlite_api.FEATURES, features_df.to_dict('records'))

        if (dbReferences_df is not None) and (len(dbReferences_df) > 0):
            cls.sqlite_api.sync_insert(
                cls.sqlite_api.DB_REFERENCES, dbReferences_df.to_dict('records'))
            if iso_df is not None:
                cls.sqlite_api.sync_insert(
                    cls.sqlite_api.ALTERNATIVE_PRODUCTS, iso_df.to_dict('records'))
            else:
                cls.logger.info(
                    f"Can't find ALTERNATIVE_PRODUCTS with {identifier}")
        else:
            cls.logger.info(
                f"Can't find (reviewed) dbReference with {identifier}")

        if len(other_dbReferences_df) > 0:
            cls.sqlite_api.sync_insert(
                cls.sqlite_api.OTHER_DB_REFERENCES, other_dbReferences_df.to_dict('records'))
        if int_df is not None and len(int_df) > 0:
            cls.sqlite_api.sync_insert(
                cls.sqlite_api.INTERACTION, int_df.to_dict('records'))

    @unsync
    async def fetch_from_ProteinsAPI_with_unp(self):
        assert self.source == 'UniProt'
        res = await ProteinsAPI.pipe_summary([await ProteinsAPI.single_retrieve(
            'proteins/',
            dict(),
            self.proteins_api_folder,
            self.proteins_api_web_semaphore,
            identifier=self.identifier
        ).then(a_load_json)])
        if res is None:
            return
        else:
            self.save_ProteinsAPI_data_to_DB(res, identifier=self.identifier)
            return res

    @unsync
    async def query_from_DB_with_unp(self, table_name:str, columns: str = '*', exists:Optional[bool]=None):
        default_tables = ('DB_REFERENCES', 'OTHER_DB_REFERENCES', 'ALTERNATIVE_PRODUCTS', 'FEATURES', 'INTERACTION', 'INFO')
        assert table_name in default_tables
        exists = (await self.sqlite_api.INFO.objects.filter(accession=self.identifier).exists()) if exists is None else exists
        if not exists:
            res = await self.fetch_from_ProteinsAPI_with_unp()
            if res is None:
                warn(repr(self), PossibleObsoletedUniProtWarning)
                return
            res = dict(zip(default_tables, res))
            return res[table_name]
        else:
            if columns == '*':
                return DataFrame(
                    await getattr(self.sqlite_api, table_name).objects.filter(accession=self.identifier).all())
            else:
                return DataFrame(
                    (await self.sqlite_api.database.fetch_all(query=f"SELECT {columns} FROM {table_name} WHERE accession = '{self.identifier}'")),
                    columns=columns.split(','))
            '''
            dbReferences_df = DataFrame(await cls.sqlite_api.DB_REFERENCES.objects.filter(accession=accession).all())
            other_dbReferences_df = DataFrame(await cls.sqlite_api.OTHER_DB_REFERENCES.objects.filter(accession=accession).all())
            iso_df = DataFrame(await cls.sqlite_api.ALTERNATIVE_PRODUCTS.objects.filter(accession=accession).all())
            features_df = DataFrame(res)
            int_df = DataFrame(await cls.sqlite_api.INTERACTION.objects.filter(accession1__contains=accession).all())
            return dbReferences_df, other_dbReferences_df, iso_df, features_df, int_df
            '''

    @unsync
    async def fetch_from_ProteinsAPI(self, reviewed='true', isoform=0):
        res = await ProteinsAPI.pipe_summary(await ProteinsAPI.single_retrieve(
            'proteins/',
            dict(offset=0, size=-1, reviewed=reviewed, isoform=isoform),
            self.proteins_api_folder,
            self.proteins_api_web_semaphore,
            identifier=f'{self.source}:{self.identifier}'
        ).then(a_load_json))
        if res is None:
            return
        else:
            self.save_ProteinsAPI_data_to_DB(res, identifier=self.identifier)
            return res

    @unsync
    async def get_all_level_identifiers(self):
        try:
            cur_id = self.raw_identifier if self.source == 'RefSeq' else self.identifier
            return dict(zip(('protein', 'transcript', 'gene'), await self.sqlite_api.database.fetch_one(
                query=f"""
                    SELECT protein,transcript,gene FROM dbReferences
                    WHERE type == '{self.source}' AND ({self.level} == '{cur_id}' OR {self.level} LIKE '{cur_id}%')""")))
        except TypeError:
            return

    @unsync
    async def map2unp_from_DB(self):
        '''
        return accession, UniProt Isoform, is_canonical
        '''
        cur_id = self.raw_identifier if self.source == 'RefSeq' else self.identifier
        res = await self.sqlite_api.database.fetch_one(
            query=f"""
                SELECT accession,isoform FROM dbReferences
                WHERE type == '{self.source}' AND ({self.level} == '{cur_id}' OR {self.level} LIKE '{cur_id}%')""")
        if res is None:
            return
        else:
            accession, isoform = res
        if isoform is not None:
            sequenceStatus = await self.sqlite_api.database.fetch_val(
                query=f"""
                    SELECT sequenceStatus FROM ALTERNATIVE_PRODUCTS 
                    WHERE accession == '{accession}' AND isoform == '{isoform}'""")
            assert sequenceStatus is not None, accession
            return self.raw_identifier, accession, isoform, sequenceStatus == 'displayed'
        else:
            sequenceStatus = await self.sqlite_api.database.fetch_one(
                query=f"""
                    SELECT sequenceStatus FROM ALTERNATIVE_PRODUCTS 
                    WHERE accession == '{accession}'""")
            if sequenceStatus is None:
                return self.raw_identifier, accession, accession, True
            else:
                return self.raw_identifier, accession, 'NaN', True

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
        if self.level == 'model':
            return self.raw_identifier, 'NaN', 'NaN', False
        elif self.source == 'UniProt':
            return self.raw_identifier, self.identifier, self.raw_identifier, (self.raw_identifier == self.identifier)
        try:
            res = await self.map2unp_from_DB()
        except AssertionError:
            res = None
        if res is None:
            val = await self.fetch_from_ProteinsAPI(**kwargs)
            if val is None:
                return self.raw_identifier, 'NaN', 'NaN', False
            else:
                res = await self.map2unp_from_DB()
        return res


class Identifiers(tuple):
    def __new__(cls, iterable: Iterable):
        return super(Identifiers, cls).__new__(cls, (Identifier(i) if isinstance(i, str) else i for i in iterable))

    def __getitem__(self, slice):
        if isinstance(slice, str):
            for i in self:
                if i.identifier == slice:
                    return i
        res = tuple.__getitem__(self, slice)
        if isinstance(res, Iterable):
            return self.__class__(res)
        else:
            return res

    def fetch(self, func: Union[Callable[[Identifier], List], str], **kwargs):
        if isinstance(func, str):
            self.tasks = [getattr(id_ob, func)(**kwargs) for id_ob in self]
        else:
            self.tasks = [func(id_ob, **kwargs) for id_ob in self]
        return self
    
    @unsync
    async def run(self, tqdm=None):
        if tqdm is None:
            return [await fob for fob in as_completed(self.tasks)]
        else:
            return [await fob for fob in tqdm(as_completed(self.tasks), total=len(self.tasks))]

    def query_from_DB_with_unps(self, table_name: str, columns: str = '*'):
        default_tables = ('DB_REFERENCES', 'OTHER_DB_REFERENCES', 'ALTERNATIVE_PRODUCTS', 'FEATURES', 'INTERACTION', 'INFO')
        assert table_name in default_tables
        obs = tuple(i for i in self if i.source == 'UniProt')
        if len(obs) == 0:
            self.tasks = []
            return self
        accessions = tuple(i.identifier for i in obs)
        if columns != '*' and table_name == 'INFO':
            task = Identifier.sqlite_api.database.fetch_all(query=f'SELECT {columns} FROM INFO WHERE accession IN {accessions}')
        else:
            task = Identifier.sqlite_api.INFO.objects.filter(accession__in=accessions).all()
        exists = unsync_run(task)
        if len(exists) == 0:
            self.tasks = [ob.query_from_DB_with_unp(table_name=table_name, columns=columns, exists=False) for ob in obs]
            return self
        else:
            exist_ids = frozenset(i.accession for i in exists)
            rest_ids = frozenset(accessions) - exist_ids
            rest_dfs = [self[accession].query_from_DB_with_unp(table_name=table_name, columns=columns, exists=False) for accession in rest_ids]
            if table_name == 'INFO':
                ap = unsync_wrap(exists)
            else:
                if columns == '*':
                    ap = unsync_wrap(getattr(Identifier.sqlite_api, table_name).objects.filter(accession__in=exist_ids).all())
                else:
                    ap = unsync_wrap(Identifier.sqlite_api.database.fetch_all(query=f'SELECT {columns} FROM {table_name} WHERE accession IN {tuple(exist_ids)}'))
            rest_dfs.append(ap)
            self.tasks = rest_dfs
            return self

# @Created Date: 2020-09-28 05:42:53 pm
# @Filename: record.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2022-07-27 05:44:55 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pdb_profiling.processors.recordbase import IdentifierBase
from pdb_profiling.processors.ensembl.api import EnsemblAPI
from pdb_profiling.processors.eutils.api import EutilsAPI
from pdb_profiling.processors.proteins.api import ProteinsAPI
from pdb_profiling.processors.proteins import ProteinsDB
from pdb_profiling.log import Abclog
from pdb_profiling.warnings import PossibleObsoletedUniProtWarning, SequenceConflictWarning
from pdb_profiling.utils import init_folder_from_suffix, init_folder_from_suffixes, a_seq_reader, a_load_json, init_semaphore, unsync_wrap, unsync_run, flatten_dict
from pathlib import Path
from typing import Union, Optional, Iterable, Callable, List
from unsync import unsync
from asyncio import as_completed
from pandas import DataFrame
from warnings import warn
from orm.models import NoMatch


class Identifier(Abclog, IdentifierBase):
    '''
    Impl EBI Proteins API
    '''
    auto_assign_when_seq_conflict = False

    @classmethod
    @unsync
    async def set_web_semaphore(cls, proteins_api=20, ensembl_api=20, eutils_api=20):
        cls.proteins_api_web_semaphore = await init_semaphore(proteins_api)
        cls.ensembl_api_web_semaphore = await init_semaphore(ensembl_api)
        cls.eutils_api_web_semaphore = await init_semaphore(eutils_api)

    @classmethod
    def set_folder(cls, folder: Union[Path, str]):
        cls.folder = Path(folder)
        cls.sqlite_api = ProteinsDB(
            "sqlite:///%s" % (init_folder_from_suffix(cls.folder, 'local_db')/"proteinsAPI.db"))
        cls.proteins_api_folder = cls.folder/'proteins/api/'
        tuple(init_folder_from_suffixes(
            cls.proteins_api_folder, ((i if i[-1] == '/' else i+'_') for i in ProteinsAPI.api_set)))
        cls.seq_folder = dict(
            RefSeq=init_folder_from_suffix(cls.folder, 'eutils/efetch'),
            Ensembl=init_folder_from_suffix(cls.folder, 'ensembl/sequence/id'))
        cls.ensembl_archive_folder = init_folder_from_suffix(
            cls.folder, 'ensembl/archive/id')

    def __post_init__(self):
        super().__post_init__()
        if not hasattr(self, 'sqlite_api'):
            raise AttributeError("Please specify class variable `folder` via set_folder() first or pass `folder` in this method!")
        self.status = None

    def __repr__(self):
        return f'<{self.source} {self.level} {self.identifier} {self.identifier_suffix}>'

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
        if features_df is not None and len(features_df) > 0:
            cls.sqlite_api.sync_insert(
                cls.sqlite_api.FEATURES, features_df.to_dict('records'))
        if dbReferences_df is not None and len(dbReferences_df) > 0:
            cls.sqlite_api.sync_insert(
                cls.sqlite_api.DB_REFERENCES, dbReferences_df.to_dict('records'))
        if iso_df is not None and len(iso_df) > 0:
            cls.sqlite_api.sync_insert(
                cls.sqlite_api.ALTERNATIVE_PRODUCTS, iso_df.to_dict('records'))
        if other_dbReferences_df is not None and len(other_dbReferences_df) > 0:
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
            ret = res[table_name]
            if ret is None:
                return DataFrame([])
            else:
                return ret
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
    async def fetch_from_proteins_api(self, api_suffix, id_suffix='', with_source:bool=False, params={}, rate=1.5):
        return await ProteinsAPI.single_retrieve(
            suffix=api_suffix, 
            params=params, 
            folder=self.proteins_api_folder,
            semaphore=self.proteins_api_web_semaphore,
            identifier=(f"{self.source}:" if with_source else '')+(self.raw_identifier if self.source in ('UniProt', 'Taxonomy') else self.identifier)+id_suffix,
            rate=rate)

    @classmethod
    def yield_mapping(cls, data):
        if isinstance(data, list):
            for sub in data:
                yield from cls.yield_mapping_unit(sub)
        elif isinstance(data, dict):
            yield from cls.yield_mapping_unit(data)

    @staticmethod
    def yield_mapping_unit(data):
        for gnCoordinate in data['gnCoordinate']:
            info = dict(accession=data['accession'],
                        ensemblGeneId=gnCoordinate['ensemblGeneId'], 
                        ensemblTranscriptId=gnCoordinate['ensemblTranscriptId'], 
                        ensemblTranslationId=gnCoordinate['ensemblTranslationId'],
                        chromosome=gnCoordinate['genomicLocation']['chromosome'],
                        #start=gnCoordinate['genomicLocation']['start'],
                        #end=gnCoordinate['genomicLocation']['end'],
                        reverseStrand=gnCoordinate['genomicLocation']['reverseStrand'])
            for record in gnCoordinate['genomicLocation']['exon']:
                to_flat = record.copy()
                flatten_dict(to_flat, 'proteinLocation')
                flatten_dict(to_flat, 'proteinLocation.begin')
                flatten_dict(to_flat, 'proteinLocation.end')
                flatten_dict(to_flat, 'genomeLocation')
                flatten_dict(to_flat, 'genomeLocation.begin')
                flatten_dict(to_flat, 'genomeLocation.end')
                to_flat.update(info)
                yield to_flat

    @unsync
    async def alignment_df(self, api_suffix='coordinates/', **kwargs):
        #assert self.source in ('Taxonomy', 'UniProt')
        return DataFrame(self.yield_mapping(
            await self.fetch_from_proteins_api(api_suffix, **kwargs).then(a_load_json))).rename(columns={'id': 'ensemblExonId'})

    @unsync
    async def fetch_proteins_from_ProteinsAPI(self, reviewed='true', **kwargs):
        res = await ProteinsAPI.pipe_summary(await self.fetch_from_proteins_api(
            'proteins/',
            with_source=True,
            params=dict(offset=0, size=-1, reviewed=reviewed),  #isoform=1
            **kwargs
        ).then(a_load_json))
        if res is None:
            """
            res = await ProteinsAPI.pipe_summary(await self.fetch_from_proteins_api(
                'proteins/',
                with_source=True,
                params=dict(offset=0, size=-1, reviewed=reviewed, isoform=0),
                **kwargs
            ).then(a_load_json))
            if res is None:
            """
            return
        else:
            self.save_ProteinsAPI_data_to_DB(res, identifier=self.identifier)
            return res

    @unsync
    async def get_isoform_ob(self):
        assert self.level == 'isoform'
        try:
            if self.identifier_suffix == '':
                query_ob = await self.get_canonical_isoform_ob()
                if query_ob is None:
                    return
                else:
                    query_id = query_ob.isoform
            else:
                query_id = self.raw_identifier
            return await self.sqlite_api.ALTERNATIVE_PRODUCTS.objects.get(isoform=query_id, sequenceStatus__in=('displayed', 'described'))
        except NoMatch:
            pass

    @unsync
    async def get_canonical_isoform_ob(self):
        assert self.level == 'isoform'
        try:
            return await self.sqlite_api.ALTERNATIVE_PRODUCTS.objects.get(accession=self.identifier, sequenceStatus='displayed')
        except NoMatch:
            pass

    @unsync
    async def get_all_ref_identifiers(self, to_dataframe:bool=True, **kwargs):
        if self.level == 'isoform':
            if self.identifier_suffix == '':
                c_ob = await self.get_canonical_isoform_ob()
                if c_ob is None:
                    query_args = {'accession': self.identifier}
                else:
                    query_args = {'isoform': c_ob.isoform}
            else:
                query_args = {'isoform': self.raw_identifier}
        else:
            query_args = {f"{self.level}__contains": self.identifier, 'type': self.source}
        query_args.update(kwargs)
        ret = await self.sqlite_api.DB_REFERENCES.objects.filter(**query_args).all()
        if len(ret) == 0:
            return None
        if to_dataframe:
            return DataFrame(ret)
        else:
            return ret

    @unsync
    async def get_all_level_identifiers(self):
        try:
            #cur_id = self.raw_identifier if self.source == 'RefSeq' else self.identifier
            cur_id = self.identifier
            return dict(zip(('protein', 'transcript', 'gene'), await self.sqlite_api.database.fetch_one(
                query=f"""
                    SELECT protein,transcript,gene FROM dbReferences
                    WHERE type == '{self.source}' AND ({self.level} == '{cur_id}' 
                    OR substr({self.level}, 0, instr({self.level}, '.')) == '{cur_id}'
                    )""")))
        except TypeError:
            return

    @unsync
    async def map2unp_from_DB(self):
        '''
        return accession, UniProt Isoform, is_canonical
        '''
        # cur_id = self.raw_identifier if self.source == 'RefSeq' else self.identifier
        cur_id = self.identifier
        res = await self.sqlite_api.database.fetch_one(
            query=f"""
                SELECT accession,isoform FROM dbReferences
                WHERE type == '{self.source}' AND ({self.level} == '{cur_id}' 
                OR substr({self.level}, 0, instr({self.level}, '.')) == '{cur_id}'
                )""")
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
                if self.auto_assign_when_seq_conflict:
                    if self.level == 'protein':
                        tseq = (await self.fetch_sequence())[1]
                    else:
                        tseq = (await Identifier((await self.get_all_level_identifiers())['protein']).fetch_sequence())[1]
                    tseqlen = len(tseq)
                    tseq_head = tseq[:50]
                    for isoform, cseqlen, cseq_head, sequenceStatus in (await self.sqlite_api.database.fetch_all(query=f"SELECT isoform,length(sequence),substr(sequence,1,50),sequenceStatus FROM ALTERNATIVE_PRODUCTS WHERE accession == '{accession}' AND sequenceStatus IN ('displayed', 'described') ORDER BY isoform")):
                        if cseqlen == tseqlen:
                            if sum(1 for x1, x2 in zip(cseq_head, tseq_head) if x1 == x2) >= 40:
                                warn(f'Exists sequence conflict between {self.raw_identifier} and {accession}(with isoforms), but still assign the {isoform} since their length of protein-seq are equal. This is a naive assignment and maybe error-prone!', SequenceConflictWarning)
                                return self.raw_identifier, accession, isoform, (sequenceStatus == 'displayed')
                    warn(f'Exists sequence conflict (also unequal length) between {self.raw_identifier} and {accession}(with isoforms), assign NaN instead.', SequenceConflictWarning)
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
                return (self.status.get('id', None), self.status.get('peptide', None))
            if not newest:
                self.logger.warning(
                    "Can't retrieve older version Ensembl Sequence via Ensembl REST API!")
            return await EnsemblAPI.single_retrieve('sequence/id/',
                self.identifier,
                dict(type='protein'),
                self.seq_folder['Ensembl'],
                self.ensembl_api_web_semaphore).then(a_seq_reader)

    @unsync
    async def unp_is_canonical(self):
        res = await self.query_from_DB_with_unp('ALTERNATIVE_PRODUCTS', columns='isoform,sequenceStatus')
        if res is None:
            return
        if self.identifier_suffix == '':
            return True
        if res.shape[0] == 0:
            if self.identifier_suffix != '1':
                self.logger.warning(f'Possbile invalid isoform identifier {self.raw_identifier}')
            return True
        else:
            focus = res[res.isoform.eq(self.raw_identifier) & res.sequenceStatus.isin(('described','displayed'))]
            if focus.shape[0] == 0:
                self.logger.error(f'Invalid isoform identifier {self.raw_identifier}')
                return False
            else:
                assert focus.shape[0] == 1
                return focus.iloc[0]['sequenceStatus'] == 'displayed'

    @unsync
    async def init(self):
        if hasattr(self, 'inited'):
            return self
        else:
            self.inited = True
            await self.map2unp()
            return self

    @unsync
    async def map2unp(self, **kwargs):
        if self.source == 'UniProt':
            is_canonical = await self.unp_is_canonical()
            if is_canonical is None:
                return self.raw_identifier, 'NaN', 'NaN', False
            return self.raw_identifier, self.identifier, (self.identifier if is_canonical else self.raw_identifier), is_canonical
        elif self.source not in ('RefSeq', 'Ensembl'):
            res = None
        try:
            res = await self.map2unp_from_DB()
        except AssertionError:
            res = None
        if res is None:
            val = await self.fetch_proteins_from_ProteinsAPI(**kwargs)
            if val is None:
                return self.raw_identifier, 'NaN', 'NaN', False
            else:
                res = await self.map2unp_from_DB()
        if res is None:
            warn(f"Unexpected None: {self.raw_identifier}")
            return self.raw_identifier, 'NaN', 'NaN', False
        else:
            return res


class Identifiers(tuple):
    '''immutable iterable class (tuple-like)'''

    def __new__(cls, iterable: Iterable=tuple()):
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

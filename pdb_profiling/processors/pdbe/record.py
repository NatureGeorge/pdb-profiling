# @Created Date: 2020-08-11 10:48:08 pm
# @Filename: record.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2023-03-25 08:08:25 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Iterable, Iterator, Union, Callable, Optional, Hashable, Dict, Coroutine, List, Tuple
from inspect import isawaitable
from functools import partial, reduce
from numpy import array, where as np_where, count_nonzero, nan, exp, square
from pathlib import Path
from pandas import isna, concat, DataFrame, Series, merge, read_csv
from unsync import unsync, Unfuture
from asyncio import as_completed
from aiofiles import open as aiofiles_open
from aiofiles.os import remove as aiofiles_os_remove
from re import compile as re_compile
import orjson as json
from io import StringIO
import zlib
from bisect import bisect_left
from collections import defaultdict, namedtuple, OrderedDict
from itertools import product, combinations_with_replacement, combinations
from operator import itemgetter
from pdb_profiling.processors.recordbase import IdentifierBase
from pdb_profiling.processors.transformer import Dict2Tabular
from pdb_profiling.exceptions import *
from pdb_profiling.cython.cyrange import to_interval, lyst22interval, lyst32interval, range_len, subtract_range, add_range, overlap_range, outside_range, trim_range, overlap_range_2, overlap_range_3
from pdb_profiling.utils import (init_semaphore, init_folder_from_suffix, 
                                 a_read_csv, split_df_by_chain, unsync_wrap,
                                 related_dataframe, slice_series, 
                                 MMCIF2DictPlus, a_load_json, SeqRangeReader,
                                 sort_2_range, flat_dict_in_df,
                                 get_diff_index, get_seq_seg, id2score,
                                 get_gap_list, get_range_diff,
                                 select_range, expand_interval,
                                 lyst2range, select_ho_max_range,
                                 select_he_range, init_folder_from_suffixes,
                                 a_seq_reader, dumpsParams, get_str_dict_len, aa_three2one)
from pdb_profiling.processors.pdbe.api import ProcessPDBe, PDBeModelServer, PDBeCoordinateServer, PDBArchive, PDBeKBAnnotations, FUNCS as API_SET, FetchBase
from pdb_profiling.processors.uniprot.api import UniProtINFO
from pdb_profiling.processors.pdbe import PDBeDB
from pdb_profiling.processors.rcsb import RCSBDB
from pdb_profiling.processors.rcsb.api import RCSBDataAPI, RCSBSearchAPI, RCSB1DCoordinatesAPI
from pdb_profiling.processors.swissmodel.api import SMR
from pdb_profiling.data import miyata_similarity_matrix
from pdb_profiling import cif_gz_stream
from pdb_profiling.processors.i3d.api import Interactome3D
from pdb_profiling.warnings import (WithoutCifKeyWarning, PISAErrorWarning, 
                                    ConflictChainIDWarning, PossibleObsoletedUniProtWarning,
                                    PossibleObsoletedPDBEntryWarning, SkipAssemblyWarning,
                                    PeptideLinkingWarning, MultiWrittenWarning, WithoutRCSBClusterMembershipWarning,
                                    PDBeKBResidueMappingErrorWarning)
from pdb_profiling.ensure import aio_file_exists_stat
from pdb_profiling.cython.py_qcprot import rmsd
from textdistance import sorensen
from warnings import warn
from cachetools import LRUCache
from sklearn.neighbors import NearestNeighbors
from random import choice
from tenacity import retry, wait_random, stop_after_attempt, RetryError
from parasail import nw_trace_scan_sat, blosum62


API_SET = {api for apiset in API_SET for api in apiset[1]}


class PropertyRegister(object):
    '''
    ClassDecorator to Register Property Method

        * take method/function name
        * specified pipeline func: `get_property_suffixes`, `__call__`
    '''

    properties = []

    @classmethod
    def get_property_suffixes(cls):
        return {key: f'api/pdb/entry/{key}/' for key in cls.properties}

    def __init__(self, func):
        self._name = func.__name__
        PropertyRegister.properties.append(self._name)
    
    def __call__(self, that):
        that.init_properties(self.get_property_suffixes())
        return getattr(that, f'_{self._name}')


class Base(IdentifierBase):
    '''
    Impl
        * PDBe Entry-Based API
        * PDBe Graph API
        * RCSB Data API
        * RCSB Search API
        * RCSB 1D-Coordinates API
        * ...
    '''
    
    folder = None
    tasks = LRUCache(maxsize=1024)

    def get_id(self):
        pass

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.get_id()}>"

    @classmethod
    def register_task(cls, key: Hashable, task: Unfuture):
        cls.tasks[key] = task

    @classmethod
    def pipe_register_task(cls, key, task_func, **kwargs):
        task = cls.tasks.get(key, None)
        if task is None:
            task = task_func(**kwargs)
            cls.tasks[key] = task
        return task

    @classmethod
    @unsync
    async def set_web_semaphore(cls, web_semaphore_value: int):
        assert web_semaphore_value > 1
        cls.web_semaphore = await init_semaphore(web_semaphore_value)
    
    @classmethod
    @unsync
    async def set_rcsb_web_semaphore(cls, web_semaphore_value: int):
        assert web_semaphore_value > 1
        cls.rcsb_semaphore = await init_semaphore(web_semaphore_value)

    @classmethod
    def get_web_semaphore(cls):
        return cls.web_semaphore

    @classmethod
    def set_folder(cls, folder: Union[Path, str]):
        """Set your folder path

        Args:
            folder (Union[Path, str]): the path to set
        """        
        folder = Path(folder)
        assert folder.exists(), "Folder not exist! Please create it or input a valid folder!"
        cls.folder = folder
        tuple(init_folder_from_suffixes(cls.folder, API_SET))
        tuple(init_folder_from_suffixes(cls.folder/'data_rcsb', RCSBDataAPI.api_set | {'graphql', 'search', '1d_coordinates'}))
        cls.sqlite_api = PDBeDB("sqlite:///%s" % (init_folder_from_suffix(cls.folder, 'local_db')/"PDBeDB.db"))
        cls.rcsb_sqlite_api = RCSBDB("sqlite:///%s" % (init_folder_from_suffix(cls.folder, 'local_db')/"RCSBDB.db"))
        init_folder_from_suffix(cls.folder, 'for_version_control')

    @classmethod
    def get_folder(cls) -> Path:
        cls.check_folder()
        return cls.folder

    @classmethod
    def check_folder(cls):
        if cls.folder is None:
            raise ValueError(f"Please set folder via {cls.__name__}.set_folder(folder: Union[Path, str])")
    
    @classmethod
    def infer_ret_from_args(cls, api_suffix: str, identifier: str):
        return cls.get_folder()/api_suffix/'{}+{}.tsv'.format(
            api_suffix.replace('/', '%'), 
            identifier.replace('/', '%'))

    @classmethod
    def r_task(cls, task, then_func, api_suffix, identifier, json):
        if then_func is not None:
            task = task.then(then_func)
        cls.register_task((cls.__class__.__name__, api_suffix, then_func, json, identifier), task)
        return task

    def fetch_from_pdbe_api(self, api_suffix: str, then_func: Optional[Callable[[Unfuture], Unfuture]] = None, json: bool = False, mask_id: str = None, infer_path: bool = True, **kwargs) -> Unfuture:
        """fetch data from PDBe API

        Args:
            api_suffix (str): the suffix of the API that you want to retrieve info.
            then_func (Optional[Callable[[Unfuture], Unfuture]], optional): function arg that pass to Unfuture.then(). Defaults to None.
            json (bool, optional): whether the data is treated and returned as JSON. Defaults to False.
            mask_id (str, optional): Defaults to None.

        Returns:
            Unfuture: Unfuture object
        """        
        assert api_suffix in API_SET, f"Invlaid API SUFFIX! Valid set:\n{API_SET}"
        identifier = self.get_id() if mask_id is None else mask_id
        task = self.tasks.get((self.__class__.__name__, api_suffix, then_func, json, identifier), None)
        if task is not None:
            return task

        if infer_path and not json:
            infer = self.infer_ret_from_args(api_suffix, identifier)
            if infer.exists():
                task = unsync_wrap(infer)
                return self.r_task(task, then_func, api_suffix, identifier, json)

        args = dict(pdb=identifier,
                    suffix=api_suffix,
                    method='get',
                    folder=self.get_folder()/api_suffix,
                    semaphore=self.get_web_semaphore())
        args = {**kwargs, **args}
        if json:
            args['to_do_func'] = None
        task = ProcessPDBe.single_retrieve(**args)
        return self.r_task(task, then_func, api_suffix, identifier, json)
    
    def fetch_from_rcsb_api(self, api_suffix: str, query=None, then_func: Optional[Callable[[Unfuture], Unfuture]] = None, json: bool = False, mask_id: str = None, **kwargs):
        task = self.tasks.get((repr(self), api_suffix, query, then_func, json, mask_id), None)
        if task is not None:
            return task
        if api_suffix in RCSBDataAPI.api_set:
            args = dict(identifier=self.get_id() if mask_id is None else mask_id,
                        suffix=api_suffix,
                        folder=self.get_folder()/'data_rcsb'/api_suffix,
                        semaphore=self.rcsb_semaphore)
            task_func = RCSBDataAPI.single_retrieve
        elif api_suffix == 'graphql':
            args = dict(query=query, folder=self.get_folder()/'data_rcsb/graphql', semaphore=self.rcsb_semaphore)
            task_func = RCSBDataAPI.graphql_retrieve
        elif api_suffix == 'search':
            args = dict(query=query, folder=self.get_folder()/'data_rcsb/search', semaphore=self.rcsb_semaphore)
            task_func = RCSBSearchAPI.single_retrieve
        elif api_suffix == '1d_coordinates':
            args = dict(query=query, folder=self.get_folder()/'data_rcsb/1d_coordinates', semaphore=self.rcsb_semaphore)
            task_func = RCSB1DCoordinatesAPI.graphql_retrieve
        else:
            raise AssertionError(f"Invlaid API SUFFIX! Valid set:\n{RCSBDataAPI.api_set} or graphql or search")
        args = {**kwargs, **args}
        if json:
            args['to_do_func'] = None
        task = task_func(**args)
        if then_func is not None:
            task = task.then(then_func)
        self.register_task((repr(self), api_suffix, query, then_func, json, mask_id), task)
        return task

    @classmethod
    @unsync
    async def to_dataframe(cls, path):
        if isawaitable(path):
            path = await path
        if path is None:
            return None
        df = await a_read_csv(path, sep="\t", converters=ProcessPDBe.converters, keep_default_na=False, na_values=['NULL', 'null', ''])
        return df

    @classmethod
    @unsync
    async def to_dataframe_with_kwargs(cls, path, **kwargs):
        path = await path
        if path is None:
            return None
        default_na_values = frozenset({'NULL', 'null', ''})
        kwargs['na_values'] = list(frozenset(kwargs.get('na_values', default_na_values)) | default_na_values)
        if len(kwargs['na_values']) > 3:
            assert ' ' in kwargs['na_values']
        if 'converters' in kwargs.keys():
            converters = {**ProcessPDBe.converters, **kwargs['converters']}
            del kwargs['converters']
        else:
            converters = ProcessPDBe.converters
        df = await a_read_csv(path, sep="\t", converters=converters, keep_default_na=False, **kwargs)
        return df
    
    @staticmethod
    @unsync
    async def result_set_to_dataframe(data):
        if isawaitable(data):
            data = await data
        if data is None:
            return
        elif isinstance(data, Dict):
            pass
        else:
            data = await a_load_json(data)
        return DataFrame(data['result_set'])


class PDB(Base):

    tasks = LRUCache(maxsize=1024)

    protein_sequence_pat = re_compile(r'([A-Z]{1}|\([A-Z0-9]+\))')
    nucleotide_sequence_pat = re_compile(r'([AUCGI]{1}|\(DA\)|\(DT\)|\(DC\)|\(DG\)|\(DI\)|\(UNK\))')
    author_residue_pat = re_compile(r'([\-]*[0-9]+)([A-Z]*)')
    binding_details_pat = re_compile(r'BINDING SITE FOR RESIDUE (.+) ([A-Z]+) ([\-]*[0-9]+)([A-Z]*)')
    StatsProteinEntitySeq = namedtuple(
        'StatsProteinEntitySeq', 'pdb_id molecule_type entity_id ca_p_only SEQRES_COUNT STD_INDEX STD_COUNT NON_INDEX NON_COUNT UNK_INDEX UNK_COUNT ARTIFACT_INDEX')
    StatsNucleotideEntitySeq = namedtuple(
        'StatsNucleotideEntitySeq', 'pdb_id molecule_type entity_id ca_p_only SEQRES_COUNT dNTP_INDEX dNTP_COUNT NTP_INDEX NTP_COUNT NON_INDEX NON_COUNT UNK_INDEX UNK_COUNT')
    StatsChainSeq = namedtuple(
        'StatsChainSeq', 'pdb_id entity_id chain_id struct_asym_id OBS_INDEX OBS_COUNT OBS_RATIO_ARRAY BINDING_LIGAND_INDEX BINDING_LIGAND_COUNT')

    fix_molecule_type_dict = {
        'Bound': 'bound',
        'Carbohydrate-polymer': 'carbohydrate polymer'
    }

    @classmethod
    def set_folder(cls, folder: Union[Path, str]):
        super().set_folder(folder)
        cls.assembly_cif_folder = cls.folder/'pdbe_assembly_cif'
        cls.assembly_cif_folder.mkdir(parents=True, exist_ok=True)
        tuple(init_folder_from_suffixes(cls.folder/'model-server', PDBeModelServer.api_set))
        tuple(init_folder_from_suffixes(cls.folder/'pdb/data/structures', PDBArchive.api_set))
        tuple(init_folder_from_suffixes(cls.folder/'coordinate-server', PDBeCoordinateServer.api_set))
        tuple(init_folder_from_suffixes(cls.folder/'pdbe-kb/annotations', PDBeKBAnnotations.api_set))

    @unsync
    async def prepare_property(self, raw_data):
        raw_data = await raw_data
        if raw_data is None:
            return None
        data = raw_data[self.pdb_id]
        assert len(data) == 1, f"{repr(self)}: Unexpected status data length\n{data}"
        return data[0]
    
    @unsync
    async def prepare_properties(self, property_suffixes):
        tasks = {
            var_property: self.pdb_ob.fetch_from_pdbe_api(
                    suffix, a_load_json, json=True
                ).then(self.prepare_property) for var_property, suffix in property_suffixes.items()}
        for var_property, task in tasks.items():
            setattr(self, f'_{var_property}', await task)

    def init_properties(self, property_suffixes):
        if not self.properties_inited:
            self.prepare_properties(property_suffixes).result()
            self.properties_inited = True

    @property
    @PropertyRegister
    def summary(self):
        pass

    @property
    @PropertyRegister
    def status(self):
        """
        Check PDB Entry Status

            * status_code: REL -> relase
            * status_code: OBS -> obsolete

        >>> if self.status['status_code'] == 'REL':
                logging.info(self.status['obsoletes'])
            elif self.status['status_code'] == 'OBS':
                logging.warning(self.status['superceded_by'])
        """
        pass

    def __post_init__(self):
        super().__post_init__()
        self.check_folder()
        self.set_id()
        self.pdb_ob = self
        self.properties_inited = False
    
    def set_id(self):
        assert self.source == 'PDB'
        self.pdb_id = self.identifier.lower()

    def get_id(self):
        return self.pdb_id

    def fetch_from_coordinateServer_api(self, api_suffix: str, then_func: Optional[Callable[[Unfuture], Unfuture]] = None, root='litemol', kwargs=dict(), **params):
        assert api_suffix in PDBeCoordinateServer.api_set, f"Invlaid API SUFFIX! Valid set:\n{PDBeCoordinateServer.api_set}"
        dparams = dumpsParams(params)
        task = self.tasks.get((repr(self), 'PDBeCoordinateServer', root, api_suffix, dparams, then_func), None)
        if task is not None:
            return task
        task = PDBeCoordinateServer(root).single_retrieve(
            pdb_id=self.pdb_id,
            suffix=api_suffix,
            params=params,
            folder=self.get_folder()/'coordinate-server'/api_suffix,
            semaphore=self.get_web_semaphore(),
            **kwargs)
        if then_func is not None:
            task = task.then(then_func)
        self.register_task((repr(self), 'PDBeCoordinateServer', root, api_suffix, dparams, then_func), task)
        return task

    def fetch_from_modelServer_api(self, api_suffix: str, method: str = 'post', data_collection=None, then_func: Optional[Callable[[Unfuture], Unfuture]] = None, filename='_subset', kwargs=dict(), **params) -> Unfuture:
        assert api_suffix in PDBeModelServer.api_set, f"Invlaid API SUFFIX! Valid set:\n{PDBeModelServer.api_set}"
        dparams = dumpsParams(params) if len(params) > 0 else None
        task = self.tasks.get((repr(self), PDBeModelServer.root, api_suffix, method, data_collection, dparams, then_func, filename), None)
        if task is not None:
            return task
        task = PDBeModelServer.single_retrieve(
            pdb=self.pdb_id,
            suffix=api_suffix,
            method=method,
            folder=self.get_folder()/'model-server'/api_suffix,
            semaphore=self.get_web_semaphore(),
            params=params,
            data_collection=data_collection,
            filename=filename,
            **kwargs)
        if then_func is not None:
            task = task.then(then_func)
        self.register_task((repr(self), PDBeModelServer.root, api_suffix, method, data_collection, dparams, then_func), task)
        return task

    def fetch_from_PDBArchive(self, api_suffix: str, then_func: Optional[Callable[[Unfuture], Unfuture]] = None, **kwargs) -> Unfuture:
        assert api_suffix in PDBArchive.api_set, f"Invlaid API SUFFIX! Valid set:\n{PDBArchive.api_set}"
        task = self.tasks.get((repr(self), PDBArchive.root, api_suffix, kwargs.get('file_suffix', None), then_func), None)
        if task is not None:
            return task
        task = PDBArchive.single_retrieve(
            pdb=self.pdb_id,
            suffix=api_suffix,
            folder=self.get_folder()/'pdb/data/structures'/api_suffix,
            semaphore=self.get_web_semaphore(),
            **kwargs)
        if then_func is not None:
            task = task.then(then_func)
        self.register_task((repr(self), PDBArchive.root, api_suffix, kwargs.get('file_suffix', None), then_func), task)
        return task

    def fetch_from_pdbekb_annotations(self, api_suffix: str, then_func=None, **kwargs):
        assert api_suffix in PDBeKBAnnotations.api_set, f"Invlaid API SUFFIX! Valid set:\n{PDBArchive.api_set}"
        task = self.tasks.get((repr(self), PDBeKBAnnotations.root, api_suffix, then_func), None)
        if task is not None:
            return task
        task = PDBeKBAnnotations.single_retrieve(
            pdb=self.pdb_id,
            suffix=api_suffix,
            folder=self.get_folder()/'pdbe-kb/annotations'/api_suffix,
            semaphore=self.get_web_semaphore(),
            **kwargs)
        if then_func is not None:
            task = task.then(then_func)
        self.register_task((repr(self), PDBeKBAnnotations.root, api_suffix, then_func), task)
        return task

    @unsync
    async def pipe_pdbekb_annotations(self, api_suffix, merge_residue_listing:bool=True, only_protein_res: bool = False, **kwargs):
        data = await self.fetch_from_pdbekb_annotations(api_suffix, **kwargs).then(a_load_json)
        if data is None:
            return
        try:
            dfrm = Dict2Tabular.pandas_io(PDBeKBAnnotations.yieldPDBeKBAnnotations(data))
        except ValueError:
            warn(f"{repr(self)}: {api_suffix}: without expected data: {data}")
            return
        #dfrm.drop(columns=['site_data'], inplace=True)
        dfrm[['author_residue_number', 'author_insertion_code']] = dfrm.pdb_res_label.astype(str).apply(
            lambda x: self.author_residue_pat.fullmatch(x).groups()).apply(Series)
        dfrm.author_residue_number = dfrm.author_residue_number.astype(int)
        if merge_residue_listing:
            dfrm = dfrm.merge(await self.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe), how='left')
            assert not dfrm.aa_type.isnull().any()
            '''
            Exception Case:  (not 'X' but raw code)
            ---------------------------------------
            data_resource            M-CSA
            pdb_id                    1auk
            chain_id                     A
            struct_asym_id               A
            residue_name               FGL <-> aa_type                    Gly
            residue_number              51
            author_residue_number       69
            author_insertion_code
            '''
            # assert all((dfrm.aa_type.str.upper() == dfrm.residue_name) | (dfrm.aa_type.eq('X')))
            #, f"{repr(self)}: {api_suffix}:\n{dfrm[dfrm.aa_type.str.upper() != dfrm.residue_name].iloc[0].T}"
            dfrm.drop(columns=['aa_type', 'pdb_res_label'],inplace=True)
            if only_protein_res:
                mo_df = await self.fetch_from_pdbe_api('api/pdb/entry/molecules/', Base.to_dataframe)
                focus_entities = mo_df[mo_df.molecule_type.str.startswith('polypeptide')].entity_id
                return dfrm[dfrm.entity_id.isin(focus_entities)]
        return dfrm

    @staticmethod
    @unsync
    async def a_load_cif(path: Union[Unfuture, Coroutine, str, Path], *args) -> MMCIF2DictPlus:
        if isawaitable(path):
            path = await path
        if not isinstance(path, Path):
            path = Path(path)
        if path.suffix == '.gz':
            async with aiofiles_open(path, 'rb') as file_io:
                handle = StringIO(zlib.decompress((await file_io.read()), zlib.MAX_WBITS | 32).decode('utf-8'))
                return MMCIF2DictPlus(handle, *args)
        else:
            async with aiofiles_open(path, 'rt') as file_io:
                handle = await file_io.readlines()
                return MMCIF2DictPlus((iline for iline in handle), *args)

    @classmethod
    @unsync
    async def cif2atom_sites_df(cls, path: Union[Unfuture, Coroutine, Path]):
        cols = ('_atom_site.type_symbol', '_atom_site.label_atom_id',
                '_atom_site.label_comp_id', '_atom_site.label_seq_id',
                '_atom_site.label_asym_id', '_atom_site.Cartn_x',
                '_atom_site.Cartn_y', '_atom_site.Cartn_z')
        mmcif_dict = await cls.a_load_cif(path, cols+('_coordinate_server_result.has_error',))
        if '_coordinate_server_result.has_error' in mmcif_dict:
            assert mmcif_dict['_coordinate_server_result.has_error'][0] == 'no', str(mmcif_dict)
        try:
            dfrm = DataFrame(zip(*[mmcif_dict[col] for col in cols]), columns=cols)
        except KeyError:
            raise AssertionError(await path if isawaitable(path) else path)
        return dfrm

    @staticmethod
    def get_coordinates(x):
        return x[['_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z']].to_numpy(dtype=float)

    @unsync
    async def get_overlap_representative_atoms(self, struct_asym_id: str, label_atom_id: str = 'CA', model_nums: int = 1, overlap_pdb_segs: Optional[Iterable] = None):
        df = (await self.fetch_from_modelServer_api('atoms',
                                                    method='get',
                                                    filename=f'.{struct_asym_id}({label_atom_id})',
                                                    label_asym_id=struct_asym_id,
                                                    label_atom_id=label_atom_id,
                                                    model_nums=model_nums,
                                                    encoding='cif',
                                                    copy_all_categories=False
                                                    ).then(PDB.cif2atom_sites_df)
              ).drop_duplicates(subset='_atom_site.label_seq_id', keep='first')
        if overlap_range is not None:
            df = df.merge(DataFrame({'_atom_site.label_seq_id': list(expand_interval(overlap_pdb_segs))}, dtype=str), how='right')
        return df.reset_index(drop=True)

    @unsync
    async def get_representative_atoms(self, struct_asym_id: str, label_atom_id: str = 'CA', model_nums: int = 1, overlap_pdb_segs: Optional[Iterable] = None):
        df = (await self.fetch_from_modelServer_api('atoms',
                                                    method='get',
                                                    filename=f'.{struct_asym_id}({label_atom_id})',
                                                    label_asym_id=struct_asym_id,
                                                    label_atom_id=label_atom_id,
                                                    model_nums=model_nums,
                                                    encoding='cif',
                                                    copy_all_categories=False
                                                    ).then(PDB.cif2atom_sites_df)
              ).drop_duplicates(subset='_atom_site.label_seq_id', keep='first')
        if overlap_range is not None:
            return self.get_coordinates(df.merge(DataFrame({'_atom_site.label_seq_id': list(expand_interval(overlap_pdb_segs))}, dtype=str)))
        else:
            return self.get_coordinates(df)

    @classmethod
    def get_coordinates_dict(cls, df):
        return df.groupby(['_atom_site.label_asym_id', '_atom_site.label_seq_id']).apply(cls.get_coordinates).to_dict()

    @classmethod
    @unsync
    async def cif2residue_listing(cls, path: Union[Unfuture, str, Path]):
        cols = ('_pdbx_poly_seq_scheme.asym_id',
                '_pdbx_poly_seq_scheme.entity_id',
                '_pdbx_poly_seq_scheme.seq_id',
                '_pdbx_poly_seq_scheme.mon_id',
                #'_pdbx_poly_seq_scheme.ndb_seq_num', # NDB residue number
                '_pdbx_poly_seq_scheme.pdb_seq_num',
                '_pdbx_poly_seq_scheme.pdb_strand_id',
                '_pdbx_poly_seq_scheme.pdb_ins_code')
        new_cols = ('struct_asym_id',
                    'entity_id',
                    'residue_number',
                    'residue_name',
                    #'residue_number?',
                    'authore_residue_number',
                    'chain_id',
                    'author_insertion_code')
        mmcif_dict = await cls.a_load_cif(path, cols)
        mmcif_dict['data_'] = mmcif_dict['data_'].lower()
        dfrm = DataFrame(mmcif_dict)
        col_dict = dict(zip(cols, new_cols))
        col_dict['data_'] = 'pdb_id'
        dfrm.rename(columns=col_dict, inplace=True)
        #assert (dfrm['residue_number'] == dfrm['residue_number?']).all(), f"Unexpectd Cases: _pdbx_poly_seq_scheme.seq_id != _pdbx_poly_seq_scheme.ndb_seq_num\n{dfrm[dfrm['residue_number'] != dfrm['residue_number?']]}"
        #dfrm.drop(columns=['residue_number?'], inplace=True)
        for col in ('residue_number', 'authore_residue_number'):
            dfrm[col] = dfrm[col].astype(int)
        dfrm.author_insertion_code = dfrm.author_insertion_code.apply(lambda x: '' if x == '.' else x)
        return dfrm

    @classmethod
    @unsync
    async def assg_oper_cif_base(cls, path):
        assg_cols = ('_pdbx_struct_assembly_gen.asym_id_list',
                    '_pdbx_struct_assembly_gen.oper_expression',
                    '_pdbx_struct_assembly_gen.assembly_id')
        oper_cols = ('_pdbx_struct_oper_list.id',
                    '_pdbx_struct_oper_list.symmetry_operation',
                    '_pdbx_struct_oper_list.name')
        mmcif_dict = await cls.a_load_cif(path, assg_cols+oper_cols)
        if len(mmcif_dict) < 7:
            return None
        assg_df = DataFrame(list(zip(*[mmcif_dict[col] for col in assg_cols])), columns=assg_cols)
        assg_df = split_df_by_chain(assg_df, assg_cols, assg_cols[0:1]).rename(columns={col: col.split(
            '.')[1] for col in assg_cols}).rename(columns={'asym_id_list': 'struct_asym_id'}).reset_index(drop=True)
        oper_dict_1 = dict(zip(*[mmcif_dict[col] for col in oper_cols[:2]]))
        oper_dict_2 = dict(zip(*[mmcif_dict[col] for col in (oper_cols[0], oper_cols[2])]))
        return assg_df, oper_dict_1, oper_dict_2
    
    @staticmethod
    @unsync
    async def assg_oper_json_base(path):
        data = await a_load_json(path)
        oper_dict_1, oper_dict_2 = dict(), dict()
        pdbx_struct_assembly_gen_lyst = []
        for ass in data['data']['assemblies']:
            for i in ass['pdbx_struct_oper_list']:
                cur_id = i['id']
                cur_so = i['symmetry_operation']
                cur_so = '?' if cur_so is None else cur_so
                cur_na = i['name']
                cur_na = '?' if cur_na is None else cur_na
                assert oper_dict_1.get(cur_id, cur_so) == cur_so
                assert oper_dict_2.get(cur_id, cur_na) == cur_na
                oper_dict_1[cur_id] = cur_so
                oper_dict_2[cur_id] = cur_na
            pdbx_struct_assembly_gen_lyst.extend(ass['pdbx_struct_assembly_gen'])
        assg_df = DataFrame(sorted(pdbx_struct_assembly_gen_lyst, key=itemgetter('ordinal'))).drop(columns=['ordinal']).rename(columns={'asym_id_list': 'struct_asym_id'})
        assg_df = split_df_by_chain(assg_df, assg_df.columns, ['struct_asym_id'], mode='list').reset_index(drop=True)
        return assg_df, oper_dict_1, oper_dict_2

    @classmethod
    @unsync
    async def to_assg_oper_df(cls, path: Union[Unfuture, Coroutine, str, Path]):
        '''
        EXAMPLE OF model_id != asym_id_rank
        
            * 3hl2
            * 3ue7
        '''
        
        def to_rank(rank_dict, assembly_id, struct_asym_id):
            var = rank_dict[assembly_id, struct_asym_id]
            var[1] += 1
            assert var[1] <= var[0]
            return var[1]
        
        if isawaitable(path):
            path = await path
        path = Path(path)
        if path.suffix == '.cif':
            assg_oper_res = await cls.assg_oper_cif_base(path)
        elif path.suffix == '.json':
            assg_oper_res = await cls.assg_oper_json_base(path)
        else:
            raise ValueError(f'Invalid file: {path}')
        if assg_oper_res is None:
            return
        else:
            assg_df, oper_dict_1, oper_dict_2 = assg_oper_res
        assg_df.oper_expression = assg_df.oper_expression.apply(
            lambda x: cls.expandOperators(cls.parseOperatorList(x)))
        assg_df = split_df_by_chain(assg_df, assg_df.columns, ('oper_expression', ), mode='list').reset_index(drop=True)
        assg_df.oper_expression = assg_df.oper_expression.apply(
            lambda x: json.dumps(x).decode('utf-8'))
        
        assg_dict = assg_df.groupby('assembly_id').oper_expression.unique().to_dict()
        rank_dict = assg_df.groupby(['assembly_id', 'struct_asym_id']).struct_asym_id.count().to_dict()
        rank_dict = {key: [value, 0] for key, value in rank_dict.items()}
        assg_df['model_id'] = assg_df.apply(lambda x: np_where(assg_dict[x['assembly_id']] == x['oper_expression'])[0][0]+1, axis=1)
        # assg_df.apply(lambda x: assg_dict[x['assembly_id']].index(x['oper_expression'])+1, axis=1)
        
        assg_df['asym_id_rank'] = assg_df.apply(lambda x: to_rank(
            rank_dict, x['assembly_id'], x['struct_asym_id']), axis=1)
        assg_df['symmetry_operation'] = assg_df.oper_expression.apply(
            lambda x: [oper_dict_1[i] for i in json.loads(x)])
        assg_df.symmetry_operation = assg_df.symmetry_operation.apply(
            lambda x: json.dumps(x).decode('utf-8'))
        assg_df['symmetry_id'] = assg_df.oper_expression.apply(
            lambda x: [oper_dict_2[i] for i in json.loads(x)])
        assg_df.symmetry_id = assg_df.symmetry_id.apply(
            lambda x: json.dumps(x).decode('utf-8'))
        assg_df.assembly_id = assg_df.assembly_id.astype(int)

        return assg_df

    @unsync
    async def pipe_assg_data_collection(self) -> str:
        '''demo_dict = {"atom_site": [{"label_asym_id": "A", "label_seq_id": 23}]}'''
        res_df = await self.pdb_ob.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe)
        # res_dict = iter_first(res_df, lambda row: row.observed_ratio > 0)
        # assert res_dict is not None, f"{self.pdb_ob}: Unexpected Cases, without observed residues?!"
        res_dict = res_df.loc[res_df.observed_ratio.gt(0).idxmax()].to_dict()
        '''
        demo_dict = dict(atom_site=[dict(
            label_asym_id=res_dict['struct_asym_id'],
            label_seq_id=int(res_dict['residue_number']))])
        return json.dumps(demo_dict).decode('utf-8')
        '''
        return {col: res_dict[col] for col in ('struct_asym_id', 'residue_number')}

    @classmethod
    @unsync
    async def orr2eec(cls, path: Unfuture):
        '''
        TODO Retrieve entry-entity-chain info via `observed_residues_ratio`
        
        NOTE only for polymeric molecules 
        
            * polypeptide(L)
            * polypeptide(D) (e.g.3ue7, struct_asym_id: A)
            * polydeoxyribonucleotide
            * polyribonucleotide
            * polydeoxyribonucleotide/polyribonucleotide hybrid (e.g. 6ozg)
        '''
        eec_df = await path.then(cls.to_dataframe)
        eec_df['struct_asym_id_in_assembly'] = eec_df.struct_asym_id
        eec_df = eec_df.drop(columns=['number_residues', 'observed_ratio', 'partial_ratio'])
        return eec_df

    @unsync
    async def set_res2eec_df(self, merge_with_molecules_info:bool=True):
        '''
        NOTE except water

        Example of chain_id != struct_asym_id:

            * 6a8r
            * 6e32
        '''
        res_df = await self.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe)
        eec_df = res_df[['pdb_id', 'entity_id', 'chain_id',
                         'struct_asym_id']].drop_duplicates().reset_index(drop=True)
        # eec_df['struct_asym_id_in_assembly'] = eec_df.struct_asym_id
        if merge_with_molecules_info:
            mol_df = await self.fetch_from_pdbe_api('api/pdb/entry/molecules/', Base.to_dataframe)
            self.res2eec_df = eec_df.merge(mol_df, how='left')
        else:
            self.res2eec_df = eec_df

    @unsync
    async def get_res2eec_df(self, **kwargs):
        if hasattr(self, 'res2eec_df'):
            return self.res2eec_df
        else:
            await self.set_res2eec_df(**kwargs)
            return self.res2eec_df

    @classmethod
    @unsync
    async def assembly2eec(cls, path: Unfuture):
        '''
        NOTE except water
        '''
        assembly_df = await path.then(cls.to_dataframe)
        # exclude_molecule_type=('bound', 'water', 'carbohydrate polymer')
        assembly_focus_col = ('in_chains', 'pdb_id',
                              'entity_id', 'molecule_type', 
                              'assembly_id', 'details')
        eec_df = split_df_by_chain(
            assembly_df[assembly_df.molecule_type.ne('water')],
            assembly_focus_col, assembly_focus_col[0:1],
            'json-list').rename(columns={"in_chains": "struct_asym_id_in_assembly"}).reset_index(drop=True)
        return eec_df

    @staticmethod
    def to_assembly_id(pdb_id, assemblys):
        for assembly_id in assemblys:
            yield f"{pdb_id}/{assembly_id}"

    @unsync
    async def set_focus_assembly(self, focus_assembly_ids:Optional[Iterable[int]]=None, discard_multimer_chains_cutoff:Optional[int]=None, discard_multimer_chains_cutoff_for_au:Optional[int]=None, omit_peptide_length:int=0):
        '''
        EXAMPLE Discard multimer contains more then 16 chains (except water | bound | carbohydrate polymer)
        '''
        ass_eec_df = await self.fetch_from_pdbe_api('api/pdb/entry/assembly/', Base.to_dataframe)
        mol_df = await self.fetch_from_pdbe_api('api/pdb/entry/molecules/', Base.to_dataframe)
        mol_df = mol_df[~mol_df.molecule_type.isin(('water', 'bound', 'carbohydrate polymer'))]
        mol_df = mol_df[~(mol_df.molecule_type.isin(('polypeptide(L)', 'polypeptide(D)')) & mol_df.length.le(omit_peptide_length))]
        try:
            ass_eec_df = ass_eec_df[ass_eec_df.details.notnull() & ass_eec_df.entity_id.isin(mol_df.entity_id)]
        except AttributeError as e:
            warn(f"{repr(self)}: None assembly")
            raise e
        check_chain_num = ass_eec_df.groupby('assembly_id').in_chains.apply(lambda x: sum(i.count(',')+1 for i in x))
        if discard_multimer_chains_cutoff is not None:
            assemblys = set(ass_eec_df.assembly_id) & set(check_chain_num[check_chain_num<discard_multimer_chains_cutoff].index)
        else:
            assemblys = set(ass_eec_df.assembly_id)
        if discard_multimer_chains_cutoff_for_au is None or mol_df.in_struct_asyms.apply(lambda x: x.count(',')+1).sum() < discard_multimer_chains_cutoff_for_au:
            assemblys |= {0}
        if focus_assembly_ids is not None:
            assemblys = sorted(assemblys & set(int(i) for i in focus_assembly_ids))
        else:
            assemblys = sorted(assemblys)
        if not hasattr(self, 'assembly'):
            await self.set_assembly()
        self.focus_assembly: Dict[int, PDBAssembly] = {ass_id: ass_ob for ass_id, ass_ob in self.assembly.items() if ass_id in assemblys}

    @unsync
    async def set_assembly(self):
        '''
        NOTE even for NMR structures (e.g 1oo9), there exists assembly 1 for that entry
        NOTE Discard `details` is NaN -> not author_defined_assembly OR software_defined_assembly
        '''        
        ass_eec_df = await self.fetch_from_pdbe_api('api/pdb/entry/assembly/', Base.to_dataframe)
        ass_eec_df = ass_eec_df[ass_eec_df.details.notnull()]
        assemblys = set(ass_eec_df.assembly_id) | {0}
        self.assembly: Dict[int, PDBAssembly] = dict(zip(
            assemblys, 
            (PDBAssembly(ass_id).add_args(pdb_ob=self) for ass_id in self.to_assembly_id(self.pdb_id, assemblys))))

    def get_assembly(self, assembly_id):
        return self.assembly[assembly_id]
    
    @staticmethod
    def parseOperatorList(value: str) -> Iterable[Iterable[str]]:
        '''
        This function has been modified from `function parseOperatorList(value: string): string[][]`

            * Copyright (c) 2017 - now, Mol* contributors
            * https://github.com/molstar/molstar/
            * src/mol-model-formats/structure/property/assembly.ts

        NOTE (original notes from assembly.ts):

            * '(X0)(1-5)' becomes[['X0'], ['1', '2', '3', '4', '5']]
            * kudos to Glen van Ginkel.
        
        >>> PDB.parseOperatorList('(X0)(1-5)')
        >>> [['X0'], ['1', '2', '3', '4', '5']]
        '''

        def unit(g:str) -> Iterable[str]:
            group: Iterable[str] = []
            for e in g.split(','):
                try:
                    dashIndex = e.index('-')
                    group.extend(str(i) for i in range(
                        int(e[0:dashIndex]), int(e[dashIndex + 1:])+1))
                except ValueError:
                    group.append(e.strip())
            return group

        oeRegex = re_compile(r"\(?([^\(\)]+)\)?]*")
        return [unit(g) for g in oeRegex.findall(value)]

    @staticmethod
    def expandOperators(operatorList: Iterable[Iterable[str]]) -> Iterable[Iterable[str]]:
        '''
        This function has been modified from `function expandOperators(operatorList: string[][])`

            * Copyright (c) 2017 - now, Mol* contributors
            * https://github.com/molstar/molstar/
            * src/mol-model-formats/structure/property/assembly.ts
        
        >>> PDB.expandOperators([['X0'], ['1', '2', '3', '4', '5']])
        >>> [['X0', '1'], ['X0', '2'], ['X0', '3'], ['X0', '4'], ['X0', '5']]
        '''
        ops: Iterable[Iterable[str]] = []
        currentOp: Iterable[str] = ['' for _ in range(len(operatorList))]
        PDB.expandOperators1(operatorList, ops, len(operatorList) - 1, currentOp)
        return ops
    
    @staticmethod
    def expandOperators1(operatorNames: Iterable[Iterable[str]], lyst: Iterable[Iterable[str]], i: int, current: Iterable[str]) -> Iterable[Iterable[str]]:
        '''
        This function has been modified from `function expandOperators1(operatorNames: string[][], list: string[][], i: number, current: string[])`

            * Copyright (c) 2017 - now, Mol* contributors
            * https://github.com/molstar/molstar/
            * src/mol-model-formats/structure/property/assembly.ts
        '''
        if i < 0:
            lyst.append(current[0:])
            return

        ops = operatorNames[i]
        for j in range(len(ops)):
            current[i] = ops[j]
            PDB.expandOperators1(operatorNames, lyst, i - 1, current)

    @staticmethod
    def dumpsOperators(ops: Iterable[Iterable[str]], sep1:str=',', sep2:str='&') -> str:
        return sep1.join(sep2.join(i) for i in ops)

    @unsync
    async def stats_protein_entity_seq(self):
        store = await self.sqlite_api.StatsProteinEntitySeq.objects.filter(pdb_id=self.pdb_id).all()
        if not store:
            muta_df = await Base.to_dataframe_with_kwargs(
                self.fetch_from_pdbe_api('api/pdb/entry/mutated_AA_or_NA/'),
                usecols=lambda x: not x.startswith('author') and not x in ('struct_asym_id', 'chain_id'))
            if muta_df is not None:
                muta_df = muta_df.drop_duplicates().reset_index(drop=True)
                muta_df = flat_dict_in_df(muta_df, muta_df.mutation_details.apply(json.loads), ['type', 'from','to'])
                muta_dict = muta_df[muta_df['mutation_details.type'].isin(
                    ('Cloning artifact', 'Expression tag', 'Initiating methionine', 'Linker'))].groupby(
                    'entity_id', group_keys=False).residue_number.apply(to_interval).to_dict()
            else:
                muta_dict = dict()
            
            def regist_info(record):
                if record['molecule_type'] not in ('polypeptide(L)', 'polypeptide(D)'):
                    return
                seq = array(list(len(i) if i != '(UNK)' else
                                -1 for i in self.protein_sequence_pat.findall(record['pdb_sequence'])))
                
                std_mask = seq == 1
                std_res = np_where(std_mask)[0]+1
                non_res = np_where(~std_mask)[0]+1
                unk_res = np_where(seq == -1)[0]+1

                if not isinstance(record['ca_p_only'], bool):
                    warn(f'{record["pdb_id"]}: molecules', MultiWrittenWarning)
                    raise AssertionError(f"{record}")
                
                store.append(self.StatsProteinEntitySeq(
                    record['pdb_id'],
                    record['molecule_type'],
                    record['entity_id'],
                    record['ca_p_only'],
                    len(seq),
                    to_interval(std_res),
                    len(std_res),
                    to_interval(non_res),
                    len(non_res),
                    to_interval(unk_res),
                    len(unk_res),
                    muta_dict.get(record['entity_id'], tuple())))
                # assert len(seq)-len(std_res) == len(non_res)
            
            mol_df = await self.fetch_from_pdbe_api('api/pdb/entry/molecules/', Base.to_dataframe)
            mol_df.apply(regist_info, axis=1)
            if store:
                await self.sqlite_api.async_insert(self.sqlite_api.StatsProteinEntitySeq, store)
        return store

    @unsync
    async def stats_nucleotide_entity_seq(self):
        store = await self.sqlite_api.StatsNucleotideEntitySeq.objects.filter(pdb_id=self.pdb_id).all()
        if not store:
            def regist_info(record):
                if record['molecule_type'] not in ('polydeoxyribonucleotide', 'polyribonucleotide', 'polydeoxyribonucleotide/polyribonucleotide hybrid'):
                    return
                seq = array(list(len(i) if i != '(UNK)' else
                                -1 for i in self.nucleotide_sequence_pat.findall(record['pdb_sequence'])))
                rna_mask = seq == 1
                dna_mask = seq == 4
                unk_mask = seq == -1
                non_mask = ~(dna_mask | rna_mask)
                rna_res = np_where(rna_mask)[0]+1
                dna_res = np_where(dna_mask)[0]+1
                non_res = np_where(non_mask)[0]+1
                unk_res = np_where(unk_mask)[0]+1

                store.append(self.StatsNucleotideEntitySeq(
                    record['pdb_id'],
                    record['molecule_type'],
                    record['entity_id'],
                    record['ca_p_only'],
                    len(seq),
                    to_interval(dna_res),
                    len(dna_res),
                    to_interval(rna_res),
                    len(rna_res),
                    to_interval(non_res),
                    len(non_res),
                    to_interval(unk_res),
                    len(unk_res)))

            mol_df = await self.fetch_from_pdbe_api('api/pdb/entry/molecules/', Base.to_dataframe)
            mol_df.apply(regist_info, axis=1)
            if store:
                await self.sqlite_api.async_insert(self.sqlite_api.StatsNucleotideEntitySeq, store)
        return store

    '''
    @unsync
    async def stats_poly_chain_obs_seq(self):
        store = defaultdict(list)
        def regist_info(record):
            store[(record['pdb_id'], record['entity_id'], record['chain_id'], record['struct_asym_id'])
            ].append(
                (json.loads(record['start'])['residue_number'], json.loads(record['end'])['residue_number']))

        c_df = await self.fetch_from_pdbe_api('api/pdb/entry/polymer_coverage/', Base.to_dataframe)
        c_df.apply(regist_info, axis=1)
        return [self.StatsChainSeq(*key, value, range_len(value)) for key, value in store.items()]
    '''
    
    @staticmethod
    def ret_obs_res(value):
        obs_index = []
        obs_array = []
        for index, ratio in value:
            if ratio > 0:
                obs_index.append(index)
            obs_array.append(ratio)
        return to_interval(obs_index), len(obs_index), zlib.compress(json.dumps(obs_array))

    @unsync
    async def stats_chain_seq(self):
        store = await self.sqlite_api.StatsChainSeq.objects.filter(pdb_id=self.pdb_id).all()
        if store:
            return store

        store = defaultdict(list)
        
        def regist_info(record):
            store[(record['pdb_id'], record['entity_id'], record['chain_id'], record['struct_asym_id'])
                  ].append((record['residue_number'], record['observed_ratio']))
        
        res_df = await self.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe)
        res_df.sort_values(by=['struct_asym_id', 'residue_number'], inplace=True)
        # mol_df = await self.fetch_from_pdbe_api('api/pdb/entry/molecules/', Base.to_dataframe)
        # res_df = res_df.merge(mol_df[mol_df.molecule_type.isin(('polypeptide(L)', 'polypeptide(D)'))][['entity_id']])
        res_df.apply(lambda x: regist_info(x), axis=1)

        summary_dict = (await self.fetch_from_pdbe_api('api/pdb/entry/summary/', a_load_json, json=True)
                            )[self.pdb_id][0]['number_of_entities']
        # TODO: make sure whether bound_count including sugar & other
        bound_count = sum(value for key, value in summary_dict.items() if key in ('ligand', 'sugar', 'other', 'carbohydrate_polymer'))
        if bound_count:
            br_df = await Base.to_dataframe_with_kwargs(
                self.fetch_from_pdbe_api('api/pdb/entry/binding_sites/'),
                usecols=lambda x: not x.startswith('author'))
        else:
            br_df = None
        if br_df is None:
            ret = [self.StatsChainSeq(*key,
                                      *self.ret_obs_res(value),
                                       tuple(),
                                       0
                                       ) for key, value in store.items()]
        else:
            br_dict = br_df[br_df.entity_id.ne(-1)].groupby('struct_asym_id').residue_number.apply(lambda x: (to_interval(x), len(frozenset(x)))).to_dict()
            ret = [self.StatsChainSeq(*key, 
                                      *self.ret_obs_res(value),
                                      *br_dict.get(key[-1], (tuple(), 0))
                                      ) for key, value in store.items()]

        await self.sqlite_api.async_insert(self.sqlite_api.StatsChainSeq, ret)
        return ret

    @unsync
    async def stats_chain(self, stats_nucleotide=False):
        pe_df = DataFrame(await self.stats_protein_entity_seq())
        c_df = DataFrame(await self.stats_chain_seq())
        try:
            pec_df = c_df.merge(pe_df)
            pec_df['OBS_STD_INDEX'] = pec_df.apply(lambda x: overlap_range(x['OBS_INDEX'], x['STD_INDEX']), axis=1)
            pec_df['OBS_STD_COUNT'] = pec_df.OBS_STD_INDEX.apply(range_len)
        except Exception:
            pec_df = None
        if stats_nucleotide:
            ne_df = DataFrame(await self.stats_nucleotide_entity_seq())
            nec_df = c_df.merge(ne_df) if len(ne_df) > 0 else None
        else:
            nec_df = None
        return pec_df, nec_df

    def deal_with_pdb_sequence_indices_with_multiple_residues(self, pdb_sequence, data, kwargs):
        '''
        Deal with pdb_sequence_indices_with_multiple_residues
            4u2v entity 1 -> {"68":{"three_letter_code":"CR2","parent_chem_comp_ids":["GLY","TYR","GLY"],"one_letter_code":"GYG"}}
        '''
        if isinstance(data, float):
            warn(f'{repr(self)}: Unexpected float datatype for `pdb_sequence_indices_with_multiple_residues`')
            assert '(' not in pdb_sequence, repr(self)
            return pdb_sequence
        elif isinstance(data, str):
            if data == '{}':
                return pdb_sequence
            else:
                data = json.loads(data)
        for val_dict in data.values():
            one_letter_code = val_dict["one_letter_code"]
            three_letter_code = val_dict["three_letter_code"]
            if len(one_letter_code) != 1:
                warn(f"Possible Peptide Linking: {repr(self)}, {kwargs}, {val_dict}; select the first code", PeptideLinkingWarning)
            pdb_sequence = pdb_sequence.replace(f'({three_letter_code})', one_letter_code[0])
        return pdb_sequence

    @unsync
    async def get_sequence(self, mode='fix_seq', **kwargs):
        '''
        Get true SEQRES Sequence via entity_id | chain_id (default protein) | struct_asym_id
        '''
        mol_df = await self.fetch_from_pdbe_api('api/pdb/entry/molecules/', Base.to_dataframe)
        if mol_df is None:
            raise PossibleObsoletedPDBEntryError(f"None dataframe: {repr(self)}, either obsoleted or API lag update")
        if 'entity_id' in kwargs:
            mask = mol_df.entity_id.eq(kwargs['entity_id'])
        elif 'struct_asym_id' in kwargs:
            struct_asym_id = kwargs['struct_asym_id']
            mask = mol_df.in_struct_asyms.apply(lambda x: f'"{struct_asym_id}"' in x)
        elif 'chain_id' in kwargs:
            chain_id = kwargs['chain_id']
            mask = mol_df.molecule_type.eq('polypeptide(L)') & mol_df.in_chains.apply(lambda x: f'"{chain_id}"' in x)
        else:
            raise ValueError(f"Cannot get sequence with specified information: {kwargs}")
        if not mask.any():
            raise ValueError(f"Cannot get sequence with specified information: {kwargs}")
        else:
            cur_record = mol_df.loc[mask.idxmax()]
        if not cur_record['molecule_type'].startswith('poly'):
            raise ValueError(f"Not polymer with specified information: {repr(self)} -> {kwargs}")
        if mode == 'fix_seq':
            return self.deal_with_pdb_sequence_indices_with_multiple_residues(
                cur_record['pdb_sequence'], cur_record['pdb_sequence_indices_with_multiple_residues'], kwargs)
        elif mode == 'raw_seq':
            return cur_record['sequence']
        elif mode == 'raw_pdb_seq':
            return cur_record['pdb_sequence']
        elif mode == 'mod_x_seq':
            return ''.join(i if len(i) == 1 else 'X' for i in self.protein_sequence_pat.findall(cur_record['pdb_sequence']))
        else:
            return

    @unsync
    @retry(wait=wait_random(max=20), stop=stop_after_attempt(5))
    async def profile_asym_id_in_assembly_unit(self, parent, assembly_id, replace):
        header = f'{self.pdb_id}-assembly-{assembly_id}'
        target_file = parent/f'{header}-pdbe_chain_remapping.cif'
        exists, _ = await aio_file_exists_stat(target_file)
        if (not exists) or replace:
            await cif_gz_stream.writer(
                cif_gz_stream.reader(f'http://www.ebi.ac.uk/pdbe/static/entry/download/{header}.cif.gz'),
                target_file,
                b'data_%s\n#\nloop_\n' % bytes(header, 'utf-8'),
                b'_pdbe_chain_remapping')
        try:
            async with aiofiles_open(target_file, 'rt') as file_io:
                handle = await file_io.readlines()
                mmcif_dict = MMCIF2DictPlus(iline for iline in handle)
                if (len(mmcif_dict) == 0) or (not handle[-1] == '#\n'):
                    raise PossibleConnectionError(target_file)
                del mmcif_dict['data_']
                cur_df = DataFrame(mmcif_dict)
        except (KeyError, ValueError, PossibleConnectionError) as e:
            warn(f"{target_file}, {e.__class__.__name__}, need retry", WithoutCifKeyWarning)
            await aiofiles_os_remove(target_file)
            raise e
        '''
        mmcif_dict = await cif_gz_stream.full_io(
            f'http://www.ebi.ac.uk/pdbe/static/entry/download/{header}.cif.gz',
            target_file)
        del mmcif_dict['data_']
        cur_df = DataFrame(mmcif_dict)
        if len(cur_df) == 0:
            raise PossibleConnectionError(target_file)
        '''
        cur_df['_pdbe_chain_remapping.assembly_id'] = assembly_id  # NOTE: exception: 6E5N-assembly-1
        return cur_df

    @unsync
    async def profile_asym_id_in_assembly(self, assembly_ids, replace:bool=False):
        parent = self.assembly_cif_folder/self.pdb_id[1:3]
        parent.mkdir(parents=True, exist_ok=True)
        dfs = [await self.pipe_register_task(
            f'{self.pdb_id}-assembly-{assembly_id}',
            self.profile_asym_id_in_assembly_unit,
            parent=parent, assembly_id=assembly_id, replace=replace) for assembly_id in assembly_ids]
        in_ass_df =  concat(dfs, ignore_index=True, sort=False)
        in_ass_df.columns = [i.replace('_pdbe_chain_remapping.', '') for i in in_ass_df.columns]
        in_ass_df = in_ass_df[['assembly_id', 'entity_id', 'orig_label_asym_id', 'new_label_asym_id', 'applied_operations']].rename(
            columns={'orig_label_asym_id': 'struct_asym_id',
                     'new_label_asym_id': 'struct_asym_id_in_assembly',
                     'applied_operations': 'oper_expression'}).copy()
        in_ass_df.oper_expression = in_ass_df.oper_expression.apply(lambda x: json.dumps([i for i in x.split('_')[::-1] if i != '']).decode('utf-8'))
        # in_ass_df.assembly_id = in_ass_df.assembly_id.astype(int)  # NOTE: exception: 6E5N-assembly-1
        in_ass_df.entity_id = in_ass_df.entity_id.astype(int)
        return in_ass_df

    @unsync
    async def rd_source_ass_oper_df(self):
        upper_id = self.pdb_id.upper()
        try:
            assembly_ids = (await self.fetch_from_rcsb_api('graphql', 
                query='{entry(entry_id:"%s"){rcsb_entry_container_identifiers{assembly_ids}}}' % upper_id, 
                then_func=a_load_json,
                json=True))['data']['entry']['rcsb_entry_container_identifiers']['assembly_ids']
        except TypeError:
            warn(repr(self), PossibleObsoletedPDBEntryWarning)
            assembly_ids = None
        if assembly_ids is None:
            return
        query = '''{
            assemblies(assembly_ids: %s) {
                pdbx_struct_assembly_gen{
                    assembly_id
                    asym_id_list
                    oper_expression
                    ordinal}
                pdbx_struct_oper_list{
                    id
                    name
                    symmetry_operation}
                rcsb_assembly_info{
                    entry_id
                    assembly_id}}}
            ''' % json.dumps([f'{upper_id}-{assembly_id}' for assembly_id in assembly_ids]).decode('utf-8')
        return await self.fetch_from_rcsb_api('graphql', query=query).then(self.to_assg_oper_df)

    @unsync
    async def ms_source_ass_oper_df(self, struct_asym_id, residue_number):
        assg_oper_df = await self.fetch_from_modelServer_api('atoms',
            data_collection=json.dumps(dict(atom_site=[dict(
                label_asym_id=struct_asym_id,
                label_seq_id=int(residue_number))])).decode('utf-8'),
            then_func=self.to_assg_oper_df)
        return assg_oper_df
    
    @unsync
    async def cs_source_ass_oper_df(self, struct_asym_id, residue_number):
        assg_oper_df = await self.fetch_from_coordinateServer_api(
            'residues', 
            then_func=self.to_assg_oper_df,
            #root='ebi',
            asymId=struct_asym_id, 
            seqNumber=int(residue_number), 
            modelId=1)
        return assg_oper_df
    
    @unsync
    async def set_subset_assembly_from_df(self, profile_id_df):
        mask = profile_id_df.assembly_id.ne(0)
        if not mask.any():
            self.subset_assembly = frozenset()
        else:
            if not hasattr(self, 'assembly'):
                await self.set_assembly()
            to_fetch_assembly = profile_id_df[
                mask & 
                profile_id_df.symmetry_id.ne('["1_555"]') &
                profile_id_df.assembly_id.isin(self.assembly.keys())
                ].assembly_id.unique()
            self.subset_assembly = frozenset(profile_id_df.assembly_id) - frozenset(to_fetch_assembly) - {0}
            self.subset_assembly = {ass: profile_id_df[profile_id_df.assembly_id.eq(ass)].struct_asym_id.tolist() for ass in self.subset_assembly}

    @unsync
    async def profile_id_from_DB(self):
        task = self.tasks.get((repr(self), 'profile_id'), None)
        if task is None:
            task = await self.sqlite_api.profile_id.objects.filter(pdb_id=self.pdb_id).all()
            if not task:
                return
            task = DataFrame(task).sort_values(['assembly_id', 'model_id', 'entity_id', 'struct_asym_id']).reset_index(drop=True)
        if not hasattr(self, 'subset_assembly'):
            await self.set_subset_assembly_from_df(task)
        self.register_task((repr(self), 'profile_id'), task)
        return task
    
    @unsync
    async def profile_id(self):
        '''
        NOTE: except water
        TODO: CLEAN CODE
        '''
        profile_id_df = await self.profile_id_from_DB()
        if profile_id_df is not None:
            return profile_id_df
        
        if choice((1, 0, 0, 0, 0)):
            assg_oper_df = await self.rd_source_ass_oper_df()
        else:
            demo_dict = await self.pipe_assg_data_collection()
            assg_oper_df = await self.ms_source_ass_oper_df(**demo_dict)  # self.cs_source_ass_oper_df
        # assg_oper_df = await self.ms_source_ass_oper_df(**demo_dict)
        res2eec_df = await self.get_res2eec_df()
        focus_res2eec_df = res2eec_df[['pdb_id', 'entity_id', 'molecule_type', 'chain_id', 'struct_asym_id']]
        add_0_assg_oper_df = focus_res2eec_df.copy()
        add_0_assg_oper_df['assembly_id'] = 0
        add_0_assg_oper_df['model_id'] = 1
        add_0_assg_oper_df['asym_id_rank'] = 1
        add_0_assg_oper_df['oper_expression'] = ''
        add_0_assg_oper_df['symmetry_operation'] = ''
        add_0_assg_oper_df['symmetry_id'] = ''

        if assg_oper_df is not None:
            focus_assg_oper_df = assg_oper_df[assg_oper_df.struct_asym_id.isin(focus_res2eec_df.struct_asym_id)]
            new_focus_assg_oper_df = concat([add_0_assg_oper_df, focus_assg_oper_df.merge(focus_res2eec_df, how='left')], ignore_index=True, sort=False)
            assert not new_focus_assg_oper_df.isnull().sum().any(), f"Unexpected Cases {new_focus_assg_oper_df}"
            if not hasattr(self, 'assembly'):
                await self.set_assembly()
            new_focus_assg_oper_df = new_focus_assg_oper_df[new_focus_assg_oper_df.assembly_id.isin(self.assembly.keys())].reset_index(drop=True)
        else:
            add_0_assg_oper_df['struct_asym_id_in_assembly'] = add_0_assg_oper_df.struct_asym_id
            add_0_assg_oper_df['details'] = 'asymmetric_unit'
            self.subset_assembly = frozenset()
            if not hasattr(self, 'assembly'):
                await self.set_assembly()
            await self.sqlite_api.async_insert(self.sqlite_api.profile_id, add_0_assg_oper_df.to_dict('records'))
            self.register_task((repr(self), 'profile_id'), add_0_assg_oper_df)
            return add_0_assg_oper_df
        
        to_fetch_assembly = new_focus_assg_oper_df[
            new_focus_assg_oper_df.assembly_id.ne(0) & new_focus_assg_oper_df.symmetry_id.ne('["1_555"]')].assembly_id.unique()
        if len(to_fetch_assembly)> 0:
            in_ass_df = await self.profile_asym_id_in_assembly(to_fetch_assembly)
            profile_id_df = new_focus_assg_oper_df.merge(in_ass_df, how='left')
            self.subset_assembly = frozenset(profile_id_df.assembly_id)-frozenset(to_fetch_assembly)
            rest = frozenset(profile_id_df[profile_id_df.struct_asym_id_in_assembly.isnull()].assembly_id)
            assert rest == self.subset_assembly, f"{repr(self)}\n{rest}\n{self.subset_assembly}\n{to_fetch_assembly}"
            a0_index = profile_id_df[profile_id_df.assembly_id.isin(self.subset_assembly)].index
            self.subset_assembly = self.subset_assembly - {0}
            profile_id_df.loc[a0_index, 'struct_asym_id_in_assembly'] = profile_id_df.loc[a0_index, 'struct_asym_id']
        else:
            self.subset_assembly = frozenset(new_focus_assg_oper_df.assembly_id) - {0}
            new_focus_assg_oper_df['struct_asym_id_in_assembly'] = new_focus_assg_oper_df.struct_asym_id
            profile_id_df = new_focus_assg_oper_df

        self.subset_assembly = {ass: profile_id_df[profile_id_df.assembly_id.eq(ass)].struct_asym_id.tolist() for ass in self.subset_assembly}
        profile_id_df['au_subset'] = False
        profile_id_df.loc[profile_id_df[profile_id_df.assembly_id.isin(self.subset_assembly)].index, 'au_subset'] = True

        # TODO: Check
        ass_df = await self.fetch_from_pdbe_api('api/pdb/entry/assembly/', Base.to_dataframe)
        a_df = profile_id_df[['assembly_id', 'entity_id']].drop_duplicates().query('assembly_id != 0')
        for assembly_id, entity_id in zip(a_df.assembly_id, a_df.entity_id):
            profile_lyst = tuple(profile_id_df[profile_id_df.assembly_id.eq(assembly_id) & profile_id_df.entity_id.eq(entity_id)].struct_asym_id_in_assembly.sort_values().tolist())
            ass_lyst = tuple(json.loads(ass_df.loc[(ass_df.assembly_id.eq(assembly_id)&ass_df.entity_id.eq(entity_id)).idxmax(), 'in_chains']))
            assert profile_lyst == ass_lyst, f"\n{self.pdb_id}\n{assembly_id}\n{entity_id}\n{profile_lyst},\n{ass_lyst}"
        profile_id_df = profile_id_df.merge(ass_df[['assembly_id', 'details']].drop_duplicates(), how='left')
        profile_id_df.loc[profile_id_df[profile_id_df.assembly_id.eq(0)].index, 'details'] = 'asymmetric_unit'
        await self.sqlite_api.async_insert(self.sqlite_api.profile_id, profile_id_df.to_dict('records'))
        self.register_task((repr(self), 'profile_id'), profile_id_df)
        return profile_id_df
    
    @unsync
    async def get_expanded_map_res_df(self, UniProt, unp_range, pdb_range, conflict_pdb_index=None, **kwargs):
        assert kwargs.keys() <= frozenset({'chain_id', 'entity_id', 'struct_asym_id', 'pdb_id'})
        res_map_df = DataFrame(zip(expand_interval(unp_range),expand_interval(pdb_range)), columns=['unp_residue_number', 'residue_number'])
        res_map_df['UniProt'] = UniProt
        res_df = await self.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe)
        res_map_df_full = res_map_df.merge(res_df.query(
            'and '.join(f'{key}=="{value}"' if isinstance(value, str) else f'{key}=={value}' for key,value in kwargs.items())))
        assert res_map_df_full.shape[0] == res_map_df.shape[0], f"{(repr(self), UniProt, unp_range, pdb_range, kwargs)}"
        if conflict_pdb_index is None:
            return res_map_df_full
        else:
            cf = self.convert_conflict_pdb_index(conflict_pdb_index)
            if len(cf) > 0:
                return res_map_df_full.merge(cf, how='left')
            else:
                return res_map_df_full
    
    @staticmethod
    def three_range2range_df(unp_range, pdb_range, auth_pdb_range):
        unp_range = DataFrame(unp_range, columns=['unp_beg', 'unp_end'])
        pdb_range = DataFrame(pdb_range, columns=['pdb_beg', 'pdb_end'])
        auth_pdb_range = DataFrame(auth_pdb_range, columns=['auth_pdb_beg', 'auth_pdb_end'])
        return concat([unp_range, pdb_range, auth_pdb_range], axis=1)

    @unsync
    async def get_ranged_map_res_df(self, UniProt, unp_range, pdb_range, conflict_pdb_index=None, **kwargs):
        assert kwargs.keys() <= frozenset({'chain_id', 'entity_id', 'struct_asym_id', 'pdb_id'})
        res_df = (await self.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe)).query(
            'and '.join(f'{key}=="{value}"' if isinstance(value, str) else f'{key}=={value}' for key,value in kwargs.items()))
        assert res_df.shape[0] > 0 # , f"Cannot find expected chain: {kwargs}"
        res_map_df = DataFrame(zip(expand_interval(unp_range),expand_interval(pdb_range)), columns=['unp_residue_number', 'residue_number'])
        res_map_df_full = res_map_df.merge(res_df)
        assert res_map_df_full.shape[0] == res_map_df.shape[0] #, f"{(repr(self), UniProt, unp_range, pdb_range, kwargs)}"
        if conflict_pdb_index is None:
            res_map_df_full['conflict_code'] = nan
        else:
            cf = self.convert_conflict_pdb_index(conflict_pdb_index)
            if len(cf) > 0:
                res_map_df_full = res_map_df_full.merge(cf, how='left')
            else:
                res_map_df_full['conflict_code'] = nan
        mask = (
            res_map_df_full.author_insertion_code.ne('') |
            res_map_df_full.multiple_conformers.notnull() |
            res_map_df_full.observed_ratio.lt(1) |
            (~res_map_df_full.residue_name.isin(aa_three2one)) |
            res_map_df_full.conflict_code.notnull())
        demo_record = res_map_df_full.iloc[0]
        if (~mask).all():
            range_df = self.three_range2range_df(*lyst32interval(
                res_map_df_full.unp_residue_number, res_map_df_full.residue_number, res_map_df_full.author_residue_number))
            for col in ('pdb_id', 'entity_id', 'struct_asym_id', 'chain_id'):
                range_df[col] = demo_record[col]
            range_df['UniProt'] = UniProt
            range_df['author_insertion_code'] = ''
            range_df['residue_name'] = ''
            range_df['observed_ratio'] = 1.0
            range_df['multiple_conformers'] = nan
            range_df['conflict_code'] = nan
            return range_df
        elif mask.all():
            res_map_df_full['UniProt'] = UniProt
            final_df = res_map_df_full.rename(
                columns={'unp_residue_number': 'unp_beg', 'residue_number': 'pdb_beg', 'author_residue_number': 'auth_pdb_beg'})
            final_df['unp_end'] = final_df.unp_beg
            final_df['pdb_end'] = final_df.pdb_beg
            final_df['auth_pdb_end'] = final_df.auth_pdb_beg
            return final_df
        else:
            to_range_df = res_map_df_full[~mask]
            range_df = self.three_range2range_df(*lyst32interval(
                to_range_df.unp_residue_number, to_range_df.residue_number, to_range_df.author_residue_number))
            one_df = res_map_df_full[mask].rename(
                columns={'unp_residue_number': 'unp_beg', 'residue_number': 'pdb_beg', 'author_residue_number': 'auth_pdb_beg'})
            one_df['unp_end'] = one_df.unp_beg
            one_df['pdb_end'] = one_df.pdb_beg
            one_df['auth_pdb_end'] = one_df.auth_pdb_beg
            final_df = concat([one_df, range_df], axis=0, ignore_index=True)
            final_df['UniProt'] = UniProt
            final_df.pdb_id.fillna(demo_record['pdb_id'], inplace=True)
            final_df.entity_id = final_df.entity_id.fillna(demo_record['entity_id']).astype(int)
            final_df.chain_id.fillna(demo_record['chain_id'], inplace=True)
            final_df.struct_asym_id.fillna(
                demo_record['struct_asym_id'], inplace=True)
            final_df.author_insertion_code.fillna('', inplace=True)
            final_df.observed_ratio.fillna(1.0, inplace=True)
            final_df.residue_name.fillna('', inplace=True)
            return final_df

    @staticmethod
    def convert_conflict_pdb_index(conflict_pdb_index):
        if isinstance(conflict_pdb_index, str):
            conflict_pdb_index = json.loads(conflict_pdb_index)
        conflict_pdb_index = OrderedDict(conflict_pdb_index)
        return DataFrame({'residue_number': (int(i) for i in conflict_pdb_index.keys()), 'conflict_code': (i for i in conflict_pdb_index.values())})

    @unsync
    async def get_map_res_df(self, UniProt, unp_range, pdb_range, your_sites, conflict_pdb_index=None, unp2pdb:bool=True, author_site:bool=False, already_merged_with_residue_listing:bool=False,**kwargs):
        assert kwargs.keys() <= frozenset({'chain_id', 'entity_id', 'struct_asym_id', 'pdb_id'})
        if not isinstance(your_sites, Iterable) or len(your_sites) == 0:
            return
        if unp2pdb:
            info = {
                'unp_residue_number': your_sites,
                'residue_number': (SIFTS.convert_index(pdb_range, unp_range, your_sites)),
                'UniProt': UniProt}
            your_sites_df = DataFrame({**info, **kwargs})
            your_sites_df = your_sites_df[~your_sites_df.residue_number.isnull()]
            if len(your_sites_df) == 0:
                return
            ret = (await self.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe)).merge(your_sites_df)
        else:
            if already_merged_with_residue_listing:
                ret = your_sites.copy()
                ret['UniProt'] = UniProt
                ret['unp_residue_number'] = SIFTS.convert_index(unp_range, pdb_range, ret.residue_number)
                ret = ret[~ret.unp_residue_number.isnull()]
            else:
                if author_site:
                    author_residue_numbers, author_insertion_codes = zip(
                        *(self.author_residue_pat.fullmatch(i).groups() for i in your_sites))
                    info = {
                        'author_residue_number': (int(i) for i in author_residue_numbers),
                        'author_insertion_code': author_insertion_codes,
                        'UniProt': UniProt}
                    ret = DataFrame({**info, **kwargs})
                    ret = (await self.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe)).merge(ret)
                    ret['unp_residue_number'] = SIFTS.convert_index(unp_range, pdb_range, ret.residue_number)
                    ret = ret[~ret.unp_residue_number.isnull()]
                else:
                    info = {
                        'residue_number': your_sites,
                        'unp_residue_number': SIFTS.convert_index(unp_range, pdb_range, your_sites),
                        'UniProt': UniProt}
                    ret = DataFrame({**info, **kwargs})
                    ret = ret[~ret.unp_residue_number.isnull()]
                    ret = (await self.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe)).merge(ret)
        if ret.shape[0] == 0:
            return None
        if conflict_pdb_index is None:
            return ret
        else:
            cf = self.convert_conflict_pdb_index(conflict_pdb_index)
            if len(cf) > 0:
                return ret.merge(cf, how='left')
            else:
                return ret

    @unsync
    async def pipe_interface_res_dict_ic(self, include_chains=None, focus_assembly_ids=None, func='set_interface', discard_multimer_chains_cutoff=21, discard_multimer_chains_cutoff_for_au=None, omit_peptide_length:int=20, css_cutoff=-1, chains_count_cutoff=2, allow_not_polymer=False, **kwargs):
        # ic: include_chains version
        if include_chains is not None:
            include_chains = include_chains[self.pdb_id]
        await self.set_focus_assembly(focus_assembly_ids=focus_assembly_ids, discard_multimer_chains_cutoff=discard_multimer_chains_cutoff, discard_multimer_chains_cutoff_for_au=discard_multimer_chains_cutoff_for_au, omit_peptide_length=omit_peptide_length)
        res = []
        for assembly in self.focus_assembly.values():
            try:
                await getattr(assembly, func)(**kwargs)
            except PossibleInvalidAssemblyError:
                warn(f"skip {repr(assembly)}", SkipAssemblyWarning)
                continue
            if len(assembly.interface) == 0:
                continue
            if include_chains is not None:
                tr = await assembly.get_assemble_eec_as_df()
                cur_include_chains = frozenset(tr[tr.struct_asym_id.isin(include_chains)].struct_asym_id_in_assembly)
            else:
                cur_include_chains = include_chains
            for interface in assembly.interface.values():
                if ((cur_include_chains is None) or bool(interface.chain_set & cur_include_chains)) and interface.info['css'] > css_cutoff and len(interface.chain_set) >= chains_count_cutoff:
                    if func == 'pipe_protein_protein_interface' or func == 'pipe_protein_nucleotide_interface':
                        if not allow_not_polymer and not (interface.info['chains'][0][1] is None and interface.info['chains'][1][1] is None):
                            continue
                    elif func == 'pipe_protein_ligand_interface':
                        if not allow_not_polymer and not (interface.info['chains'][0][1] is None or interface.info['chains'][1][1] is None):
                            continue
                    res.append(interface)
        return res

    @unsync
    async def pipe_interface_res_info_for_rcsb(self, chain_pairs=None, include_chains=None, focus_assembly_ids=None, discard_multimer_chains_cutoff=21, discard_multimer_chains_cutoff_for_au=None, omit_peptide_length:int=20, **kwargs):
        '''
        NOTE: `chain_pairs` would suppress `include_chains`
        '''
        if chain_pairs is not None:
            chain_pairs = chain_pairs[self.pdb_id]
        elif include_chains is not None:
            include_chains = include_chains[self.pdb_id]
        await self.set_focus_assembly(focus_assembly_ids=focus_assembly_ids, discard_multimer_chains_cutoff=discard_multimer_chains_cutoff, discard_multimer_chains_cutoff_for_au=discard_multimer_chains_cutoff_for_au, omit_peptide_length=omit_peptide_length)
        interface_df_list = []
        for assembly in self.focus_assembly.values():
            if assembly.assembly_id == 0:
                continue
            interface_df = await assembly.get_interface_info_from_rcsb_data_api()
            if interface_df is None:
                continue
            if chain_pairs is not None:
                # NOTE: discard include_chains anyway
                mask = array([(frozenset(item) in chain_pairs) for item in zip(interface_df.struct_asym_id_1,interface_df.struct_asym_id_2)])
                if mask.any():
                    interface_df_list.append(interface_df[mask])
            else:
                if include_chains is None:
                    interface_df_list.append(interface_df)
                else:
                    mask = interface_df.struct_asym_id_1.isin(include_chains) | interface_df.struct_asym_id_2.isin(include_chains)
                    if mask.any():
                        interface_df_list.append(interface_df[mask])
        if interface_df_list:
            return concat(interface_df_list, axis=0, ignore_index=True, sort=False)
        else:
            return

    @unsync
    async def pipe_interface_res_dict(self, chain_pairs=None, focus_assembly_ids=None, func='set_interface', discard_multimer_chains_cutoff=21, discard_multimer_chains_cutoff_for_au=None, omit_peptide_length:int=20, css_cutoff=-1, chains_count_cutoff=2, allow_not_polymer=False, **kwargs):
        # maybe the name `au2bu` should be changed since its actual behavior is to use copied chains
        if chain_pairs is not None:
            chain_pairs = chain_pairs[self.pdb_id]
        await self.set_focus_assembly(focus_assembly_ids=focus_assembly_ids, discard_multimer_chains_cutoff=discard_multimer_chains_cutoff, discard_multimer_chains_cutoff_for_au=discard_multimer_chains_cutoff_for_au, omit_peptide_length=omit_peptide_length)
        res = []
        for assembly in self.focus_assembly.values():
            try:
                await getattr(assembly, func)(**kwargs)
            except PossibleInvalidAssemblyError:
                warn(f"skip {repr(assembly)}", SkipAssemblyWarning)
                continue
            if len(assembly.interface) == 0:
                continue
            if chain_pairs is not None:
                tr = await assembly.get_assemble_eec_as_df()
                tr_info = tr.groupby('struct_asym_id').struct_asym_id_in_assembly.apply(frozenset).to_dict()
                tr_info = defaultdict(frozenset, tr_info)
                '''
                for a, b in chain_pairs:
                    if a == b:
                        [frozenset(i) for i in combinations(tr_info[a])]
                    else:
                        [frozenset(i) for i in product(tr_info[a], tr_info[b])]
                '''
                cur_chain_pairs = frozenset.union(*(frozenset(frozenset(i) for i in combinations(tr_info[a], 2)) if a == b else frozenset(frozenset(i) for i in product(tr_info[a], tr_info[b])) for a, b in chain_pairs))
            else:
                cur_chain_pairs = chain_pairs
            for interface in assembly.interface.values():
                if ((cur_chain_pairs is None) or (interface.chain_set in cur_chain_pairs)) and interface.info['css'] > css_cutoff and len(interface.chain_set) >= chains_count_cutoff:
                    if func == 'pipe_protein_protein_interface' or func == 'pipe_protein_nucleotide_interface':
                        if not allow_not_polymer and not (interface.info['chains'][0][1] is None and interface.info['chains'][1][1] is None):
                            continue
                    elif func == 'pipe_protein_ligand_interface':
                        if not allow_not_polymer and not (interface.info['chains'][0][1] is None or interface.info['chains'][1][1] is None):
                            continue
                    res.append(interface)
        return res
    
    @unsync
    async def get_binding_sites(self):
        br_df = await self.to_dataframe_with_kwargs(
            self.fetch_from_pdbe_api('api/pdb/entry/binding_sites/'),
            usecols=lambda x: not x.startswith('author'))
        res_df = await self.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe)
        if br_df is None:
            return
        else:
            '''
            try:
                br_df = br_df.merge(res_df)
                br_df[['details.chem_comp_id',
                    'details.chain_id',
                    'details.author_residue_number',
                    'details.author_insertion_code']] = br_df.details.apply(
                        lambda x: self.binding_details_pat.fullmatch(x).groups()).apply(Series)
                br_df['details.author_residue_number'] = br_df['details.author_residue_number'].astype(int)
                br_df['details.author_insertion_code'] = br_df['details.author_insertion_code'].fillna('')
            except Exception:
                # raise AssertionError(f"{repr(self)}: Unexpected Case for get_binding_sites")
                return br_df
            return br_df.drop(columns=['details'])
            '''
            return br_df.merge(res_df)

    @unsync
    async def get_bound_molecules(self, excluding_branched=True, molecules_col=['pdb_id', 'entity_id', 'molecule_name']):
        suffix = 'bound_excluding_branched' if excluding_branched else 'bound_molecules'
        bm_df = await self.fetch_from_pdbe_api(f'graph-api/pdb/{suffix}/', Base.to_dataframe)
        if bm_df is None:
            return
        if isinstance(bm_df.composition.iloc[0], str):
            bm_df.composition = bm_df.composition.apply(lambda x: json.loads(x)['ligands'])
        bm_df = split_df_by_chain(bm_df, bm_df.columns, ['composition'], 'list')
        bm_df = concat([bm_df, bm_df.composition.apply(Series)], axis=1).drop(columns=['composition']).rename(columns={'entity': 'entity_id'})
        bm_df[['chain_id', 'chain_id_suffix']] = bm_df.chain_id.str.split('_').apply(Series)
        bm_df.author_insertion_code = bm_df.author_insertion_code.str.strip()
        bm_df.entity_id = bm_df.entity_id.astype(int)
        bm_df.molecule_type = bm_df.molecule_type.apply(lambda x: self.fix_molecule_type_dict[x])
        res_df = await self.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe)
        mol_df = (await self.fetch_from_pdbe_api('api/pdb/entry/molecules/', Base.to_dataframe))[molecules_col]
        new_bm_df = bm_df.merge(res_df).merge(mol_df)
        assert new_bm_df.shape[0] == bm_df.shape[0]
        return new_bm_df
    
    @staticmethod
    def get_fixed_site_for_bm(info):
        if isinstance(info ,str):
            info = json.loads(info)
        rets = info['chain_id'].split('_')
        return rets[0], rets[1], info['chem_comp_id'], int(info['author_residue_number']), info['author_insertion_code'].strip()

    @unsync
    async def get_bound_molecule_interaction(self, bm_id, molecules_col=['pdb_id', 'entity_id', 'molecule_type', 'molecule_name']):
        br_df = await self.fetch_from_pdbe_api('graph-api/pdb/bound_molecule_interactions/', Base.to_dataframe, mask_id=f"{self.pdb_id}/{bm_id}")
        if br_df is None:
            return
        bind_res_df = br_df.end.apply(self.get_fixed_site_for_bm).apply(Series)
        bind_res_df.columns = ['chain_id', 'chain_id_suffix', 'residue_name', 'author_residue_number', 'author_insertion_code']
        bind_res_df['bm_id'] = br_df.bm_id
        res_df = await self.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe)
        mol_df = (await self.fetch_from_pdbe_api('api/pdb/entry/molecules/', Base.to_dataframe))[molecules_col]
        return bind_res_df.merge(res_df.merge(mol_df))  # except waters

    @unsync
    async def rcsb_cluster_membership(self, entity_id, identity_cutoff:int=100):
        assert identity_cutoff in (100, 95, 90, 70, 50, 30)
        res = await self.rcsb_sqlite_api.rcsb_cluster_membership.objects.filter(pdb_id=self.pdb_id, entity_id=entity_id, identity=identity_cutoff).all()
        if res:
            assert len(res) == 1, f"{repr(self)} {(entity_id, identity_cutoff)}: Unexpected number of cluster hits!"
            return DataFrame(await self.rcsb_sqlite_api.rcsb_cluster_membership.objects.filter(cluster_id=res[0].cluster_id, identity=res[0].identity).all())
        query = '''
        {
            polymer_entity(entry_id: "%s", entity_id:"%s") {
                rcsb_cluster_membership {
                    cluster_id
                    identity
                }
            }
        }
        ''' % (self.pdb_id.upper(), entity_id)
        res = await self.fetch_from_rcsb_api('graphql', query=query, json=True, then_func=a_load_json)
        if res is None:
            return
        search_template = '''
        {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_cluster_membership.cluster_id",
                            "operator": "equals",
                            "value": %s
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_cluster_membership.identity",
                            "operator": "equals",
                            "value": %s
                        }
                    }
                ]
            },
            "request_options": {
                "return_all_hits": true
            },
            "return_type": "polymer_entity"
        }
        '''
        dfs = []
        try:
            assert res['data']['polymer_entity']['rcsb_cluster_membership'] is not None
        except Exception as e:
            info = f"polymer_entity(entry_id: \"{self.pdb_id}\", entity_id: \"{entity_id}\") -> {res}"
            if isinstance(e, AssertionError):
                warn(info, WithoutRCSBClusterMembershipWarning)
                return
            else:
                raise ValueError(info)
        for i in res['data']['polymer_entity']['rcsb_cluster_membership']:
            if i['identity'] != identity_cutoff:
                continue
            cur = await self.fetch_from_rcsb_api('search', query=search_template % (i['cluster_id'], i['identity']), then_func=Base.result_set_to_dataframe)
            cur['cluster_id'] = i['cluster_id']
            cur['identity'] = i['identity']
            dfs.append(cur.drop(columns=['services']))
        if not dfs:
            return
        df = concat(dfs, ignore_index=True, sort=False)
        df[['pdb_id', 'entity_id']] = df.identifier.str.split('_').apply(Series)
        df.rename(columns={'identifier': 'rcsb_id'}, inplace=True)
        df.pdb_id = df.pdb_id.str.lower()
        df.entity_id = df.entity_id.astype(int)
        await self.rcsb_sqlite_api.async_insert(self.rcsb_sqlite_api.rcsb_cluster_membership, df.to_dict('records'))
        return df


class PDBAssembly(PDB):

    tasks = LRUCache(maxsize=1024)

    struct_range_pattern = re_compile(r"\[.+\]([A-Z]+[_0-9]*):(-?[0-9]+[\?A-Z]*)")  # e.g. [FMN]B:149 [C2E]A:301 [ACE]H:-8? [BR]BA:957A
    #rare_pat = re_compile(r"([A-Z]+)_([0-9]+)")  # e.g. 2rde assembly 1 A_1, B_1...
    interface_structures_pat = re_compile(r"(\[.+\])?([A-Z]+):?(-?[0-9]+[\?A-Z]*)?\+(\[.+\])?([A-Z]+):?(-?[0-9]+[\?A-Z]*)?")  # [4CA]BB:170+AB [ZN]D:154A+[CU]C:154 [ACE]B:0?+B	
    seq_pat = re_compile(r"(-?[0-9]+)([\?A-Z]*)?")

    @property
    def assembly_summary(self) -> Dict:
        for ass in self.summary['assemblies']:
            if int(ass['assembly_id']) == self.assembly_id:
                return ass
        raise ValueError(f"{repr(self)}: Without expected assemble info\n{self.summary['assemblies']}")

    def add_args(self, pdb_ob=None):
        if pdb_ob is None:
            self.pdb_ob = PDB(self.pdb_id)
        else:
            self.pdb_ob = pdb_ob
        return self

    def __post_init__(self):
        super().__post_init__()
        '''
        NOTE: reference: <https://www.ebi.ac.uk/training/online/course/pdbepisa-identifying-and-interpreting-likely-biolo/1555-special-code-doing-nothing-structure>
        '''
        self.interface_filters = {
            'symmetry_operator': ('isin', ('1_555', '1555', 1555))  # 1555 for api%pisa%asiscomponent%+6e4h%0%interfaces
        }  # 'structure_2.symmetry_id': ('eq', '1_555'),'css': ('ge', 0)

    def set_id(self):
        assert self.level == 'entry_like' and self.raw_identifier.count('/') == 1
        self.pdb_ass_id = self.raw_identifier.lower()
        self.pdb_id, self.assembly_id = self.pdb_ass_id.split('/')
        self.assembly_id = int(self.assembly_id)
    
    def get_id(self):
        return self.pdb_ass_id

    """@classmethod
    def transform(cls, x):
        res = cls.rare_pat.search(x)
        assert bool(res), f"Unexpected Case: {x}"
        chain, num = res.groups()
        num = int(num)
        if num == 1:
            return chain
        else:
            return chain+chr(63+num)"""

    @classmethod
    def fix_auth_seq_id(cls, seq_id):
        if seq_id is None:
            return None,None
        else:
            auth_res_num, insertion_code = cls.seq_pat.fullmatch(seq_id).groups()
            return auth_res_num, '' if insertion_code == '?' else insertion_code

    @classmethod
    async def to_asiscomponent_interfaces_df(cls, path: Unfuture):
        interfacelist_df = await path.then(cls.to_dataframe)
        if interfacelist_df is None:
            return None
        interfacelist_df.rename(columns={"interface_number": "interface_id"}, inplace=True)
        try:
            interfacelist_df[[
                'struct_asym_id_in_assembly_1', 'auth_seq_id_1',
                'struct_asym_id_in_assembly_2', 'auth_seq_id_2']
                ] = interfacelist_df.interface_structures.apply(
                            lambda x: cls.interface_structures_pat.fullmatch(x).group(2,3,5,6)).apply(Series)
            interfacelist_df[['author_residue_number_1', 'author_insertion_code_1']] = interfacelist_df.auth_seq_id_1.apply(cls.fix_auth_seq_id).apply(Series)
            interfacelist_df[['author_residue_number_2', 'author_insertion_code_2']] = interfacelist_df.auth_seq_id_2.apply(cls.fix_auth_seq_id).apply(Series)
        except AttributeError as e:
            check = interfacelist_df.interface_structures.apply(lambda x: bool(cls.interface_structures_pat.fullmatch(x)))
            warn(str(interfacelist_df[~check]))
            raise e
        return interfacelist_df.drop(columns=['auth_seq_id_1', 'auth_seq_id_2'])

    @classmethod
    async def to_interfacelist_df(cls, path: Unfuture):
        interfacelist_df = await path.then(cls.to_dataframe)
        if interfacelist_df is None:
            return None
        interfacelist_df.rename(columns={"id": "interface_id"}, inplace=True)
        assert ((interfacelist_df['structure_1.range'].notnull()) & (interfacelist_df['structure_2.range'].notnull())).all(), str(interfacelist_df[interfacelist_df['structure_1.range'].isnull() | interfacelist_df['structure_2.range'].isnull()])
        for col in ('struct_asym_id_in_assembly_1', 'auth_seq_id_1', 'struct_asym_id_in_assembly_2', 'auth_seq_id_2'):
            interfacelist_df[col] = None
        if (interfacelist_df['structure_1.range'].str.contains('_') | interfacelist_df['structure_2.range'].str.contains('_')).any():
            raise AssertionError(f'Error structure.range in \n{interfacelist_df}')
        else:
            range_check = interfacelist_df['structure_1.range'].apply(cls.struct_range_pattern.fullmatch)
            range_check_yes = range_check.apply(bool)
            check_m = range_check_yes | (interfacelist_df['structure_1.original_range']=='{-}')
            assert check_m.all(), f"{interfacelist_df[~check_m].T}"
            if range_check_yes.any():
                to_df_index = interfacelist_df[range_check_yes].index
                #try:
                interfacelist_df.loc[to_df_index, ['struct_asym_id_in_assembly_1', 'auth_seq_id_1']] = range_check.loc[to_df_index].apply(lambda x: x.groups()).tolist()
                #except ValueError:
                #    raise ValueError("%s\n%s" % (range_check.loc[to_df_index].apply(lambda x: x.groups()), interfacelist_df.loc[to_df_index, 'structure_1.range']))
            range_check_no = ~range_check_yes
            if range_check_no.any():
                no_index = interfacelist_df[range_check_no].index
                interfacelist_df.loc[no_index, 'struct_asym_id_in_assembly_1'] = interfacelist_df.loc[no_index, 'structure_1.range']
            #
            range_check = interfacelist_df['structure_2.range'].apply(cls.struct_range_pattern.fullmatch)
            range_check_yes = range_check.apply(bool)
            check_m = range_check_yes | (interfacelist_df['structure_2.original_range']=='{-}')
            assert check_m.all(), f"{interfacelist_df[~check_m].T}"
            if range_check_yes.any():
                to_df_index = interfacelist_df[range_check_yes].index
                interfacelist_df.loc[to_df_index, ['struct_asym_id_in_assembly_2', 'auth_seq_id_2']] = range_check.loc[to_df_index].apply(lambda x: x.groups()).tolist()
            range_check_no = ~range_check_yes
            if range_check_no.any():
                no_index = interfacelist_df[range_check_no].index
                interfacelist_df.loc[no_index, 'struct_asym_id_in_assembly_2'] = interfacelist_df.loc[no_index, 'structure_2.range']
        interfacelist_df[['author_residue_number_1', 'author_insertion_code_1']] = interfacelist_df.auth_seq_id_1.apply(cls.fix_auth_seq_id).apply(Series)
        interfacelist_df[['author_residue_number_2', 'author_insertion_code_2']] = interfacelist_df.auth_seq_id_2.apply(cls.fix_auth_seq_id).apply(Series)
        return interfacelist_df.drop(columns=['auth_seq_id_1', 'auth_seq_id_2'])

    @staticmethod
    def extract_partner_detail_pipeline(partner_detail):
        segment_feat_dict_0 = {}
        segment_feat_dict_1 = {}
        for segment_feat in partner_detail['interface_partner_feature'][0]['feature_positions']:
            segment_feat_dict_0.update(dict(zip(range(segment_feat['beg_seq_id'], segment_feat['end_seq_id']+1), segment_feat['values'])))
        for segment_feat in partner_detail['interface_partner_feature'][1]['feature_positions']:
            segment_feat_dict_1.update(dict(zip(range(segment_feat['beg_seq_id'], segment_feat['end_seq_id']+1), segment_feat['values'])))
        assert segment_feat_dict_0.keys() == segment_feat_dict_1.keys()
        partner_interface_residue_numbers = [key for key,val in segment_feat_dict_0.items() if abs(val - segment_feat_dict_1[key]) > 0]
        return partner_detail['interface_partner_identifier']['asym_id'], partner_detail['interface_partner_identifier']['entity_id'], partner_interface_residue_numbers

    @unsync
    async def get_interface_info_from_rcsb_data_api(self):
        # assemblies(assembly_ids: ["4HHB-1", "1A01-1", "3PQR-1"]) {}
        # interface(entry_id: "1A01", assembly_id: "1", interface_id: "1"){}
        
        query = '''
        {
            assembly(entry_id: "%s", assembly_id: "%s") {
                interfaces {
                    rcsb_interface_operator
                    rcsb_interface_info {
                        interface_area
                        interface_character
                        num_core_interface_residues
                        num_interface_residues
                        polymer_composition
                        self_jaccard_contact_score
                    }
                    rcsb_interface_container_identifiers {
                        interface_id
                    }
                    rcsb_interface_partner {
                        interface_partner_identifier {
                            asym_id
                            entity_id
                        }
                        interface_partner_feature {
                            type
                            feature_positions {
                                beg_seq_id
                                end_seq_id
                                values
                            }
                        }
                    }
                }
            }
        }
        ''' % (self.pdb_id, self.assembly_id)
        interfaces_data = (await self.fetch_from_rcsb_api('graphql', query=query, then_func=a_load_json, json=True))['data']['assembly']
        if interfaces_data is None or interfaces_data['interfaces'] is None:
            return
        else:
            interfaces_data = interfaces_data['interfaces']
        records = []
        
        profile_df = (await PDB(self.pdb_id).profile_id()).query(f'assembly_id == {self.assembly_id}')
        
        for interface_data in interfaces_data:
            meta = interface_data['rcsb_interface_info']

            interface_id = interface_data['rcsb_interface_container_identifiers']['interface_id']
            partner_1, partner_2 = tuple(self.extract_partner_detail_pipeline(partner_detail) for partner_detail in interface_data['rcsb_interface_partner'])

            operator = json.dumps(interface_data["rcsb_interface_operator"]).decode("utf-8")
            
            operator_arr = array(interface_data["rcsb_interface_operator"]).transpose(1,0,2)
            for chain_1_operator, chain_2_operator in operator_arr:
                chain_1_operator, chain_2_operator = json.dumps(chain_1_operator.tolist()).decode("utf-8"), json.dumps(chain_2_operator.tolist()).decode("utf-8")
                chain_1_record = profile_df[profile_df.oper_expression.eq(chain_1_operator) & profile_df.struct_asym_id.eq(partner_1[0])]
                assert chain_1_record.shape[0] == 1
                chain_1_record = chain_1_record.iloc[0]
                chain_2_record = profile_df[profile_df.oper_expression.eq(chain_2_operator) & profile_df.struct_asym_id.eq(partner_2[0])]
                assert chain_2_record.shape[0] == 1
                chain_2_record = chain_2_record.iloc[0]
                break # just get the first pair

            record = dict(pdb_id=self.pdb_id, assembly_id=self.assembly_id, interface_id=interface_id,
                          struct_asym_id_1=partner_1[0], entity_id_1=int(partner_1[1]), interface_range_1=to_interval(partner_1[2]), model_id_1=chain_1_record.model_id, asym_id_rank_1=chain_1_record.asym_id_rank, struct_asym_id_in_assembly_1=chain_1_record.struct_asym_id_in_assembly, oper_expression_1=chain_1_operator,
                          struct_asym_id_2=partner_2[0], entity_id_2=int(partner_2[1]), interface_range_2=to_interval(partner_2[2]), model_id_2=chain_2_record.model_id, asym_id_rank_2=chain_2_record.asym_id_rank, struct_asym_id_in_assembly_2=chain_2_record.struct_asym_id_in_assembly, oper_expression_2=chain_2_operator,
                          operator=operator)
            record.update(meta)
            records.append(record)
        return DataFrame(records)

    @unsync
    async def get_interfacelist_df(self, api_suffix, func):
        use_au = self.assembly_id in self.pdb_ob.subset_assembly
        if api_suffix == 'api/pisa/asiscomponent/':
            temp1 = '{}/0/interfaces'
            temp2 = '{}/interfaces'
        elif api_suffix == 'api/pisa/interfacelist/':
            temp1 = '{}/0'
            temp2 = '{}'
        else:
            raise ValueError(f"Invalid suffix: {api_suffix}")
        try:
            interfacelist_df = await self.fetch_from_pdbe_api(api_suffix, func, mask_id=temp1.format(self.pdb_id) if use_au else temp2.format(self.get_id()))
        except WithoutExpectedKeyError:
            try:
                interfacelist_df = await self.fetch_from_pdbe_api(api_suffix, func, mask_id=temp2.format(self.get_id()))
                use_au = False
            except WithoutExpectedKeyError:
                return None, use_au
        return interfacelist_df, use_au

    @unsync
    async def set_interface(self, obligated_class_chains: Optional[Iterable[str]] = None, allow_same_class_interaction:bool=True):

        def to_interface_id(pdb_assembly_id, focus_interface_ids):
            for interface_id in focus_interface_ids:
                yield f"{pdb_assembly_id}/{interface_id}"
        
        interfacelist_df, use_au = await self.get_interfacelist_df(
            'api/pisa/asiscomponent/', PDBAssembly.to_asiscomponent_interfaces_df)
        if interfacelist_df is None:
            interfacelist_df, use_au = await self.get_interfacelist_df(
                'api/pisa/interfacelist/', PDBAssembly.to_interfacelist_df)
            self.interface_filters['structure_2.symmetry_id'] = ('isin', ('1_555', '1555', 1555))
            del self.interface_filters['symmetry_operator']
        else:
            interfacelist_df = interfacelist_df.rename(columns={'complex_formation_score': 'css'})
        
        if interfacelist_df is None:
            warn(f"cannot get interfacelist data from PISA: {self.get_id()}", PISAErrorWarning)
            self.interface = dict()
            return
        else:
            interfacelist_df = interfacelist_df.copy()
            interfacelist_df['assembly_id'] = self.assembly_id
        if use_au:
            a_chains = self.pdb_ob.subset_assembly[self.assembly_id]
            interfacelist_df = interfacelist_df[
                interfacelist_df.struct_asym_id_in_assembly_1.isin(a_chains) &
                interfacelist_df.struct_asym_id_in_assembly_2.isin(a_chains)]

        if obligated_class_chains is not None:
            bool_1 = interfacelist_df.struct_asym_id_in_assembly_1.isin(obligated_class_chains)
            bool_2 = interfacelist_df.struct_asym_id_in_assembly_2.isin(obligated_class_chains)
            interfacelist_df = interfacelist_df[bool_1 | bool_2]
            if not allow_same_class_interaction:
                bool_1 = interfacelist_df.struct_asym_id_in_assembly_1.isin(obligated_class_chains)
                bool_2 = interfacelist_df.struct_asym_id_in_assembly_2.isin(obligated_class_chains)
                interfacelist_df = interfacelist_df[~(bool_1 & bool_2)]
        
        focus_interface_df = related_dataframe(
            self.interface_filters, interfacelist_df)
        focus_interface_ids = focus_interface_df.interface_id
        focus_interface_chains = zip(
            focus_interface_df.struct_asym_id_in_assembly_1,
            focus_interface_df.author_residue_number_1,
            focus_interface_df.author_insertion_code_1,
            focus_interface_df.struct_asym_id_in_assembly_2,
            focus_interface_df.author_residue_number_2,
            focus_interface_df.author_insertion_code_2)

        self.interface: Dict[int, PDBInterface] = dict(zip(
            focus_interface_ids, (PDBInterface(if_id).add_args(PDBAssembly_ob=self, use_au=use_au).store(chains=(chains[:3], chains[3:]), css=css) for if_id, chains, css in zip(to_interface_id(self.get_id(), focus_interface_ids), focus_interface_chains, focus_interface_df.css))))

    def get_interface(self, interface_id):
        return self.interface[interface_id]

    @unsync
    async def set_assemble_eec_as_df(self):
        eec_as_df = await self.pdb_ob.profile_id()
        if self.assembly_id not in eec_as_df.assembly_id:
            raise PossibleInvalidAssemblyError(f"{repr(self)}: Invalid assembly_id!")
        self.assemble_eec_as_df = eec_as_df[eec_as_df.assembly_id.eq(self.assembly_id)]
    
    @unsync
    async def get_assemble_eec_as_df(self):
        if hasattr(self, 'assemble_eec_as_df'):
            return self.assemble_eec_as_df
        else:
            await self.set_assemble_eec_as_df()
            return self.assemble_eec_as_df

    @unsync
    async def pipe_protein_protein_interface(self):
        '''
        DEMO

        >>> eec_as_df = self.pdb_ob.get_eec_as_df().result()
        >>> protein_type_asym = eec_as_df[
            eec_as_df.assembly_id.eq(self.assembly_id) &
            eec_as_df.molecule_type.eq('polypeptide(L)')].struct_asym_id_in_assembly
        >>> self.interface_filters['struct_asym_id_in_assembly_1'] = ('isin', protein_type_asym)
        >>> self.interface_filters['struct_asym_id_in_assembly_2'] = ('isin', protein_type_asym)
        >>> self.set_interface().result()
        '''
        
        eec_as_df = await self.get_assemble_eec_as_df()
        eec_as_df = eec_as_df[eec_as_df.molecule_type.isin(('polypeptide(L)', 'polypeptide(D)'))]
        if len(eec_as_df) == 1:
            self.interface = dict()
            return
        protein_type_asym = eec_as_df.struct_asym_id_in_assembly
        self.interface_filters['struct_asym_id_in_assembly_1'] = ('isin', protein_type_asym)
        self.interface_filters['struct_asym_id_in_assembly_2'] = ('isin', protein_type_asym)
        await self.set_interface()

    @unsync
    async def pipe_protein_else_interface(self, molecule_types, allow_same_class_interaction=False):
        eec_as_df = await self.get_assemble_eec_as_df()
        eec_as_df = eec_as_df[eec_as_df.molecule_type.isin(frozenset({'polypeptide(L)', 'polypeptide(D)'}) | frozenset(molecule_types))]
        if len(eec_as_df) == 1:
            self.interface = dict()
            return
        molType_dict = eec_as_df.groupby('molecule_type').struct_asym_id_in_assembly.apply(frozenset).to_dict()

        protein_type_asym = molType_dict.get('polypeptide(L)', frozenset()) | molType_dict.get('polypeptide(D)', frozenset())
        else_type_asym = frozenset.union(*(molType_dict.get(molecule_type, frozenset()) for molecule_type in molecule_types))
        target_type_asym = protein_type_asym | else_type_asym
        
        self.interface_filters['struct_asym_id_in_assembly_1'] = ('isin', target_type_asym)
        self.interface_filters['struct_asym_id_in_assembly_2'] = ('isin', target_type_asym)
        await self.set_interface(obligated_class_chains=protein_type_asym, allow_same_class_interaction=allow_same_class_interaction)

    @unsync
    async def pipe_protein_ligand_interface(self):
        '''
        DEMO

        >>> eec_as_df = self.pdb_ob.get_eec_as_df().result()
        >>> molType_dict = eec_as_df[eec_as_df.assembly_id.eq(self.assembly_id)].groupby('molecule_type').struct_asym_id_in_assembly.apply(set).to_dict()

        >>> protein_type_asym = molType_dict.get('polypeptide(L)', set())
        >>> ligand_type_asym = molType_dict.get('bound', set())
        >>> target_type_asym = protein_type_asym | ligand_type_asym

        >>> self.interface_filters['struct_asym_id_in_assembly_1'] = ('isin', target_type_asym)
        >>> self.interface_filters['struct_asym_id_in_assembly_2'] = ('isin', target_type_asym)
        >>> self.set_interface(obligated_class_chains=protein_type_asym).result()
        '''
        await self.pipe_protein_else_interface(('bound',), False)

    @unsync
    async def pipe_protein_nucleotide_interface(self, molecule_types=None):
        await self.pipe_protein_else_interface((
            'polydeoxyribonucleotide', 
            'polyribonucleotide', 
            'polydeoxyribonucleotide/polyribonucleotide hybrid') if molecule_types is None else molecule_types, False)


class PDBInterface(PDBAssembly):
 
    tasks = LRUCache(maxsize=1024)

    def add_args(self, PDBAssembly_ob: Optional[PDBAssembly]=None, use_au:bool=False):
        self.use_au = use_au
        if PDBAssembly_ob is None:
            self.PDBAssembly_ob = PDBAssembly(
                f"{self.pdb_id}/{self.assembly_id}").add_args()
        else:
            self.PDBAssembly_ob = PDBAssembly_ob
        return self

    def __repr__(self):
        if hasattr(self, 'info'):
            return "<{} {} {}>".format(self.__class__.__name__, self.get_id(), self.info).replace('frozenset', '')
        else:
            return f"<{self.__class__.__name__} {self.get_id()}>"

    def set_id(self):
        assert self.level == 'entry_like' and self.raw_identifier.count('/') == 2
        self.pdb_ass_int_id = self.raw_identifier.lower()
        self.pdb_id, self.assembly_id, self.interface_id = self.raw_identifier.split('/')
        self.assembly_id = int(self.assembly_id)
        self.interface_id = int(self.interface_id)
    
    def get_id(self):
        return self.pdb_ass_int_id

    def store(self, **kwargs):
        self.info = kwargs
        self.chain_set = frozenset(i[0] for i in self.info['chains'])
        return self

    @classmethod
    async def to_interfacedetail_df(cls, path: Unfuture):

        """def check_struct_selection(interfacedetail_df, colName):
            sele = next(iter(interfacedetail_df[colName]))
            sele_m = cls.struct_range_pattern.fullmatch(sele)
            if bool(sele_m):
                interfacedetail_df[colName] = sele_m.group(1)"""

        interfacedetail_df = await cls.to_dataframe_with_kwargs(path,
            #usecols=['pdb_code', 'assemble_code', 'interface_number', 'chain_id',
            #         'residue', 'sequence', 'insertion_code',
            #         'buried_surface_area', 'solvent_accessible_area',
            #         'interface_detail.interface_structure_1.structure.selection',
            #         'interface_detail.interface_structure_2.structure.selection'],
            na_values=[' ', '?'])
        if interfacedetail_df is None:
            return None
        interfacedetail_df.rename(
            columns={"pdb_code": "pdb_id",
                     "sequence": "author_residue_number",
                     "insertion_code": "author_insertion_code",
                     "residue": "residue_name",
                     "chain_id": "struct_asym_id_in_assembly",
                     "assemble_code": "assembly_id",
                     "interface_number": "interface_id",
                     'interface_detail.interface_structure_1.structure.selection': "s1_selection",
                     'interface_detail.interface_structure_2.structure.selection': "s2_selection"},
            inplace=True)
        interfacedetail_df.author_insertion_code.fillna('', inplace=True)
        return interfacedetail_df

    @staticmethod
    def molecule_type_score(m_type):
        if m_type == 'polypeptide(L)':
            return 1
        elif m_type == 'polypeptide(D)':
            return 2
        else:
            return 3

    @unsync
    async def set_interface_res(self, keep_interface_res_df:bool=False):
        # [ ]X-X_mod interaction
        # [ ]X-Y_mod interaction
        # [ ]X-X interaction -> error?
        # [-]X-Y interaction
        run_following = False
        task = self.tasks.get((repr(self), 'PISAInterfaceDict'), None)
        if task is None:
            res = await self.sqlite_api.PISAInterfaceDict.objects.filter(pdb_id=self.pdb_id, assembly_id=self.assembly_id, interface_id=self.interface_id).all()
            if len(res) == 1:
                self.interface_res_dict = res[0]
                self.register_task((repr(self), 'PISAInterfaceDict'), self.interface_res_dict)
            elif len(res) > 1:
                raise AssertionError(f"{repr(self)}: Unexpected length of {res}")
            else:
                run_following = True
        else:
            self.interface_res_dict = task
        if not run_following and not keep_interface_res_df:
            return
        try:
            interfacedetail_df = await self.fetch_from_pdbe_api('api/pisa/interfacedetail/', PDBInterface.to_interfacedetail_df, mask_id=f'{self.pdb_id}/0/{self.interface_id}' if self.use_au else None)
            if interfacedetail_df is None:
                raise WithoutExpectedContentError('')
            else:
                interfacedetail_df = interfacedetail_df.copy()
        except (WithoutExpectedKeyError, WithoutExpectedContentError, RetryError) as e:
            warn(f"cannot get interfacedetail data from PISA: {repr(self)}, {e.__class__.__name__}", PISAErrorWarning)
            return
        except Exception as e:
            raise AssertionError(e)
        interfacedetail_df['assembly_id'] = self.assembly_id

        assert not interfacedetail_df.struct_asym_id_in_assembly.str.contains('_').any(), f"Error chain_id in {repr(self)}"
        s1_selection = interfacedetail_df.loc[0, 's1_selection']
        s2_selection = interfacedetail_df.loc[0, 's2_selection']
        s1_match = self.struct_range_pattern.fullmatch(s1_selection)
        if bool(s1_match):
            struct_asym_id_in_assembly_1, auth_seq_id_1 = s1_match.groups()
            author_residue_number_1, author_insertion_code_1 = self.fix_auth_seq_id(auth_seq_id_1)
        else:
            struct_asym_id_in_assembly_1 = s1_selection
            assert '_' not in struct_asym_id_in_assembly_1, f"Error s1_selection in {repr(self)}"
            author_residue_number_1, author_insertion_code_1 = None, None
        s2_match = self.struct_range_pattern.fullmatch(s2_selection)
        if bool(s2_match):
            struct_asym_id_in_assembly_2, auth_seq_id_2 = s2_match.groups()
            author_residue_number_2, author_insertion_code_2 = self.fix_auth_seq_id(auth_seq_id_2)
        else:
            struct_asym_id_in_assembly_2 = s2_selection
            assert '_' not in struct_asym_id_in_assembly_2, f"Error s2_selection in {repr(self)}"
            author_residue_number_2, author_insertion_code_2 = None, None
        map_dict = {s1_selection: struct_asym_id_in_assembly_1, s2_selection: struct_asym_id_in_assembly_2}
        interfacedetail_df.struct_asym_id_in_assembly = interfacedetail_df.struct_asym_id_in_assembly.map(map_dict)
        #if len(struct_sele_set) != 2:
            # NOTE: Exception example: 2beq assembly_id 1 interface_id 32
            #warn(f"skip possible modified_residue <-> standard_residues interaction within the same chain:\n{interfacedetail_df.head(1)[['pdb_id', 'assembly_id', 'interface_id', 's1_selection']].to_dict('records')[0]}")
            #return
        if not hasattr(self, 'info'):
            pass
        elif self.info['chains'] != ((struct_asym_id_in_assembly_1, author_residue_number_1, author_insertion_code_1), (struct_asym_id_in_assembly_2, author_residue_number_2, author_insertion_code_2)):
            # NOTE: Exception example: 2beq assembly_id 1 interface_id 32
            warn(f"{repr(self)}: interfacedetail({((struct_asym_id_in_assembly_1, author_residue_number_1, author_insertion_code_1), (struct_asym_id_in_assembly_2, author_residue_number_2, author_insertion_code_2))}) inconsistent with interfacelist({self.info['chains']}) ! May miss some data.", PISAErrorWarning)
            return
        eec_as_df = await self.PDBAssembly_ob.get_assemble_eec_as_df()
        res_df = await self.PDBAssembly_ob.pdb_ob.fetch_from_pdbe_api('api/pdb/entry/residue_listing/', Base.to_dataframe)
        interfacedetail_df = interfacedetail_df.merge(eec_as_df, how="left")
        interfacedetail_df = interfacedetail_df.merge(res_df, how="left")
        if author_residue_number_1 is not None:
            try:
                struct_asym_id_1 = eec_as_df[eec_as_df.struct_asym_id_in_assembly.eq(struct_asym_id_in_assembly_1)].struct_asym_id.iloc[0]
                residue_number_1 = res_df[res_df.struct_asym_id.eq(struct_asym_id_1) &
                    res_df.author_residue_number.eq(int(author_residue_number_1)) &
                    res_df.author_insertion_code.eq(author_insertion_code_1)
                ].residue_number.iloc[0]
            except IndexError:
                raise PISAOutOfDateError(f"{repr(self)}: {struct_asym_id_in_assembly_1}, {author_residue_number_1}, {author_insertion_code_1}")
                #raise e
        else:
            residue_number_1 = None
        if author_residue_number_2 is not None:
            try:
                struct_asym_id_2 = eec_as_df[eec_as_df.struct_asym_id_in_assembly.eq(struct_asym_id_in_assembly_2)].struct_asym_id.iloc[0]
                residue_number_2 = res_df[res_df.struct_asym_id.eq(struct_asym_id_2) &
                    res_df.author_residue_number.eq(int(author_residue_number_2)) &
                    res_df.author_insertion_code.eq(author_insertion_code_2)
                ].residue_number.iloc[0]
            except IndexError:
                raise PISAOutOfDateError(f"{repr(self)}: {struct_asym_id_in_assembly_2}, {author_residue_number_2}, {author_insertion_code_2}")
                #raise e
        else:
            residue_number_2 = None
        if keep_interface_res_df:
            self.interface_res_df = interfacedetail_df
            if hasattr(self, 'interface_res_dict'):
                return
        # NOTE: Not rigorous filter for multiple_conformers
        check_merge = interfacedetail_df.residue_number.isnull()
        check_mc = interfacedetail_df.author_residue_number < 0
        if check_merge.any():
            # raise ValueError(f"Unexpected Data in Residue DataFrame: {check.head(1).to_dict('records')[0]}")
            example = interfacedetail_df[check_merge & check_mc].head(3)
            if example.shape[0]:
                warn(f"Unexpected Data in Residue DataFrame: \n{example[['pdb_id', 'assembly_id', 'interface_id', 'struct_asym_id_in_assembly', 'residue_name','author_residue_number']].to_dict('index')}", PISAErrorWarning)
            to_drop = interfacedetail_df[check_merge].index
            if len(interfacedetail_df) == len(to_drop):
                # ['assembly_id', 'pdb_id', 'struct_asym_id_in_assembly']
                # ['author_insertion_code', 'author_residue_number', 'chain_id', 'entity_id', 'pdb_id', 'residue_name', 'struct_asym_id']
                raise ValueError(f"{repr(self)}: Unexpected null in Residue DataFrame\n{interfacedetail_df.head().T}")
            interfacedetail_df = interfacedetail_df.drop(index=to_drop)
        
        focus_cols = ['pdb_id', 'entity_id', 'chain_id', 'struct_asym_id', 
                      'struct_asym_id_in_assembly', 'asym_id_rank', 'model_id',
                      'assembly_id', 'interface_id', 'molecule_type',
                      'residue_number', 'buried_surface_area', 'solvent_accessible_area']
        asa_col = 12 # focus_cols.index('solvent_accessible_area')
        bsa_col = 11 # focus_cols.index('buried_surface_area')
        res_col = 10 # focus_cols.index('residue_number')
        nda = interfacedetail_df[focus_cols].to_numpy()

        def yield_record():
            for _, (start, end) in slice_series(interfacedetail_df.struct_asym_id_in_assembly).items():
                asa_index = start+np_where(nda[start:end, asa_col] > 0)[0]
                asa_res = nda[asa_index, res_col]
                bsa_index = asa_index[np_where(nda[asa_index, bsa_col] > 0)]
                bsa_res = nda[bsa_index, res_col]
                info = nda[start, 0:len(focus_cols)-3].tolist() + \
                    [json.dumps(to_interval(asa_res)).decode('utf-8'),
                    json.dumps(to_interval(bsa_res)).decode('utf-8')]
                yield dict(zip(focus_cols[:-3]+['surface_range', 'interface_range'], info))
        
        records = sorted(yield_record(), key=lambda x: (self.molecule_type_score(x['molecule_type']), -id2score(x['struct_asym_id']), x['model_id']))

        common_keys = ('pdb_id', 'assembly_id', 'interface_id')
        record1 = {f"{key}_1": value for key,
                value in records[0].items() if key not in common_keys}
        assert 'interface_range_1' in record1, f"Error in {repr(self)}"
        record1['is_polymer_1'] = True
        if record1['struct_asym_id_in_assembly_1'] == struct_asym_id_in_assembly_1:
            the_other_sele = struct_asym_id_in_assembly_2
            the_other_interface_range = f'[[{residue_number_2},{residue_number_2}]]' if residue_number_2 is not None else None
        else:
            the_other_sele = struct_asym_id_in_assembly_1
            the_other_interface_range = f'[[{residue_number_1},{residue_number_1}]]'if residue_number_1 is not None else None
        if len(records) == 2:
            record2 = {f"{key}_2": value for key,
                    value in records[1].items() if key not in common_keys}
            assert record2['struct_asym_id_in_assembly_2'] == the_other_sele, f"Error in {repr(self)}"
            assert 'interface_range_2' in record2, f"Error in {repr(self)}"
            record2['is_polymer_2'] = True
        else:
            record2 = {'struct_asym_id_in_assembly_2': the_other_sele}
            cur_keys = list(set(focus_cols[:-3])-set(common_keys)-{'struct_asym_id_in_assembly'})
            try:
                cur_record = eec_as_df[eec_as_df.struct_asym_id_in_assembly.eq(the_other_sele)][cur_keys].to_dict('records')[0]
            except Exception:
                raise ValueError(f"\n{self.get_id()},\n{the_other_sele},\n{eec_as_df}")
            for key, value in cur_record.items():
                record2[f"{key}_2"] = value
            assert the_other_interface_range is not None, f"Error in {repr(self)}"
            record2['interface_range_2'] = the_other_interface_range
            record2['is_polymer_2'] = False

        record_dict = {**record1, **record2}
        for key in common_keys:
            record_dict[key] = records[0][key]
        record_dict['use_au'] = self.use_au
        if hasattr(self, 'info'):
            record_dict['css'] =  self.info['css']
        self.interface_res_dict = record_dict
        self.register_task((repr(self), 'PISAInterfaceDict'), record_dict)
        # self.sqlite_api.sync_insert(self.sqlite_api.PISAInterfaceDict, [record_dict])
        # await self.sqlite_api.PISAInterfaceDict.objects.create(**record_dict)

    @unsync
    async def get_interface_res_df(self):
        if hasattr(self, 'interface_res_df'):
            return self.interface_res_df
        else:
            await self.set_interface_res(True)
            if hasattr(self, 'interface_res_df'):
                return self.interface_res_df
            else:
                return
    
    @unsync
    async def get_interface_res_dict(self, **kwargs):
        if hasattr(self, 'interface_res_dict'):
            return self.interface_res_dict
        else:
            await self.set_interface_res(**kwargs)
            if hasattr(self, 'interface_res_dict'):
                return self.interface_res_dict
            else:
                return


class RCSB1DCoordinates(Base):
    '''RCSB 1D-Coordinates'''
    tasks = LRUCache(maxsize=1024)

    sequence_reference = {
        ('RefSeq', 'genome'): 'NCBI_GENOME',
        ('RefSeq', 'protein'): 'NCBI_PROTEIN',
        ('RefSeq', 'model_protein'): 'NCBI_PROTEIN',
        ('UniProt', 'isoform'): 'UNIPROT',
        ('PDB', 'entity'): 'PDB_ENTITY',
        ('PDB', 'instance'): 'PDB_INSTANCE',
    }

    def __post_init__(self):
        super().__post_init__()
        self.check_folder()
        self.seq_ref_type = self.sequence_reference[self.source, self.level]
    
    @property
    def seq_ref_id(self):
        if self.source in ('RefSeq', 'UniProt'):
            return self.identifier
        elif self.source == 'PDB':
            return self.raw_identifier.upper()
    
    def get_id(self):
        return self.seq_ref_id

    def alignment(self, seq_ref_type, with_seq:bool=False, **kwargs):
        assert seq_ref_type in self.sequence_reference.values()
        args = (self.seq_ref_type, seq_ref_type, self.seq_ref_id, 'query_sequence', 'target_sequence') if with_seq else (self.seq_ref_type, seq_ref_type, self.seq_ref_id, '', '')
        return self.fetch_from_rcsb_api(api_suffix='1d_coordinates', query='''
        {
            alignment(from:%s, to:%s, queryId:"%s"){
                %s
                target_alignment {
                    target_id
                    %s
                    coverage{
                        query_length
                        target_length
                    }
                    aligned_regions {
                        query_begin
                        query_end
                        target_begin
                        target_end
                    }
                    orientation
                }
            }
        }
        ''' % args, **kwargs)
    
    def yield_mapping(self, data):
        for mapping in data['data']['alignment']['target_alignment']:
            info = dict(query_id=self.seq_ref_id,
                        target_id=mapping['target_id'], 
                        query_length=mapping['coverage']['query_length'], 
                        target_length=mapping['coverage']['target_length'],
                        orientation=mapping['orientation'])
            for region in mapping['aligned_regions']:
                yield {**region,**info}
    
    @unsync
    async def alignment_df(self, seq_ref_type, **kwargs):
        try:
            return DataFrame(self.yield_mapping(await self.alignment(seq_ref_type, **kwargs).then(a_load_json)))
        except TypeError:
            pass


class SIFTS(PDB):
    '''
    TODO
    
    1. Better OligoState
      * RAW (both from wwPDB and self assigned)
      * FILTERED 
    ~~2. Define Best Isoform~~
    ~~3. UniProt Isoform Interaction~~
    ~~4. PDBChain Instance Interaction (Biological Relevance)~~
    '''

    tasks = LRUCache(maxsize=1024)

    EntityChain = namedtuple('EntityChain', 'pdb_id entity_chain_info entity_count chain_count')
    UniProtEntity = namedtuple('UniProtEntity', 'pdb_id unp_entity_info entity_unp_info entity_with_unp_count min_unp_count')
    OligoState = namedtuple('OligoState', 'pdb_id oligo_state has_unmapped_protein')
    MappingMeta = namedtuple('MappingMeta', 'UniProt species identity')

    chain_filter = 'UNK_COUNT < SEQRES_COUNT and ca_p_only == False and new_identity >=0.9 and repeated == False and reversed == False and OBS_STD_COUNT >= 20'
    entry_filter = '(experimental_method in ["X-ray diffraction", "Electron Microscopy"] and resolution <= 3) or experimental_method == "Solution NMR"'

    # complete_chains_run_as_completed = False

    # weight = array([1, -1, -1, -1.79072623, -2.95685934, -4.6231746])

    deletion_part_kwargs = dict()

    unp_head = re_compile(r'>sp\|(.+)\|')

    UniProtFASTA = UniProtINFO('fasta')

    def set_id(self):
        if self.source == 'PDB':
            self.for_get_id = self.identifier.lower()
            self.pdb_id = self.for_get_id
        elif self.source == 'UniProt':
            self.for_get_id = self.raw_identifier.upper()
        else:
            raise AssertionError(f'Invalid identifier: {self.raw_identifier} for {self.__class__.__name__}!')

    def get_id(self):
        return self.for_get_id

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.source} {self.level} {self.get_id()}>"

    @classmethod
    def fetch_unp_fasta(cls, identifier):
        task = cls.tasks.get((identifier, 'UniProtFASTA.single_retrieve(identifier).then(a_seq_reader)'), None)
        if task is None:    
            task = cls.UniProtFASTA.single_retrieve(identifier).then(a_seq_reader)
            cls.register_task((identifier, 'UniProtFASTA.single_retrieve(identifier).then(a_seq_reader)'), task)
        return task
    
    @unsync
    async def get_sequence(self, **kwargs):
        if self.source == 'UniProt':
            return (await SIFTS.fetch_unp_fasta(self.get_id()))[1]
        else:
            return await super().get_sequence(**kwargs)

    @classmethod
    @unsync
    async def complete_chains(cls, dfrm: Union[DataFrame, Unfuture, Coroutine]):
        if isawaitable(dfrm):
            dfrm = await dfrm
        #if cls.complete_chains_run_as_completed:
        #    res = await SIFTSs(dfrm.pdb_id.unique()).fetch('fetch_from_pdbe_api', 
        #                    api_suffix='api/mappings/all_isoforms/',
        #                    then_func=Base.to_dataframe).run()
        #else:
        res = [await task for task in SIFTSs(frozenset(dfrm.pdb_id)).fetch('fetch_from_pdbe_api', 
                    api_suffix='api/mappings/all_isoforms/',
                    then_func=Base.to_dataframe).tasks]
        return concat(res, sort=False, ignore_index=True)

    @staticmethod
    @unsync
    async def check_whether_pdbe_api_lag(dfrm):
        if isawaitable(dfrm):
            dfrm = await dfrm
        if dfrm is None:
            return
        elif isinstance(dfrm, DataFrame):
            pass
        else:
            return
        pdbs = PDBs(frozenset(dfrm.pdb_id))
        tasks = [await task for task in pdbs.fetch('fetch_from_pdbe_api', api_suffix='api/pdb/entry/molecules/').tasks]
        pass_pdbs = [i.name[-8:-4] for i in tasks if i is not None]
        mask = dfrm.pdb_id.isin(pass_pdbs)
        if mask.any():
            if (~mask).any():
                warn(f"{dfrm[~mask].pdb_id.tolist()}: either obsoleted or API lag update", PossibleObsoletedPDBEntryWarning)
            return dfrm[mask].reset_index(drop=True)
        else:
            return

    """@staticmethod
    @unsync
    async def check_pdb_status(dfrm):
        if isinstance(dfrm, Unfuture):
            dfrm = await dfrm
        elif dfrm is None:
            return
        elif isinstance(dfrm, DataFrame):
            pass
        else:
            return
        pdbs = PDBs(dfrm.pdb_id.unique())
        tasks = [await task for task in pdbs.fetch('fetch_from_pdbe_api', api_suffix='api/pdb/entry/status/', then_func=a_load_json, json=True).tasks]
        pass_pdbs = [next(iter(i)) for i in tasks if next(iter(i.values()))[0]['status_code'] == 'REL']
        res = dfrm[dfrm.pdb_id.isin(pass_pdbs)]
        if len(res) > 0:
            return res.reset_index(drop=True)
        else:
            return"""

    @staticmethod
    @unsync
    def generate_new_identity(dfrm):
        '''
        new_identity = #identical characters / #aligned columns including columns containing a gap in either sequence
        '''
        if isinstance(dfrm, Unfuture):
            dfrm = dfrm.result()
        if 'new_pdb_range_raw' in dfrm.columns:
            pdb_range_col = 'new_pdb_range_raw' 
            unp_range_col = 'new_unp_range_raw'
        else:
            pdb_range_col = 'new_pdb_range' 
            unp_range_col = 'new_unp_range'
        new_pdb_range_len = dfrm[pdb_range_col].apply(range_len)
        new_unp_range_len = dfrm[unp_range_col].apply(range_len) # TODO: drop
        assert (new_pdb_range_len==new_unp_range_len).all() # TODO: drop
        conflict_range_len = dfrm.conflict_pdb_range.apply(range_len)
        pdb_gaps = dfrm[pdb_range_col].apply(get_gap_list)
        unp_gaps = dfrm[unp_range_col].apply(get_gap_list)
        indel_len = pdb_gaps.apply(lambda x: sum(i for i in x if i > 0)) + unp_gaps.apply(lambda x: sum(i for i in x if i > 0))
        dfrm['new_identity'] = round((new_pdb_range_len-conflict_range_len)/(indel_len+new_pdb_range_len), 2)
        return dfrm

    @staticmethod
    @unsync
    def check_identity(dfrm: Union[DataFrame, Unfuture]):
        if isinstance(dfrm, Unfuture):
            dfrm = dfrm.result()
        check = dfrm.groupby(['UniProt','pdb_id','struct_asym_id']).identity.unique().apply(len).gt(1)
        res = check[check==True].to_frame().reset_index().rename(columns={'identity': 'multi_identity'})
        df = dfrm.merge(res, how='left')
        df.multi_identity = df.multi_identity.fillna(False)
        return df

    @staticmethod
    @unsync
    def reformat(dfrm: Union[DataFrame, Unfuture], drop_non_sequencial:bool=True) -> DataFrame:
        if isinstance(dfrm, Unfuture):
            dfrm = dfrm.result()
        if 'pdb_start' not in dfrm.columns:
            flat_dict_in_df(dfrm, dfrm.start.apply(json.loads), ('residue_number',))
            flat_dict_in_df(dfrm, dfrm.end.apply(json.loads), ('residue_number',))
            dfrm.rename(columns={
                'start.residue_number': 'pdb_start', 
                'end.residue_number': 'pdb_end'}, inplace=True)
        
        if drop_non_sequencial:
            """
            NOTE:
                {'pdb_id': '5o32', 'entity_id': 3, 'UniProt': 'P08603'}
                [[1, 246], [260, 383]]
                [[19, 264], [1107, 387]]
            """
            pass_mask = (dfrm.pdb_start <= dfrm.pdb_end) & (dfrm.unp_start <= dfrm.unp_end)
            if not pass_mask.all():
                warn(f"Drop:\n{dfrm[~pass_mask][['UniProt','pdb_id','struct_asym_id','pdb_start','pdb_end','unp_start','unp_end']]}")
                dfrm = dfrm[pass_mask].reset_index(drop=True)
        '''
        NOTE: sort for cases like P00441 5j0c (chain A,B)
        NOTE: hasn't handle multiple identity values!
        '''
        dfrm.sort_values(by=['UniProt', 'pdb_id', 'struct_asym_id', 'pdb_start'], inplace=True)
        group_info_col = ['pdb_id', 'chain_id', 'UniProt']
        range_info_col = ['pdb_start', 'pdb_end', 'unp_start', 'unp_end']
        reader = SeqRangeReader(group_info_col)
        dfrm[['pdb_range', 'unp_range']] = DataFrame(dfrm.apply(
            lambda x: reader.check(tuple(x[i] for i in group_info_col), tuple(
                x[i] for i in range_info_col)),
            axis=1).values.tolist(), index=dfrm.index)
        dfrm = dfrm.drop(columns=range_info_col).drop_duplicates(
            subset=group_info_col, keep='last').reset_index(drop=True)
        dfrm["Entry"] = dfrm["UniProt"].apply(lambda x: x.split('-')[0])
        return dfrm
    
    @staticmethod
    def add_tage_to_range(df: DataFrame, tage_name: str):
        # ADD TAGE FOR SIFTS
        df[tage_name] = 'Safe'
        # No Insertion But Deletion[Pure Deletion]
        df.loc[df[(df['group_info'] == 1) & (
            df['diff+'] > 0)].index, tage_name] = 'Deletion'
        # Insertion & No Deletion
        df.loc[df[
            (df['group_info'] == 1) &
            (df['diff-'] > 0)].index, tage_name] = 'Insertion_Undivided'
        df.loc[df[
            (df['group_info'] > 1) &
            (df['diff0'] == df['group_info']) &
            (df['unp_gaps0'] == (df['group_info'] - 1))].index, tage_name] = 'Insertion'
        # Insertion & Deletion
        df.loc[df[
            (df['group_info'] > 1) &
            (df['diff0'] == df['group_info']) &
            (df['unp_gaps0'] != (df['group_info'] - 1))].index, tage_name] = 'InDel_1'
        df.loc[df[
            (df['group_info'] > 1) &
            (df['diff0'] != df['group_info']) &
            (df['unp_gaps0'] != (df['group_info'] - 1))].index, tage_name] = 'InDel_2'
        df.loc[df[
            (df['group_info'] > 1) &
            (df['diff0'] != df['group_info']) &
            (df['unp_gaps0'] == (df['group_info'] - 1))].index, tage_name] = 'InDel_3'

    @classmethod
    @unsync
    def dealWithInDel(cls, dfrm: Union[DataFrame, Unfuture], sort_by_unp: bool = True) -> DataFrame:
        if isinstance(dfrm, Unfuture):
            dfrm = dfrm.result()

        dfrm.pdb_range = dfrm.pdb_range.apply(json.loads)
        dfrm.unp_range = dfrm.unp_range.apply(json.loads)
        dfrm['group_info'] = dfrm.apply(lambda x: len(
            x['pdb_range']), axis=1)

        focus_index = dfrm[dfrm.group_info.gt(1)].index
        if sort_by_unp and (len(focus_index) > 0):
            focus_df = dfrm.loc[focus_index].apply(lambda x: sort_2_range(x['unp_range'], x['pdb_range']), axis=1).apply(Series)
            focus_df.index = focus_index
            focus_df.columns = ['unp_range', 'pdb_range']
            dfrm.loc[focus_index, ['unp_range', 'pdb_range']] = focus_df

        dfrm['pdb_gaps'] = dfrm.pdb_range.apply(get_gap_list)
        dfrm['unp_gaps'] = dfrm.unp_range.apply(get_gap_list)
        dfrm['range_diff'] = dfrm.apply(lambda x: get_range_diff(
            x['unp_range'], x['pdb_range']), axis=1)
        dfrm['diff0'] = dfrm.range_diff.apply(
            lambda x: count_nonzero(x == 0))
        dfrm['diff+'] = dfrm.range_diff.apply(
            lambda x: count_nonzero(x > 0))
        dfrm['diff-'] = dfrm.range_diff.apply(
            lambda x: count_nonzero(x < 0))
        dfrm['unp_gaps0'] = dfrm.unp_gaps.apply(lambda x: x.count(0))
        cls.add_tage_to_range(dfrm, tage_name='sifts_range_tag')
        dfrm['repeated'] = dfrm.apply(
            lambda x: x['diff-'] > 0 and x['sifts_range_tag'] != 'Insertion_Undivided', axis=1)
        dfrm['repeated'] = dfrm.apply(
            lambda x: True if any(i < 0 for i in x['unp_gaps']) else x['repeated'], axis=1)
        dfrm['reversed'] = dfrm.pdb_gaps.apply(lambda x: any(i < 0 for i in x))
        dfrm.pdb_range = dfrm.pdb_range.apply(
            lambda x: json.dumps(x).decode('utf-8'))
        dfrm.unp_range = dfrm.unp_range.apply(
            lambda x: json.dumps(x).decode('utf-8'))
        #dfrm['InDel_sum'] = dfrm.pdb_gaps.apply(sum) + dfrm.unp_gaps.apply(sum) + dfrm.range_diff.apply(sum)
        temp_cols = ['start', 'end', 'group_info', 'pdb_gaps', 'unp_gaps',
                     'diff0', 'diff+', 'diff-', 'unp_gaps0']
        return dfrm.drop(columns=temp_cols), dfrm[temp_cols]

    @unsync
    async def summary_ec_ue(self):
        def regist_info(store, record):
            store[self.MappingMeta(record['UniProt'], record['name'].split('_')[-1], record['identity'])][record['entity_id']] = record['unp_range']

        def reverse_dict(info_dict):
            res = defaultdict(list)
            for key, valueDict in info_dict.items():
                for innerKey in valueDict.keys():
                    res[innerKey].append(key)
            return res

        def min_unp(info_dict):
            return min(len(frozenset(res)) for res in product(*info_dict.values()))
        
        mol_df = await self.fetch_from_pdbe_api('api/pdb/entry/molecules/', Base.to_dataframe)
        mol_df = mol_df[mol_df.molecule_type.eq('polypeptide(L)')]
        info_dict = dict(zip(mol_df.entity_id,mol_df.in_chains))
        entity_count = len(info_dict)
        chain_count = sum(i.count(',')+1 for i in info_dict.values())

        try:
            best_iso = await self.fetch_from_pdbe_api('api/mappings/isoforms/', Base.to_dataframe).then(self.reformat)
            # best_iso = best_iso.sort_values(by='identity', ascending=False).drop_duplicates(subset='entity_id')
            unp_entity_info = defaultdict(dict)
            best_iso.apply(lambda x: regist_info(unp_entity_info, x), axis=1)
            entity_unp_info = reverse_dict(unp_entity_info)
            entity_with_unp_count = len(entity_unp_info)
            min_unp_count = min_unp(entity_unp_info)
        except AttributeError:
            unp_entity_info, entity_unp_info, entity_with_unp_count, min_unp_count = nan, nan, 0, 0

        return self.EntityChain(self.pdb_id, info_dict, entity_count, chain_count), self.UniProtEntity(self.pdb_id, unp_entity_info, entity_unp_info, entity_with_unp_count, min_unp_count)

    @unsync
    async def get_oligo_state(self):
        ec, ue = await self.summary_ec_ue()
        assert ec.entity_count >= ue.entity_with_unp_count
        has_unmapped_protein = ec.entity_count > ue.entity_with_unp_count
        if ((ec.entity_count > ue.entity_with_unp_count) and (ue.min_unp_count == 1)) or (ue.min_unp_count > 1):
            return self.OligoState(self.pdb_id, 'he', has_unmapped_protein)
        elif (ec.entity_count == ue.entity_with_unp_count) and (ue.min_unp_count == 1) and (ec.chain_count > 1):
            return self.OligoState(self.pdb_id, 'ho', has_unmapped_protein)
        elif ec.chain_count == 1:
            return self.OligoState(self.pdb_id, 'mo', has_unmapped_protein)
        else:
            raise ValueError('Unexpected cases for get_oligo_state!')

    @staticmethod
    def convert_index(lrange: Union[List, Tuple, str], rrange: Union[List,Tuple, str], sites: Iterable[int], nan_value=None):
        '''convert from rrange to lrange'''
        def unit(lrange, rrange, site):
            for (lbeg, lend), (rbeg, rend) in zip(lrange, rrange):
                assert (lend - lbeg) == (rend-rbeg), "convert_index(): Invalid range {} {}".format(lrange, rrange)
                if (site >= rbeg) and (site <= rend):
                    return int(site + lbeg - rbeg)
                else:
                    continue
            return nan_value
        lrange = json.loads(lrange) if isinstance(lrange, str) else lrange
        rrange = json.loads(rrange) if isinstance(rrange, str) else rrange
        assert len(lrange) == len(rrange), "{} {}".format(lrange, rrange)
        return tuple(unit(lrange, rrange, site) for site in sites)

    @classmethod
    @unsync
    async def get_residue_conflict(cls, pdb_id, entity_id, pdb_range, Entry, UniProt, unp_range):
        pdb_seq = await PDB(pdb_id).get_sequence(entity_id=entity_id)
        try:
            _, unp_seq = await cls.fetch_unp_fasta(UniProt)
        except TypeError:
            raise PossibleObsoletedUniProtError(UniProt)
        indexes = tuple(get_diff_index(pdb_seq, pdb_range, unp_seq, unp_range))
        if len(indexes) > 0:
            pdb_diff_index, unp_diff_index = zip(*indexes)
            return (
                json.dumps(dict(zip((str(i) for i in pdb_diff_index), (unp_seq[i-1] for i in unp_diff_index)))).decode('utf-8'),
                json.dumps(dict(zip((str(i) for i in pdb_diff_index), (pdb_seq[i-1] for i in pdb_diff_index)))).decode('utf-8'), 
                json.dumps(to_interval(pdb_diff_index)).decode('utf-8'), 
                json.dumps(to_interval(unp_diff_index)).decode('utf-8'), 
                len(unp_seq))
        else:
            return ('{}', '{}', '[]', '[]', len(unp_seq)) 

    @classmethod
    @unsync
    async def add_residue_conflict(cls, dfrm: Union[DataFrame, Tuple, Unfuture, Coroutine]):
        '''
        TODO: optimization
        '''
        if isawaitable(dfrm):
            dfrm = await dfrm
        if isinstance(dfrm, Tuple):
            dfrm = dfrm[0]
        pdb_range_col = 'new_pdb_range'# if 'new_pdb_range' in dfrm.columns else 'pdb_range'
        unp_range_col = 'new_unp_range'# if 'new_unp_range' in dfrm.columns else 'unp_range'
        focus = ['UniProt', 'entity_id', 'pdb_id', 'pdb_range', 'unp_range']
        '''
                (UniProt	chain_id	entity_id	identifier	identity	is_canonical	name	pdb_id	struct_asym_id	pdb_range	unp_range)
        NOTE: add pdb_range because of (P00720	B	1	ENLYS_BPT4	0.94	True	ENLYS_BPT4	2b7x	B	[[1,24],[31,170]]	[[1,24],[25,164]])
                                       (P00720	A	1	ENLYS_BPT4	0.94	True	ENLYS_BPT4	2b7x	A	[[1,30],[37,170]]	[[1,30],[31,164]])
        
              same entity_id with different mapped range
        '''
        f_dfrm = dfrm[focus+[pdb_range_col, 'Entry', unp_range_col]].drop_duplicates(subset=focus).reset_index(drop=True)
        tasks = f_dfrm.apply(lambda x: cls.get_residue_conflict(
            x['pdb_id'], x['entity_id'], x[pdb_range_col], x['Entry'], x['UniProt'], x[unp_range_col]), axis=1)
        f_dfrm['conflict_pdb_index'], f_dfrm['raw_pdb_index'], f_dfrm['conflict_pdb_range'], f_dfrm['conflict_unp_range'], f_dfrm['unp_len'] = zip(*[await i for i in tasks])
        dfrm_ed = merge(dfrm, f_dfrm)
        assert dfrm_ed.shape[0] == dfrm.shape[0], f"\n{dfrm_ed.shape}\n{dfrm.shape}\n{f_dfrm.shape}"
        return dfrm_ed

    @staticmethod
    @unsync
    async def rsmfga_unit(UniProt, pdb_id, entity_id, start, end):
        '''rsmfga: short for renew_sifts_mapping_from_graph_api'''
        df = await SIFTS(pdb_id).fetch_residue_mapping(entity_id=entity_id, start=start, end=end, UniProt=UniProt)
        if df is None:
            raise AttributeError(f"{(UniProt, pdb_id, entity_id, start, end)}")
        return lyst22interval(df.unp_residue_number, df.residue_number)

    @classmethod
    @unsync
    async def renew_sifts_mapping_from_graph_api(cls, UniProt, pdb_id, entity_id, pdb_range, unp_range, range_diff):
        if isinstance(pdb_range, str):
            pdb_range = json.loads(pdb_range)
        if isinstance(unp_range, str):
            unp_range = json.loads(unp_range)
        new_unp_range, new_pdb_range = [], []
        for diff, urange, prange in zip(range_diff, unp_range, pdb_range):
            if diff == 0:
                new_unp_range.append(urange)
                new_pdb_range.append(prange)
            else:
                try:
                    newurange, newprange = await cls.rsmfga_unit(UniProt, pdb_id, entity_id, *prange)
                    new_unp_range.extend(newurange)
                    new_pdb_range.extend(newprange)
                except AttributeError as e:
                    warn(str(e), PDBeKBResidueMappingErrorWarning)
        return json.dumps(new_unp_range).decode('utf-8'), json.dumps(new_pdb_range).decode('utf-8')

    @staticmethod
    @unsync
    async def deal_with_identical_entity_seq(dfrm):
        if isawaitable(dfrm):
            dfrm = await dfrm
        #already = set()
        #cluster_dfs = []
        dfrm = dfrm.copy()
        dfrm['pdb_sequence'] = b''
        dfrm_nr = dfrm[['pdb_id', 'entity_id']].drop_duplicates()
        for pdb_id, entity_id in zip(dfrm_nr.pdb_id, dfrm_nr.entity_id):
            dfrm.loc[dfrm[dfrm.pdb_sequence.eq(b'') & dfrm.pdb_id.eq(pdb_id) & dfrm.entity_id.eq(entity_id)].index, 'pdb_sequence'] = zlib.compress(bytes(await PDB(pdb_id).get_sequence(entity_id=entity_id, mode='raw_pdb_seq'), encoding='utf-8'))
            """if (pdb_id, entity_id) in already:
                continue
            cur_cluster_df = await PDB(pdb_id).rcsb_cluster_membership(entity_id=entity_id, identity_cutoff=100)
            try:
                assert cur_cluster_df is not None
                already |= set(zip(cur_cluster_df.pdb_id, cur_cluster_df.entity_id))
            except AssertionError:
                cur_cluster_df = DataFrame([dict(pdb_id=pdb_id, entity_id=entity_id, cluster_id=-1)])
            cluster_dfs.append(cur_cluster_df)"""

        #cluster_df = concat(cluster_dfs, sort=False, ignore_index=True)
        #assert not any(cluster_df.duplicated())
        #dfrm = dfrm.merge(cluster_df[['pdb_id','entity_id','cluster_id']], how='left')
        #assert not any(dfrm.cluster_id.isnull()), f"{dfrm[dfrm.cluster_id.isnull()]}"
        #dfrm['fix_cluster_id'] = dfrm.groupby(['cluster_id', 'pdb_sequence']).ngroup().astype(str) + '_' + dfrm.cluster_id.astype(str)
        dfrm['fix_cluster_id'] = dfrm.groupby('pdb_sequence').ngroup()
        # ignore/overried cases like (P00720,2b7x,B v.s P00720,2b7x,A)
        return dfrm.drop(columns=['pdb_sequence'])

    @classmethod
    @unsync
    async def double_check_conflict_and_range(cls, dfrm: Union[DataFrame, Unfuture, Coroutine]):
        if isawaitable(dfrm):
            dfrm = await dfrm
        focus_part = dfrm[
            dfrm.sifts_range_tag.isin(('Deletion', 'Insertion_Undivided', 'InDel_2', 'InDel_3')) &
            (dfrm.conflict_pdb_index.apply(get_str_dict_len)/dfrm.new_pdb_range.apply(range_len)).ge(0.1)]
        if len(focus_part) == 0:
            return dfrm
        focus_part_iden = await cls.deal_with_identical_entity_seq(focus_part)
        focus_part_iden_dd = focus_part_iden.drop_duplicates(subset=['UniProt', 'fix_cluster_id']).copy()
        tasks = tuple(map(cls.renew_sifts_mapping_from_graph_api, focus_part_iden_dd.UniProt, focus_part_iden_dd.pdb_id, focus_part_iden_dd.entity_id, focus_part_iden_dd.pdb_range, focus_part_iden_dd.unp_range, focus_part_iden_dd.range_diff))
        focus_part_iden_dd[['new_unp_range', 'new_pdb_range']] = [await task for task in tasks]
        focus_part_iden_dd = await cls.add_residue_conflict(focus_part_iden_dd.drop(columns=['conflict_pdb_index', 'raw_pdb_index', 'conflict_pdb_range', 'conflict_unp_range']))
        focus_part_iden_dd.reversed = focus_part_iden_dd.new_pdb_range.apply(get_gap_list).apply(lambda x: any(i < 0 for i in x))
        focus_cols = ['UniProt', 'fix_cluster_id', 'new_unp_range', 'new_pdb_range', 'reversed',
                      'conflict_pdb_index', 'raw_pdb_index', 'conflict_pdb_range', 'conflict_unp_range']
        res = focus_part_iden.drop(columns=focus_cols[2:]).merge(focus_part_iden_dd[focus_cols], how='left')
        assert res.isnull().sum().sum() == 0
        #res = res.drop(columns=['fix_cluster_id', 'cluster_id'])
        res = res.drop(columns=['fix_cluster_id'])
        assert res.shape == focus_part.shape, f"{res.shape}, {focus_part.shape}"
        res.index = focus_part.index
        dfrm.loc[focus_part.index] = res
        return dfrm

    @unsync
    async def fetch_new_residue_mapping(self, entity_id, start, end, chunksize=50):
        tasks = []
        for index in range(0, end-start+1, chunksize):
            cur_start = start+index
            cur_end = cur_start+chunksize-1
            if cur_end > end:
                cur_end = end
            task = self.fetch_from_pdbe_api(
                'graph-api/residue_mapping/', 
                Base.to_dataframe, 
                mask_id=f'{self.pdb_id}/{entity_id}/{cur_start}/{cur_end}')
            tasks.append(task)
        dfs = [await i for i in tasks]
        dfs = [i for i in dfs if i is not None]
        if dfs:
            df = concat(dfs, sort=False, ignore_index=True)
            await self.sqlite_api.async_insert_chunk(self.sqlite_api.ResidueMapping, df.to_dict('records'))
            return df

    @unsync
    async def fetch_residue_mapping(self, entity_id:int, start:int, end:int, columns:str='UniProt,residue_number,unp_residue_number', UniProt:Optional[str]=None):
        unp = "" if UniProt is None else f"UniProt == '{UniProt}' AND"
        query=f"""
            SELECT DISTINCT {columns} FROM ResidueMapping
            WHERE {unp} pdb_id == '{self.pdb_id}' AND entity_id == {entity_id}
                AND residue_number >= {start} AND residue_number <= {end}"""
        res = await self.sqlite_api.database.fetch_all(query=query)
        columns = columns.split(',')

        if len(res) > 0:
            df = DataFrame(res, columns=columns)
            start_diff = start - df.residue_number.min()
            end_diff = end - df.residue_number.max()
            new_range = []
            dfs = [df]
            if start_diff < 0:
                new_range.append((start, start-start_diff-1))
            if end_diff > 0:
                new_range.append((end-end_diff+1, end))
            for s, e in new_range:
                new_df = await self.fetch_new_residue_mapping(entity_id, s, e)
                new_df = new_df[new_df.UniProt.eq(UniProt)] if ((new_df is not None) and (UniProt is not None)) else new_df
                dfs.append(new_df[columns])
            return concat(dfs, ignore_index=True)
        else:
            df = await self.fetch_new_residue_mapping(entity_id, start, end)
            if df is not None:
                df = df[df.UniProt.eq(UniProt)] if UniProt is not None else df
                return df[columns]
            else:
                return

    @classmethod
    @unsync
    async def fix_range(cls, dfrm: Union[DataFrame, Tuple, Unfuture, Coroutine]):
        if isawaitable(dfrm):
            dfrm = await dfrm
        if isinstance(dfrm, Tuple):
            dfrm = dfrm[0]
        '''
        NOTE: same entity with different range
            UniProt         P42262	P42262
            chain_id	    A	B
            entity_id	    1	1
            identity	    0.66	0.66
            is_canonical	True	True
            pdb_id	        2xhd	2xhd
            struct_asym_id	A	B
            pdb_range	    [[3,263]]	[[2,263]]
            unp_range	    [[413,796]]	[[412,796]]
        '''
        focus = ['UniProt', 'entity_id', 'pdb_id', 'pdb_range', 'unp_range']
        f_dfrm = dfrm[focus+['Entry', 'range_diff', 'sifts_range_tag']].drop_duplicates(subset=focus)
        f_dfrm = f_dfrm[f_dfrm.sifts_range_tag.isin(('Deletion', 'Insertion_Undivided', 'InDel_2', 'InDel_3'))].reset_index(drop=True)

        if len(f_dfrm) > 0:
            tasks = f_dfrm.apply(lambda x: cls.a_re_align(
                x['range_diff'], x['pdb_range'], x['unp_range'], x['pdb_id'], x['entity_id'], x['UniProt']), axis=1)
            res = [await i for i in tasks]
            f_dfrm[['new_pdb_range', 'new_unp_range']] = DataFrame(list(zip(*i)) for i in res)
            #assert all(cls.check_range_tail(*args) for args in zip(f_dfrm.new_pdb_range, f_dfrm.new_unp_range, f_dfrm.pdb_range))
            #f_dfrm[['new_pdb_range', 'new_unp_range']] = DataFrame([cls.check_range_tail(*args) for args in zip(f_dfrm.new_pdb_range, f_dfrm.new_unp_range, f_dfrm.pdb_range)])
            dfrm_ed = merge(dfrm, f_dfrm.drop(columns=['range_diff']), how='left')
            assert dfrm_ed.shape[0] == dfrm.shape[0]
            dfrm_ed.new_pdb_range = dfrm_ed.apply(lambda x: x['pdb_range'] if isna(x['new_pdb_range']) else x['new_pdb_range'],axis=1)
            dfrm_ed.new_unp_range = dfrm_ed.apply(lambda x: x['unp_range'] if isna(x['new_unp_range']) else x['new_unp_range'],axis=1)
            return dfrm_ed
        else:
            dfrm_ed = dfrm.copy()
            dfrm_ed['new_pdb_range'] = dfrm.pdb_range
            dfrm_ed['new_unp_range'] = dfrm.unp_range
            return dfrm_ed

    @staticmethod
    def re_alignment(range_diff, pdb_seq, pdb_range, unp_seq, unp_range, **kwargs):
        def get_optimal_range(lseq, rseq, lbeg, lend, rbeg, rend):
            '''
            NOTE:
                UniProt	chain_id	entity_id	pdb_id	struct_asym_id	pdb_range	unp_range	range_diff	new_pdb_range	new_unp_range	
                P36022	A	1	4w8f	A	[[2,1676],[1683,1684],[1849,2649]]	[[1364,3038],[3064,3076],[3292,4092]]	[0, 11, 0]	((2, 1676), (1683, 1684), (1685, 1684), (1849,...	((1364, 3038), (3064, 3065), (3077, 3076), (32...
            '''
            aln = nw_trace_scan_sat(lseq, rseq, 10, 1, blosum62)
            cur_lbeg = lbeg
            cur_rbeg = rbeg
            cur_lend = lend
            cur_rend = rend
            for v in aln.cigar.seq:
                seg_len = aln.cigar.decode_len(v)
                seg_type = aln.cigar.decode_op(v)
                if seg_type == b'=' or seg_type == b'X':
                    cur_lend = cur_lbeg+seg_len-1
                    cur_rend = cur_rbeg+seg_len-1
                    yield (cur_lbeg, cur_lend), (cur_rbeg, cur_rend)
                    cur_lbeg = cur_lend + 1
                    cur_rbeg = cur_rend + 1
                elif seg_type == b'I':
                    cur_lbeg += seg_len
                elif seg_type == b'D':
                    cur_rbeg += seg_len
                else:
                    raise ValueError(f'Unexpected type: {seg_type}')
        """
        NOTE:
            {'pdb_id': '5o32', 'entity_id': 3, 'UniProt': 'P08603'}
            [[1, 246], [260, 383]]
            [[19, 264], [1107, 387]]
        """
        for diff, (lbeg, lseg), (rbeg, rseg) in zip(range_diff, get_seq_seg(pdb_seq, pdb_range, **kwargs), get_seq_seg(unp_seq, unp_range, **kwargs)):

            # assert bool(lseg) and bool(rseg), f"{pdb_range}\n{unp_range}"

            lend = lbeg + len(lseg) - 1
            rend = rbeg + len(rseg) - 1

            if diff != 0:
                yield from get_optimal_range(lseg, rseg, lbeg, lend, rbeg, rend)
            else:
                yield (lbeg, lend), (rbeg, rend)

    @classmethod
    @unsync
    async def a_re_align(cls, range_diff, pdb_range, unp_range, pdb_id, entity_id, UniProt):
        pdb_seq = await PDB(pdb_id).get_sequence(entity_id=entity_id)
        try:
            _, unp_seq = await cls.fetch_unp_fasta(UniProt)
        except TypeError:
            raise PossibleObsoletedUniProtError(UniProt)
        pdb_range = json.loads(pdb_range) if isinstance(pdb_range, str) else pdb_range
        unp_range = json.loads(unp_range) if isinstance(unp_range, str) else unp_range
        return cls.re_alignment(range_diff, pdb_seq, pdb_range, unp_seq, unp_range, 
                pdb_id=pdb_id,entity_id=entity_id,UniProt=UniProt)

    @unsync
    async def pipe_base(self, complete_chains:bool=False, skip_pdbs=None, only_canonical:bool=False, skip_carbohydrate_polymer:bool=False):
        init_task = await self.fetch_from_pdbe_api('api/mappings/all_isoforms/', Base.to_dataframe).then(self.check_whether_pdbe_api_lag)
        if init_task is None:
            return
        elif skip_pdbs is not None:
            init_task = init_task[~init_task.pdb_id.isin(skip_pdbs)].reset_index(drop=True)
            if len(init_task) == 0:
                return
        if skip_carbohydrate_polymer:
            pdbs = PDBs(frozenset(init_task.pdb_id))
            tasks = [await task for task in pdbs.fetch(partial(PDBs.fetch_number_of_entities, 'carbohydrate_polymer')).tasks]
            cpdf = DataFrame(tasks, columns=['pdb_id', 'num_cp'])
            pass_pdbs = cpdf[cpdf.num_cp.eq(0)][['pdb_id']]
            init_task = init_task.merge(pass_pdbs).reset_index(drop=True)
            if len(init_task) == 0:
                return
        if complete_chains:
            init_task = await self.complete_chains(init_task)
            if len(init_task) == 0:
                return
        if only_canonical:
            init_task = init_task[init_task.is_canonical.eq(True)].reset_index(drop=True)
            if len(init_task) == 0:
                return
        try:
            sifts_df = await self.reformat(init_task
                ).then(self.dealWithInDel
                ).then(self.fix_range
                ).then(self.add_residue_conflict)
                #).then(self.double_check_conflict_and_range)
        except PossibleObsoletedUniProtError as e:
            warn(f'{repr(self)}, {e}', PossibleObsoletedUniProtWarning)
            return None
        return sifts_df
    
    @unsync
    async def pipe_score(self, sifts_df=None, complete_chains:bool=False, skip_pdbs=None, skip_carbohydrate_polymer=False):
        if sifts_df is None:
            sifts_df = await self.pipe_base(complete_chains=complete_chains, skip_pdbs=skip_pdbs, skip_carbohydrate_polymer=skip_carbohydrate_polymer)
        if sifts_df is None:
            return
        exp_cols = ['pdb_id', 'resolution', 'experimental_method_class',
                    'experimental_method', 'multi_method', '-r_factor', '-r_free']

        if self.source == 'PDB':
            return await self.pipe_score_for_pdb_entry(sifts_df, exp_cols)
        else:
            return await self.pipe_score_for_unp_isoform(sifts_df, exp_cols)

    @staticmethod
    def bs_score_aligned_part(new_pdb_range, conflict_pdb_range, conflict_pdb_index, raw_pdb_index, NON_INDEX, OBS_RATIO_ARRAY):
        assert new_pdb_range is not None
        L_aligned = range_len(new_pdb_range)
        un_range = add_range(conflict_pdb_range, NON_INDEX)
        aligned_equal_range = subtract_range(new_pdb_range, un_range)
        assert aligned_equal_range is not None
        aligned_unequal_range = overlap_range(new_pdb_range, un_range)

        aligned_unequal_score = 0
        conflict_pdb_index = json.loads(conflict_pdb_index) if isinstance(conflict_pdb_index, str) else conflict_pdb_index
        raw_pdb_index = json.loads(raw_pdb_index) if isinstance(raw_pdb_index, str) else raw_pdb_index
        non_set = frozenset(expand_interval(NON_INDEX))
        for i in expand_interval(aligned_unequal_range):
            unp_aa = conflict_pdb_index.get(i, None)
            pdb_aa = raw_pdb_index.get(i, None)
            if (i not in non_set) and (unp_aa is not None) and (unp_aa != pdb_aa):
                # NOT Modified & Conflict Residues are fall into here
                theta = miyata_similarity_matrix.get((unp_aa, pdb_aa), -3.104)
            else:
                # UNK | Modified Residue
                theta = -3.104
            aligned_unequal_score += OBS_RATIO_ARRAY[i-1]*(theta+1)
        
        return -L_aligned + \
               2*OBS_RATIO_ARRAY[[i-1 for i in expand_interval(aligned_equal_range)]].sum() + \
               aligned_unequal_score
    
    @staticmethod
    def bs_score_CN_terminal_part(new_pdb_range, ARTIFACT_INDEX, SEQRES_COUNT, OBS_RATIO_ARRAY):
        outside_range_ignore_artifact = subtract_range(outside_range(new_pdb_range, SEQRES_COUNT), ARTIFACT_INDEX)
        return -OBS_RATIO_ARRAY[[i-1 for i in expand_interval(outside_range_ignore_artifact)]].sum(), outside_range_ignore_artifact
    
    @classmethod
    def bs_score_insertion_part(cls, new_pdb_range, OBS_RATIO_ARRAY):
        insertion_range = cls.get_InDel_part(new_pdb_range)
        if len(insertion_range) == 0:
            return 0
        else:
            return -(OBS_RATIO_ARRAY[[i-1 for i in expand_interval(insertion_range)]].sum())-range_len(insertion_range)

    @classmethod
    @unsync
    async def bs_score_deletion_part(cls, pdb_id, struct_asym_id, new_pdb_range, new_unp_range, OBS_INDEX):
        new_unp_range = json.loads(new_unp_range) if isinstance(new_unp_range, str) else new_unp_range
        deletion_range = cls.get_InDel_part(new_unp_range)
        if len(deletion_range) == 0:
            return 0
        del_range_len = range_len(deletion_range)
        deletion_edge = cls.convert_index(new_pdb_range, new_unp_range, [i-1 for i, j in deletion_range])
        obs_lyst = tuple(expand_interval(OBS_INDEX))
        assert len(obs_lyst) > 0, f"{pdb_id} {struct_asym_id}, {new_pdb_range}"
        iedge = (bisect_left(obs_lyst, edge) for edge in deletion_edge)
        score = 0
        for i, j in zip(iedge, deletion_edge):
            if i == 0:
                edge = obs_lyst[0]
            elif i == len(obs_lyst):
                edge = obs_lyst[-1]
            else:
                this = obs_lyst[i]
                that = obs_lyst[i-1]
                if this - j <= j - that:
                    edge = this
                else:
                    edge = that
            score += await cls.clostest_distance_sum(pdb_id, struct_asym_id, edge, **cls.deletion_part_kwargs)
        return -del_range_len*score
        
    @classmethod
    @unsync
    async def clostest_distance_sum(cls, pdb_id, struct_asym_id, residue_number, radius=5, modelId=1, d=10):
        try:
            df = await PDB(pdb_id).fetch_from_coordinateServer_api(
                'ambientResidues', then_func=PDB.cif2atom_sites_df,
                asymId=struct_asym_id, seqNumber=residue_number, radius=radius, atomSitesOnly=1, modelId=modelId)
        except TypeError:
            df = None
        if df is None:
            try:
                df = await PDB(pdb_id).fetch_from_coordinateServer_api(
                    'ambientResidues', then_func=PDB.cif2atom_sites_df, root='ebi',
                    asymId=struct_asym_id, seqNumber=residue_number, radius=radius, atomSitesOnly=1, modelId=modelId)
            except TypeError:
                raise AssertionError(f"{pdb_id} {struct_asym_id} {residue_number}")
        df = df[df['_atom_site.label_comp_id'].ne('HOH') & df['_atom_site.type_symbol'].ne('H')].reset_index(drop=True)
        '''
        # for some ligand | carbohydrate
        df['_atom_site.label_seq_id'] = df['_atom_site.label_seq_id'].astype(int)
        assert all(~df['_atom_site.label_seq_id'].eq('.')), f"{pdb_id} {struct_asym_id} {residue_number}"
        '''
        coordinates_dict = cls.get_coordinates_dict(df)
        nbrs = NearestNeighbors(n_neighbors=1).fit(coordinates_dict[(struct_asym_id, str(residue_number))])
        total = 0
        for key in coordinates_dict.keys() - {(struct_asym_id, str(residue_number))}:
            distances, _ = nbrs.kneighbors(coordinates_dict[key])
            total += cls.exp_score(distances.min(), d)
        return total
    
    @staticmethod
    def exp_score(x, d):
        return exp(-square(x)/d)

    @staticmethod
    def get_InDel_part(li):
        '''
        when li is new_pdb_range, the return value is insertion_range
        when li is new_unp_range, the return value is deletion_range
        '''
        res = ((li[i][1]+1, li[i+1][0]-1) for i in range(len(li)-1))
        return [[i, j] for i, j in res if i<=j]

    @classmethod
    @unsync
    async def bs_score_base(cls, *args):
        if len(args) == 1:
            pdb_id, struct_asym_id, new_pdb_range, new_unp_range, conflict_pdb_range, conflict_pdb_index, raw_pdb_index, SEQRES_COUNT, ARTIFACT_INDEX, OBS_INDEX, NON_INDEX, OBS_RATIO_ARRAY = args[0]
        elif len(args) == 12:
            pdb_id, struct_asym_id, new_pdb_range, new_unp_range, conflict_pdb_range, conflict_pdb_index, raw_pdb_index, SEQRES_COUNT, ARTIFACT_INDEX, OBS_INDEX, NON_INDEX, OBS_RATIO_ARRAY = args
        else:
            raise TypeError("bs_score_base() with invalid positional arguments")
        new_pdb_range = json.loads(new_pdb_range) if isinstance(new_pdb_range, str) else new_pdb_range
        OBS_RATIO_ARRAY = array(json.loads(zlib.decompress(OBS_RATIO_ARRAY))) if isinstance(OBS_RATIO_ARRAY, bytes) else OBS_RATIO_ARRAY
        aligned_part = cls.bs_score_aligned_part(new_pdb_range, conflict_pdb_range, conflict_pdb_index, raw_pdb_index, NON_INDEX, OBS_RATIO_ARRAY)
        CN_terminal_part, outside_range_ignore_artifact = cls.bs_score_CN_terminal_part(new_pdb_range, ARTIFACT_INDEX, SEQRES_COUNT, OBS_RATIO_ARRAY)
        insertion_part = cls.bs_score_insertion_part(new_pdb_range, OBS_RATIO_ARRAY)
        deletion_part = await cls.bs_score_deletion_part(pdb_id, struct_asym_id, new_pdb_range, new_unp_range, OBS_INDEX)
        return aligned_part, CN_terminal_part, insertion_part, deletion_part, outside_range_ignore_artifact

    @staticmethod
    def wrap_trim_range(iobs_range, ipdb_range, iunp_range, irepeated, ireversed):
        if irepeated or ireversed:
            return ipdb_range, iunp_range
        else:
            return trim_range(iobs_range, ipdb_range, iunp_range)

    @unsync
    async def pipe_score_base(self, sifts_df, pec_df):
        full_df = sifts_df.merge(pec_df)
        try:
            assert sifts_df.shape[0] == full_df.shape[0]
            # {full_df.duplicated(subset=['UniProt','pdb_id','struct_asym_id'], keep=False)}
            # f"{pec_df[pec_df.pdb_id.eq('6mzl')]}"
            '''
            FOR O00330 1zy8 struct_asym_id O
            '''
        except AssertionError:
            warn(f"{repr(self)}: {sifts_df.shape}, {full_df.shape}", ConflictChainIDWarning)
            full_df = sifts_df.merge(pec_df.drop(columns=['chain_id']), how='left')
            for_check = full_df[full_df.OBS_STD_COUNT.isnull()][['UniProt','pdb_id','struct_asym_id']]
            assert sifts_df.shape[0] == full_df.shape[0], f"\n{sifts_df.shape[0]} != {full_df.shape[0]}\n{for_check}"
        
        full_df[['new_pdb_range_raw', 'new_unp_range_raw']] = full_df[['new_pdb_range', 'new_unp_range']]
        full_df[['new_pdb_range', 'new_unp_range']] = DataFrame([self.wrap_trim_range(
            iobs_range, ipdb_range, iunp_range, irepeated, ireversed) for iobs_range, ipdb_range, iunp_range, irepeated, ireversed in zip(
                full_df.OBS_INDEX, full_df.new_pdb_range, full_df.new_unp_range, full_df.repeated, full_df.reversed)])
        grouped_len = full_df.groupby('UniProt').new_unp_range.apply(partial(reduce, add_range)).apply(range_len)
        grouped_len.name = 'c_unp_len'
        grouped_len = grouped_len.to_frame().reset_index()
        # assert frozenset(grouped_len.columns) == frozenset({'UniProt', 'c_unp_len'}), str(grouped_len.columns)
        full_df = full_df.merge(grouped_len, how='left')
        assert full_df.c_unp_len.isnull().sum() == 0, f"{full_df[full_df.c_unp_len.isnull()].UniProt}"

        cols = ['pdb_id', 'struct_asym_id', 'new_pdb_range', 'new_unp_range', 'conflict_pdb_range', 'conflict_pdb_index',
                'raw_pdb_index', 'SEQRES_COUNT', 'ARTIFACT_INDEX', 'OBS_INDEX', 'NON_INDEX', 'OBS_RATIO_ARRAY']
        # tasks = full_df[cols].agg(self.bs_score_base, axis=1).tolist()
        tasks = tuple(map(self.bs_score_base, *tuple(full_df[col] for col in cols)))
        bs_score_df = DataFrame([await task for task in tasks],
            columns=['aligned_p', 'CN_terminal_p', 'insertion_p', 'deletion_p', 'mapped_out_range'])
        assert bs_score_df.shape[0] == full_df.shape[0]
        bs_score_df['bs_score'] = bs_score_df.aligned_p + bs_score_df.CN_terminal_p + bs_score_df.insertion_p + bs_score_df.deletion_p
        ret = concat([full_df, bs_score_df], axis=1)
        ret.bs_score = ret.bs_score/ret.c_unp_len
        ret = await self.generate_new_identity(ret)
        return ret

    @unsync
    async def pipe_score_for_pdb_entry(self, sifts_df, exp_cols):
        exp_df = DataFrame(await PDBs.fetch_exp_pipe(self), columns=exp_cols)
        summary_dict = await self.fetch_from_pdbe_api('api/pdb/entry/summary/', a_load_json, json=True)
        d = summary_dict[self.pdb_id][0]
        exp_df['revision_date'] = d['revision_date']
        exp_df['deposition_date'] = d['deposition_date']
        exp_df['release_date'] = d['release_date']

        pec_df, _ = await self.stats_chain()

        full_df = await self.pipe_score_base(sifts_df, pec_df)

        return full_df, exp_df

    @unsync
    async def pipe_score_for_unp_isoform(self, sifts_df, exp_cols):
        cur_pdbs = sifts_df.pdb_id.unique()

        pdbs = PDBs(cur_pdbs)
        # tasks = await pdbs.fetch(PDBs.fetch_exp_pipe).run()
        tasks = [await task for task in pdbs.fetch(PDBs.fetch_exp_pipe).tasks]
        exp_df = DataFrame((j for i in tasks for j in i), columns=exp_cols)
        r_len = exp_df.shape[0]
        # tasks = await pdbs.fetch(PDBs.fetch_date).run()
        tasks = [await task for task in pdbs.fetch(PDBs.fetch_date).tasks]
        exp_df = exp_df.merge(
            DataFrame(tasks, columns=['pdb_id', 'revision_date', 'deposition_date', 'release_date']))
        assert r_len == exp_df.shape[0]
        if self.entry_filter:
            exp_df = exp_df.query(self.entry_filter)
            if exp_df.shape[0] == 0:
                return
            cur_pdbs = exp_df.pdb_id.unique()
            sifts_df = sifts_df[sifts_df.pdb_id.isin(cur_pdbs)]
        # tasks = await PDBs(cur_pdbs).fetch('stats_chain').run()
        # tasks = [await task for task in PDBs(cur_pdbs).fetch('stats_chain').tasks]
        pec_df, _ = await PDBs(cur_pdbs).stats_chain()
        # pec_df, _ = zip(*tasks)
        # pec_df = concat(pec_df, sort=False, ignore_index=True)
        
        full_df = await self.pipe_score_base(sifts_df, pec_df)
        if self.chain_filter:
            full_df = full_df.query(self.chain_filter)
            if full_df.shape[0] == 0:
                return

        return full_df, exp_df

    @unsync
    async def pipe_select_base(self, exclude_pdbs=frozenset(), **kwargs):    
        res = await self.pipe_score(**kwargs)
        if res is None: return
        full_df, exp_df = res
        if self.source == 'UniProt':
            m_df = full_df[~full_df.pdb_id.isin(exclude_pdbs)]
            sele_df = merge(m_df, exp_df)
        else:
            sele_df = merge(full_df, exp_df, how='left')
        if len(sele_df) == 0:
            return
        else:
            assert sele_df.experimental_method_class.isnull().sum() == 0
        sele_df['1/resolution'] = 1 / sele_df.resolution
        sele_df['id_score'] = sele_df.chain_id.apply(id2score)
        return sele_df

    @staticmethod
    def select_mo(sele_df, OC_cutoff=0.2, sort_cols=['bs_score', '1/resolution', 'revision_date', 'id_score'], infer_new_col:bool=False, ascending=False, allow_mask=None):
        sele_df['select_tag'] = False
        sele_df['select_rank'] = -1
        if infer_new_col:
            for col in sort_cols:
                if (col not in sele_df.columns) and (col[1:] in sele_df.columns) and (col[0] == '-'):
                    sele_df[col] = -sele_df[col[1:]]
        sele_df = sele_df.sort_values(by=sort_cols, ascending=ascending).reset_index(drop=True)
        if allow_mask is True:
            allow_sele_df = sele_df
        elif allow_mask is None:
            allow_sele_df = sele_df[sele_df.bs_score > 0]
        else:
            allow_sele_df = sele_df[allow_mask]

        def sele_func(dfrm):
            rank_index = dfrm.index
            sele_df.loc[rank_index, 'select_rank'] = range(1, len(rank_index)+1)
            return select_range(dfrm.new_unp_range, rank_index, cutoff=OC_cutoff)

        if len(allow_sele_df.UniProt.unique()) > 1:
            sele_indexes = allow_sele_df.groupby('UniProt').apply(sele_func)
            sele_df.loc[[j for i in sele_indexes for j in i], 'select_tag'] = True
        else:
            sele_df.loc[sele_func(allow_sele_df), 'select_tag'] = True
        sele_df['select_rank_score'] = 1/sele_df.select_rank
        sele_df['after_select_rank'] = sele_df[['select_tag', 'select_rank_score']].apply(tuple, axis=1).rank(method='dense', ascending=False).astype(int)
        return sele_df.drop(columns=['select_rank_score'])

    @unsync
    async def pipe_select_mo(self, exclude_pdbs=frozenset(), OC_cutoff=0.2, select_mo_kwargs={}, **kwargs):        
        sele_df = await self.pipe_select_base(exclude_pdbs, **kwargs)
        if sele_df is None:
            return
        else:
            return self.select_mo(sele_df, OC_cutoff, **select_mo_kwargs)

    @unsync
    async def pipe_cluster3d_mo(self, exclude_pdbs=frozenset(), oc_cutoff: float = 0.8, rmsd_cutoff: float = 2.5, sort_cols=['bs_score', '1/resolution', 'revision_date', 'id_score'], infer_new_col: bool = False, ascending: bool = False, **kwargs):
        assert oc_cutoff >= 0 and oc_cutoff <= 1 and rmsd_cutoff >= 0
        sifts_df = await self.pipe_select_base(exclude_pdbs, **kwargs)
        if sifts_df is None:
            return
        else:
            sifts_df['oc_rmsd_select_tag'] = False
            sifts_df['oc_select_tag'] = False
            sifts_df['select_rank'] = -1
            if infer_new_col:
                for col in sort_cols:
                    if (col not in sifts_df.columns) and (col[1:] in sifts_df.columns) and (col[0] == '-'):
                        sifts_df[col] = -sifts_df[col[1:]]
            sifts_df = sifts_df.sort_values(by=sort_cols, ascending=ascending).reset_index(drop=True)
            sifts_df['select_rank'] = range(1, sifts_df.shape[0]+1)
        pdbekb_cluster = await self.fetch_from_pdbe_api('graph-api/uniprot/superposition/', Base.to_dataframe)
        if pdbekb_cluster is not None:
            pdbekb_cluster = pdbekb_cluster.filter(items=('UniProt', 'pdb_id', 'entity_id', 'struct_asym_id', 'pdbekb_cluster'))  # , 'is_representative'
            sifts_df = sifts_df.merge(pdbekb_cluster, how='left')
            fill_index = sifts_df[sifts_df.pdbekb_cluster.isnull()].index
            # NOTE: for cases like 2l53 2lqp
            sifts_df.loc[fill_index, 'pdbekb_cluster'] = list(f'multiple_model_{iclu}' for iclu in range(len(fill_index)))
            focus = sifts_df.loc[sifts_df.groupby('pdbekb_cluster').select_rank.idxmin().sort_values()].filter(items=('pdb_id', 'struct_asym_id', 'OBS_INDEX', 'OBS_COUNT', 'new_pdb_range_raw', 'new_unp_range_raw', 'select_rank', 'pdbekb_cluster'))
        else:
            sifts_df['pdbekb_cluster'] = list(f'none_assigned_{iclu}' for iclu in range(sifts_df.shape[0]))
            focus = sifts_df.filter(items=('pdb_id', 'struct_asym_id', 'OBS_INDEX', 'OBS_COUNT', 'new_pdb_range_raw', 'new_unp_range_raw', 'select_rank', 'pdbekb_cluster'))
        assert focus.shape[0] > 0
        focus[['new_pdb_range_obs', 'new_unp_range_obs']] = focus[['OBS_INDEX', 'new_pdb_range_raw', 'new_unp_range_raw']].apply(
            lambda x: overlap_range_2(*x), axis=1).apply(Series)
        cluster = await self.cluster3d_mo_algorithm(focus.drop(columns=['OBS_INDEX', 'new_pdb_range_raw', 'new_unp_range_raw']).itertuples(), oc_cutoff, rmsd_cutoff)
        selected_indexes = cluster.groupby(['oc_cluster', 'rmsd_cluster']).Index.min().sort_values()
        sifts_df.loc[selected_indexes, 'oc_rmsd_select_tag'] = True
        cluster = cluster.merge(focus[['pdbekb_cluster', 'new_unp_range_obs']].reset_index().rename(columns={'index': 'Index'}), how='left')
        assert not cluster.pdbekb_cluster.isnull().any()
        focus_oc_cluster = cluster.drop_duplicates(subset='oc_cluster', keep='first')
        selected_oc_indexes = cluster.loc[select_range(focus_oc_cluster.new_unp_range_obs, focus_oc_cluster.index, cutoff=1-oc_cutoff)].Index
        sifts_df.loc[selected_oc_indexes, 'oc_select_tag'] = True
        sifts_df = sifts_df.merge(cluster.drop(columns=['Index', 'new_unp_range_obs']), how='left')
        assert (not sifts_df.oc_cluster.isnull().any()) and (not sifts_df.rmsd_cluster.isnull().any())
        return sifts_df
    
    @staticmethod
    @unsync
    async def cluster3d_mo_algorithm(rows: Iterator, oc_cutoff: float, rmsd_cutoff: float) -> DataFrame:
        """
        # Monomer 3D Cluster (impl Greedy Algorithm)
        
        ## Data Structure

            [
                [
                    # one oc cluster
                    
                    [next(rows), ], # one rmsd cluster
                    
                    [],           # another rmsd cluster
                ],
                [
                    # another oc cluster
                ]
            ]
        """
        cluster = [[[next(rows)]]]
        for cur_row in rows:
            stop_sub_for_loop = False
            to_append_oc_index = -1
            to_append_rmsd_index = -1
            for oc_cluster_index, oc_cluster in enumerate(cluster):
                if stop_sub_for_loop: break
                for rmsd_cluster_index, rmsd_cluster in enumerate(oc_cluster):
                    if stop_sub_for_loop: break
                    for sele_row in rmsd_cluster:
                        overlap_unp_segs, overlap_sele_pdb_segs, overlap_cur_pdb_segs = overlap_range_3(
                            sele_row.new_unp_range_obs, 
                            cur_row.new_unp_range_obs,
                            sele_row.new_pdb_range_obs,
                            cur_row.new_pdb_range_obs)
                        oc_score = range_len(overlap_unp_segs) / min(sele_row.OBS_COUNT, cur_row.OBS_COUNT)
                        if oc_score >= oc_cutoff:
                            # same oc_cluster
                            sele_coords = await PDB(sele_row.pdb_id).get_representative_atoms(sele_row.struct_asym_id, overlap_pdb_segs=overlap_sele_pdb_segs)
                            cur_coords = await PDB(cur_row.pdb_id).get_representative_atoms(cur_row.struct_asym_id, overlap_pdb_segs=overlap_cur_pdb_segs)
                            if sele_coords.shape == cur_coords.shape:
                                rmsd_score = rmsd(sele_coords, cur_coords)
                            else:
                                #print('did not pass coords')
                                sele_coords_df = await PDB(sele_row.pdb_id).get_overlap_representative_atoms(sele_row.struct_asym_id, overlap_pdb_segs=overlap_sele_pdb_segs)
                                cur_coords_df = await PDB(cur_row.pdb_id).get_overlap_representative_atoms(cur_row.struct_asym_id, overlap_pdb_segs=overlap_cur_pdb_segs)
                                assert sele_coords_df.shape == cur_coords_df.shape
                                mask_1 = sele_coords_df['_atom_site.Cartn_x'].notnull()
                                mask_2 = cur_coords_df['_atom_site.Cartn_x'].notnull()
                                sele_coords = PDB.get_coordinates(sele_coords_df[mask_1 & mask_2])
                                cur_coords = PDB.get_coordinates(cur_coords_df[mask_1 & mask_2])
                                rmsd_score = rmsd(sele_coords, cur_coords)
                            if rmsd_score <= rmsd_cutoff:
                                # same oc_cluster & same rmsd_cluster
                                to_append_oc_index = oc_cluster_index
                                to_append_rmsd_index = rmsd_cluster_index
                                stop_sub_for_loop = True
                                break
                            else:
                                # same oc_cluster & possible different rmsd_cluster
                                to_append_oc_index = oc_cluster_index
                        else:
                            # possible different oc_cluster
                            pass
            # to append
            if to_append_oc_index < 0:
                cluster.append([[cur_row]])
            else:
                if to_append_rmsd_index < 0:
                    cluster[to_append_oc_index].append([cur_row])
                else:
                    cluster[to_append_oc_index][to_append_rmsd_index].append(cur_row)
        # to index
        return DataFrame((
            {'Index': sele_row.Index, 'oc_cluster': oc_cluster_index, 'rmsd_cluster': rmsd_cluster_index
            } for oc_cluster_index, oc_cluster in enumerate(cluster) for rmsd_cluster_index, rmsd_cluster in enumerate(oc_cluster) for sele_row in rmsd_cluster
        ))
                    

    @staticmethod
    @unsync
    async def search_partner_from_i3d(Entry, interaction_types, columns='pdb_id,Entry_1,Entry_2,assembly_id,model_id_1,chain_id_1,model_id_2,chain_id_2,organism,interaction_type'):
        if isinstance(interaction_types, str):
            query=f"""
                        SELECT {columns} FROM InteractionMeta
                        WHERE (Entry_1 = '{Entry}' or Entry_2 = '{Entry}') and interaction_type = '{interaction_types}'
            """
        elif isinstance(interaction_types, Tuple):
            if len(interaction_types) > 1:
                query=f"""
                        SELECT {columns} FROM InteractionMeta
                        WHERE (Entry_1 = '{Entry}' or Entry_2 = '{Entry}') and (interaction_type in {interaction_types})
                """
            elif len(interaction_types) == 1:
                query = f"""
                        SELECT {columns} FROM InteractionMeta
                        WHERE (Entry_1 = '{Entry}' or Entry_2 = '{Entry}') and interaction_type = '{interaction_types[0]}'
                """
            else:
                raise ValueError("Invalid interaction_types")
        else:
            raise ValueError("Invalid interaction_types")
            
        res = await Interactome3D.sqlite_api.database.fetch_all(query=query)
        return DataFrame(res, columns=columns.split(','))

    @staticmethod
    def parallel_interact_df(sifts_df, i3d_df, common_cols=('revision_date', 'deposition_date', 'release_date', 'experimental_method', 'experimental_method_class', 'multi_method', '-r_factor', '-r_free','resolution', '1/resolution')):
        rename_dict = dict(zip((f'{i}_1' for i in common_cols), common_cols))
        rename_dict['pdb_id_1'] = 'pdb_id'
        sifts_df_ = sifts_df.add_suffix('_1').rename(columns=rename_dict)
        i3d_df = i3d_df.merge(sifts_df_)
        sifts_df_ = sifts_df.drop(columns=sifts_df.columns.intersection(common_cols)).add_suffix('_2').rename(columns={'pdb_id_2': 'pdb_id'})
        i3d_df = i3d_df.merge(sifts_df_)
        swap_index = i3d_df[
            (i3d_df.struct_asym_id_1 > i3d_df.struct_asym_id_2) | 
            ((i3d_df.struct_asym_id_1 == i3d_df.struct_asym_id_2) & (i3d_df.model_id_1 > i3d_df.model_id_2))].index
        cols_1 = [col for col in i3d_df.columns if '_1' in col]
        cols_2 = [col for col in i3d_df.columns if '_2' in col]
        store_1 = i3d_df.loc[swap_index, cols_1].rename(columns=dict(zip(cols_1, [col.replace('_1', '_2') for col in cols_1])))
        store_2 = i3d_df.loc[swap_index, cols_2].rename(columns=dict(zip(cols_2, [col.replace('_2', '_1') for col in cols_2])))
        i3d_df.loc[swap_index, cols_1] = store_2
        i3d_df.loc[swap_index, cols_2] = store_1
        assert ((i3d_df.struct_asym_id_1 < i3d_df.struct_asym_id_2) | ((i3d_df.struct_asym_id_1 == i3d_df.struct_asym_id_2) & (i3d_df.model_id_1 < i3d_df.model_id_2))).all()
        return i3d_df

    '''
    @staticmethod
    async def pipe_interface_res_dict(p_df, pdb_id):
        chain_pairs = p_df[p_df.pdb_id.eq(pdb_id)][['struct_asym_id_1', 'struct_asym_id_2']].apply(frozenset, axis=1).to_list()
        pdb_ob = PDB(pdb_id)
        await pdb_ob.set_assembly(focus_assembly_ids=p_df[p_df.pdb_id.eq(pdb_id)].assembly_id.unique())
        for assembly in pdb_ob.assembly.values():
            try:
                await assembly.pipe_protein_protein_interface()
                for interface in assembly.interface.values():
                    if interface.info['chains'] in chain_pairs:
                        yield await interface.get_interface_res_dict()
            except Exception:
                pass
    '''

    @unsync
    async def pipe_select_ho_base(self, exclude_pdbs=frozenset(), run_as_completed: bool=False, progress_bar=None, skip_pdbs=None, select_mo_kwargs={}, skip_carbohydrate_polymer=False, **kwargs):
        sele_df = await self.pipe_select_mo(exclude_pdbs=exclude_pdbs, skip_pdbs=skip_pdbs, select_mo_kwargs=select_mo_kwargs, skip_carbohydrate_polymer=skip_carbohydrate_polymer)
        if sele_df is None:
            return
        chain_pairs = sele_df.groupby('pdb_id').struct_asym_id.apply(lambda x: frozenset(frozenset(item) for item in combinations_with_replacement(x, 2)))
        
        #p_df = await self.pisa_interact_integrate_with_i3d(sele_df, chain_pairs, 'ho', run_as_completed, progress_bar, **kwargs)
        p_df = await self.retrieve_rcsb_interface('ho', sele_df, chain_pairs=chain_pairs, **kwargs)
        if p_df is None:
            return None
        p_df['unp_range_DSC'] = p_df.apply(lambda x: sorensen.similarity(
            lyst2range(x['new_unp_range_1']), lyst2range(x['new_unp_range_2'])), axis=1)
        p_df['id_score_1'] = list(map(self.get_id_score_for_assembly, zip(p_df.struct_asym_id_1, p_df.asym_id_rank_1, p_df.assembly_id)))
        p_df['id_score_2'] = list(map(self.get_id_score_for_assembly, zip(p_df.struct_asym_id_2, p_df.asym_id_rank_2, p_df.assembly_id)))
        return self.add_interact_common_cols(p_df)

    @staticmethod
    def select_ho(p_df, interface_mapped_cov_cutoff=0.8, unp_range_DSC_cutoff=0.8, DSC_cutoff=0.2, sort_cols=['best_select_rank_score', 'second_select_rank_score'], ascending=False, allow_mask=None):
        p_df = p_df.sort_values(by=sort_cols, ascending=ascending).reset_index(drop=True)
        p_df['i_select_tag'] = False
        p_df['i_select_rank'] = -1
        if allow_mask is None:
            allow_p_df = p_df[
                (p_df.best_select_rank_score > 0) &
                (p_df.unp_range_DSC >= unp_range_DSC_cutoff) &
                ((p_df.unp_interface_range_1.apply(range_len)/p_df.interface_range_1.apply(range_len)) >= interface_mapped_cov_cutoff) &
                ((p_df.unp_interface_range_2.apply(range_len)/p_df.interface_range_2.apply(range_len)) >= interface_mapped_cov_cutoff)]
        elif allow_mask is True:
            allow_p_df = p_df[
                (p_df.unp_range_DSC >= unp_range_DSC_cutoff) &
                ((p_df.unp_interface_range_1.apply(range_len)/p_df.interface_range_1.apply(range_len)) >= interface_mapped_cov_cutoff) &
                ((p_df.unp_interface_range_2.apply(range_len)/p_df.interface_range_2.apply(range_len)) >= interface_mapped_cov_cutoff)]
        else:
            allow_p_df = p_df[
                allow_mask & 
                (p_df.unp_range_DSC >= unp_range_DSC_cutoff) &
                ((p_df.unp_interface_range_1.apply(range_len)/p_df.interface_range_1.apply(range_len)) >= interface_mapped_cov_cutoff) &
                ((p_df.unp_interface_range_2.apply(range_len)/p_df.interface_range_2.apply(range_len)) >= interface_mapped_cov_cutoff)]

        if allow_p_df.shape[0] == 0:
            return p_df

        def sele_func(dfrm):
            rank_index = dfrm.index
            p_df.loc[rank_index, 'i_select_rank'] = range(1, len(rank_index)+1)
            return select_ho_max_range(dfrm.unp_interface_range_1, dfrm.unp_interface_range_2, rank_index, cutoff=DSC_cutoff)

        if len(p_df.i_group.unique()) == 1:
            p_df.loc[sele_func(allow_p_df), 'i_select_tag'] = True
        else:
            sele_indexes = allow_p_df.groupby('i_group').apply(sele_func)
            p_df.loc[[j for i in sele_indexes for j in i], 'i_select_tag'] = True

        return p_df

    @unsync
    async def pipe_select_ho(self, interface_mapped_cov_cutoff=0.8, unp_range_DSC_cutoff=0.8, DSC_cutoff=0.2, **kwargs):
        p_df = await self.pipe_select_ho_base(**kwargs)
        if p_df is None:
            return
        else:
            return self.select_ho(p_df, interface_mapped_cov_cutoff, unp_range_DSC_cutoff, DSC_cutoff)

    @unsync
    async def pipe_select_ho_iso_base(self, exclude_pdbs=frozenset(), run_as_completed:bool=False, progress_bar=None, skip_pdbs=None, select_mo_kwargs={}, skip_carbohydrate_polymer=False, **kwargs):
        sele_df = await self.pipe_select_mo(exclude_pdbs=exclude_pdbs, complete_chains=True, skip_pdbs=skip_pdbs, select_mo_kwargs=select_mo_kwargs, skip_carbohydrate_polymer=skip_carbohydrate_polymer)
        if sele_df is None:
            return
        sele_df = sele_df[sele_df.Entry.eq(self.get_id().split('-')[0])]
        chain_pairs = sele_df.groupby('pdb_id').struct_asym_id.apply(lambda x: frozenset(frozenset(item) for item in combinations_with_replacement(x, 2)))
        
        #p_df = await self.pisa_interact_integrate_with_i3d(sele_df, chain_pairs, 'ho_iso', run_as_completed, progress_bar, **kwargs)
        p_df = await self.retrieve_rcsb_interface('ho_iso', sele_df, chain_pairs=chain_pairs, **kwargs)
        if p_df is None:
            return
        else:
            p_df['id_score_1'] = list(map(self.get_id_score_for_assembly, zip(p_df.struct_asym_id_1, p_df.asym_id_rank_1, p_df.assembly_id)))
            p_df['id_score_2'] = list(map(self.get_id_score_for_assembly, zip(p_df.struct_asym_id_2, p_df.asym_id_rank_2, p_df.assembly_id)))
            return self.add_interact_common_cols(p_df)
    
    @classmethod
    def select_ho_iso(cls, p_df, interface_mapped_cov_cutoff=0.8, DSC_cutoff=0.2, sort_cols=['bs_score_1', 'bs_score_2', '1/resolution', 'revision_date', 'id_score_1', 'id_score_2'], ascending=False, allow_mask=None):
        return cls.select_he(p_df, interface_mapped_cov_cutoff, DSC_cutoff, sort_cols, ascending, allow_mask)

    @unsync
    async def pipe_select_ho_iso(self, interface_mapped_cov_cutoff=0.8, DSC_cutoff=0.2, **kwargs):
        p_df = await self.pipe_select_ho_iso_base(**kwargs)
        if p_df is None:
            return
        else:
            return self.select_ho_iso(self.sort_interact_cols(p_df), interface_mapped_cov_cutoff, DSC_cutoff)   

    @unsync
    async def pipe_select_else_base(self, func:str = "Protein/NA", exclude_pdbs=frozenset(), skip_pdbs=None, select_mo_kwargs={}, skip_carbohydrate_polymer=False, **kwargs):
        assert func != 'pipe_protein_protein_interface'
        sele_df = await self.pipe_select_mo(exclude_pdbs=exclude_pdbs, skip_pdbs=skip_pdbs, select_mo_kwargs=select_mo_kwargs, skip_carbohydrate_polymer=skip_carbohydrate_polymer)
        if sele_df is None:
            return
        include_chains = sele_df.groupby('pdb_id').struct_asym_id.apply(frozenset)
        #return await self.pisa_interact_protein_else(sele_df, include_chains, func, run_as_completed, progress_bar, **kwargs)
        return await self.retrieve_rcsb_interface(func, sele_df, include_chains=include_chains, **kwargs)

    @unsync
    async def pipe_select_he_base(self, exclude_pdbs=frozenset(), run_as_completed:bool=False, progress_bar=None, skip_pdbs=None, select_mo_kwargs={}, skip_carbohydrate_polymer=False, **kwargs):
        sele_df = await self.pipe_select_mo(exclude_pdbs=exclude_pdbs, complete_chains=True, skip_pdbs=skip_pdbs, select_mo_kwargs=select_mo_kwargs, skip_carbohydrate_polymer=skip_carbohydrate_polymer)
        if sele_df is None:
            return
        if len(sele_df.Entry.unique()) == 1:
            return
        chain_pairs = sele_df.groupby(['pdb_id', 'entity_id']).struct_asym_id.apply(frozenset).groupby('pdb_id').apply(tuple).apply(lambda x: frozenset(res for i,j in combinations(range(len(x)), 2) for res in product(x[i], x[j])))
        chain_pairs = chain_pairs[chain_pairs.apply(len) > 0]
        chain_pairs = chain_pairs.apply(lambda x: frozenset(frozenset(item) for item in x))
        
        #p_df = await self.pisa_interact_integrate_with_i3d(sele_df, chain_pairs, 'he', run_as_completed, progress_bar, **kwargs)
        p_df = await self.retrieve_rcsb_interface('he', sele_df, chain_pairs=chain_pairs, **kwargs)
        if p_df is None:
            return
        else:
            p_df['id_score_1'] = list(map(self.get_id_score_for_assembly, zip(p_df.struct_asym_id_1, p_df.asym_id_rank_1, p_df.assembly_id)))
            p_df['id_score_2'] = list(map(self.get_id_score_for_assembly, zip(p_df.struct_asym_id_2, p_df.asym_id_rank_2, p_df.assembly_id)))
            return self.add_interact_common_cols(p_df)

    @staticmethod
    def select_else(p_df, interface_mapped_cov_cutoff=0.8, DSC_cutoff=0.2, sort_cols=['bs_score_1', '1/resolution', 'revision_date', 'id_score_1'], ascending=False, allow_mask=None):
        p_df = p_df.sort_values(by=sort_cols, ascending=ascending).reset_index(drop=True)
        p_df['i_select_tag'] = False
        p_df['i_select_rank'] = -1
        if allow_mask is None:
            allow_mask = (p_df.bs_score_1 > 0)
        elif allow_mask is True:
            allow_mask = p_df.bs_score_1.notnull()
        mask = allow_mask & (
            (p_df.unp_interface_range_1.apply(range_len)/p_df.interface_range_1.apply(range_len)) >= interface_mapped_cov_cutoff)
        if not mask.any():
            return p_df
        allow_p_df = p_df[mask]
        
        def sele_func(dfrm):
            rank_index = dfrm.index
            p_df.loc[rank_index, 'i_select_rank'] = range(1, len(rank_index)+1)
            return select_range(dfrm.unp_interface_range_1, rank_index, cutoff=DSC_cutoff, similarity_func=sorensen.similarity)
        
        sele_indexes = allow_p_df.groupby('UniProt_1').apply(sele_func)
        try:
            p_df.loc[[j for i in sele_indexes for j in i], 'i_select_tag'] = True
        except Exception:
            raise ValueError("{sele_indexes}\n{p_df}")
        return p_df

    @staticmethod
    def select_he(p_df, interface_mapped_cov_cutoff=0.8, DSC_cutoff=0.2, sort_cols=['bs_score_1', 'bs_score_2', '1/resolution', 'revision_date', 'id_score_1', 'id_score_2'], ascending=False, allow_mask=None):
        p_df = p_df.sort_values(by=['i_group']+sort_cols, ascending=ascending).reset_index(drop=True)
        # Groupby preserves the order of rows within each group. <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.groupby.html>
        p_df['i_select_tag'] = False
        p_df['i_select_rank'] = -1
        
        if allow_mask is None:
            allow_p_df = p_df[
                (p_df.best_select_rank_score > 0) &
                ((p_df.unp_interface_range_1.apply(range_len)/p_df.interface_range_1.apply(range_len)) >= interface_mapped_cov_cutoff) &
                ((p_df.unp_interface_range_2.apply(range_len)/p_df.interface_range_2.apply(range_len)) >= interface_mapped_cov_cutoff)]
        elif allow_mask is True:
            allow_p_df = p_df[
                ((p_df.unp_interface_range_1.apply(range_len)/p_df.interface_range_1.apply(range_len)) >= interface_mapped_cov_cutoff) &
                ((p_df.unp_interface_range_2.apply(range_len)/p_df.interface_range_2.apply(range_len)) >= interface_mapped_cov_cutoff)]
        else:
            allow_p_df = p_df[
                allow_mask &
                ((p_df.unp_interface_range_1.apply(range_len)/p_df.interface_range_1.apply(range_len)) >= interface_mapped_cov_cutoff) &
                ((p_df.unp_interface_range_2.apply(range_len)/p_df.interface_range_2.apply(range_len)) >= interface_mapped_cov_cutoff)]

        if allow_p_df.shape[0] == 0:
            return p_df

        def sele_func(dfrm):
            rank_index = dfrm.index
            p_df.loc[rank_index, 'i_select_rank'] = range(1, len(rank_index)+1)
            return select_he_range(dfrm.UniProt_1, dfrm.UniProt_2, dfrm.unp_interface_range_1, dfrm.unp_interface_range_2, rank_index, cutoff=DSC_cutoff)

        sele_indexes = allow_p_df.groupby('i_group').apply(sele_func)
        p_df.loc[[j for i in sele_indexes for j in i], 'i_select_tag'] = True
        return p_df

    @unsync
    async def pipe_select_he(self, interface_mapped_cov_cutoff=0.8, DSC_cutoff=0.2, **kwargs):
        p_df = await self.pipe_select_he_base(**kwargs)
        if p_df is None:
            return
        else:
            return self.select_he(self.sort_interact_cols(p_df), interface_mapped_cov_cutoff, DSC_cutoff)

    @unsync
    async def pipe_select_else(self, interface_mapped_cov_cutoff=0.8, DSC_cutoff=0.2, **kwargs):
        p_df = await self.pipe_select_else_base(**kwargs)
        if p_df is None:
            return
        else:
            return self.select_else(p_df, interface_mapped_cov_cutoff, DSC_cutoff)

    #@unsync
    def sort_interact_cols(self, dfrm):
        assert self.source == 'UniProt'
        if isinstance(dfrm, Unfuture):
            dfrm = dfrm.result()
        swap_index = dfrm[dfrm.UniProt_1.ne(self.get_id())].index
        cols_1 = [col for col in dfrm.columns if '_1' in col]
        cols_2 = [col for col in dfrm.columns if '_2' in col]
        store_1 = dfrm.loc[swap_index, cols_1].rename(columns=dict(zip(cols_1, [col.replace('_1', '_2') for col in cols_1])))
        store_2 = dfrm.loc[swap_index, cols_2].rename(columns=dict(zip(cols_2, [col.replace('_2', '_1') for col in cols_2])))
        dfrm.loc[swap_index, cols_1] = store_2
        dfrm.loc[swap_index, cols_2] = store_1
        dfrm['i_group'] = dfrm.apply(lambda x: "('%s','%s')" % (x['UniProt_1'], x['UniProt_2']), axis=1)
        return dfrm

    @staticmethod
    def select_smr_mo(smr_df, allow_oligo_state=None, selected_sifts_unp_ranges=list(), smr_sort_cols=None, ascending=False, OC_cutoff=0.2):
        smr_df['select_rank'] = -1
        smr_df['select_tag'] = False
        if allow_oligo_state is not None:
            allow_smr_df = smr_df[smr_df['oligo-state'].isin(allow_oligo_state)] # allow_oligo_state=('monomer',)
        else:
            allow_smr_df = smr_df
        
        if smr_sort_cols is None:
            rank_index = allow_smr_df.index
        else:
            rank_index = allow_smr_df.sort_values(
                by=smr_sort_cols, ascending=ascending).index
        
        smr_df.loc[rank_index, 'select_rank'] = range(1, len(rank_index)+1)
        smr_df.loc[select_range(smr_df.unp_range, rank_index, cutoff=OC_cutoff, selected_ranges=list(selected_sifts_unp_ranges)), 'select_tag'] = True
        return smr_df
    
    @unsync
    async def pipe_select_smr_mo(self, smr_df=None, **kwargs):
        if 'sifts_mo_df' in kwargs:
            sifts_mo_df = kwargs['sifts_mo_df']
        else:
            sifts_mo_df = await self.pipe_select_mo(**kwargs)
        smr_df = (await SMR.single_retrieve(self.get_id()).then(SMR.to_dataframe)) if smr_df is None else smr_df
        if smr_df is None or len(smr_df) == 0:
            return
        return self.select_smr_mo(
            smr_df, 
            kwargs.get('allow_oligo_state', None),
            sifts_mo_df[sifts_mo_df.select_tag.eq(True)].new_unp_range if sifts_mo_df is not None else [], 
            kwargs.get('smr_sort_cols', None),
            kwargs.get('ascending', False),
            kwargs.get('OC_cutoff', 0.2))

    @classmethod
    def add_interact_common_cols(cls, p_df):
        p_df['best_select_rank_score'] = p_df[['select_rank_1',
                                               'select_rank_2']].apply(lambda x: 1/min(x), axis=1)
        p_df['second_select_rank_score'] = p_df[['select_rank_1',
                                                 'select_rank_2']].apply(lambda x: 1/max(x), axis=1)
        p_df['unp_interface_range_1'] = p_df.apply(lambda x: to_interval(cls.convert_index(x['new_unp_range_1'], x['new_pdb_range_1'], expand_interval(x['interface_range_1']))), axis=1)
        p_df['unp_interface_range_2'] = p_df.apply(lambda x: to_interval(cls.convert_index(x['new_unp_range_2'], x['new_pdb_range_2'], expand_interval(x['interface_range_2']))), axis=1)
        p_df['i_group'] = p_df.apply(lambda x: "('%s','%s')" % tuple(sorted((x['UniProt_1'], x['UniProt_2']))), axis=1)
        return p_df

    @unsync
    async def pipe_scheduled_interface_res_dict_for_pdb(self, **kwargs):
        assert self.source == 'PDB'
        return await self.schedule_interface_tasks((await self.pipe_interface_res_dict(**kwargs)),False,None)

    @staticmethod
    async def schedule_interface_tasks(ob, run_as_completed, progress_bar):
        if run_as_completed:
            assert isinstance(ob, (SIFTSs,PDBs))
            res = await ob.run(tqdm=progress_bar)
            ob.tasks = [i.get_interface_res_dict() for interfaces in res for i in interfaces]
            res = await ob.run(tqdm=progress_bar)
        else:
            if isinstance(ob, (SIFTSs,PDBs)):
                res = [await i for i in ob.tasks]
                inteface_lyst = [i for interfaces in res for i in interfaces if not i.use_au]  # TODO: check whether 'not use_au' will affect selected results
            else:
                assert isinstance(ob, list)
                inteface_lyst = [i for i in ob if not i.use_au]
            res = []
            for index in range(0, len(inteface_lyst), 100):
                tasks = [i.get_interface_res_dict() for i in inteface_lyst[index:index+100]]
                res.extend([await i for i in tasks])
        return DataFrame(j for j in res if j is not None)

    @unsync
    async def pipe_scheduled_ranged_map_res_df(self, chunksize=100, func_for_this='pipe_select_mo', default_mask=True, with_sele_cols_for_unp=False, **kwargs):
        if self.source == 'UniProt':
            df = await getattr(self, func_for_this)(**kwargs)
            if default_mask:
                records = df[df.select_rank.ne(-1)].to_records()
            else:
                records = df.to_records()
        elif self.source == 'PDB':
            df = await self.pipe_select_base(**kwargs)
            records = df.to_records()
        res = []
        for index in range(0, len(records), chunksize):
            tasks = [PDB(row.pdb_id).get_ranged_map_res_df(
                        UniProt=row.UniProt,
                        unp_range=row.new_unp_range_raw,
                        pdb_range=row.new_pdb_range_raw,
                        conflict_pdb_index=row.conflict_pdb_index,
                        struct_asym_id=row.struct_asym_id) for row in records[index:index+chunksize]]
            res.extend([await i for i in tasks])
        ret = concat(res, sort=False, ignore_index=True)
        if with_sele_cols_for_unp and (self.source == 'UniProt') and (func_for_this == 'pipe_select_mo'):
            ret = ret.merge(df[['UniProt','pdb_id','struct_asym_id','bs_score','select_tag','select_rank','after_select_rank']], how='left')
            assert not ret.select_tag.isnull().any()
        return ret

    @staticmethod
    def get_id_score_for_assembly(args):
        struct_asym_id, asym_id_rank, assembly_id = args
        return id2score(struct_asym_id) - asym_id_rank - assembly_id

    @unsync
    async def pisa_interact_protein_else(self, sele_df, include_chains, func:str, run_as_completed:bool=False, progress_bar=None, **kwargs):
        # TODO: check
        ob = PDBs(include_chains.index).fetch('pipe_interface_res_dict_ic', include_chains=include_chains, func=func, **kwargs)
        interact_df = await self.schedule_interface_tasks(ob, run_as_completed, progress_bar)
        if len(interact_df) == 0:
            return
        await self.sqlite_api.async_insert(self.sqlite_api.PISAInterfaceDict, interact_df.to_dict('records'))
        for col in ('surface_range_2', 'interface_range_2'):
            if col not in interact_df.columns:
                interact_df[col] = nan
        check_mask = interact_df.molecule_type_1.isin(('polypeptide(L)', 'polypeptide(D)'))
        if not check_mask.all():
            # EXAMPLE: 5b0y/0/78
            warn('Outdated PISA chain identifier! Current data could be ligand related: ' +
                 str(interact_df[~check_mask].head(1).to_dict('records')[0]), PISAErrorWarning)
        interact_df = interact_df[check_mask & interact_df.interface_range_1.notnull()]
        if len(interact_df) == 0:
            return
        common_cols = ('pdb_id', 'revision_date', 'deposition_date', 'release_date', 'experimental_method', 'experimental_method_class', 'multi_method', '-r_factor', '-r_free','resolution', '1/resolution')
        rename_dict = dict(zip((f'{i}_1' for i in common_cols), common_cols))
        sele_df = sele_df.add_suffix('_1').rename(columns=rename_dict)
        p_df = interact_df.merge(sele_df)
        same_cols = interact_df.columns.intersection(sele_df.columns)
        assert p_df.shape[0] > 0, f"{interact_df[same_cols]}\n{sele_df[same_cols]}"
        p_df['unp_interface_range_1'] = p_df.apply(lambda x: to_interval(self.convert_index(x['new_unp_range_1'], x['new_pdb_range_1'], expand_interval(x['interface_range_1']))), axis=1)
        p_df['id_score_1'] = list(map(self.get_id_score_for_assembly, zip(p_df.struct_asym_id_1, p_df.asym_id_rank_1, p_df.assembly_id)))
        return p_df

    @unsync
    async def retrieve_rcsb_interface(self, interaction_type, sele_df, chain_pairs=None, include_chains=None, **kwargs):
        if chain_pairs is not None:
            ob = PDBs(chain_pairs.index).fetch('pipe_interface_res_info_for_rcsb', chain_pairs=chain_pairs, include_chains=include_chains, **kwargs)
        else:
            ob = PDBs(include_chains.index).fetch('pipe_interface_res_info_for_rcsb', chain_pairs=chain_pairs, include_chains=include_chains, **kwargs)
        interfaces_dfs = [await i for i in ob.tasks]
        interfaces_dfs = [i for i in interfaces_dfs if i is not None]
        if interfaces_dfs:
            interfaces_dfs = concat(interfaces_dfs, axis=0, ignore_index=True, sort=False)
            assert interfaces_dfs.shape[0] > 0
            if interaction_type in ("ho", "he", "ho_iso"):
                p_a_df = self.parallel_interact_df(sele_df, interfaces_dfs)
                if interaction_type == 'ho':
                    assert len(p_a_df) > 0
                elif interaction_type == 'ho_iso':
                    p_a_df = p_a_df[(p_a_df.UniProt_1.eq(self.get_id()) | p_a_df.UniProt_2.eq(self.get_id())) & (p_a_df.UniProt_1 != p_a_df.UniProt_2) & (p_a_df.struct_asym_id_in_assembly_1 != p_a_df.struct_asym_id_in_assembly_2)].reset_index(drop=True)
                elif interaction_type == 'he':
                    p_a_df = p_a_df[(p_a_df.UniProt_1.eq(self.get_id()) | p_a_df.UniProt_2.eq(self.get_id())) & (p_a_df.Entry_1 != p_a_df.Entry_2)].reset_index(drop=True)
            else:
                assert (chain_pairs is None) and (include_chains is not None)
                interfaces_dfs = interfaces_dfs[interfaces_dfs.polymer_composition.eq(interaction_type)] # e.g. "Protein/NA"
                if len(interfaces_dfs) == 0:
                    return
                
                common_cols = ('pdb_id', 'revision_date', 'deposition_date', 'release_date', 'experimental_method', 'experimental_method_class', 'multi_method', '-r_factor', '-r_free','resolution', '1/resolution')
                rename_dict = dict(zip((f'{i}_1' for i in common_cols), common_cols))
                sele_df = sele_df.add_suffix('_1').rename(columns=rename_dict)
                p_a_df = interfaces_dfs.merge(sele_df)
                same_cols = interfaces_dfs.columns.intersection(sele_df.columns)
                assert p_a_df.shape[0] > 0, f"{interfaces_dfs[same_cols]}\n{sele_df[same_cols]}"
                p_a_df['unp_interface_range_1'] = p_a_df.apply(lambda x: to_interval(self.convert_index(x['new_unp_range_1'], x['new_pdb_range_1'], expand_interval(x['interface_range_1']))), axis=1)
                p_a_df['id_score_1'] = list(map(self.get_id_score_for_assembly, zip(p_a_df.struct_asym_id_1, p_a_df.asym_id_rank_1, p_a_df.assembly_id)))

            if len(p_a_df) == 0:
                return
            return p_a_df
        else:
            return

    @unsync
    async def pisa_interact_integrate_with_i3d(self, sele_df, chain_pairs, interaction_type:str, run_as_completed:bool=False, progress_bar=None, **kwargs):
        ob = PDBs(chain_pairs.index).fetch('pipe_interface_res_dict', chain_pairs=chain_pairs, func='pipe_protein_protein_interface', **kwargs)
        interact_df = await self.schedule_interface_tasks(ob, run_as_completed, progress_bar)
        if len(interact_df) == 0:
            return
        await self.sqlite_api.async_insert(self.sqlite_api.PISAInterfaceDict, interact_df.to_dict('records'))
        '''
        NOTE: exception case: 5b0y (see api%pisa%interfacelist%+5b0y%0.tsv)

        # assert all((~interact_df.interface_range_1.isnull()) & (~interact_df.interface_range_2.isnull())), f"{interact_df[interact_df.interface_range_1.isnull() | interact_df.interface_range_2.isnull()]}"
        '''
        interact_df = interact_df[(~interact_df.interface_range_1.isnull()) & (~interact_df.interface_range_2.isnull())]
        if len(interact_df) == 0:
            return
        p_a_df = self.parallel_interact_df(sele_df, interact_df)
        if interaction_type == 'ho':
            assert len(p_a_df) > 0
        elif interaction_type == 'ho_iso':
            p_a_df = p_a_df[(p_a_df.UniProt_1.eq(self.get_id()) | p_a_df.UniProt_2.eq(self.get_id())) & (p_a_df.UniProt_1 != p_a_df.UniProt_2) & (p_a_df.struct_asym_id_in_assembly_1 != p_a_df.struct_asym_id_in_assembly_2)].reset_index(drop=True)
        elif interaction_type == 'he':
            p_a_df = p_a_df[(p_a_df.UniProt_1.eq(self.get_id()) | p_a_df.UniProt_2.eq(self.get_id())) & (p_a_df.Entry_1 != p_a_df.Entry_2)].reset_index(drop=True)
        else:
            raise ValueError(f"Invalid interaction_type: {interaction_type}!")
        if len(p_a_df) == 0:
            return
        i3d_df = await self.search_partner_from_i3d(self.get_id().split('-')[0], interaction_type[:2])

        if len(i3d_df) == 0:
            p_df = p_a_df
            p_df['in_i3d'] = False
        else:
            p_b_df = self.parallel_interact_df(sele_df[['UniProt', 'pdb_id', 'chain_id', 'struct_asym_id']], i3d_df)
            p_df = p_a_df.merge(p_b_df, how='left')
            p_df['in_i3d'] = p_df.organism.apply(lambda x: False if isna(x) else True)
            p_df.drop(columns=['organism', 'interaction_type'], inplace=True)
        return p_df
    
    @unsync
    async def unp_is_canonical(self):
        if self.source != 'UniProt':
            return None
        if '-' not in self.get_id():
            return True
        try:
            header = (await self.fetch_unp_fasta(self.get_id()))[0]
        except TypeError:
            warn(self.get_id(), PossibleObsoletedUniProtWarning)
            return None
        return self.get_id() != self.unp_head.match(header).group(1)
    
    @unsync
    async def unp_is_canonical_with_id(self):
        return self.get_id(), (await self.unp_is_canonical())

    @classmethod
    @unsync
    def meta_sifts_annotation(cls, path):
        df = cls.to_dataframe(path).result()
        if df is None or len(df) == 0:
            return None
        df = df.copy()
        df.start = df.start.apply(json.loads).apply(itemgetter('residue_number'))
        df.end = df.end.apply(json.loads).apply(itemgetter('residue_number'))
        if 'secondary_structure' in df.columns:
            resource_col = 'secondary_structure'
            df['secondary_structure'] = df.apply(lambda x: '{}_{}'.format(x['secondary_structure'], x['sheet_id'] if x['secondary_structure'] == 'strands' else x['start']), axis=1)
        else:
            resource_col = df.columns[0]
        df['resource'] = resource_col
        return df[['pdb_id', 'entity_id', 'struct_asym_id', 'chain_id', 'resource', resource_col, 'start', 'end']].rename(columns={
            resource_col: 'resource_id',
            'start': 'pdb_beg',
            'end': 'pdb_end'})

    @staticmethod
    async def get_mapped_pdbekb_annotaions_task_unit(pdb_ob, record, res_df):
        sites_df = res_df[res_df.struct_asym_id.eq(record.struct_asym_id)]
        if sites_df.shape[0] == 0:
            return None
        return await pdb_ob.get_map_res_df(
            record.UniProt,
            record.new_unp_range,
            record.new_pdb_range,
            your_sites=sites_df,
            conflict_pdb_index=record.conflict_pdb_index,
            unp2pdb=False,
            author_site=True,
            already_merged_with_residue_listing=True,
            struct_asym_id=record.struct_asym_id)

    @classmethod
    @unsync
    async def get_mapped_pdbekb_annotaions_task(cls, pdb_id, sub_sifts_df, api_suffix, **kwargs):
        rets = []
        pdb_ob = PDB(pdb_id)
        res_df = await pdb_ob.pipe_pdbekb_annotations(api_suffix, **kwargs)
        if res_df is None:
            return rets
        for record in sub_sifts_df.to_records():
            res = await cls.get_mapped_pdbekb_annotaions_task_unit(pdb_ob, record, res_df)
            if res is not None:
                rets.append(res)
        return rets

    @unsync
    async def get_mapped_pdbekb_annotaions(self, api_suffix, query_filter:str='bs_score > 0', pipe_select_mo_kwargs={}, **kwargs):
        sifts_df = await self.pipe_select_mo(**pipe_select_mo_kwargs)
        if sifts_df is None:
            return
        if query_filter != '' and query_filter is not None:
            sifts_df = sifts_df.query(query_filter)
            if sifts_df.shape[0] == 0:
                return None
        tasks = [self.get_mapped_pdbekb_annotaions_task(pdb_id, sub_sifts_df, api_suffix, **kwargs) for pdb_id, sub_sifts_df in sifts_df.groupby('pdb_id')]
        dfs = []
        for task in tasks:
            dfs.extend(await task)
        if len(dfs) == 0:
            return None
        else:
            return concat(dfs, sort=False, ignore_index=True)

    sifts_version_pat = re_compile(r'# (.+) \| PDB: (.+) \| UniProt: (.+)')

    @classmethod
    def check_version(cls):
        sifts_df = read_csv(
            FetchBase('EBI', 'pub/databases/msd/sifts/flatfiles/tsv/uniprot_pdb.tsv.gz').download(
            cls.web_semaphore, cls.get_folder()/'for_version_control/uniprot_pdb.tsv.gz'
        ).result(), sep='\t')
        # 2021/08/22 - 03:28 | PDB: 33.21 | UniProt: 2021.03
        sifts_version, sifts_pdb_version, sifts_unp_version = cls.sifts_version_pat.fullmatch(sifts_df.columns[0]).groups()
        return sifts_df


class Compounds(Base):

    tasks = LRUCache(maxsize=1024)

    def __post_init__(self, **kwargs):
        super().__post_init__(**kwargs)
        assert self.level == 'compounds'
        self.check_folder()
        self.raw_identifier = self.raw_identifier.upper()


class PDBs(tuple):
    '''immutable iterable class (tuple-like)'''

    def __new__(cls, iterable:Iterable=tuple()):
        return super(PDBs, cls).__new__(cls, (PDB(i) if isinstance(i, str) else i for i in iterable))
    
    def __getitem__(self, slice):
        res = tuple.__getitem__(self, slice)
        if isinstance(res, Iterable):
            return self.__class__(res)
        else:
            return res

    @staticmethod
    @unsync
    async def fetch_exp_pipe(pdb_ob:PDB):
        exp_data = await pdb_ob.fetch_from_pdbe_api('api/pdb/entry/experiment/', a_load_json, json=True)
        assert exp_data is not None, f"{pdb_ob}: no data for 'api/pdb/entry/experiment/'"
        methods = set(lyst['experimental_method'] for lyst in exp_data[pdb_ob.pdb_id])
        multi_method = len(methods) > 1
        ret = []
        for i in exp_data[pdb_ob.pdb_id]:
            if i['experimental_method'] in methods:
                methods.remove(i['experimental_method'])
                ret.append((
                    pdb_ob.pdb_id, 
                    i.get('resolution', None), 
                    i['experimental_method_class'], 
                    i['experimental_method'],
                    multi_method,
                    -i.get('r_factor', nan) if i.get('r_factor', nan) is not None else nan,
                    -i.get('r_free', nan) if i.get('r_free', nan) is not None else nan
                ))
        return ret
    
    @staticmethod
    @unsync
    async def fetch_number_of_entities(key: str, pdb_ob: PDB):
        summary_data = await pdb_ob.fetch_from_pdbe_api('api/pdb/entry/summary/', a_load_json, json=True)
        if key:
            try:
                return pdb_ob.pdb_id, summary_data[pdb_ob.pdb_id][0]['number_of_entities'][key]
            except KeyError:
                warn(f"{repr(pdb_ob)}: without '{key}' key")
                return pdb_ob.pdb_id, 0
        else:
            return pdb_ob.pdb_id, summary_data[pdb_ob.pdb_id][0]['number_of_entities']

    @staticmethod
    @unsync
    async def fetch_date(pdb_ob: PDB):
        summary_data = await pdb_ob.fetch_from_pdbe_api('api/pdb/entry/summary/', a_load_json, json=True)
        data = summary_data[pdb_ob.pdb_id][0]
        return pdb_ob.pdb_id, data['revision_date'], data['deposition_date'], data['release_date']
    
    def fetch(self, fetch_func:Union[Callable[[PDB], List], str], **kwargs):
        if isinstance(fetch_func, str):
            self.tasks = [getattr(pdb_ob, fetch_func)(**kwargs) for pdb_ob in self]
        else:
            self.tasks = [fetch_func(pdb_ob, **kwargs) for pdb_ob in self]
        return self

    @unsync
    async def run(self, tqdm=None):
        if tqdm is None:
            return [await fob for fob in as_completed(self.tasks)]
        else:
            return [await fob for fob in tqdm(as_completed(self.tasks), total=len(self.tasks))]

    @unsync
    async def stats_protein_entity_seq(self):
        all_pdbs = frozenset(pdb_ob.pdb_id for pdb_ob in self)
        stored = await PDB.sqlite_api.StatsProteinEntitySeq.objects.filter(pdb_id__in=all_pdbs).all()
        ed_pdbs = frozenset(i.pdb_id for i in stored)
        new = [await task for task in PDBs(all_pdbs-ed_pdbs).fetch('stats_protein_entity_seq').tasks]
        return stored + [j._asdict() if isinstance(j, PDB.StatsProteinEntitySeq) else j for i in new for j in i]

    @unsync
    async def stats_nucleotide_entity_seq(self):
        all_pdbs = frozenset(pdb_ob.pdb_id for pdb_ob in self)
        stored = await PDB.sqlite_api.StatsNucleotideEntitySeq.objects.filter(pdb_id__in=all_pdbs).all()
        ed_pdbs = frozenset(i.pdb_id for i in stored)
        new = [await task for task in PDBs(all_pdbs-ed_pdbs).fetch('stats_nucleotide_entity_seq').tasks]
        return stored + [j._asdict() if isinstance(j, PDB.StatsNucleotideEntitySeq) else j for i in new for j in i]

    @unsync
    async def stats_chain_seq(self):
        all_pdbs = frozenset(pdb_ob.pdb_id for pdb_ob in self)
        stored = await PDB.sqlite_api.StatsChainSeq.objects.filter(pdb_id__in=all_pdbs).all()
        ed_pdbs = frozenset(i.pdb_id for i in stored)
        new = [await task for task in PDBs(all_pdbs-ed_pdbs).fetch('stats_chain_seq').tasks]
        return stored + [j._asdict() if isinstance(j, PDB.StatsChainSeq) else j for i in new for j in i]

    @unsync
    async def stats_chain(self, stats_nucleotide=False):
        pe_df = DataFrame(await self.stats_protein_entity_seq())
        c_df = DataFrame(await self.stats_chain_seq())
        try:
            pec_df = c_df.merge(pe_df)
            pec_df['OBS_STD_INDEX'] = pec_df.apply(lambda x: overlap_range(x['OBS_INDEX'], x['STD_INDEX']), axis=1)
            pec_df['OBS_STD_COUNT'] = pec_df.OBS_STD_INDEX.apply(range_len)
        except Exception:
            pec_df = None
        if stats_nucleotide:
            ne_df = DataFrame(await self.stats_nucleotide_entity_seq())
            nec_df = c_df.merge(ne_df) if len(ne_df) > 0 else None
        else:
            nec_df = None
        return pec_df, nec_df


class SIFTSs(PDBs):
    def __new__(cls, iterable: Iterable=tuple()):
        return super(SIFTSs, cls).__new__(cls, (SIFTS(i) if isinstance(i, str) else i for i in iterable))

'''
TODO: Deal with carbohydrate polymer in PISA

    * possible change in struct_asym_id
'''

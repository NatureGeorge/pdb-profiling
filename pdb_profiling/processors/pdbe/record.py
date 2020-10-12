# @Created Date: 2020-08-11 10:48:08 pm
# @Filename: record.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-08-11 10:48:11 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Iterable, Union, Callable, Optional, Hashable, Dict, Coroutine, List, Tuple
from numpy import array, where as np_where, count_nonzero, nan, dot
from pathlib import Path
from pandas import isna, concat, DataFrame, Series, merge
from unsync import unsync, Unfuture
from asyncio import as_completed
from aiofiles import open as aiofiles_open
from smart_open import open as smart_open
from re import compile as re_compile
import orjson as json
from collections import defaultdict, namedtuple
from itertools import product
from pdb_profiling import default_id_tag
from pdb_profiling.utils import (init_semaphore, init_folder_from_suffix, 
                                 a_read_csv, split_df_by_chain, 
                                 related_dataframe, slice_series, 
                                 to_interval, MMCIF2DictPlus, 
                                 a_load_json, SeqRangeReader,
                                 sort_2_range, range_len,
                                 overlap_range, flat_dict_in_df,
                                 a_seq_parser, get_seq_from_parser,
                                 get_diff_index, get_seq_seg,
                                 get_gap_list,get_range_diff,
                                 outside_range_len,add_range,
                                 subtract_range, select_range,
                                 interval2set)
from pdb_profiling.processors.pdbe.api import ProcessPDBe, PDBeModelServer, PDBArchive, FUNCS as API_SET
from pdb_profiling.processors.uniprot.api import UniProtFASTA
from pdb_profiling.processors.pdbe import PDBeDB
from pdb_profiling.data import blosum62


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


class Base(object):
    
    folder = None

    def __init__(self):
        self.tasks = dict()

    def get_id(self):
        pass

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.get_id()}>"

    def register_task(self, key: Hashable, task: Unfuture):
        self.tasks[key] = task

    def set_neo4j_connection(self, api):
        pass

    def set_sqlite_connection(self, api):
        pass

    @classmethod
    @unsync
    async def set_web_semaphore(cls, web_semaphore_value):
        cls.web_semaphore = await init_semaphore(web_semaphore_value)

    @classmethod
    def get_web_semaphore(cls):
        return cls.web_semaphore

    @classmethod
    @unsync
    async def set_db_semaphore(cls, db_semaphore_value):
        cls.db_semaphore = await init_semaphore(db_semaphore_value)

    @classmethod
    def get_db_semaphore(cls):
        return cls.db_semaphore

    @classmethod
    def set_folder(cls, folder: Union[Path, str]):
        folder = Path(folder)
        assert folder.exists(), "Folder not exist! Please create it or input a valid folder!"
        cls.folder = folder

    @classmethod
    def get_folder(cls) -> Path:
        cls.check_folder()
        return cls.folder

    @classmethod
    def check_folder(cls):
        if cls.folder is None:
            raise ValueError("Please set folder via PDB.set_folder(folder: Union[Path, str])")
    
    def fetch_from_web_api(self, api_suffix: str, then_func: Optional[Callable[[Unfuture], Unfuture]] = None, json: bool = False, mask_id: str = None) -> Unfuture:
        assert api_suffix in API_SET, f"Invlaid API SUFFIX! Valid set:\n{API_SET}"
        task = self.tasks.get((api_suffix, then_func, json, mask_id), None)
        if task is not None:
            return task

        args = dict(pdb=self.get_id() if mask_id is None else mask_id,
                    suffix=api_suffix,
                    method='get',
                    folder=init_folder_from_suffix(
                        self.get_folder(), api_suffix),
                    semaphore=self.get_web_semaphore())
        if json:
            args['to_do_func'] = None
        task = ProcessPDBe.single_retrieve(**args)
        if then_func is not None:
            task = task.then(then_func)
        self.register_task((api_suffix, then_func, json, mask_id), task)
        return task

    @classmethod
    @unsync
    async def to_dataframe(cls, path):
        path = await path
        if path is None:
            return None
        df = await a_read_csv(path, sep="\t", converters=ProcessPDBe.converters)
        return df

    @classmethod
    @unsync
    async def to_dataframe_with_kwargs(cls, path):
        path = await path
        if path is None:
            return None
        df = await a_read_csv(path, sep="\t", converters=ProcessPDBe.converters, **cls.get_to_df_kwargs())
        return df

    @classmethod
    def set_to_df_kwargs(cls, **kwargs):
        cls.to_df_kwargs = kwargs

    @classmethod
    def get_to_df_kwargs(cls):
        return cls.to_df_kwargs


class PDB(Base):

    protein_sequence_pat = re_compile(r'([A-Z]{1}|\([A-Z0-9]+\))')
    nucleotide_sequence_pat = re_compile(r'([AUCG]{1}|\(DA\)|\(DT\)|\(DC\)|\(DG\)|\(UNK\))')
    StatsProteinEntitySeq = namedtuple(
        'StatsProteinEntitySeq', 'pdb_id molecule_type entity_id ca_p_only SEQRES_COUNT STD_INDEX STD_COUNT NON_INDEX NON_COUNT UNK_INDEX UNK_COUNT ARTIFACT_INDEX')
    StatsNucleotideEntitySeq = namedtuple(
        'StatsNucleotideEntitySeq', 'pdb_id molecule_type entity_id ca_p_only SEQRES_COUNT dNTP_INDEX dNTP_COUNT NTP_INDEX NTP_COUNT NON_INDEX NON_COUNT UNK_INDEX UNK_COUNT')
    StatsChainSeq = namedtuple(
        'StatsChainSeq', 'pdb_id entity_id chain_id struct_asym_id OBS_INDEX OBS_COUNT OBS_RATIO_SUM BINDING_LIGAND_INDEX BINDING_LIGAND_COUNT')

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
            var_property: self.pdb_ob.fetch_from_web_api(
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

    def __init__(self, pdb_id: str):
        self.check_folder()
        self.set_id(pdb_id)
        self.tasks = dict()
        self.pdb_ob = self
        self.properties_inited = False
    
    def set_id(self, pdb_id: str):
        assert default_id_tag(pdb_id) == 'pdb_id', "Invalid PDB ID!"
        self.pdb_id = pdb_id.lower()

    def get_id(self):
        return self.pdb_id

    def fetch_from_modelServer_api(self, api_suffix: str, method: str = 'post', data_collection=None, params=None, then_func: Optional[Callable[[Unfuture], Unfuture]] = None) -> Unfuture:
        assert api_suffix in PDBeModelServer.api_set, f"Invlaid API SUFFIX! Valid set:\n{PDBeModelServer.api_set}"
        task = self.tasks.get((PDBeModelServer.root, api_suffix, method, data_collection, params, then_func), None)
        if task is not None:
            return task
        task = PDBeModelServer.single_retrieve(
            pdb=self.pdb_id,
            suffix=api_suffix,
            method=method,
            folder=init_folder_from_suffix(self.get_folder()/'model-server', api_suffix),
            semaphore=self.get_web_semaphore(),
            data_collection=data_collection,
            params=params)
        if then_func is not None:
            task = task.then(then_func)
        self.register_task((PDBeModelServer.root, api_suffix, method, data_collection, params, then_func), task)
        return task

    def fetch_from_PDBArchive(self, api_suffix: str, then_func: Optional[Callable[[Unfuture], Unfuture]] = None, **kwargs) -> Unfuture:
        assert api_suffix in PDBArchive.api_set, f"Invlaid API SUFFIX! Valid set:\n{PDBArchive.api_set}"
        task = self.tasks.get((PDBArchive.root, api_suffix, then_func), None)
        if task is not None:
            return task
        task = PDBArchive.single_retrieve(
            pdb=self.pdb_id,
            suffix=api_suffix,
            folder=init_folder_from_suffix(self.get_folder()/'pdb/data/structures', api_suffix),
            semaphore=self.get_web_semaphore(),
            **kwargs)
        if then_func is not None:
            task = task.then(then_func)
        self.register_task((PDBArchive.root, api_suffix, then_func), task)
        return task

    @classmethod
    @unsync
    async def cif2residue_listing(cls, path: Unfuture):
        cols = ('_pdbx_poly_seq_scheme.asym_id',
                '_pdbx_poly_seq_scheme.entity_id',
                '_pdbx_poly_seq_scheme.seq_id',
                '_pdbx_poly_seq_scheme.mon_id',
                '_pdbx_poly_seq_scheme.ndb_seq_num',
                '_pdbx_poly_seq_scheme.pdb_seq_num',
                '_pdbx_poly_seq_scheme.pdb_strand_id',
                '_pdbx_poly_seq_scheme.pdb_ins_code')
        new_cols = ('struct_asym_id',
                    'entity_id',
                    'residue_number',
                    'residue_name',
                    'residue_number?',
                    'authore_residue_number',
                    'chain_id',
                    'author_insertion_code')
        with smart_open(await path) as handle:
            mmcif_dict = MMCIF2DictPlus(handle, cols)
            mmcif_dict['data_'] = mmcif_dict['data_'].lower()
            dfrm = DataFrame(mmcif_dict)
        col_dict = dict(zip(cols, new_cols))
        col_dict['data_'] = 'pdb_id'
        dfrm.rename(columns=col_dict, inplace=True)
        assert all(dfrm['residue_number'] == dfrm['residue_number?']
                   ), f"Unexpectd Cases: _pdbx_poly_seq_scheme.seq_id != _pdbx_poly_seq_scheme.ndb_seq_num\n{dfrm[dfrm['residue_number'] != dfrm['residue_number?']]}"
        dfrm.drop(columns=['residue_number?'], inplace=True)
        for col in ('residue_number', 'authore_residue_number'):
            dfrm[col] = dfrm[col].astype(int)
        dfrm.author_insertion_code = dfrm.author_insertion_code.apply(lambda x: '' if x == '.' else x)
        return dfrm

    @classmethod
    @unsync
    async def to_assg_oper_df(cls, path: Unfuture):
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
        
        assg_cols = ('_pdbx_struct_assembly_gen.asym_id_list',
                     '_pdbx_struct_assembly_gen.oper_expression',
                     '_pdbx_struct_assembly_gen.assembly_id')
        oper_cols = ('_pdbx_struct_oper_list.id', 
                     '_pdbx_struct_oper_list.symmetry_operation')
        async with aiofiles_open(await path, 'rt') as file_io:
            handle = await file_io.read()
            handle = (i+'\n' for i in handle.split('\n'))
            mmcif_dict = MMCIF2DictPlus(handle, assg_cols+oper_cols)
        
        if len(mmcif_dict) < 2:
            return None
        assg_df = DataFrame(list(zip(*[mmcif_dict[col] for col in assg_cols])), columns=assg_cols)
        assg_df = split_df_by_chain(assg_df, assg_cols, assg_cols[0:1]).rename(columns={col: col.split(
            '.')[1] for col in assg_cols}).rename(columns={'asym_id_list': 'struct_asym_id'})
        
        assg_df.oper_expression = assg_df.oper_expression.apply(
            lambda x: cls.expandOperators(cls.parseOperatorList(x)))
        assg_df = split_df_by_chain(
            assg_df, assg_df.columns, ('oper_expression', ), mode='list')
        assg_df.oper_expression = assg_df.oper_expression.apply(
            lambda x: json.dumps(x).decode('utf-8'))
        
        assg_dict = assg_df.groupby('assembly_id').oper_expression.unique().to_dict()
        rank_dict = assg_df.groupby(['assembly_id', 'struct_asym_id']).struct_asym_id.count().to_dict()
        rank_dict = {key: [value, 0] for key, value in rank_dict.items()}
        assg_df['model_id'] = assg_df.apply(lambda x: np_where(assg_dict[x['assembly_id']] == x['oper_expression'])[0][0]+1, axis=1)
        # assg_df.apply(lambda x: assg_dict[x['assembly_id']].index(x['oper_expression'])+1, axis=1)
        
        assg_df['asym_id_rank'] = assg_df.apply(lambda x: to_rank(
            rank_dict, x['assembly_id'], x['struct_asym_id']), axis=1)
        oper_dict = dict(zip(*[mmcif_dict[col] for col in oper_cols]))
        assg_df['symmetry_operation'] = assg_df.oper_expression.apply(
            lambda x: [oper_dict[i] for i in json.loads(x)])
        assg_df.symmetry_operation = assg_df.symmetry_operation.apply(
            lambda x: json.dumps(x).decode('utf-8'))
        assg_df.assembly_id = assg_df.assembly_id.astype(int)

        return assg_df

    @unsync
    async def pipe_assg_data_collection(self) -> str:
        '''demo_dict = {"atom_site": [{"label_asym_id": "A", "label_seq_id": 23}]}'''
        res_df = await self.pdb_ob.fetch_from_web_api('api/pdb/entry/residue_listing/', PDB.to_dataframe)
        # res_dict = iter_first(res_df, lambda row: row.observed_ratio > 0)
        # assert res_dict is not None, f"{self.pdb_ob}: Unexpected Cases, without observed residues?!"
        res_dict = res_df.loc[res_df.observed_ratio.gt(0).idxmax()].to_dict()
        demo_dict = dict(atom_site=[dict(
            label_asym_id=res_dict['struct_asym_id'],
            label_seq_id=int(res_dict['residue_number']))])
        return json.dumps(demo_dict).decode('utf-8')

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
        res_df = await self.fetch_from_web_api('api/pdb/entry/residue_listing/', self.to_dataframe)
        eec_df = res_df[['pdb_id', 'entity_id', 'chain_id',
                         'struct_asym_id']].drop_duplicates().reset_index(drop=True)
        eec_df['struct_asym_id_in_assembly'] = eec_df.struct_asym_id
        if merge_with_molecules_info:
            mol_df = await self.fetch_from_web_api('api/pdb/entry/molecules/', self.to_dataframe)
            self.res2eec_df = eec_df.merge(mol_df, how='left')
        else:
            self.res2eec_df = eec_df

    @unsync
    async def get_res2eec_df(self):
        try:
            return self.res2eec_df
        except AttributeError:
            await self.set_res2eec_df()
            return self.res2eec_df

    @classmethod
    @unsync
    async def assembly2eec(cls, path: Unfuture):
        '''
        NOTE except waters
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

    @unsync
    async def set_assembly(self, focus_assembly_ids:Optional[Iterable[int]]=None):
        '''
        NOTE even for NMR structures (e.g 1oo9), there exists assembly 1 for that entry
        NOTE Discard `details` is NaN -> not author_defined_assembly OR software_defined_assembly
        '''
        
        def to_assembly_id(pdb_id, assemblys):
            for assembly_id in assemblys:
                yield f"{pdb_id}/{assembly_id}"
        
        ass_eec_df = await self.fetch_from_web_api('api/pdb/entry/assembly/', PDB.to_dataframe)
        ass_eec_df = ass_eec_df[ass_eec_df.details.notnull()]
        assemblys = set(ass_eec_df.assembly_id) | {0}
        if focus_assembly_ids is not None:
            assemblys = sorted(assemblys & set(int(i) for i in focus_assembly_ids))
        else:
            assemblys = sorted(assemblys)
        self.assembly: Dict[int, PDBAssemble] = dict(zip(
            assemblys, 
            (PDBAssemble(ass_id, self) for ass_id in to_assembly_id(self.pdb_id, assemblys))))

    def get_assembly(self, assembly_id):
        return self.assembly[assembly_id]

    @unsync
    async def set_eec_as_df(self):
        def convert(struct_asym_id, asym_id_in_assembly):
            if not isna(struct_asym_id):
                return struct_asym_id, 1
            else:
                assert len(asym_id_in_assembly) > 1, f"id: {asym_id_in_assembly}"
                asym_id_rank = ord(asym_id_in_assembly[1])-63
                assert asym_id_rank > 1
                return asym_id_in_assembly[0], asym_id_rank

        ed_eec_df = await self.get_res2eec_df()
        ed_as_df = await self.fetch_from_web_api('api/pdb/entry/assembly/', PDB.assembly2eec) # NOTE may exists null pdb assembly
        if (ed_eec_df is None) or (ed_as_df is None):
            return None
        
        ed_eec_df = ed_eec_df[['pdb_id', 'entity_id', 'chain_id',
                               'struct_asym_id', 'struct_asym_id_in_assembly',
                               'molecule_type', 'molecule_name']]
        eec_as_df = ed_eec_df.merge(ed_as_df, how='outer')

        eec_as_df['asym_id_rank'] = 1
        eec_as_df[['struct_asym_id', 'asym_id_rank']] = array([convert(*i)
                                                                for i in zip(eec_as_df.struct_asym_id, eec_as_df.struct_asym_id_in_assembly)])
        chain_info = {(pdb_id, struct_asym_id): chain_id for pdb_id, chain_id, struct_asym_id in zip(
            eec_as_df.pdb_id, eec_as_df.chain_id, eec_as_df.struct_asym_id) if not isna(chain_id)}
        eec_as_df['chain_id'] = eec_as_df.apply(
            lambda x: chain_info[(x['pdb_id'], x['struct_asym_id'])], axis=1)
        eec_as_df.asym_id_rank = eec_as_df.asym_id_rank.astype(int)

        ed_eec_df_0 = ed_eec_df.copy()
        ed_eec_df_0['assembly_id'] = 0
        ed_eec_df_0['asym_id_rank'] = 1
        ed_eec_df_0['details'] = 'asymmetric_unit'
        eec_as_df = concat([eec_as_df, ed_eec_df_0]).sort_values(
            ['pdb_id', 'assembly_id', 'entity_id', 'struct_asym_id', 'asym_id_rank'])

        self.eec_as_df = eec_as_df

    @unsync
    async def get_eec_as_df(self):
        try:
            return self.eec_as_df
        except AttributeError:
            await self.set_eec_as_df()
            return self.eec_as_df
    
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
        store = []
        self.set_to_df_kwargs(usecols=lambda x: not x.startswith('author') and not x in ('struct_asym_id', 'chain_id'))
        muta_df = await self.fetch_from_web_api('api/pdb/entry/mutated_AA_or_NA/', self.to_dataframe_with_kwargs)
        if muta_df is not None:
            muta_df = muta_df.drop_duplicates().reset_index(drop=True)
            muta_df = flat_dict_in_df(muta_df, muta_df.mutation_details.apply(json.loads), ['type', 'from','to'])
            muta_dict = muta_df[muta_df['mutation_details.type'].isin(
                ('Cloning artifact', 'Expression tag', 'Initiating methionine', 'Linker'))].groupby(
                'entity_id').residue_number.apply(to_interval).to_dict()
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
        
        mol_df = await self.fetch_from_web_api('api/pdb/entry/molecules/', self.to_dataframe)
        mol_df.apply(regist_info, axis=1)
        return store

    @unsync
    async def stats_nucleotide_entity_seq(self):
        store = []

        def regist_info(record):
            if record['molecule_type'] not in ('polydeoxyribonucleotide', 'polyribonucleotide', 'polydeoxyribonucleotide/polyribonucleotide'):
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

        mol_df = await self.fetch_from_web_api('api/pdb/entry/molecules/', self.to_dataframe)
        mol_df.apply(regist_info, axis=1)
        return store

    @unsync
    async def stats_poly_chain_obs_seq(self):
        store = defaultdict(list)
        def regist_info(record):
            store[(record['pdb_id'], record['entity_id'], record['chain_id'], record['struct_asym_id'])
            ].append(
                (json.loads(record['start'])['residue_number'], json.loads(record['end'])['residue_number']))

        c_df = await self.fetch_from_web_api('api/pdb/entry/polymer_coverage/', self.to_dataframe)
        c_df.apply(regist_info, axis=1)
        return [self.StatsChainSeq(*key, value, range_len(value)) for key, value in store.items()]
    
    @unsync
    async def stats_chain_seq(self):
        store = defaultdict(list)
        
        def regist_info(record):
            store[(record['pdb_id'], record['entity_id'], record['chain_id'], record['struct_asym_id'])
                  ].append((record['residue_number'], record['observed_ratio']))
        
        res_df = await self.fetch_from_web_api('api/pdb/entry/residue_listing/', self.to_dataframe)
        # mol_df = await self.fetch_from_web_api('api/pdb/entry/molecules/', self.to_dataframe)
        # res_df = res_df.merge(mol_df[mol_df.molecule_type.isin(('polypeptide(L)', 'polypeptide(D)'))][['entity_id']])
        res_df[res_df.observed_ratio.gt(0)].apply(lambda x: regist_info(x), axis=1)

        self.set_to_df_kwargs(usecols=lambda x: not x.startswith('author'))
        br_df = await self.fetch_from_web_api('api/pdb/entry/binding_sites/', self.to_dataframe_with_kwargs)
        if br_df is None:
            return [self.StatsChainSeq(*key,
                                       to_interval(i[0] for i in value),
                                       len(value),
                                       sum(i[1] for i in value),
                                       tuple(),
                                       0
                                       ) for key, value in store.items()]
        else:
            br_dict = br_df[br_df.entity_id.ne(-1)].groupby('struct_asym_id').residue_number.apply(lambda x: (to_interval(x), len(frozenset(x)))).to_dict()
            return [self.StatsChainSeq(*key, 
                                    to_interval(i[0] for i in value), 
                                    len(value), 
                                    sum(i[1] for i in value),
                                    *br_dict.get(key[-1], (tuple(), 0))
                                    ) for key, value in store.items()]

    @unsync
    async def stats_chain(self):
        pe_df = DataFrame(await self.stats_protein_entity_seq())
        ne_df = DataFrame(await self.stats_nucleotide_entity_seq())
        c_df = DataFrame(await self.stats_chain_seq())
        try:
            pec_df = c_df.merge(pe_df)
            pec_df['OBS_STD_INDEX'] = pec_df.apply(lambda x: overlap_range(x['OBS_INDEX'], x['STD_INDEX']), axis=1)
            pec_df['OBS_STD_COUNT'] = pec_df.OBS_STD_INDEX.apply(range_len)
        except Exception:
            pec_df = None
        nec_df = c_df.merge(ne_df) if len(ne_df) > 0 else None
        return pec_df, nec_df

    @unsync
    async def get_sequence(self, **kwargs):
        '''
        Get SEQRES Sequence via entity_id | chain_id (default protein) | struct_asym_id
        '''
        sequence_col = 'sequence' if kwargs.get('one_letter_code', True) else 'pdb_sequence'
        mol_df = await self.fetch_from_web_api('api/pdb/entry/molecules/', self.to_dataframe)
        if 'entity_id' in kwargs:
            return mol_df.loc[mol_df.entity_id.eq(kwargs['entity_id']).idxmax(), sequence_col]
        elif 'struct_asym_id' in kwargs:
            struct_asym_id = kwargs['struct_asym_id']
            return mol_df.loc[mol_df.in_struct_asyms.apply(lambda x: f'"{struct_asym_id}"' in x).idxmax(), sequence_col]
        elif 'chain_id' in kwargs:
            chain_id = kwargs['chain_id']
            return mol_df.loc[(mol_df.molecule_type.eq('polypeptide(L)') & mol_df.in_chains.apply(lambda x: f'"{chain_id}"' in x)).idxmax(), sequence_col]
        else:
            raise ValueError(f"Cannot get sequence with specified information: {kwargs}")

    @unsync
    async def profile_id(self):
        '''except water'''
        assg_oper_df = await self.fetch_from_modelServer_api(
                'atoms',
                data_collection=await self.pipe_assg_data_collection(), 
                then_func=PDB.to_assg_oper_df)
        res2eec_df = await self.get_res2eec_df()
        focus_res2eec_df = res2eec_df[['pdb_id', 'entity_id', 'molecule_type', 'chain_id', 'struct_asym_id']]
        add_0_assg_oper_df = focus_res2eec_df.copy()
        add_0_assg_oper_df['assembly_id'] = 0
        add_0_assg_oper_df['model_id'] = 1
        add_0_assg_oper_df['asym_id_rank'] = 1
        add_0_assg_oper_df['oper_expression'] = ''
        add_0_assg_oper_df['symmetry_operation'] = ''

        if assg_oper_df is not None:
            focus_assg_oper_df = assg_oper_df[assg_oper_df.struct_asym_id.isin(focus_res2eec_df.struct_asym_id)]
            new_focus_assg_oper_df = concat([add_0_assg_oper_df, focus_assg_oper_df.merge(focus_res2eec_df, how='outer')])
            assert any(new_focus_assg_oper_df.isnull().sum()) is False, f"Unexpected Cases {new_focus_assg_oper_df}"
        else:
            new_focus_assg_oper_df = add_0_assg_oper_df
        return new_focus_assg_oper_df


class PDBAssemble(PDB):

    id_pattern = re_compile(r"([a-z0-9]{4})/([0-9]+)")
    struct_range_pattern = re_compile(r"\[.+\]([A-Z]+):[0-9]+")  # e.g. [FMN]B:149 [C2E]A:301
    rare_pat = re_compile(r"([A-Z]+)_([0-9]+)")  # e.g. 2rde assembly 1 A_1, B_1...

    @property
    def assemble_summary(self) -> Dict:
        for ass in self.summary['assemblies']:
            if int(ass['assembly_id']) == self.assembly_id:
                return ass
        raise ValueError(f"{repr(self)}: Without expected assemble info\n{self.summary['assemblies']}")

    def __init__(self, pdb_ass_id, pdb_ob: Optional[PDB]=None):
        super().__init__(pdb_ass_id)
        if pdb_ob is None:
            self.pdb_ob = PDB(self.pdb_id)
        else:
            self.pdb_ob = pdb_ob
        self.interface_filters = {
            'structure_2.symmetry_operator': ('eq', 'x,y,z'), 
            'css': ('ge', 0)}

    def set_id(self, pdb_ass_id: str):
        self.pdb_ass_id = pdb_ass_id.lower()
        try:
            self.pdb_id, self.assembly_id = self.id_pattern.fullmatch(self.pdb_ass_id).groups()
        except AttributeError:
            raise ValueError(f"Invalid ID: {self.pdb_ass_id}")
        self.assembly_id = int(self.assembly_id)
    
    def get_id(self):
        return self.pdb_ass_id

    @classmethod
    async def to_interfacelist_df(cls, path: Unfuture):
        def transform(x):
            res = cls.rare_pat.search(x)
            assert bool(res), "Unexpected Case"
            chain, num = res.groups()
            num = int(num)
            if num == 1:
                return chain
            else:
                return chain+chr(63+num)

        interfacelist_df = await path.then(cls.to_dataframe)
        if interfacelist_df is None:
            return None
        interfacelist_df.rename(columns={
                                "id": "interface_id", 
                                "pdb_code": "pdb_id", 
                                "assemble_code": "assembly_id"
                                }, inplace=True)
        if any('_' in i for i in interfacelist_df['structure_1.range']):
            interfacelist_df['structure_1.range'] = interfacelist_df['structure_1.range'].apply(lambda x: transform(x) if '_' in x else x)
            interfacelist_df['struct_asym_id_in_assembly_1'] = interfacelist_df['structure_1.range']
        else:
            interfacelist_df['struct_asym_id_in_assembly_1'] = interfacelist_df.apply(
                lambda x: cls.struct_range_pattern.match(x['structure_1.range']).group(1) if x['structure_1.original_range'] != '{-}' else x['structure_1.range'], axis=1)
        if any('_' in i for i in interfacelist_df['structure_2.range']):
            interfacelist_df['structure_2.range'] = interfacelist_df['structure_2.range'].apply(lambda x: transform(x) if '_' in x else x)
            interfacelist_df['struct_asym_id_in_assembly_2'] = interfacelist_df['structure_2.range']
        else:
            interfacelist_df['struct_asym_id_in_assembly_2'] = interfacelist_df.apply(
                lambda x: cls.struct_range_pattern.match(x['structure_2.range']).group(1) if x['structure_2.original_range'] != '{-}' else x['structure_2.range'], axis=1)
        return interfacelist_df

    @unsync
    async def set_interface(self, obligated_class_chains: Optional[Iterable[str]] = None, allow_same_class_interaction:bool=True):

        def to_interface_id(pdb_assembly_id, focus_interface_ids):
            for interface_id in focus_interface_ids:
                yield f"{pdb_assembly_id}/{interface_id}"

        interfacelist_df = await self.fetch_from_web_api('api/pisa/interfacelist/', self.to_interfacelist_df)
        
        if interfacelist_df is None:
            self.interface = dict()
            return
        
        if obligated_class_chains is not None:
            bool_1 = interfacelist_df.struct_asym_id_in_assembly_1.isin(obligated_class_chains)
            bool_2 = interfacelist_df.struct_asym_id_in_assembly_2.isin(obligated_class_chains)
            interfacelist_df = interfacelist_df[bool_1 | bool_2]
            if not allow_same_class_interaction:
                interfacelist_df = interfacelist_df[~(bool_1 & bool_2)]
        
        focus_interface_ids = related_dataframe(
            self.interface_filters, interfacelist_df).interface_id.unique()

        self.interface: Dict[int, PDBInterface] = dict(zip(
            focus_interface_ids, (PDBInterface(if_id, self) for if_id in to_interface_id(self.get_id(), focus_interface_ids))))

    def get_interface(self, interface_id):
        return self.interface[interface_id]

    @unsync
    async def set_assemble_eec_as_df(self):
        eec_as_df = await self.pdb_ob.get_eec_as_df()
        assert self.assembly_id in eec_as_df.assembly_id, f"{repr(self)}: Invalid assembly_id!"
        self.assemble_eec_as_df = eec_as_df[eec_as_df.assembly_id.eq(self.assembly_id)]
    
    @unsync
    async def get_assemble_eec_as_df(self):
        try:
            return self.assemble_eec_as_df
        except AttributeError:
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
        if len(eec_as_df) == 1:
            self.interface = dict()
            return
        protein_type_asym = eec_as_df[
            eec_as_df.molecule_type.eq('polypeptide(L)')].struct_asym_id_in_assembly
        self.interface_filters['struct_asym_id_in_assembly_1'] = ('isin', protein_type_asym)
        self.interface_filters['struct_asym_id_in_assembly_2'] = ('isin', protein_type_asym)
        await self.set_interface()

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
        eec_as_df = await self.get_assemble_eec_as_df()
        if len(eec_as_df) == 1:
            self.interface = dict()
            return
        molType_dict = eec_as_df.groupby('molecule_type').struct_asym_id_in_assembly.apply(set).to_dict()

        protein_type_asym = molType_dict.get('polypeptide(L)', set())
        ligand_type_asym = molType_dict.get('bound', set())
        target_type_asym = protein_type_asym | ligand_type_asym

        self.interface_filters['struct_asym_id_in_assembly_1'] = ('isin', target_type_asym)
        self.interface_filters['struct_asym_id_in_assembly_2'] = ('isin', target_type_asym)
        await self.set_interface(obligated_class_chains=protein_type_asym, allow_same_class_interaction=False)


class PDBInterface(PDBAssemble):
 
    id_pattern = re_compile(r"([a-z0-9]{4})/([0-9]+)/([0-9]+)")

    def __init__(self, pdb_ass_int_id, pdbAssemble_ob: Optional[PDBAssemble]=None):
        super().__init__(pdb_ass_int_id)
        if pdbAssemble_ob is None:
            self.pdbAssemble_ob = PDBAssemble(f"{self.pdb_id}/{self.assembly_id}")
        else:
            self.pdbAssemble_ob = pdbAssemble_ob

    def set_id(self, pdb_ass_int_id: str):
        self.pdb_ass_int_id = pdb_ass_int_id.lower()
        try:
            self.pdb_id, self.assembly_id, self.interface_id = self.id_pattern.fullmatch(
                self.pdb_ass_int_id).groups()
        except AttributeError:
            raise ValueError(f"Invalid ID: {self.pdb_ass_int_id}")
        self.assembly_id = int(self.assembly_id)
        self.interface_id = int(self.interface_id)
    
    def get_id(self):
        return self.pdb_ass_int_id

    @classmethod
    async def to_interfacedetail_df(cls, path: Unfuture):

        def check_struct_selection(interfacedetail_df, colName):
            sele = next(iter(interfacedetail_df[colName]))
            sele_m = cls.struct_range_pattern.fullmatch(sele)
            if bool(sele_m):
                interfacedetail_df[colName] = sele_m.group(1)

        cls.set_to_df_kwargs(usecols=['pdb_code', 'assemble_code', 'interface_number', 'chain_id',
                                      'residue', 'sequence', 'insertion_code',
                                      'buried_surface_area', 'solvent_accessible_area',
                                      'interface_detail.interface_structure_1.structure.selection',
                                      'interface_detail.interface_structure_2.structure.selection'
                                      ],
                             na_values=[' ', '?'])
        interfacedetail_df = await path.then(cls.to_dataframe_with_kwargs)
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
        check_struct_selection(interfacedetail_df, 's1_selection')
        check_struct_selection(interfacedetail_df, 's2_selection')
        return interfacedetail_df

    @unsync
    async def set_interface_res(self):
        interfacedetail_df = await self.fetch_from_web_api('api/pisa/interfacedetail/', self.to_interfacedetail_df)
        if interfacedetail_df is None:
            return None
        else:
            struct_sele_set = set(interfacedetail_df.head(1)[['s1_selection', 's2_selection']].to_records(index=False)[0])
        eec_as_df = await self.pdbAssemble_ob.get_assemble_eec_as_df()
        res_df = await self.pdbAssemble_ob.pdb_ob.fetch_from_web_api('api/pdb/entry/residue_listing/', PDB.to_dataframe)
        interfacedetail_df = interfacedetail_df.merge(eec_as_df, how="left")
        interfacedetail_df = interfacedetail_df.merge(res_df, how="left")
        check = interfacedetail_df[interfacedetail_df.residue_number.isnull()]
        self.interface_res_df = interfacedetail_df
        if len(check) != 0:
            raise ValueError(f"Unexcepted Data in Residue DataFrame: {check.head(1).to_dict('record')[0]}")
        
        focus_cols = ['pdb_id', 'entity_id', 'chain_id', 'struct_asym_id', 
                      'struct_asym_id_in_assembly', 'asym_id_rank',
                      'assembly_id', 'interface_id', 'molecule_type', 'molecule_name',
                      'residue_number', 'buried_surface_area', 'solvent_accessible_area']
        nda = interfacedetail_df[focus_cols].to_numpy()

        def yield_record():
            for _, (start, end) in slice_series(interfacedetail_df.struct_asym_id_in_assembly).items():
                asa_col = focus_cols.index('solvent_accessible_area')
                bsa_col = focus_cols.index('buried_surface_area')
                res_col = focus_cols.index('residue_number')
                asa_index = start+np_where(nda[start:end, asa_col] > 0)[0]
                asa_res = nda[asa_index, res_col]
                bsa_index = asa_index[np_where(nda[asa_index, bsa_col] > 0)]
                bsa_res = nda[bsa_index, res_col]
                info = nda[start, 0:len(focus_cols)-3].tolist() + \
                    [json.dumps(to_interval(asa_res)).decode('utf-8'),
                    json.dumps(to_interval(bsa_res)).decode('utf-8')]
                yield dict(zip(focus_cols[:-3]+['surface_range', 'interface_range'], info))
        
        records = sorted(yield_record(), key=lambda x: x['struct_asym_id_in_assembly'])

        common_keys = ('pdb_id', 'assembly_id', 'interface_id')
        record1 = {f"{key}_1": value for key,
                value in records[0].items() if key not in common_keys}
        struct_sele_set = struct_sele_set - {record1['struct_asym_id_in_assembly_1']}
        if len(records) == 2:
            record2 = {f"{key}_2": value for key,
                    value in records[1].items() if key not in common_keys}
        else:
            saiia2 = struct_sele_set.pop()
            record2 = {'struct_asym_id_in_assembly_2': saiia2}
            cur_keys = list(set(focus_cols[:-3])-set(common_keys)-{'struct_asym_id_in_assembly'})
            try:
                cur_record = eec_as_df[eec_as_df.struct_asym_id_in_assembly.eq(saiia2)][cur_keys].to_dict('record')[0]
            except Exception:
                raise ValueError(f"\n{self.get_id()},\n{saiia2},\n{eec_as_df}")
            for key, value in cur_record.items():
                record2[f"{key}_2"] = value

        record_dict = {**record1, **record2}
        for key in common_keys:
            record_dict[key] = records[0][key]
        self.interface_res_dict = record_dict

    @unsync
    async def get_interface_res_df(self):
        try:
            return self.interface_res_df
        except AttributeError:
            await self.set_interface_res()
            return self.interface_res_df
    
    @unsync
    async def get_interface_res_dict(self):
        try:
            return self.interface_res_dict
        except AttributeError:
            await self.set_interface_res()
            return self.interface_res_dict


class SIFTS(PDB):

    EntityChain = namedtuple('EntityChain', 'pdb_id entity_chain_info entity_count chain_count')
    UniProtEntity = namedtuple('UniProtEntity', 'pdb_id unp_entity_info entity_unp_info entity_with_unp_count min_unp_count')
    OligoState = namedtuple('OligoState', 'pdb_id oligo_state has_unmapped_protein')
    MappingMeta = namedtuple('MappingMeta', 'UniProt species identity')

    @classmethod
    def set_folder(cls, folder: Union[Path, str]):
        super().set_folder(folder)
        cls.sqlite_api = PDBeDB(
            "sqlite:///%s" % (init_folder_from_suffix(cls.folder, 'local_db')/"PDBeDB.db"))

    def set_id(self, identifier: str):
        tag = default_id_tag(identifier, None)
        if tag == 'pdb_id':
            self.level = 'PDB Entry'
            self.identifier = identifier.lower()
            self.pdb_id = self.identifier
        elif tag == 'UniProt':
            self.level = tag
            self.identifier = identifier.upper()
        else:
            raise AssertionError(f"Invalid identifier: <{identifier}, {tag}>")

    def get_id(self):
        return self.identifier

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.level} {self.get_id()}>"

    @staticmethod
    @unsync
    async def reformat(dfrm: Union[DataFrame, Unfuture, Coroutine]) -> DataFrame:
        if isinstance(dfrm, (Coroutine, Unfuture)):
            dfrm = await dfrm
        if 'pdb_start' not in dfrm.columns:
            flat_dict_in_df(dfrm, dfrm.start.apply(json.loads), ('residue_number',))
            flat_dict_in_df(dfrm, dfrm.end.apply(json.loads), ('residue_number',))
            dfrm.rename(columns={
                'start.residue_number': 'pdb_start', 
                'end.residue_number': 'pdb_end'}, inplace=True)
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
    @unsync
    async def dealWithInDel(dfrm: Union[DataFrame, Unfuture, Coroutine], sort_by_unp: bool = True) -> DataFrame:
        if isinstance(dfrm, (Coroutine, Unfuture)):
            dfrm = await dfrm

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

        dfrm.pdb_range = dfrm.pdb_range.apply(json.loads)
        dfrm.unp_range = dfrm.unp_range.apply(json.loads)
        dfrm['group_info'] = dfrm.apply(lambda x: len(
            x['pdb_range']), axis=1)

        focus_index = dfrm[dfrm.group_info.gt(1)].index
        if sort_by_unp and (len(focus_index) > 0):
            focus_df = dfrm.loc[focus_index].apply(lambda x: sort_2_range(
                x['unp_range'], x['pdb_range']), axis=1, result_type='expand')
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
        add_tage_to_range(dfrm, tage_name='sifts_range_tag')
        dfrm['repeated'] = dfrm.apply(
            lambda x: x['diff-'] > 0 and x['sifts_range_tag'] != 'Insertion (Specail Case)', axis=1)
        dfrm['repeated'] = dfrm.apply(
            lambda x: True if any(i < 0 for i in x['unp_gaps']) else x['repeated'], axis=1)
        dfrm['reversed'] = dfrm.pdb_gaps.apply(lambda x: any(i < 0 for i in x))
        dfrm.pdb_range = dfrm.pdb_range.apply(
            lambda x: json.dumps(x).decode('utf-8'))
        dfrm.unp_range = dfrm.unp_range.apply(
            lambda x: json.dumps(x).decode('utf-8'))
        dfrm['InDel_sum'] = dfrm.pdb_gaps.apply(sum) + dfrm.unp_gaps.apply(sum) + dfrm.range_diff.apply(sum)
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
            return min(len(set(res)) for res in product(*info_dict.values()))
        
        mol_df = await self.fetch_from_web_api('api/pdb/entry/molecules/', self.to_dataframe)
        mol_df = mol_df[mol_df.molecule_type.eq('polypeptide(L)')]
        info_dict = dict(zip(mol_df.entity_id,mol_df.in_chains))
        entity_count = len(info_dict)
        chain_count = sum(i.count(',')+1 for i in info_dict.values())

        try:
            best_iso = await self.fetch_from_web_api('api/mappings/isoforms/', self.to_dataframe).then(self.reformat)
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
            for (lstart, lend), (rstart, rend) in zip(lrange, rrange):
                assert (lend - lstart) == (rend-rstart), "convert_index(): Invalid range"
                if (site >= rstart) and (site <= rend):
                    return int(site + lstart - rstart)
                else:
                    continue
            return nan_value
        lrange = json.loads(lrange) if isinstance(lrange, str) else lrange
        rrange = json.loads(rrange) if isinstance(rrange, str) else rrange
        return tuple(unit(lrange, rrange, site) for site in sites)

    @staticmethod
    @unsync
    async def get_residue_conflict(pdb_id, entity_id, pdb_range, Entry, UniProt, unp_range, on_pdb: bool = True):
        pdb_seq = await PDB(pdb_id).get_sequence(entity_id=entity_id)
        unp_seq = await get_seq_from_parser(UniProtFASTA.single_retrieve(Entry).then(a_seq_parser), UniProt)
        return to_interval(get_diff_index(pdb_seq, pdb_range, unp_seq, unp_range, on_pdb))

    @classmethod
    @unsync
    async def add_residue_conflict(cls, dfrm: Union[DataFrame, Tuple, Unfuture, Coroutine], on_pdb: bool = True):
        if isinstance(dfrm, (Coroutine, Unfuture)):
            dfrm = await dfrm
        if isinstance(dfrm, Tuple):
            dfrm = dfrm[0]
        pdb_range_col = 'new_pdb_range' if 'new_pdb_range' in dfrm.columns else 'pdb_range'
        unp_range_col = 'new_unp_range' if 'new_unp_range' in dfrm.columns else 'unp_range'
        focus = ['UniProt', 'entity_id', 'pdb_id']
        f_dfrm = dfrm[focus+[pdb_range_col, 'Entry', unp_range_col]].drop_duplicates(subset=focus).reset_index(drop=True)
        tasks = f_dfrm.apply(lambda x: cls.get_residue_conflict(
            x['pdb_id'], x['entity_id'], x[pdb_range_col], x['Entry'], x['UniProt'], x[unp_range_col]), axis=1)
        f_dfrm['conflict_pdb_range'] = [await i for i in tasks]
        dfrm_ed = merge(dfrm, f_dfrm)
        assert dfrm_ed.shape[0] == dfrm.shape[0]
        return dfrm_ed

    @unsync
    async def fetch_new_residue_mapping(self, entity_id, start, end):
        df = await self.fetch_from_web_api(
            'graph-api/residue_mapping/', 
            self.to_dataframe, 
            mask_id=f'{self.pdb_id}/{entity_id}/{start}/{end}')
        if df is not None:
            await self.sqlite_api.async_insert(self.sqlite_api.ResidueMapping, df.to_dict('records'))
        return df

    @unsync
    async def fetch_residue_mapping(self, entity_id:int, start:int, end:int, columns:str='UniProt,residue_number,unp_residue_number', usecols:bool=True, UniProt:Optional[str]=None):
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
                dfs.append(new_df)
            return concat(dfs, ignore_index=True)
        else:
            df = await self.fetch_new_residue_mapping(entity_id, start, end)
            if df is not None:
                df = df[df.UniProt.eq(UniProt)] if UniProt is not None else df
                return df[columns] if usecols else df
            else:
                return
    
    @classmethod
    @unsync
    async def fix_range(cls, dfrm: Union[DataFrame, Tuple, Unfuture, Coroutine]):
        if isinstance(dfrm, (Coroutine, Unfuture)):
            dfrm = await dfrm
        if isinstance(dfrm, Tuple):
            dfrm = dfrm[0]
        
        focus = ['UniProt', 'entity_id', 'pdb_id']
        f_dfrm = dfrm[focus+['pdb_range', 'unp_range', 'Entry', 'range_diff', 'sifts_range_tag']].drop_duplicates(subset=focus)
        f_dfrm = f_dfrm[f_dfrm.sifts_range_tag.isin(('Deletion', 'Insertion_Undivided', 'InDel_2', 'InDel_3'))].reset_index(drop=True)

        if len(f_dfrm) > 0:
            tasks = f_dfrm.apply(lambda x: cls.a_sliding_alignment(
                x['range_diff'], x['pdb_range'], x['unp_range'], x['pdb_id'], x['entity_id'], x['UniProt'], x['Entry']), axis=1)
            res = [await i for i in tasks]
            f_dfrm[['new_pdb_range', 'new_unp_range']] = DataFrame(list(zip(*i)) for i in res)
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
    def sliding_alignment_score(range_diff, pdb_seq, pdb_range, unp_seq, unp_range, **kwargs):
        def get_optimal_range(abs_diff, seg_to_add, seg_to_ori, lstart, lend, rstart, rend, on_left):
            gap_seg = '-' * abs_diff
            res = tuple(sum(blosum62.get((l, r), blosum62.get((r, l), 0)
                                        ) for l, r in zip(seg_to_add[:i] + gap_seg + seg_to_add[i:], seg_to_ori)
                            ) for i in range(len(seg_to_add)+1))
            max_val = max(res)
            index = res.index(max_val)
            assert index > 0
            # assert max_val not in res[index+1:]
            '''
            if max_val in res[index+1:]:
                import matplotlib.pyplot as plt
                plt.figure(figsize=(10,8))
                plt.ylim(max_val-100, max_val+30)
                plt.xlim(index-50, index+50)
                plt.scatter(tuple(range(len(res))),res, label=f"{kwargs}")
                plt.plot(tuple(range(len(res))), res)
                plt.legend()
                plt.show()
            '''
            yield (lstart, lstart+index), (rstart, rstart+index)
            if on_left:
                yield (lstart+index+1, lend), (rstart+index+1+abs_diff, rend)
            else:
                yield (lstart+index+1+abs_diff, lend), (rstart+index+1, rend)

        for diff, (lstart, lseg), (rstart, rseg) in zip(range_diff, get_seq_seg(pdb_seq, pdb_range), get_seq_seg(unp_seq, unp_range)):

            lend = lstart + len(lseg) - 1
            rend = rstart + len(rseg) - 1

            if diff > 0:
                yield from get_optimal_range(diff, lseg, rseg, lstart, lend, rstart, rend, True)
            elif diff < 0:
                yield from get_optimal_range(-diff, rseg, lseg, lstart, lend, rstart, rend, False)
            else:
                yield (lstart, lend), (rstart, rend)
                continue

    @classmethod
    @unsync
    async def a_sliding_alignment(cls, range_diff, pdb_range, unp_range, pdb_id, entity_id, UniProt, Entry):
        pdb_seq = await PDB(pdb_id).get_sequence(entity_id=entity_id)
        unp_seq = await get_seq_from_parser(UniProtFASTA.single_retrieve(Entry).then(a_seq_parser), UniProt)
        pdb_range = json.loads(pdb_range) if isinstance(pdb_range, str) else pdb_range
        unp_range = json.loads(unp_range) if isinstance(unp_range, str) else unp_range
        return cls.sliding_alignment_score(range_diff, pdb_seq, pdb_range, unp_seq, unp_range, 
                pdb_id=pdb_id,entity_id=entity_id,UniProt=UniProt)

    @classmethod
    async def add_unp_len(cls, dfrm: Union[DataFrame, Tuple, Unfuture, Coroutine]):
        if isinstance(dfrm, (Coroutine, Unfuture)):
            dfrm = await dfrm
        if isinstance(dfrm, Tuple):
            dfrm = dfrm[0]
        f_dfrm = dfrm[['UniProt', 'Entry']].drop_duplicates().reset_index(drop=True)
        tasks = f_dfrm.apply(lambda x: get_seq_from_parser(UniProtFASTA.single_retrieve(x['Entry']).then(a_seq_parser), x['UniProt']), axis=1)
        f_dfrm['unp_len'] = [len(await i) for i in tasks]
        dfrm_ed = merge(dfrm, f_dfrm)
        assert dfrm_ed.shape[0] == dfrm.shape[0]
        return dfrm_ed

    @unsync
    async def pipe_score(self):
        weight = array([1, -1, -1, -1.79072623, -2.95685934, -4.6231746])
        def raw_score(vec): return dot(vec, weight)
        
        try:
            sifts_df = await self.fetch_from_web_api('api/mappings/all_isoforms/', self.to_dataframe
                ).then(self.reformat
                ).then(self.dealWithInDel
                ).then(self.fix_range
                ).then(self.add_residue_conflict
                ).then(self.add_unp_len)
        except Exception:
            return
        if self.level == 'PDB Entry':
            pec_df, _ = await self.stats_chain()
        else:
            tasks = await PDBs(sifts_df.pdb_id.unique()).fetch('stats_chain').run()
            pec_df, _ = zip(*tasks)
            pec_df = concat(pec_df, sort=False, ignore_index=True)
        
        full_df = sifts_df.merge(pec_df)
        assert sifts_df.shape[0] == full_df.shape[0] # f"{pec_df[pec_df.pdb_id.eq('6mzl')]}"
        
        s1 = full_df.apply(lambda x: range_len(overlap_range(x['OBS_STD_INDEX'], x['new_pdb_range'])), axis=1)
        s2 = full_df.apply(lambda x: outside_range_len(subtract_range(x['new_pdb_range'], x['ARTIFACT_INDEX']), x['SEQRES_COUNT']), axis=1)
        s3 = full_df.BINDING_LIGAND_COUNT
        s4 = full_df.new_pdb_range.apply(range_len) - full_df.apply(
            lambda x: range_len(overlap_range(x['OBS_INDEX'], x['new_pdb_range'])), axis=1)
        s5 = full_df.apply(lambda x: range_len(overlap_range(
            overlap_range(x['OBS_INDEX'], x['new_pdb_range']), 
            add_range(x['conflict_pdb_range'], x['NON_INDEX']))), axis=1)
        s6 = full_df.InDel_sum
        s_df = DataFrame(dict(s1=s1,s2=s2,s3=s3,s4=s4,s5=s5,s6=s6))
        s_df['RAW_BS'] = s_df.apply(raw_score, axis=1) / full_df.unp_len

        exp_cols = ['pdb_id', 'resolution', 'experimental_method_class', 'experimental_method', 'multi_method']
        if self.level == 'PDB Entry':
            exp_df = DataFrame(await PDBs.fetch_exp_pipe(self), columns=exp_cols)
            summary_dict = await self.fetch_from_web_api('api/pdb/entry/summary/', a_load_json, json=True)
            d = summary_dict[self.pdb_id][0]
            exp_df['revision_date'] = d['revision_date']
            exp_df['deposition_date'] = d['deposition_date']
        else:
            pdbs = PDBs(sifts_df.pdb_id.unique())
            tasks = await pdbs.fetch(PDBs.fetch_exp_pipe).run()
            exp_df = DataFrame((j for i in tasks for j in i), columns=exp_cols)
            r_len = exp_df.shape[0]
            tasks = await pdbs.fetch(PDBs.fetch_date).run()
            exp_df = exp_df.merge(DataFrame(tasks, columns=['pdb_id', 'revision_date','deposition_date']))
            assert r_len == exp_df.shape[0]

        return full_df, s_df, exp_df

    @unsync
    async def pipe_select_mo(self):
        assert self.level == 'UniProt'
        full_df, s_df, exp_df = await self.pipe_score()
        sele_df = merge(
            concat([full_df, s_df[['RAW_BS']]], axis=1).query(
                'UNK_COUNT == 0 and ca_p_only == False and identity >=0.9 and repeated == False and reversed == False and OBS_COUNT > 20'),
            exp_df.query(
                '(experimental_method == "X-ray diffraction" and resolution <= 3) or experimental_method == "Solution NMR"'))
        sele_df['1/resolution'] = 1 / sele_df.resolution
        sele_df['id_score'] = sele_df.chain_id.apply(lambda x: -sum(ord(i) for i in x))
        sele_df['COV_SCORE'] = sele_df.RAW_BS * sele_df.unp_len / len(frozenset.union(*sele_df.new_unp_range.apply(interval2set).to_list()))
        sele_df = sele_df[sele_df.RAW_BS>0]
        sele_index = select_range(sele_df.new_unp_range, sele_df.sort_values(by=['RAW_BS', '1/resolution', 'revision_date', 'id_score'], ascending=False).index)
        sele_df['select_tag'] = False
        sele_df.loc[sele_index, 'select_tag'] = True
        return sele_df


class Compounds(Base):

    def set_id(self, identifier: str):
        assert len(identifier) > 0, "Empty string is not a valid identifier!"
        self.identifier = identifier.upper()

    def get_id(self):
        return self.identifier

    def __init__(self, identifier:str):
        self.check_folder()
        self.set_id(identifier)
        self.tasks = dict()


class PDBs(list):
    
    def __init__(self, iterable:Iterable):
        super().__init__((PDB(i) if isinstance(i, str) else i for i in iterable))

    @staticmethod
    @unsync
    async def fetch_exp_pipe(pdb_ob:PDB):
        exp_data = await pdb_ob.fetch_from_web_api('api/pdb/entry/experiment/', a_load_json, json=True)
        multi_method = len(exp_data[pdb_ob.pdb_id]) > 1
        return ((pdb_ob.pdb_id, 
                i.get('resolution', None), 
                i['experimental_method_class'], 
                i['experimental_method'],
                multi_method) for i in exp_data[pdb_ob.pdb_id])
    
    @staticmethod
    @unsync
    async def fetch_date(pdb_ob: PDB):
        summary_data = await pdb_ob.fetch_from_web_api('api/pdb/entry/summary/', a_load_json, json=True)
        data = summary_data[pdb_ob.pdb_id][0]
        return pdb_ob.pdb_id, data['revision_date'], data['deposition_date']
    
    def fetch(self, func:Union[Callable[[PDB], List], str]):
        if isinstance(func, str):
            self.tasks = [getattr(pdb_ob, func)() for pdb_ob in self]
        else:
            self.tasks = [func(pdb_ob) for pdb_ob in self]
        return self

    @unsync
    async def run(self, tqdm=None):
        if tqdm is None:
            return [await fob for fob in as_completed(self.tasks)]
        else:
            return [await fob for fob in tqdm(as_completed(self.tasks), total=len(self.tasks))]


class SIFTSs(PDBs):
    def __init__(self, iterable: Iterable):
        super().__init__((SIFTS(i) if isinstance(i, str) else i for i in iterable))

'''
TODO: Deal with carbohydrate polymer in PISA

    * possible change in struct_asym_id
'''

# @Created Date: 2020-08-11 10:48:08 pm
# @Filename: record.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-08-11 10:48:11 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Iterable, Union, Callable, Optional, Hashable, Dict
from numpy import array, where as np_where
from pathlib import Path
from pandas import isna, concat, DataFrame
from unsync import unsync, Unfuture
from aiofiles import open as aiofiles_open
from smart_open import open as smart_open
from re import compile as re_compile
import orjson as json
from pdb_profiling.utils import init_semaphore, init_folder_from_suffix, a_read_csv, split_df_by_chain, related_dataframe, slice_series, to_interval, MMCIF2DictPlus, a_load_json
from pdb_profiling.processers.pdbe.api import ProcessPDBe, PDBeModelServer, PDBArchive, FUNCS as API_SET


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


class PDB(object):

    folder = None

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

    def register_task(self, key: Hashable, task: Unfuture):
        self.tasks[key] = task

    @classmethod
    def set_web_semaphore(cls, web_semaphore_value):
        cls.web_semaphore = init_semaphore(web_semaphore_value).result()

    @classmethod
    def get_web_semaphore(cls):
        return cls.web_semaphore

    @classmethod
    def set_db_semaphore(cls, db_semaphore_value):
        cls.db_semaphore = init_semaphore(db_semaphore_value).result()

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

    def __init__(self, pdb_id: str):
        self.check_folder()
        self.set_id(pdb_id)
        self.tasks = dict()
        self.pdb_ob = self
        self.properties_inited = False
    
    def __repr__(self):
        return f"<{self.__class__.__name__} {self.get_id()}>"
    
    def set_id(self, pdb_id: str):
        assert len(pdb_id) == 4, "Invalid PDB ID!"
        self.pdb_id = pdb_id.lower()

    def get_id(self):
        return self.pdb_id

    def set_neo4j_connection(self, api):
        pass

    def set_sqlite_connection(self, api):
        pass

    def fetch_from_web_api(self, api_suffix: str, then_func: Optional[Callable[[Unfuture], Unfuture]] = None, json: bool = False) -> Unfuture:
        assert api_suffix in API_SET, f"Invlaid API SUFFIX! Valid set:\n{API_SET}"
        task = self.tasks.get((api_suffix, then_func), None)
        if task is not None:
            return task
        
        args = dict(pdb=self.get_id(),
                    suffix=api_suffix,
                    method='get',
                    folder=next(init_folder_from_suffix(self.get_folder(), (api_suffix, ))),
                    semaphore=self.get_web_semaphore())
        if json:
            args['to_do_func'] = None
        task = ProcessPDBe.single_retrieve(**args)
        if then_func is not None:
            task = task.then(then_func)
        self.register_task((api_suffix, then_func), task)
        return task

    def fetch_from_modelServer_api(self, api_suffix: str, method: str = 'post', data_collection=None, params=None, then_func: Optional[Callable[[Unfuture], Unfuture]] = None) -> Unfuture:
        assert api_suffix in PDBeModelServer.api_sets, f"Invlaid API SUFFIX! Valid set:\n{PDBeModelServer.api_sets}"
        task = self.tasks.get((PDBeModelServer.root, api_suffix, method, data_collection, params, then_func), None)
        if task is not None:
            return task
        task = PDBeModelServer.single_retrieve(
            pdb=self.pdb_id,
            suffix=api_suffix,
            method=method,
            folder=next(init_folder_from_suffix(
                self.get_folder()/'model-server', (api_suffix, ))),
            semaphore=self.get_web_semaphore(),
            data_collection=data_collection,
            params=params)
        if then_func is not None:
            task = task.then(then_func)
        self.register_task((PDBeModelServer.root, api_suffix, method, data_collection, params, then_func), task)
        return task

    def fetch_from_PDBArchive(self, api_suffix: str, then_func: Optional[Callable[[Unfuture], Unfuture]] = None, **kwargs) -> Unfuture:
        assert api_suffix in PDBArchive.api_sets, f"Invlaid API SUFFIX! Valid set:\n{PDBArchive.api_sets}"
        task = self.tasks.get((PDBArchive.root, api_suffix, then_func), None)
        if task is not None:
            return task
        task = PDBArchive.single_retrieve(
            pdb=self.pdb_id,
            suffix=api_suffix,
            folder=next(init_folder_from_suffix(
                self.get_folder()/'pdb/data/structures', (api_suffix, ))),
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
        NOTE except waters

        Example of chain_id != struct_asym_id:

            * 6a8r
            * 6e32
        '''
        res_df = await self.fetch_from_web_api('api/pdb/entry/residue_listing/', PDB.to_dataframe)
        eec_df = res_df[['pdb_id', 'entity_id', 'chain_id',
                         'struct_asym_id']].drop_duplicates().reset_index(drop=True)
        eec_df['struct_asym_id_in_assembly'] = eec_df.struct_asym_id
        if merge_with_molecules_info:
            mol_df = await self.fetch_from_web_api('api/pdb/entry/molecules/', PDB.to_dataframe)
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

'''
TODO: Deal with carbohydrate polymer in PISA

    * possible change in struct_asym_id
'''

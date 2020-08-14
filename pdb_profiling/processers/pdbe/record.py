# @Created Date: 2020-08-11 10:48:08 pm
# @Filename: record.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-08-11 10:48:11 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Iterable, Union, Callable, Optional, Hashable
from numpy import array, where as np_where
from pathlib import Path
from pandas import isna, concat
from unsync import unsync, Unfuture
from re import compile as re_compile
import orjson as json
from pdb_profiling.utils import init_semaphore, init_folder_from_suffix, a_read_csv, split_df_by_chain, related_dataframe, slice_series, to_interval
from pdb_profiling.processers.pdbe.api import ProcessPDBe, FUNCS as API_SET


API_SET = {api for apiset in API_SET for api in apiset[1]}


class PDB(object):

    folder = None

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

    def fetch_from_web_api(self, api_suffix: str, then_func: Optional[Callable[[Unfuture], Unfuture]] = None) -> Unfuture:
        assert api_suffix in API_SET, f"Invlaid API SUFFIX! Valid set:\n{API_SET}"
        task = self.tasks.get((api_suffix, then_func), None)
        if task is not None:
            return task
        task = ProcessPDBe.single_retrieve(
            pdb=self.get_id(),
            suffix=api_suffix,
            method='get', 
            folder=next(init_folder_from_suffix(self.get_folder(), (api_suffix, ))),
            semaphore=self.get_web_semaphore())
        if then_func is not None:
            task = task.then(then_func)
        self.register_task((api_suffix, then_func), task)
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
        '''
        
        def to_assembly_id(pdb_id, assemblys):
            for assembly_id in assemblys:
                yield f"{pdb_id}/{assembly_id}"
        
        ass_eec_df = await self.fetch_from_web_api('api/pdb/entry/assembly/', self.assembly2eec)
        ass_eec_df = ass_eec_df[ass_eec_df.details.notnull()]
        assemblys = set(ass_eec_df.assembly_id) | {0}
        if focus_assembly_ids is not None:
            assemblys = sorted(assemblys & set(int(i) for i in focus_assembly_ids))
        else:
            assemblys = sorted(assemblys)
        self.assembly = dict(zip(
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
        

class PDBAssemble(PDB):

    id_pattern = re_compile(r"([a-z0-9]{4})/([0-9]+)")
    struct_range_pattern = re_compile(r"\[.+\]([A-Z]+):[0-9]+")  # e.g. [FMN]B:149

    def __init__(self, pdb_ass_id, pdb_ob: Optional[PDB]=None, interface_filters={'structure_2.symmetry_operator': ('eq', 'x,y,z'), 'css': ('gt', 0)}):
        super().__init__(pdb_ass_id)
        self.pdb_ob = pdb_ob
        self.interface_filters = interface_filters

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
        interfacelist_df = await path.then(cls.to_dataframe)
        if interfacelist_df is None:
            return None
        interfacelist_df.rename(columns={
                                "id": "interface_id", 
                                "pdb_code": "pdb_id", 
                                "assemble_code": "assembly_id"
                                }, inplace=True)
        interfacelist_df['struct_asym_id_in_assembly_1'] = interfacelist_df.apply(
            lambda x: cls.struct_range_pattern.match(x['structure_1.range']).group(1) if x['structure_1.original_range'] != '{-}' else x['structure_1.range'], axis=1)
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

        self.interface = dict(zip(
            focus_interface_ids, (PDBInterface(if_id, self) for if_id in to_interface_id(self.get_id(), focus_interface_ids))))

    def get_interface(self, interface_id):
        return self.interface[interface_id]

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
        
        eec_as_df = await self.pdb_ob.get_eec_as_df()
        protein_type_asym = eec_as_df[
            eec_as_df.assembly_id.eq(self.assembly_id) &
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

        eec_as_df = await self.pdb_ob.get_eec_as_df()
        molType_dict = eec_as_df[eec_as_df.assembly_id.eq(self.assembly_id)].groupby('molecule_type').struct_asym_id_in_assembly.apply(set).to_dict()

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
        eec_as_df = await self.pdbAssemble_ob.pdb_ob.get_eec_as_df()
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
            cur_record = eec_as_df[eec_as_df.struct_asym_id_in_assembly.eq(saiia2) &
                      eec_as_df.assembly_id.eq(self.assembly_id)][cur_keys].to_dict('record')[0]
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

# @Created Date: 2020-08-11 10:48:08 pm
# @Filename: record.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-08-11 10:48:11 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Iterable, Union, Callable, Optional, Hashable
from pathlib import Path
from unsync import unsync, Unfuture
from re import compile as re_compile
from pdb_profiling.utils import init_semaphore, init_folder_from_suffix, a_read_csv, split_df_by_chain, related_dataframe
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
        residue_df = await a_read_csv(path, sep="\t", converters=ProcessPDBe.converters)
        return residue_df

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
    async def res2eec(self, merge_with_molecules_info:bool=True):
        '''
        NOTE except waters
        '''
        res_df = await self.fetch_from_web_api('api/pdb/entry/residue_listing/', PDB.to_dataframe)
        eec_df = res_df[['pdb_id', 'entity_id', 'chain_id',
                         'struct_asym_id']].drop_duplicates().reset_index(drop=True)
        if merge_with_molecules_info:
            mol_df = await self.fetch_from_web_api('api/pdb/entry/molecules/', PDB.to_dataframe)
            return eec_df.merge(mol_df, how='left')
        else:
            return eec_df

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
    async def set_interface(self, obligated_class_chains: Optional[Iterable[str]] = None):

        def to_interface_id(pdb_assembly_id, focus_interface_ids):
            for interface_id in focus_interface_ids:
                yield f"{pdb_assembly_id}/{interface_id}"

        interfacelist_df = await self.fetch_from_web_api('api/pisa/interfacelist/', self.to_interfacelist_df)
        
        if interfacelist_df is None:
            self.interface = dict()
            return
        
        if obligated_class_chains is not None:
            interfacelist_df = interfacelist_df[
                interfacelist_df.struct_asym_id_in_assembly_1.isin(obligated_class_chains) |
                interfacelist_df.struct_asym_id_in_assembly_2.isin(obligated_class_chains)]
        
        focus_interface_ids = related_dataframe(
            self.interface_filters, interfacelist_df).interface_id.unique()

        self.interface = dict(zip(
            focus_interface_ids, (PDBInterface(if_id, self) for if_id in to_interface_id(self.get_id(), focus_interface_ids))))

    def get_interface(self, interface_id):
        return self.interface[interface_id]

    @property
    def demo_interface_filters(self):
        return '''
        ass_eec_df = await self.fetch_from_web_api('api/pdb/entry/assembly/', PDB.assembly2eec)
        molType_label_dict = ass_eec_df[ass_eec_df.assembly_id.eq(self.assembly_id)].groupby(
            'molecule_type').struct_asym_id_in_assembly.apply(tuple).to_dict()
        '''
        


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





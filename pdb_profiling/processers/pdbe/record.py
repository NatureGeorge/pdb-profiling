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
        return cls.folder

    def __init__(self, pdb_id: str, folder: Union[Path, str]):
        self.set_pdb_id(pdb_id)
        self.set_folder(folder)
        self.tasks = dict()
    
    def __repr__(self):
        return f"<{self.__class__.__name__} {self.pdb_id}>"
    
    def set_pdb_id(self, pdb_id: str):
        assert len(pdb_id) == 4, "Invalid PDB ID!"
        self.pdb_id = pdb_id.lower()

    def get_pdb_id(self):
        return (self.pdb_id, )

    def set_neo4j_connection(self, api):
        pass

    def set_sqlite_connection(self, api):
        pass

    def fetch_from_web_api(self, api_suffix: str, then_func: Optional[Callable[[Unfuture], Unfuture]] = None) -> Unfuture:
        assert api_suffix in API_SET, f"Invlaid API SUFFIX! Valid set:\n{API_SET}"
        task = self.tasks.get((api_suffix, then_func), None)
        if task is not None:
            return task
        task = ProcessPDBe.retrieve(
            pdbs=self.get_pdb_id(),
            suffix=api_suffix,
            method='get', 
            folder=next(init_folder_from_suffix(self.get_folder(), (api_suffix, ))),
            ret_res=False, 
            semaphore=self.get_web_semaphore())[0]
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
            * polydeoxyribonucleotide
            * polyribonucleotide
            * polydeoxyribonucleotide/polyribonucleotide hybrid (e.g. 6ozg)
        '''
        eec_df = await path.then(cls.to_dataframe)
        eec_df['struct_asym_id_in_assembly'] = eec_df.struct_asym_id
        eec_df = eec_df.drop(columns=['number_residues', 'observed_ratio', 'partial_ratio'])
        return eec_df

    @unsync
    async def res2eec(self, merge_with_molecules_info:bool=False):
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
                              'assembly_id')
        eec_df = split_df_by_chain(
            assembly_df[assembly_df.molecule_type.ne('water')],
            assembly_focus_col, assembly_focus_col[0:1],
            'json-list').rename(columns={"in_chains": "struct_asym_id_in_assembly"}).reset_index(drop=True)
        return eec_df

    @unsync
    async def set_assembly(self, focus_assembly_ids:Optional[Iterable[int]]=None):
        
        def to_assembly_id(pdb_id, assemblys):
            for assembly_id in assemblys:
                yield f"{pdb_id}/{assembly_id}"
        
        ass_eec_df = await self.fetch_from_web_api('api/pdb/entry/assembly/', self.assembly2eec)
        assemblys = set(ass_eec_df.assembly_id) | {0}
        if focus_assembly_ids is not None:
            assemblys = sorted(assemblys & set(int(i) for i in focus_assembly_ids))
        else:
            assemblys = sorted(assemblys)
        self.assembly = dict(zip(
            assemblys, 
            (PDBAssemble(ass_id, self.get_folder()) for ass_id in to_assembly_id(self.pdb_id, assemblys))))

    def get_assembly(self, assembly_id):
        return self.assembly[assembly_id]


class PDBAssemble(PDB):

    id_pattern = re_compile(r"[a-z0-9]{4}/[0-9]+")
    interface_filters = {
        'structure_2.symmetry_operator': ('eq', 'x,y,z'),
        'css': ('gt', 0),
        }

    def set_pdb_id(self, pdb_id: str):
        self.pdb_id = pdb_id.lower()
        assert bool(self.id_pattern.fullmatch(self.pdb_id)), f"Invalid ID: {self.pdb_id}"

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
        return interfacelist_df

    @unsync
    async def set_interface(self):
        
        def to_interface_id(pdb_assembly_id, interfaces):
            for interface_id in interfaces:
                yield f"{pdb_assembly_id}/{interface_id}"

        interfacelist_df = await self.fetch_from_web_api('api/pisa/interfacelist/', self.to_interfacelist_df)
        
        if interfacelist_df is None:
            self.interface = dict()
            return
        
        interfaces = related_dataframe(
            self.interface_filters, interfacelist_df).interface_id.unique()
        self.interface = dict(zip(
            interfaces,
            (PDBInterface(if_id, self.get_folder()) for if_id in to_interface_id(self.pdb_id, interfaces))))


class PDBInterface(PDBAssemble):
 
    id_pattern = re_compile(r"[a-z0-9]{4}/[0-9]+/[0-9]+")







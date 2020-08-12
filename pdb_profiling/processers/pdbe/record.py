# @Created Date: 2020-08-11 10:48:08 pm
# @Filename: record.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-08-11 10:48:11 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Iterable, Union
from unsync import unsync, Unfuture
from re import compile
from pdb_profiling.utils import init_semaphore

class PDB(object):

    @classmethod
    def set_web_semaphore(cls, web_semaphore_value):
        cls.web_semaphore = init_semaphore(web_semaphore_value).result()
    
    @classmethod
    def set_db_semaphore(cls, db_semaphore_value):
        cls.db_semaphore = init_semaphore(db_semaphore_value).result()

    def __init__(self, pdb_id: str):
        self.set_pdb_id(pdb_id)
    
    def __repr__(self):
        return f"<{self.__class__.__name__} {self.pdb_id}>"
    
    def set_pdb_id(self, pdb_id: str):
        assert len(pdb_id) == 4, "Invalid PDB ID!"
        self.pdb_id = pdb_id.lower()

    def set_neo4j_connection(self, api):
        pass

    def set_sqlite_connection(self, api):
        pass

    @staticmethod
    def submit_task(task_func: Unfuture, **task_kwargs) -> Unfuture:
        return task_func(**task_kwargs)


class PDBAssemble(PDB):

    id_pattern = compile("[a-z0-9]{4}/[0-9]+")

    def set_pdb_id(self, pdb_id: str):
        self.pdb_id = pdb_id.lower()
        assert bool(self.id_pattern.fullmatch(self.pdb_id)), f"Invalid ID: {self.pdb_id}"


class PDBInterface(PDB):
    
    id_pattern = compile("[a-z0-9]{4}/[0-9]+/[0-9]+")

    def set_pdb_id(self, pdb_id: str):
        self.pdb_id = pdb_id.lower()
        assert bool(self.id_pattern.fullmatch(self.pdb_id)), f"Invalid ID: {self.pdb_id}"


class PDBCollection(PDB):
    def __init__(self, pdbs: Iterable[PDB]):
        self.pdbs = pdbs

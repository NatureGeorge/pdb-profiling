# @Created Date: 2021-04-22 07:23:06 pm
# @Filename: recordbase.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2021-04-22 09:21:16 pm
# @Copyright (c) 2021 MinghuiGroup, Soochow University
from dataclasses import dataclass
from collections import OrderedDict
from re import compile as re_compile

ID_SUFFIX = r'[0-9]+)[\.]*([0-9]*)'
PDB_COMMON_PREFIX = r'^((?=.*[A-Za-z])(?=.*\d)[A-Za-z\d]'
#RANGE_SUFFIX = r':([0-9]+)-*([0-9]*)$'


@dataclass
class IdentifierBase:
    PATS = OrderedDict({
        re_compile(f'(NC_{ID_SUFFIX}'): ('RefSeq', 'genome'),
        re_compile(r'([NX]{1}M_'+ID_SUFFIX): ('RefSeq', 'transcript'),
        re_compile(r'([NX]{1}P_'+ID_SUFFIX): ('RefSeq', 'protein'),
        re_compile(f'(ENS[A-Z]*E{ID_SUFFIX}'): ('Ensembl', 'exon'),
        re_compile(f'(ENS[A-Z]*G{ID_SUFFIX}'): ('Ensembl', 'gene'),
        re_compile(f'(ENS[A-Z]*T{ID_SUFFIX}'): ('Ensembl', 'transcript'),
        re_compile(f'(ENS[A-Z]*P{ID_SUFFIX}'): ('Ensembl', 'protein'),
        re_compile(r'(CCDS[0-9]+)'): ('CCDS', 'CCDS'),
        re_compile(r'(rs[0-9]+)'): ('dbSNP', 'mutation'),
        re_compile(r'^((?=.*[A-Za-z])(?=.*\d)[A-Za-z\d]{6,})[\-]*([0-9]*)$'): ('UniProt', 'isoform'),
        re_compile(PDB_COMMON_PREFIX+r'{4})$'): ('PDB', 'entry'),
        re_compile(PDB_COMMON_PREFIX+r'{4})-([0-9]+)$'): ('PDB', 'assembly'),
        re_compile(PDB_COMMON_PREFIX+r'{4})_([0-9]+)$'): ('PDB', 'entity'),
        re_compile(PDB_COMMON_PREFIX+r'{4})\.([A-Z]+)$'): ('PDB', 'instance'),
        re_compile(PDB_COMMON_PREFIX+r'{4})/([0-9/]+)$'): ('PDB', 'entry_like'),
        re_compile(r'(PDB-CPX-[0-9]+)'): ('PDB', 'complex'),
        re_compile(r'([0-9]+)/([0-9]+)'): ('Taxonomy', 'genome'),
        re_compile(r'([0-9]+)'): ('HGNC', 'HGNC'),
        #re_compile(r''): ('EC', 'EC'),
        re_compile(r'([A-z0-9\-]+)'): ('PDB', 'compounds'),
    })

    raw_identifier: str
    renew:bool = True
    source: str = ''
    level: str = ''
    identifier: str = ''
    identifier_suffix: str = ''
    
    @classmethod
    def get_type(cls, raw_identifier: str):
        for pat, group in cls.PATS.items():
            res = pat.fullmatch(raw_identifier)
            if bool(res):
                return group, res.groups()
        raise AssertionError(f"Unexpected identifier type: {raw_identifier}")

    def __post_init__(self):
        if self.renew:
            (self.source, self.level), identifier_tp = self.get_type(self.raw_identifier)
            if len(identifier_tp) == 2:
                self.identifier, self.identifier_suffix = identifier_tp
            else:
                self.identifier, = identifier_tp

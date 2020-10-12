# @Created Date: 2020-05-13 08:54:03 pm
# @Filename: __init__.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-05-13 08:54:09 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from re import compile as re_compile

__version__ = '0.1.9'


common_pat = r'^(?=.*[A-Za-z])(?=.*\d)[A-Za-z\d]'


pats = dict(pdb_id=re_compile(common_pat+r'{4}$'),
            pdb_entity_id=re_compile(common_pat+r'{4}_[0-9]+$'),
            UniProt=re_compile(common_pat+r'{6,}[\-]*[0-9]*$'),
            pdb_complex_id=re_compile(r'PDB-CPX-[0-9]+'))


def default_id_tag(identifier:str, default:str='', raise_error:bool=False):
    for pat_name, pat in pats.items():
        if bool(pat.fullmatch(identifier)):
            return pat_name
    if raise_error:
        raise ValueError(f'Unexcepted Identifiers: {identifier}')
    else:
        return default


def default_config(folder='./'):
    from pdb_profiling.log import Abclog
    from pdb_profiling.fetcher.webfetch import UnsyncFetch
    from pdb_profiling.processors.pdbe.record import Base, SIFTS
    from pdb_profiling.processors.pdbe.api import ProcessPDBe
    from pdb_profiling.processors.proteins.record import Identifier
    from pdb_profiling.processors import UniProtFASTA
    # Use Existing Handled PDBe API Results (e.g. tsv format results)
    ProcessPDBe.use_existing = True
    # Use Existing API Results (e.g. json format results downloaded from web)
    UnsyncFetch.use_existing = True
    # Init Abclog Logger
    Abclog.init_logger(logName='PDB-Profiling')
    # Set WebFetcher's Semaphore
    Base.set_web_semaphore(30).result()
    Identifier.set_web_semaphore(30).result()
    UniProtFASTA.set_web_semaphore(30).result()
    # Set Folder that store downloaded and handled files
    Base.set_folder(folder)
    SIFTS.set_folder(folder)
    Identifier.set_folder(folder)
    UniProtFASTA.set_folder(folder)

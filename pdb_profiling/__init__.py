# @Created Date: 2020-05-13 08:54:03 pm
# @Filename: __init__.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-05-13 08:54:09 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from re import compile as re_compile

__version__ = '0.1.6'


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


def default_config():
    from pdb_profiling.log import Abclog
    from pdb_profiling.fetcher.webfetch import UnsyncFetch
    from pdb_profiling.processors.pdbe.api import ProcessPDBe
    from pdb_profiling.processors.pdbe.record import Base, PDBeModelServer, PDBArchive
    from pdb_profiling.processors.uniprot.api import UniProtFASTA
    # Use Existing Handled PDBe API Results (e.g. tsv format results)
    ProcessPDBe.use_existing = True
    # Init Abclog Logger
    Abclog.init_logger(logName='PDB-Profiling')
    # Init ProcessPDBe's Logger (pass it with Abclog Logger)
    ProcessPDBe.init_logger(logger=Abclog.logger)
    # Use Existing API Results (e.g. json format results downloaded from web)
    UnsyncFetch.use_existing = True
    # Init WebFetcher's Logger (pass it with Abclog Logger)
    UnsyncFetch.init_setting(ProcessPDBe.logger)
    # Set WebFetcher's Semaphore
    Base.set_web_semaphore(30)
    # Set Folder that store downloaded and handled files
    Base.set_folder('./')
    # Init ModelServer API's Logger (pass it with Abclog Logger)
    PDBeModelServer.init_logger(logger=ProcessPDBe.logger)
    # Init PDBArchive API's Logger (pass it with Abclog Logger)
    PDBArchive.init_logger(logger=ProcessPDBe.logger)
    # Init UniProtFASTA API's Logger (pass it with Abclog Logger)
    UniProtFASTA.init_logger(logger=ProcessPDBe.logger)

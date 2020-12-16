# @Created Date: 2020-05-13 08:54:03 pm
# @Filename: __init__.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-05-13 08:54:09 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
__version__ = '0.2.5'


def default_config(folder='./'):
    from pdb_profiling.log import Abclog
    from pdb_profiling.fetcher.webfetch import UnsyncFetch
    from pdb_profiling.processors.pdbe.record import Base, PDB
    from pdb_profiling.processors.pdbe.api import ProcessPDBe
    from pdb_profiling.processors.proteins.record import Identifier
    from pdb_profiling.processors.uniprot.api import UniProtFASTA, UniProtAPI
    from pdb_profiling.processors.uniprot.record import UniProts
    from pdb_profiling.processors.i3d.api import Interactome3D
    from pdb_profiling.processors.swissmodel.api import SMR
    # Use Existing Handled PDBe API Results (e.g. tsv format results)
    ProcessPDBe.use_existing = True
    # Use Existing API Results (e.g. json format results downloaded from web)
    UnsyncFetch.use_existing = True
    # Init Abclog Logger
    Abclog.init_logger(logName='PDB-Profiling')
    # Set WebFetcher's Semaphore
    Base.set_web_semaphore(50).result()
    Identifier.set_web_semaphore(30).result()
    UniProtFASTA.set_web_semaphore(30).result()
    UniProtAPI.set_web_semaphore(30).result()
    Interactome3D.set_web_semaphore(30).result()
    SMR.set_web_semaphore(30).result()
    # Set Folder that store downloaded and handled files
    Base.set_folder(folder)
    PDB.set_folder(folder)
    Identifier.set_folder(folder)
    UniProtFASTA.set_folder(folder)
    UniProtAPI.set_folder(folder)
    UniProts.set_folder(folder)
    Interactome3D.set_folder(folder)
    SMR.set_folder(folder)

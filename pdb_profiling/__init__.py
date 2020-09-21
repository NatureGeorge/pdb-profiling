# @Created Date: 2020-05-13 08:54:03 pm
# @Filename: __init__.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-05-13 08:54:09 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
__version__ = '0.1.5'


def default_config():
    from pdb_profiling.fetcher.webfetch import UnsyncFetch
    from pdb_profiling.processers.pdbe.api import ProcessPDBe
    from pdb_profiling.processers.pdbe.record import PDB, PDBeModelServer, PDBArchive
    # Use Existing Handled PDBe API Results (e.g. tsv format results)
    ProcessPDBe.use_existing = True
    # Init PDBe API Logger
    ProcessPDBe.init_logger()
    # Use Existing API Results (e.g. json format results downloaded from web)
    UnsyncFetch.use_existing = True
    # Init WebFetcher's Logger (pass it with PDBe API Logger)
    UnsyncFetch.init_setting(ProcessPDBe.logger)
    # Set WebFetcher's Semaphore
    PDB.set_web_semaphore(30)
    # Set Folder that store downloaded and handled files
    PDB.set_folder('./')
    # Init ModelServer API's Logger (pass it with PDBe API Logger)
    PDBeModelServer.init_logger(logger=ProcessPDBe.logger)
    # Init PDBArchive API's Logger (pass it with PDBe API Logger)
    PDBArchive.init_logger(logger=ProcessPDBe.logger)

from pdb_profiling import default_config
from tqdm import tqdm

default_config()

def test_init():
    from pdb_profiling.processors.i3d.api import Interactome3D
    Interactome3D.pipe_init_interaction_meta().result()

def test_single_select():
    from pdb_profiling.processors import SIFTS
    # SIFTS.chain_filter, SIFTS.entry_filter = '', ''
    demo = SIFTS('P21359')
    demo.pipe_select_mo().result()
    demo.pipe_select_ho(run_as_completed=True, progress_bar=tqdm).result()
    demo.pipe_select_he(run_as_completed=True, progress_bar=tqdm).result()

def test_identifiers():
    from pdb_profiling.processors import Identifiers
    demo = Identifiers([
        'ENSP00000491589', 'ENST00000379268',
        'ENSP00000427757', 'ENSP00000266732',
        'NP_001291289.1', 'ENST00000335295',
        'NP_001165602.1', 'ENST00000402254',
        'ENST00000371100'])
    demo.fetch('map2unp').run().result()

def test_uniprots_alt():
    from pdb_profiling.processors import UniProts, Identifiers
    from pdb_profiling.utils import a_concat
    UniProts.fetch_VAR_SEQ_from_DB(('Q5VST9', 'Q5JWF2', 'P08631', 'O92972')).result()
    
    demo_unps = ('Q5VST9', 'Q5JWF2', 'P21359', 'P68871', 'P63092')
    Identifiers(demo_unps).query_from_DB_with_unps('ALTERNATIVE_PRODUCTS').run().then(a_concat).result()

def test_command():
    from os import system
    system('pdb_profiling --help')
    system('pdb_profiling insert-mutation --help')
    system('pdb_profiling id-mapping --help')
    system('pdb_profiling sifts-mapping --help')

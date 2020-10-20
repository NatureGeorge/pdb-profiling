from pdb_profiling import default_config

default_config()

def test_init():
    from pdb_profiling.processors.i3d.api import Interactome3D
    Interactome3D.pipe_init_interaction_meta().result()

def test_single_select():
    from pdb_profiling.processors import SIFTS
    SIFTS.chain_filter, SIFTS.entry_filter = '', ''
    demo = SIFTS('O15350')
    demo.pipe_select_mo().result()
    demo.pipe_select_ho().result()
    demo.pipe_select_he().result()

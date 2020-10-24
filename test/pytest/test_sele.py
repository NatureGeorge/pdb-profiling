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

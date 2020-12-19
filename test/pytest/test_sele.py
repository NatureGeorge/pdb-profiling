from pdb_profiling import default_config
from rich.progress import track

default_config()

def test_init():
    from pdb_profiling.processors.i3d.api import Interactome3D
    Interactome3D.pipe_init_interaction_meta().result()

def test_single_select():
    from pdb_profiling.processors import SIFTS
    # SIFTS.chain_filter, SIFTS.entry_filter = '', ''
    demo = SIFTS('P21359')
    demo.pipe_select_mo().result()
    demo.pipe_select_ho(run_as_completed=True, progress_bar=track).result()
    demo.pipe_select_he(run_as_completed=True, progress_bar=track).result()

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
    
    demo_unps = ('Q5VST9', 'Q5JWF2', 'P21359', 'P68871', 'P63092', 'Q29960')
    Identifiers(demo_unps).query_from_DB_with_unps('ALTERNATIVE_PRODUCTS').run().then(a_concat).result()

def test_command():
    from os import system
    system('pdb_profiling --help')
    system('pdb_profiling insert-mutation --help')
    system('pdb_profiling id-mapping --help')
    system('pdb_profiling sifts-mapping --help')

def test_other_api():
    from pdb_profiling.processors import PDB
    from pdb_profiling.processors.pdbe.api import PDBVersioned
    pdb_ob = PDB('1a01')
    pdb_ob.fetch_from_pdbe_api('api/pdb/entry/secondary_structure/').result()
    pdb_ob.fetch_from_pdbe_api('api/pdb/entry/files/').result()
    pdb_ob.fetch_from_pdbe_api('graph-api/pdb/funpdbe_annotation/').result()
    pdb_ob.fetch_from_pdbe_api('graph-api/pdb/sequence_conservation/').result()
    PDB('4zai').fetch_from_PDBArchive('obsolete/mmCIF/', PDB.cif2residue_listing).result()
    pv_path = PDB.get_folder()/'pdb-versioned/entries'
    pv_path.mkdir(parents=True, exist_ok=True)
    PDBVersioned.single_retrieve(('4fc3', '_v1-2'), 'entries/', pv_path, PDB.get_web_semaphore(), file_suffix='.cif.gz').result()
    bm_df = pdb_ob.get_bound_molecules().result()
    [pdb_ob.get_bound_molecule_interaction(bm_id).result() for bm_id in bm_df.bm_id.unique()[:2]]

def test_fetch_residue_mapping():
    from pdb_profiling.processors import SIFTS
    pdb_ob = SIFTS('1a01')
    pdb_ob.fetch_residue_mapping(entity_id=1, start=20, end=25).result()
    pdb_ob.fetch_residue_mapping(entity_id=1, start=24, end=27).result()

def test_rcsb_data_api():
    from pdb_profiling.processors import PDB, PDBAssemble
    from pdb_profiling.utils import a_load_json
    pdb_id = '3hl2'
    ob = PDB(pdb_id)
    assembly_ids = ob.fetch_from_rcsb_data_api(
        'graphql',
        query='{entry(entry_id:"%s"){rcsb_entry_container_identifiers{assembly_ids}}}' % pdb_id,
        then_func=a_load_json,
        json=True).result()['data']['entry']['rcsb_entry_container_identifiers']['assembly_ids']

    for assembly_id in assembly_ids:
        data = PDBAssemble(f'{pdb_id}/{assembly_id}').fetch_from_rcsb_data_api('assembly/', then_func=a_load_json, json=True).result()
        data['pdbx_struct_assembly_gen']
        data['pdbx_struct_oper_list']

    df1 = ob.rd_source_ass_oper_df().result()
    df2 = ob.ms_source_ass_oper_df(**ob.pipe_assg_data_collection().result()).result()
    assert df1.shape == df2.shape
    assert (df1.merge(df2).shape) == df1.shape

def test_get_sequence():
    from pdb_profiling.processors import PDB
    ob = PDB('4u2v')
    ob.get_sequence(entity_id=1).result()
    ob.get_sequence(mode='raw_seq', entity_id=1).result()
    ob.get_sequence(mode='raw_pdb_seq', entity_id=1).result()
    ob.get_sequence(mode='mod_x_seq', entity_id=1).result()
    ob.get_sequence(struct_asym_id='A').result()
    ob.get_sequence(chain_id='A').result()

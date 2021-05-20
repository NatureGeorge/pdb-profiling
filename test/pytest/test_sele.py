# @Created Date: 2020-10-20 10:00:38 am
# @Filename: test_sele.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2021-03-15 09:07:11 pm
# @Copyright (c) 2021 MinghuiGroup, Soochow University
from pdb_profiling import default_config
from pdb_profiling.utils import a_load_json, a_concat
from pdb_profiling.processors import *
from rich.progress import track
from pandas import DataFrame
import pytest

default_config('test/pytest/demo_dir')


@pytest.mark.timeout(60)
def test_init():
    from pdb_profiling.processors.i3d.api import Interactome3D
    Interactome3D.pipe_init_interaction_meta().result()


@pytest.mark.timeout(300)
def test_single_select():
    # SIFTS.chain_filter, SIFTS.entry_filter = '', ''
    demo = SIFTS('P21359-2')
    demo.pipe_base().then(SIFTS.double_check_conflict_and_range).result()
    demo.pipe_scheduled_ranged_map_res_df().result()
    demo.pipe_select_he(run_as_completed=True, progress_bar=track).result()
    demo.pipe_select_ho_iso(run_as_completed=True).result()
    demo.pipe_select_else(func='pipe_protein_ligand_interface', css_cutoff=0.5, run_as_completed=True).result()


@pytest.mark.timeout(90)
def test_identifiers():
    demo = Identifiers([
        'ENSP00000491589', 'ENST00000379268',
        'ENSP00000427757', 'ENSP00000266732',
        'NP_001291289.1', 'ENST00000335295',
        'NP_001165602.1', 'ENST00000402254',
        'ENST00000371100', 'ENST00000401731',
        'ENSP00000387612'])
    demo.fetch('map2unp').run().result()
    Identifier('P21359-3').fetch_from_proteins_api('coordinates/location/', ':550').result()
    Identifier('P21359-2').init().result().get_isoform_ob().result()
    Identifier('P21359').get_all_ref_identifiers().result()
    Identifier('P21359').alignment_df().result()


@pytest.mark.timeout(60)
def test_uniprots_alt():
    UniProts.fetch_VAR_SEQ_from_DB(('Q5VST9', 'Q5JWF2', 'P08631', 'O92972'), via_txt=True).result()
    
    demo_unps = ('Q5VST9', 'Q5JWF2', 'P21359', 'P68871', 'P63092', 'Q29960')
    Identifiers(demo_unps).query_from_DB_with_unps('ALTERNATIVE_PRODUCTS').run().then(a_concat).result()


@pytest.mark.timeout(120)
def test_other_api():
    from pdb_profiling.processors.pdbe.api import PDBVersioned, PDBeKBAnnotations
    pdb_ob = PDB('1a01')
    pdb_ob.status
    pdb_ob.summary
    pdb_ob.stats_chain().result()
    SIFTS('1a01').get_oligo_state().result()
    pdb_ob.fetch_from_pdbe_api('api/pdb/entry/secondary_structure/').result()
    pdb_ob.fetch_from_pdbe_api('api/pdb/entry/files/').result()
    pdb_ob.fetch_from_pdbe_api('graph-api/pdb/funpdbe_annotation/').result()
    pdb_ob.fetch_from_pdbe_api('graph-api/pdb/sequence_conservation/', mask_id='1cbs/1').result()
    pdb_ob.fetch_from_pdbe_api('api/validation/RNA_pucker_suite_outliers/entry/').result()
    pdb_ob.fetch_from_pdbe_api('api/validation/rama_sidechain_listing/entry/').result()
    PDB('4zai').fetch_from_PDBArchive('obsolete/mmCIF/', PDB.cif2residue_listing).result()
    pv_path = PDB.get_folder()/'pdb-versioned/entries'
    pv_path.mkdir(parents=True, exist_ok=True)
    PDBVersioned.single_retrieve(('4fc3', '_v1-2'), 'entries/', pv_path, PDB.get_web_semaphore(), file_suffix='.cif.gz').result()
    bm_df = pdb_ob.get_bound_molecules().result()
    [pdb_ob.get_bound_molecule_interaction(bm_id).result() for bm_id in bm_df.bm_id.unique()[:2]]
    SIFTS('P21359-2').fetch_from_pdbe_api('graph-api/uniprot/superposition/', SIFTS.to_dataframe).result()
    PDBAssembly('1a01/1').add_args().assembly_summary


@pytest.mark.timeout(80)
def test_pdbekdb_self_annotation():
    """from pdb_profiling.processors.pdbe.api import PDBeKBAnnotations
    PDBeKBAnnotations.root = PDBeKBAnnotations.ftp_root
    assert PDB('12ca').pipe_pdbekb_annotations('MetalPDB/').result() is not None
    PDBeKBAnnotations.root = PDBeKBAnnotations.http_root"""
    assert SIFTS('P69327-2').get_mapped_pdbekb_annotaions('M-CSA/', only_protein_res=True).result() is not None


@pytest.mark.timeout(60)
def test_fetch_residue_mapping():
    pdb_ob = SIFTS('3pg7')
    pdb_ob.fetch_residue_mapping(entity_id=1, start=251, end=256).result()
    pdb_ob.fetch_residue_mapping(entity_id=1, start=252, end=255).result()


@pytest.mark.timeout(80)
def test_rcsb_data_api():
    pdb_id = '3hl2'
    ob = PDB(pdb_id)
    assembly_ids = ob.fetch_from_rcsb_api(
        'graphql',
        query='{entry(entry_id:"%s"){rcsb_entry_container_identifiers{assembly_ids}}}' % pdb_id,
        then_func=a_load_json,
        json=True).result()['data']['entry']['rcsb_entry_container_identifiers']['assembly_ids']

    for assembly_id in assembly_ids:
        data = PDBAssembly(f'{pdb_id}/{assembly_id}').fetch_from_rcsb_api('assembly/', then_func=a_load_json, json=True).result()
        data['pdbx_struct_assembly_gen']
        data['pdbx_struct_oper_list']

    df1 = ob.rd_source_ass_oper_df().result()
    df2 = ob.ms_source_ass_oper_df(**ob.pipe_assg_data_collection().result()).result()
    assert df1.shape == df2.shape
    assert (df1.merge(df2).shape) == df1.shape


@pytest.mark.timeout(40)
def test_rcsb_cluster_membership():
    PDB('2d4q').rcsb_cluster_membership(entity_id=1, identity_cutoff=100).result()
    PDB('2e2x').rcsb_cluster_membership(entity_id=1, identity_cutoff=100).result()


@pytest.mark.timeout(40)
def test_other_SIFTS_func():
    try:
        SIFTS('P21359').fetch_from_pdbe_api('api/mappings/all_isoforms/'
            ).then(SIFTS.to_dataframe
            #).then(SIFTS.check_pdb_status
            ).then(SIFTS.check_identity
            ).then(SIFTS.reformat
            ).then(SIFTS.deal_with_identical_entity_seq).result()
    except Exception:
        pass


@pytest.mark.timeout(30)
def test_get_sequence():
    ob = PDB('4u2v')
    ob.get_sequence(entity_id=1).result()
    ob.get_sequence(mode='raw_seq', entity_id=1).result()
    ob.get_sequence(mode='raw_pdb_seq', entity_id=1).result()
    ob.get_sequence(mode='mod_x_seq', entity_id=1).result()
    ob.get_sequence(struct_asym_id='A').result()
    ob.get_sequence(chain_id='A').result()


@pytest.mark.timeout(20)
def test_show_rcsb_error():
    #assert RCSB1DCoordinates('6OB3.B').alignment_df('NCBI_GENOME').result() is not None
    assert RCSB1DCoordinates('P21359').alignment_df('PDB_INSTANCE').result() is not None

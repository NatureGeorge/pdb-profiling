# @Created Date: 2021-03-13 05:52:56 pm
# @Filename: test_command.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2021-03-13 06:43:06 pm
# @Copyright (c) 2021 MinghuiGroup, Soochow University
from click.testing import CliRunner
from pdb_profiling.commands.command import Interface

def test_command1():
    runner = CliRunner()
    dargs = ['--folder', 'test/pytest/demo_dir']
    for task in ('insert-mutation --input test/pytest/data/mutation.tsv --usecols Alt,Pos,Ref,ftId',
                 'id-mapping --auto_assign',
                 'check-muta-conflict',
                 'sifts-mapping --chunksize 15',
                 'insert-sele-mapping --input test/pytest/demo_dir/pipe_select_mo.tsv',
                 'residue-mapping --input test/pytest/demo_dir/pipe_select_mo.tsv',
                 'export-mutation-mapping -o e_resmap.tsv --sele',
                 'insert-sele-mutation-mapping -i test/pytest/demo_dir/e_resmap.tsv',
                 'sifts-mapping --func pipe_select_smr_mo --chunksize 10',
                 'insert-smr-mapping -i test/pytest/demo_dir/pipe_select_smr_mo.tsv',
                 'export-smr-mutation-mapping -o e_smr_resmap.tsv --sele',
                 'insert-sifts-meta --input test/pytest/data/pdb_demo.tsv --api_suffix api/mappings/pfam/',
                 'insert-isoform-range',
                 ['sifts-mapping', '--func', 'unp_is_canonical_with_id', '--input', 'test/pytest/data/unp_demo.tsv', '--iteroutput'],
                 ['fetch1pdb', '-i', '4hho', '-a', 'atoms', '-t', 'B_ONLY',
                  '-p', 'label_asym_id=B', '-p', 'copy_all_categories=false'],
                 ['fetch1pdb', '-i', '3hl2', '-a', 'atoms',
                 '-p', 'encoding=bcif',
                 '-d', 'label_asym_id=A,label_seq_id=90',
                 '-d', 'label_asym_id=B,label_seq_id=100',
                 '-m', 'post',
                 '-t', 'A_90_B_10'],
                 'sifts-mapping --func pipe_select_ho --chunksize 5',
                 'insert-interaction -i test/pytest/demo_dir/pipe_select_ho.tsv',
                 'export-interaction-mapping -o e_interaction_resmap.tsv',
                 ):
        result = runner.invoke(Interface, dargs+task.split(' ') if not isinstance(task, list) else dargs+task)
        assert result.exit_code == 0, str(task)


def test_command2():
    runner = CliRunner()
    dargs = ['--folder', 'test/pytest/demo_dir', '--custom_db', 'pdb2unp_unp2pdb.db']
    for task in ("insert-pdb-mutation -i test/pytest/data/pdb2unp.csv --sep , --usecols pdb_id,chain_id,Ref,Alt,author_residue_number,author_insertion_code --auth",
                 "sifts-mapping --func pipe_base --autotype from_PDBAuthMutation -o task_20210602_pdb2unp_SIFTS-MAPPING.tsv",
                 "residue-mapping -i test/pytest/demo_dir/task_20210602_pdb2unp_SIFTS-MAPPING.tsv",
                 "export-pdb-mutation-mapping --auth -o task_20210602_pdb2unp_RESIDUE-MAPPING.tsv",
                 "insert-mutation -i test/pytest/data/unp2pdb.csv --sep , --usecols ftId,Ref,Pos,Alt id-mapping check-muta-conflict",
                 "sifts-mapping --func pipe_base -o task_20210602_unp2pdb_SIFTS-MAPPING.tsv residue-mapping -i test/pytest/demo_dir/task_20210602_unp2pdb_SIFTS-MAPPING.tsv",
                 #pypath=$(python -c 'from distutils.sysconfig import get_python_lib; print(get_python_lib())')
                 #sqlite3 - header - cmd ".mode tabs" ./local_db/task_20210602.db < $pypath/pdb_profiling/sql/export_mutation_mapping_all.sql > task_20210602_unp2pdb_RESIDUE-MAPPING.tsv
                 ):
        result = runner.invoke(Interface, dargs+task.split(' '))
        assert result.exit_code == 0, str(task)

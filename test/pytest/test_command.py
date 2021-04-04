# @Created Date: 2021-03-13 05:52:56 pm
# @Filename: test_command.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2021-03-13 06:43:06 pm
# @Copyright (c) 2021 MinghuiGroup, Soochow University
from click.testing import CliRunner
from pdb_profiling.commands.command import Interface

def test_command():
    runner = CliRunner()
    for task in ('insert-mutation --input test/pytest/data/mutation.tsv --usecols Alt,Pos,Ref,ftId',
                 'id-mapping',
                 'check-muta-conflict',
                 'sifts-mapping --chunksize 15',
                 'insert-sele-mapping --input pipe_select_mo.tsv',
                 'residue-mapping --input pipe_select_mo.tsv',
                 'export-residue-mapping -o e_resmap.tsv --sele',
                 'insert-sele-mutation-mapping -i e_resmap.tsv',
                 'sifts-mapping --func pipe_select_smr_mo --chunksize 10',
                 'insert-smr-mapping -i pipe_select_smr_mo.tsv',
                 'export-smr-residue-mapping -o e_smr_resmap.tsv --sele'
                 'insert-sifts-meta --input test/pytest/data/pdb_demo.tsv --api_suffix api/mappings/pfam/',
                 'insert-isoform-range',
                 ['sifts-mapping', '--func', 'unp_is_canonical_with_id', '--input', 'test/pytest/data/unp_demo.tsv', '--iteroutput', '--skip_pdbs', ''],
                 ['fetch1pdb', '-i', '4hho', '-a', 'atoms', '-t', 'B_ONLY',
                  '-p', 'label_asym_id=B', '-p', 'copy_all_categories=false'],
                 ['fetch1pdb', '-i', '3hl2', '-a', 'atoms',
                 '-p', 'encoding=bcif',
                 '-d', 'label_asym_id=A,label_seq_id=90',
                 '-d', 'label_asym_id=B,label_seq_id=100',
                 '-m', 'post',
                 '-t', 'A_90_B_10'],
                 ):
        result = runner.invoke(Interface, task.split(' ') if not isinstance(task, list) else task)
        assert result.exit_code == 0, str(task)

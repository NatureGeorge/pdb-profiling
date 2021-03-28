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
                 'sifts-mapping --func pipe_base --chunksize 15',
                 'residue-mapping --input pipe_base.tsv',
                 'insert-sifts-meta --input test/pytest/data/pdb_demo.tsv --api_suffix api/mappings/pfam/',
                 'insert-isoform-range',
                 ['sifts-mapping', '--func', 'unp_is_canonical_with_id', '--input', 'test/pytest/data/unp_demo.tsv', '--iteroutput', '--skip_pdbs', ''],
                 #['pisa-range', '--func', 'pipe_protein_ligand_interface', '--input', 'test/pytest/data/pdb_demo.tsv'],
                 'insert-sele-mapping --input test/pytest/data/pipe_select_mo.tsv --tag',
                 ):
        result = runner.invoke(Interface, task.split(' ') if not isinstance(task, list) else task)
        assert result.exit_code == 0, str(task)

from os import system

def test_command():
    system('pdb_profiling insert-mutation --input ./data/mutation.tsv --usecols Alt,Pos,Ref,ftId')
    system('pdb_profiling id-mapping')
    system('pdb_profiling sifts-mapping --func pipe_base --chunksize 15')
    system('pdb_profiling residue-mapping --input pipe_base.tsv')
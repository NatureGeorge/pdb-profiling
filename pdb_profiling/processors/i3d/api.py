# @Created Date: 2020-06-07 06:39:15 pm
# @Filename: api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-06-07 06:39:23 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pdb_profiling.utils import sort_sub_cols
from pathlib import Path
from pandas import read_csv, concat
from unsync import unsync
from pdb_profiling.log import Abclog
from pdb_profiling.fetcher.webfetch import UnsyncFetch
from pdb_profiling.utils import init_semaphore, init_folder_from_suffixes
from typing import Union
from pathlib import Path

class I3D_API(object):

    @staticmethod
    def concat_dat(files, outpath):
        outpath = Path(outpath)
        with outpath.open('w+') as out:
            for index, file in enumerate(files):
                with Path(file).open('r+') as cur:
                    if index:
                        next(cur)
                    for line in cur:
                        out.write(line)


    @staticmethod
    def download(folder: Path, orgnisim):
        pass
        

    def __init__(self, path, identity_cutoff: int = 90):
        self.identity_cutoff = identity_cutoff
        self.all = read_csv(path, sep="\t")
        
        self.all_unp_pairs = self.all[['PROT1', 'PROT2']].drop_duplicates().reset_index(drop=True)
        self.all_unp_pairs['from_i3d_unp_pairs'] = True
        self.all_unp_pairs_1 = self.all_unp_pairs
        self.all_unp_pairs_2 = self.all_unp_pairs.rename(columns={"PROT1":"PROT2", "PROT2":"PROT1"})
        self.all_unp_pairs_12 = concat([self.all_unp_pairs_1, self.all_unp_pairs_2]).drop_duplicates()

        self.chain_pairs = self.all[
            self.all.TYPE.eq('Structure') &
            self.all.SEQ_IDENT1.ge(identity_cutoff) &
            self.all.SEQ_IDENT2.ge(identity_cutoff) &
            (self.all.MODEL1 == self.all.MODEL2)
        ][['PDB_ID', 'CHAIN1', 'CHAIN2']].drop_duplicates()
        self.chain_pairs = self.pipe_i3d(self.chain_pairs)
        
        self.unp_chain_pairs = self.all[
            self.all.TYPE.eq('Structure') &
            self.all.SEQ_IDENT1.ge(identity_cutoff) &
            self.all.SEQ_IDENT2.ge(identity_cutoff) &
            (self.all.MODEL1 == self.all.MODEL2)
        ][['PROT1', 'PROT2', 'PDB_ID', 'CHAIN1', 'CHAIN2']].drop_duplicates()
        self.unp_chain_pairs['partner_1'] = self.unp_chain_pairs.CHAIN1+'_'+self.unp_chain_pairs.PROT1
        self.unp_chain_pairs['partner_2'] = self.unp_chain_pairs.CHAIN2+'_'+self.unp_chain_pairs.PROT2
        self.unp_chain_pairs = self.unp_chain_pairs[['PDB_ID', 'partner_1', 'partner_2']]
        self.unp_chain_pairs = self.pipe_i3d(self.unp_chain_pairs, ['partner_1', 'partner_2'], "from_i3d_unp_chain_pair")
        self.unp_chain_pairs[['chain_id_1', 'PROT1']] = self.unp_chain_pairs.partner_1.str.split('_', expand=True)
        self.unp_chain_pairs[['chain_id_2', 'PROT2']] = self.unp_chain_pairs.partner_2.str.split('_', expand=True)
        self.unp_chain_pairs.drop(columns=['partner_1', 'partner_2'], inplace=True)

    
    @staticmethod
    def pipe_i3d(df, focus_cols=['chain_id_1', 'chain_id_2'], newcol="from_i3d_chain_pair"):
        df = df.rename(columns={
                        "PDB_ID": "pdb_id", 
                        "CHAIN1": "chain_id_1",
                        "CHAIN2": "chain_id_2"}).reset_index(drop=True)
        df[newcol] = True
        return sort_sub_cols(df, focus_cols)
    
    def chain_pairs_p(self, pdbs):
        return self.chain_pairs[self.chain_pairs.pdb_id.isin(pdbs)]
    
    def unp_chain_pairs_p(self, pdbs):
        return self.unp_chain_pairs[self.unp_chain_pairs.pdb_id.isin(pdbs)]


class Interactome3D(Abclog):
    '''
    Implement Interactome3D API
    
    DEMO URL: https://interactome3d.irbbarcelona.org/user_data/bsubtilis/download/complete/interactions.dat
    '''
    root = 'https://interactome3d.irbbarcelona.org/'
    headers = {'accept': 'text/plain'}
    organisms = frozenset(
        {'human', 'athaliana', 'bsubtilis', 'btaurus', 'celegans', 
         'cjejuni', 'fly', 'ecoli', 'hpylori', 'mouse',
         'mtuberculosis', 'mpneumoniae', 'pfalciparum', 
         'rat', 'yeast', 'spombe', 'ssp', 'tpallidum'})
    api_set = frozenset({'api/getProteinStructures', 'api/getInteractionStructures', 'api/getUniprotACs',
                         'api/getPdbFile', 'api/getVersion'})

    @classmethod
    def set_folder(cls, folder: Union[Path, str]):
        folder = Path(folder)
        assert folder.exists(), "Folder not exist! Please create it or input a valid folder!"
        cls.folder, cls.split_folder = tuple(init_folder_from_suffixes(folder, ('interactome3d', 'interactome3d/split')))

    @classmethod
    @unsync
    async def set_web_semaphore(cls, web_semaphore_value):
        cls.web_semaphore = await init_semaphore(web_semaphore_value)

    @classmethod
    def get_web_semaphore(cls):
        return cls.web_semaphore

    @classmethod
    def task_unit(cls, folder, suffix, file_suffix, file_name, params=None):
        args = dict(
            url=f'{cls.root}{suffix}',
            headers=cls.headers)
        if params is not None:
            args['params'] = params
        return 'get', args, folder/f'{file_name}{file_suffix}'

    @classmethod
    def retrieve_metadata(cls, organism: str, dataset: str = 'complete', suffix: str = 'user_data/{organism}/download/{dataset}/'):
        assert organism in cls.organisms, f"Invalid organism!\nValid set:{cls.organisms}"
        return UnsyncFetch.single_task(cls.task_unit(cls.folder,
            f'{suffix.format(organism=organism, dataset=dataset)}interactions.dat', '.tsv', f'{organism}_{dataset}_interactions'), 
            cls.get_web_semaphore()).result()

    @classmethod
    def retrieve_all_meta(cls, dataset: str = 'complete', suffix: str = 'user_data/{organism}/download/{dataset}/'):
        tasks = [cls.task_unit(f'{suffix.format(organism=organism, dataset=dataset)}interactions.dat', '.tsv', f'{organism}_{dataset}_interactions') for organism in cls.organisms]
        return UnsyncFetch.multi_tasks(tasks, semaphore=cls.get_web_semaphore())

    @classmethod
    @unsync
    async def fetch(cls, suffix: str, file_suffix: str = '.xml', **kwargs):
        assert suffix in cls.api_set, f"Invalid api!\nValid set:{cls.api_set}"
        return await UnsyncFetch.single_task(cls.task_unit(cls.split_folder,
            suffix, file_suffix, '&'.join(f'{key}={value}' for key, value in kwargs.items()), params=kwargs),
            cls.get_web_semaphore())

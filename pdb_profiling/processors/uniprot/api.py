# @Created Date: 2019-12-08 06:46:49 pm
# @Filename: api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-16 10:54:32 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Iterable, Iterator, Optional, Union, Generator, Dict, List
import re
from time import perf_counter
import pandas as pd
import numpy as np
import logging
from pathlib import Path
import asyncio
from unsync import unsync, Unfuture
from copy import deepcopy
from collections import Counter
from pdb_profiling.log import Abclog
from pdb_profiling.fetcher.webfetch import UnsyncFetch
from pdb_profiling.processors.uniprot.process import ExtractIsoAlt
from pdb_profiling.utils import init_semaphore, init_folder_from_suffix


QUERY_COLUMNS: List[str] = [
    'id', 'length', 'reviewed',
    'comment(ALTERNATIVE%20PRODUCTS)',
    'feature(ALTERNATIVE%20SEQUENCE)',
    'genes', 'organism', 'protein%20names']


RESULT_COLUMNS: List[str] = [
    'Entry', 'Length', 'Status',
    'Alternative products (isoforms)',
    'Alternative sequence',
    'Gene names', 'Organism', 'Protein names']


COLUMNS_DICT: Dict = dict(zip(QUERY_COLUMNS, RESULT_COLUMNS))


RESULT_NEW_COLUMN: List[str] = ['yourlist', 'isomap']


BASE_URL: str = 'https://www.uniprot.org'


PARAMS: Dict = {
    # 'fil': 'organism%3A"Homo+sapiens+(Human)+[9606]"+AND+reviewed%3Ayes',
    # reviewed:yes+AND+organism:9606
    'columns': None,
    'query': None,
    'from': None,
    'to': 'ACC',
    'format': 'tab'}


class MapUniProtID(Abclog):
    '''
    Implement UniProt Retrieve/ID Mapping API
    '''

    def __init__(self, id_col: str, id_type: str,
                 dfrm: Optional[pd.DataFrame],
                 ids: Optional[Iterable] = None,
                 sites: Optional[Iterable] = None,
                 genes: Optional[Iterable] = None,
                 usecols: Optional[Iterable] = QUERY_COLUMNS,
                 site_col: Optional[str] = None,
                 gene_col: Optional[str] = None,
                 logger: Optional[logging.Logger] = None,
                 loggingPath: Optional[str] = None):

        self.init_logger(self.__class__.__name__, logger)
        if dfrm is not None:
            self.dfrm = dfrm.drop_duplicates().reset_index(drop=True)
        else:
            '''
            the length of dataframe is based on:
            
                * the num of `ids` if there is more than one id
                * the num of `sites` if there is just one id with specified `sites`
            '''
            if isinstance(ids, str):
                if sites is not None and not isinstance(sites, str):
                    index_len = len(sites)
                else:
                    index_len = 1
            else:
                index_len = len(ids)

            self.dfrm = pd.DataFrame(dict(zip(
                (col for col in (id_col, site_col, gene_col) if col is not None),
                (value for value in (ids, sites, genes) if value is not None))),
                index=list(range(index_len)))

        self.index = dfrm.index
        self.id_col = id_col
        self.id_type = id_type
        self.site_col = site_col
        self.gene_col = gene_col
        self.loggingPath = loggingPath
        if isinstance(usecols, str):
            PARAMS['columns'] = usecols
            usecols = usecols.split(',')
        elif isinstance(usecols, (Iterable, Iterator)):
            PARAMS['columns'] = ','.join(usecols)
        else:
            raise ValueError('Invalid usecols')
        self.usecols = usecols
        PARAMS['from'] = id_type
        if isinstance(loggingPath, (str, Path)):
            self.set_logging_fileHandler(loggingPath)

    @property
    def sites(self) -> Generator:
        if self.site_col is not None:
            for name, group in self.dfrm.groupby(by=self.id_col, sort=False):
                yield name, group[self.site_col]
        else:
            yield None

    @staticmethod
    def split_df(dfrm, colName, sep):
        """Split DataFrame"""
        df = dfrm.copy()
        return df.drop([colName], axis=1).join(df[colName].str.split(sep, expand=True).stack().reset_index(level=1, drop=True).rename(colName))

    def yieldTasks(self, lyst: Iterable, chunksize: int = 100, sep: str = ',') -> Generator:
        fileName = self.outputPath.stem
        for i in range(0, len(lyst), chunksize):
            cur_fileName = f'{fileName}+{i}'
            cur_params = deepcopy(PARAMS)
            cur_params['query'] = sep.join(lyst[i:i+chunksize]) # self.outputPath.suffix
            yield ('get', {'url': f'{BASE_URL}/uploadlists/', 'params': cur_params}, str(Path(self.outputPath.parent, cur_fileName+'.tsv')))

    def retrieve(self, outputPath: Union[str, Path], 
                       finishedPath: Optional[str] = None, 
                       sep: str = '\t', 
                       chunksize: int = 100, 
                       concur_req: int = 20, 
                       rate: float = 1.5,
                       ret_res: bool = True,
                       semaphore = None):
        finish_id = list()
        self.outputPath = Path(outputPath)
        self.result_cols = [COLUMNS_DICT.get(
            i, i) for i in self.usecols] + RESULT_NEW_COLUMN
        if finishedPath is not None:
            try:
                target_col = RESULT_NEW_COLUMN[0]
                finish: pd.Series = pd.read_csv(
                    finishedPath,
                    sep=sep,
                    usecols=[target_col],
                    names=self.result_cols,
                    skiprows=1,
                    header=None)[target_col]
            except Exception as e:
                col_to_add = RESULT_NEW_COLUMN[1]
                self.logger.warning(
                    f"{e}\nSomething wrong with finished raw file, probably without '{col_to_add}' column.")
                finish_df = pd.read_csv(
                    finishedPath, sep=sep, names=self.result_cols[:-1], skiprows=1, header=None)
                finish_df[col_to_add] = np.nan
                finish_df.to_csv(finishedPath, sep=sep, index=False)
                finish: pd.Series = finish_df[target_col]

            for query_id in finish:
                if ',' in query_id:
                    finish_id.extend(query_id.split(','))
                else:
                    finish_id.append(query_id)

        query_id: pd.Series = self.dfrm[self.id_col]
        if finish_id:
            rest_id = list(set(query_id) - set(finish_id))
        else:
            rest_id = query_id.unique()

        self.logger.info(
            f"Have finished {len(finish_id)} ids, {len(rest_id)} ids left.")
        res = UnsyncFetch.multi_tasks(
            tasks=self.yieldTasks(rest_id, chunksize),
            to_do_func=self.process,
            concur_req=concur_req, 
            rate=rate, 
            logger=self.logger,
            ret_res=ret_res,
            semaphore=semaphore)
        return res

    def getCanonicalInfo(self, dfrm: pd.DataFrame):
        """
        Will Change the dfrm

        * Add new column (canonical_isoform)
        * Change the content of column (UniProt)
        """
        # Get info from Alt Product file
        if self.altProPath is None:
            dfrm['canonical_isoform'] = np.nan
            return dfrm
        else:
            usecols = ["IsoId", "Sequence", "Entry", "UniProt"]
            altPro_df = pd.read_csv(self.altProPath, sep="\t", usecols=usecols)
            altPro_df = altPro_df[altPro_df["Sequence"]
                                  == "Displayed"].reset_index(drop=True)
            altPro_df.rename(
                columns={"IsoId": "canonical_isoform"}, inplace=True)
            # Modify dfrm
            dfrm = pd.merge(
                dfrm, altPro_df[["canonical_isoform", "Entry"]], how="left")
            return dfrm

    def getGeneStatus(self, handled_df: pd.DataFrame, colName: str = 'GENE_status'):
        """
        Will Change the dfrm, add Gene Status

        * Add new column (GENE) # if id_col != gene_col
        * Add new column (GENE_status)

        **About GENE_status**

        * ``False`` : First element of Gene names is not correspond with refSeq's GENE (e.g)
        * others(corresponding GENE)

        """
        self.gene_status_col = colName
        if self.id_type != 'GENENAME':
            if self.gene_col is None:
                handled_df[colName] = True
                return None
            gene_map = self.dfrm[[self.id_col,
                                  self.gene_col]].drop_duplicates()
            gene_map = gene_map.groupby(self.id_col)[self.gene_col].apply(
                lambda x: np.array(x) if len(x) > 1 else list(x)[0])
            handled_df['GENE'] = handled_df.apply(
                lambda z: gene_map[z['yourlist']], axis=1)
            handled_df[colName] = handled_df.apply(lambda x: x['GENE'] == x['Gene names'].split(
                ' ')[0] if not isinstance(x['Gene names'], float) else False, axis=1)
            handled_df['GENE'] = handled_df['GENE'].apply(
                lambda x: ','.join(x) if not isinstance(x, str) else x)
        else:
            handled_df[colName] = handled_df.apply(lambda x: x['yourlist'] == x['Gene names'].split(
                ' ')[0] if not isinstance(x['Gene names'], float) else False, axis=1)

    def label_mapping_status(self, dfrm: pd.DataFrame, colName: str = 'Mapping_status'):
        self.mapping_status_col = colName
        gene_status_col = self.gene_status_col
        dfrm[colName] = 'No'
        dfrm[gene_status_col] = dfrm[gene_status_col].apply(
            lambda x: x.any() if isinstance(x, Iterable) else x)

        if self.id_col == 'GENENAME':
            pass_df = dfrm[
                (dfrm[gene_status_col] == True) &
                (dfrm['Status'] == 'reviewed') &
                (dfrm['unp_map_tage'] != 'Untrusted & No Isoform')]
        else:
            pass_df = dfrm[
                (dfrm['Status'] == 'reviewed') &
                (dfrm['unp_map_tage'] != 'Untrusted & No Isoform')]
        pass_index = pass_df.index
        dfrm.loc[pass_index, colName] = 'Yes'

        # Deal with 'one to many' situation
        multipleCounter = Counter(dfrm.loc[pass_index, 'yourlist'])
        err_li = [i for i, j in multipleCounter.items() if j > 1]
        err_index = pass_df[pass_df['yourlist'].isin(err_li)].index
        dfrm.loc[err_index, colName] = 'Error'

    @unsync
    async def process(self, path: Union[str, Path, Unfuture], sep: str = '\t'):
        self.logger.debug("Start to handle id mapping result")
        if not isinstance(path, (Path, str)):
            path = await path  # .result()
        if not Path(path).stat().st_size:
            return None
        self.altSeqPath, self.altProPath = ExtractIsoAlt.main(path=path)
        try:
            df = pd.read_csv(
                path, sep='\t', names=self.result_cols, skiprows=1, header=None)
        except ValueError:
            df = pd.read_csv(
                path, sep='\t', names=self.result_cols[:-1], skiprows=1, header=None)

        # Add New Column: canonical_isoform
        df = self.getCanonicalInfo(df)
        # Add New Column: unp_map_tage
        df['unp_map_tage'] = np.nan
        # Classification
        df_with_no_isomap = df[df['isomap'].isnull()]  # Class A
        df_with_isomap = df[df['isomap'].notnull()]  # Class B
        # ----------------------------------------------------------------------
        # In Class A
        # ----------------------------------------------------------------------
        if len(df_with_no_isomap) > 0:
            df_wni_split = self.split_df(df_with_no_isomap, 'yourlist', ',')
            df_wni_split.drop(columns=['isomap'], inplace=True)
            # [yourlist <-> UniProt]
            df_wni_split['UniProt'] = df_wni_split['Entry']
            df_wni_split['unp_map_tage'] = 'Trusted & No Isoform'
            # Find out special cases 1
            df_wni_split_warn = df_wni_split[df_wni_split['Alternative products (isoforms)'].notnull(
            )].index
            df_wni_split.loc[df_wni_split_warn,
                             'unp_map_tage'] = 'Untrusted & No Isoform'
            # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt'
        # ----------------------------------------------------------------------
        # In Class B
        # ----------------------------------------------------------------------
        if len(df_with_isomap) > 0:
            wi_yourlist_count = df_with_isomap.apply(
                lambda x: x['yourlist'].count(','), axis=1)
            wi_isomap_count = df_with_isomap.apply(
                lambda x: x['isomap'].count(','), axis=1)
            # In subClass 1
            df_wi_eq = df_with_isomap.loc[wi_yourlist_count[wi_yourlist_count ==
                                                            wi_isomap_count].index]
            if len(df_wi_eq) > 0:
                df_wi_eq_split = self.split_df(
                    df_wi_eq.drop(columns=['yourlist']), 'isomap', ',')
                df_wi_eq_split[['yourlist', 'UniProt']] = df_wi_eq_split['isomap'].str.split(
                    ' -> ', expand=True)
                # [yourlist <-> UniProt]
                df_wi_eq_split.drop(columns=['isomap'], inplace=True)
                df_wi_eq_split['unp_map_tage'] = 'Trusted & Isoform'
                # # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt'

            # In subClass 2
            df_wi_ne = df_with_isomap.loc[wi_yourlist_count[wi_yourlist_count !=
                                                            wi_isomap_count].index]
            if len(df_wi_ne) > 0:
                df_wi_ne_split = self.split_df(df_wi_ne, 'isomap', ',')
                df_wi_ne_split.rename(
                    columns={'yourlist': 'checkinglist'}, inplace=True)
                df_wi_ne_split[['yourlist', 'UniProt']] = df_wi_ne_split['isomap'].str.split(
                    ' -> ', expand=True)
                df_wi_ne_split.drop(columns=['isomap'], inplace=True)
                df_wi_ne_split['unp_map_tage'] = 'Trusted & Isoform & Contain Warnings'
                # 'Entry', 'Gene names', 'Status', 'Alternative products (isoforms)', 'Organism', 'yourlist', 'UniProt', 'checkinglist'
                # Find out special cases 2
                usecols = pd.Index(set(df_wi_ne_split.columns) -
                                   {'yourlist', 'UniProt'})
                df_wi_ne_warn = self.split_df(
                    df_wi_ne_split[usecols].drop_duplicates(), 'checkinglist', ',')
                df_wi_ne_warn = df_wi_ne_warn[~df_wi_ne_warn['checkinglist'].isin(
                    df_wi_ne_split['yourlist'])].rename(columns={'checkinglist': 'yourlist'})
                df_wi_ne_warn['UniProt'] = df_wi_ne_warn['Entry']
                # sequence conflict
                df_wi_ne_warn['unp_map_tage'] = 'Untrusted & No Isoform'
                df_wi_ne_split.drop(columns=['checkinglist'], inplace=True)

        # Concat Dfrm
        variables = ["df_wni_split", "df_wi_eq_split",
                     "df_wi_ne_split", "df_wi_ne_warn"]
        lvs = locals()
        varLyst = [lvs[variable] for variable in variables if variable in lvs]
        final_df = pd.concat(varLyst, sort=False).reset_index(drop=True)
        cano_index = final_df[final_df["canonical_isoform"].notnull()].index
        if len(cano_index) > 0:
            final_df.loc[cano_index, "UniProt"] = final_df.loc[cano_index, ].apply(
                lambda x: x["Entry"] if x["UniProt"] in x["canonical_isoform"] else x["UniProt"], axis=1)

        # Add Gene Status
        self.getGeneStatus(final_df)
        # Label Mapping Status
        self.label_mapping_status(final_df)

        pathOb = Path(path)
        edPath = str(Path(pathOb.parent, f'{pathOb.stem}_ed.tsv'))  # {pathOb.suffix}
        final_df.to_csv(edPath, sep=sep, index=False)
        self.logger.debug(f"Handled id mapping result saved in {edPath}")
        return edPath


class UniProtFASTA(Abclog):
    '''
    Download UniProt Fasta Sequences

    >>> UniProtFASTA.retrieve(
        ('Q6NZ36', 'P12755'),
        init_folder_from_suffix(yourfolder, 'UniProt/fasta/'))
    '''

    params = {'include': 'yes'}
    obj = {}

    @classmethod
    @unsync
    async def set_web_semaphore(cls, web_semaphore_value:int):
        cls.web_semaphore = await init_semaphore(web_semaphore_value)
    
    @classmethod
    def set_folder(cls, folder: Union[Path, str]):
        cls.folder = init_folder_from_suffix(folder, 'UniProt/fasta/')

    @classmethod
    @unsync
    async def process(cls, path: Union[str, Path, Unfuture]):
        if not isinstance(path, (Path, str)):
            path = await path  # .result()
        path = Path(path)
        if not path.stat().st_size:
            return None
        if not cls.obj:
            folder = path.parent
            kwargs = {'ret_res': False}
        else:
            folder = cls.obj['fasta_folder']
            kwargs = {'concur_req': cls.obj['unp_concurreq'], 'rate': cls.obj['unp_concurrate'],
                      'ret_res': False, 'semaphore': cls.obj['semaphore']}
        unps = pd.read_csv(path, sep='\t', usecols=['UniProt']).UniProt.drop_duplicates()
        for fob in cls.retrieve(unps, folder, **kwargs):
            await fob
        return path

    @classmethod
    def task_unit(cls, unp:str, folder: Union[str, Path]):
        cur_fileName = f'{unp}.fasta'
        cur_filePath = str(Path(folder, cur_fileName))
        return ('get', {'url': f'{BASE_URL}/uniprot/{cur_fileName}', 'params': cls.params}, cur_filePath)

    @classmethod
    def yieldTasks(cls, lyst: Iterable, folder: Union[str, Path]) -> Generator:
        for unp in lyst:
            return cls.task_unit(unp, folder)

    @classmethod
    def retrieve(cls, lyst: Iterable, folder: Union[str, Path], concur_req: int = 20, rate: float = 1.5, ret_res: bool = True, semaphore=None):
        return UnsyncFetch.multi_tasks(
            cls.yieldTasks(lyst, folder), 
            concur_req=concur_req, 
            rate=rate, 
            logger=cls.logger,
            ret_res=ret_res,
            semaphore=semaphore)
    
    @classmethod
    def single_retrieve(cls, identifier: str, folder: Union[str, Path], semaphore, rate: float = 1.5):
        return UnsyncFetch.single_task(
            task=cls.task_unit(identifier, folder),
            semaphore=semaphore,
            rate=rate)

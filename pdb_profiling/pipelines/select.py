# @Created Date: 2020-06-07 06:38:00 pm
# @Filename: select.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-06-07 06:38:04 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pandas import DataFrame, read_csv, isna, merge
import numpy as np
from unsync import unsync
import asyncio
import aiofiles
import ujson as json
import orjson
from pathlib import Path
import networkx as nx
from tablib import Dataset
from typing import Dict
from textdistance import overlap
from pdb_profiling.processers.pdbe.sqlite_api import converters
from pdb_profiling.processers.pdbe.neo4j_api import SIFTS, slice_series, lyst2range, related_dataframe
import logging

def str2ord(string: str):
    if len(string) > 1:
        num = 0
        for i in string:
            num += ord(i)
        return num
    else:
        return ord(string)


async def pipe_out(df, path):
    path = Path(path)
    if isinstance(df, DataFrame):
        sorted_col = sorted(df.columns)
        if path.exists():
            headers = None
        else:
            headers = sorted_col
        async with aiofiles.open(path, 'a') as fileOb:
            dataset = Dataset(headers=headers)
            dataset.extend(df[sorted_col].to_records(index=False))
            await fileOb.write(dataset.export('tsv'))
    elif isinstance(df, Dataset):
        async with aiofiles.open(path, 'a') as fileOb:
            await fileOb.write(df.export('tsv'))
    else:
        raise TypeError("Invalid Object for pipe_out()")


class Select_API(object):

    json_cols = ('entity_chain_map', 'unp_entity_info',
            'entity_unp_info', 'del_entity_chain')
    usecols = ['UniProt', 'pdb_id', 'chain_id',
               'new_unp_range', 'entity_id', 'RAW_BS']

    def __init__(self, sqlite, neo4j, folder, oscutoff, logger, **kwargs):
        self.sqlite = sqlite
        self.neo4j = neo4j
        self.oligo_path = folder/'oligo.tsv'
        self.moho_path = folder/'select_moho.tsv'
        self.he_path = folder/'select_he.tsv'
        self.oscutoff = oscutoff
        self.logger = logger
        self.kwargs = kwargs

    @staticmethod
    def select_range(ranges, indexes, cutoff=0.2, skip_index=[]):
        select_index = []
        def unit(cur_index):
            if cur_index in skip_index:
                return
            cur_range = orjson.loads(ranges[cur_index])
            for selected in select_index:
                selected_range = orjson.loads(ranges[selected])
                score = overlap.similarity(lyst2range(cur_range),
                                   lyst2range(selected_range))
                if score > cutoff:
                    return
            select_index.append(cur_index)

        for index in indexes:
            unit(index)
        return select_index

    @staticmethod
    def check_range(focus_sifts_nda, cutoff):
        skip_index = []
        # checked_id = []
        for checking_index, record in enumerate(focus_sifts_nda):
            pdb_id = record[1]
            # if pdb_id not in checked_id:
            #     checked_id.append(pdb_id)
            # else:
            #     continue
            related_records = focus_sifts_nda[focus_sifts_nda[:, 1] == pdb_id]
            if len(related_records) > 1:
                cur_range = orjson.loads(record[3])
                cur_chain = record[2]
                count = 0
                for other_record in related_records:
                    other_chain = other_record[2]
                    if other_chain == cur_chain:
                        continue
                    other_range = orjson.loads(other_record[3])
                    score = overlap.similarity(lyst2range(cur_range),
                                                lyst2range(other_range))
                    if score < 1 - cutoff:
                        continue
                    else:
                        count += 1
                if not count:
                    skip_index.append(checking_index)
            else:
                skip_index.append(checking_index)
        return skip_index

    def select_o(self, sifts_df, entry_info, cutoff=0.2, check_range:bool=False, tag:str='mo') -> Dataset:
        sifts_dict = slice_series(sifts_df.UniProt.to_numpy())
        sifts_nda = sifts_df.to_numpy()
        dataob = Dataset(headers=self.usecols+['1/resolution', 'revision_date_score', 'chain_id_score', 'select_tag', 'oligo_type'])
        for start, end in sifts_dict.values():
            focus_sifts_nda = sifts_nda[start:end]
            # TODO: check range
            if check_range:
                skip_index = self.check_range(focus_sifts_nda, cutoff)
            else:
                skip_index = []
            # TODO: Sort ndarray
            # Score, 1/R, Date, chain_id
            bs_chainid = focus_sifts_nda[:, (-1,2)]
            rs_date = [entry_info[pdb_id] for pdb_id in focus_sifts_nda[:, 1]]
            lyst = np.array([(bs, rs, str2ord(date), -str2ord(chain_id)) for (bs, chain_id), (rs, date) in zip(bs_chainid, rs_date)])
            sorted_index = np.lexsort(np.transpose(lyst)[::-1])[::-1]
            # TODO: select range
            selected_index = self.select_range(focus_sifts_nda[:, 3], sorted_index, cutoff, skip_index)
            new_sifts_nda = np.concatenate(
                (focus_sifts_nda[sorted_index, :-1], lyst[sorted_index]), axis=1)
            select_tag = np.array([(0, tag) for _ in range(len(sorted_index))])
            if len(selected_index):
                select_tag[selected_index] = np.array([(1, tag) for _ in range(len(selected_index))])
            new_sifts_nda = np.concatenate((new_sifts_nda, select_tag[sorted_index]), axis=1)
            for i in new_sifts_nda:
                dataob.append(list(i))
        return dataob

    @unsync
    async def select_e(self, sifts_df, entry_info, oligo_info, cutoff=0.2) -> Dataset:
        sifts_dict = slice_series(sifts_df.UniProt.to_numpy())
        sifts_nda = sifts_df.to_numpy()
        pass

    @staticmethod
    def yieldallchains(dictob):
        for lyst in dictob.values():
            yield from lyst
    
    @unsync
    async def process(self, sifts_path):
        sifts_df = read_csv(sifts_path, sep="\t", converters=converters, usecols=self.usecols)[self.usecols]
        sifts_df.drop_duplicates(inplace=True)
        sifts_df.sort_values(by='UniProt', inplace=True)
        pdbs = sifts_df.pdb_id.unique()
        oligo_df = await SIFTS.pipe_oligo(pdbs, self.neo4j.neo4j_api, self.sqlite, **self.kwargs)
        valid_chains = dict(zip(oligo_df.pdb_id, oligo_df.entity_chain_map.apply(
            lambda x: list(self.yieldallchains(x)))))
        check_chains = sifts_df[['pdb_id', 'chain_id']].apply(lambda x: x['chain_id'] in valid_chains[x['pdb_id']], axis=1)
        sifts_df.drop(index=check_chains[check_chains.eq(False)].index, inplace=True)
        pdbs = sifts_df.pdb_id.unique()
        entry_info = await self.sqlite.Entry_Info.objects.filter(pdb_id__in=pdbs).all()
        entry_info = {value.pdb_id: (1/value.resolution, value.REVISION_DATE) for value in entry_info}
        # TODO: select
        # dropna
        pdbs_focus = oligo_df.dropna(subset=['oligo_state']).pdb_id
        sifts_df = sifts_df[sifts_df.pdb_id.isin(pdbs_focus)]
        # mo
        pdbs_mo = oligo_df[oligo_df.oligo_state.eq('mo')].pdb_id
        sifts_mo_df = sifts_df[sifts_df.pdb_id.isin(pdbs_mo)]
        mo_res = self.select_o(sifts_mo_df, entry_info, self.oscutoff)
        # ho (+he)
        sifts_others_df = sifts_df.loc[sifts_df.index.difference(sifts_mo_df.index)]
        ho_res = self.select_o(sifts_others_df, entry_info, self.oscutoff, True, 'ho')
        # he
        pdbs_he = oligo_df[oligo_df.oligo_state.eq('he')].pdb_id.to_numpy()
        sifts_he_df = sifts_others_df[sifts_others_df.pdb_id.isin(pdbs_he)]
        ## TODO: Partner Pipe
        oligo_info = {key: value for key,value in zip(oligo_df.pdb_id, oligo_df.entity_unp_info) if key in pdbs_he}
        ### -----------------------------------------------------------------
        new_sifts_df = await self.neo4j.pipe_new_sifts(pdbs_he, 'pdb')
        ### -----------------------------------------------------------------
        # TODO: output
        for col in self.json_cols:
            oligo_df[col] = oligo_df[col].apply(lambda x: json.dumps(x) if not isna(x) else x)
        await pipe_out(oligo_df, self.oligo_path)
        await pipe_out(mo_res.stack(ho_res), self.moho_path)
        await pipe_out(new_sifts_df, self.he_path)

        


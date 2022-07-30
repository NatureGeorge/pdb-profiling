# @Created Date: 2020-11-13 03:02:36 pm
# @Filename: record.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2022-07-27 04:48:15 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from unsync import unsync, Unfuture
from typing import Union
from pathlib import Path
from pandas import DataFrame, concat, isna, Series
from numpy import nan
from re import compile as re_compile
from cachetools import LRUCache
from aiofiles import open as aiofiles_open
from pdb_profiling.utils import init_folder_from_suffix, a_concat, a_read_csv
from pdb_profiling.processors.uniprot import UniProtDB
from pdb_profiling.processors.uniprot.api import UniProtAPI, UniProtINFO

class UniProts(object):
    
    pattern_sep = re_compile(r"([A-z0-9]+):VAR_SEQ\s<?([0-9\.]+);\s+")
    pattern_iso_Value = re_compile(r'/([^=]+)="([^;]+)"')
    pattern_inIso = re_compile(r"([A-z\s\->]+)\s\(in\s([^\)]+)\)")
    pattern_iso_keyValue = re_compile(r"([^=]+)=([^;]+);\s+")
    pattern_inDel = re_compile(r"([A-z\s]+)\s->\s([A-z\s]+)")
    
    tasks = LRUCache(maxsize=1024)

    UniProtTXT = UniProtINFO('txt')

    @classmethod
    def register_task(cls, key, task):
        cls.tasks[key] = task

    @classmethod
    def set_folder(cls, folder: Union[Path, str]):
        cls.sqlite_api = UniProtDB("sqlite:///%s" % (init_folder_from_suffix(folder, 'local_db')/"uniprot.db"))
    
    @classmethod
    @unsync
    async def as_txt2tsv(cls, identifier:str, **kwargs):
        path = await cls.UniProtTXT.stream_retrieve_txt(identifier, **kwargs)
        if path is None:
            return identifier, None
        else:
            async with aiofiles_open(path, 'rt') as handle:
                txt = await handle.read()
            raw = txt.replace('FT                  ', '').replace('\n','').split('FT   VAR_SEQ         ')[1:]
            return identifier, '; '.join('VAR_SEQ '+i.replace(' /', ';  /') for i in raw)

    @staticmethod
    def getAltProInfo(inputLyst, groupCol="isoforms", flagCol="Name"):
        def toFrame(dict):
            df = DataFrame(dict[groupCol])
            for key, value in dict.items():
                if key == groupCol:
                    continue
                else:
                    df[key] = value
            return df

        outputDict = {groupCol: []}
        flag = 0
        for key, value in inputLyst:
            if not flag and key != flagCol:
                outputDict[key] = value
            if key == flagCol:
                outputDict[groupCol].append({key: value})
                flag += 1
            else:
                if flag:
                    outputDict[groupCol][flag-1][key] = value

        return toFrame(outputDict)

    @staticmethod
    def get_affected_interval(info, ac_len, ac_seq):
        if isinstance(info, float):
            return nan, nan, nan
        n_data = []
        delta_len = 0
        iso_seq = ''
        start_site = 0
        for index in range(len(info)):
            (cur_site, cur_old_end), (cur_old_len, cur_len), after_seq = info[index]
            delta_len += cur_len - cur_old_len
            iso_seq += ac_seq[start_site:cur_site-1] + after_seq
            start_site = cur_old_end
            if index == 0:
                n_data.append((cur_site, cur_site+cur_len-1 if cur_len > 0 else None))
            else:
                (b_site, b_old_end), (b_old_len, _), _ = info[index-1]
                b_start, b_end = n_data[index-1]
                cur_start = (b_end + cur_site - b_old_end) if b_end is not None else (b_start - b_site + cur_site - b_old_len)
                n_data.append((cur_start, cur_start+cur_len-1 if cur_len > 0 else None))
        n_data = [i for i in n_data if i[1] is not None]
        iso_seq += ac_seq[start_site:]
        return (n_data if n_data else nan), ac_len + delta_len, iso_seq

    '''
    @classmethod
    def get_alt_interval(cls, alt_id, alt_seq_dict):
        if alt_id in ("Displayed", "External", "Not described"):
            return nan, nan
        elif not alt_id.startswith("VSP"):
            raise ValueError("Unexpected alt_id: %s" % alt_id)
        else:
            alt_id_lyst = alt_id.split(", ")
            return cls.get_alt_interval_base(alt_id_lyst, alt_seq_dict)
    '''

    @staticmethod
    def get_alt_interval_base(alt_id_lyst, alt_seq_dict):
        info = []
        for altID in alt_id_lyst:
            try:
                index = alt_seq_dict["ftId"].index(altID)
            except ValueError:
                return []
            info.append((
                (alt_seq_dict["begin"][index], alt_seq_dict["end"][index]), 
                (alt_seq_dict["before_len"][index], alt_seq_dict["after_len"][index]),
                alt_seq_dict["after"][index]))
        assert len(info) > 0
        return sorted(info, key=lambda x: x[0][0])
    '''
    @classmethod
    def get_iso_len(cls, ac_len, sequenceStatus, info):
        if sequenceStatus == 'displayed':
            return ac_len
        elif sequenceStatus == 'described':
            _, alt_lens, _ = zip(*info)
            return ac_len + sum(after-before for before, after in alt_lens)
        else:
            return nan
    '''

    '''
    @classmethod
    def deal_with_alternative_pro(cls, dfrm, alt_seq_dict, ac_len, ac_seq):
        assert frozenset(dfrm.columns) >= frozenset(("Entry", "Alternative products (isoforms)", "Sequence", "Length"))
        altPro_df = concat((cls.getAltProInfo(cls.pattern_iso_keyValue.findall(i)) for i in dfrm["Alternative products (isoforms)"].dropna()), sort=False)
        altPro_df["AltInterval"] = altPro_df["Sequence"].apply(lambda x: cls.get_affected_interval(cls.get_alt_interval(x, alt_seq_dict)))
        return altPro_df
    '''

    @staticmethod
    def split_dot_range(x):
        res = x.split("..")
        if len(res) == 2:
            return [int(i) for i in res]
        elif len(res) == 1:
            val = int(res[0])
            return [val, val]
        else:
            raise AssertionError(f"Unexpected range: {x}")

    @classmethod
    @unsync
    def deal_with_alternative_seq(cls, dfrm):
        if isinstance(dfrm, Unfuture):
            dfrm = dfrm.result()
        assert frozenset(dfrm.columns) >= frozenset(("Entry", "Alternative sequence")), str(dfrm.columns)
        if all(dfrm['Alternative sequence'].isnull()):
            return
        dfrm['Alternative sequence'] = dfrm.apply(lambda x: x['Alternative sequence'].replace("VAR_SEQ", f"{x['Entry']}:VAR_SEQ") if not isna(x['Alternative sequence']) else nan, axis=1)
        altSeq_li = []
        for content in dfrm['Alternative sequence'].dropna().to_numpy():
            result = cls.pattern_sep.split(content)
            for i in range(1, len(result)-1, 3):
                altSeq_li.append(result[i+2] + '; /Entry="%s"; /AltRange="%s"; ' % tuple(result[i:i+2]))
        altSeq_df = DataFrame(dict(i.groups() for i in cls.pattern_iso_Value.finditer(content)) for content in altSeq_li)
        altSeq_df.rename(columns={"id": "ftId"}, inplace=True)
        try:
            altSeq_df[["AltInfo", "description"]] = altSeq_df.note.apply(lambda x: cls.pattern_inIso.search(x).groups()).apply(Series) 
        except AttributeError:
            from warnings import warn
            warn(str(dfrm))
            raise
        altSeq_df[['begin', 'end']] = altSeq_df.AltRange.apply(cls.split_dot_range).apply(Series)
        altSeq_df[["before", "after"]] = altSeq_df.AltInfo.apply(lambda x: cls.pattern_inDel.search(x).groups() if x != "Missing" else (nan, '')).apply(Series)
        altSeq_df.before = altSeq_df.before.apply(lambda x: x.replace(' ', '') if isinstance(x, str) else x)
        altSeq_df.after = altSeq_df.after.apply(lambda x: x.replace(' ', '') if isinstance(x, str) else x)
        altSeq_df['before_len'] = altSeq_df.apply(lambda x: len(x['before']) if isinstance(x['before'], str) else len(range(x['begin'], x['end']))+1, axis=1)
        altSeq_df['after_len'] = altSeq_df.after.apply(len)

        sup = DataFrame({
            'Entry': dfrm[dfrm['Alternative sequence'].isnull()].Entry.unique(),
            'ftId': '',
            'begin': nan,
            'end': nan,
            'before': nan,
            'after': nan,
            'description': nan,
            'before_len': nan,
            'after_len': nan,
            'evidence': nan})
        return concat((altSeq_df.drop(columns=['note', 'AltRange', 'AltInfo']), sup), ignore_index=True, sort=False)
    
    @classmethod
    @unsync
    async def fetch_VAR_SEQ_to_DB(cls, accessions, name='VAR_SEQ', via_txt:bool=True, **kwargs):
        if via_txt:
            tasks = [cls.as_txt2tsv(ac, name_suffix=name) for ac in accessions]
            dfrm = await cls.deal_with_alternative_seq(DataFrame([await task for task in tasks], columns=['Entry', 'Alternative sequence']))
            if dfrm is not None:
                await cls.sqlite_api.async_insert(cls.sqlite_api.VAR_SEQ, dfrm.to_dict('records'))
        else:
            UniProtAPI.params['columns'] = kwargs.get('columns', 'id,feature(ALTERNATIVE%20SEQUENCE)')
            UniProtAPI.params['from'] = kwargs.get('from', 'ACC+ID')
            UniProtAPI.with_name_suffix = kwargs.get('with_name_suffix', True)
            for res in UniProtAPI.retrieve(accessions, name, **kwargs):
                dfrm = await a_read_csv(res, sep='\t').then(cls.deal_with_alternative_seq)
                if dfrm is None:
                    continue
                await cls.sqlite_api.async_insert(cls.sqlite_api.VAR_SEQ, dfrm.to_dict('records'))
    
    @classmethod
    @unsync
    async def fetch_VAR_SEQ_from_DB(cls, accessions, name='VAR_SEQ', **kwargs):
        task = cls.tasks.get(('fetch_VAR_SEQ_from_DB', tuple(accessions)), None)
        if task is not None:
            return task
        res = await cls.sqlite_api.VAR_SEQ.objects.filter(Entry__in=accessions).all()
        if len(res) == 0:
            await cls.fetch_VAR_SEQ_to_DB(accessions, name, **kwargs)
            dfrm = DataFrame(await cls.sqlite_api.VAR_SEQ.objects.filter(Entry__in=accessions).all())
        else:
            dfrm = DataFrame(res)
            rest = tuple(frozenset(accessions) - frozenset(dfrm.Entry))
            if len(rest) > 0:
                await cls.fetch_VAR_SEQ_to_DB(rest, name, **kwargs)
                rest_df = DataFrame(await cls.sqlite_api.VAR_SEQ.objects.filter(Entry__in=rest).all())
                dfrm = concat((dfrm, rest_df), ignore_index=True, sort=False)
        cls.register_task(('fetch_VAR_SEQ_from_DB', tuple(accessions)), dfrm)
        return dfrm
        
        
    

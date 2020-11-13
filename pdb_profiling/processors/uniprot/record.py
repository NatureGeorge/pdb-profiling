# @Created Date: 2020-11-13 03:02:36 pm
# @Filename: record.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-11-13 03:02:39 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from unsync import unsync, Unfuture
from typing import Union
from pathlib import Path
from pandas import DataFrame, concat, isna, Series
from numpy import nan
from re import compile as re_compile
from pdb_profiling.utils import init_folder_from_suffix, a_concat, a_read_csv
from pdb_profiling.processors.uniprot import UniProtDB
from pdb_profiling.processors.uniprot.api import UniProtAPI

class UniProts(object):
    
    pattern_sep = re_compile(r"([A-z0-9]+):VAR_SEQ\s([0-9\.]+);\s+")
    pattern_iso_Value = re_compile(r'/([^=]+)="([^;]+)"')
    pattern_inIso = re_compile(r"([A-z\s\->]+)\s\(in\s([^\)]+)\)")
    pattern_iso_keyValue = re_compile(r"([^=]+)=([^;]+);\s+")
    pattern_inDel = re_compile(r"([A-z]+)\s->\s([A-z]+)")
    
    @classmethod
    def set_folder(cls, folder: Union[Path, str]):
        cls.sqlite_api = UniProtDB("sqlite:///%s" % (init_folder_from_suffix(folder, 'local_db')/"uniprot.db"))
    
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
    def getAffectedInterval(mis_info, inDe_info):
        if isinstance(mis_info, float):
            return nan
        affect_Interval = []
        for inDeRange, inDeLen in inDe_info:
            start = inDeRange[0]
            affect_Interval.append((start, start+inDeLen[1]-1))

        for index in range(len(affect_Interval)):
            raw = affect_Interval[index]
            interval = list(raw)
            for misRange, misLen in mis_info:
                if misRange[-1] < raw[0]:
                    interval = [x-misLen for x in interval]
            affect_Interval[index] = interval

        if affect_Interval:
            return affect_Interval
        else:
            return nan

    @classmethod
    def getAltInterval(cls, alt_id, altSeq_dict):
        if alt_id in ("Displayed", "External", "Not described"):
            return nan, nan
        elif not alt_id.startswith("VSP"):
            raise ValueError("Unexpected alt_id: %s" % alt_id)
        else:
            alt_id_lyst = alt_id.split(", ")
            return cls.getAltInterval_base(alt_id_lyst, altSeq_dict)

    @staticmethod
    def getAltInterval_base(alt_id_lyst, altSeq_dict):
        mis_info = []
        inDe_info = []
        for altID in alt_id_lyst:
            index = altSeq_dict["AltID"].index(altID)
            len_info = altSeq_dict["AltLen"][index]
            range_info = altSeq_dict["AltRange"][index]
            if isinstance(len_info, int):
                mis_info.append((range_info, len_info))
            else:
                inDe_info.append((range_info, len_info))
                if len_info[0] > len_info[1]:
                    mis_info.append(([range_info[0], range_info[1]-1], len_info[0]-len_info[1]))
        return mis_info, inDe_info

    @classmethod
    def deal_with_alternative_pro(cls, dfrm, altSeq_dict):
        assert frozenset(dfrm.columns) >= frozenset(("Entry", "Alternative products (isoforms)"))
        altPro_df = concat((cls.getAltProInfo(cls.pattern_iso_keyValue.findall(i)) for i in dfrm["Alternative products (isoforms)"].dropna()), sort=False)
        altPro_df["AltInterval"] = altPro_df["Sequence"].apply(lambda x: cls.getAffectedInterval(*cls.getAltInterval(x, altSeq_dict)))
        return altPro_df

    @classmethod
    @unsync
    def deal_with_alternative_seq(cls, dfrm):
        if isinstance(dfrm, Unfuture):
            dfrm = dfrm.result()
        assert frozenset(dfrm.columns) >= frozenset(("Entry", "Alternative sequence"))
        dfrm['Alternative sequence'] = dfrm.apply(lambda x: x['Alternative sequence'].replace("VAR_SEQ", f"{x['Entry']}:VAR_SEQ") if not isna(x['Alternative sequence']) else nan, axis=1)
        altSeq_li = []
        for content in dfrm['Alternative sequence'].dropna().to_numpy():
            result = cls.pattern_sep.split(content)
            for i in range(1, len(result)-1, 3):
                altSeq_li.append(result[i+2] + '; /Entry="%s"; /AltRange="%s"; ' % tuple(result[i:i+2]))
        altSeq_df = DataFrame(dict(i.groups() for i in cls.pattern_iso_Value.finditer(content)) for content in altSeq_li)
        altSeq_df.rename(columns={"id": "AltID"}, inplace=True)
        altSeq_df[["AltInfo", "AltIso"]] = altSeq_df.note.apply(lambda x: cls.pattern_inIso.search(x).groups()).apply(Series)
        altSeq_df["AltRange"] = altSeq_df.AltRange.apply(lambda x: [int(i) for i in x.split("..")])
        altSeq_df["AltLen"] = altSeq_df.apply(lambda x: [len(i) for i in cls.pattern_inDel.search(x["AltInfo"]).groups()] if x["AltInfo"] != "Missing" else len(range(*x["AltRange"]))+1, axis=1)
        # altSeq_dict = altSeq_df[["AltID", "AltLen", "AltRange"]].to_dict("list")
        sup = DataFrame({
            'Entry': dfrm[dfrm['Alternative sequence'].isnull()].Entry.unique(),
            'AltID': '',
            'AltRange': nan,
            'AltInfo': nan,
            'AltIso': nan,
            'AltLen': nan,
            'evidence': nan})
        return concat((altSeq_df.drop(columns=['note']), sup), ignore_index=True, sort=False)
    
    @classmethod
    @unsync
    async def fetch_VAR_SEQ_to_localDB(cls, accessions, name='VAR_SEQ', **kwargs):
        UniProtAPI.params['columns'] = kwargs.get('columns', 'id,feature(ALTERNATIVE%20SEQUENCE)')
        UniProtAPI.params['from'] = kwargs.get('from', 'ACC+ID')
        UniProtAPI.with_name_suffix = kwargs.get('with_name_suffix', True)
        for res in UniProtAPI.retrieve(accessions, name, **kwargs):
            dfrm = await a_read_csv(res, sep='\t').then(cls.deal_with_alternative_seq)
            await cls.sqlite_api.async_insert(cls.sqlite_api.VAR_SEQ, dfrm.to_dict('records'))
    
    @classmethod
    @unsync
    async def fetch_VAR_SEQ_from_localDB(cls, accessions, name='VAR_SEQ', **kwargs):
        res = await cls.sqlite_api.VAR_SEQ.objects.filter(Entry__in=accessions).all()
        if len(res) == 0:
            await cls.fetch_VAR_SEQ_to_localDB(accessions, name, **kwargs)
            dfrm = DataFrame(await cls.sqlite_api.VAR_SEQ.objects.filter(Entry__in=accessions).all())
        else:
            dfrm = DataFrame(res)
            rest = tuple(frozenset(accessions) - frozenset(dfrm.Entry))
            if len(rest) > 0:
                await cls.fetch_VAR_SEQ_to_localDB(rest, name, **kwargs)
                rest_df = DataFrame(await cls.sqlite_api.VAR_SEQ.objects.filter(Entry__in=rest).all())
                dfrm = concat((dfrm, rest_df), ignore_index=True, sort=False)
        return dfrm
        
        
    

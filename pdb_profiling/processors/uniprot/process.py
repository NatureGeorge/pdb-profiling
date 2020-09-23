# @Created Date: 2019-12-08 06:46:49 pm
# @Filename: process.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-16 10:54:32 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import re
from pathlib import Path
import pandas as pd
from numpy import nan


class ExtractIsoAlt:
    # pattern_all = re.compile(r"([A-z0-9]+):VAR_SEQ\s([0-9\s]+)\s([A-z\s\->]+)\s\(in\s([^\)]+)\)[^=]+FTId=([A-z0-9_,\s]+)")
    pattern_sep = re.compile(r"([A-z0-9]+):VAR_SEQ\s([0-9\.]+);\s+")
    pattern_iso_Value = re.compile(r'/([^=]+)="([^;]+)"')
    pattern_inIso = re.compile(r"([A-z\s\->]+)\s\(in\s([^\)]+)\)")
    pattern_iso_keyValue = re.compile(r"([^=]+)=([^;]+);\s+")
    pattern_inDe = re.compile(r"([A-z]+)\s->\s([A-z]+)")
    usecols = ["Entry", "Alternative sequence",
               "Alternative products (isoforms)"]

    def __init__(self, path: str, usecols=usecols, sep: str = '\t'):
        pathOb = Path(path)
        prefix = pathOb.stem
        self.dfrm = pd.read_csv(path, sep=sep, usecols=usecols)
        self.dfrm.dropna(inplace=True)
        self.dfrm.drop_duplicates(inplace=True)
        self.dfrm.reset_index(drop=True, inplace=True)
        self.altSeqPath = str(Path(pathOb.parent, f"{prefix}_UNP_AltSeq.tsv"))
        self.altProPath = str(Path(pathOb.parent, f"{prefix}_UNP_AltPro.tsv"))

    def __len__(self):
        return len(self.dfrm)

    def extractAltSeq(self):
        target_col_1 = "Alternative sequence"
        self.dfrm[target_col_1] = self.dfrm.apply(lambda x: x[target_col_1].replace(
            "VAR_SEQ", "{}:VAR_SEQ".format(x["Entry"])), axis=1)
        target_col_2 = "Alternative products (isoforms)"
        self.dfrm[target_col_2] = self.dfrm.apply(
            lambda x: "Entry={}; {}".format(x["Entry"], x[target_col_2]), axis=1)

        altSeq_li = []
        for content in self.dfrm[target_col_1].dropna():
            result = self.pattern_sep.split(content)
            for i in range(1, len(result)-1, 3):
                altSeq_li.append(
                    result[i+2] + '; /Entry="%s"; /AltRange="%s"; ' % tuple(result[i:i+2]))

        altSeq_df = pd.DataFrame(dict(i.groups() for i in self.pattern_iso_Value.finditer(
            content)) for content in altSeq_li)
        altSeq_df.rename(columns={"id": "AltID"}, inplace=True)
        altSeq_df[["AltInfo", "AltIso"]] = altSeq_df.apply(
            lambda x: self.pattern_inIso.search(x["note"]).groups(), result_type='expand', axis=1)
        altSeq_df["AltRange"] = altSeq_df["AltRange"].apply(
            lambda x: [int(i) for i in x.split("..")])

        altSeq_df["AltLen"] = altSeq_df.apply(lambda x: [len(i) for i in self.pattern_inDe.search(
            x["AltInfo"]).groups()] if x["AltInfo"] != "Missing" else len(range(*x["AltRange"]))+1, axis=1)
        altSeq_df.to_csv(self.altSeqPath, sep="\t", index=False)
        self.altSeq_dict = altSeq_df[[
            "AltID", "AltLen", "AltRange"]].to_dict("list")

    @staticmethod
    def getAltProInfo(inputLyst, groupCol="isoforms", flagCol="Name"):
        def toFrame(dict):
            df = pd.DataFrame(dict[groupCol])
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
    def getAltInterval(alt_id, altSeq_dict):
        if alt_id in ["Displayed", "External", "Not described"]:
            return nan, nan
        elif not alt_id.startswith("VSP"):
            raise ValueError("Unexcepted alt_id: %s" % alt_id)
        else:
            mis_info = []
            inDe_info = []
            alt_id_lyst = alt_id.split(", ")
            for altID in alt_id_lyst:
                index = altSeq_dict["AltID"].index(altID)
                len_info = altSeq_dict["AltLen"][index]
                range_info = altSeq_dict["AltRange"][index]
                if isinstance(len_info, int):
                    mis_info.append((range_info, len_info))
                else:
                    inDe_info.append((range_info, len_info))
                    if len_info[0] > len_info[1]:
                        mis_info.append(
                            ([range_info[0], range_info[1]-1], len_info[0]-len_info[1]))

        return mis_info, inDe_info

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

    def extractAltPro(self):
        altPro_df = pd.concat((self.getAltProInfo(self.pattern_iso_keyValue.findall(
            i)) for i in self.dfrm["Alternative products (isoforms)"].dropna()), sort=False)
        altPro_df["AltInterval"] = altPro_df["Sequence"].apply(
            lambda x: self.getAffectedInterval(*self.getAltInterval(x, self.altSeq_dict)))
        altPro_df["UniProt"] = altPro_df.apply(lambda x: x["IsoId"].split(
            ",")[0] if x["Sequence"] != "Displayed" else x["Entry"], axis=1)
        altPro_df.to_csv(self.altProPath, sep="\t", index=False)

    @classmethod
    def main(cls, **kwargs):
        demo = cls(**kwargs)
        if len(demo.dfrm) > 0:
            demo.extractAltSeq()
            demo.extractAltPro()
            return demo.altSeqPath, demo.altProPath
        else:
            return None, None

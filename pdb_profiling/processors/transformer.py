# @Created Date: 2020-01-12 01:27:18 pm
# @Filename: transformer.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-11 04:22:22 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pyexcel import get_sheet, Sheet
import tablib
from collections import OrderedDict
from typing import Dict, Iterable, Tuple


class Dict2Tabular(object):
    """
    Convert Json-Dict to Tabular format
    """

    @staticmethod
    def sync_with_pyexcel(*args) -> Sheet:
        records, *remain = args
        if not len(records):
            return None
        keys = {key for record in records for key in record.keys()}
        head = records[0]
        records[0] = {key: head.get(key, None) for key in keys}
        sheet = get_sheet(records=records)
        if len(remain) > 1:
            append_header, append_value = remain
            append_value = tuple(i if i is not None else '' for i in append_value)
            append_data = [append_header] + [append_value for _ in range(len(records))]
            sheet.column += Sheet(append_data)
        sheet.name_columns_by_row(0)
        return sheet

    @classmethod
    def pyexcel_io(cls, data: Iterable[Tuple]) -> Sheet:
        cur_sheet = None
        for res in data:  # traverseSuffixes(suffix, data)
            try:
                cur_cols = set(cur_sheet.colnames)
                new_sheet = cls.sync_with_pyexcel(*res)
                new_cols = set(new_sheet.colnames)
                # assert cur_cols >= new_cols, f"Unexpected new columns: {new_cols-cur_cols}"
                if not cur_cols >= new_cols:
                    append_header = tuple(new_cols-cur_cols)
                    append_data = [append_header] + [
                        tuple('' for _ in range(len(append_header))) for _ in range(len(cur_sheet))]
                    # NOTE: cur_sheet and new_sheet have colnames while the new Sheet does not
                    cur_sheet.column += Sheet(append_data)
                cur_sheet = get_sheet(
                    records=list(cur_sheet.records)+list(new_sheet.records), 
                    name_columns_by_row=0)
            except AttributeError:
                cur_sheet = cls.sync_with_pyexcel(*res)
            # except TypeError:
            #     continue
        return cur_sheet

    @staticmethod
    def sync_with_tablib(*args) -> tablib.Dataset:
        records, *remain = args
        if not len(records):
            return None
        keys = {key for record in records for key in record.keys()}
        ob = tablib.Dataset()
        records = [OrderedDict(sorted((key, record.get(key, None)) for key in keys))
                   for record in records]
        ob.dict = records
        if len(remain) > 1:
            append_header, append_value = remain
            for i in range(len(append_value)):
                ob.append_col([append_value[i]]*len(records), append_header[i])
        return ob

    @classmethod
    def tablib_io(cls, data: Iterable[Tuple]) -> tablib.Dataset:
        cur_ob = None
        for res in data:
            try:
                cur_cols = set(cur_ob.headers)
                new = cls.sync_with_tablib(*res)
                new_cols = set(new.headers)
                assert cur_cols == new_cols, f"Unequal columns: new_cols-cur_cols={new_cols-cur_cols}, cur_cols-new_cols={cur_cols-new_cols}"
                cur_ob.dict += new.dict
            except AttributeError:
                cur_ob = cls.sync_with_tablib(*res)
        return cur_ob

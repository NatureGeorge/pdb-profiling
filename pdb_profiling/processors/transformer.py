# @Created Date: 2020-01-12 01:27:18 pm
# @Filename: transformer.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-11 04:22:22 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pyexcel import get_sheet, Sheet
from pandas import DataFrame, concat
from collections import OrderedDict
from typing import Dict, Iterable, Tuple
# from tablib import Dataset


class Dict2Tabular(object):
    """
    Convert Json-Dict to Tabular format
    """

    @staticmethod
    def sync_with_pandas(*args) -> DataFrame:
        records, *remain = args
        if not len(records):
            return None
        df = DataFrame(records)
        if len(remain) > 1:
            for append_header, append_value in zip(*remain):
                df[append_header] = append_value
        return df

    @classmethod
    def pandas_io(cls, data: Iterable[Tuple]):
        return concat((cls.sync_with_pandas(*res) for res in data), ignore_index=True, sort=False)

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
        new_sheets = []
        for res in data:  # traverseSuffixes(suffix, data)
            if cur_sheet is None:
                cur_sheet = cls.sync_with_pyexcel(*res)
                if cur_sheet is None:
                    continue
                cur_cols = set(cur_sheet.colnames)
            else:
                new_sheet = cls.sync_with_pyexcel(*res)
                if new_sheet is None:
                    continue
                new_cols = set(new_sheet.colnames)
                # assert cur_cols >= new_cols, f"Unexpected new columns: {new_cols-cur_cols}"
                if not cur_cols >= new_cols:
                    append_header = tuple(new_cols-cur_cols)
                    append_data = [append_header] + [
                        tuple('' for _ in range(len(append_header))) for _ in range(len(cur_sheet))]
                    # NOTE: cur_sheet and new_sheet have colnames while the new Sheet does not
                    cur_sheet.column += Sheet(append_data)
                    cur_cols |= new_cols
                new_sheets.append(new_sheet)
                for sheet in new_sheets:
                    e_cols = set(sheet.colnames)
                    x_cols = cur_cols-e_cols
                    if len(x_cols) > 0:
                        append_header = tuple(x_cols)
                        append_data = [append_header] + [
                            tuple('' for _ in range(len(append_header))) for _ in range(len(sheet))]
                        sheet.column += Sheet(append_data)
        yield cur_sheet
        yield from new_sheets

    '''
    @staticmethod
    def sync_with_tablib(*args) -> Dataset:
        records, *remain = args
        if not len(records):
            return None
        keys = {key for record in records for key in record.keys()}
        ob = Dataset()
        records = [OrderedDict(sorted((key, record.get(key, None)) for key in keys))
                   for record in records]
        ob.dict = records
        if len(remain) > 1:
            append_header, append_value = remain
            for i in range(len(append_value)):
                ob.append_col([append_value[i]]*len(records), append_header[i])
        return ob

    @classmethod
    def tablib_io(cls, data: Iterable[Tuple]) -> Dataset:
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
    '''

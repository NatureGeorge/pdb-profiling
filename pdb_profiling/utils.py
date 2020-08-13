# @Created Date: 2019-11-20 06:33:41 pm
# @Filename: utils.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-17 10:25:37 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import os
import gzip
import shutil
from typing import Optional, Union, Dict, Tuple, Iterable, Iterator, List
from logging import Logger
from pandas import read_csv, DataFrame, isna
import numpy as np
from pathlib import Path
import aiofiles
from tablib import Dataset
from pyexcel import Sheet
import asyncio
from unsync import unsync, Unfuture
from itertools import chain
import orjson as json
import logging
from io import StringIO


def to_interval(lyst: Union[Iterable, Iterator]) -> List:
    def pass_check(lyst):
        try:
            if not lyst or isna(lyst):
                return False
            else:
                return True
        except ValueError:
            if isinstance(lyst, float):
                return False
            else:
                return True
    if not pass_check(lyst):
        return []
    else:
        lyst = set(int(i) for i in lyst)
        if not pass_check(lyst):
            return []
        start = []
        interval_lyst = []
        true_interval_lyst = []
        max_edge = max(lyst)
        min_edge = min(lyst)
        if len(lyst) == (max_edge + 1 - min_edge):
            true_interval_lyst.append((min_edge, max_edge))
        else:
            lyst_list = sorted(lyst)
            for j in lyst_list:
                if not start:
                    i = j
                    start.append(j)
                    i += 1
                else:
                    if (i != j) or (j == max(lyst_list)):
                        if j == max(lyst_list):
                            if (i != j):
                                interval_lyst.append(start)
                                interval_lyst.append([j])
                                break
                            else:
                                start.append(j)
                        interval_lyst.append(start)
                        start = [j]
                        i = j + 1
                    else:
                        start.append(j)
                        i += 1
            for li in interval_lyst:
                max_edge = max(li)
                min_edge = min(li)
                true_interval_lyst.append((min_edge, max_edge))
        return true_interval_lyst


@unsync
async def init_semaphore(concurreq) -> Unfuture:
    """
    `semaphore` initiated in the `unsync` event loop
    """
    await asyncio.sleep(.01)
    return asyncio.Semaphore(concurreq)


def decompression(path: str, extension: str =".gz", remove: bool =True, outputPath: Optional[str] = None, logger: Optional[Logger] = None):
    """
    Decompress gz file

    :param path: The file path of the compressed file
    :param extension: Compressed file tail, default value: `.gz`
    :param remove: Whether remove the compressed file, default value: `True`
    :param outputPath: File safe path, default value: `None`
    :param logger: logger Object
    """

    """
    with gzip.GzipFile(mode="rb", fileobj=open(path, 'rb')) as raw:
        with open(path[:-len(extension)], "wb") as file:
            file.write(raw.read())
    """
    if outputPath is None:
        outputPath = path[:-len(extension)]

    with gzip.open(path, 'rb') as raw:
        with open(outputPath, 'wb') as file:
            shutil.copyfileobj(raw, file)
    try:
        if remove:
            os.remove(path)
    except Exception as e:
        if isinstance(logger, Logger):
            logger.error(e)

    return outputPath


def related_dataframe(filters: Optional[Union[Dict, Iterable[Tuple]]] = None, dfrm: Optional[DataFrame] = None, path: Union[str, Path, None] = None, sep: str = '\t', **kwargs):
    '''
    valid symbol: `eq, ne, le, lt, ge, gt, isin, isnull`
    '''
    
    if dfrm is None:
        if path is not None:
            dfrm = read_csv(path, sep=sep)
        else:
            raise ValueError('path should not be None')
    elif not isinstance(dfrm, DataFrame):
        raise ValueError('dfrm should be a pandas.DataFrame')

    if filters is None:
        return dfrm
    elif isinstance(filters, Dict):
        filters = filters.items()

    for col, (symbol, value) in filters:
        dfrm = dfrm[getattr(getattr(dfrm, col), symbol)(value)]
    return dfrm


def sort_sub_cols(dfrm, cols):
    if set(dfrm.columns) >= set(cols):
        dfrm[cols] = np.sort(dfrm[cols].to_numpy())
        return dfrm.drop_duplicates()
    else:
        return dfrm


async def a_read_csv(path, read_mode='r',**kwargs):
    '''
    only suitable for small dataframe
    '''
    async with aiofiles.open(path, read_mode) as file_io:
        with StringIO(await file_io.read()) as text_io:
            return read_csv(text_io, **kwargs)


async def pipe_out(df, path, **kwargs):
    if not isinstance(df, (DataFrame, Dataset, Sheet)):
        raise TypeError(f"Invalid Object for pipe_out(): {type(df)}")
    if len(df) == 0:
        raise ValueError("Zero record!")
    path = Path(path)
    var_format = kwargs.get('format', 'tsv')
    async with aiofiles.open(path, kwargs.get('mode', 'a')) as file_io:
        if isinstance(df, DataFrame):
            sorted_col = sorted(df.columns)
            if path.exists():
                headers = None
            else:
                headers = sorted_col
            dataset = Dataset(headers=headers)
            dataset.extend(df[sorted_col].to_records(index=False))
            to_write = dataset.export(var_format)
        elif isinstance(df, Dataset):
            to_write = df.export(var_format)
        elif isinstance(df, Sheet):
            to_write = getattr(df, var_format)
        else:
            pass
        await file_io.write(to_write)


def flatten_dict(data: Dict, root: str, with_root: bool = True):
    if with_root:
        iterator = yield_flatten_dict(data[root], root)
    else:
        iterator = yield_flatten_dict(data[root])
    for key, value in iterator:
        data[key] = value
    del data[root]


def yield_flatten_dict(data: Dict, root: Optional[str] = None):
    if root is None:
        yield from data.items()
    else:
        for key, value in data.items():
            yield f'{root}.{key}', value


def slice_series(se: Iterable) -> Dict:
    '''
    For Sorted Series
    '''
    data = {}
    cur = next(iter(se))
    start = 0
    try:
        for index, i in enumerate(se):
            if i != cur:
                assert cur not in data, "Invalid Series"
                data[cur] = (start, index)
                cur = i
                start = index
        assert cur not in data, "Invalid Series"
        data[cur] = (start, index+1)
    except AssertionError as e:
        logging.error(e)
        raise e
    return data


def split_df_by_chain(df, all_cols, cols_to_split, mode='sep', sep=','):
    '''
    Reference: <https://stackoverflow.com/a/50731258/12876491>
    '''
    def chainer_sep(s):
        return list(chain.from_iterable(s.str.split(sep)))

    def chainer_json(s):
        return list(chain.from_iterable(s.apply(json.loads)))
    
    if mode == 'sep':
        chainer = chainer_sep
        lens = df[cols_to_split[0]].str.split(sep).map(len)
    elif mode == 'json-list':
        chainer = chainer_json
        lens = df[cols_to_split[0]].apply(json.loads).map(len)
    else:
        raise ValueError("Invalid mode!")
    
    repeat_dict = {col: np.repeat(df[col], lens)
                   for col in set(all_cols)-set(cols_to_split)}
    chain_dict = {col: chainer(df[col]) for col in cols_to_split}
    return DataFrame({**repeat_dict, **chain_dict})


def init_folder_from_suffix(folder: Union[Path, str], suffixes: Iterable) -> Iterable[Path]:
    folder = Path(folder)
    for suffix in suffixes:
        new_path = folder/suffix
        new_path.mkdir(parents=True, exist_ok=True)
        yield new_path

# @Created Date: 2019-11-20 06:33:41 pm
# @Filename: utils.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-17 10:25:37 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import os
import gzip
import shutil
from typing import Optional, Union, Dict, Tuple, Iterable
from logging import Logger
from pandas import read_csv, DataFrame
import numpy as np
from pathlib import Path
import aiofiles
from tablib import Dataset
import asyncio
from unsync import unsync, Unfuture
from itertools import chain
import orjson as json


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


def split_df_by_chain(df, all_cols, cols_to_split, mode='sep'):
    '''
    Reference: <https://stackoverflow.com/a/50731258/12876491>
    '''
    def chainer_sep(s):
        return list(chain.from_iterable(s.str.split(',')))

    def chainer_json(s):
        return list(chain.from_iterable(s.apply(json.loads)))
    
    if mode == 'sep':
        chainer = chainer_sep
        lens = df[cols_to_split[0]].str.split(',').map(len)
    elif mode == 'json-list':
        chainer = chainer_json
        lens = df[cols_to_split[0]].apply(json.loads).map(len)
    else:
        raise ValueError("Invalid mode!")
    
    repeat_dict = {col: np.repeat(df[col], lens)
                   for col in set(all_cols)-set(cols_to_split)}
    chain_dict = {col: chainer(df[col]) for col in cols_to_split}
    return DataFrame({**repeat_dict, **chain_dict})



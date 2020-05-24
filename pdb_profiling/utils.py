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
from pathlib import Path


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

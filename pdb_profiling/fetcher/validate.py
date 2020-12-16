# @Created Date: 2020-12-14 11:43:18 pm
# @Filename: validate.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-12-14 11:43:20 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import orjson as json
from aiofiles import open as aiofiles_open
from typing import Coroutine
from unsync import Unfuture
from pathlib import Path
from pdb_profiling.exceptions import InvalidFileContentError
from re import compile as re_compile


class ValidateBase(object):

    fasta_pat = re_compile(r'(>.+)\n([A-Z\*\n]+)')

    @classmethod
    async def validate(cls, path):
        dispatch = {
            '.json': cls.json_load,
            '.cif': cls.cif_load,
            '.fasta': cls.fasta_load
        }
        if isinstance(path, (Coroutine, Unfuture)):
            path = await path
        path = Path(path)
        func = dispatch.get(path.suffix, None)
        if func is not None:
            try:
                await func(path)
            except (json.JSONDecodeError, AssertionError):
                raise InvalidFileContentError(path)
        else:
            return

    @staticmethod
    async def json_load(path):
        async with aiofiles_open(path, 'rt') as handle:
            data = await handle.read()
        json.loads(data)
    
    @staticmethod
    async def cif_load(path):
        async with aiofiles_open(path, 'rt') as handle:
            data = await handle.read()
        assert data.endswith('\n#\n')
    
    @classmethod
    async def fasta_load(cls, path):
        async with aiofiles_open(path, 'rt') as handle:
            data = await handle.read()
        assert bool(cls.fasta_pat.fullmatch(data))
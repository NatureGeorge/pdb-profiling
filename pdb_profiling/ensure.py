# @Created Date: 2020-12-20 07:56:06 pm
# @Filename: ensure.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-12-20 07:56:14 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Dict
from uuid import uuid4
import asyncio
from unsync import unsync
from tenacity import retry
import aiofiles
import aiofiles.os
from pathlib import Path
from warnings import warn
from pdb_profiling.exceptions import InvalidFileContentError
from pdb_profiling.warnings import ZeroSizeWarning, InvalidFileContentWarning, FileExistsWarning
from pdb_profiling.validate import ValidateBase


class EnsureBase(object):

    def __init__(self, use_existing: bool = True):
        self.use_existing = use_existing

    def set_use_existing(self, status:bool):
        self.use_existing = status

    @staticmethod
    def create_tmp_path_from_kw(kwargs: Dict):
        path = kwargs['path']
        temp_path = str(path)+f'.{uuid4().hex}.tmp'  # still possible to meet collision in high freq
        kwargs['path'] = temp_path
        return kwargs

    def make_sure_complete(self, *retry_args, **retry_kwargs):
        """
        Decorator to wrap a function | coroutinefunction | unsync object 
        to ensure the output file has been completed
        """
        def decorator(func):
            @retry(*retry_args, **retry_kwargs)
            async def wrapper(*args, **kwargs):
                raw_path = kwargs['path']
                if raw_path is None:
                    return
                exists, stat_result = await aio_file_exists_stat(raw_path)
                if self.use_existing and exists:
                    if not (stat_result.st_size > 0):
                        warn(str(raw_path), ZeroSizeWarning)
                    try:
                        await ValidateBase.validate(raw_path)
                        return raw_path
                    except InvalidFileContentError as e:
                        warn(f"InvalidFileContentError for '{raw_path}', will retry", InvalidFileContentWarning)
                        await aiofiles.os.remove(raw_path)
                        raise e
                if asyncio.iscoroutinefunction(func) or isinstance(func, unsync):
                    path = await func(*args, **self.create_tmp_path_from_kw(kwargs))
                else:
                    path = func(*args, **self.create_tmp_path_from_kw(kwargs))
                if path is None:
                    return
                try:
                    await ValidateBase.validate(path, suffix=Path(raw_path).suffix)
                except InvalidFileContentError as e:
                    warn(f"InvalidFileContentError for '{path}', will retry", InvalidFileContentWarning)
                    await aiofiles.os.remove(path)
                    raise e
                try:
                    await aiofiles.os.rename(path, raw_path)
                except FileExistsError:
                    if self.use_existing:
                        await aiofiles.os.remove(path)
                        warn(f"remove {str(path)}", FileExistsWarning)
                    else:
                        await aiofiles.os.remove(raw_path)
                        await aiofiles.os.rename(path, raw_path)
                return raw_path
            return wrapper
        return decorator


async def aio_file_exists_stat(path):
    try:
        stat_result = await aiofiles.os.stat(path)
    except (OSError, ValueError):
        return False, None
    return True, stat_result

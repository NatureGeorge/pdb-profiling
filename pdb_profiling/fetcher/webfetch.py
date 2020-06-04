# @Created Date: 2020-01-18 10:53:07 am
# @Filename: webfetch.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-09 08:50:23 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import os
from time import perf_counter
import asyncio
import aiohttp
import aioftp
import aiofiles
from unsync import unsync, Unfuture
from tenacity import retry, wait_random, stop_after_attempt, after_log, RetryError
import logging
from tqdm import tqdm
from typing import Iterable, Iterator, Union, Any, Optional, List, Dict, Coroutine, Callable
import ujson as json
from furl import furl
from pdb_profiling.log import Abclog
import re


class UnsyncFetch(Abclog):
    '''
    Fetch files through API in unsync/async way

    Note:

    * Implement both GET and POST method
    * Load all the available data(stream) in memory as soon as it is received
    * Since the methods in this class would not load the entire file in memory,
      the response data could be a large file
    * Parameters of `tenacity.retry` is built-in
    
    following packages provide me with a lot of insights and inspiration

    Reference 
    ---------

        `aiohttp`

        * https://github.com/aio-libs/aiohttp
        * https://docs.aiohttp.org/en/stable/client_quickstart.html
        * https://docs.aiohttp.org/en/stable/streams.html

        `unsync`

        * https://github.com/alex-sherman/unsync

        `tenacity`

        * https://github.com/jd/tenacity
        * https://tenacity.readthedocs.io/en/latest/
    '''

    retry_kwargs = {
        'wait': wait_random(max=10),
        'stop': stop_after_attempt(3),
        # 'after': after_log(logger, logging.WARNING)
        }
    use_existing: bool = False

    @classmethod
    async def http_download(cls, method: str, info: Dict, path: str):
        if cls.use_existing is True and os.path.exists(path):
            return path
        cls.logger.debug(f"Start to download file: {info['url']}")
        async with aiohttp.ClientSession() as session:
            async_func = getattr(session, method)
            async with async_func(**info) as resp:
                if resp.status == 200:
                    async with aiofiles.open(path, 'wb') as fileOb:
                        # Asynchronous iterator implementation of readany()
                        async for chunk in resp.content.iter_any():
                            await fileOb.write(chunk)
                    cls.logger.debug(f"File has been saved in: {path}")
                    return path
                elif resp.status in (404, 405):
                    cls.logger.warning(f"404/405 for: {info}")
                    return None
                else:
                    mes = "code={resp.status}, message={resp.reason}, headers={resp.headers}".format(resp=resp)
                    cls.logger.error(mes)
                    cls.logger.error(info)
                    raise Exception(mes)

    @classmethod
    async def ftp_download(cls, method: str, info: Dict, path: str):
        url = furl(info['url'])
        fileName = url.path.segments[-1]
        filePath = os.path.join(path, fileName)
        if cls.use_existing is True and os.path.exists(filePath):
            return filePath
        cls.logger.debug(f"Start to download file: {url}")  # info
        async with aioftp.ClientSession(url.host) as session:
            await session.change_directory('/'.join(url.path.segments[:-1]))
            await session.download(fileName, path) # , write_info=True
        cls.logger.debug(f"File has been saved in: {filePath}")
        return filePath

    @classmethod
    def download_func_dispatch(cls, method: str):
        method = method.lower()
        if method in ('get', 'post'):
            return cls.http_download
        elif method == 'ftp':
            return cls.ftp_download
        else:
            raise ValueError(f'Invalid method: {method}, valid method should be get, post or ftp')

    @classmethod
    @unsync
    async def save_file(cls, path: str, data: bytes):
        '''Deprecated'''
        async with aiofiles.open(path, 'wb') as fileOb:
            await fileOb.write(data)

    @classmethod
    @unsync
    async def fetch_file(cls, semaphore: asyncio.Semaphore, method: str, info: Dict, path: str, rate: float):
        download_func = cls.download_func_dispatch(method)
        try:
            async with semaphore:
                res = await download_func(method, info, path)
                if res is not None:
                    await asyncio.sleep(rate)
                return res
        except RetryError:
            cls.logger.error(f"Retry failed for: {info}")

    @classmethod
    @unsync
    async def unsync_tasks(cls, tasks):
        return [await fob for fob in tqdm(asyncio.as_completed(tasks), total=len(tasks))]

    @classmethod
    def multi_tasks(cls, tasks: Union[Iterable, Iterator],
              to_do_func: Optional[Callable] = None,
              concur_req: int = 4, rate: float = 1.5,
              logger: Optional[logging.Logger] = None,
              run_tasks: bool = True,
              semaphore = None
              ):

        cls.init_logger(cls.__name__, logger)
        cls.retry_kwargs['after'] = after_log(cls.logger, logging.WARNING)
        cls.http_download = retry(cls.http_download, **cls.retry_kwargs)
        cls.ftp_download = retry(cls.ftp_download, **cls.retry_kwargs)
        if semaphore is None:
            semaphore = asyncio.Semaphore(concur_req)
        else:
            cls.logger.info(f'{cls.multi_tasks.__qualname__}: pass asyncio.Semaphore')
        if to_do_func is None:
            tasks = [cls.fetch_file(semaphore, method, info, path, rate)
                     for method, info, path in tasks]
        else:
            tasks = [cls.fetch_file(semaphore, method, info, path, rate).then(
                to_do_func) for method, info, path in tasks]
        if run_tasks:
            t0 = perf_counter()
            res = cls.unsync_tasks(tasks).result()
            elapsed = perf_counter() - t0
            cls.logger.info('{} chunks downloaded in {:.2f}s'.format(len(res), elapsed))
            return res
        else:
            return tasks
        

    @classmethod
    def main(cls, workdir: str, data: Union[Iterable, Iterator], concur_req: int = 4, rate: float = 1.5, logName: str = 'UnsyncFetch'):
        cls.set_logging_fileHandler(os.path.join(workdir, f'{logName}.log'), logName=cls.__name__)
        t0 = perf_counter()
        # res = asyncio.run(cls.multi_tasks(data, concur_req=concur_req, rate=rate))
        res = cls.multi_tasks(data, concur_req=concur_req, rate=rate).result()
        elapsed = perf_counter() - t0
        cls.logger.info(f'downloaded in {elapsed}s')
        return res

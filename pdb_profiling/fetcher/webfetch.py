# @Created Date: 2020-01-18 10:53:07 am
# @Filename: webfetch.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-09 08:50:23 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from time import perf_counter
import asyncio
import aiohttp
from aiofiles import open as aiofiles_open
from unsync import unsync, Unfuture
from tenacity import retry, wait_random, stop_after_attempt, RetryError, retry_if_exception_type
from rich.progress import track
from typing import Iterable, Iterator, Union, Optional, Dict, Callable
from urllib.parse import urlparse
from aioftp import Client as aioftp_Client
from pathlib import Path
#from email.utils import parsedate_to_datetime
#from time import mktime
#from os import utime
from pdb_profiling.log import Abclog
from pdb_profiling.utils import init_semaphore
from pdb_profiling.exceptions import InvalidFileContentError, RemoteServerError
from pdb_profiling.ensure import EnsureBase


ensure = EnsureBase()
msc_rt_kw = dict(wait=wait_random(max=10), stop=stop_after_attempt(5), retry=retry_if_exception_type(InvalidFileContentError))
rt_kw = dict(wait=wait_random(max=20), stop=stop_after_attempt(6))


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

    @staticmethod
    @unsync
    @retry(**rt_kw)
    async def get_headers_only(url: str) -> Dict:
        async with aiohttp.ClientSession() as session:
            async with session.head(url) as resp:
                return resp.headers
        # resp.headers.get('Last-Modified', None), resp.headers.get('Etag', None)

    @classmethod
    @retry(**rt_kw)
    async def http_download(cls, semaphore, method: str, info: Dict, path: Union[str, Path], rate: float):
        '''
        possible exceptions:
            RemoteServerError, 
            aiohttp.client_exceptions.ServerDisconnectedError, 
            aiohttp.client_exceptions.ClientConnectorError, 
            aiohttp.client_exceptions.ClientPayloadError
        '''
        cls.logger.debug(f"http_download: Start to download file: {info['url']}")
        try:
            async with semaphore:
                async with aiohttp.ClientSession(connector=aiohttp.TCPConnector(ssl=False), trust_env=True) as session:
                    async_func = getattr(session, method)
                    async with async_func(**info) as resp:
                        if resp.status == 200:
                            async with aiofiles_open(path, 'wb') as file_out:
                                # Asynchronous iterator implementation of readany()
                                async for chunk in resp.content.iter_any():
                                    await file_out.write(chunk)
                            cls.logger.debug(f"File has been saved in: '{path}'")
                            await asyncio.sleep(rate)
                            #mod_time = resp.headers.get('Last-Modified', None)
                            #if mod_time is not None:
                            #    mod_time = mktime(parsedate_to_datetime(mod_time).timetuple())
                            #    utime(path, (mod_time, mod_time))
                            return path
                        elif resp.status in (204, 300, 400, 403, 404, 405, 406):
                            cls.logger.debug(f"{resp.status} for: {info}")
                            return None
                        else:
                            mes = "code={resp.status}, message={resp.reason}, headers={resp.headers}".format(resp=resp)
                            cls.logger.error(f"{info} -> {mes}")
                            raise RemoteServerError(mes)
        except Exception as e:
            cls.logger.error(f"{info} -> {e}")
            raise e

    @classmethod
    @retry(**rt_kw)
    async def ftp_download(cls, semaphore, method: str, info: Dict, path, rate: float):
        url = urlparse(info['url'])
        cls.logger.debug("ftp_download: Start to download file: {}".format(info['url']))
        async with semaphore:
            async with aioftp_Client.context(url.netloc) as client:
                if await client.exists(url.path):
                    async with aiofiles_open(path, 'wb') as file_out, client.download_stream(url.path) as stream:
                        async for block in stream.iter_by_block():
                            await file_out.write(block)
                    cls.logger.debug(f"File has been saved in: {path}")
                    await asyncio.sleep(rate)
                    return path
        return None

    @classmethod
    def download_func_dispatch(cls, method: str):
        method = method.lower()
        if method in ('get', 'post'):
            return cls.http_download
        elif method == 'ftp':
            return cls.ftp_download
        else:
            raise NotImplementedError(
                f'Invalid method: {method}, valid method should be get, post or ftp')

    @classmethod
    @unsync
    @ensure.make_sure_complete(**msc_rt_kw)
    async def fetch_file(cls, semaphore, method: str, info: Dict, path: Union[str, Path], rate: float):
        download_func = cls.download_func_dispatch(method)
        try:
            return await download_func(semaphore=semaphore, method=method, info=info, path=path, rate=rate)
        except RetryError:
            cls.logger.error(f"Retry failed for: {info}")
            return
        except Exception:
            cls.logger.error(f"Unexpected Exception for: {info}")
            raise

    @classmethod
    @unsync
    async def unsync_tasks(cls, tasks):
        return [await fob for fob in track(asyncio.as_completed(tasks), total=len(tasks))]

    @classmethod
    def single_task(cls, task, semaphore, to_do_func: Optional[Callable] = None, rate: float = 1.5) -> Unfuture:
        method, info, path = task
        task = cls.fetch_file(semaphore=semaphore, method=method, info=info, path=path, rate=rate)
        if to_do_func is not None:
            task = task.then(to_do_func)
        return task

    @classmethod
    def multi_tasks(cls, tasks: Union[Iterable, Iterator],
                    to_do_func: Optional[Callable] = None,
                    concur_req: int = 4, rate: float = 1.5,
                    ret_res: bool = True,
                    semaphore=None
                    ):

        if semaphore is None:
            # asyncio.Semaphore(concur_req)
            semaphore = init_semaphore(concur_req).result()
        else:
            cls.logger.debug(
                f'{cls.multi_tasks.__qualname__}: pass {repr(semaphore)}')
        if to_do_func is None:
            tasks = [cls.fetch_file(semaphore, method, info, path, rate)
                     for method, info, path in tasks]
        else:
            tasks = [cls.fetch_file(semaphore, method, info, path, rate).then(
                to_do_func) for method, info, path in tasks]
        if ret_res:
            t0 = perf_counter()
            res = cls.unsync_tasks(tasks).result()
            elapsed = perf_counter() - t0
            cls.logger.info(
                '{} chunks downloaded in {:.2f}s'.format(len(res), elapsed))
            return res
        else:
            return tasks

    '''
    @classmethod
    def main(cls, workdir: str, data: Union[Iterable, Iterator], concur_req: int = 4, rate: float = 1.5, logName: str = 'UnsyncFetch'):
        cls.set_logging_fileHandler(os.path.join(
            workdir, f'{logName}.log'), logName=cls.__name__)
        t0 = perf_counter()
        # res = asyncio.run(cls.multi_tasks(data, concur_req=concur_req, rate=rate))
        res = cls.multi_tasks(data, concur_req=concur_req, rate=rate).result()
        elapsed = perf_counter() - t0
        cls.logger.info(f'downloaded in {elapsed}s')
        return res
    '''

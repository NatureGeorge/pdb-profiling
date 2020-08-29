# @Created Date: 2020-08-27 10:35:34 am
# @Filename: api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-08-27 10:35:41 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pdb_profiling.fetcher.webfetch import UnsyncFetch
from pdb_profiling.utils import pipe_out
from pdb_profiling.processers.pdbe.api import PDBeDecoder
from typing import Union, Dict, Generator
from pathlib import Path
import logging
from unsync import unsync
from aiofiles import open as aiofiles_open
from orjson import loads as orjson_loads

BASE_URL: str = 'https://swissmodel.expasy.org/'


class SMR(object):
    '''
    Implement SWISS-MODEL Repository API

        * <https://swissmodel.expasy.org/docs/smr_openapi>
        * <https://swissmodel.expasy.org/docs/repository_help#smr_api>
    '''

    root = 'repository/uniprot/'
    headers = {'accept': 'text/plain'}
    use_existing: bool = False

    @classmethod
    def yieldTasks(cls, unps, params: Dict, file_format: str, folder: Union[Path, str]) -> Generator:
        for unp in unps:
            args = dict(
                url=f'{BASE_URL}{cls.root}{unp}.{file_format}',
                headers=cls.headers,
                params=params)
            yield 'get', args, Path(folder)/f'{unp}.{file_format}'

    @classmethod
    def retrieve(cls, unps, params: Dict, folder: Union[Path, str], concur_req: int = 20, rate: float = 1.5, file_format: str = 'json', ret_res: bool = True, **kwargs):
        assert file_format in ('json', 'pdb'), "Invalid file format"
        res = UnsyncFetch.multi_tasks(
            cls.yieldTasks(unps, params, file_format, folder),
            cls.process,
            concur_req=concur_req,
            rate=rate,
            logger=logging,
            ret_res=ret_res,
            semaphore=kwargs.get('semaphore', None))
        return res

    @classmethod
    def single_retrieve(cls, unp: str, params: Dict, folder: Union[Path, str], semaphore, rate: float = 1.5, file_format: str = 'json'):
        assert file_format in ('json', 'pdb'), "Invalid file format"

        return UnsyncFetch.single_task(
            task=next(cls.yieldTasks((unp, ), params, file_format, folder)),
            semaphore=semaphore,
            to_do_func=cls.process,
            rate=rate)

    @classmethod
    @unsync
    async def process(cls, path):
        logging.debug('Start to decode SMR JSON')
        if not isinstance(path, (str, Path)):
            path = await path
        if path is None or str(path).endswith('.pdb'):
            return path
        new_path = str(path).replace('.json', '.tsv')
        if Path(new_path).exists() and cls.use_existing:
            return new_path
        async with aiofiles_open(path) as inFile:
            try:
                data = orjson_loads(await inFile.read())
            except Exception as e:
                logging.error(f"Error in {path}")
                raise e
        res = PDBeDecoder.pyexcel_io(suffix='swissmodel/repository/uniprot/', data=data)
        if res is not None:
            await pipe_out(
                df=res, path=new_path,
                format='tsv', mode='w')
            logging.debug(f'Decoded file in {new_path}')
            return new_path
        else:
            logging.warning(f"Without Expected Data (swissmodel/repository/uniprot/): {data}")
            return None

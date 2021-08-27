# @Created Date: 2020-09-26 02:51:49 pm
# @Filename: api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-09-26 02:51:59 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from typing import Dict, Tuple, Iterable, Generator, Union
from pathlib import Path
from pdb_profiling.log import Abclog
from pdb_profiling.fetcher.webfetch import UnsyncFetch

BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


class EutilsAPI(Abclog):
    '''
    Implement The Entrez Programming Utilities (E-utilities) API

    * Entrez Programming Utilities Help: <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
    * DEMO URL 1: <https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=NP_001291289&rettype=fasta>
    * DEMO URL 2: <https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id=NM_001304360&rettype=fasta>
    '''
    headers = {"Content-Type": "text/plain"}
    api_set = frozenset(('efetch.fcgi', 'einfo.fcgi', 'esearch.fcgi',
                         'epost.fcgi', 'esummary.fcgi', 'egquery.fcgi'))

    @classmethod
    def dumpsParams(cls, params: Dict) -> str:
        return '&'.join(f'{key}={value}' for key, value in params.items())

    @classmethod
    def task_unit(cls, suffix: str, params: Dict, folder: Path) -> Tuple:
        args = dict(
            url=f'{BASE_URL}{suffix}',
            headers=cls.headers,
            params=params)
        return 'get', args, folder/f'{cls.dumpsParams(params)}.{params.get("retmode", params.get("rettype", "xml"))}'

    @classmethod
    def yieldTasks(cls, suffix: str, params_collection: Iterable[Dict], folder: Path) -> Generator:
        for params in params_collection:
            yield cls.task_unit(suffix, params, folder)

    """@classmethod
    def retrieve(cls, suffix: str, params_collection: Iterable[Dict], folder: Union[Path, str], concur_req: int = 20, rate: float = 1.5, ret_res: bool = True, **kwargs):
        assert suffix in cls.api_set, f"Invalid suffix! Valid set is \n{cls.api_set}"
        folder = Path(folder)
        res = UnsyncFetch.multi_tasks(
            cls.yieldTasks(suffix, params_collection, folder),
            concur_req=concur_req,
            rate=rate,
            ret_res=ret_res,
            semaphore=kwargs.get('semaphore', None))
        return res"""

    @classmethod
    def single_retrieve(cls, suffix: str, params: Dict, folder: Union[Path, str], semaphore, rate: float = 1.5):
        assert suffix in cls.api_set, f"Invalid suffix! Valid set is \n{cls.api_set}"
        folder = Path(folder)
        return UnsyncFetch.single_task(
            task=cls.task_unit(suffix, params, folder),
            semaphore=semaphore,
            rate=rate)

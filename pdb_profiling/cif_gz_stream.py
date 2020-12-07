# @Created Date: 2020-10-13 04:38:16 pm
# @Filename: cif_gz_stream.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-10-13 04:38:24 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from aiohttp import ClientSession
from aiofiles import open as aiofiles_open
from unsync import unsync
import zlib
from typing import Iterable
from pdb_profiling.utils import MMCIF2DictPlus, init_semaphore, unsync_run
from pdb_profiling.fetcher.webfetch import UnsyncFetch
from gzip import open as gzip_open

'''
>>> pdb_id, assembly_id = '2o2q', 2
>>> header = f'{pdb_id}-assembly-{assembly_id}'

>>> writer(
        reader(f'http://www.ebi.ac.uk/pdbe/static/entry/download/{header}.cif.gz'), 
        f'{header}-pdbe_chain_remapping.cif',
        b'data_%s\n#\nloop_\n' % bytes(header, 'utf-8'),
        b'_pdbe_chain_remapping'
        ).result()

>>> parser(
        reader(f'http://www.ebi.ac.uk/pdbe/static/entry/download/{header}.cif.gz'),
        ('data_%s\n' % header, '#\n', 'loop_\n'),
        b'_pdbe_chain_remapping'
        )
>>> 
'''

semaphore = unsync_run(init_semaphore(10))


def iter_index(text, target, add):
    '''
    >>> text = b'sdgfsd\nsdgsdg\nfdsg\nd'
    >>> index = (None, *iter_index(text), None)
    >>> print(index)
    >>> tuple(text[start:end] for start,end in zip(index, index[1:]))
    '''
    text_len = len(text)
    start = -1
    while True:
        try:
            res = text.index(target, start+1) + add
            yield res
            start = res
            if start >= text_len-1:
                break
        except ValueError:
            break


@unsync
def full_io(url, path, remove=True):
    path = UnsyncFetch.fetch_file(semaphore, 'get', dict(url=url), path, 1).result()
    with gzip_open(path, 'rt') as handle:
        mmcif_dict = MMCIF2DictPlus(handle, ('_pdbe_chain_remapping.',))
    if remove:
        path.unlink()
    return mmcif_dict


async def reader(url):
    dec = zlib.decompressobj(32 + zlib.MAX_WBITS)
    remain_part = b''
    async with semaphore:
        async with ClientSession() as session:
            async with session.get(url=url, timeout=3600) as resp:
                if resp.status == 200:
                    async for chunk in resp.content.iter_any():
                        rv = dec.decompress(chunk)
                        if rv:
                            # yield rv
                            index = (None, *iter_index(rv, b'\n', 1), None)
                            if len(index) == 2:
                                remain_part += rv
                                continue
                            if remain_part:
                                yield remain_part + rv[:index[1]]
                                remain_part = b''
                                for start, end in zip(index[1:-1], index[2:-1]):
                                    yield rv[start:end]
                            else:
                                for start, end in zip(index[:-1], index[1:-1]):
                                    yield rv[start:end]
                            if index[-2] != len(rv):
                                remain_part = rv[index[-2]:]
                    if remain_part:
                        yield remain_part
                else:
                    raise Exception(
                        "code={resp.status}, message={resp.reason}, headers={resp.headers}".format(resp=resp) + \
                        f"\nurl={url}")


@unsync
async def writer(handle, path, header:bytes, start_key:bytes):
    start = False
    async with aiofiles_open(path, 'wb') as fileOb:
        if header:
            await fileOb.write(header)
        async for line in handle:
            if line.startswith(start_key):
                start = True
            elif start and line.startswith(b'#'):
                await fileOb.write(line)
                return path
            if start:
                await fileOb.write(line)


def parser(handle, header_lines: Iterable[bytes], start_key: bytes):
    @unsync
    async def a_handle(handle):
        start = False
        res = []
        async for line in handle:
            if line.startswith(start_key):
                start = True
            elif start and line.startswith(b'#'):
                res.append(line.decode('utf-8'))
                return res
            if start:
                res.append(line.decode('utf-8'))

    return MMCIF2DictPlus(tuple(header_lines)+tuple(a_handle(handle).result()))

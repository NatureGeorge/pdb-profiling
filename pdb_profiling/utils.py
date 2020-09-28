# @Created Date: 2019-11-20 06:33:41 pm
# @Filename: utils.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-02-17 10:25:37 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import os
import gzip
import shutil
from typing import Optional, Union, Dict, Tuple, Iterable, Iterator, List, Coroutine, NamedTuple, Callable
from logging import Logger
from pandas import read_csv, DataFrame, isna
import numpy as np
from pathlib import Path
from aiofiles import open as aiofiles_open
from re import compile as re_compile
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
    async with aiofiles_open(path, read_mode) as file_io:
        with StringIO(await file_io.read()) as text_io:
            return read_csv(text_io, **kwargs)


async def a_load_json(path):
    if isinstance(path, (Coroutine, Unfuture)):
        path = await path
    if path is None:
        return None
    async with aiofiles_open(path) as inFile:
        return json.loads(await inFile.read())


async def pipe_out(df, path, **kwargs):
    if not isinstance(df, (DataFrame, Dataset, Sheet)):
        raise TypeError(f"Invalid Object for pipe_out(): {type(df)}")
    if len(df) == 0:
        raise ValueError("Zero record!")
    path = Path(path)
    mode = kwargs.get('mode', 'a')
    clear_headers:bool = path.exists() and mode.startswith('a')
    var_format = kwargs.get('format', 'tsv').lower()
    async with aiofiles_open(path, mode) as file_io:
        if isinstance(df, DataFrame):
            sorted_col = sorted(df.columns)
            if clear_headers:
                headers = None
            else:
                headers = sorted_col
            dataset = Dataset(headers=headers)
            dataset.extend(df[sorted_col].to_records(index=False))
            to_write = dataset.export(var_format, lineterminator='\n')
        elif isinstance(df, Dataset):
            if clear_headers:
                df.headers = None
            to_write = df.export(var_format, lineterminator='\n')
        elif isinstance(df, Sheet):
            if clear_headers:
                df.colnames = []
            to_write = getattr(df, f"get_{var_format}")(lineterminator='\n')
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
    
    def chainer_list(s):
        return list(chain.from_iterable(s))
    
    if mode == 'sep':
        chainer = chainer_sep
        lens = df[cols_to_split[0]].str.split(sep).map(len)
    elif mode == 'json-list':
        chainer = chainer_json
        lens = df[cols_to_split[0]].apply(json.loads).map(len)
    elif mode == 'list':
        chainer = chainer_list
        lens = df[cols_to_split[0]].map(len)
    else:
        raise ValueError("Invalid mode!")
    
    repeat_dict = {col: np.repeat(df[col], lens)
                   for col in set(all_cols)-set(cols_to_split)}
    chain_dict = {col: chainer(df[col]) for col in cols_to_split}
    return DataFrame({**repeat_dict, **chain_dict})


def init_folder_from_suffix(folder: Union[Path, str], suffix: str):
    folder = Path(folder)
    new_path = folder/suffix
    new_path.mkdir(parents=True, exist_ok=True)
    return new_path


def init_folder_from_suffixes(folder: Union[Path, str], suffixes: Iterable) -> Iterable[Path]:
    folder = Path(folder)
    for suffix in suffixes:
        new_path = folder/suffix
        new_path.mkdir(parents=True, exist_ok=True)
        yield new_path


def iter_first(df: DataFrame, criteria: Callable[[NamedTuple], bool], **kwargs) -> Optional[NamedTuple]:
    '''
    Implement pandas.DataFrame.itertuples
    
    Returns the value as soon as you find the first row/record 
    that meets the requirements and NOT iterating other rows

    Originated from: https://stackoverflow.com/a/63826677/12876491

    >>> iter_first(df, lambda row: row.A > 4 and row.B > 3)
    '''
    for row in df.itertuples(**kwargs):
        if criteria(row):
            return row

class MMCIF2DictPlus(dict):
    """
    Parse a mmCIF file and return a dictionary
        
    NOTE: Override methods based on Biopython's `Bio.PDB.MMCIF2Dict.MMCIF2Dict`
    """

    def _check_token_with_focus_keys(self, token: Iterable[str]) -> bool:
        is_key, key_value = token
        return is_key == 0 and ((key_value in self.focus_keys) or any(key in key_value for key in self.focus_keys))

    def __init__(self, handle, focus_keys: Iterable[str]=['']):
        self.focus_keys = set(focus_keys)
        self.quote_chars = ["'", '"']
        self.whitespace_chars = [" ", "\t"]
        # TODO: init first loop
        loop_flag = False
        key = None
        tokens = self._tokenize(handle)
        try:
            token = next(tokens)
        except StopIteration:
            return  # NOTE: annotation from biopython: for Python 3.7 and PEP 479
        self[token[1][0:5]] = token[1][5:]
        i = 0
        n = 0
        use = []
        # TODO: loops
        for token in tokens:
            if token[1].lower() == "loop_":
                loop_flag = True
                keys = []
                i = 0
                n = 0
                use = []
                continue
            elif loop_flag:
                '''
                NOTE: annotation from biopython:
                # The second condition checks we are in the first column
                # Some mmCIF files (e.g. 4q9r) have values in later columns
                # starting with an underscore and we don't want to read
                # these as keys
                '''
                if token[1].startswith("_") and (n == 0 or i % n == 0):
                    if i > 0:
                        loop_flag = False
                    else:
                        if self._check_token_with_focus_keys(token):
                            use.append(n)
                            self[token[1]] = []
                        keys.append(token[1])
                        n += 1
                        continue
                else:
                    key_index = i % n
                    try:
                        if key_index in use:
                            self[keys[key_index]].append(token[1])
                    except Exception:
                        raise ValueError(f"{keys}, {key_index}, {use}")
                    i += 1
                    continue
            if key is None:
                if self._check_token_with_focus_keys(token):
                    key = token[1]
            else:
                # Always returns a list
                self[key] = [token[1]]
                key = None
            if self.keys() >= self.focus_keys:
                break

    def _splitline(self, line: str):
        # NOTE: annotation from biopython: See https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for the syntax
        in_token = False
        # NOTE: annotation from biopython: quote character of the currently open quote, or None if no quote open
        quote_open_char = None
        start_i = 0
        for (i, c) in enumerate(line):
            if c in self.whitespace_chars:
                if in_token and not quote_open_char:
                    in_token = False
                    yield start_i, line[start_i:i]
            elif c in self.quote_chars:
                if not quote_open_char and not in_token:
                    # raise ValueError(f"{self['data_']}: Opening quote in middle of word: " + line)
                    quote_open_char = c
                    in_token = True
                    start_i = i + 1
                elif c == quote_open_char and (i + 1 == len(line) or line[i + 1] in self.whitespace_chars):
                    quote_open_char = None
                    in_token = False
                    yield start_i, line[start_i:i]
            elif c == "#" and not in_token:
                ''' NOTE: annotation from biopython:
                # Skip comments. "#" is a valid non-comment char inside of a
                # quote and inside of an unquoted token (!?!?), so we need to
                # check that the current char is not in a token.
                '''
                return
            elif not in_token:
                in_token = True
                start_i = i
        if in_token:
            yield start_i, line[start_i:]
        if quote_open_char:
            raise ValueError("Line ended with quote open: " + line)

    def _tokenize(self, handle):
        empty = True
        for line in handle:
            empty = False
            if line.startswith("#"):
                continue
            elif line.startswith(";"):
                '''
                NOTE: annotation from biopython:
                # The spec says that leading whitespace on each line must be
                # preserved while trailing whitespace may be stripped.  The
                # trailing newline must be stripped.
                '''
                token_buffer = [line[1:].rstrip()]
                for line in handle:
                    line = line.rstrip()
                    if line.startswith(";"):
                        yield 1, "\n".join(token_buffer)
                        line = line[1:]
                        if line and not line[0] in self.whitespace_chars:
                            raise ValueError("Missing whitespace")
                        break
                    token_buffer.append(line)
                else:
                    raise ValueError("Missing closing semicolon")
            yield from self._splitline(line.strip())
        if empty:
            raise ValueError("Empty file.")


class DisplayPDB(object):

    a_name = 'Asymmetric unit'

    b_name = 'Biological assembly {assembly_id}'

    a_code = 'model-1'

    b_code = 'assembly-{assembly_id}'

    header_unit = '''
                <td>
                    <b>{name}</b> of {pdb_id}
                </td>
    '''

    content_unit = '''
                <td>
                    <img width="300em" src="https://cdn.rcsb.org/images/structures/{in_code}/{pdb_id}/{pdb_id}_{code}.jpeg"/>
                </td>
    '''

    template = '''
        <table align="center">
            <tr>
            {headers}
            </tr>
            <tr>
            {content}
            </tr>
        </table>
        '''

    @classmethod
    def setting(cls, pdb_id, assemblies):
        headers = [cls.header_unit.format(name=cls.a_name, pdb_id=pdb_id)]
        content = [cls.content_unit.format(pdb_id=pdb_id, in_code=pdb_id[1:3], code=cls.a_code)]
        for assembly_id in assemblies:
            headers.append(cls.header_unit.format(
                name=cls.b_name.format(assembly_id=assembly_id),
                pdb_id=pdb_id))
            content.append(cls.content_unit.format(
                pdb_id=pdb_id, 
                in_code=pdb_id[1:3], 
                code=cls.b_code.format(assembly_id=assembly_id)
            ))
        return ''.join(headers), ''.join(content)

    def display(self, pdb_id, assemblies: Iterable[int]= [1]):
        from IPython.display import display, HTML
        assemblies = sorted(int(i) for i in assemblies if int(i) > 0)
        headers, content = self.setting(pdb_id, assemblies)
        self.table = self.template.format(headers=headers, content=content)
        display(HTML(self.table))
    
    def __init__(self, *args):
        if len(args) > 0:
            self.display(*args)


fasta_pat = re_compile(r'(>.+)\n([A-Z\*\n]+)')
unp_header_pat = re_compile(r'>sp\|(.+)\|')

async def a_seq_reader(path: Union[Unfuture, Union[Path, str]]):
    if isinstance(path, Unfuture):
        path = await path
    async with aiofiles_open(path, 'rt') as handle:
        header, content = fasta_pat.match(await handle.read()).groups()
        return header, content.replace('\n', '')


@unsync
async def get_seq_from_parser(res, identifier, seq_only:bool = True):
    async for header, content in await res:
        if identifier in header:
            return content if seq_only else (header, content)


@unsync
async def get_seqs_from_parser(res, identifiers:Optional[Iterable[str]]=None):
    ret = []
    async for header, content in await res:
        header = unp_header_pat.match(header).group(1)
        if identifiers is None or header in identifiers:
            ret.append((header, content))
    return ret


async def a_seq_parser(path: Union[Unfuture, Union[Path, str]]):
    if isinstance(path, Unfuture):
        path = await path
    async with aiofiles_open(path, 'rt') as handle:
        header, content = None, ''
        async for line in handle:
            if line.startswith('>'):
                if header is not None:
                    yield header.strip(), content.replace('\n', '')
                header, content = line, ''
            else:
                content += line
        yield header.strip(), content.replace('\n', '')

nu2aa_dict = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

def translate2aa(seq:str, check:bool=False):
    assert len(seq) % 3 == 0, "Invalid length of dna OR rna sequence!"
    seq = seq.replace('U', 'T')
    p_seq = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        p_seq += nu2aa_dict[codon]
    if check:
        assert "_" not in p_seq, "Invalid Sequence!"
    return p_seq


def unsync_run(arg):
    @unsync
    async def unsync_wrap(arg):
        return await arg 
    return unsync_wrap(arg).result()

# @Created Date: 2019-11-20 06:33:41 pm
# @Filename: utils.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2021-12-28 07:07:34 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import os
import gzip
import shutil
from inspect import isawaitable
from typing import Optional, Union, Dict, Tuple, Iterable, Iterator, List, Coroutine, NamedTuple, Callable, Generator
from logging import Logger
from pandas import read_csv, DataFrame, isna, Series, concat
import numpy as np
from pathlib import Path
from aiofiles import open as aiofiles_open
from re import compile as re_compile
from pyexcel import Sheet
import asyncio
from unsync import unsync, Unfuture
from itertools import chain
import orjson as json
import logging
from io import StringIO
from operator import itemgetter
from textdistance import overlap, sorensen
from collections import Counter, OrderedDict
from warnings import warn


"""def to_interval(lyst: Union[Iterable, Iterator]) -> List:
    def pass_check(lyst):
        try:
            if (not isinstance(lyst, Generator) and  (len(lyst) == 0)) or isna(lyst):
                return False
            else:
                return True
        except ValueError:
            if isinstance(lyst, float):
                return False
            else:
                return True
    if not pass_check(lyst):
        return ()
    else:
        if not isinstance(lyst, set):
            lyst = frozenset(int(i) for i in lyst if i is not None)
        if not pass_check(lyst):
            return ()
        start = []
        interval_lyst = []
        max_edge = max(lyst)
        min_edge = min(lyst)
        if len(lyst) == (max_edge + 1 - min_edge):
            return ((min_edge, max_edge),)
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
            return tuple((min(li), max(li)) for li in interval_lyst)


def lyst22interval(x, y):
    # x, y = sorted(x), sorted(y)
    data = frozenset({i for i in zip(x,y)})
    x, y = zip(*sorted(data, key=itemgetter(0)))
    start_x, start_y = x[0], y[0]
    index_x, index_y = x[0]-1, y[0]-1
    interval_x, interval_y = [], []
    for i, j in zip(x, y):
        pre_x = index_x + 1
        pre_y = index_y + 1
        if pre_x == i and pre_y == j:
            index_x, index_y = i, j
        else:
            interval_x.append((start_x, index_x))
            interval_y.append((start_y, index_y))
            start_x, start_y = i, j
            index_x, index_y = i, j
    interval_x.append((start_x, index_x))
    interval_y.append((start_y, index_y))
    return interval_x, interval_y"""


@unsync
async def init_semaphore(concurreq) -> Unfuture:
    """
    `semaphore` initiated in the `unsync` event loop
    """
    await asyncio.sleep(.01)
    return asyncio.Semaphore(concurreq)


'''
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
'''


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


@unsync
async def a_read_csv(path, read_mode='r',**kwargs):
    '''
    only suitable for small dataframe
    '''
    try:
        if isinstance(path, (Coroutine, Unfuture)):
            path = await path
        if isinstance(path, (Path, str)):
            async with aiofiles_open(path, read_mode) as file_io:
                with StringIO(await file_io.read()) as text_io:
                    return read_csv(text_io, **kwargs)
        elif isinstance(path, DataFrame):
            return path
        else:
            return DataFrame(path, columns=kwargs.get('columns', None))
    except Exception:
        raise ValueError(f'{path}')


async def a_load_json(path):
    try:
        if isawaitable(path):
            path = await path
        if path is None:
            return None
        async with aiofiles_open(path) as inFile:
            return json.loads(await inFile.read())
    except Exception as e:
        raise ValueError(f"Exception {e} ({path})")


async def pipe_out(df, path, **kwargs):
    if not isinstance(df, Sheet):
        raise TypeError(f"Invalid Object for pipe_out(): {type(df)}")
    if len(df) == 0:
        raise ValueError("Zero record!")
    path = Path(path)
    mode = kwargs.get('mode', 'a')
    clear_headers:bool = path.exists() and mode.startswith('a')
    var_format = kwargs.get('format', 'tsv').lower()
    async with aiofiles_open(path, mode) as file_io:
        """
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
        """
        df = df.project(sorted(df.colnames))
        if clear_headers:
            df.colnames = []
        to_write = getattr(df, f"get_{var_format}")(lineterminator='\n')
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
                    <img class="display" width="300em" src="https://cdn.rcsb.org/images/structures/{in_code}/{pdb_id}/{pdb_id}_{code}.jpeg"/>
                </td>
    '''

    css = '''
        <style>
            img.display {
                -webkit-filter: invert(1);
                filter: invert(1);
                }
        </style>
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

    def show(self, pdb_id, assemblies: Iterable[int]= [1]):
        from IPython.display import display, HTML
        assemblies = sorted(int(i) for i in assemblies if int(i) > 0)
        headers, content = self.setting(pdb_id, assemblies)
        self.table = self.template.format(headers=headers, content=content)
        if self.dark:
            self.table = self.css + self.table
        display(HTML(self.table))
    
    def __init__(self, dark:bool=False):
        self.dark = dark


fasta_pat = re_compile(r'(>.+)\n([A-Z\*\n]+)')
unp_header_pat = re_compile(r'>sp\|(.+)\|')

async def a_seq_reader(path: Union[Unfuture, Union[Path, str]]):
    if isinstance(path, Unfuture):
        path = await path
    async with aiofiles_open(path, 'rt') as handle:
        header, content = fasta_pat.match(await handle.read()).groups()
        content = content.replace('\n', '')
        assert content != '', str(path)
        return header, content


@unsync
async def get_seq_from_parser(res, identifier, seq_only:bool = True):
    async for header, content in await res:
        if (identifier in header) and (f'{identifier}-' not in header):
            return content if seq_only else (header, content)


@unsync
async def get_seqs_from_parser(res, identifiers:Optional[Iterable[str]]=None):
    ret = []
    async for header, content in await res:
        header = unp_header_pat.match(header).group(1)
        if identifiers is None or header in identifiers:
            ret.append((header, content))
    return ret


async def a_seq_parser(path: Union[Unfuture, Coroutine, Path, str]):
    if isawaitable(path):
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
        assert header is not None, f"\npath: {path}\ncur_content: {content}"
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
    async def s_unsync_wrap(arg):
        return await arg 
    return s_unsync_wrap(arg).result()


@unsync
def unsync_wrap(var):
    return var


class SeqRangeReader(object):
    def __init__(self, name_group):
        self.name = name_group  # ('pdb_id', 'chain_id', 'UniProt')
        self.pdb_range = []
        self.unp_range = []

    def output(self):
        if self.pdb_range:
            pdb_range = json.dumps(self.pdb_range).decode('utf-8')
            unp_range = json.dumps(self.unp_range).decode('utf-8')
            return pdb_range, unp_range
        else:
            return self.default_pdb_range, self.default_unp_range

    def check(self, name_group_to_check, data_group):
        self.default_pdb_range = '[[%s, %s]]' % data_group[:2]
        self.default_unp_range = '[[%s, %s]]' % data_group[2:4]

        if self.name == name_group_to_check:
            self.pdb_range.append([int(data_group[0]), int(data_group[1])])
            self.unp_range.append([int(data_group[2]), int(data_group[3])])
        else:
            self.name = name_group_to_check
            self.pdb_range = [[int(data_group[0]), int(data_group[1])]]
            self.unp_range = [[int(data_group[2]), int(data_group[3])]]

        return self.output()


def sort_2_range(unp_range: List, pdb_range: List):
    unp_range, pdb_range = zip(
        *sorted(zip(unp_range, pdb_range), key=lambda x: x[0][0]))
    return unp_range, pdb_range


def flat_dict_in_df(dfrm:DataFrame, targetCol:Union[str, Series], cols:List):
    try:
        new_cols = list(f'{targetCol.name}.{col}' for col in cols)
        dfrm[new_cols] = DataFrame(
            targetCol.apply(lambda x: list(x[col] for col in cols)).to_list(),
            columns=new_cols)
        return dfrm.drop(columns=[targetCol.name])
    except AttributeError:
        assert isinstance(targetCol, str)
        new_cols = list(f'{targetCol}.{col}' for col in cols)
        dfrm[cols] = DataFrame(
            dfrm[targetCol].apply(lambda x: list(x[col] for col in cols)).to_list(),
            columns=new_cols)
        return dfrm.drop(columns=[targetCol])


class AAThree2One(dict):
    def __missing__(self, key):
        return 'X'

aa_three2one = AAThree2One({
    "GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C", "VAL": "V", "LEU": "L",
    "ILE": "I", "MET": "M", "PRO": "P", "PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D",
    "GLU": "E", "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R"})


standardAA = list(aa_three2one.keys())

standardNu = ['DA', 'DT', 'DC', 'DG', 'DI', 'A', 'U', 'C', 'G', 'I']


"""def range_len(lyst: Union[List, str, float]) -> int:
    if isinstance(lyst, float) or lyst is None:
        return 0
    elif isinstance(lyst, str):
        lyst = json.loads(lyst)
    length = 0
    for left, right in lyst:
        assert right >= left, f"\n{lyst}"
        length += right - left + 1
    return length


def interval2set(lyst: Union[Iterable, Iterator, str]):
    if isinstance(lyst, str):
        lyst = json.loads(lyst)
    range_set = frozenset()
    for left, right in lyst:
        range_set |= frozenset(range(left, right+1))
    return range_set"""


def expand_interval(lyst: Union[Iterable, Iterator, str]):
    lyst = json.loads(lyst) if isinstance(lyst, str) else lyst
    yield from (i for start, end in lyst for i in range(start, end+1))


def lyst2range(lyst, add_end=1):
    if isinstance(lyst, str):
        lyst = json.loads(lyst)
    for start, end in lyst:
        yield from range(int(start), int(end)+add_end)


"""def subtract_range(pdb_range: Union[str, Iterable], mis_range: Union[str, Iterable]) -> List:
    if isinstance(mis_range, float) or mis_range is None:
        return pdb_range
    elif len(pdb_range) == 0:
        return ()
    elif len(mis_range) == 0:
        return pdb_range
    pdb_range_set = interval2set(pdb_range)
    mis_range_set = interval2set(mis_range)
    return to_interval(pdb_range_set - mis_range_set)


def check_range(i) -> bool:
    if isinstance(i, float) or (i is None) or (len(i) == 0) or (i == 'nan'):
        return False
    return True


def add_range(left: Union[str, Iterable], right: Union[str, Iterable]) -> List:
    check_left = check_range(left)
    check_right = check_range(right)
    if check_left and not check_right:
        return left
    elif not check_left and check_right:
        return right
    elif not check_left and not check_right:
        return np.nan
    try:
        left_range_set = interval2set(left)
        right_range_set = interval2set(right)
        return to_interval(left_range_set | right_range_set)
    except Exception as e:
        print(left, right)
        print(type(left), type(right))
        raise e


def overlap_range(obs_range:Union[str, Iterable], unk_range: Union[str, Iterable]) -> List:
    if isinstance(unk_range, float) or unk_range is None:
        return ()
    '''
    obs_range_set = interval2set(obs_range)
    unk_range_set = interval2set(unk_range)
    return to_interval(obs_range_set.intersection(unk_range_set))
    '''
    obs_range = json.loads(obs_range) if isinstance(obs_range, str) else obs_range
    unk_range = json.loads(unk_range) if isinstance(unk_range, str) else unk_range

    def unit(i1,i2):
        for start1, end1 in i1:
            for start2, end2 in i2:
                sl = start2 >= start1
                sr = start2 <= end1
                el = end2 >= start1
                er = end2 <= end1
                s_in = sl and sr
                e_in = el and er
                ini = s_in or e_in
                # out = (sl and el) or (sr and er)
                cov = (not sl) and (not er)
                start = start2 if s_in else start1
                end = end2 if e_in else end1
                if ini or cov:
                    yield start, end

    return tuple(unit(obs_range, unk_range))"""


def get_seq_seg(seq, ranges, **kwargs):
    for start,end in ranges:
        if end >= start:
            yield start, seq[start-1:end]
        else:
            warn(f"{kwargs} -> Invalid Order: {ranges}, skip")


def get_diff_index(lseq, lrange, rseq, rrange):
    if isinstance(lrange, str):
        lrange = json.loads(lrange)
    if isinstance(rrange, str):
        rrange = json.loads(rrange)
    for (lstart, lseg), (rstart, rseg) in zip(get_seq_seg(lseq, lrange), get_seq_seg(rseq, rrange)):
        yield from ((lstart+index, rstart+index) for index, (r1, r2) in enumerate(zip(lseg, rseg)) if r1 != r2)


def red_seq_seg(seq, ranges):
    edge = 0
    for start, end in ranges:
        yield f"{seq[edge:start-1]}\x1b[31m{seq[start-1:end]}\x1b[0m"
        edge = end
    yield seq[end:]


"""def outside_range(pdb_range: Union[str, Iterable], seqres_len: int):
    pdb_range = json.loads(pdb_range) if isinstance(pdb_range, str) else pdb_range
    out_head = pdb_range[0][0] - 1
    out_tail = pdb_range[-1][-1] + 1
    ret = [[1, out_head], [out_tail, seqres_len]]
    return [i for i in ret if i[0]<=i[1]]"""


def outside_range_len(pdb_range: Union[str, Iterable], seqres_len: int, omit: int = 5) -> int:
    if isinstance(pdb_range, str):
        lyst = json.loads(pdb_range)
    else:
        lyst = pdb_range
    out_head = lyst[0][0]-1
    out_tail = seqres_len - lyst[-1][-1]
    if out_head <= omit:
        out_head = 0
    else:
        out_head -= omit
    if out_tail <= omit:
        out_tail = 0
    else:
        out_tail -= omit
    return out_head + out_tail


def get_gap_list(li: Union[str,List,Tuple]):
    if isinstance(li, str):
        li = json.loads(li)
    return [li[i+1][0] - li[i][1] - 1 for i in range(len(li)-1)]


def get_range_diff(lyst_a: Union[str, List, Tuple], lyst_b: Union[str, List, Tuple]):
    lyst_a = json.loads(lyst_a) if isinstance(lyst_a, str) else lyst_a
    lyst_b = json.loads(lyst_b) if isinstance(lyst_b, str) else lyst_b
    array_a = np.array([right - left + 1 for left, right in lyst_a])
    array_b = np.array([right - left + 1 for left, right in lyst_b])
    return array_a - array_b


def select_range(ranges, indexes, cutoff=0.2, skip_index=[], selected_ranges=None, similarity_func=overlap.similarity):
    select_index = []
    selected_ranges = [] if selected_ranges is None else selected_ranges
    def unit(cur_index):
        if cur_index in skip_index:
            return
        cur_range = ranges[cur_index]
        cur_range = json.loads(cur_range) if isinstance(cur_range, str) else cur_range
        for selected_range in selected_ranges:
            selected_range = json.loads(selected_range) if isinstance(selected_range, str) else selected_range
            if len(cur_range) == 0:
                return
            score = similarity_func(lyst2range(cur_range),
                               lyst2range(selected_range))
            if score > cutoff:
                return
        select_index.append(cur_index)
        selected_ranges.append(cur_range)

    for index in indexes:
        unit(index)
    return select_index


"""def select_ho_range(ranges1, ranges2, indexes, cutoff=0.2, skip_index=[]):
    from scipy.stats import wasserstein_distance
    select_index = []

    def unit(cur_index):
        if cur_index in skip_index:
            return
        cur_range1, cur_range2 = ranges1[cur_index], ranges2[cur_index]
        c1_1 = Counter(expand_interval(cur_range1))
        c1_2 = Counter(expand_interval(cur_range2))
        if len(c1_1) == 0 or len(c1_2) == 0:
            return
        c1 = c1_1 + c1_2
        for selected in select_index:
            selected_range1,selected_range2 = ranges1[selected], ranges2[selected]
            c_c1 = c1.copy()
            c2 = Counter(expand_interval(selected_range1))+Counter(expand_interval(selected_range2))
            for key in c_c1.keys() | c2.keys():
                if key not in c_c1:
                    c_c1[key] = 0
                if key not in c2:
                    c2[key] = 0
            oc2 = OrderedDict(sorted((item for item in c2.items()), key=lambda x: x[0]))
            oc1 = OrderedDict(sorted((item for item in c_c1.items()), key=lambda x: x[0]))
            score = wasserstein_distance(tuple(oc1.values()), tuple(oc2.values()))
            if score < cutoff:
                return
        select_index.append(cur_index)

    for index in indexes:
        unit(index)
    return select_index"""


def select_ho_max_range(ranges1, ranges2, indexes, cutoff=0.2, skip_index=[]):
    select_range_set = OrderedDict()

    def unit(cur_index):
        if cur_index in skip_index:
            return
        cur_range1, cur_range2 = ranges1[cur_index], ranges2[cur_index]

        c1_1 = frozenset((1, i) for i in expand_interval(cur_range1))
        c1_2 = frozenset((2, i) for i in expand_interval(cur_range2))
        if len(c1_1) == 0 or len(c1_2) == 0:
            return
        c1 = c1_1 | c1_2
        for c2s in select_range_set.values():
            for c2 in c2s:
                score = sorensen.similarity(c1, c2)
                if score > cutoff:
                    return
        select_range_set[cur_index] = (
                frozenset((1, i) for i in expand_interval(cur_range1)) | frozenset((2, i) for i in expand_interval(cur_range2)),
                frozenset((2, i) for i in expand_interval(cur_range1)) | frozenset((1, i) for i in expand_interval(cur_range2)))

    for index in indexes:
        unit(index)
    return list(select_range_set.keys())


def select_he_range(Entry_1, Entry_2, ranges1, ranges2, indexes, cutoff=0.2, skip_index=[]):
    select_index = []

    def unit(cur_index):
        if cur_index in skip_index:
            return
        cur_range1, cur_range2 = ranges1[cur_index], ranges2[cur_index]
        cur_e1, cur_e2 = Entry_1[cur_index], Entry_2[cur_index]
        (cur_e1, cur_range1), (cur_e2, cur_range2) = sorted(((cur_e1, cur_range1), (cur_e2, cur_range2)), key=lambda x: x[0])

        c1_1 = frozenset((1, i) for i in expand_interval(cur_range1))
        c1_2 = frozenset((2, i) for i in expand_interval(cur_range2))
        if len(c1_1) == 0 or len(c1_2) == 0:
            return
        c1 = c1_1 | c1_2
        for selected in select_index:
            selected_range1, selected_range2 = ranges1[selected], ranges2[selected]
            selected_e1, selected_e2 = Entry_1[selected], Entry_2[selected]
            (selected_e1, selected_range1), (selected_e2, selected_range2) = sorted(((selected_e1, selected_range1), (selected_e2, selected_range2)), key=lambda x: x[0])
            c2 = frozenset((1, i) for i in expand_interval(selected_range1)) | frozenset((2, i) for i in expand_interval(selected_range2))
            score = sorensen.similarity(c1, c2)
            if score > cutoff:
                return
        select_index.append(cur_index)

    for index in indexes:
        unit(index)
    return select_index


def dumpsParams(params: Dict) -> str:
    return '&'.join(f'{key}={value}' for key, value in params.items())


@unsync
async def a_concat(pathes, sep='\t', sort=False, ignore_index=True, columns=None):
    if isinstance(pathes, (Unfuture, Coroutine)):
        pathes = await pathes
    res = [await a_read_csv((await path) if isinstance(Unfuture, Coroutine) else path, sep=sep, columns=columns) for path in pathes]
    return concat((i for i in res if i is not None), sort=sort, ignore_index=ignore_index)


def get_str_dict_len(x):
    if isinstance(x, str):
        return x.count(':')
    else:
        return len(x)


def id2score(identifier):
    len_id = len(identifier)
    return -sum(ord(i)*2*(len_id-level) for level, i in enumerate(identifier))

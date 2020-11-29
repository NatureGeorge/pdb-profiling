# @Created Date: 2020-11-23 10:29:17 am
# @Filename: command.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-11-23 10:29:36 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pdb_profiling import default_config
from pdb_profiling.commands import CustomDB
from pdb_profiling.processors import Identifiers, Identifier, SIFTSs
from pdb_profiling.utils import unsync_run
import click
from pandas import read_csv
from importlib import util as imp_util
from pathlib import Path
from math import ceil
from tqdm.auto import tqdm
from time import sleep
from random import uniform


def colorClick(name: str, template: str = "Initializing %s", fg: str = "green"):
    return click.style(template % name, fg=fg)


@click.group(chain=True, invoke_without_command=False)
@click.option("--folder", default="./", help="The output folder.", type=click.Path())
@click.option('--dropall/--no-dropall', help="whether to use existing custom DB", default=False, is_flag=True)
@click.pass_context
def Interface(ctx, folder, dropall):
    folder = Path(folder)
    click.echo(colorClick(f"Folder: {folder}"))
    ctx.ensure_object(dict)
    ctx.obj['folder'] = folder
    default_config(folder)
    ctx.obj['custom_db'] = CustomDB(
        "sqlite:///%s" % (folder/'local_db'/'custom.db'), dropall)


@Interface.command("insert-mutation")
@click.option("--input", help="the file that contains sites info", type=click.Path())
@click.option("--sep", default="\t", help="the seperator of input file", type=str)
@click.option("--usecols", default='from_id,Ref,Pos,Alt', help="The comma-sep columns of site info", type=str)
@click.option('--readchunk', type=int, help="the chunksize parameter of pandas.read_csv", default=100000)
@click.option('--nrows', type=int, help="the nrows parameter of pandas.read_csv", default=None)
@click.option('--skiprows', type=int, help="the skiprows parameter of pandas.read_csv", default=None)
@click.option('--functionfile', help="the py file that contains custom function", default=None, type=click.Path())
@click.option('--functionname', default='do_something', type=str)
@click.pass_context
def insert_sites(ctx, input, sep, usecols, readchunk, nrows, skiprows, functionfile, functionname):
    def do_nothing(dfrm):
        return dfrm.to_dict('records')

    click.echo(colorClick("DB Mutation Insertion"))
    usecols = usecols.split(',')
    if functionfile is not None:
        spec = imp_util.spec_from_file_location("CustomFunc", functionfile)
        CustomFunc = imp_util.module_from_spec(spec)
        spec.loader.exec_module(CustomFunc)
        deal = getattr(CustomFunc, functionname)
        click.echo(f"Success: load func: {functionname} from {functionfile}")
    else:
        deal = do_nothing
    df = read_csv(input, sep=sep, usecols=usecols, chunksize=readchunk,
                  nrows=nrows, skiprows=skiprows)
    sqlite_api = ctx.obj['custom_db']
    start = 0
    for index, dfrm in enumerate(df):
        end = readchunk*(index+1)
        click.echo(f"Try to insert: {start}-{end}")
        start = end+1
        sqlite_api.sync_insert(sqlite_api.Mutation, deal(dfrm))


@Interface.command("id-mapping")
@click.option('--chunksize', type=int, help="the chunksize parameter", default=500)
@click.pass_context
def id_mapping(ctx, chunksize):
    sqlite_api = ctx.obj['custom_db']
    total = unsync_run(sqlite_api.database.fetch_one(
        query="SELECT COUNT(DISTINCT ftId) FROM Mutation WHERE ftId NOT IN (SELECT DISTINCT ftId FROM IDMapping)"))[0]
    cols = ('ftId', 'Entry', 'isoform', 'is_canonical')
    click.echo(f"Total {total} to query")
    query = f"""
            SELECT DISTINCT ftId FROM Mutation
            WHERE ftId NOT IN (SELECT DISTINCT ftId FROM IDMapping)
            LIMIT {chunksize}
            """
    for _ in range(ceil(total/chunksize)):
        res = unsync_run(sqlite_api.database.fetch_all(query=query))
        if len(res) == 0:
            break
        res = Identifiers(i[0] for i in res).fetch(
            'map2unp').run(tqdm).result()
        values = [dict(zip(cols, i)) for i in res]
        if values:
            sqlite_api.sync_insert(sqlite_api.IDMapping, values)
        sleep(uniform(1, 10))


@Interface.command("sifts-mapping")
@click.option('--func', type=str, default='pipe_select_mo')
@click.option('--chunksize', type=int, help="the chunksize parameter", default=200)
@click.option('--output', type=str, default='')
@click.pass_context
def sifts_mapping(ctx, func, chunksize, output):
    def get_unp_id(args):
        Entry, isoform, is_canonical = args
        return Entry if is_canonical else isoform

    sqlite_api = ctx.obj['custom_db']
    total = unsync_run(sqlite_api.database.fetch_one(
        query="SELECT COUNT(DISTINCT isoform) FROM IDMapping WHERE isoform != 'NaN'"))[0]
    click.echo(f"Total {total} to query")
    output = f'{func}.tsv' if output == '' else output
    for i in range(ceil(total/chunksize)):
        res = unsync_run(sqlite_api.database.fetch_all(
            query=f"""
            SELECT DISTINCT Entry,isoform,is_canonical FROM IDMapping
            WHERE isoform != 'NaN'
            LIMIT {chunksize} OFFSET {chunksize*i}
            """))
        res = SIFTSs(map(get_unp_id, res)).fetch(func, skip_pdbs=(
            '6vnn', '2i6l', '4zai', '5jn1', '6bj0', '6yth')+('4fc3', '7acu', '6lsd', '6llc', '6xoz')).run(tqdm).result()
        output_path = ctx.obj['folder']/output
        for dfrm in res:
            if dfrm is None:
                continue
            dfrm.to_csv(output_path, sep='\t', index=False,
                        header=not output_path.exists(), mode='a+')
        sleep(uniform(1, 10))


'''
@Interface.command("residue-mapping")
@click.option('--chunksize', type=int, help="the chunksize parameter", default=200)
@click.option('--output', type=str)
@click.pass_context
def sifts_mapping(ctx, func, chunksize, output):
'''


if __name__ == '__main__':
    Interface(obj={})

# @Created Date: 2020-11-23 10:29:17 am
# @Filename: command.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-11-23 10:29:36 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pdb_profiling import default_config
from pdb_profiling.commands import CustomDB
from pdb_profiling.processors import Identifiers, Identifier, SIFTSs, SIFTS, PDB, PDBs
from pdb_profiling.utils import unsync_run
import click
from pandas import read_csv, concat, DataFrame
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
@click.option('--input', type=click.Path(), default=None)
@click.option('--column', type=str, default=None)
@click.option('--sep', type=str, default='\t')
@click.option('--chunksize', type=int, help="the chunksize parameter", default=500)
@click.pass_context
def id_mapping(ctx, input, column, sep, chunksize):
    sqlite_api = ctx.obj['custom_db']
    cols = ('ftId', 'Entry', 'isoform', 'is_canonical')
    if input is None:
        total = unsync_run(sqlite_api.database.fetch_one(
            query="SELECT COUNT(DISTINCT ftId) FROM Mutation WHERE ftId NOT IN (SELECT DISTINCT ftId FROM IDMapping)"))[0]
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
    else:
        if column is None:
            ids = read_csv(input, sep=sep, header=None)[0].unique()
        else:
            ids = read_csv(input, sep=sep, usecols=[column])[column].unique()
        total = len(ids)
        click.echo(f"Total {total} to query")
        for i in range(0, total, chunksize):
            res = Identifiers(ids[i:i+chunksize]).fetch(
                'map2unp').run(tqdm).result()
            values = [dict(zip(cols, i)) for i in res]
            if values:
                sqlite_api.sync_insert(sqlite_api.IDMapping, values)
            sleep(uniform(1, 10))


@Interface.command("sifts-mapping")
@click.option('--input', type=click.Path(), default=None)
@click.option('--column', type=str, default=None)
@click.option('--sep', type=str, default='\t')
@click.option('--func', type=str, default='pipe_select_mo')
@click.option('--kwargs', type=str, default='{}')
@click.option('--chunksize', type=int, help="the chunksize parameter", default=200)
@click.option('--entry_filter', type=str, default='(release_date < "20201020") and ((experimental_method in ["X-ray diffraction", "Electron Microscopy"] and resolution <= 3) or experimental_method == "Solution NMR")')
@click.option('--chain_filter', type=str, default="UNK_COUNT < SEQRES_COUNT and ca_p_only == False and identity >=0.9 and repeated == False and reversed == False and OBS_COUNT > 20")
@click.option('--skip_pdbs', type=str, default='6wrg,5jm5,6vnn,2i6l,4zai,5jn1,6bj0,6yth,4fc3,7acu,6lsd,6llc,6xoz,6xp0,6xp1,6xp2,6xp3,6xp4,6xp5,6xp6,6xp7,6xp8,6xpa,6zqz,6t5h,6xwd,6xxc')
@click.option('--omit', type=int, default=0)
@click.option('--output', type=str, default='')
@click.pass_context
def sifts_mapping(ctx, input, column, sep, func, kwargs, chunksize, entry_filter, chain_filter, skip_pdbs, omit, output):
    def get_unp_id(args):
        Entry, isoform, is_canonical = args
        return Entry if is_canonical else isoform

    SIFTS.entry_filter = entry_filter
    SIFTS.chain_filter = chain_filter
    sqlite_api = ctx.obj['custom_db']
    output = f'{func}.tsv' if output == '' else output
    output_path = ctx.obj['folder']/output
    skip_pdbs = skip_pdbs.split(',')
    kwargs = eval(kwargs)
    
    if input is None:
        total = unsync_run(sqlite_api.database.fetch_one(
            query="SELECT COUNT(DISTINCT isoform) FROM IDMapping WHERE isoform != 'NaN'"))[0]
        click.echo(f"Total {total} to query")
        for i in range(ceil(total/chunksize)):
            res = unsync_run(sqlite_api.database.fetch_all(
                query=f"""
                SELECT DISTINCT Entry,isoform,is_canonical FROM IDMapping
                WHERE isoform != 'NaN'
                LIMIT {chunksize} OFFSET {omit+chunksize*i}
                """))
            res = SIFTSs(map(get_unp_id, res)).fetch(func, skip_pdbs=skip_pdbs, **kwargs).run(tqdm).result()
            for dfrm in res:
                if dfrm is None:
                    continue
                dfrm[sorted(dfrm.columns)].to_csv(output_path, sep='\t', index=False,
                            header=not output_path.exists(), mode='a+')
            click.echo(f'Done: {len(res)+chunksize*i}')
            if len(res) < chunksize:
                break
            sleep(uniform(1, 10))
    else:
        if column is None:
            ids = read_csv(input, sep=sep, header=None, skiprows=omit if omit > 0 else None)[0].unique()
        else:
            ids = read_csv(input, sep=sep, usecols=[column], skiprows=omit if omit > 0 else None)[column].unique()
        total = len(ids)
        click.echo(f"Total {total} to query")
        for i in range(0, total, chunksize):
            res = SIFTSs(ids[i:i+chunksize]).fetch(func, skip_pdbs=skip_pdbs, **kwargs).run(tqdm).result()
            for dfrm in res:
                if dfrm is None:
                    continue
                elif isinstance(dfrm, DataFrame):
                    dfrm[sorted(dfrm.columns)].to_csv(output_path, sep='\t', index=False, header=not output_path.exists(), mode='a+')
                else:
                    pass
            click.echo(f'Done: {i+len(res)}')
            if len(res) < chunksize:
                break
            sleep(uniform(1, 10))


'''
('5jm5', '6vnn', '2i6l', '4zai', '5jn1', '6bj0', '6yth', '6wrg') + 
('4fc3', '7acu', '6lsd', '6llc', '6xoz', '6xp0', '6xp1', '6xp2', '6xp3', 
 '6xp4', '6xp5', '6xp6', '6xp7', '6xp8', '6xpa', '6zqz', '6t5h', '6xwd', 
 '6xxc')
'''


@Interface.command("residue-mapping")
@click.option('--input', type=click.Path())
@click.option('--chunksize', type=int, help="the chunksize parameter", default=10000)
@click.option('--output', type=str)
def residue_mapping(input, chunksize, output):
    output = Path(output)
    dfs = read_csv(input, sep='\t', keep_default_na=False,
                   na_values=['NULL', 'null'], chunksize=chunksize)
    for df in dfs:
        for col in ('new_pdb_range', 'new_unp_range', 'conflict_pdb_index'):
            df[col] = df[col].apply(eval)
        ob = PDBs(())
        ob.tasks = [PDB(row.pdb_id).get_expanded_map_res_df(
                    row.UniProt,
                    row.new_unp_range,
                    row.new_pdb_range,
                    conflict_pdb_index=row.conflict_pdb_index,
                    struct_asym_id=row.struct_asym_id) for _, row in df.iterrows()]
        res = ob.run(tqdm).result()
        res_mapping_df = concat(res, sort=False, ignore_index=True)
        res_mapping_df[sorted(res_mapping_df.columns)].to_csv(output, sep='\t', mode='a+', index=False, header=not output.exists())
        sleep(uniform(0, 1))


if __name__ == '__main__':
    Interface(obj={})

# @Created Date: 2020-11-23 10:29:17 am
# @Filename: command.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2022-09-03 05:09:14 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pdb_profiling import default_config
from pdb_profiling.commands import CustomDB
from pdb_profiling.processors import Identifiers, Identifier, SIFTSs, SIFTS, PDB, PDBs
from pdb_profiling.utils import unsync_run, aa_three2one, a_concat
from pdb_profiling.fetcher.webfetch import ensure
import click
from unsync import unsync
from pandas import read_csv, concat, DataFrame, read_sql_query
from importlib import util as imp_util
from pathlib import Path
from math import ceil
from rich.progress import Progress, BarColumn, TimeRemainingColumn
from rich.console import Console
from time import sleep as tsleep
import orjson as json
from random import uniform


console = Console()


progress_bar_args = ("[progress.percentage]{task.percentage:>3.1f}%", BarColumn(), "[magenta]{task.completed} of {task.total}", "[", TimeRemainingColumn(), "{task.elapsed:>3.2f}s", "]")


def format_info(name: str, template: str = "[green]Initializing %s"):
    return template % name


@click.group(chain=True, invoke_without_command=False)
@click.option("--folder", default="./", help="The output folder.", type=click.Path())
@click.option("--custom_db", default="custom.db", type=str)
@click.option('--dropall/--no-dropall', help="whether to use existing custom DB", default=False, is_flag=True)
@click.option('--initaa/--no-initaa', default=True, is_flag=True)
@click.pass_context
def Interface(ctx, folder, custom_db, dropall, initaa):
    folder = Path(folder)
    console.log(format_info(f"Folder: {folder.absolute()}"))
    ctx.ensure_object(dict)
    ctx.obj['folder'] = folder
    default_config(folder)
    ctx.obj['custom_db'] = CustomDB("sqlite:///%s" % (folder/'local_db'/custom_db), dropall)
    if initaa:
        ctx.obj['custom_db'].sync_insert(ctx.obj['custom_db'].AAThree2one, [dict(three_letter_code=three, one_letter_code=one) for three, one in aa_three2one.items()])


@Interface.command("init")
def init_folder():
    pass


@Interface.command("insert-mutation")
@click.option("-i", "--input", help="the file that contains sites info", type=click.Path())
@click.option("--sep", default="\t", help="the seperator of input file", type=str)
@click.option("--usecols", default='ftId,Ref,Pos,Alt', help="The comma-sep columns of site info", type=str)
@click.option("--headers/--no-headers", default=True, is_flag=True)
@click.option('--readchunk', type=int, help="the chunksize parameter of pandas.read_csv", default=100000)
@click.option('--nrows', type=int, help="the nrows parameter of pandas.read_csv", default=None)
@click.option('--skiprows', type=int, help="the skiprows parameter of pandas.read_csv", default=None)
@click.option('--functionfile', help="the py file that contains custom function", default=None, type=click.Path())
@click.option('--functionname', default='do_something', type=str)
@click.pass_context
def insert_sites(ctx, input, sep, usecols, headers, readchunk, nrows, skiprows, functionfile, functionname):
    def do_nothing(dfrm):
        return dfrm.to_dict('records')

    console.log(format_info("DB: Mutation Insertion"))
    usecols = usecols.split(',')
    if functionfile is not None:
        spec = imp_util.spec_from_file_location("CustomFunc", functionfile)
        CustomFunc = imp_util.module_from_spec(spec)
        spec.loader.exec_module(CustomFunc)
        deal = getattr(CustomFunc, functionname)
        console.log(f"Success: load func: {functionname} from {functionfile}")
    else:
        deal = do_nothing
    if headers:
        df = read_csv(input, sep=sep, usecols=usecols, chunksize=readchunk,
                    nrows=nrows, skiprows=skiprows)
    else:
        df = read_csv(input, sep=sep, header=None, names=usecols, chunksize=readchunk,
                      nrows=nrows, skiprows=skiprows)
    sqlite_api = ctx.obj['custom_db']
    start = 0
    with console.status("[bold green]Trying to insert..."):
        for index, dfrm in enumerate(df):
            end = readchunk*(index+1)
            console.log(f"{start}-{end}")
            start = end+1
            sqlite_api.sync_insert(sqlite_api.Mutation, deal(dfrm))


@Interface.command("id-mapping")
@click.option('-i', '--input', type=click.Path(), default=None)
@click.option('--column', type=str, default=None)
@click.option('--sep', type=str, default='\t')
@click.option('--chunksize', type=int, help="the chunksize parameter", default=50)
@click.option('--auto_assign/--no-auto_assign', default=False, is_flag=True)
@click.option('--sleep/--no-sleep', default=True, is_flag=True)
@click.pass_context
def id_mapping(ctx, input, column, sep, chunksize, auto_assign, sleep):
    sqlite_api = ctx.obj['custom_db']
    cols = ('ftId', 'Entry', 'isoform', 'is_canonical')
    Identifier.auto_assign_when_seq_conflict = auto_assign
    if input is None:
        total = unsync_run(sqlite_api.database.fetch_one(
            query="SELECT COUNT(DISTINCT ftId) FROM Mutation WHERE ftId NOT IN (SELECT DISTINCT ftId FROM IDMapping)"))[0]
        console.log(f"Total {total} to query")
        query = f"""
                SELECT DISTINCT ftId FROM Mutation
                WHERE ftId NOT IN (SELECT DISTINCT ftId FROM IDMapping)
                LIMIT {chunksize}
                """
        for index in range(ceil(total/chunksize)):
            res = unsync_run(sqlite_api.database.fetch_all(query=query))
            if len(res) == 0:
                break
            with Progress(*progress_bar_args) as p:
                res = Identifiers(i[0] for i in res).fetch('map2unp').run(p.track).result()
            values = [dict(zip(cols, i)) for i in res]
            if values:
                sqlite_api.sync_insert(sqlite_api.IDMapping, values)
            console.log(f'Done: {len(res)+chunksize*index}')
            if sleep:
                tsleep(uniform(1, 10))
    else:
        if column is None:
            ids = read_csv(input, sep=sep, header=None)[0].unique()
        else:
            ids = read_csv(input, sep=sep, usecols=[column])[column].unique()
        total = len(ids)
        console.log(f"Total {total} to query")
        for index in range(0, total, chunksize):
            with Progress(*progress_bar_args) as p:
                res = Identifiers(ids[index:index+chunksize]).fetch('map2unp').run(p.track).result()
            values = [dict(zip(cols, i)) for i in res]
            if values:
                sqlite_api.sync_insert(sqlite_api.IDMapping, values)
            console.log(f'Done: {len(res)+index}')
            if sleep:
                tsleep(uniform(1, 10))


@Interface.command('check-muta-conflict')
@click.option('--chunksize', type=int, default=100000)
@click.pass_context
def check_muta_conflict(ctx, chunksize):

    def get_seq(seq_dict, iso, pos):
        try:
            return seq_dict[iso][pos-1]
        except IndexError:
            return 'X'
        except TypeError as e:
            print(str(iso), str({key for key,val in seq_dict.items() if val is None}))
            raise e

    custom_db = ctx.obj['custom_db']
    root_query = """SELECT DISTINCT isoform, Pos FROM IDMapping, Mutation
        WHERE IDMapping.ftId = Mutation.ftId AND isoform != 'NaN'"""
    fetch_iso_seq_query = "SELECT isoform, sequence FROM ALTERNATIVE_PRODUCTS WHERE isoform IN ({}) AND (sequenceStatus = 'displayed' OR sequenceStatus = 'described');"
    fetch_can_seq_query = "SELECT accession, sequence FROM INFO WHERE accession IN ({}) ;"
    total = unsync_run(custom_db.database.fetch_val(query=f"SELECT COUNT(*) FROM ({root_query});"))
    console.log(f"Total {total} to query")
    with console.status("[bold green]checking..."):
        for i in range(ceil(total/chunksize)):
            unp_pos = DataFrame(unsync_run(custom_db.database.fetch_all(query=f"{root_query} LIMIT {chunksize} OFFSET {chunksize*i};")), columns=['isoform', 'Pos'])
            mask = unp_pos.isoform.str.contains('-')
            seq_dict = dict(unsync_run(Identifier.sqlite_api.database.fetch_all(
                                query=fetch_iso_seq_query.format(','.join(f"'{ix}'" for ix in set(unp_pos[mask].isoform))
                                ))) +
                            unsync_run(Identifier.sqlite_api.database.fetch_all(
                                query=fetch_can_seq_query.format(','.join(f"'{ix}'" for ix in set(unp_pos[~mask].isoform)))
                                )))
            unp_pos['Ref'] = [get_seq(seq_dict, iso, pos) for iso, pos in zip(unp_pos.isoform, unp_pos.Pos)]
            custom_db.sync_insert(custom_db.UniProtSeq, unp_pos.to_dict('records'))
            console.log(f'Done: {len(unp_pos)+chunksize*i}')


@Interface.command("sifts-mapping")
@click.option('-i', '--input', type=click.Path(), default=None)
@click.option('--column', type=str, default=None)
@click.option('--sep', type=str, default='\t')
@click.option('--func', type=str, default='pipe_select_mo')
@click.option('-k', '--kwargs', multiple=True, type=str)
@click.option('--chunksize', type=int, help="the chunksize parameter", default=50)
@click.option('--entry_filter', type=str, default=None)
@click.option('--chain_filter', type=str, default=None)
@click.option('--skip_pdbs', multiple=True, type=str)
@click.option('--skip_carbohydrate_polymer/--no-skip_carbohydrate_polymer', default=False, is_flag=True)
@click.option('--omit', type=int, default=0)
@click.option('-o', '--output', type=str, default='')
@click.option('--iteroutput/--no-iteroutput', default=True, is_flag=True)
@click.option('--sleep/--no-sleep', default=True, is_flag=True)
@click.option('--autotype', type=str, default='from_IDMapping')
@click.pass_context
def sifts_mapping(ctx, input, column, sep, func, kwargs, chunksize, entry_filter, chain_filter, skip_pdbs, skip_carbohydrate_polymer, omit, output, iteroutput, sleep, autotype):
    def get_unp_id(args):
        Entry, isoform, is_canonical = args
        return Entry if is_canonical else isoform

    kwargs = dict(sub.split('=') for item in kwargs for sub in item.split(';'))
    if len(kwargs) > 0:
        for key,value in kwargs.items():
            kwargs[key] = eval(value)
        console.log(f"take args: {kwargs}")
    kwargs['skip_carbohydrate_polymer'] = skip_carbohydrate_polymer
    skip_pdbs = [pdbi for item in skip_pdbs for pdbi in item.split(',')]
    if skip_pdbs:
        kwargs['skip_pdbs'] = skip_pdbs

    if entry_filter is not None:
        SIFTS.entry_filter = entry_filter
    if chain_filter is not None:
        SIFTS.chain_filter = chain_filter
    sqlite_api = ctx.obj['custom_db']
    output = f'{func}.tsv' if output == '' else output
    output_path = ctx.obj['folder']/output
    
    if input is None:
        if autotype == 'from_IDMapping':
            total = unsync_run(sqlite_api.database.fetch_one(
                query="SELECT COUNT(DISTINCT isoform) FROM IDMapping WHERE isoform != 'NaN'"))[0] - omit
            console.log(f"Total {total} to query")
            for i in range(ceil(total/chunksize)):
                res = unsync_run(sqlite_api.database.fetch_all(
                    query=f"""
                    SELECT DISTINCT Entry,isoform,is_canonical FROM IDMapping
                    WHERE isoform != 'NaN'
                    LIMIT {chunksize} OFFSET {omit+chunksize*i}
                    """))
                with Progress(*progress_bar_args) as p:
                    res = SIFTSs(map(get_unp_id, res)).fetch(func, **kwargs).run(p.track).result()
                for dfrm in res:
                    if dfrm is None:
                        continue
                    dfrm[sorted(dfrm.columns)].to_csv(output_path, sep='\t', index=False,
                                header=not output_path.exists(), mode='a+')
                console.log(f'Done: {len(res)+chunksize*i}')
                #if len(res) < chunksize:
                #    break
                if sleep and len(res) == chunksize:
                    tsleep(uniform(1, 10))
        elif autotype.startswith('from_PDB'):
            cur_db_table = autotype[5:]
            total = unsync_run(sqlite_api.database.fetch_one(query=f"SELECT COUNT(DISTINCT pdb_id) FROM {cur_db_table}"))[0] - omit
            console.log(f"Total {total} to query")
            for i in range(ceil(total/chunksize)):
                res = unsync_run(sqlite_api.database.fetch_all(query=f"SELECT DISTINCT pdb_id FROM {cur_db_table} LIMIT {chunksize} OFFSET {omit+chunksize*i}"))
                with Progress(*progress_bar_args) as p:
                    res = SIFTSs(pdbi[0] for pdbi in res).fetch(func, **kwargs).run(p.track).result()
                for dfrm in res:
                    if dfrm is None:
                        continue
                    dfrm[sorted(dfrm.columns)].to_csv(output_path, sep='\t', index=False,
                                                      header=not output_path.exists(), mode='a+')
                console.log(f'Done: {len(res)+chunksize*i}')
        else:
            console.log('Unknown autotype')
            return
    else:
        if column is None:
            ids = read_csv(input, sep=sep, header=None, skiprows=omit if omit > 0 else None)[0].unique()
        else:
            ids = read_csv(input, sep=sep, usecols=[column], skiprows=omit if omit > 0 else None)[column].unique()
        total = len(ids)
        console.log(f"Total {total} to query")
        for i in range(0, total, chunksize):
            with Progress(*progress_bar_args) as p:
                res = SIFTSs(ids[i:i+chunksize]).fetch(func, **kwargs).run(p.track).result()
            if iteroutput:
                for dfrm in res:
                    if dfrm is None:
                        continue
                    elif isinstance(dfrm, DataFrame):
                        dfrm[sorted(dfrm.columns)].to_csv(output_path, sep='\t', index=False, header=not output_path.exists(), mode='a+')
                    else:
                        pass
            else:
                DataFrame(res).to_csv(output_path, sep='\t', index=False, header=False, mode='a+')
            console.log(f'Done: {i+len(res)}')
            #if len(res) < chunksize:
            #    break
            if sleep and len(res) == chunksize:
                tsleep(uniform(1, 10))


'''
('5jm5', '6vnn', '2i6l', '4zai', '5jn1', '6bj0', '6yth', '6wrg') + 
('4fc3', '7acu', '6lsd', '6llc', '6xoz', '6xp0', '6xp1', '6xp2', '6xp3', 
 '6xp4', '6xp5', '6xp6', '6xp7', '6xp8', '6xpa', '6zqz', '6t5h', '6xwd', 
 '6xxc', '1fc2')
'''


@Interface.command("residue-mapping")
@click.option('-i', '--input', type=click.Path())
@click.option('--chunksize', type=int, help="the chunksize parameter", default=500)
@click.option('-o', '--output', type=str, default=None)
@click.option('--sleep/--no-sleep', default=True, is_flag=True)
@click.pass_context
def residue_mapping(ctx, input, chunksize, output, sleep):
    dfs = read_csv(input, sep='\t', keep_default_na=False,
                   na_values=['NULL', 'null', ''], chunksize=chunksize)
    sqlite_api = ctx.obj['custom_db']
    if output is not None:
        output = ctx.obj['folder']/output
    done = 0
    for df in dfs:
        if 'new_pdb_range_raw' in df.columns:
            cur_columns = ('new_pdb_range_raw', 'new_unp_range_raw', 'conflict_pdb_index')
        else:
            cur_columns = ('new_pdb_range', 'new_unp_range', 'conflict_pdb_index')
        for col in cur_columns:
            df[col] = df[col].apply(eval)
        ob = PDBs(())
        ob.tasks = [PDB(row.pdb_id).get_ranged_map_res_df(
                    row.UniProt,
                    getattr(row, cur_columns[1]),
                    getattr(row, cur_columns[0]),
                    conflict_pdb_index=row.conflict_pdb_index,
                    struct_asym_id=row.struct_asym_id) for row in df.to_records()]
        with Progress(*progress_bar_args) as p:
            res = ob.run(p.track).result()
        res_mapping_df = concat(res, sort=False, ignore_index=True)
        done += len(df)
        if output is not None:
            res_mapping_df[sorted(res_mapping_df.columns)].to_csv(output, sep='\t', mode='a+', index=False, header=not output.exists())
        else:
            sqlite_api.sync_insert(sqlite_api.ResidueMappingRange, res_mapping_df.to_dict('records'))
        console.log(f'Done: {done}')
        if sleep:
            tsleep(uniform(0, 2))


@Interface.command('insert-sele-mapping')
@click.option('-i', '--input', type=click.Path())
@click.option('--chunksize', type=int, help="the chunksize parameter", default=10000)
@click.pass_context
def sele_mapping(ctx, input, chunksize):
    usecols = ['UniProt', 'pdb_id', 'entity_id', 'struct_asym_id', 'chain_id', 'bs_score', 'select_rank', 'select_tag', 'after_select_rank']
    dfs = read_csv(input, sep='\t', keep_default_na=False,
                   na_values=[''], 
                   chunksize=chunksize,
                   usecols=usecols)
    custom_db = ctx.obj['custom_db']
    done = 0
    with console.status("[bold green]inserting..."):
        for df in dfs:
            custom_db.sync_insert(
                custom_db.SelectedMappingMeta, 
                df.to_dict('records'))
            done += df.shape[0]
    console.log(f'Done: {done}')


@Interface.command('insert-sifts-meta')
@click.option('-i', '--input', type=click.Path())
@click.option('--chunksize', type=int, help="the chunksize parameter", default=500)
@click.option('--func', type=str, default='fetch_from_pdbe_api')
@click.option('--api_suffix', type=str)
@click.option('--then_func', type=str, default='meta_sifts_annotation')
@click.option('--sleep/--no-sleep', default=True, is_flag=True)
@click.pass_context
def insert_sifts_meta(ctx, input, chunksize, func, api_suffix, then_func, sleep):
    custom_db = ctx.obj['custom_db']
    
    @unsync
    async def insert_meta(pdb):
        df = await getattr(pdb, func)(api_suffix).then(getattr(SIFTS, then_func))
        if df is not None:
            await custom_db.async_insert(custom_db.ResidueAnnotation, df.to_dict('records'))
    
    df = read_csv(input, header=None, chunksize=chunksize, keep_default_na=False, na_values=[''])
    done = 0
    for ids in df:
        pdbs = PDBs(ids[0].unique())
        with Progress(*progress_bar_args) as p:
            pdbs.fetch(insert_meta).run(p.track).result()
        done += len(pdbs)
        console.log(f'Done: {done}')
        if sleep:
            tsleep(uniform(0, 3))


@Interface.command('insert-isoform-range')
@click.option('--chunksize', type=int, help="the chunksize parameter", default=500)
@click.pass_context
def insert_iso_range(ctx, chunksize):
    def expand_iso_range(res):
        for UniProt, iso_range in res:
            iso_range = json.loads(iso_range)
            for start, end in iso_range:
                yield dict(UniProt=UniProt, unp_start=start, unp_end=end, resource='iso_range', resource_id=str(start))
    
    custom_db = ctx.obj['custom_db']
    proteins_db = Identifier.sqlite_api
    total = unsync_run(proteins_db.database.fetch_one(query="SELECT COUNT(*) FROM ALTERNATIVE_PRODUCTS WHERE sequenceStatus='described' AND iso_range != 'NaN'"))[0]
    console.log(f"Total {total} to query")
    for i in range(ceil(total/chunksize)):
        res = unsync_run(proteins_db.database.fetch_all(
            query=f"""
            SELECT isoform, iso_range FROM ALTERNATIVE_PRODUCTS
                WHERE sequenceStatus = 'described' AND iso_range != 'NaN'
            LIMIT {chunksize} OFFSET {chunksize*i}
            """))
        custom_db.sync_insert(custom_db.UniProtAnnotation, tuple(expand_iso_range(res)))
        console.log(f'Done: {len(res)+chunksize*i}')


def pi2records(dfrm: DataFrame, usecols: list, pair_cols: list):
    yield from yield_interact_records(dfrm[usecols[:13]].rename(columns=dict(zip(usecols[6:13], pair_cols))))
    yield from yield_interact_records(dfrm[usecols[:6]+usecols[13:]].rename(columns=dict(zip(usecols[13:], pair_cols))))

def yield_interact_records(dfrm: DataFrame):
    if 'UniProt' in dfrm.columns:
        for row in dfrm.itertuples(index=False):
            for beg, end in eval(row.interface_range):
                yield dict(UniProt=row.UniProt, pdb_id=row.pdb_id, entity_id=row.entity_id, 
                           struct_asym_id=row.struct_asym_id, chain_id=row.chain_id, 
                           assembly_id=row.assembly_id, model_id=row.model_id,
                           struct_asym_id_in_assembly=row.struct_asym_id_in_assembly,
                           interface_id=row.interface_id, css=row.css,
                           i_select_tag=row.i_select_tag, i_select_rank=row.i_select_rank,
                           pdb_beg=beg, pdb_end=end)
    else:
        for row in dfrm.itertuples(index=False):
            for beg, end in eval(row.interface_range):
                yield dict(UniProt='NaN', pdb_id=row.pdb_id, entity_id=row.entity_id,
                           struct_asym_id=row.struct_asym_id, chain_id=row.chain_id,
                           assembly_id=row.assembly_id, model_id=row.model_id,
                           struct_asym_id_in_assembly=row.struct_asym_id_in_assembly,
                           interface_id=row.interface_id, css=row.css,
                           i_select_tag=row.i_select_tag, i_select_rank=row.i_select_rank,
                           pdb_beg=beg, pdb_end=end)

@Interface.command('insert-interaction')
@click.option('-i', '--input', type=click.Path())
@click.option('--chunksize', type=int, help="the chunksize parameter", default=5000)
@click.option('--ppi/--no-ppi', is_flag=True, default=True)
@click.pass_context
def insert_interaction(ctx, input, chunksize, ppi):
    custom_db = ctx.obj['custom_db']
    common_cols = ['pdb_id', 'assembly_id', 'interface_id', 'css', 'i_select_tag', 'i_select_rank']
    pair_cols = ['entity_id', 'struct_asym_id', 'chain_id', 'model_id', 'struct_asym_id_in_assembly', 'interface_range', 'UniProt']
    usecols = common_cols + [col+'_1' for col in pair_cols] + [col+'_2' for col in pair_cols]
    df_usecols = usecols if ppi else usecols[:-1]
    dfs = read_csv(input, sep='\t', keep_default_na=False, na_values=[''], chunksize=chunksize, usecols=df_usecols)
    done: int = 0
    with console.status("[bold green]inserting..."):
        for df in dfs:
            custom_db.sync_insert(custom_db.PI, list(pi2records(df[df.i_select_rank.ne(-1)], usecols, pair_cols)))
            done += df.shape[0]
    console.log(f'Done: {done}')


@Interface.command('export-mutation-mapping')
@click.option('--with_id/--no-with_id', is_flag=True, default=False)
@click.option('--sele/--no-sele', is_flag=True, default=True)
@click.option('-o', '--output', type=str, help='filename of output file')
@click.pass_context
def export_residue_remapping(ctx, with_id, sele, output):
    output_path = ctx.obj['folder']/output
    query = """
        SELECT DISTINCT 
                    %s
                    CASE IDMapping.is_canonical
                        WHEN 1
                        THEN IDMapping.Entry
                        ELSE IDMapping.isoform
                    END edUniProt, Mutation.Ref, Mutation.Pos, Mutation.Alt,
                    Mutation.Pos - ResidueMappingRange.unp_beg + ResidueMappingRange.pdb_beg AS residue_number,
                    Mutation.Pos - ResidueMappingRange.unp_beg + ResidueMappingRange.auth_pdb_beg AS author_residue_number,
                    ResidueMappingRange.author_insertion_code,
                    ResidueMappingRange.observed_ratio,
                    ResidueMappingRange.pdb_id,
                    ResidueMappingRange.entity_id,
                    ResidueMappingRange.struct_asym_id,
                    ResidueMappingRange.chain_id,
                    {}
        FROM Mutation,ResidueMappingRange
            INNER JOIN IDMapping ON Mutation.ftId = IDMapping.ftId
            INNER JOIN UniProtSeq ON UniProtSeq.isoform = IDMapping.isoform 
                                AND UniProtSeq.Pos = Mutation.Pos 
                                AND UniProtSeq.Ref = Mutation.Ref
            INNER JOIN SelectedMappingMeta ON SelectedMappingMeta.UniProt = ResidueMappingRange.UniProt
                                        AND SelectedMappingMeta.pdb_id = ResidueMappingRange.pdb_id
                                        AND SelectedMappingMeta.entity_id = ResidueMappingRange.entity_id
                                        AND SelectedMappingMeta.struct_asym_id = ResidueMappingRange.struct_asym_id
        WHERE ResidueMappingRange.UniProt = edUniProt
        AND Mutation.Pos >= ResidueMappingRange.unp_beg
        AND Mutation.Pos <= ResidueMappingRange.unp_end
        AND ResidueMappingRange.conflict_code IS NULL
        AND ResidueMappingRange.observed_ratio > 0
        AND (ResidueMappingRange.residue_name = '' OR ResidueMappingRange.residue_name IN (SELECT three_letter_code FROM AAThree2one))
        AND SelectedMappingMeta.select_rank != -1
        {} ;"""
    if with_id:
        query = query % 'Mutation.ftId,'
        if sele:
            query = query.format('MIN(SelectedMappingMeta.after_select_rank)', 'GROUP BY Mutation.ftId, ResidueMappingRange.UniProt, Mutation.Pos, Mutation.Alt')
        else:
            query = query.format('SelectedMappingMeta.after_select_rank', '')
    else:
        query = query % ''
        if sele:
            query = query.format('MIN(SelectedMappingMeta.after_select_rank)', 'GROUP BY ResidueMappingRange.UniProt, Mutation.Pos, Mutation.Alt')
        else:
            query = query.format('SelectedMappingMeta.after_select_rank', '')
    with console.status("[bold green]query..."):
        dfs = read_sql_query(query, ctx.obj['custom_db'].engine, chunksize=10000)
        for df in dfs:
            if df.shape[0] == 0:
                continue
            df.rename(columns={'edUniProt': 'UniProt'}).to_csv(
                output_path, index=False, mode='a+', sep='\t', header=not output_path.exists())
    console.log(f'result saved in {output_path}')


@Interface.command('export-pdb-mutation-mapping')
@click.option('--auth/--no-auth', is_flag=True, default=True)
@click.option('-o', '--output', type=str, help='filename of output file')
@click.pass_context
def export_residue_remapping(ctx, auth, output):
    output_path = ctx.obj['folder']/output
    if auth:
        query = """
        SELECT DISTINCT ResidueMappingRange.UniProt,
                        PDBAuthMutation.pdb_id,
                        ResidueMappingRange.entity_id,
                        ResidueMappingRange.struct_asym_id,
                        ResidueMappingRange.chain_id,
                        PDBAuthMutation.Ref,
                        PDBAuthMutation.Alt,
                        PDBAuthMutation.author_residue_number,
                        PDBAuthMutation.author_insertion_code,
                        ResidueMappingRange.observed_ratio,
                        ResidueMappingRange.conflict_code,
                        PDBAuthMutation.author_residue_number - ResidueMappingRange.auth_pdb_beg + ResidueMappingRange.pdb_beg AS residue_number,
                        PDBAuthMutation.author_residue_number - ResidueMappingRange.auth_pdb_beg + ResidueMappingRange.unp_beg AS unp_residue_number
        FROM PDBAuthMutation,ResidueMappingRange
        WHERE PDBAuthMutation.pdb_id == ResidueMappingRange.pdb_id
        AND PDBAuthMutation.author_residue_number >= ResidueMappingRange.auth_pdb_beg
        AND PDBAuthMutation.author_residue_number <= ResidueMappingRange.auth_pdb_end
        AND PDBAuthMutation.author_insertion_code == ResidueMappingRange.author_insertion_code
        AND 
            CASE 
                WHEN PDBAuthMutation.struct_asym_id IS NOT NULL
                THEN PDBAuthMutation.struct_asym_id == ResidueMappingRange.struct_asym_id
                WHEN PDBAuthMutation.chain_id IS NOT NULL
                THEN PDBAuthMutation.chain_id == ResidueMappingRange.chain_id
            ELSE
                1
            END
        ;
        """
    else:
        query = """
        SELECT DISTINCT ResidueMappingRange.UniProt,
                        PDBMutation.pdb_id,
                        ResidueMappingRange.entity_id,
                        ResidueMappingRange.struct_asym_id,
                        ResidueMappingRange.chain_id,
                        PDBMutation.Ref,
                        PDBMutation.Alt,
                        PDBMutation.residue_number,
                        ResidueMappingRange.observed_ratio,
                        ResidueMappingRange.conflict_code,
                        PDBMutation.residue_number - ResidueMappingRange.pdb_beg + ResidueMappingRange.auth_pdb_beg AS author_residue_number,
                        ResidueMappingRange.author_insertion_code,
                        PDBMutation.residue_number - ResidueMappingRange.pdb_beg + ResidueMappingRange.unp_beg AS unp_residue_number
        FROM PDBMutation,ResidueMappingRange
        WHERE PDBMutation.pdb_id == ResidueMappingRange.pdb_id
        AND PDBMutation.residue_number >= ResidueMappingRange.pdb_beg
        AND PDBMutation.residue_number <= ResidueMappingRange.pdb_end
        AND 
            CASE 
                WHEN PDBMutation.struct_asym_id IS NOT NULL
                THEN PDBMutation.struct_asym_id == ResidueMappingRange.struct_asym_id
                WHEN PDBMutation.chain_id IS NOT NULL
                THEN PDBMutation.chain_id == ResidueMappingRange.chain_id
            ELSE
                1
            END
        ;
        """
    with console.status("[bold green]query..."):
        dfs = read_sql_query(query, ctx.obj['custom_db'].engine, chunksize=10000)
        for df in dfs:
            if df.shape[0] == 0:
                continue
            reslist = PDBs(frozenset(df.pdb_id)).fetch(
                'fetch_from_pdbe_api',
                api_suffix='api/pdb/entry/residue_listing/',
                then_func=PDB.to_dataframe).run().then(a_concat).result()
            newdf = df.merge(reslist, how='left')
            newdf['pdb_one_letter_code'] = newdf.residue_name.map(aa_three2one)
            newdf['unp_one_letter_code'] = newdf.apply(lambda x: x['conflict_code'] if isinstance(x['conflict_code'], str) and x['conflict_code'] != '' else x['pdb_one_letter_code'], axis=1)
            newdf.drop(columns=['conflict_code']).to_csv(
                output_path, index=False, mode='a+', sep='\t', header=not output_path.exists())
    console.log(f'result saved in {output_path}')

@Interface.command('export-interaction-mapping')
@click.option('--with_id/--no-with_id', is_flag=True, default=False)
@click.option('-o', '--output', type=str, help='filename of output file')
@click.pass_context
def export_interface_mapping(ctx, with_id, output):
    output_path = ctx.obj['folder']/output
    query = """
        SELECT DISTINCT 
                    %s
                    CASE IDMapping.is_canonical
                        WHEN 1
                        THEN IDMapping.Entry
                        ELSE IDMapping.isoform
                    END edUniProt, Mutation.Ref, Mutation.Pos, Mutation.Alt,
                    Mutation.Pos - ResidueMappingRange.unp_beg + ResidueMappingRange.pdb_beg AS residue_number,
                    Mutation.Pos - ResidueMappingRange.unp_beg + ResidueMappingRange.auth_pdb_beg AS author_residue_number,
                    ResidueMappingRange.author_insertion_code,
                    ResidueMappingRange.observed_ratio,
                    ResidueMappingRange.pdb_id,
                    ResidueMappingRange.entity_id,
                    ResidueMappingRange.struct_asym_id,
                    ResidueMappingRange.chain_id,
                    PI.assembly_id,
                    PI.model_id,
                    PI.struct_asym_id_in_assembly,
                    PI.interface_id,
                    PI.css,
                    PI.i_select_tag,
                    PI.i_select_rank,
                    (Mutation.Pos - ResidueMappingRange.unp_beg + ResidueMappingRange.pdb_beg >= PI.pdb_beg AND Mutation.Pos - ResidueMappingRange.unp_beg + ResidueMappingRange.pdb_beg <= PI.pdb_end) AS is_interface_residue
        FROM Mutation,ResidueMappingRange
            INNER JOIN IDMapping ON Mutation.ftId = IDMapping.ftId
            INNER JOIN UniProtSeq ON UniProtSeq.isoform = IDMapping.isoform 
                                AND UniProtSeq.Pos = Mutation.Pos 
                                AND UniProtSeq.Ref = Mutation.Ref
            INNER JOIN SelectedMappingMeta ON SelectedMappingMeta.UniProt = ResidueMappingRange.UniProt
                                        AND SelectedMappingMeta.pdb_id = ResidueMappingRange.pdb_id
                                        AND SelectedMappingMeta.struct_asym_id = ResidueMappingRange.struct_asym_id
            INNER JOIN PI ON PI.UniProt = ResidueMappingRange.UniProt
                         AND PI.pdb_id = ResidueMappingRange.pdb_id
                         AND PI.struct_asym_id = ResidueMappingRange.struct_asym_id
        WHERE ResidueMappingRange.UniProt = edUniProt
        AND Mutation.Pos >= ResidueMappingRange.unp_beg
        AND Mutation.Pos <= ResidueMappingRange.unp_end
        AND ResidueMappingRange.conflict_code IS NULL
        AND ResidueMappingRange.observed_ratio > 0
        AND (ResidueMappingRange.residue_name = '' OR ResidueMappingRange.residue_name IN (SELECT three_letter_code FROM AAThree2one))
        AND SelectedMappingMeta.select_rank != -1
        ;
    """
    query = query % ('Mutation.ftId,' if with_id else '')
    with console.status("[bold green]query..."):
        dfs = read_sql_query(query, ctx.obj['custom_db'].engine, chunksize=10000)
        for df in dfs:
            if df.shape[0] == 0:
                continue
            #df.sort()
            #df.drop_duplicates(subset=df.columns[:-1], inplace=True)
            df.rename(columns={'edUniProt': 'UniProt'}).to_csv(
                output_path, index=False, mode='a+', sep='\t', header=not output_path.exists())
    console.log(f'result saved in {output_path}')


@Interface.command('insert-sele-mutation-mapping')
@click.option('-i', '--input', type=click.Path())
@click.option('--chunksize', type=int, help="the chunksize parameter", default=10000)
@click.pass_context
def insert_mapped_resmap(ctx, input, chunksize):
    custom_db = ctx.obj['custom_db']
    dfs = read_csv(input, sep='\t', keep_default_na=False,
                   na_values=[''], chunksize=chunksize,
                   usecols=['UniProt', 'Ref', 'Pos', 'Alt'])
    done = 0
    for df in dfs:
        custom_db.sync_insert(custom_db.MappedMutation, df.to_dict('records'))
        done += df.shape[0]
        console.log(f'Done: {done}')


@Interface.command('export-smr-mutation-mapping')
@click.option('--identity_cutoff', type=float, default=0)
@click.option('--length_cutoff', type=int, default=0)
@click.option('--with_id/--no-with_id', is_flag=True, default=False)
@click.option('--sele/--no-sele', is_flag=True, default=True)
@click.option('--allow_oligo_state', type=str, default=None)
@click.option('-o', '--output', type=str, help='filename of output file')
@click.pass_context
def export_smr_residue_remapping(ctx, identity_cutoff, length_cutoff, with_id, sele, allow_oligo_state, output):
    output_path = ctx.obj['folder']/output
    # sele_o_path = ctx.obj['folder']/(output_path.name.replace(output_path.suffix,'')+'.sele'+output_path.suffix)
    query = """
    SELECT DISTINCT
        %s
        CASE IDMapping.is_canonical
                    WHEN 1
                    THEN IDMapping.Entry
                    ELSE IDMapping.isoform
        END edUniProt, Mutation.Ref, Mutation.Pos, Mutation.Alt,
        SMRModel.oligo_state,SMRModel.select_tag,SMRModel.coordinates,
        {}
    FROM Mutation, SMRModel
        INNER JOIN IDMapping ON Mutation.ftId = IDMapping.ftId
        INNER JOIN UniProtSeq ON UniProtSeq.isoform = IDMapping.isoform 
                            AND UniProtSeq.Pos = Mutation.Pos 
                            AND UniProtSeq.Ref = Mutation.Ref
    WHERE SMRModel.UniProt = edUniProt
    AND Mutation.Pos >= SMRModel.unp_beg
    AND Mutation.Pos <= SMRModel.unp_end
    AND SMRModel.identity >= %s
    AND SMRModel.select_rank > 0
    AND SMRModel.unp_end - SMRModel.unp_beg + 1 >= %s
    %s
    AND NOT EXISTS (SELECT * FROM MappedMutation 
                  WHERE edUniProt = MappedMutation.UniProt 
                    AND MappedMutation.Pos = Mutation.Pos 
                    AND MappedMutation.Alt = Mutation.Alt LIMIT 1)
    {};
    """
    if with_id:
        if allow_oligo_state is None:
            query = query % ('Mutation.ftId,', identity_cutoff, length_cutoff, '')
        else:
            query = query % ('Mutation.ftId,', identity_cutoff, length_cutoff, f"AND SMRModel.oligo_state IN {allow_oligo_state}")
    else:
        if allow_oligo_state is None:
            query = query % ('', identity_cutoff, length_cutoff, '')
        else:
            query = query % ('', identity_cutoff, length_cutoff, f"AND SMRModel.oligo_state IN {allow_oligo_state}")
    if sele:
        if with_id:
            query = query.format('MIN(SMRModel.select_rank)', 'GROUP BY Mutation.ftId, SMRModel.UniProt, Mutation.Pos, Mutation.Alt')
        else:
            query = query.format('MIN(SMRModel.select_rank)', 'GROUP BY SMRModel.UniProt, Mutation.Pos, Mutation.Alt')
    else:
        query = query.format('SMRModel.select_rank', '')
    with console.status("[bold green]query..."):
        dfs = read_sql_query(query, ctx.obj['custom_db'].engine, chunksize=10000)
        for df in dfs:
            if df.shape[0] == 0:
                continue
            df.rename(columns={'edUniProt': 'UniProt'}).to_csv(
                output_path, index=False, mode='a+', sep='\t', header=not output_path.exists())
    console.log(f'result saved in {output_path}')
    #full_df = read_csv(output_path, sep='\t', keep_default_na=False)
    #best_indexes = full_df.groupby(['UniProt','Pos', 'Alt']).select_rank.idxmin()
    #full_df.loc[best_indexes].to_csv(sele_o_path, sep='\t', index=False)
    #console.log(f'sele result saved in {sele_o_path}')


@Interface.command('insert-smr-mapping')
@click.option('-i', '--input', type=click.Path())
@click.option('--chunksize', type=int, help="the chunksize parameter", default=10000)
@click.pass_context
def insert_smr_mapping(ctx, input, chunksize):
    custom_db = ctx.obj['custom_db']
    dfs = read_csv(input, sep='\t', keep_default_na=False, 
                   na_values=[''], chunksize=chunksize, 
                   usecols=['UniProt', 'coordinates', 'from', 'to', 'identity', 'similarity', 'coverage', 'oligo-state', 'ligand_chains', 'select_rank', 'select_tag'])
    done = 0
    for df in dfs:
        df['with_ligand'] = df.ligand_chains.notnull()
        df = df.drop(columns=['ligand_chains']).rename(columns={'oligo-state': 'oligo_state', 'from': 'unp_beg', 'to': 'unp_end'})
        custom_db.sync_insert(custom_db.SMRModel, df.to_dict('records'))
        done += df.shape[0]
        console.log(f'Done: {done}')


@Interface.command('fetch1pdb')
@click.option('-i', '--pdb', type=str, help="PDB Identifier")
@click.option('-a', '--api', type=str, help="API Name")
@click.option('-p', '--params', multiple=True, type=str, default=None)
@click.option('-d', '--data_collection', multiple=True, type=str, default=None)
@click.option('-m', '--method', type=str, default='get')
@click.option('-t', '--tag', type=str, default='subset')
@click.option('--use_existing/--no-use_existing', is_flag=True, default=True)
def fetch1pdb(pdb, api, params, data_collection, method, tag, use_existing):
    if params is not None:
        params = dict(sub.split('=') for item in params for sub in item.split(','))
        if len(params)> 0:
            console.log(f"take args: {params}")
    else:
        params = {}
    if method == 'get':
        if len(data_collection) > 0:
            console.log(f"current impl '{method}' method, would omit data_collection: {str(data_collection)[:100]}!")
        data = None
    elif method == 'post':
        data = {'atom_site': []}
        for item in data_collection:
            data['atom_site'].append(dict(sub.split('=') for sub in item.split(',')))
        data = json.dumps(data)
        console.log(f"take data_collection: {data[:100]}")
    else:
        return
    pdb = PDB(pdb)
    ensure.set_use_existing(use_existing)
    with console.status("[bold green]download..."):
        res = pdb.fetch_from_modelServer_api(
            api_suffix=api, method=method, data_collection=data, filename=tag, **params).result()
    console.log(f"Result saved in {res}")


@Interface.command("insert-pdb-mutation")
@click.option("-i", "--input", help="the file that contains sites info", type=click.Path())
@click.option("--sep", default="\t", help="the seperator of input file", type=str)
@click.option("--usecols", help="The comma-sep columns of site info.", type=str)
@click.option("--headers/--no-headers", default=True, is_flag=True)
@click.option('--readchunk', type=int, help="the chunksize parameter of pandas.read_csv", default=100000)
@click.option('--nrows', type=int, help="the nrows parameter of pandas.read_csv", default=None)
@click.option('--skiprows', type=int, help="the skiprows parameter of pandas.read_csv", default=None)
@click.option("--auth/--no-auth", default=True, is_flag=True)
@click.pass_context
def insert_pdb_mutation(ctx, input, sep, usecols, headers, readchunk, nrows, skiprows, auth):
    console.log(format_info("DB: PDB Mutation Insertion"))
    usecols = usecols.split(',')
    if headers:
        df = read_csv(input, sep=sep, usecols=usecols, chunksize=readchunk,
                      nrows=nrows, skiprows=skiprows, keep_default_na=False)
    else:
        df = read_csv(input, sep=sep, header=None, names=usecols, chunksize=readchunk,
                      nrows=nrows, skiprows=skiprows, keep_default_na=False)
    sqlite_api = ctx.obj['custom_db']
    db_table = sqlite_api.PDBAuthMutation if auth else sqlite_api.PDBMutation
    start = 0
    with console.status("[bold green]Trying to insert..."):
        for index, dfrm in enumerate(df):
            end = readchunk*(index+1)
            console.log(f"{start}-{end}")
            start = end+1
            sqlite_api.sync_insert(db_table, dfrm.to_dict('records'))


if __name__ == '__main__':
    Interface(obj={})

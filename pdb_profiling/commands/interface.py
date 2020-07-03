# @Created Date: 2019-11-24 11:03:59 pm
# @Filename: interface.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-05-23 02:24:53 pm
# @Copyright (c) 2019 MinghuiGroup, Soochow University
from __future__ import absolute_import
import click
import json
from pandas import read_csv, merge, Series, DataFrame
from numpy import nan
from pathlib import Path
from typing import Iterator
from unsync import Unfuture, unsync
import asyncio
import sys
from importlib import util as imp_util
from tqdm import tqdm


if "\\" in __file__:
    # Windows
    _SEP = "\\"
    sys.path.append("C:\GitWorks\pdb-profiling")
else:
    # Linux
    _SEP = "/"
    sys.path.append("/data/zzf/2020/src/pdb-profiling")

try:
    from pdb_profiling.processers.uniprot.api import MapUniProtID, UniProtFASTA
    from pdb_profiling.processers.pdbe.neo4j_api import Neo4j_API
    from pdb_profiling.processers.pdbe.sqlite_api import Sqlite_API
    from pdb_profiling.fetcher.webfetch import UnsyncFetch
    from pdb_profiling.fetcher.dbfetch import Neo4j
    from pdb_profiling.log import Abclog
    from pdb_profiling.utils import related_dataframe
    from pdb_profiling.pipelines.score import Score_API
    from pdb_profiling.pipelines.select import Select_API
except Exception as e:
    raise e

def colorClick(name: str, template: str = "Initializing %s", fg: str = "green"):
    return click.style(template % name, fg=fg)


@unsync
async def init_semaphore(concurreq) -> Unfuture:
    """
    `semaphore` initiated in the `unsync` event loop
    """
    await asyncio.sleep(.01)
    return asyncio.Semaphore(concurreq)


@click.group(chain=True, invoke_without_command=False)
@click.option("--folder", default="", help="The file folder of new files.", type=click.Path())
@click.option("--loggingpath", default=None, help="The file path of logging.", type=click.Path())
@click.option('--useexisting/--no-useexisting', help="whether to use existing result files (For UnsyncFetch)", default=False, is_flag=True)
@click.pass_context
def Interface(ctx, folder, loggingpath, useexisting):
    if not folder:
        return
    # Init Folder Setting
    click.echo(colorClick("Folder"))
    folder = Path(folder)
    ctx.ensure_object(dict)
    ctx.obj['folder'] = folder
    ctx.obj['loggingpath'] = loggingpath
    for cur_path in (folder/'UniProt'/'mapping',
                     folder/'UniProt'/'fasta',
                     folder/'DB',
                     folder/'I3D'):
        cur_path.mkdir(parents=True, exist_ok=True)
        ctx.obj[f'{cur_path.stem}_folder'] = cur_path
    # Init Logging Setting
    click.echo(colorClick("Logger"))
    ctx.obj['logger'] = Abclog.set_logging_fileHandler(loggingpath, logName='CommandLine')
    # For TEST
    UnsyncFetch.use_existing = useexisting
  

@Interface.resultcallback()
@click.pass_context
def process_pipeline(ctx, processors, folder, loggingpath, useexisting):    
    def iter_task(task: Unfuture):
        if isinstance(task, Unfuture):
            return iter_task(task.result())
        else:
            return task
    
    iterator = ctx.obj.get('iterator', None)
    if iterator is None:
        return
    for processor in processors:
        if processor is not None:
            iterator = processor(iterator)
    tasks = list(iterator)
    for task in tqdm(tasks, total=len(tasks)):
        iter_task(task)
    # click.echo(res)
    # UnsyncFetch.unsync_tasks(list(iterator)).result()
    # [await fob for fob in tqdm(asyncio.as_completed(iterator), total=len(iterator))]
    # sleep(30)


@Interface.command("UniProt.init")
@click.option('--concurreq', type=int, help="the number of concurent requests (For UniProt API)", default=20)
@click.option('--concurrate', type=float, help="the rate of concurent requests (For UniProt API)", default=2)
@click.pass_context
def init_unp(ctx, concurreq, concurrate):
    click.echo(colorClick("UniProt API"))
    ctx.obj['unp_concurreq'] = concurreq
    ctx.obj['unp_concurrate'] = concurrate
    # Set Semaphore
    click.echo("Set Semaphore (For UniProt API)")
    ctx.obj['semaphore'] = init_semaphore(concurreq).result()


@Interface.command("UniProt.id-mapping")
@click.option("--input", 
              default="",
              help="The reference file of IDs that need to map via UniProt RESTful API.", 
              type=click.Path(exists=True))
@click.option("--sep", 
              default="\t", 
              help="The seperator of referenceFile.",
              type=str)
@click.option("--idcol",
              help="The column name of IDs in referenceFile.", 
              type=str)
@click.option("--idtype", 
              help="ID Abbreviation that stands for the type of ID. (e.g P_REFSEQ_AC)",
              type=str)
@click.option("--usecols", 
              # feature(CHAIN)
              default="id,genes,reviewed,comment(ALTERNATIVE%20PRODUCTS),feature(ALTERNATIVE%20SEQUENCE),organism,protein%20names",
              help="Comma-separated list of the column names for programmatic access to the UniProtKB search results.", 
              type=str)
@click.option("--genecol",
              default=None, 
              help="The column name of gene info in referenceFile.",
              type=str)
@click.option('--querychunk', type=int, help="the chunksize of query ids", default=50)
@click.option('--readchunk', type=int, help="the chunksize parameter of pandas.read_csv", default=None)
@click.option('--nrows', type=int, help="the nrows parameter of pandas.read_csv", default=None)
@click.option('--skiprows', type=int, help="the skiprows parameter of pandas.read_csv", default=None)
@click.option('--outname', type=str, help="the filename stem of output files", default="unp_yourlist")
@click.pass_context
def idMappingViaUnpApi(ctx, input, sep, idcol, idtype, usecols, genecol, querychunk, readchunk, nrows, skiprows, outname):
    def yieldTasks(df):
        for dfrm in df:
            demo = MapUniProtID(id_col=idcol,
                                id_type=idtype,
                                dfrm=dfrm,
                                usecols=usecols,
                                gene_col=genecol,
                                logger=ctx.obj['logger'])
            unsync_tasks = demo.retrieve(ctx.obj['mapping_folder']/outname,
                                        chunksize=querychunk,
                                        concur_req=ctx.obj['unp_concurreq'],
                                        rate=ctx.obj['unp_concurrate'],
                                        run_tasks=False,
                                        semaphore = ctx.obj['semaphore'])
            for task in unsync_tasks:
                yield task

    def processor(iterator: Iterator[Unfuture]):
        '''This processor do nothing'''
        for task in iterator:
            yield task

    click.echo(colorClick("UniProt Retrieve/ID Mapping"))
    df = read_csv(input, sep=sep, chunksize=readchunk, nrows=nrows, skiprows=skiprows, usecols=[
                  col for col in (idcol, genecol) if col is not None])
    if isinstance(df, DataFrame):
        df = (df,)
    click.echo(colorClick("Set Task Iterator", "%s"))
    ctx.obj['iterator'] = yieldTasks(df)
    return processor
    

@Interface.command("UniProt.download-seq")
@click.option("--input", default="", help="the file of UniProt IDs", type=click.Path())
@click.option("--idcol", default=None, help="the column name of UniProt IDs (If not specified, use the first col)", type=str)
@click.option("--sep", default="\t", help="the seperator of input file", type=str)
@click.option("--include", default="no", help="whether to include isoform sequences in one single fasta file", type=click.Choice(['yes', 'no']))
@click.pass_context
def downloadUnpSeq(ctx, input, idcol, sep, include):
    def processor(iterator: Iterator[Unfuture]):
        '''Continuations for Retrieving UniProt Fasta Sequence'''
        for task in iterator:
            yield task.then(UniProtFASTA.process)

    click.echo(colorClick("UniProt Sequence Download"))
    UniProtFASTA.logger = ctx.obj['logger']
    UniProtFASTA.params['include'] = include
    if input:
        path = Path(input)
        if idcol is not None:
            unps = read_csv(path, sep=sep, usecols=[idcol])[idcol].drop_duplicates()
        else:
            with path.open('rt') as infile:
                unps = [line.strip() for line in infile]
        UniProtFASTA.retrieve(unps, ctx.obj['fasta_folder'], ctx.obj['unp_concurreq'],
                              ctx.obj['unp_concurrate'], semaphore=ctx.obj['semaphore'])
        return
    else:
        UniProtFASTA.obj = {key: ctx.obj[key]
                            for key in ('fasta_folder', 'unp_concurreq', 'unp_concurrate', 'semaphore')}
        return processor


'''
TODO:

1. Init Local Database
2. Init Site Info [Insert] For residue mapping 
3. Seperate sifts mapping/filtering and residue mapping
'''

@Interface.command("DB.init")
@click.option('--db', help="the name of local sqlite database file in the folder", default="local_sqlite.db", type=str)
@click.option('--dropall/--no-dropall', help="whether to drop all the tables in the local sqlite database", default=False, is_flag=True)
@click.option('--remotedburl', default=None, help="the url of remote neo4j database", type=str)
@click.option('--remotedbuser', default=None, help="the user-name that accesses to remote neo4j database", type=str)
@click.option('--remotedbpass', default=None, help="the password that accesses to remote neo4j database", prompt=True, hide_input=True)
@click.option('--concurreq', type=int, help="the number of concurent requests", default=3)
@click.option('--insertsleep', type=float, help="the sleep duration when encounter database lock error", default=45.5)
@click.pass_context
def init_db(ctx, db, dropall, remotedburl, remotedbuser, remotedbpass, concurreq, insertsleep):
    # Init Local DataBase Setting
    click.echo(colorClick(f"Local DB (dropall: {dropall})"))
    db_path = ctx.obj['DB_folder']/db
    sqlite_api = Sqlite_API("sqlite:///%s" % str(db_path), dropall, insertsleep)
    Neo4j_API.sqlite_api = sqlite_api
    ctx.obj['sqlite_api'] = sqlite_api
    if remotedburl is not None:
        click.echo(colorClick("Remote DB (Neo4j)"))
        config = {'user': remotedbuser, 'pass': remotedbpass,
                'url': remotedburl}
        click.echo("Set Semaphore (For Neo4j API)")
        Neo4j_API.neo4j_api = Neo4j(
            config, concurreq,
            init_semaphore(concurreq).result(),
            log_func=ctx.obj['logger'].info).connnect().result()
        ctx.obj['neo4j_api'] = Neo4j_API.neo4j_api
        ctx.obj['Neo4j_API'] = Neo4j_API


@Interface.command("DB.insert-sites-info")
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
    
    click.echo(colorClick("DB Insertions"))
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
    sqlite_api = ctx.obj['sqlite_api']
    start = 0
    for index, dfrm in enumerate(df):
        end = readchunk*(index+1)
        click.echo(f"Try to insert: {start}-{end}")
        start = end+1
        sqlite_api.sync_insert(sqlite_api.Site_Info, deal(dfrm))


@Interface.command("DB.GraphDB.unp-to-pdb")
@click.option("--input", default="", help="the file of UniProt Mapping result", type=click.Path())
@click.option("--sep", default="\t", help="the seperator of input file", type=str)
@click.option('--readchunk', type=int, help="the chunksize parameter of pandas.read_csv", default=None)
@click.option('--nrows', type=int, help="the nrows parameter of pandas.read_csv", default=None)
@click.option('--skiprows', type=int, help="the skiprows parameter of pandas.read_csv", default=None)
@click.option("--unpmapfilter",
              default='{"Mapping_status":["eq","Yes"],"Organism":["eq","Homo sapiens (Human)"]}',
              help="the filter(JSON-Format Dict) of UniProt Mapping result", type=str)
@click.option("--siftsfilter",
              default='{"identity":["ge",0.9],"repeated":["eq",false]}',
              help="the filter(JSON-Format Dict) of SIFTS Mapping result", type=str)
@click.option("--entryfilter",
              default='{"METHOD_CLASS":["isin",["x-ray","nmr"]],"resolution":["le",3],"has_hybrid_nucleotides":["eq",false]}',
              help="the filter(JSON-Format Dict) of PDB Entry Info result", type=str)
@click.option("--outname", default=None, help="the output file name of filtered SIFTS Mapping result", type=str)
@click.pass_context
def neo4j_unp2pdb(ctx, input, sep, readchunk, nrows, skiprows, unpmapfilter, siftsfilter, entryfilter, outname):
    def yieldTasks(df):
        for dfrm in df:
            yield Neo4j_API.process_unp2pdb(dfrm)
    
    def processor_backup(iterator: Iterator[Unfuture]):
        '''This processor do nothing'''
        for task in iterator:
            yield task

    def processor(iterator: Iterator[Unfuture]):
        '''
        Continuations for Mapping from unp to pdb
        '''
        for task in iterator:
            yield task.then(Neo4j_API.process_unp2pdb)

    if input:
        df = read_csv(input, sep=sep, chunksize=readchunk,
                      nrows=nrows, skiprows=skiprows)
        if isinstance(df, DataFrame):
            df = (df,)
        click.echo(colorClick("Set Task Iterator", "%s"))
        ctx.obj['iterator'] = yieldTasks(df)
        ret = processor_backup
    else:
        ret = processor
    
    if ctx.obj.get('iterator', None) is not None:
        click.echo(colorClick("Set Filters", "%s"))
        click.echo(f"unpmap_filter:\n\t{unpmapfilter}")
        click.echo(f"sifts_filter:\n\t{siftsfilter}")
        click.echo(f"entry_filter:\n\t{entryfilter}")
        Neo4j_API.unpmap_filter = json.loads(unpmapfilter)
        Neo4j_API.sifts_filter = json.loads(siftsfilter)
        Neo4j_API.entry_filter = json.loads(entryfilter)
        resolution_cut_off = Neo4j_API.entry_filter.get('resolution', None)
        if resolution_cut_off is not None:
            Neo4j_API.entry_info_add['resolution'] = resolution_cut_off[1]
        if outname is not None:
            Neo4j_API.filtered_sifts_path = ctx.obj['folder']/outname
    return ret


@Interface.command("DB.GraphDB.unpres-to-pdbres")
@click.option("--outname", help="the output file name of residue-level mapping", type=click.Path())
@click.pass_context
def neo4j_res2pdb(ctx, outname):
    def processor(iterator: Iterator[Unfuture]):
        '''
        Continuations for Mapping from unp to pdb
        '''
        for task in iterator:
            yield task.then(Neo4j_API.process_mapres2pdb)
    Neo4j_API.res2pdbpath = ctx.obj['folder']/outname
    return processor


@Interface.command("DB.GraphDB.to-unpres")
@click.option("--input", type=click.Path())
@click.option("--sep", default="\t", help="the seperator of input file", type=str)
@click.option("--observedonly/--no-observedonly", default=True, is_flag=True)
@click.option('--readchunk', type=int, help="the chunksize parameter of pandas.read_csv", default=None)
@click.pass_context
def neo4j_res2unp(ctx, input, sep, observedonly, readchunk):
    respath = ctx.obj['folder']/'pdb2unp_rm.tsv'
    siftspath = ctx.obj['folder']/'pdb2unp_fs.tsv'
    def yieldTasks(df):
        for dfrm in df:
            dfrm = dfrm.drop_duplicates()
            yield Neo4j_API.process_map2unp(dfrm, respath, siftspath, observedonly)

    def processor(iterator: Iterator[Unfuture]):
        '''This processor do nothing'''
        for task in iterator:
            yield task

    df = read_csv(input, sep=sep, chunksize=readchunk)  # usecols=['pdb_id', 'entity_id', 'chain_id']
    if isinstance(df, DataFrame):
        df = (df,)
    click.echo(colorClick("Set Task Iterator", "%s"))
    ctx.obj['iterator'] = yieldTasks(df)
    return processor


@Interface.command("Stat.score-sifts")
@click.option("--input", default="", help="the file of SIFTS Mapping result", type=click.Path())
@click.option("--outname", help="the output file name of result file", type=click.Path())
@click.option("--omit", help="omit the number of unmapped range", type=int, default=5)
@click.option("--score/--no-score", help="whether to score", default=True, is_flag=True)
@click.pass_context
def score_sifts(ctx, input, outname, omit, score):
    score_api = Score_API(
        ctx.obj['sqlite_api'],
        ctx.obj['neo4j_api'],
        outpath=ctx.obj['folder']/outname,
        logger=ctx.obj['logger'],
        omit=omit,
        add_score=score)
    ctx.obj['score_api'] = score_api
    
    def processor(iterator: Iterator[Unfuture]):
        '''
        Continuations for Mapping from unp to pdb
        '''
        for task in iterator:
            yield task.then(score_api.process)
    
    if input:
        score_api.process(input).result()
        ctx.obj['logger'].info("Finish stat-sifts")
    else:
        return processor
    

@Interface.command("Stat.select-sifts")
@click.option("--input", default="", help="the file of SIFTS Mapping result", type=click.Path())
@click.option("--omitcutoff", help="the length of omitted chains", type=int, default=50)
@click.option("--omitcol", help="the column that apply omit-cutoff", type=str, default="ATOM_RECORD_COUNT")
@click.option("--oscutoff", help="the cutoff of overlap coefficient", type=float, default=0.2)
@click.option("--siftsfilter",
              default='{"identity":["ge",0.9],"repeated":["eq",false]}',
              help="the filter(JSON-Format Dict) of SIFTS Mapping result", type=str)
@click.pass_context
def oligo_sifts(ctx, input, omitcutoff, omitcol, oscutoff, siftsfilter):
    def processor(iterator: Iterator[Unfuture]):
        '''
        Continuations for Mapping from unp to pdb
        '''
        for task in iterator:
            yield task.then(ctx['select_api'].process)
    
    select_api = Select_API(
        ctx.obj['sqlite_api'],
        ctx.obj['Neo4j_API'],
        folder=ctx.obj['folder'],
        oscutoff=oscutoff,
        logger=ctx.obj['logger'],
        cutoff=omitcutoff,
        omit_col=omitcol)
    ctx.obj['select_api'] = select_api
    if ctx.obj['Neo4j_API'].sifts_filter is None:
        ctx.obj['Neo4j_API'].sifts_filter = json.loads(siftsfilter)

    if input:
        select_api.process(input).result()
        ctx.obj['logger'].info("Finish select")
    else:
        return processor


if __name__ == '__main__':
    Interface(obj={})

# @Created Date: 2019-11-24 11:03:59 pm
# @Filename: interface.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-05-23 02:24:53 pm
# @Copyright (c) 2019 MinghuiGroup, Soochow University
import click
import configparser
import os
import json
from pandas import read_csv, merge, Series, DataFrame
from numpy import nan
from pathlib import Path
from typing import Iterator
from unsync import Unfuture
import ujson as json
import sys; sys.path.append("C:/GitWorks/pdb-profiling")
from pdb_profiling.processers.uniprot.api import MapUniProtID, UniProtFASTA
from pdb_profiling.processers.pdbe.neo4j_api import Neo4j_API
from pdb_profiling.fetcher.webfetch import UnsyncFetch
from pdb_profiling.log import Abclog
from pdb_profiling.utils import related_dataframe


if "\\" in __file__:
    # Windows
    _SEP = "\\"
else:
    # Linux
    _SEP = "/"


def colorClick(name: str, template: str = "Initializing %s", fg: str = "green"):
    return click.style(template % name, fg=fg)


@click.group(chain=True, invoke_without_command=False)
@click.option("--folder", default="", help="The file folder of new files.", type=click.Path())
@click.option("--loggingpath", default=None, help="The file path of logging.", type=click.Path())
@click.option('--concurReq', type=int, help="the number of concurent requests", default=100)
@click.option('--concurRate', type=float, help="the rate of concurent requests", default=1.5)
@click.option('--useexisting/--no-useexisting', help="whether to use existing result files", default=False, is_flag=True)
@click.pass_context
def Interface(ctx, folder, loggingpath, concurreq, concurrate, useexisting):
    if not folder:
        return
    click.echo(colorClick("Folder"))
    folder = Path(folder)
    ctx.ensure_object(dict)
    ctx.obj['folder'] = folder
    ctx.obj['loggingpath'] = loggingpath
    ctx.obj['concurreq'] = concurreq
    ctx.obj['concurrate'] = concurrate
    for sub in ('mapping', 'fasta'):
        cur_path = folder/'UniProt'/sub
        cur_path.mkdir(parents=True, exist_ok=True)
        ctx.obj[f'{sub}_folder'] = cur_path
    click.echo(colorClick("Logger"))
    ctx.obj['logger'] = Abclog.set_logging_fileHandler(loggingpath, logName='UniProt Module')
    UnsyncFetch.use_existing = useexisting
    

@Interface.resultcallback()
@click.pass_context
def process_pipeline(ctx, processors, folder, loggingpath, concurreq, concurrate, useexisting):
    iterator = ctx.obj.get('iterator', None)
    if iterator is None:
        return
    for processor in processors:
        if processor is not None:
            iterator = processor(iterator)
    UnsyncFetch.unsync_tasks(list(iterator)).result()


@Interface.command("UniProt.id-mapping")
@click.option("--input", 
              default="",
              help="The reference file of IDs that need to map via UniProt RESTful API.", 
              type=click.Path(exists=True))
@click.option("--sep", 
              default="\t", 
              help="The seperator of referenceFile.",
              type=str)
@click.option("--idCol",
              default="RefSeq_protein", 
              help="The column name of IDs in referenceFile.", 
              type=str)
@click.option("--idType", 
              default="P_REFSEQ_AC", 
              help="ID Abbreviation that stands for the type of ID.", 
              type=str)
@click.option("--usecols", 
              default="id,genes,reviewed,comment(ALTERNATIVE%20PRODUCTS),feature(ALTERNATIVE%20SEQUENCE),organism,protein%20names,feature(CHAIN)",
              help="Comma-separated list of the column names for programmatic access to the UniProtKB search results.", 
              type=str)
@click.option("--geneCol",
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
                                        concur_req=ctx.obj['concurreq'],
                                        rate=ctx.obj['concurrate'],
                                        run_tasks=False)
            for task in unsync_tasks:
                yield task

    def processor(iterator: Iterator[Unfuture]):
        '''This processor do nothing'''
        for task in iterator:
            yield task

    click.echo(colorClick("UniProt Retrieve/ID Mapping"))
    df = read_csv(input, sep=sep, chunksize=readchunk, nrows=nrows, skiprows=skiprows)
    if isinstance(df, DataFrame):
        df = (df,)
    click.echo(colorClick("Set Task Iterator", "%s"))
    ctx.obj['iterator'] = yieldTasks(df)
    return processor
    

@Interface.command("UniProt.download-seq")
@click.option("--input", default="", help="the file of UniProt IDs", type=click.Path())
@click.option("--idCol", default=None, help="the column name of UniProt IDs (If not specified, use the first col)", type=str)
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
        UniProtFASTA.retrieve(unps, ctx.obj['fasta_folder'], ctx.obj['concurreq'], ctx.obj['concurrate'])
        return
    else:
        UniProtFASTA.obj = {key: ctx.obj[key]
                            for key in ('fasta_folder', 'concurreq', 'concurrate')}
        return processor


@Interface.command("PDB.SIFTS")
@click.option("--input", default="", help="the file of UniProt Mapping result", type=click.Path())
@click.option("--sep", default="\t", help="the seperator of input file", type=str)
@click.option("--filter",
              default='{"Mapping_status":["eq","Yes"],"Organism":["eq","Homo sapiens (Human)"]}',
              help="the filter(JSON-Format Dict) of UniProt Mapping result", type=str)
@click.pass_context
def sifts(ctx, input, sep, filter):
    def processor(iterator: Iterator[Unfuture]):
        '''Continuations for Retrieving UniProt Fasta Sequence'''
        for task in iterator:
            yield task.then(Neo4j_API.process)

    filter = json.loads(filter)
    if input:
        unp_map_df = read_csv(input, sep, filter)
    else:
        pass
    unp_map_df = related_dataframe(filter, unp_map_df)


if __name__ == '__main__':
    Interface(obj={})

# @Created Date: 2020-04-08 09:44:01 am
# @Filename: neo4j_api.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-04-08 09:44:05 am
# @Copyright (c) 2020 MinghuiGroup, Soochow University
"""
from typing import List, Iterable, Iterator, Union, Dict, Optional, Tuple
from functools import partial, reduce
import pandas as pd
from pathlib import Path
from collections import defaultdict
from itertools import combinations, product
from numpy import nan
import numpy as np
import orjson as json
from functools import lru_cache
from textdistance import jaccard, overlap
from Bio import Align
from Bio.SubsMat import MatrixInfo as matlist
from unsync import unsync, Unfuture
import asyncio
import aiofiles
from tablib import Dataset
import traceback
from pdb_profiling.utils import (pipe_out, sort_sub_cols, slice_series, to_interval, 
                                 lyst22interval, aa_three2one, standardAA, standardNu, range_len,
                                 interval2set, lyst2range, subtract_range,
                                 add_range, overlap_range, outside_range_len, related_dataframe)
from pdb_profiling.log import Abclog
from pdb_profiling.fetcher.dbfetch import Neo4j
from pdb_profiling.processors.pdbe.sqlite_api import Sqlite_API
import logging
# logging.basicConfig(level=logging.INFO)


def lyst2dict(lyst: List) -> Dict:
    try:
        res = dict(lyst)
    except TypeError:
        res = dict((tuple(key), value)for key,value in lyst)
    for key in res.keys():
        res[key] = dict(res[key])
    return res


def sub_index(init_index, subtract_index) -> pd.Index:
    if len(subtract_index) == 0:
        return init_index
    else:
        return pd.Index(set(init_index)-set(subtract_index))


def subtract_dict(left: Dict, right: Dict, copy:bool=True) -> Dict:
    if copy:
        left = left.copy()
    for key in left.keys() & right.keys():
        res = set(left[key]) - set(right[key])
        if res:
            left[key] = list(res)
        else:
            del left[key]
    return left


def tidy_na(dfrm: pd.DataFrame, colName: str, fill, astype):
    dfrm[colName] = dfrm[colName].fillna(fill).apply(astype)


class Entry(object):

    session = None
    
    @classmethod
    def set_session(cls, session):
        cls.session = session
    
    @staticmethod
    def to_data_frame(res):
        return pd.DataFrame(dict(i) for i in res)
    
    @classmethod
    def summary_method(cls, pdbs, session=None):
        query = '''
            MATCH (entry:Entry)-[:EXPERIMENT]->(method:Method)
            WHERE entry.ID in $lyst
            RETURN entry.ID as pdb_id,
                   entry.PDB_REV_DATE_ORIGINAL as PDB_REV_DATE_ORIGINAL, 
                   entry.FIRST_REV_DATE as FIRST_REV_DATE,
                   entry.PDB_REV_DATE as PDB_REV_DATE,
                   entry.REVISION_DATE as REVISION_DATE,
                   tofloat(entry.RESOLUTION) as resolution,
                   method.METHOD_CLASS as METHOD_CLASS
        '''
        if session is None:
            if cls.session is not None:
                session = cls.session
        try:
            return session.run(query, lyst=list(pdbs))
        except AttributeError:
            return query, dict(lyst=list(pdbs))

    @classmethod
    def summary_ligand(cls, pdbs, session=None):
        query = '''
            MATCH (entry:Entry)-[:HAS_BOUND_MOLECULE]->(bmol:BoundMolecule)<-[:IS_PART_OF]-(bli:BoundLigand) 
            WHERE entry.ID in $lyst
            RETURN entry.ID as pdb_id, COUNT(distinct bli) as BOUND_LIGAND_COUNT, COUNT(distinct bmol) as BOUND_MOL_COUNT
        '''
        if session is None:
            if cls.session is not None:
                session = cls.session
        try:
            return session.run(query, lyst=list(pdbs))
        except AttributeError:
            return query, dict(lyst=list(pdbs))
    
    @classmethod
    def summary_nucleotides(cls, pdbs, session=None):
        query = '''
            MATCH (entry:Entry)-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdbRes:PDBResidue)
            WHERE entity.POLYMER_TYPE IN ['D','R','D/R'] AND entry.ID in $lyst
            RETURN 
                entry.ID as pdb_id, 
                entity.ID as entity_id, 
                entity.POLYMER_TYPE as polymer_type,
                size([x in collect(pdbRes.CHEM_COMP_ID) where x in ['DA','DT','DC','DG']]) as DNA_COUNT,
                size([x in collect(pdbRes.CHEM_COMP_ID) where x in ['A','U','C','G']]) as RNA_COUNT,
                size([x in collect(pdbRes.CHEM_COMP_ID) where not x in $standardNu]) as OTHER_COUNT
        '''
        if session is None:
            if cls.session is not None:
                session = cls.session
        try:
            return session.run(query, lyst=list(pdbs), standardNu=standardNu)
        except AttributeError:
            return query, dict(lyst=list(pdbs), standardNu=standardNu)
    
    @staticmethod
    def get_polymer_type(dna, rna):
        if dna and rna:
            return 'D/R'
        elif dna and not rna:
            return 'D'
        elif not dna and rna:
            return 'R'
        elif not dna and not rna:
            return nan

    @classmethod
    def deal_nucleotides(cls, res):
        if isinstance(res, pd.DataFrame):
            dfrm = res
        elif res is None:
            return None
        else:
            dfrm = cls.to_data_frame(res)
        if len(dfrm) == 0:
            return None
        # dfrm['polymer_type'] = dfrm.apply(lambda x: cls.get_polymer_type(x['DNA_COUNT'], x['RNA_COUNT']), axis=1)
        dfrm['polymer_type'] = [cls.get_polymer_type(x, y) for x, y in zip(dfrm['DNA_COUNT'], dfrm['RNA_COUNT'])]
        res = pd.DataFrame(
            ((pdb_id, json.dumps(dict(zip(data.entity_id, data.polymer_type))).decode('utf-8'))
             for pdb_id, data in dfrm.groupby('pdb_id')),
            columns=('pdb_id', 'nucleotides_entity_type'))
        # res['has_hybrid_nucleotides'] = res.nucleotides_entity_type.apply(lambda x: 'D\/R' in x)
        res['has_hybrid_nucleotides'] = ['D/R' in x for x in res.nucleotides_entity_type.values]
        return res

    @classmethod
    def summary_seq(cls, pdbs, session=None):
        query = '''
            MATCH (entry:Entry)-[:HAS_ENTITY]->(entity:Entity{POLYMER_TYPE:'P'})-[:CONTAINS_CHAIN]->(chain:Chain)-[inChain:IS_IN_CHAIN]-(res:PDBResidue)
            WHERE entry.ID in $lyst
            WITH 
                entry.ID as pdb_id,
                entity.ID as entity_id,
                chain.AUTH_ASYM_ID as chain_id,
                COUNT(res) as SEQRES_COUNT,
                avg(tofloat(inChain.OBSERVED_RATIO)) as AVG_OBS_RATIO,
                [x in COLLECT([res.ID, inChain.AUTH_COMP_ID]) WHERE NOT x[1] IN $standardAA] as NON_INDEX,
                [x in COLLECT([res.ID, inChain.AUTH_COMP_ID]) WHERE x[1] = 'UNK'] as UNK_INDEX,
                [x in COLLECT([res.ID, inChain.OBSERVED]) WHERE x[1]='N'] as MIS_INDEX,
                chain
            MATCH (chain)-[obs_inChain:IS_IN_CHAIN]-()
           	WHERE tofloat(obs_inChain.OBSERVED_RATIO) > 0
            RETURN pdb_id, entity_id, chain_id, SEQRES_COUNT, AVG_OBS_RATIO, avg(tofloat(obs_inChain.OBSERVED_RATIO)) as AVG_OBS_OBS_RATIO, NON_INDEX, UNK_INDEX, MIS_INDEX
        '''
        if session is None:
            if cls.session is not None:
                session = cls.session
        try:
            return session.run(query, lyst=list(pdbs), standardAA=standardAA)
        except AttributeError:
            return query, dict(lyst=list(pdbs), standardAA=standardAA)

    @classmethod
    def summary_eec(cls, pdbs, protein_only:bool=True):
        query = '''
            MATCH (entry:Entry)-[:HAS_ENTITY]-(entity:Entity{})-[:CONTAINS_CHAIN|:IS_AN_INSTANCE_OF]-(chain)
            WHERE entry.ID IN $pdbs
            RETURN entry.ID as pdb_id, entity.ID as entity_id, chain.AUTH_ASYM_ID as chain_id, chain.STRUCT_ASYM_ID as struct_asym_id, entity.POLYMER_TYPE as POLYMER_TYPE
        '''
        if protein_only:
            query = query.format("{POLYMER_TYPE:'P'}")
        else:
            query = query.format("")
        return query, dict(pdbs=list(pdbs))

    @staticmethod
    def deal_seq_index(dfrm: pd.DataFrame) -> pd.DataFrame:
        def index2range(lyst):
            res = to_interval(int(i[0]) for i in lyst)
            if isinstance(res, List):
                return json.dumps(res).decode('utf-8')
            else:
                return res
        
        def get_obs_rel_count(seqres_count, mis_range, unk_range, non_range):
            seq_range = [(1, seqres_count)]
            obs_range = subtract_range(seq_range, mis_range)
            obs_unk_range = overlap_range(obs_range, unk_range)
            obs_pure_range = subtract_range(obs_range, non_range)
            return range_len(obs_unk_range), range_len(obs_pure_range)

        dfrm['UNK_COUNT'] = dfrm.UNK_INDEX.apply(len)
        dfrm['PURE_SEQRES_COUNT'] = dfrm.SEQRES_COUNT - dfrm.NON_INDEX.apply(len)
        dfrm['OBS_RECORD_COUNT'] = dfrm.SEQRES_COUNT - dfrm.MIS_INDEX.apply(len)
        dfrm.NON_INDEX = [index2range(x) for x in dfrm.NON_INDEX]
        dfrm.UNK_INDEX = [index2range(x) for x in dfrm.UNK_INDEX]
        dfrm.MIS_INDEX = [index2range(x) for x in dfrm.MIS_INDEX]
        dfrm[['OBS_UNK_COUNT', 'ATOM_RECORD_COUNT']] = dfrm.apply(lambda x: get_obs_rel_count(
            x['SEQRES_COUNT'], x['MIS_INDEX'], x['UNK_INDEX'], x['NON_INDEX']), axis=1, result_type='expand')
        return dfrm

    @classmethod
    @lru_cache()
    def get_seqres(cls, pdb_id: str, entity_id: str = None, chain_id: str = None, three2one: bool = True, session=None):
        '''
        TODO: Different behavior
        '''
        if session is None:
            if cls.session is not None:
                session = cls.session
        if entity_id is not None:
            query = '''
                MATCH (entity:Entity)-[:HAS_PDB_RESIDUE]->(pdbRes:PDBResidue)
                WHERE entity.UNIQID = $uniqid
                RETURN pdbRes.CHEM_COMP_ID AS residue_name ORDER BY toInteger(pdbRes.ID)
            '''
            # res = session.run(query, uniqid=f'{pdb_id}_{entity_id}')
            kwargs = dict(uniqid=f'{pdb_id}_{entity_id}')
        elif chain_id is not None:
            query = '''
                MATCH (entry:Entry)-[:HAS_ENTITY]-(:Entity)-[:CONTAINS_CHAIN]-(chain:Chain)-[:IS_IN_CHAIN]-(pdbRes:PDBResidue)
                WHERE entry.ID = $pdb_id AND chain.AUTH_ASYM_ID = $chain_id
                RETURN pdbRes.CHEM_COMP_ID AS residue_name ORDER BY toInteger(pdbRes.ID)
            '''
            # res = session.run(query, pdb_id=pdb_id, chain_id=chain_id)
            kwargs = dict(pdb_id=pdb_id, chain_id=chain_id)
        else:
            raise ValueError('please specify entity_id or chain_id')
        
        try:
            res = session.run(query, **kwargs)
            if three2one:
                return ''.join(aa_three2one.get(r['residue_name'], 'X') for r in res)
            else:
                return res
        except AttributeError:
            return query, kwargs

    @classmethod
    def get_residues(cls, pdb_id: str, entity_id: str, chain_id: str, res_ids=None, observed_only:bool=False, session=None):
        if session is None:
            if cls.session is not None:
                session = cls.session
        query = '''
            MATCH (entry:Entry)-[:HAS_ENTITY]->(entity:Entity)-[:CONTAINS_CHAIN]->(chain:Chain)-[inChain:IS_IN_CHAIN]-(res:PDBResidue)
            WHERE entry.ID = $pdb_id AND entity.ID = $entity_id AND chain.AUTH_ASYM_ID = $chain_id {}
            RETURN entry.ID as pdb_id, entity.ID as entity_id, chain.AUTH_ASYM_ID as chain_id, res.CHEM_COMP_ID as residue_name, toInteger(res.ID) as residue_number, tofloat(inChain.OBSERVED_RATIO) as obs_ratio, toInteger(inChain.AUTH_SEQ_ID) as author_residue_number, inChain.PDB_INS_CODE as author_insertion_code ORDER BY residue_number
        '''
        if res_ids is None:
            if observed_only:
                query = query.format("AND inChain.OBSERVED = 'Y'")
            else:
                query = query.format('')
            # return session.run(query, pdb_id=pdb_id, entity_id=str(entity_id), chain_id=chain_id)
            kwargs = dict(pdb_id=pdb_id, entity_id=str(entity_id), chain_id=chain_id)
        else:
            if observed_only:
                query = query.format("AND res.ID IN $res_ids AND inChain.OBSERVED = 'Y'")
            else:
                query = query.format('AND res.ID IN $res_ids')
            res_ids = [str(i) for i in res_ids]
            # return session.run(query, pdb_id=pdb_id, entity_id=str(entity_id), chain_id=chain_id, res_ids=res_ids)
            kwargs = dict(pdb_id=pdb_id, entity_id=str(entity_id), chain_id=chain_id, res_ids=res_ids)
        try:
            return session.run(query, **kwargs)
        except AttributeError:
            return query, kwargs

    @classmethod
    def get_binding_redidues(cls, pdb_id: str, entity_id: str, chain_id: str, count_only: bool = True, session=None):
        if session is None:
            if cls.session is not None:
                session = cls.session
        if count_only:
            query = '''
                MATCH (entry:Entry)-[:HAS_BOUND_MOLECULE]->(:BoundMolecule)<-[:IS_PART_OF]-(:BoundLigand)-[:HAS_ARP_CONTACT]-(res:PDBResidue)-[:IS_IN_CHAIN]-(chain:Chain)-[:CONTAINS_CHAIN]-(entity:Entity)
                WHERE entry.ID = $pdb_id AND entity.ID = $entity_id AND chain.AUTH_ASYM_ID = $chain_id
                RETURN COUNT(DISTINCT res) as LIGAND_BINDING_RES_COUNT
            '''
        else:
            query = '''
                MATCH (entry:Entry)-[:HAS_BOUND_MOLECULE]->(:BoundMolecule)<-[:IS_PART_OF]-(:BoundLigand)-[:HAS_ARP_CONTACT]-(res:PDBResidue)-[inChain:IS_IN_CHAIN]-(chain:Chain)-[:CONTAINS_CHAIN]-(entity:Entity)
                WHERE entry.ID = $pdb_id AND entity.ID = $entity_id AND chain.AUTH_ASYM_ID = $chain_id
                RETURN distinct entry.ID as pdb_id, 
                                entity.ID as entity_id, 
                                chain.AUTH_ASYM_ID as chain_id, 
                                res.CHEM_COMP_ID as residue_name, 
                                toInteger(res.ID) as residue_number, 
                                tofloat(inChain.OBSERVED_RATIO) as obs_ratio, 
                                toInteger(inChain.AUTH_SEQ_ID) as author_residue_number, 
                                inChain.PDB_INS_CODE as author_insertion_code
                '''
        kwargs = dict(
            pdb_id=pdb_id, 
            entity_id=str(entity_id), 
            chain_id=chain_id)
        try:
            return session.run(query, **kwargs)
        except AttributeError:
            return query, kwargs

    @classmethod
    def summary_entity_chain(cls, pdbs, session=None):
        query = '''
                MATCH (entry:Entry)-[:HAS_ENTITY]-(entity:Entity{POLYMER_TYPE:'P'})-[:CONTAINS_CHAIN]-(chain:Chain)
                WHERE entry.ID IN $pdbs
                WITH entry, entity, COLLECT(chain.AUTH_ASYM_ID) as chain_ids
                RETURN entry.ID as pdb_id, COLLECT([entity.ID, chain_ids]) as entity_chain_map
            '''
        if session is None:
            if cls.session is not None:
                session = cls.session
        try:
            return session.run(query, pdbs=list(pdbs))
        except AttributeError:
            return query, dict(pdbs=list(pdbs))
    
    @classmethod
    def deal_entity_chain(cls, res):
        if isinstance(res, pd.DataFrame):
            dfrm = res
        else:
            dfrm = cls.to_data_frame(res)
            dfrm.entity_chain_map = dfrm.entity_chain_map.apply(dict)
        dfrm['entity_count'] = dfrm.entity_chain_map.apply(len)
        dfrm['chain_count'] = dfrm.entity_chain_map.apply(lambda x: sum(len(i) for i in x.values()))
        return dfrm
 
    @classmethod
    def summary_entity_interaction(cls, pdbs, session=None):
        '''
        Target: 
        
        * Whether two chains have interaction, if yes, find the interacting residues;
        * Then check whether the interacting residues are in the mapped range of the corresponding uniprot-isoform

        Possible Solution:

        1. Find the interacting entity
        2. Fing the Interface Node
        3. Directly go to the residues with contact
        4. Use InterfaceBindingResidue Node
        5. Use MAP_TO_UNIPROT_RESIDUE Relationship

        Problem to solve:

        1. Difference between ARP_CONTACT and PISA_CONTACT
        2. 
        '''
        query = '''
                MATCH (entry:Entry)-[:HAS_ENTITY]-(entityA:Entity{POLYMER_TYPE:'P'})-[:INTERACTS_WITH_ARP|:INTERACTS_WITH_PISA]-(entityB:Entity{POLYMER_TYPE:'P'})
                WHERE entry.ID IN $pdbs
                RETURN entityA, entityB
            '''
        if session is None:
            if cls.session is not None:
                session = cls.session
        try:
            return session.run(query, pdbs=list(pdbs))
        except AttributeError:
            return query, dict(lyst=list(pdbs))

    @classmethod
    def summary_chain_interaction(cls, pdbs, source:str):
        query = '''
        MATCH (entry:Entry)-[:HAS_ENTITY]-(:Entity{POLYMER_TYPE:'P'})-[:HAS_PDB_RESIDUE]-(:PDBResidue)-[has_contact:%s]-(:PDBResidue)
        WHERE entry.ID IN $pdbs AND has_contact.AUTH_ASYM_ID_1 <> has_contact.AUTH_ASYM_ID_2
        RETURN distinct entry.ID as pdb_id, has_contact.AUTH_ASYM_ID_1 as chain_id_1, has_contact.AUTH_ASYM_ID_2 as chain_id_2, true as from_%s
        '''
        # //COLLECT(has_contact.TYPE) as type
        source = source.lower()
        if source == 'pisa_bond':
            query = query % ("HAS_PISA_BOND", source)
        elif source == 'arp_contact':
            query = query % ("HAS_ARP_CONTACT", source)
        else:
            raise ValueError("Invalid source! Please input specified source within ('pisa_bond', 'arp_contact')")
        return query, dict(pdbs=list(pdbs))

    @classmethod
    def summary_assembly(cls, pdbs, protein_only:bool=True):
        query = '''
        MATCH (entry:Entry)-[:HAS_ENTITY]-(entity:Entity{})-[:CONTAINS_CHAIN]-(chain:Chain)-[:IS_PART_OF_ASSEMBLY]-(ass:Assembly)
        WHERE entry.ID IN $pdbs
        WITH entry,entity,chain,ass
        MATCH (entity:Entity)-[inass:IS_PART_OF_ASSEMBLY]-(ass:Assembly)
        RETURN
            entry.ID as pdb_id,
            entity.ID as entity_id,
            entity.POLYMER_TYPE as POLYMER_TYPE,
            chain.AUTH_ASYM_ID as chain_id,
            ass.UNIQID as biounit_id, 
            ass.PREFERED as PREFERED,
            inass.NUMBER_OF_CHAINS as NUMBER_OF_CHAINS
        '''
        if protein_only:
            query = query.format("{POLYMER_TYPE:'P'}")
        else:
            query = query.format("")
        return query, dict(pdbs=list(pdbs))
    
    @classmethod
    def summary_assembly_interface(cls, pdbs, protein_only:bool=True):
        query = '''
        MATCH (entry:Entry)-[:HAS_ENTITY]-(entity:Entity{})-[:CONTAINS_CHAIN]-(chain:Chain)-[:IS_PART_OF_ASSEMBLY]-(ass:Assembly)-[:HAS_INTERFACE]-(interface:Interface)
        WHERE entry.ID IN $pdbs
        WITH entry,entity,chain,ass,interface
        MATCH (interface:Interface)-[:HAS_ENTITY]-(entity:Entity)-[inass:IS_PART_OF_ASSEMBLY]-(ass:Assembly)
        RETURN
            entry.ID as pdb_id,
            entity.ID as entity_id,
            entity.POLYMER_TYPE as POLYMER_TYPE,
            chain.AUTH_ASYM_ID as chain_id,
            ass.UNIQID as biounit_id, 
            ass.PREFERED as PREFERED,
            inass.NUMBER_OF_CHAINS as NUMBER_OF_CHAINS,
            interface.UNIQID as interface_id
        '''
        if protein_only:
            query = query.format("{POLYMER_TYPE:'P'}")
        else:
            query = query.format("")
        return query, dict(pdbs=list(pdbs))


class SIFTS(Entry):
    @classmethod
    def summary_entity_unp(cls, pdbs, tax_id: Optional[str]=None, session=None):
        query = '''
                MATCH (entry:Entry)-[:HAS_ENTITY]-(entity:Entity{POLYMER_TYPE:'P'})-[unp_rel:HAS_UNIPROT{BEST_MAPPING:"1"}]-(unp:UniProt)-[:HAS_TAXONOMY]-(tax:Taxonomy)
                WHERE entry.ID in $pdbs %s
                WITH entry, entity, unp, tax
                MATCH (unp)-[seg:HAS_UNIPROT_SEGMENT]-(entity)
                WITH COLLECT(DISTINCT [toInteger(seg.UNP_START), toInteger(seg.UNP_END)]) as range_info, unp, entry, entity, tax
                WITH COLLECT([entity.ID, range_info]) as eneitiy_unp_info, entry, unp, tax
                RETURN entry.ID as pdb_id, COLLECT([[unp.ACCESSION,tax.TAX_ID], eneitiy_unp_info]) as unp_entity_info
            '''
        if tax_id is None:
            query = query % ''
        else:
            query = query % f'AND tax.ID = "{tax_id}"'
        if session is None:
            if cls.session is not None:
                session = cls.session
        try:
            return session.run(query, pdbs=list(pdbs))
        except AttributeError:
            return query, dict(pdbs=list(pdbs))
    
    @classmethod
    def deal_entity_unp(cls, res):
        def reverse_dict(info_dict: Dict) -> Dict:
            res = defaultdict(list)
            for key, valueDict in info_dict.items():
                for innerKey in valueDict.keys():
                    res[innerKey].append(key)
            return res

        def min_unp(info_dict: Dict) -> int:
            return min(len(set(res)) for res in product(*info_dict.values()))
        
        if isinstance(res, pd.DataFrame):
            dfrm = res
        else:
            dfrm = cls.to_data_frame(res)
            dfrm.unp_entity_info = dfrm.unp_entity_info.apply(lambda x: lyst2dict(x) if not isinstance(x, float) else nan)
        # !!!
        dfrm['entity_unp_info'] = dfrm.unp_entity_info.apply(lambda x: reverse_dict(x) if not isinstance(x, float) else nan)
        dfrm['entity_with_unp_count'] = dfrm.entity_unp_info.apply(lambda x: len(x) if not isinstance(x, float) else nan)
        dfrm['min_unp_count'] = dfrm.entity_unp_info.apply(lambda x: min_unp(x) if not isinstance(x, float) else nan)
        return dfrm
    
    @staticmethod
    def get_oligo_state(df: pd.DataFrame):
        state_name, unmapped_name = 'oligo_state', 'has_unmapped_protein'
        assert len(df[df.entity_count < df.entity_with_unp_count]) == 0
        cur_df = df
        # He
        he = cur_df[
            ((cur_df.entity_count > cur_df.entity_with_unp_count) &
             (cur_df.min_unp_count.eq(1))) | (cur_df.min_unp_count.gt(1))
        ].index
        cur_df = cur_df.loc[sub_index(cur_df.index, he)]
        # Ho
        ho = cur_df[
            (cur_df.entity_count == cur_df.entity_with_unp_count)
            & (cur_df.min_unp_count.eq(1))
            & (cur_df.chain_count.gt(1))
        ].index
        cur_df = cur_df.loc[sub_index(cur_df.index, ho)]
        # Mo
        mo = df[df.chain_count.eq(1)].index
        cur_df_index = sub_index(cur_df.index, mo)
        # Tag
        df[state_name] = nan
        df.loc[mo, state_name] = 'mo'
        df.loc[ho, state_name] = 'ho'
        df.loc[he, state_name] = 'he'
        # Tag
        non_ho = df.loc[sub_index(df.index, ho)]
        unmapped = non_ho[non_ho.entity_count > non_ho.entity_with_unp_count].index
        df[unmapped_name] = False
        df.loc[unmapped, unmapped_name] = True
        # return
        return df, (mo, ho, he, cur_df_index)

    @classmethod
    def summary_oligo_state(cls, pdbs):
        entity_chain = cls.deal_entity_chain(cls.summary_entity_chain(pdbs))
        unp_entity = cls.deal_entity_unp(cls.summary_entity_unp(pdbs))
        unp_entity_chain = pd.merge(entity_chain, unp_entity, how='left')
        tidy_na(unp_entity_chain, 'entity_with_unp_count', 0, int)
        tidy_na(unp_entity_chain, 'min_unp_count', 0, int)
        # pass_result_oli, indexes_oli
        return cls.get_oligo_state(unp_entity_chain)
    
    @classmethod
    @unsync
    async def asummary_oligo_state(cls, pdbs, neo4j):
        query, kwargs = cls.summary_entity_chain(pdbs)
        res = await neo4j.afetch(query, **kwargs)
        entity_chain = cls.deal_entity_chain(res)
        query, kwargs = cls.summary_entity_unp(pdbs)
        res = await neo4j.afetch(query, **kwargs)
        unp_entity = cls.deal_entity_unp(res)
        unp_entity_chain = pd.merge(entity_chain, unp_entity, how='left')
        tidy_na(unp_entity_chain, 'entity_with_unp_count', 0, int)
        tidy_na(unp_entity_chain, 'min_unp_count', 0, int)
        return cls.get_oligo_state(unp_entity_chain)

    @classmethod
    @unsync
    async def pipe_oligo(cls, pdbs, neo4j, sqlite, **kwargs):
        seqres_df = await sqlite.SEQRES_Info.objects.filter(pdb_id__in=pdbs).all()
        seqres_df = pd.DataFrame(seqres_df)
        omit_df = cls.omit_chains(seqres_df, **kwargs)
        oligo_df, _ = await cls.asummary_oligo_state(pdbs, neo4j)
        oligo_df, _ = cls.update_oligo_state(oligo_df, omit_df)
        return oligo_df

    @staticmethod
    def omit_chains(seq_df: pd.DataFrame, cutoff: int = 50, omit_col: str = 'ATOM_RECORD_COUNT'):
        focus = seq_df[seq_df[omit_col].lt(cutoff)]
        return pd.DataFrame(
            ((pdb_id, dict((entity_id, entityData.chain_id.to_list()) for entity_id,
                           entityData in pdbData.groupby('entity_id'))) for pdb_id, pdbData in focus.groupby('pdb_id')),
            columns=('pdb_id', 'del_entity_chain')
        )

    @classmethod
    def update_oligo_state(cls, oligo_df: pd.DataFrame, omit_df: pd.DataFrame):
        def new_entity_chain_map(entity_chain_map, del_entity_chain) -> Dict:
            if isinstance(del_entity_chain, Dict):
                return subtract_dict(entity_chain_map, del_entity_chain)
            else:
                return entity_chain_map

        def rest_unp_entity(unp_entity_info, entity_chain_map: Dict):
            if not isinstance(unp_entity_info, Dict):
                return nan
            new_unp_entity_info = defaultdict(dict)
            for unp, valuedict in unp_entity_info.items():
                entities = valuedict.keys() & entity_chain_map.keys()
                if not entities:
                    continue
                for entity in entities:
                    new_unp_entity_info[unp][entity] = unp_entity_info[unp][entity]
            return new_unp_entity_info
        
        oligo_df = pd.merge(oligo_df, omit_df, how='left')
        oligo_df.entity_chain_map = oligo_df.apply(lambda x: new_entity_chain_map(
            x['entity_chain_map'], x['del_entity_chain']), axis=1)
        oligo_df.unp_entity_info = oligo_df.apply(lambda x: rest_unp_entity(
            x['unp_entity_info'], x['entity_chain_map']), axis=1)
        oligo_df = cls.deal_entity_unp(cls.deal_entity_chain(oligo_df))
        tidy_na(oligo_df, 'entity_with_unp_count', 0, int)
        tidy_na(oligo_df, 'min_unp_count', 0, int)
        return cls.get_oligo_state(oligo_df)

    @classmethod
    def set_from(cls, id_type: str):
        id_type = id_type.lower()
        if id_type in ('unp', 'uniprot', 'unp.ACCESSION'):
            cls.from_str = 'unp.ACCESSION'
        elif id_type in ('pdb', 'pdb_id', 'entry.ID'):
            cls.from_str = 'entry.ID'
        else:
            raise ValueError('unknown id type')

    @classmethod
    def summary_mapping(cls, lyst, id_type: str, session=None):
        cls.set_from(id_type)
        query = '''
            MATCH (entry:Entry)-[:HAS_ENTITY]-(entity:Entity)-[seg:HAS_UNIPROT_SEGMENT]-(unp:UniProt)
            WHERE {} in $lyst
            RETURN unp.ACCESSION as UniProt, entry.ID as pdb_id, entity.ID as entity_id, seg.AUTH_ASYM_ID as chain_id, tofloat(seg.IDENTITY) as identity, COLLECT([toInteger(seg.PDB_START), toInteger(seg.PDB_END)]) as pdb_range, COLLECT([toInteger(seg.UNP_START), toInteger(seg.UNP_END)]) as unp_range
        '''.format(cls.from_str)
        if session is None:
            if cls.session is not None:
                session = cls.session
        try:
            return session.run(query, lyst=list(lyst))
        except AttributeError:
            return query, dict(lyst=list(lyst))
    
    @classmethod
    def summary_mapping_detail(cls, pdb_id, entity_id, unp):
        query = '''
        MATCH (entry:Entry)-[:HAS_ENTITY]-(entity:Entity)-[:HAS_PDB_RESIDUE]-(pdbres:PDBResidue)-[:MAP_TO_UNIPROT_RESIDUE]-(unpres:UNPResidue)-[:HAS_UNP_RESIDUE]-(unp:UniProt)
        WHERE unp.ACCESSION = $unp AND entity.ID = $entity_id AND entry.ID = $pdb_id
        RETURN toInteger(pdbres.ID) as residue_number,
               toInteger(SPLIT(unpres.UNIQID, '_')[1]) as unp_index
        '''
        return query, dict(pdb_id=pdb_id, entity_id=entity_id, unp=unp)

    @classmethod
    @unsync
    async def deal_mapping(cls, res, neo4j_api):
        dfrm = cls.to_data_frame(res)
        if len(dfrm) == 0:
            return
        dfrm = cls.deal_InDe(dfrm)
        dfrm.pdb_range = dfrm.pdb_range.apply(
            lambda x: json.dumps(x).decode('utf-8'))
        dfrm.unp_range = dfrm.unp_range.apply(
            lambda x: json.dumps(x).decode('utf-8'))
        dfrm = await cls.update_range(dfrm, neo4j_api)
        # dfrm = await cls.update_range_local_align(dfrm, neo4j_api)
        return dfrm

    @classmethod
    def extra_mapping(cls, entity_uniqids, unps) -> pd.DataFrame:
        '''
        TODO: Different behavior
        '''
        query = '''
            MATCH (entity:Entity)-[seg:HAS_UNIPROT_SEGMENT]->(unp:UniProt)
            WHERE entity.UNIQID IN $entity_uniqids AND unp.ACCESSION IN $unps
            RETURN unp.ACCESSION as UniProt,
                   substring(entity.UNIQID,0,4) as pdb_id, 
                   entity.ID as entity_id, 
                   seg.AUTH_ASYM_ID as chain_id, 
                   tofloat(seg.IDENTITY) as identity, 
                   COLLECT([toInteger(seg.PDB_START), 
                   toInteger(seg.PDB_END)]) as pdb_range, 
                   COLLECT([toInteger(seg.UNP_START), 
                   toInteger(seg.UNP_END)]) as unp_range
        '''
        '''
        if session is None:
            if cls.session is not None:
                session = cls.session
        try:
            # TODO: this part of code is wrong right now, need to be fix
            res = session.run(query,
                              entity_uniqids=list(entity_uniqids),
                              unps=list(unps))
            return cls.deal_mapping(res)
        except AttributeError:
        '''
        return query, dict(entity_uniqids=list(entity_uniqids), unps=list(unps))

    @classmethod
    def pdbchain2unpres(cls, lyst, observed_only: bool = False):
        query = '''
        MATCH (chain:Chain)-[inChain:IS_IN_CHAIN]-(pdbres:PDBResidue)-[:MAP_TO_UNIPROT_RESIDUE]-(unpres:UNPResidue)-[:HAS_UNP_RESIDUE]-(unp:UniProt)
        WHERE chain.UNIQID in $lyst {}
        RETURN chain.UNIQID as UNIQID,
               chain.AUTH_ASYM_ID as chain_id,
               pdbres.CHEM_COMP_ID as residue_name, 
               toInteger(pdbres.ID) as residue_number,
               tofloat(inChain.OBSERVED_RATIO) as obs_ratio,
               toInteger(inChain.AUTH_SEQ_ID) as author_residue_number,
               inChain.PDB_INS_CODE as author_insertion_code,
               unpres.UNIQID as unp_index,
               unpres.ONE_LETTER_CODE as unp_res
        '''
        if observed_only:
            query = query.format("AND inChain.OBSERVED = 'Y'")
        else:
            query = query.format('')

        return query, dict(lyst=list(lyst))

    @staticmethod
    def sort_2_range(unp_range: List, pdb_range: List):
        unp_range, pdb_range = zip(*sorted(zip(unp_range, pdb_range), key=lambda x: x[0][0]))
        return unp_range, pdb_range

    @classmethod
    def deal_InDe(cls, dfrm: pd.DataFrame) -> pd.DataFrame:
        def get_gap_list(li: List):
            return [li[i+1][0] - li[i][1] - 1 for i in range(len(li)-1)]

        def get_range_diff(lyst_a: List, lyst_b: List):
            array_a = np.array([right - left + 1 for left, right in lyst_a])
            array_b = np.array([right - left + 1 for left, right in lyst_b])
            return (array_a - array_b).tolist()

        def add_tage_to_range(df: pd.DataFrame, tage_name: str):
            # ADD TAGE FOR SIFTS
            df[tage_name] = 'Safe'
            # No Insertion But Deletion[Pure Deletion]
            df.loc[df[(df['group_info'] == 1) & (
            df['unp_pdb_var'] > 0)].index, tage_name] = 'Deletion'
            # Insertion & No Deletion
            df.loc[df[
                (df['group_info'] != 1) &
                (df['var_0_count'] == df['group_info']) &
                (df['unp_gap_0_count'] == (df['group_info'] - 1))].index, tage_name] = 'Insertion'
            # Insertion & Deletion
            df.loc[df[
                (df['group_info'] != 1) &
                ((df['var_0_count'] != df['group_info']) |
                (df['unp_gap_0_count'] != (df['group_info'] - 1)))].index, tage_name] = 'Insertion & Deletion'

        dfrm['group_info'] = dfrm.apply(lambda x: len(
            x['pdb_range']), axis=1)
        
        focus_index = dfrm[dfrm.group_info.gt(1)].index
        if len(focus_index) > 0: 
            focus_df = dfrm.loc[focus_index].apply(lambda x: cls.sort_2_range(x['unp_range'], x['pdb_range']), axis=1, result_type='expand')
            focus_df.index = focus_index
            focus_df.columns = ['unp_range', 'pdb_range']
            dfrm.loc[focus_index, ['unp_range', 'pdb_range']] = focus_df
        
        dfrm['pdb_gap_list'] = dfrm.apply(lambda x: json.dumps(
            get_gap_list(x['pdb_range'])).decode('utf-8'), axis=1)
        dfrm['unp_gap_list'] = dfrm.apply(lambda x: json.dumps(
            get_gap_list(x['unp_range'])).decode('utf-8'), axis=1)
        dfrm['var_list'] = dfrm.apply(lambda x: json.dumps(get_range_diff(
            x['unp_range'], x['pdb_range'])).decode('utf-8'), axis=1)
        dfrm['repeated'] = dfrm.apply(
            lambda x: '-' in x['var_list'], axis=1)
        dfrm['repeated'] = dfrm.apply(
            lambda x: True if '-' in x['unp_gap_list'] else x['repeated'], axis=1)
        dfrm['var_0_count'] = dfrm.apply(
            lambda x: json.loads(x['var_list']).count(0), axis=1)
        dfrm['unp_gap_0_count'] = dfrm.apply(
            lambda x: json.loads(x['unp_gap_list']).count(0), axis=1)
        dfrm['unp_pdb_var'] = dfrm.apply(
            lambda x: json.loads(x['var_list'])[0], axis=1)
        add_tage_to_range(dfrm, tage_name='sifts_range_tag')
        return dfrm

    @classmethod
    @unsync
    async def update_range(cls, dfrm: pd.DataFrame, neo4j_api) -> pd.DataFrame:
        new_unp_range, new_pdb_range = 'new_unp_range', 'new_pdb_range'
        focus_df = dfrm[
            (dfrm.sifts_range_tag.isin(('Deletion', 'Insertion & Deletion')))
            & (dfrm.repeated.eq(False))]
        focus = focus_df[['pdb_id', 'entity_id', 'UniProt']].to_numpy()
        focus_index = focus_df.index
        updated_pdb_range, updated_unp_range = list(), list()
        for record in focus:
            query, kwargs = cls.summary_mapping_detail(*record)
            res = await neo4j_api.afetch(query, **kwargs)
            res = cls.to_data_frame(res)
            res_unp_intervel, res_pdb_intervel = lyst22interval(res.unp_index, res.residue_number)
            updated_unp_range.append(json.dumps(res_unp_intervel).decode('utf-8'))
            updated_pdb_range.append(json.dumps(res_pdb_intervel).decode('utf-8'))
        updated_range_df = pd.DataFrame(
            {new_unp_range: updated_unp_range, new_pdb_range: updated_pdb_range}, index=focus_index)
        dfrm = pd.merge(dfrm, updated_range_df, left_index=True,
                        right_index=True, how='left')
        dfrm[new_unp_range] = dfrm.apply(lambda x: x['unp_range'] if pd.isna(
            x[new_unp_range]) else x[new_unp_range], axis=1)
        dfrm[new_pdb_range] = dfrm.apply(lambda x: x['pdb_range'] if pd.isna(
            x[new_pdb_range]) else x[new_pdb_range], axis=1)
        return dfrm

    @classmethod
    @unsync
    async def update_range_local_align(cls, dfrm: pd.DataFrame, neo4j_api) -> pd.DataFrame:
        new_unp_range, new_pdb_range = 'new_unp_range_align', 'new_pdb_range_align'
        focus_df = dfrm[
            (dfrm.sifts_range_tag.isin(('Deletion', 'Insertion & Deletion')))
            & (dfrm.repeated.eq(False))]
        focus = focus_df[['pdb_id', 'entity_id', 'UniProt']].to_numpy()
        focus_index = focus_df.index
        updated_pdb_range, updated_unp_range = list(), list()
        seqAligner = SeqPairwiseAlign()
        for record in focus:
            # pdbSeq = cls.get_seqres(record['pdb_id'], record['entity_id'])
            query, kwargs = cls.get_seqres(record[0], record[1])
            res = await neo4j_api.afetch(query, **kwargs)
            pdbSeq = ''.join(aa_three2one.get(r['residue_name'], 'X') for r in res)
            # unpSeq = cls.get_unp_seq(record["UniProt"])
            query, kwargs = cls.get_unp_seq(record[2])
            res = await neo4j_api.afetch(query, **kwargs)
            unpSeq = ''.join(r['aa'] for r in res)
            # ---
            res = seqAligner.makeAlignment(unpSeq, pdbSeq)
            updated_unp_range.append(res[0])
            updated_pdb_range.append(res[1])

        updated_range_df = pd.DataFrame(
            {new_unp_range: updated_unp_range, new_pdb_range: updated_pdb_range}, index=focus_index)
        dfrm = pd.merge(dfrm, updated_range_df, left_index=True,
                        right_index=True, how='left')
        dfrm[new_unp_range] = dfrm.apply(lambda x: x['unp_range'] if pd.isna(
            x[new_unp_range]) else x[new_unp_range], axis=1)
        dfrm[new_pdb_range] = dfrm.apply(lambda x: x['pdb_range'] if pd.isna(
            x[new_pdb_range]) else x[new_pdb_range], axis=1)
        return dfrm

    @classmethod
    @lru_cache()
    def get_unp_seq(cls, unp: str, session=None):
        '''
        TODO: Different behavior
        '''
        query = '''
            MATCH (unp:UniProt{ACCESSION: "%s"})-[:HAS_UNP_RESIDUE]-(unpRes:UNPResidue)
            RETURN unpRes.ONE_LETTER_CODE AS aa ORDER BY toInteger(unpRes.ID)
        ''' % unp
        if session is None:
            if cls.session is not None:
                session = cls.session
        try:
            return ''.join(r['aa'] for r in session.run(query))
        except AttributeError:
            return query, {}
    
    @classmethod
    def get_unp_seq_len(cls, unp: str, session=None):
        query = '''
            MATCH (unp:UniProt{ACCESSION: "%s"})-[:HAS_UNP_RESIDUE]-(unpRes:UNPResidue)
            RETURN COUNT(unpRes) as unp_len
        ''' % unp
        if session is None:
            if cls.session is not None:
                session = cls.session
        try:
            return session.run(query)
        except AttributeError:
            return query, {}

    @classmethod
    def res_conflict_then(cls, res):
        def get_mutaRes(dfrm: pd.DataFrame):
            sites = defaultdict(list)
            def storeSites(lyst, sitelyst):
                sitelyst.extend([i.split('|')[1] for i in lyst if '|' in i])# !!!
            dfrm.apply(lambda x: storeSites(x['resCon.DETAILS'], sites[(x['pdb_id'], x['entity_id'], x['resCon.ID'])]), axis=1)
            return sites

        def yield_dfrm_of_mutaRes(sites: Dict):
            for (pdb, entity, tag), value in sites.items():
                cur_df = pd.DataFrame(set(value), columns=['residue_number'])
                cur_df['pdb_id'] = pdb
                cur_df['entity_id'] = entity
                cur_df['tag'] = tag
                yield cur_df
        
        result = cls.to_data_frame(res)
        if len(result) == 0:
            return None
        try:
            sites = get_mutaRes(result)
        except IndexError:
            logging.error(result.to_string())
            return None
        return pd.concat(yield_dfrm_of_mutaRes(sites), sort=False, ignore_index=True)

    @classmethod
    def summary_res_conflict(cls, lyst, session=None):
        muta_query = '''
            MATCH (entry:Entry)-[:HAS_ENTITY]->(entity:Entity{POLYMER_TYPE:'P'})-[:HAS_RESIDUE_CONFLICT]->(resCon:ResidueConflict)
            WHERE entry.ID in $lyst
            RETURN entry.ID as pdb_id, entity.ID as entity_id, resCon.DETAILS, resCon.ID
        '''
        if session is None:
            if cls.session is not None:
                session = cls.session
        try:
            res = session.run(muta_query, lyst=list(lyst))
            return cls.res_conflict_then(res)
        except AttributeError:
            return muta_query, dict(lyst=list(lyst))

    @classmethod
    def deal_res_conflict(cls, res: pd.DataFrame):
        res.tag = res.tag.apply(lambda x: x.split('_')[0])
        group_cols = ['pdb_id', 'entity_id', 'tag']
        return pd.DataFrame(
            ((pdb_id, entity_id, tag, 
              json.dumps(to_interval(groupData.residue_number)).decode('utf-8')
              ) for (pdb_id, entity_id, tag), groupData in res.groupby(group_cols)),
            columns=group_cols+['conflict_range'])

    @staticmethod
    def yieldPureHo(entity_chain_map: Dict, unp_entity_info: Dict):
        for unp, entities in unp_entity_info.items():
            for entity in entities:
                chains = entity_chain_map[entity]
                if len(chains) > 1:
                    yield unp, tuple(combinations(chains, 2))

    @staticmethod
    def pass_range_check(left, right):
        print(
            jaccard.similarity(lyst2range(left), lyst2range(right)),
            overlap.similarity(lyst2range(left), lyst2range(right)))
        return True

    @classmethod
    def yieldDetectHo(cls, entity_chain_map: Dict, unp_entity_info: Dict):
        for unp, entities in unp_entity_info.items():
            if len(entities) > 1:
                for entity_a, entity_b in combinations(entities, 2):
                    if cls.pass_range_check(entities[entity_a], entities[entity_b]):
                        yield unp, tuple(product(entity_chain_map[entity_a], entity_chain_map[entity_b]))
    
    @staticmethod
    def yieldHe(entity_chain_map: Dict, unp_entity_info: Dict):
        for unp_a, unp_b in combinations(unp_entity_info, 2):
            # Check Taxonomy
            if unp_a[1] != unp_b[1]:
                continue
            for entity_a, entity_b in product(unp_entity_info[unp_a], unp_entity_info[unp_b]):
                assert entity_a != entity_b
                yield (unp_a, unp_b), tuple(product(entity_chain_map[entity_a], entity_chain_map[entity_b]))


class SeqPairwiseAlign(object):
    def __init__(self):
        self.seqa = None
        self.seqb = None
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.open_gap_score = -10
        self.aligner.extend_gap_score = -0.5
        self.aligner.substitution_matrix = matlist.blosum62
        self.alignment_count = 0

    @lru_cache()
    def makeAlignment(self, seqa, seqb):
        if seqa is None or seqb is None:
            return np.nan, np.nan
        self.seqa = seqa
        self.seqb = seqb
        alignments = self.aligner.align(seqa, seqb)
        for alignment in alignments:
            result = self.getAlignmentSegment(alignment)
            self.alignment_count += 1
            return json.dumps(result[0]).decode('utf-8'), json.dumps(result[1]).decode('utf-8')

    @staticmethod
    def getAlignmentSegment(alignment):
        segments1 = []
        segments2 = []
        i1, i2 = alignment.path[0]
        for node in alignment.path[1:]:
            j1, j2 = node
            if j1 > i1 and j2 > i2:
                segment1 = (i1 + 1, j1)
                segment2 = (i2 + 1, j2)
                segments1.append(segment1)
                segments2.append(segment2)
            i1, i2 = j1, j2
        return segments1, segments2


@lru_cache()
def convert_index(lrange_str: str, rrange_str: str, site: int):
    lrange = json.loads(lrange_str)
    rrange = json.loads(rrange_str)
    try:
        for (lstart, lend), (rstart, rend) in zip(lrange, rrange):
            assert (lend - lstart) == (rend-rstart), "convert_index(): Invalid range"
            if (site >= rstart) and (site <= rend):
                return int(site + lstart - rstart)
            else:
                continue
        return -1
    except AssertionError as e:
        logging.error(e)
        raise e


class Neo4j_API(Abclog):
    
    output_kwargs = {'sep': '\t', 'index': False, 'header': False, 'mode': 'a'}
    filtered_sifts_path = None
    res2pdbpath = None
    neo4j_api = None
    sqlite_api = None
    unpmap_filter:Dict = None
    sifts_filter:Dict = None
    entry_filter:Dict = None
    entry_info_all = ('method', 'ligand', 'nucleotides')
    entry_info_funcs = dict(zip(
        entry_info_all,
        (Entry.summary_method,Entry.summary_ligand,Entry.summary_nucleotides)))
    entry_info_add = {
        'has_hybrid_nucleotides': False,
        'nucleotides_entity_type': '',
        'BOUND_LIGAND_COUNT': 0,
        'BOUND_MOL_COUNT': 0
        }

    @classmethod
    @unsync
    async def pipe_new_sifts(cls, lyst, id_type):
        if len(lyst) > 0:
            # Assume that the elements of the lyst are all valid and are covered by the remote db
            query, kwargs = SIFTS.summary_mapping(lyst, id_type)
            res = await cls.neo4j_api.afetch(query, **kwargs)
            sifts_df = await SIFTS.deal_mapping(res, cls.neo4j_api)
            if sifts_df is not None:
                sifts_df = related_dataframe(cls.sifts_filter, sifts_df)
                if len(sifts_df) > 0:
                    # TODO: ADD Residue Conflict Info
                    query, kwargs = SIFTS.summary_res_conflict(
                        sifts_df.pdb_id.unique())
                    res = await cls.neo4j_api.afetch(query, **kwargs)
                    res = SIFTS.res_conflict_then(res)
                    if res is not None:
                        conflict_df = SIFTS.deal_res_conflict(res)
                        sifts_df = sifts_df.merge(conflict_df.rename(
                            columns={'tag': 'UniProt'}), how='left')
                    else:
                        sifts_df['conflict_range'] = None
                    # TODO: Export New Data [Sync]
                    values = sifts_df.to_dict('records')
                    await cls.sqlite_api.async_insert(cls.sqlite_api.SIFTS_Info, values)
                    return sifts_df
                else:
                    return None
            else:
                return None
        else:
            return None

    @classmethod
    @unsync
    async def pdb2unp(cls, lyst: Union[Iterable, Unfuture]):
        if not isinstance(lyst, Iterable):
            lyst = await lyst
        # Can be optimized
        sifts_df = await cls.pipe_new_sifts(lyst, 'pdb')
        return sifts_df

    @classmethod
    @unsync
    async def amapres2unp(cls, site_df: pd.DataFrame, unp:str, pdb_id: str, entity_id: str, chain_id: str, pdb_range: str, unp_range: str):
        sqlite_api = cls.sqlite_api
        sites = [convert_index(unp_range, pdb_range, x) for x in site_df.Pos]
        site_df['unp_index'] = sites
        res = await sqlite_api.PDBRes_Info.objects.filter(
            pdb_id=pdb_id, entity_id=entity_id,
            chain_id=chain_id, author_residue_number__in=site_df.Pos.unique()).all()
        pass

    @classmethod
    @unsync
    async def process_map2unp(cls, df, respath, siftspath, observed_only):
        df.pdb_id = df.pdb_id.str.lower()
        if 'entity_id' in df.columns:
            df.entity_id = df.entity_id.astype(str)
        query, kwargs = Entry.summary_eec(df.pdb_id.unique())
        res = await cls.neo4j_api.afetch(query, **kwargs)
        df = df.merge(Entry.to_data_frame(res), how="left")
        pdb_entity_struct_asym_id = df.pdb_id + '_' + df.entity_id + '_' + df.struct_asym_id
        query, kwargs = SIFTS.pdbchain2unpres(pdb_entity_struct_asym_id.unique(), observed_only)
        res = await cls.neo4j_api.afetch(query, **kwargs)
        res = SIFTS.to_data_frame(res)
        if len(res) > 0:
            res[['UniProt', 'unp_index']] = res.unp_index.str.split('_', expand=True)
            res.UNIQID = res.apply(lambda x: '_'.join(x['UNIQID'].split('_')[0:2]+[x['chain_id']]), axis=1)
            res.drop(columns=['chain_id'], inplace=True)
            res.UNIQID = res.UNIQID + '_' + res.UniProt
            res['sifts_check'] = 0
            sifts_df = await cls.pdb2unp(df.pdb_id.unique())
            sifts_df['UNIQID'] = sifts_df.pdb_id + '_' + sifts_df.entity_id.astype(str) + '_' + sifts_df.chain_id + '_' + sifts_df.UniProt
            focus_index = res[res.UNIQID.isin(set(res.UNIQID) & set(sifts_df.UNIQID))].index
            res.loc[focus_index, 'sifts_check'] = 1
            await pipe_out(res, respath)
            await pipe_out(sifts_df, siftspath)

    @classmethod
    @unsync
    async def pipe_new_entry_info(cls, pdbs):
        entry_info_func_res: Dict = {}
        for key, func in cls.entry_info_funcs.items():
            query, kwargs = func(pdbs)
            res = await cls.neo4j_api.afetch(query, **kwargs)
            entry_info_func_res[key] = Entry.to_data_frame(res)
        res = Entry.deal_nucleotides(entry_info_func_res.get('nucleotides', None))
        if res is not None:
            entry_info_func_res['nucleotides'] = res
        merge = partial(pd.merge, on=['pdb_id'], how='outer')
        pdb_entry_df = reduce(merge, (df for df in entry_info_func_res.values() if len(df) > 0))
        pdb_entry_df.rename(columns={'index': 'pdb_id'}, inplace=True)
        for col, default in cls.entry_info_add.items():
            if col not in pdb_entry_df:
                pdb_entry_df[col] = default
            else:
                pdb_entry_df[col].fillna(default, inplace=True)
        pdb_entry_df.BOUND_LIGAND_COUNT = pdb_entry_df.BOUND_LIGAND_COUNT.apply(int)
        pdb_entry_df.BOUND_MOL_COUNT = pdb_entry_df.BOUND_MOL_COUNT.apply(int)
        return pdb_entry_df

    @classmethod
    @unsync
    async def unp2pdb(cls, lyst: Union[Iterable, Unfuture]):
        '''
        Map from unp to pdb [DB]
        
        Principle
        ---------
        * Keep Consistency
        * Only Save What We Want [Filtered Data]
        '''
        if not isinstance(lyst, Iterable):
            lyst = await lyst  # .result()
        sqlite_api = cls.sqlite_api
        
        # TODO: Load Existing Data
        exists = await sqlite_api.SIFTS_Info.objects.filter(UniProt__in=lyst).all()
        if exists:
            e_sifts_df = pd.DataFrame(exists)
            lyst = list(set(lyst) - set(e_sifts_df.UniProt))
            # Get pdb lyst that already been saved in local db
            # e_pdbs = e_sifts_df.pdb_id.unique()
        else:
            e_sifts_df = None
        # TODO: Fetch New Data
        sifts_df = await cls.pipe_new_sifts(lyst, 'unp')
        sifts = [df for df in (sifts_df, e_sifts_df) if df is not None]
        if sifts:
            sifts_df = pd.concat(sifts, sort=False, ignore_index=True)
        else:
            return None
        # Following the principle of KEEP CONSISTENCY, we does not need the code below
        # e_sifts_df = related_dataframe(cls.sifts_filter, e_sifts_df)
        # Renew the pdb lyst. Now the elements in the lyst are all new to the local db
        # pdbs = list(set(pdbs)-set(e_pdbs))
        # TODO: Fetch Existing Data
        pdbs = sifts_df.pdb_id.unique()
        exists = await sqlite_api.Entry_Info.objects.filter(pdb_id__in=pdbs).all()
        if exists:
            e_pdb_entry_df = pd.DataFrame(exists)
            pdbs = list(set(pdbs) - set(e_pdb_entry_df.pdb_id))
        else:
            e_pdb_entry_df = None
        # TODO: Fetch New Data
        if len(pdbs) > 0:
            pdb_entry_df = await cls.pipe_new_entry_info(pdbs)
            pdb_entry_df = related_dataframe(cls.entry_filter, pdb_entry_df)
            # TODO: Export New Data [Sync]
            values = pdb_entry_df.to_dict('records')
            if values:
                await sqlite_api.async_insert(
                    sqlite_api.Entry_Info, values)
                # pdbs = pdb_entry_df.pdb_id.unique()
            else:
                pdb_entry_df = None
                # pdbs = []
            # del values
        else:
            pdb_entry_df = None
        
        pdb_entries = [df for df in (pdb_entry_df, e_pdb_entry_df) if df is not None]
        if pdb_entries:
            pdb_entry_df = pd.concat(pdb_entries, sort=False, ignore_index=True)
        else:
            return None
        
        # TODO: Load Existing Data
        pdbs = pdb_entry_df.pdb_id.unique()
        exists = await sqlite_api.SEQRES_Info.objects.filter(pdb_id__in=pdbs).all()
        if exists:
            e_seq_res_df = pd.DataFrame(exists)
            pdbs = list(set(pdbs) - set(e_seq_res_df.pdb_id))
        else:
            e_seq_res_df = None
        # TODO: Fetch New Data
        if len(pdbs) > 0:
            query, kwargs = Entry.summary_seq(pdbs)
            res = await cls.neo4j_api.afetch(query, **kwargs)
            seq_res_df = SIFTS.deal_seq_index(SIFTS.to_data_frame(res))
            # WARNING: NEED TO BE FIXED, SHOULD BE MORE FLEXIBLE
            del_pdbs = seq_res_df[
                seq_res_df.AVG_OBS_OBS_RATIO.le(0.25) & seq_res_df.SEQRES_COUNT.ge(50)
            ].pdb_id.unique()
            seq_res_df = seq_res_df[~seq_res_df.pdb_id.isin(del_pdbs)]
            values = seq_res_df.to_dict('records')
            # TODO: Export New Data [Sync]
            if values:
                await sqlite_api.async_insert(sqlite_api.SEQRES_Info, values)
            else:
                seq_res_df = None
        else:
            seq_res_df = None
        
        seq_ress = [df for df in (seq_res_df, e_seq_res_df) if df is not None]
        if seq_ress:
            seq_res_df = pd.concat(seq_ress, sort=False, ignore_index=True)
        else:
            return None
        
        return sifts_df[sifts_df.pdb_id.isin(seq_res_df.pdb_id.unique())]

    @classmethod
    @unsync
    async def amapres2pdb(cls, site_df: pd.DataFrame, unp:str, pdb_id: str, entity_id: str, chain_id: str, pdb_range: str, unp_range: str):
        sqlite_api = cls.sqlite_api
        sites = [convert_index(pdb_range, unp_range, x) for x in site_df.Pos]
        site_df['residue_number'] = sites
        sites = list(set(sites) | set({1}))
        res = await sqlite_api.PDBRes_Info.objects.filter(
            pdb_id=pdb_id, entity_id=entity_id,
            chain_id=chain_id, residue_number__in=sites).all()
        if not res:
            query, kwargs = Entry.get_residues(
                pdb_id, entity_id,
                chain_id)  # , sites
            res = await cls.neo4j_api.afetch(query, **kwargs)
            res = Entry.to_data_frame(res)
            if len(res) > 0:
                res.author_insertion_code.fillna('', inplace=True)
                await sqlite_api.async_insert(
                    sqlite_api.PDBRes_Info, res.to_dict('records'))
            else:
                return None
        else:
            res = pd.DataFrame(res)
        merge_df = pd.merge(site_df, res)
        merge_df['UniProt'] = unp
        return merge_df

    @classmethod
    @unsync
    async def mapres2pdb(cls, unp_map_df: pd.DataFrame, sifts_df: pd.DataFrame) -> Union[Unfuture, str, Path]:
        if unp_map_df is None or sifts_df is None:
            return
        unp_map_df.sort_values(by='UniProt', inplace=True)
        sifts_df.sort_values(by='UniProt', inplace=True)
        sqlite_api = cls.sqlite_api
        focus_cols = ['pdb_id', 'entity_id', 'chain_id',
                      'new_pdb_range', 'new_unp_range']
        sifts_nda = sifts_df[focus_cols].to_numpy()
        yourlist = unp_map_df.yourlist.to_numpy()
        unp_map_i_dict = slice_series(unp_map_df.UniProt.to_numpy())
        sifts_dict = slice_series(sifts_df.UniProt.to_numpy())        
        try:
            res_lyst = []
            for unp, (start, end) in sifts_dict.items():
                focus_sifts_nda = sifts_nda[start:end]
                start, end = unp_map_i_dict[unp]
                from_lyst = np.unique(yourlist[start:end])
                site_df = await sqlite_api.Site_Info.objects.filter(from_id__in=from_lyst).all()
                site_df = pd.DataFrame(site_df)
                for record in focus_sifts_nda:
                    res = await cls.amapres2pdb(site_df, unp, *record)
                    if res is not None and len(res) > 0:
                        res_lyst.append(res)
            if res_lyst:
                res = pd.concat(res_lyst, sort=False, ignore_index=True)
                cols = sorted(res.columns)
                res = res[cols]
                if cls.res2pdbpath.exists():
                    headers = None
                else:
                    headers = cols
                async with aiofiles.open(cls.res2pdbpath, 'a') as fileOb:
                    dataset = Dataset(headers=headers)
                    dataset.extend(res.to_records(index=False))
                    await fileOb.write(dataset.export('tsv'))
                return 1
            return 0
        except Exception as e:
            logging.error("res_map exception: "+str(e))
            traceback.print_exc()
            raise e

    @classmethod
    @unsync
    async def process_unp2pdb(cls, indata: Union[str, Path, pd.DataFrame, Unfuture], sep='\t') -> Union[Unfuture, Tuple]:
        '''
        Default chain process
        Default file seperator: tab
        '''
        # TODO: Add Converters
        if not isinstance(indata, (str, Path, pd.DataFrame)):
            indata = await indata  # .result()
        if isinstance(indata, pd.DataFrame):
            unp_map_df = indata
        elif indata is None:
            return None, None
        else:
            unp_map_df = pd.read_csv(indata, sep=sep)
        unp_map_df = related_dataframe(cls.unpmap_filter, unp_map_df)
        related_unp = unp_map_df.UniProt.unique()
        sifts_df = await cls.unp2pdb(related_unp)
        if sifts_df is not None:
            unp_map_df = unp_map_df[unp_map_df.UniProt.isin(sifts_df.UniProt.to_numpy())].reset_index(drop=True)
            sifts_df = sifts_df.reset_index(drop=True)
            if cls.filtered_sifts_path is not None:
                cols = sorted(sifts_df.columns)
                sifts_df = sifts_df[cols]
                if cls.filtered_sifts_path.exists():
                    headers = None
                else:
                    headers = cols
                async with aiofiles.open(cls.filtered_sifts_path, 'a') as fileOb:
                    dataset = Dataset(headers=headers)
                    dataset.extend(sifts_df.to_records(index=False))
                    await fileOb.write(dataset.export('tsv'))
                    # await fileOb.flush()
            return unp_map_df, sifts_df
        else:
            return None, None

    @classmethod
    @unsync
    async def process_mapres2pdb(cls, dfs: Union[Tuple[pd.DataFrame], Unfuture]) -> Union[Unfuture, str, Path]:
        if isinstance(dfs, Tuple):
            unp_map_df, sifts_df = dfs
        else:
            unp_map_df, sifts_df = await dfs  # .result()
        return cls.mapres2pdb(unp_map_df, sifts_df)

    @classmethod
    @unsync
    async def pipe_api_info(cls, func, **kwargs):
        query, newkwargs = func(**kwargs)
        res = await cls.neo4j_api.afetch(query, **newkwargs)
        return Entry.to_data_frame(res)

    @classmethod
    @unsync
    async def load_interact_chain_info(cls, i3d, pdbs):
        pre_cols = ['from_pisa_bond', 'from_arp_contact', 'from_i3d_chain_pair']
        pisa_df = await cls.pipe_api_info(Entry.summary_chain_interaction, pdbs=pdbs, source='pisa_bond')
        arp_df = await cls.pipe_api_info(Entry.summary_chain_interaction, pdbs=pdbs, source='arp_contact')
        cols = ['chain_id_1', 'chain_id_2']
        dfs = [df for df in (sort_sub_cols(pisa_df, cols), sort_sub_cols(arp_df, cols), i3d) if len(df) > 0]
        if len(dfs) > 0:
            merge = partial(pd.merge, how='outer')
            res = reduce(merge, dfs)
            for col in pre_cols:
                if col not in res.columns:
                    res[col] = False
                else:
                    res[col] = res[col].fillna(False)
            return res[['pdb_id']+cols+pre_cols]
        else:
            return None

    @staticmethod
    def pipe_sifts2interact(sifts_df):
        sifts_df['Entry'] = sifts_df.UniProt.str.split('-').str[0]
        # sifts_df = sifts_df[['Entry', 'UniProt', 'pdb_id', 'chain_id']]
        cols = [i for i in sifts_df.columns if i != 'pdb_id']
        sifts_df1 = sifts_df.rename(columns=dict(zip(cols, (i+'_1' for i in cols))))
        sifts_df2 = sifts_df.rename(columns=dict(zip(cols, (i+'_2' for i in cols))))
        sifts_df1 = sifts_df1.rename(columns={"Entry_1": "PROT1"})
        sifts_df2 = sifts_df2.rename(columns={"Entry_2": "PROT2"})
        sifts_df12 = pd.merge(sifts_df1, sifts_df2)
        return sifts_df12[sifts_df12.chain_id_1 != sifts_df12.chain_id_2].reset_index(drop=True)

    @classmethod
    @unsync
    async def pipe_interaction(cls, pdbs, i3d):
        int_df = await cls.load_interact_chain_info(i3d.chain_pairs_p(pdbs), pdbs)
        sifts_df = await cls.pipe_new_sifts(pdbs, 'pdb') # ! FILTERED & Can be optimized
        sifts_df = cls.pipe_sifts2interact(sifts_df)
        return int_df.merge(sifts_df, how='left').merge(i3d.unp_chain_pairs_p(pdbs), how='left').merge(i3d.all_unp_pairs_12, how='left')
"""


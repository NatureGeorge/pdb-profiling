# @Created Date: 2020-06-07 06:37:16 pm
# @Filename: score.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-06-07 06:37:18 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from __future__ import absolute_import
from itertools import combinations
from pandas import DataFrame, read_csv
import numpy as np
from typing import Dict, List, Iterable, Optional
from unsync import unsync
import asyncio
import orjson as json
from pdb_profiling.processors.pdbe.neo4j_api import SIFTS, slice_series, subtract_range, interval2set, range_len, outside_range_len, add_range, overlap_range
from pdb_profiling.processors.pdbe.sqlite_api import converters


def geometric_mean(array: Iterable):
    try:
        assert len(array[array<=0]) == 0
        array = np.log(array)
        return np.exp(array.sum()/len(array))
    except AssertionError:
        return 0


def harmonic_mean(array: Iterable):
    try:
        assert len(array[array==0]) == 0
        if not isinstance(array, np.ndarray):
            array = np.array(array)
        return len(array) / np.sum(1.0/array)
    except AssertionError:
        return 0


def minkowski_distance(array: Iterable, best: Optional[Iterable] = None, ord=None):
    '''
    Minkowski Distance

    * Hamming Distance: ord=1
    * Euclidean Distance: ord=2
    * Chebychev Distance: ord=np.inf
    * ...
    '''
    if not isinstance(array, np.ndarray):
        array = np.array(array)
    if best is None:
        best = np.array([1]*len(array))
    return np.linalg.norm(best-array, ord=ord)


max_ham_val = minkowski_distance([0, 0, 0, 0, 0, 0], ord=1)
max_eu_val = minkowski_distance([0, 0, 0, 0, 0, 0], ord=2)


class AHP(object):

    ELE_LIST = [
        # pdb_map_range - mis_index - non_index / UniPort Isoform Seq Len
        "COUNT(observed mapped standard residues)/LEN(UNP ISOFORM SEQ)",
        # pdb_map_range, seqres_count
        "1-COUNT(unmapped residues in C/N terminal)/LEN(SEQRES)",  # IF _ > 5 ELSE 0
        # ligand binding residues, seqres_count-mis_index
        "1-COUNT(ligand binding residues)/LEN(observed residues)",
        # pdb_map_range&mis_index, pdb_map_range
        "1-COUNT(missing mapped residues)/COUNT(mapped residues)",  # del Head Tail
        # (pdb_map_range - mis_index)&(non_index | mutated residues), pdb_map_range
        "1-COUNT(observed mapped modified | mutated residues)/COUNT(mapped residues)",
        # (insertion range, deletion range), seqres_count
        "1-MIN(COUNT(residues in insertion | deletion range)/LEN(SEQRES), 1)"]
    INI_LIST = [1, 1, 2, 3, 4, 1, 2, 3, 4, 2, 3, 4, 2, 3, 2]
    factor = -1
    
    def __init__(self, weight=None):
        self._weight = weight

    @property
    def weight(self) -> List:
        com_list = list(combinations(self.ELE_LIST, 2))
        com_dict = dict(zip(com_list, self.INI_LIST))
        mat = DataFrame([[1]*len(self.ELE_LIST)]*len(self.ELE_LIST))
        mat.columns = self.ELE_LIST
        mat.index = self.ELE_LIST
        for i in self.ELE_LIST:
            for j in self.ELE_LIST:
                mat[i][j] = com_dict.get((i, j), 1)
        x = mat.to_numpy()
        y = 1/x
        value, vector = np.linalg.eig(np.triu(x.T).T + np.triu(y.T) - np.diag(x.diagonal()))
        max_val = np.max(value)
        index = list(value).index(max_val)
        max_vector = vector[:, index]
        select_vector = np.real(max_vector/np.linalg.norm(max_vector, ord=1))
        return select_vector * self.factor
    
    def raw_score(self, array):
        if self._weight is not None:
            weight = self._weight
        else:
            weight = self.weight
            weight = weight/weight[0]
            weight[1:] = -weight[1:]
        return np.dot(array, weight)

    def score(self, array):
        '''
        Assume that the elements in the array are sorted
        '''
        return min(np.dot(array, self.weight), 1.0)


class Score_API(object):
    def __init__(self, sqlite, neo4j, outpath, logger, omit=5, add_score: bool=True):
        self.sqlite = sqlite
        self.neo4j = neo4j
        self.outpath = outpath
        self.omit = omit
        self.logger = logger
        self.add_score = add_score

    @unsync
    async def stat(self, focus_cols, record, unp, unp_len):
        kwargs = dict(zip(focus_cols[:3], record[:3]))
        seqres = await self.sqlite.SEQRES_Info.objects.filter(**kwargs).all()
        seqres = seqres[0]
        query, kwargs = SIFTS.get_binding_redidues(*record[:3])
        binding_res = await self.neo4j.afetch(query, **kwargs)
        binding_res = dict(next(iter(binding_res)))[
            'LIGAND_BINDING_RES_COUNT']
        pdb_range, unp_range, conflict_range, pdb_GAP_list, unp_GAP_list, var_list = record[3:]
        # None/NAN/'': conflict_range, seqres.MIS_INDEX, seqres.NON_INDEX

        # pdb_map_range - mis_index - non_index / UniPort Isoform Seq Len
        # COUNT(observed mapped standard residues)/LEN(UNP ISOFORM SEQ)
        obs_res_map_range: List = subtract_range(
            pdb_range, seqres.MIS_INDEX)
        obs_res_map_count = range_len(obs_res_map_range)
        obs_std_res_map_range: List = subtract_range(
            obs_res_map_range, seqres.NON_INDEX)
        obs_std_res_map_count = range_len(obs_std_res_map_range)
        v1 = obs_std_res_map_count/unp_len
        # "1-COUNT(unmapped residues in C/N terminal)/LEN(SEQRES)"
        out_res_count = outside_range_len(
            pdb_range, seqres.SEQRES_COUNT, omit=self.omit)
        v2 = 1 - out_res_count/seqres.SEQRES_COUNT
        # 1-COUNT(ligand binding residues)/LEN(observed residues)
        obs_res_range = subtract_range(
            [[1, seqres.SEQRES_COUNT]], seqres.MIS_INDEX)
        obs_res_count = range_len(obs_res_range)
        v3 = 1 - binding_res/obs_res_count
        # pdb_map_range&mis_index, pdb_map_range  del Head Tail
        # 1-COUNT(missing mapped residues)/COUNT(mapped residues)
        res_map_count = range_len(pdb_range)
        mis_res_map_count = res_map_count - obs_res_map_count
        v4 = 1 - mis_res_map_count/res_map_count
        # (pdb_map_range - mis_index)&(non_index | mutated residues), pdb_map_range
        # 1-COUNT(observed mapped modified | mutated residues)/COUNT(mapped residues)
        obs_non_muta_res_map_range = overlap_range(
            obs_res_map_range, add_range(seqres.NON_INDEX, conflict_range))
        obs_non_muta_res_map_count = range_len(
            obs_non_muta_res_map_range)
        v5 = 1 - obs_non_muta_res_map_count/res_map_count
        # (insertion range, deletion range), seqres_count
        # 1-MIN(COUNT(residues in insertion | deletion range)/LEN(SEQRES), 1)
        indel_count = \
            sum(json.loads(pdb_GAP_list)) + \
            sum(json.loads(unp_GAP_list)) + \
            sum(json.loads(var_list))
        v6 = 1 - min(indel_count/seqres.SEQRES_COUNT, 1)
        return (unp, record[0], record[1], record[2],
                      pdb_range, unp_range, seqres.SEQRES_COUNT,
                      unp_len, binding_res, obs_res_map_range, obs_res_map_count,
                      obs_std_res_map_range, obs_std_res_map_count,
                      out_res_count, obs_res_range, obs_res_count,
                      res_map_count, mis_res_map_count,
                      obs_non_muta_res_map_range, obs_non_muta_res_map_count,
                      indel_count,
                      v1, v2, v3, v4, v5, v6
                      )

    @unsync
    async def process(self, sifts_path):
        '''
        TODO: Make it real async
        '''
        sifts_df = read_csv(sifts_path, sep="\t", converters=converters)
        sifts_df.drop_duplicates(inplace=True)
        sifts_df.sort_values(by='UniProt', inplace=True)
        focus_cols = ['pdb_id', 'entity_id', 'chain_id',
                      'new_pdb_range', 'new_unp_range', 
                      'conflict_range', 'pdb_gap_list', 
                      'unp_gap_list', 'var_list']
        sifts_nda = sifts_df[focus_cols].to_numpy()
        sifts_dict = slice_series(sifts_df.UniProt.to_numpy())
        res_record_lyst = []
        for unp, (start, end) in sifts_dict.items():
            focus_sifts_nda = sifts_nda[start:end]
            query, kwargs = SIFTS.get_unp_seq_len(unp)
            unp_len = await self.neo4j.afetch(query, **kwargs)
            unp_len = dict(next(iter(unp_len)))['unp_len']
            res_record_lyst.extend(
                [await self.stat(focus_cols, record, unp, unp_len) for record in focus_sifts_nda])
        res_record_df = DataFrame(
            res_record_lyst,
            columns=('UniProt', 'pdb_id', 'entity_id', 'chain_id',
                     'new_pdb_range', 'new_unp_range', 'SEQRES_COUNT', 
                     'LEN(UniProt)',
                     'LIGAND_BINDING_RES_COUNT', 'obs_res_map_range',
                     'obs_res_map_count',
                     'obs_std_res_map_range', 'obs_std_res_map_count',
                     'out_res_count', 'obs_res_range', 'obs_res_count',
                     'res_map_count', 'mis_res_map_count',
                     'obs_non_muta_res_map_range', 'obs_non_muta_res_map_count',
                     'indel_count', 'v1', 'v2', 'v3', 'v4', 'v5', 'v6'
                     )
            )
        if self.add_score:
            res_record_df = self.score(res_record_df)
        res_record_df.to_csv(self.outpath, sep="\t", index=False)
        return self.outpath

    def score_pipe(self, array):
        return (self.ahp.score(array),
                geometric_mean(array),
                harmonic_mean(array),
                1-minkowski_distance(array, ord=1)/max_ham_val,
                1-minkowski_distance(array, ord=2)/max_eu_val,
                1-minkowski_distance(array, ord=np.inf))

    def score(self, dfrm):
        # arrays = dfrm[[f'v{i}' for i in range(1,7)]].to_numpy()
        self.ahp = AHP()
        dfrm['RAW_BS'] = dfrm.apply(lambda x: self.ahp.raw_score([
            x['obs_std_res_map_count'],
            x['out_res_count'],
            x['LIGAND_BINDING_RES_COUNT'],
            x['mis_res_map_count'],
            x['obs_non_muta_res_map_count'],
            x['indel_count']
        ])/x['LEN(UniProt)'],
            axis=1)
        dfrm[['BS', 'GM', 'HM', 'HS', 'ES', 'LS']] = dfrm[[f'v{i}' for i in range(
            1, 7)]].apply(lambda x: self.score_pipe(np.array([x[f'v{i}'] for i in range(1,7)])), axis=1, result_type='expand')
        return dfrm


if __name__ == "__main__":
    pass
    # '''
    array0 = [0, 0, 0, 0, 0, 0]
    array1 = [1, 1, 1, 1, 1, 1]
    array2 = [.5, .1, .1, .1, .1, .1]
    array3 = [.5, .5, .5, .5, .5, .5]
    array4 = [.3, 1, 1, 1, 1, 1]
    array5 = [0.402631579,	0.906735751,	1,	0.9,	1,	1]
    demo = AHP()
    
    print("DEMO ARRAY:")

    print(demo.score(array1), 
          geometric_mean(array1),
          harmonic_mean(array1),
          1-minkowski_distance(array1, ord=1)/max_ham_val,
          1-minkowski_distance(array1, ord=2)/max_eu_val,
          1-minkowski_distance(array1, ord=np.inf)
    )
    
    print(demo.score(array2),
          geometric_mean(array2),
          harmonic_mean(array2),
          1-minkowski_distance(array2, ord=1)/max_ham_val,
          1-minkowski_distance(array2, ord=2)/max_eu_val,
          1-minkowski_distance(array2, ord=np.inf)
    )

    print(demo.score(array3),
          geometric_mean(array3),
          harmonic_mean(array3),
          1-minkowski_distance(array3, ord=1)/max_ham_val,
          1-minkowski_distance(array3, ord=2)/max_eu_val,
          1-minkowski_distance(array3, ord=np.inf)
    )

    print(demo.score(array4),
          geometric_mean(array4),
          harmonic_mean(array4),
          1-minkowski_distance(array4, ord=1)/max_ham_val,
          1-minkowski_distance(array4, ord=2)/max_eu_val,
          1-minkowski_distance(array4, ord=np.inf)
    )

    print(demo.score(array5),
          geometric_mean(array5),
          harmonic_mean(array5),
          1-minkowski_distance(array5, ord=1)/max_ham_val,
          1-minkowski_distance(array5, ord=2)/max_eu_val,
          1-minkowski_distance(array5, ord=np.inf)
    )

    print("TEST ARRAY:")
    test_mat = np.eye(6)
    for cur_array in test_mat:
        print(cur_array)
        print(demo.score(cur_array),
              geometric_mean(cur_array),
              harmonic_mean(cur_array),
              1-minkowski_distance(cur_array, ord=1)/max_ham_val,
              1-minkowski_distance(cur_array, ord=2)/max_eu_val,
              1-minkowski_distance(cur_array, ord=np.inf)
              )
        print("++++")
    # '''

# @Created Date: 2021-03-10 11:36:31 am
# @Filename: cyrange.pyx
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2021-03-10 11:36:34 am
# @Copyright (c) 2021 MinghuiGroup, Soochow University
# cython: language_level=3
cimport cython
from operator import itemgetter
import orjson as json

def to_interval(object lyst):
    if isinstance(lyst, float):
        return tuple()
    if not isinstance(lyst, (set, frozenset)):
        lyst = frozenset(int(i) for i in lyst if i is not None)
    if len(lyst) == 0:
        return tuple()
    
    cdef:
        list start = []
        list interval_lyst = []
        int max_edge = max(lyst)
        int min_edge = min(lyst)
    
    if len(lyst) == (max_edge + 1 - min_edge):
        return ((min_edge, max_edge),)
    
    cdef:
        list lyst_list = sorted(lyst)
        int i, j, li
    
    for j in lyst_list:
        if len(start) == 0:
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


def lyst22interval(object x, object y):
    cdef:
        int i, j, start_x, start_y, index_x, index_y, pre_x, pre_y
        list interval_x, interval_y
        frozenset data = frozenset({item for item in zip(x,y)})
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
    return interval_x, interval_y


def lyst32interval(object x, object y, object z):
    cdef:
        int i, j, k, start_x, start_y, start_z, index_x, index_y, index_z, pre_x, pre_y, pre_z
        list interval_x, interval_y, interval_z
        frozenset data = frozenset({item for item in zip(x,y,z)})
    x, y, z = zip(*sorted(data, key=itemgetter(0)))
    start_x, start_y, start_z = x[0], y[0], z[0]
    index_x, index_y, index_z = x[0]-1, y[0]-1, z[0]-1
    interval_x, interval_y, interval_z = [], [], []
    for i, j, k in zip(x, y, z):
        pre_x = index_x + 1
        pre_y = index_y + 1
        pre_z = index_z + 1
        if pre_x == i and pre_y == j and pre_z == k:
            index_x, index_y, index_z = i, j, k
        else:
            interval_x.append((start_x, index_x))
            interval_y.append((start_y, index_y))
            interval_z.append((start_z, index_z))
            start_x, start_y, start_z = i, j, k
            index_x, index_y, index_z = i, j, k
    interval_x.append((start_x, index_x))
    interval_y.append((start_y, index_y))
    interval_z.append((start_z, index_z))
    return interval_x, interval_y, interval_z


cpdef int range_len(object lyst):
    if isinstance(lyst, float) or lyst is None:
        return 0
    lyst = convert_range(lyst)
    cdef:
        int length = 0
        int left, right
    for left, right in lyst:
        assert right >= left, f"\n{lyst}"
        length += right - left + 1
    return length


cpdef frozenset interval2set(object lyst):
    lyst = convert_range(lyst)
    cdef:
        frozenset range_set = frozenset()
        int left, right
    for left, right in lyst:
        range_set |= frozenset(range(left, right+1))
    return range_set


def subtract_range(object source_range, object target_range):
    if isinstance(target_range, float) or target_range is None:
        return source_range
    elif len(source_range) == 0:
        return tuple()
    elif len(target_range) == 0:
        return source_range
    cdef:
        frozenset source_range_set = interval2set(source_range)
        frozenset target_range_set = interval2set(target_range)
    return to_interval(source_range_set - target_range_set)


cdef bint check_range(object i):
    if isinstance(i, float) or (i is None) or (i == 'nan'):
        return False
    return True


def add_range(object left, object right):
    cdef:
        bint check_left = check_range(left)
        bint check_right = check_range(right)
    if check_left and not check_right:
        return left
    elif not check_left and check_right:
        return right
    elif not check_left and not check_right:
        return float('nan')
    cdef:
        frozenset left_range_set, right_range_set
    try:
        left_range_set = interval2set(left)
        right_range_set = interval2set(right)
        return to_interval(left_range_set | right_range_set)
    except Exception as e:
        print("{left}:{tleft}, {right}:{tright}".formath(left=left, tleft=type(left), right=right, tright=type(right)))
        raise e


cpdef list overlap_range(object obs_range, object unk_range):
    """
    * assumption: segment does not overlap with each other given a segment set
    * feature: overlap algorithm between two segment sets require no order in each segment set
    """
    if isinstance(unk_range, float) or unk_range is None:
        return []
    cdef:
        size_t index_obs, index_unk
        int[2] seg_obs, seg_unk
        int start1, start2, end1, end2, start, end
        bint sl, sr, el, er, s_in, e_in, ini, cov
        list ret = []
    obs_range = convert_range(obs_range)
    unk_range = convert_range(unk_range)
    cdef size_t unk_len = len(unk_range)
    
    for index_obs in range(len(obs_range)):
        seg_obs = obs_range[index_obs]
        start1 = seg_obs[0]
        end1 = seg_obs[1]
        for index_unk in range(unk_len):
            seg_unk = unk_range[index_unk]
            start2 = seg_unk[0]
            end2 = seg_unk[1]
            sl = start2 >= start1
            sr = start2 <= end1
            el = end2 >= start1
            er = end2 <= end1
            s_in = sl and sr
            e_in = el and er
            ini = s_in or e_in
            cov = (not sl) and (not er)
            start = start2 if s_in else start1
            end = end2 if e_in else end1
            if ini or cov:
                ret.append((start, end))
    return sorted(ret, key=itemgetter(0))


cpdef tuple overlap_range_2(object obs_range, object pdb_range, object unp_range):
    cdef:
        int start1, start2, start3, end1, end2, end3, start_a, start_b, end_a, end_b
        bint sl, sr, el, er, s_in, e_in, ini, cov
        list ret_a = []
        list ret_b = []
        size_t index_pdb_unp, index_obs
        int[2] seg_pdb, seg_unp, seg_obs
    obs_range = convert_range(obs_range)
    pdb_range = convert_range(pdb_range)
    unp_range = convert_range(unp_range)
    cdef size_t obs_len = len(obs_range)

    for index_pdb_unp in range(len(pdb_range)):
        seg_pdb = pdb_range[index_pdb_unp]
        seg_unp = unp_range[index_pdb_unp]
        start1 = seg_pdb[0]
        end1 = seg_pdb[1]
        start3 = seg_unp[0]
        end3 = seg_unp[1]
        for index_obs in range(obs_len):
            seg_obs = obs_range[index_obs]
            start2 = seg_obs[0]
            end2 = seg_obs[1]
            sl = start2 >= start1
            sr = start2 <= end1
            el = end2 >= start1
            er = end2 <= end1
            s_in = sl and sr
            e_in = el and er
            ini = s_in or e_in
            cov = (not sl) and (not er)
            if s_in:
                start_a = start2
                start_b = start2 - start1 + start3
            else:
                start_a = start1
                start_b = start3
            if e_in:
                end_a = end2
                end_b = end2 - end1 + end3
            else:
                end_a = end1
                end_b = end3
            if ini or cov:
                ret_a.append((start_a, end_a))
                ret_b.append((start_b, end_b))
    return ret_a, ret_b


cpdef tuple overlap_range_3(object unp_range_1, object unp_range_2, object pdb_range_1, object pdb_range_2):
    cdef:
        int start1, start2, start3, start4, end1, end2, end3, end4, start_a, start_b, start_c, end_a, end_b, end_c
        bint sl, sr, el, er, s_in, e_in, ini, cov
        list ret_unp = []
        list ret_pdb_1 = []
        list ret_pdb_2 = []
        size_t index_pdb_unp_1, index_pdb_unp_2
        int[2] seg_unp_1, seg_unp_2, seg_pdb_1, seg_pdb_2
    unp_range_1 = convert_range(unp_range_1)
    unp_range_2 = convert_range(unp_range_2)
    pdb_range_1 = convert_range(pdb_range_1)
    pdb_range_2 = convert_range(pdb_range_2)
    cdef size_t pdb_unp_len_1 = len(unp_range_1)
    cdef size_t pdb_unp_len_2 = len(unp_range_2)

    for index_pdb_unp_1 in range(pdb_unp_len_1):
        seg_unp_1 = unp_range_1[index_pdb_unp_1]
        seg_pdb_1 = pdb_range_1[index_pdb_unp_1]
        start1 = seg_unp_1[0]
        end1 = seg_unp_1[1]
        start3 = seg_pdb_1[0]
        end3 = seg_pdb_1[1]
        for index_pdb_unp_2 in range(pdb_unp_len_2):
            seg_unp_2 = unp_range_2[index_pdb_unp_2]
            seg_pdb_2 = pdb_range_2[index_pdb_unp_2]
            start2 = seg_unp_2[0]
            end2 = seg_unp_2[1]
            start4 = seg_pdb_2[0]
            end4 = seg_pdb_2[1]
            sl = start2 >= start1
            sr = start2 <= end1
            el = end2 >= start1
            er = end2 <= end1
            s_in = sl and sr
            e_in = el and er
            ini = s_in or e_in
            cov = (not sl) and (not er)
            if s_in:
                start_a = start2
                start_b = start2 - start1 + start3
                start_c = start4
            else:
                start_a = start1
                start_b = start3
                start_c = start1 - start2 + start4
            if e_in:
                end_a = end2
                end_b = end2 - end1 + end3
                end_c = end4
            else:
                end_a = end1
                end_b = end3
                end_c = end1 - end2 + end4
            if ini or cov:
                ret_unp.append((start_a, end_a))
                ret_pdb_1.append((start_b, end_b))
                ret_pdb_2.append((start_c, end_c))
    return ret_unp, ret_pdb_1, ret_pdb_2


cpdef tuple outside_range(object pdb_range, int seqres_len):
    pdb_range = convert_range(pdb_range)
    cdef:
        int out_head = pdb_range[0][0] - 1
        int out_tail = pdb_range[-1][-1] + 1
        tuple ret
        bint pass1 = 1 <= out_head
        bint pass2 = out_tail <= seqres_len
    if pass1 and pass2:
        ret = ((1, out_head), (out_tail, seqres_len))
    elif pass1 and not pass2:
        ret = ((1, out_head),)
    elif pass2 and not pass1:
        ret = ((out_tail, seqres_len),)
    else:
        ret = tuple()
    return ret


cpdef object convert_range(object input_range):
    cdef object ret
    if isinstance(input_range, str):
        ret = json.loads(input_range)
    else:
        ret = input_range
    return ret


cpdef bint isin_range(object input_range, int value):
    input_range = convert_range(input_range)
    cdef size_t size = len(input_range)
    cdef size_t i
    cdef int x, y
    for i in range(size):
        x, y = input_range[i]
        if value >= x and value <= y:
            return True
    return False


cdef int convert_index(object lrange, object rrange, int site) except *:
    # convert from rrange to lrange
    cdef:
        int lstart, rstart, lend, rend, index
        int[2] lseg, rseg
    for index in range(len(lrange)):
        lseg = lrange[index]
        rseg = rrange[index]
        lstart = lseg[0]
        lend = lseg[1]
        rstart = rseg[0]
        rend = rseg[1]
        assert lstart - lend == rstart - rend
        if (site >= rstart) and (site <= rend):
            return site + lstart - rstart
        else:
            continue
    return -1


cdef tuple new_tp_range(int start, int end):
    cdef:
        tuple r = (start, end)
        tuple ret = (r,)
    return ret



cpdef tuple trim_range(object obs_range, object pdb_range, object unp_range):
    obs_range = convert_range(obs_range)
    pdb_range = convert_range(pdb_range)
    unp_range = convert_range(unp_range)
    cdef:
        int pdb_head = pdb_range[0][0]
        int pdb_tail = pdb_range[-1][-1]
        bint pdb_head_obs = isin_range(obs_range, pdb_head)
        bint pdb_tail_obs = isin_range(obs_range, pdb_tail)
        list candidate_new_range, pdb_new_range, unp_new_range
        int pdb_obs_head, pdb_obs_tail, unp_obs_head, unp_obs_tail
        tuple tp_pdb_range, tp_unp_range

    if (not pdb_head_obs) or (not pdb_tail_obs):
        candidate_new_range = overlap_range(obs_range, pdb_range)
        if len(candidate_new_range) > 0:
            pdb_obs_head = candidate_new_range[0][0]
            pdb_obs_tail = candidate_new_range[-1][-1]
            unp_obs_head = convert_index(unp_range, pdb_range, pdb_obs_head)
            unp_obs_tail = convert_index(unp_range, pdb_range, pdb_obs_tail)
            tp_pdb_range = new_tp_range(pdb_obs_head, pdb_obs_tail)
            pdb_new_range = overlap_range(pdb_range, tp_pdb_range)
            tp_unp_range = new_tp_range(unp_obs_head, unp_obs_tail)
            unp_new_range = overlap_range(unp_range, tp_unp_range)
            # assert len(pdb_new_range) == len(unp_new_range), "{}, {}, {}\n{} {}".format(obs_range, pdb_range, unp_range, pdb_new_range, unp_new_range)
            return (pdb_new_range, unp_new_range)
    return (pdb_range, unp_range)


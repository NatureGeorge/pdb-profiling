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


cpdef int range_len(object lyst):
    if isinstance(lyst, float) or lyst is None:
        return 0
    elif isinstance(lyst, str):
        lyst = json.loads(lyst)
    cdef:
        int length = 0
        int left, right
    for left, right in lyst:
        assert right >= left, f"\n{lyst}"
        length += right - left + 1
    return length


cpdef frozenset interval2set(object lyst):
    if isinstance(lyst, str):
        lyst = json.loads(lyst)
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
    if isinstance(i, float) or (i is None) or (len(i) == 0) or (i == 'nan'):
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


def overlap_range(object obs_range, object unk_range):
    if isinstance(unk_range, float) or unk_range is None:
        return tuple()
    cdef:
        int start1, start2, end1, end2, start, end
        bint sl, sr, el, er, s_in, e_in, ini, cov
        list ret = []
    obs_range = json.loads(obs_range) if isinstance(obs_range, str) else obs_range
    unk_range = json.loads(unk_range) if isinstance(unk_range, str) else unk_range
    
    for start1, end1 in obs_range:
        for start2, end2 in unk_range:
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
    return ret


cpdef tuple outside_range(object pdb_range, int seqres_len):
    pdb_range = json.loads(pdb_range) if isinstance(pdb_range, str) else pdb_range
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



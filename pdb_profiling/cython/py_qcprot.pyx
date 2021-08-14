__author__ = "Bernhard Thiel"
__email__ = "thiel@tbi.univie.ac.at"
__license__ = "BSD 3-clause"


import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
cdef extern from "qcprot.h":
    double CalcRMSDRotationalMatrix(double **coords1, double **coords2, const int length, double *rot_out, const double *weight, int is_centered);

def rmsd(coords1, coords2, is_centered=False):
    if len(coords1) != len(coords2):
        raise ValueError("Both structures need to have the same number of coordinates")
    if np.all(coords1 == coords2):
        #qcprot would give nan after >50 iterations
        #This only covers the case of identical structures. If the centered structures are identical, nan is still returned (Should happen rarely)
        return 0.
    cdef np.ndarray[double, ndim=2, mode="fortran"] coords1_c = np.asarray(coords1, dtype = float, order="F")
    cdef np.ndarray[double, ndim=2, mode="fortran"] coords2_c = np.asarray(coords2, dtype = float, order="F")
    if coords1_c.shape[0] != coords2_c.shape[0] or coords1_c.shape[1] != coords2_c.shape[1] or coords1_c.shape[1] != 3: 
        raise ValueError("Both structures need to have the same shape (Nx3), found {}x{}, {}x{}".format(coords1_c.shape[0], coords1_c.shape[1], coords2_c.shape[0], coords2_c.shape[1]))
    cdef double[::1] rot_out = np.zeros((9))
    cdef double** point_to_coords1 = <double **>malloc(3 * sizeof(double*))
    cdef double** point_to_coords2 = <double **>malloc(3 * sizeof(double*))
    if not point_to_coords1 or not point_to_coords2:
        raise MemoryError
    try:
        for i in range(3):
            point_to_coords1[i] = &coords1_c[0,i]
        for i in range(3):
            point_to_coords2[i] = &coords2_c[0,i] 
        rmsd = CalcRMSDRotationalMatrix(&point_to_coords1[0], &point_to_coords2[0], len(coords1), &rot_out[0], NULL, is_centered)
    finally:
        free(point_to_coords1)
        free(point_to_coords2)
    return rmsd


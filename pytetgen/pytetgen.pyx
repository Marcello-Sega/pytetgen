#!/usr/bin/python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
# distutils: language=c++
# distutils: sources = pytetgen/tetgen.cxx pytetgen/predicates.cxx
cimport numpy as np
import numpy as np
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.string cimport memcpy


class Delaunay(object):
    """ Calculate a Delaunay triangulation in 3D from a set of points
    
        Parameters
        ----------
        points : ndarray of floats, shape (npoints, 3)
            Coordinates of points to triangulate


        References
        ----------
        [tetgen](http://wias-berlin.de/software/tetgen/)




        Examples
        --------

        Triangulation of a set of points:

        >>> import pytetgen
        >>> import numpy as np
        >>> N = 4
        >>> points = np.random.random(3*N).reshape(N,3)
        >>> tri = pytetgen.Delaunay(points)
        >>> np.sort(tri.simplices)
        array([[0, 1, 2, 3]], dtype=int32)


        Attributes
        ----------

        simplices   :    (ndarray of ints, shape (nsimplex, 4)) 
            Indices of the points forming the simplices in the triangulation.

        
    """
    def __init__(self,points):

        b = Tetgenbehavior() 
        b.quiet = True

        datain = Tetgenio()
        dataout= Tetgenio()

        datain.pointlist = points

        Tetrahedralize(b,datain,dataout,None,None)

        self.simplices = np.copy(dataout.tetrahedronlist)




cdef extern from "tetgen.h":

    cdef cppclass tetgenio:
        tetgenio() except+
        int numberofpoints, numberoftetrahedra, numberofcorners,mesh_dim
        double * pointlist
        int * tetrahedronlist

    cdef cppclass tetgenbehavior:
        tetgenbehavior() except+
        int quiet

    cdef void tetrahedralize(tetgenbehavior *b, tetgenio *data_in, tetgenio *data_out,tetgenio *addin, tetgenio *bgmin)


def Tetrahedralize(Tetgenbehavior behavior, Tetgenio data_in,Tetgenio data_out,Tetgenio addin, Tetgenio bgmin):
    tetrahedralize(&(behavior.c_tetgenbehavior),
                   &(data_in.c_tetgenio),
                   &(data_out.c_tetgenio),
                   &(addin.c_tetgenio),
                   &(bgmin.c_tetgenio))


cdef class Tetgenio:
    cdef tetgenio c_tetgenio      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.c_tetgenio = tetgenio()

    @property
    def numberofpoints(self):
        return self.c_tetgenio.numberofpoints

    @property
    def numberoftetrahedra(self):
        return self.c_tetgenio.numberoftetrahedra

    @property
    def mesh_dim(self):
        return self.c_tetgenio.mesh_dim



    @property
    def numberofcorners(self):
        return self.c_tetgenio.numberofcorners

    @property
    def pointlist(self):
        return np.asarray(<np.float64_t[:self.numberofpoints]> self.c_tetgenio.pointlist)

    @property
    def tetrahedronlist(self):
        cdef int nt = self.c_tetgenio.numberoftetrahedra
        cdef int nc = self.c_tetgenio.numberofcorners
        return np.asarray(<int[:nt*nc]> self.c_tetgenio.tetrahedronlist).reshape(nt,nc)

    @pointlist.setter
    def pointlist(self,val):
        # the numpy end 
        flatval = np.ascontiguousarray(val.flatten())
        cdef int npoints = flatval.shape[0]
        self.c_tetgenio.numberofpoints = npoints / 3
        cdef int size = 3*sizeof(double)*self.c_tetgenio.numberofpoints
        cdef np.float64_t[:] view = flatval
        try:
            PyMem_Free(self.c_tetgenio.pointlist)
        except:
            pass
        self.c_tetgenio.pointlist = <double*>PyMem_Malloc(size)
        memcpy(self.c_tetgenio.pointlist, <double*> (&view[0]), size)


cdef class Tetgenbehavior:
    cdef tetgenbehavior c_tetgenbehavior # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.c_tetgenbehavior = tetgenbehavior()

    @property
    def quiet(self):
        return self.c_tetgenbehavior.quiet
    @quiet.setter
    def quiet(self,val):
        self.c_tetgenbehavior.quiet=val 



#cdef tetrahedralize(behavior,data_in,data_out,addin,bgmin):


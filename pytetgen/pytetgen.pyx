# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
# distutils: language=c++
# distutils: sources=pytetgen/tetgen.cxx pytetgen/predicates.cxx

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


        Searching for neighboring tetrahedra (-1 means no neighbor).
        In the following example, the first tetrahedron has two neighbors,
        with index 1 and 2.
        
        >>> points = np.array([1,1,1, 1,-1,-1, -1,1,-1, -1,-1,1,-1,-1,-1],dtype=np.float64).reshape(5,3)
        >>> tri = pytetgen.Delaunay(points)
        >>> print tri.neighbors[0]
        [-1  2  1 -1]
    
        Searching for the tetrahedra containing a set of points. Two 
        are withing the triangulation, and the third is out

        >>> p = np.array([[0.,0.,0.],[0.,0,-0.5],[50.,0.,1.]])
        >>> print tri.find_simplex(p)
        [ 2  1 -1]

        Attributes
        ----------

        points    :   ndarray of double, shape (npoints, ndim)
            Coordinates of input points.

        simplices :   (ndarray of ints, shape (nsimplex, 4)) 
            Indices of the points forming the simplices in the triangulation.

        vertices : deprecated alias for simplices
    
        neighbors :   ndarray of ints, shape (nsimplex, ndim+1)
            Indices of neighbor simplices for each simplex.
            The kth neighbor is opposite to the kth vertex.
            For simplices at the boundary, -1 denotes no neighbor.

    """

    def __init__(self,points,neighbors=True):

        self._b = Tetgenbehavior() 
        self._b.quiet = True
        self._b.neighout = neighbors

        self._datain = Tetgenio()
        self._dataout= Tetgenio()
        self._m = Tetgenmesh()

        self._datain.pointlist = np.ascontiguousarray(points,dtype=np.float64)

        Tetrahedralize(self._b,self._datain,self._dataout,None,None,self._m)

    def find_simplex(self,xi,bruteforce=False,tol=None):
        """ Find the simplices containing the given points.

            Parameters
            ----------
            xi: ndarray of doubles, shape (npoints, 3)
                Coordinates of points to locate
            bruteforce: bool, optional
                Does nothing, for compatibility with scipy.spatial.Delaunay
            tol: float, optional
                Does nothing, for compatibility with scipy.spatial.Delaunay 

            Returns
            -------
            i : ndarray of int, same shape as xi
                Indices of simplices containing each point. Points outside the triangulation get the value -1.
        """    
        return self._m.locate(xi)

    @property
    def simplices(self):
        return self._dataout.tetrahedronlist

    @property
    def vertices(self):
        return self._dataout.tetrahedronlist

    @property
    def points(self):
        return self._datain.pointlist

    @property
    def neighbors(self):
        return self._dataout.neighborlist


cdef extern from "tetgen.h":

    cdef cppclass tetgenio:
        tetgenio() except+
        
        int numberofpoints, numberoftetrahedra, numberofcorners,mesh_dim
        double * pointlist
        int * tetrahedronlist
        int * neighborlist

    cdef cppclass tetgenbehavior:
        tetgenbehavior() except+
        int quiet
        int neighout

    cdef cppclass tetgenmesh:
        ctypedef double **tetrahedron 

        tetgenmesh() except+
        cppclass triface:
            triface() except+
            tetrahedron *tet
        int elemindex(tetrahedron* ptr)
        int locate(double *searchpt, tetgenmesh.triface* searchtet)


    cdef void tetrahedralize(tetgenbehavior *b, tetgenio *data_in, tetgenio *data_out,tetgenio *addin, tetgenio *bgmin, tetgenmesh *m)


def Tetrahedralize(Tetgenbehavior behavior, Tetgenio data_in,Tetgenio data_out,Tetgenio addin, Tetgenio bgmin, Tetgenmesh m):
    tetrahedralize(&(behavior.c_tetgenbehavior),
                   &(data_in.c_tetgenio),
                   &(data_out.c_tetgenio),
                   &(addin.c_tetgenio),
                   &(bgmin.c_tetgenio),
                   &(m.c_tetgenmesh))

cdef class Tetgenmesh:
    cdef tetgenmesh c_tetgenmesh

    cdef int UNKNOWN
    cdef int OUTSIDE
    cdef int INTETRAHEDRON
    cdef int ONFACE
    cdef int ONEDGE
    cdef int ONVERTEX
    cdef int ENCVERTEX
    cdef int ENCSEGMENT
    cdef int ENCSUBFACE
    cdef int NEARVERTEX
    cdef int NONREGULAR
    cdef int INSTAR
    cdef int BADELEMENT


    def __cinit__(self):
        self.c_tetgenmesh = tetgenmesh()
        self.UNKNOWN = 0
        self.OUTSIDE = 1
        self.INTETRAHEDRON = 2
        self.ONFACE = 3
        self.ONEDGE = 4 
        self.ONVERTEX = 5
        self.ENCVERTEX = 6 
        self.ENCSEGMENT = 7
        self.ENCSUBFACE = 8 
        self.NEARVERTEX = 9 
        self.NONREGULAR = 10
        self.INSTAR = 11
        self.BADELEMENT = 12 

    def locate(self,searchpt):
        cdef int npoints = searchpt.shape[0]
        val = np.ascontiguousarray(searchpt.flatten())
        res = np.empty(npoints,dtype=np.int32)
        cdef np.float64_t[:] val_v = val
        cdef np.int32_t[:] res_v = res
        cdef tetgenmesh.triface searchtet 
        searchtet.tet = NULL
        cdef int i,ret

        for i in range(npoints):
            ret = <int>self.c_tetgenmesh.locate(<double*> (&val_v[3*i]),&searchtet)
            if ret == self.OUTSIDE:
                res_v[i] = <int>-1
            else:
                res_v[i] = self.c_tetgenmesh.elemindex(searchtet.tet)
        return res

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
        return np.asarray(<np.float64_t[:3*self.numberofpoints]> self.c_tetgenio.pointlist).reshape(self.numberofpoints,3)

    @property
    def tetrahedronlist(self):
        cdef int nt = self.c_tetgenio.numberoftetrahedra
        cdef int nc = self.c_tetgenio.numberofcorners
        return np.asarray(<int[:nt*nc]> self.c_tetgenio.tetrahedronlist).reshape(nt,nc)

    @property
    def neighborlist(self):
        cdef int nt = self.c_tetgenio.numberoftetrahedra
        cdef int nc = 4
        return np.asarray(<int[:nt*nc]> self.c_tetgenio.neighborlist).reshape(nt,nc)


    @pointlist.setter
    def pointlist(self,val):
        # the numpy end 
        flatval = np.ascontiguousarray(val.flatten())
        cdef int npoints = flatval.shape[0]
        self.c_tetgenio.numberofpoints = npoints / <int>3
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

    @property
    def neighout(self):
        return self.c_tetgenbehavior.neighout
    @neighout.setter
    def neighout(self,val):
        self.c_tetgenbehavior.neighout=val 



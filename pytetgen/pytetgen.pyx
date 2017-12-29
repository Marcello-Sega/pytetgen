# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
# distutils: language=c++
# distutils: sources=pytetgen/tetgen.cxx pytetgen/predicates.cxx

cimport numpy as np
import numpy as np
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.string cimport memcpy

ctypedef struct voroedge:
    int v1
    int v2
    double vnormal[3]

ctypedef struct vorofacet:
    int c1
    int c2
    int *elist

class Voronoi(object):
    def __init__(self,points,_use_predicates=False,_quiet=True):

        self._ridge_dict = None
        self._regions = None
        self._b = Tetgenbehavior() 
        if _use_predicates == True:
            self._b.noexact = False
            self._b.nostaticfilter = True
        self._b.quiet = _quiet 
        self._b.nonodewritten = True
        self._b.weighted = False
        self._b.neighout = False

        self._b.voroout = True

        self._datain = Tetgenio()
        self._dataout = Tetgenio()
        self._m = Tetgenmesh()
        
        self._datain.pointlist = np.ascontiguousarray(points,dtype=np.float64)
        if points.shape[1]!=3:
            raise ValueError("Only 3d arrays are supported so far")
        Tetrahedralize(self._b,self._datain,self._dataout,None,None,self._m,points.shape[1])

    @property
    def nvertices(self):
        return self._dataout.numberofvpoints

    @property
    def vertices(self):
        return self._dataout.vpointlist

    @property
    def nridges(self):
        return self._dataout.numberofvedges

    # Note that in scipy.spatial.Voronoi (qhull) ridges in 3d are facets
    def _compute_ridge_dict(self):
        if self._ridge_dict is None:
            self._ridge_dict , self._ridge_list = self._dataout.ridge_dict()

    @property
    def regions(self):
        if self._regions is not None:
            return self._regions
        self._compute_ridge_dict()
        self._regions = self._dataout.regions(self._ridge_list)
        return self._regions

    @property
    def ridge_dict(self):
        self._compute_ridge_dict()
        return self._ridge_dict

    @property
    def ridge_points(self):
        self._compute_ridge_dict()
        return np.asarray(self._ridge_dict.keys())

    @property
    def ridge_vertices(self):
        self._compute_ridge_dict()
        return self._ridge_dict.values()


    @property
    def npoints(self):
        return self._datain.numberofpoints

    @property
    def points(self):
        return self._datain.pointlist


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

        Write the triangulation to an unstructured mesh in vtk format
        >>> tri.writevtk('tri.vtk')

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

        weights :   ndarray of double, shape (npoints, )
            Value of the weights associated to each point, or None. 
            If not None, a regular triangulation (weighted Delaunay) 
            is performed, instead of the plain Delaunay triangulation

    """

    def __init__(self,points,neighbors=True,weights=None,_use_predicates=False,_quiet=True):

        self._b = Tetgenbehavior() 
        if _use_predicates == True:
            self._b.noexact = False
            self._b.nostaticfilter = True
        self._b.quiet = _quiet 
        self._b.nonodewritten = True
        self._b.weighted = False
        self._b.neighout = neighbors

        self._datain = Tetgenio()
        self._dataout = Tetgenio()
        self._m = Tetgenmesh()

        self._datain.pointlist = np.ascontiguousarray(points,dtype=np.float64)
        if weights is not None:
            self._datain.numberofpointattributes = 1
            if len(weights) != self._datain.numberofpoints:
                raise ValueError("weights must be in the same number of points")
            self._datain.pointattributelist = np.ascontiguousarray(weights,dtype=np.float64)
            self._b.weighted = True
        else:
            self._b.weighted = False

        self._datain.pointlist = np.ascontiguousarray(points,dtype=np.float64)
        if points.shape[1]!=3:
            raise ValueError("Only 3d arrays are supported so far")
        Tetrahedralize(self._b,self._datain,self._dataout,None,None,self._m,points.shape[1])

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

    def writevtk(self,filename):
        self._m.Outmesh2vtk(filename.encode())
        """ Write the mesh to a vtk file

            Parameters
            ----------

            filename: string

        """

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
    def weights(self):
        if self._b.weighted == True:
            return self._datain.pointattributelist
        return None

    @property
    def neighbors(self):
        return self._dataout.neighborlist

    @property
    def npoints(self):
        return self._datain.numberofpoints

    @property
    def dim(self):
        return self._datain.mesh_dim


cdef extern from "tetgen.h":

    ctypedef double **tetrahedron 

    cdef cppclass tetgenio:
        tetgenio() except+

        ctypedef voroedge voroedge
        ctypedef vorofacet vorofacet

        int numberofpoints, numberoftetrahedra, numberofcorners,mesh_dim
        int numberofvpoints
        int numberofvedges
        int numberofvfacets
        int numberofvcells
        int numberofpointattributes
        voroedge * vedgelist
        vorofacet * vfacetlist
        int ** vcelllist
        double * pointlist
        double * vpointlist
        double * pointattributelist
        int * tetrahedronlist
        int * neighborlist

    cdef cppclass tetgenbehavior:
        tetgenbehavior() except+
        int quiet
        int nonodewritten
        int noelewritten
        int nofacewritten
        int noexact
        int nostaticfilter
        char * outfilename
        int weighted
        int neighout
        int voroout

    cdef cppclass tetgenmesh:

        tetgenmesh() except+
        cppclass triface:
            triface() except+
            tetrahedron *tet
        void outmesh2vtk(char* ofilename)
        int elemindex(tetrahedron* ptr)
        int locate(double *searchpt, tetgenmesh.triface* searchtet)


    cdef void tetrahedralize(tetgenbehavior *b, tetgenio *data_in, tetgenio *data_out,tetgenio *addin, tetgenio *bgmin, tetgenmesh *m)


def Tetrahedralize(Tetgenbehavior behavior, Tetgenio data_in,Tetgenio data_out, addin, bgmin, Tetgenmesh m, dim):
    tetrahedralize(&(behavior.c_tetgenbehavior),
		   &(data_in.c_tetgenio),
                   &(data_out.c_tetgenio),
                   NULL,NULL,&(m.c_tetgenmesh))


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
        if len(searchpt.shape) == 1: # a single point
            searchpt = searchpt.reshape((1,searchpt.shape[0]))
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

    def Outmesh2vtk(self, filename):
        if filename[-4:]=='.vtk':
            filename = filename[:-4]
        self.c_tetgenmesh.outmesh2vtk(filename)


cdef class Tetgenio:
    cdef tetgenio c_tetgenio      # hold a C++ instance which we're wrapping

    def __cinit__(self):
        self.c_tetgenio = tetgenio()

    @property
    def numberofpoints(self):
        return self.c_tetgenio.numberofpoints

    @property
    def numberofvpoints(self):
        return self.c_tetgenio.numberofvpoints

    @property
    def numberofvedges(self):
        return self.c_tetgenio.numberofvedges

    @property
    def numberofvfacets(self):
        return self.c_tetgenio.numberofvfacets

    @property
    def numberofvcells(self):
        return self.c_tetgenio.numberofvcells

    @property
    def numberofpointattributes(self):
        return self.c_tetgenio.numberofpointattributes

    @numberofpointattributes.setter
    def numberofpointattributes(self,val):
        self.c_tetgenio.numberofpointattributes = <int>val

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
    def vpointlist(self):
        return np.asarray(<np.float64_t[:3*self.numberofvpoints]> self.c_tetgenio.vpointlist).reshape(self.numberofvpoints,3)

    @property
    def vedgelist(self):
        return np.asarray(<voroedge[:self.numberofvedges]> self.c_tetgenio.vedgelist)


    def regions(self,ridge_list):
        cdef int i 
        cdef int j
        cdef int k
        cdef int nfacets
        cdef int facet

        _regions = []
        for i in range(self.numberofvcells):
            nfacets =  self.c_tetgenio.vcelllist[i][0]
            region_dict = {}
            for j in range(1,nfacets+1):
                facet =  self.c_tetgenio.vcelllist[i][j]
                if facet > -1:
                    for vertex in ridge_list[facet]:
                        region_dict[vertex]=''
                else:
                    region_dict[-1]=''
            _regions.append(region_dict.keys())
        return _regions

    def ridge_dict(self):
        cdef int i 
        cdef int j
        cdef int nridges
        _ridge_dict  = {}
        _ridge_list = []
        for i in range(self.numberofvfacets):
            vf = self.c_tetgenio.vfacetlist[i]
            nridges = vf.elist[0]
            L = {}
            for j in range(1,nridges+1):
                if vf.elist[j] > -1: # no bound check here, if elist[j]==-1, vedgelist will point nowhere
                    ve = self.c_tetgenio.vedgelist[vf.elist[j]]
                    L[ve.v1]=''
                    L[ve.v2]=''
                else:
                    L[-1]=''
            _ridge_dict[(vf.c1,vf.c2)] =  L.keys()
            _ridge_list.append(L.keys())
        return _ridge_dict,_ridge_list
        

    @property
    def pointattributelist(self):
        if self.c_tetgenio.pointattributelist == NULL:
            return None
        return np.asarray(<np.float64_t[:self.numberofpoints]> self.c_tetgenio.pointattributelist)

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

    @pointattributelist.setter
    def pointattributelist(self,val):
        flatval = np.ascontiguousarray(val.flatten())
        cdef int npoints = flatval.shape[0]
        cdef int size = sizeof(double)*npoints
        cdef np.float64_t[:] view = flatval
        try:
            PyMem_Free(self.c_tetgenio.pointattributelist)
        except:
            pass
        self.c_tetgenio.pointattributelist = <double*>PyMem_Malloc(size)
        memcpy(self.c_tetgenio.pointattributelist, <double*> (&view[0]), size)




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
    def nonodewritten(self):
        return self.c_tetgenbehavior.nonodewritten

    @nonodewritten.setter
    def nonodewritten(self,val):
        self.c_tetgenbehavior.nonodewritten=val 

    @property
    def noelewritten(self):
        return self.c_tetgenbehavior.noelewritten

    @noelewritten.setter
    def noelewritten(self,val):
        self.c_tetgenbehavior.noelewritten=val 

    @property
    def nofacewritten(self):
        return self.c_tetgenbehavior.nofacewritten

    @nofacewritten.setter
    def nofacewritten(self,val):
        self.c_tetgenbehavior.nofacewritten=val 

    @property
    def noexact(self):
        return self.c_tetgenbehavior.noexact

    @noexact.setter
    def noexact(self,val):
        self.c_tetgenbehavior.noexact=val 

    @property
    def nostaticfilter(self):
        return self.c_tetgenbehavior.nostaticfilter

    @nostaticfilter.setter
    def nostaticfilter(self,val):
        self.c_tetgenbehavior.noexact=val 


    @property
    def weighted(self):
        return self.c_tetgenbehavior.weighted

    @weighted.setter
    def weighted(self,val):
        self.c_tetgenbehavior.weighted=val 

    @property
    def neighout(self):
        return self.c_tetgenbehavior.neighout

    @neighout.setter
    def neighout(self,val):
        self.c_tetgenbehavior.neighout=val 

    @property
    def voroout(self):
        return self.c_tetgenbehavior.voroout

    @voroout.setter
    def voroout(self,val):
        self.c_tetgenbehavior.voroout=val 


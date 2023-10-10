"""
A simple test for pytetgen
"""
def test_simple():
    import pytetgen
    import numpy as np
    import scipy.spatial
    N=100
    np.random.seed(1)
    a=np.random.random(3*N).reshape(N,3)
    del1=pytetgen.Delaunay(a)
    del2=scipy.spatial.Delaunay(a)
    T1=np.sort(np.sort(del1.simplices,axis=1),axis=0)
    T2=np.sort(np.sort(del2.simplices,axis=1),axis=0)
    assert(np.all(T1==T2))

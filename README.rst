========
pytetgen
========

This is a python interface to tetgen, a powerful and fast mesh generator (http://wias-berlin.de/software/tetgen/)

This package includes the tetgen source, v.1.5, and provides (so far) minimal bindings to be ables to generate 
meshes in python. The basic interface follows that of scipy.spatial.Delaunay

tetgen (and this project) are distributed under the terms of the  GNU Affero General Public Licence (https://www.gnu.org/licenses/agpl-3.0.en.html). Quoting the preamble of the licence:

	"A secondary benefit of defending all users' freedom is that improvements made in alternate versions of the program, if they receive widespread use, become available for other developers to incorporate.[...]  However, in the case of software used on network servers, this result may fail to come about. The GNU General Public License permits making a modified version and letting the public access it on a server without ever releasing its source code to the public. The GNU Affero General Public License is designed specifically to ensure that, in such cases, the modified source code becomes available to the community. It requires the operator of a network server to provide the source code of the modified version running there to the users of that server. Therefore, public use of a modified version, on a publicly accessible server, gives the public access to the source code of the modified version."


References
==========
A technical paper about TetGen is available at 
Hang Si. 2015. "TetGen, a Delaunay-Based Quality Tetrahedral Mesh Generator". ACM Trans. on Mathematical Software. 41 (2), Article 11 (February 2015), 36 pages. DOI=10.1145/2629697 http://doi.acm.org/10.1145/2629697 

I guess the author would be glad if you cite him, if you use this wrapper. 



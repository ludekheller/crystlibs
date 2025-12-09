# CRYSTLIB - Quick Reference

**Module**: `crystlib.py`  
**Functions**: 76  
**Last Updated**: December 08, 2025

---

**B19p_B2_lattice_correspondence**
```python
def B19p_B2_lattice_correspondence(notation='Miyazaki'):
```
*Generate B19'→B2 lattice correspondence matrix for NiTi*

**B19p_B2_lattice_correspondence_ini**
```python
def B19p_B2_lattice_correspondence_ini():
```
*Initialize B19'→B2 correspondence for specific transformation variant*

**Rp_B2_lattice_correspondence**
```python
def Rp_B2_lattice_correspondence():
```
*Generate R-phase→B2 lattice correspondence for NiTi*

**an_between_vecs**
```python
def an_between_vecs(v1,v2,deg=True,full2pi=False):
```
*Calculate angle between two vectors*

**array2tuple**
```python
def array2tuple(arr, decimals=2):
```
*Convert a numpy array to a tuple with rounded elements*

**cubic2tetragonal_lattice_correspondence**
```python
def cubic2tetragonal_lattice_correspondence():
```
*Generate cubic→tetragonal lattice correspondence*

**cubic_lattice_vec**
```python
def cubic_lattice_vec(a):
```
*Generate cubic lattice vectors*

**def_gradient**
```python
def def_gradient(Cd,LA, LM,StressT=np.zeros((3,3)),STA=np.zeros((3,3,3,3)),STM=np.zeros((3,3,3,3)),CId=None):
```
*def_gradient - Crystallographic function for materials analysis*

**def_gradient_ini**
```python
def def_gradient_ini(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,StressT=np.zeros((3,3)),Parent_ST=np.zeros((3,3,3,3)),Product_ST=np.zeros((3,3,3,3))):
```
*def_gradient_ini - Crystallographic function for materials analysis*

**def_gradient_ini2**
```python
def def_gradient_ini2(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,StressT=np.zeros((3,3)),Parent_ST=np.zeros((3,3,3,3)),Product_ST=np.zeros((3,3,3,3))):
```
*def_gradient_ini2 - Crystallographic function for materials analysis*

**def_gradient_stressfree**
```python
def def_gradient_stressfree(Cd,LA, LM,CId=None):
```
*Calculate stress-free transformation strain (deformation gradient)*

**def_gradient_stressfree_ini**
```python
def def_gradient_stressfree_ini(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,Ci_d=None):
```
*Initialize stress-free deformation gradient from lattice parameters*

**dir2string**
```python
def dir2string(v, digits=2):
```
*Format direction as string representation [u,v,w] with Miller notation*

**equivalent_elements**
```python
def equivalent_elements(element,lattice):
```
*Find symmetrically equivalent elements*

**find_gcd**
```python
def find_gcd(x, y):
```
*Find greatest common divisor using recursive Euclidean algorithm*

**flipvector**
```python
def flipvector(v, Tol=1e-9):
```
*Flip vector to ensure positive first non-zero component*

**flipvector2negative**
```python
def flipvector2negative(v, Tol=1e-9):
```
*Flip vector to ensure negative first non-zero component*

**gen_twinned_lattice_points**
```python
def gen_twinned_lattice_points(ParentLatticePoints,eta1,shear_angle,K1,shift=0.0,dK1=None,bvr=None,deta1=None):
```
*Generate lattice points for parent and twinned regions*

**genallHexSys**
```python
def genallHexSys():
```
*Generate all hexagonal slip systems including <c+a> pyramidal*

**generate_hkls**
```python
def generate_hkls(hklmax, syms, hkls=[]):
```
*Generate unique Miller indices (hkl) considering crystal symmetry operations*

**generate_hkls01**
```python
def generate_hkls01(hklmax, syms, hkls=[]):
```
*Generate unique Miller indices (hkl) considering crystal symmetry operations*

**generate_hkls02**
```python
def generate_hkls02(hklmax, syms, G, hkls=[]):
```
*Generate unique Miller indices (hkl) with metric tensor transformation*

**generate_lattice_faces**
```python
def generate_lattice_faces(uvw2xyz,basal_dirs):
```
*Generate face polygons for unit cell visualization*

**generate_lattice_points**
```python
def generate_lattice_points(uvw2xyz,basal_dirs):
```
*Generate lattice points within specified unit cell range*

**generate_lattice_vectors**
```python
def generate_lattice_vectors(Points,uvw2xyz,S=1,Q=np.eye(3),xlim=[],ylim=[],zlim=[],fitpoints=False,shift=np.zeros(3)):
```
*generate_lattice_vectors - Crystallographic function for materials analysis*

**generate_lattite_atom_positions**
```python
def generate_lattite_atom_positions(atoms_xyz_position,uvw2xyz,S=1,Q=np.eye(3),R=np.eye(3),shift=np.zeros(3),xlim=[],ylim=[],zlim=[]):
```
*generate_lattite_atom_positions - Crystallographic function for materials ana...*

**generate_plane_vertices**
```python
def generate_plane_vertices(PlanePoints,normal,Q=np.eye(3),move=np.zeros(3)):
```
*generate_plane_vertices - Crystallographic function for materials analysis*

**generate_product_lattice_faces**
```python
def generate_product_lattice_faces(F,Parentlattices):
```
*Generate faces for two lattices*

**generate_product_lattice_points**
```python
def generate_product_lattice_points(F,Parentlattice_points,Q=np.eye(3)):
```
*generate_product_lattice_points - Crystallographic function for materials ana...*

**gensystemsHex**
```python
def gensystemsHex(eta1,K1,L, Lr,sm=None,eta2=None, K2=None):
```
*Generate hexagonal slip systems in real space*

**gensystemsHexIni**
```python
def gensystemsHexIni(eta1,K1,L, Lr,sm=None,eta2=None):
```
*Generate initial hexagonal slip system templates*

**get_interface2d**
```python
def get_interface2d(pointOutproj,normal,horizontalproj,verticalproj):
```
*Extract interface atoms from two lattices*

**get_twinning_dislocation**
```python
def get_twinning_dislocation(K1,eta1,eta2,L,G=None,Gr=None):
```
*Calculate twinning dislocation Burgers vector*

**get_twinning_plane_points**
```python
def get_twinning_plane_points(K1,Pointsout,horizontal,vertical):
```
*Get atoms on twinning plane*

**get_twinningdata**
```python
def get_twinningdata(orim,eus,Ldir_css,twin_systems,twt,phase, tension=True):
```
*Extract complete twinning data from parent and twin lattices*

**get_unique_families**
```python
def get_unique_families(hkls):
```
*Returns unique families of Miller indices based on permutation symmetry*

**habitplane_equation_solution**
```python
def habitplane_equation_solution(Uj,Ui,Qj,n,a,tol=1e-10):
```
*Solve habit plane equation for phase transformation*

**hkil2hkl**
```python
def hkil2hkl(hkil):
```
*Convert 4-index hexagonal plane (hkil) to 3-index (hkl)*

**hkl2hkil**
```python
def hkl2hkil(hkl):
```
*Convert 3-index hexagonal plane (hkl) to 4-index (hkil)*

**kronecker**
```python
def kronecker():
```
*Compute Kronecker delta δ_{ij}*

**lattice_correspondence**
```python
def lattice_correspondence(LatCorr,parent_symops,product_symops):
```
*Calculate lattice correspondence matrix between two crystal structures*

**lattice_vec**
```python
def lattice_vec(lattice_param):
```
*Calculate lattice vectors (a1, a2, a3) for various crystal systems*

**miller2fractional**
```python
def miller2fractional(uvw,frac=10,eps2=1e-2,decimals=5):
```
*Reduce Miller indices to lowest integer form*

**mohr_circles**
```python
def mohr_circles(tensor):
```
*Calculate Mohr's circles from strain or stress tensor*

**monoclinic_lattice_vec**
```python
def monoclinic_lattice_vec(a,b,c,beta):
```
*Generate monoclinic lattice vectors (B19' martensite)*

**niti_twinning**
```python
def niti_twinning(B2_symops,B2_recsymops,B19p_recsymops,B19p_symops,Uv,Parent_uvw2xyz,Parent_hkl2xyz,Product_uvw2xyz, Product_hkl2xyz,Parent_uvw_2_Product_uvw_all, Parent_hkl_2_Product_hkl_all, Parent_uvw_2_Product_uvw_all_norm,miller='greaterthanone',Qv=None):
```
*Get NiTi twinning parameters for specific variant*

**normArrayColumns**
```python
def normArrayColumns(arr):
```
*Normalize each column of matrix to unit length*

**np_kronecker**
```python
def np_kronecker():
```
*Compute Kronecker delta (NumPy version)*

**np_permut_tensor3**
```python
def np_permut_tensor3():
```
*Compute Levi-Civita symbol (NumPy-compatible version)*

**permut_tensor3**
```python
def permut_tensor3():
```
*Compute 3D Levi-Civita permutation symbol ε_{ijk}*

**perpendicular_vector**
```python
def perpendicular_vector(v):
```
*Find arbitrary unit vector perpendicular to input 3D vector*

**plane2string**
```python
def plane2string(v, digits=2):
```
*Format plane normal as string representation (h,k,l) with Miller notation*

**plane_line_intersection**
```python
def plane_line_intersection(n,V0,P0,P1):
```

**print_correspondence**
```python
def print_correspondence(Mcorr,VecA,latticeA, latticeB,planes=False,returnB=False):
```
*Print lattice correspondence matrix in readable format*

**read_txt**
```python
def read_txt(filename,delimiter='\t',skiprows=1):
```
*Read numerical data from text file*

**reciprocal_basis**
```python
def reciprocal_basis(a1,a2,a3):
```
*Calculate reciprocal lattice basis vectors from real space lattice vectors*

**select_atomic_plane**
```python
def select_atomic_plane(LatticePoints,normal,eps=1e-1,shift=0.,eps2=None):
```
*Select atoms lying on or near a crystallographic plane*

**select_atomic_region**
```python
def select_atomic_region(LatticePoints,normal,side='lower',eps=1e-1,shift=0.):
```
*Select atoms within spherical region*

**select_crystal_planes**
```python
def select_crystal_planes(NormalsOnCircle,ShearsOnCircle,maxhkl):
```
*Generate list of low-index crystallographic planes*

**select_plane**
```python
def select_plane(LatticeVectors,normal,eps=1e-1,shift=0.,Q=np.eye(3)):
```
*select_plane - Crystallographic function for materials analysis*

**strains_along_13mohrcirle**
```python
def strains_along_13mohrcirle(Strain,VV,normdiri,phi_around_V2,Parent_xyz2hkl):
```
*Calculate strains along maximum Mohr's circle*

**tetragonal_lattice_vec**
```python
def tetragonal_lattice_vec(a,b,c):
```
*Generate tetragonal lattice vectors*

**twin_equation_solution**
```python
def twin_equation_solution(Uj,Ui,L_A,Lr_A,L_M,Lr_M, R_AM, Ci_d,Ci_p,tol=1e-10,miller='greaterthanone',printlambda=False, Qj=None,Qi=None):
```
*Solve complete twin equation system*

**twin_equation_solution_ini**
```python
def twin_equation_solution_ini(Uj,Ui,Parent_uvw2xyz,Parent_hkl2xyz,Product_uvw2xyz,Product_hkl2xyz, Parent_uvw_2_Product_uvw_rot, Parent_uvw_2_Product_uvw,Parent_hkl_2_Product_hkl,tol=1e-10,miller='greaterthanone',printlambda=False, Qj=None,Qi=None):
```
*Initialize twin equation solver with default parameters*

**twinnedhabitplane**
```python
def twinnedhabitplane(Ui,Uj,Qij,a1,n1,hbplanes=[],addondata={},method='bhata'):
```
*Calculate twinned habit plane from lattice parameters*

**uvtw2uvw**
```python
def uvtw2uvw(uvtw):
```
*Convert 4-index hexagonal direction [uvtw] to 3-index [uvw]*

**uvw2uvtw**
```python
def uvw2uvtw(uvw):
```
*Convert 3-index hexagonal direction [uvw] to 4-index [uvtw]*

**vec2string**
```python
def vec2string(v, digits=2):
```
*Format vector as string representation [x,y,z] with specified precision*

**vector2miller_ini**
```python
def vector2miller_ini(v, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3):
```
*Convert Cartesian vector to Miller indices (initial guess)*

**vectors2miller**
```python
def vectors2miller(V, MIN=True, Tol=1e-9,tol=1e5,text=False):
```
*Convert multiple vectors to Miller indices*

**write_lattice_correspondence**
```python
def write_lattice_correspondence(ax,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,Product_lattice,Parent_lattice,FontSize=10,Fontweight="bold"):
```
*Write lattice correspondence analysis to file*

**write_mohr_planes**
```python
def write_mohr_planes(ax,Upperhalftext,Lowerhalftext,colors,markersize=8,markeredgewidth=2):
```
*Write Mohr circle analysis results to file*

**write_txt**
```python
def write_txt(filename,Header,DATA):
```
*Write numerical data to text file*

**xyz2fractional**
```python
def xyz2fractional(Txyz2uvw,V,frac=10,eps2=1e-2,decimals=5):
```
*Convert Cartesian coordinates to fractional Miller indices with reduction*

**xyz2fractional02**
```python
def xyz2fractional02(Txyz2uvw,V):
```
*Simple Cartesian to fractional coordinate transformation without reduction*

**zero_normal_strains**
```python
def zero_normal_strains(Strain, mcircles,VV,normdiri,phi_around_normdiri,Parent_xyz2hkl):
```
*Find planes with zero normal strain in given strain state*


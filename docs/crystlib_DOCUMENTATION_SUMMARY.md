# CRYSTLIB - Documentation Summary

**Module**: `crystlib.py`  
**Purpose**: Crystal Structure and Lattice Operations  
**Total Functions**: 76  
**Last Updated**: December 08, 2025

---

## Function Overview

### B19p_B2_lattice_correspondence

`def B19p_B2_lattice_correspondence(notation='Miyazaki'):`

Generate B19'→B2 lattice correspondence matrix for NiTi.

### B19p_B2_lattice_correspondence_ini

`def B19p_B2_lattice_correspondence_ini():`

Initialize B19'→B2 correspondence for specific transformation variant.

### Rp_B2_lattice_correspondence

`def Rp_B2_lattice_correspondence():`

Generate R-phase→B2 lattice correspondence for NiTi.

### an_between_vecs

`def an_between_vecs(v1,v2,deg=True,full2pi=False):`

Calculate angle between two vectors.

### array2tuple

`def array2tuple(arr, decimals=2):`

Convert a numpy array to a tuple with rounded elements.

### cubic2tetragonal_lattice_correspondence

`def cubic2tetragonal_lattice_correspondence():`

Generate cubic→tetragonal lattice correspondence.

### cubic_lattice_vec

`def cubic_lattice_vec(a):`

Generate cubic lattice vectors.

### def_gradient

`def def_gradient(Cd,LA, LM,StressT=np.zeros((3,3)),STA=np.zeros((3,3,3,3)),STM=np.zeros((3,3,3,3)),CId=None):`

def_gradient - Crystallographic function for materials analysis.

### def_gradient_ini

`def def_gradient_ini(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,StressT=np.zeros((3,3)),Parent_ST=np.zeros((3,3,3,3)),Product_ST=np.zeros((3,3,3,3))):`

def_gradient_ini - Crystallographic function for materials analysis.

### def_gradient_ini2

`def def_gradient_ini2(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,StressT=np.zeros((3,3)),Parent_ST=np.zeros((3,3,3,3)),Product_ST=np.zeros((3,3,3,3))):`

def_gradient_ini2 - Crystallographic function for materials analysis.

### def_gradient_stressfree

`def def_gradient_stressfree(Cd,LA, LM,CId=None):`

Calculate stress-free transformation strain (deformation gradient).

### def_gradient_stressfree_ini

`def def_gradient_stressfree_ini(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,Ci_d=None):`

Initialize stress-free deformation gradient from lattice parameters.

### dir2string

`def dir2string(v, digits=2):`

Format direction as string representation [u,v,w] with Miller notation.

### equivalent_elements

`def equivalent_elements(element,lattice):`

Find symmetrically equivalent elements.

### find_gcd

`def find_gcd(x, y):`

Find greatest common divisor using recursive Euclidean algorithm.

### flipvector

`def flipvector(v, Tol=1e-9):`

Flip vector to ensure positive first non-zero component.

### flipvector2negative

`def flipvector2negative(v, Tol=1e-9):`

Flip vector to ensure negative first non-zero component.

### gen_twinned_lattice_points

`def gen_twinned_lattice_points(ParentLatticePoints,eta1,shear_angle,K1,shift=0.0,dK1=None,bvr=None,deta1=None):`

Generate lattice points for parent and twinned regions.

### genallHexSys

`def genallHexSys():`

Generate all hexagonal slip systems including <c+a> pyramidal.

### generate_hkls

`def generate_hkls(hklmax, syms, hkls=[]):`

Generate unique Miller indices (hkl) considering crystal symmetry operations.

### generate_hkls01

`def generate_hkls01(hklmax, syms, hkls=[]):`

Generate unique Miller indices (hkl) considering crystal symmetry operations.

### generate_hkls02

`def generate_hkls02(hklmax, syms, G, hkls=[]):`

Generate unique Miller indices (hkl) with metric tensor transformation.

### generate_lattice_faces

`def generate_lattice_faces(uvw2xyz,basal_dirs):`

Generate face polygons for unit cell visualization.

### generate_lattice_points

`def generate_lattice_points(uvw2xyz,basal_dirs):`

Generate lattice points within specified unit cell range.

### generate_lattice_vectors

`def generate_lattice_vectors(Points,uvw2xyz,S=1,Q=np.eye(3),xlim=[],ylim=[],zlim=[],fitpoints=False,shift=np.zeros(3)):`

generate_lattice_vectors - Crystallographic function for materials analysis.

### generate_lattite_atom_positions

`def generate_lattite_atom_positions(atoms_xyz_position,uvw2xyz,S=1,Q=np.eye(3),R=np.eye(3),shift=np.zeros(3),xlim=[],ylim=[],zlim=[]):`

generate_lattite_atom_positions - Crystallographic function for materials analysis.

### generate_plane_vertices

`def generate_plane_vertices(PlanePoints,normal,Q=np.eye(3),move=np.zeros(3)):`

generate_plane_vertices - Crystallographic function for materials analysis.

### generate_product_lattice_faces

`def generate_product_lattice_faces(F,Parentlattices):`

Generate faces for two lattices.

### generate_product_lattice_points

`def generate_product_lattice_points(F,Parentlattice_points,Q=np.eye(3)):`

generate_product_lattice_points - Crystallographic function for materials analysis.

### gensystemsHex

`def gensystemsHex(eta1,K1,L, Lr,sm=None,eta2=None, K2=None):`

Generate hexagonal slip systems in real space.

### gensystemsHexIni

`def gensystemsHexIni(eta1,K1,L, Lr,sm=None,eta2=None):`

Generate initial hexagonal slip system templates.

### get_interface2d

`def get_interface2d(pointOutproj,normal,horizontalproj,verticalproj):`

Extract interface atoms from two lattices.

### get_twinning_dislocation

`def get_twinning_dislocation(K1,eta1,eta2,L,G=None,Gr=None):`

Calculate twinning dislocation Burgers vector.

### get_twinning_plane_points

`def get_twinning_plane_points(K1,Pointsout,horizontal,vertical):`

Get atoms on twinning plane.

### get_twinningdata

`def get_twinningdata(orim,eus,Ldir_css,twin_systems,twt,phase, tension=True):`

Extract complete twinning data from parent and twin lattices.

### get_unique_families

`def get_unique_families(hkls):`

Returns unique families of Miller indices based on permutation symmetry.

### habitplane_equation_solution

`def habitplane_equation_solution(Uj,Ui,Qj,n,a,tol=1e-10):`

Solve habit plane equation for phase transformation.

### hkil2hkl

`def hkil2hkl(hkil):`

Convert 4-index hexagonal plane (hkil) to 3-index (hkl).

### hkl2hkil

`def hkl2hkil(hkl):`

Convert 3-index hexagonal plane (hkl) to 4-index (hkil).

### kronecker

`def kronecker():`

Compute Kronecker delta δ_{ij}.

### lattice_correspondence

`def lattice_correspondence(LatCorr,parent_symops,product_symops):`

Calculate lattice correspondence matrix between two crystal structures.

### lattice_vec

`def lattice_vec(lattice_param):`

Calculate lattice vectors (a1, a2, a3) for various crystal systems.

### miller2fractional

`def miller2fractional(uvw,frac=10,eps2=1e-2,decimals=5):`

Reduce Miller indices to lowest integer form.

### mohr_circles

`def mohr_circles(tensor):`

Calculate Mohr's circles from strain or stress tensor.

### monoclinic_lattice_vec

`def monoclinic_lattice_vec(a,b,c,beta):`

Generate monoclinic lattice vectors (B19' martensite).

### niti_twinning

`def niti_twinning(B2_symops,B2_recsymops,B19p_recsymops,B19p_symops,Uv,Parent_uvw2xyz,Parent_hkl2xyz,Product_uvw2xyz, Product_hkl2xyz,Parent_uvw_2_Product_uvw_all, Parent_hkl_2_Product_hkl_all, Parent_uvw_2_Product_uvw_all_norm,miller='greaterthanone',Qv=None):`

Get NiTi twinning parameters for specific variant.

### normArrayColumns

`def normArrayColumns(arr):`

Normalize each column of matrix to unit length.

### np_kronecker

`def np_kronecker():`

Compute Kronecker delta (NumPy version).

### np_permut_tensor3

`def np_permut_tensor3():`

Compute Levi-Civita symbol (NumPy-compatible version).

### permut_tensor3

`def permut_tensor3():`

Compute 3D Levi-Civita permutation symbol ε_{ijk}.

### perpendicular_vector

`def perpendicular_vector(v):`

Find arbitrary unit vector perpendicular to input 3D vector.

### plane2string

`def plane2string(v, digits=2):`

Format plane normal as string representation (h,k,l) with Miller notation.

### plane_line_intersection

`def plane_line_intersection(n,V0,P0,P1):`



### print_correspondence

`def print_correspondence(Mcorr,VecA,latticeA, latticeB,planes=False,returnB=False):`

Print lattice correspondence matrix in readable format.

### read_txt

`def read_txt(filename,delimiter='\t',skiprows=1):`

Read numerical data from text file.

### reciprocal_basis

`def reciprocal_basis(a1,a2,a3):`

Calculate reciprocal lattice basis vectors from real space lattice vectors.

### select_atomic_plane

`def select_atomic_plane(LatticePoints,normal,eps=1e-1,shift=0.,eps2=None):`

Select atoms lying on or near a crystallographic plane.

### select_atomic_region

`def select_atomic_region(LatticePoints,normal,side='lower',eps=1e-1,shift=0.):`

Select atoms within spherical region.

### select_crystal_planes

`def select_crystal_planes(NormalsOnCircle,ShearsOnCircle,maxhkl):`

Generate list of low-index crystallographic planes.

### select_plane

`def select_plane(LatticeVectors,normal,eps=1e-1,shift=0.,Q=np.eye(3)):`

select_plane - Crystallographic function for materials analysis.

### strains_along_13mohrcirle

`def strains_along_13mohrcirle(Strain,VV,normdiri,phi_around_V2,Parent_xyz2hkl):`

Calculate strains along maximum Mohr's circle.

### tetragonal_lattice_vec

`def tetragonal_lattice_vec(a,b,c):`

Generate tetragonal lattice vectors.

### twin_equation_solution

`def twin_equation_solution(Uj,Ui,L_A,Lr_A,L_M,Lr_M, R_AM, Ci_d,Ci_p,tol=1e-10,miller='greaterthanone',printlambda=False, Qj=None,Qi=None):`

Solve complete twin equation system.

### twin_equation_solution_ini

`def twin_equation_solution_ini(Uj,Ui,Parent_uvw2xyz,Parent_hkl2xyz,Product_uvw2xyz,Product_hkl2xyz, Parent_uvw_2_Product_uvw_rot, Parent_uvw_2_Product_uvw,Parent_hkl_2_Product_hkl,tol=1e-10,miller='greaterthanone',printlambda=False, Qj=None,Qi=None):`

Initialize twin equation solver with default parameters.

### twinnedhabitplane

`def twinnedhabitplane(Ui,Uj,Qij,a1,n1,hbplanes=[],addondata={},method='bhata'):`

Calculate twinned habit plane from lattice parameters.

### uvtw2uvw

`def uvtw2uvw(uvtw):`

Convert 4-index hexagonal direction [uvtw] to 3-index [uvw].

### uvw2uvtw

`def uvw2uvtw(uvw):`

Convert 3-index hexagonal direction [uvw] to 4-index [uvtw].

### vec2string

`def vec2string(v, digits=2):`

Format vector as string representation [x,y,z] with specified precision.

### vector2miller_ini

`def vector2miller_ini(v, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3):`

Convert Cartesian vector to Miller indices (initial guess).

### vectors2miller

`def vectors2miller(V, MIN=True, Tol=1e-9,tol=1e5,text=False):`

Convert multiple vectors to Miller indices.

### write_lattice_correspondence

`def write_lattice_correspondence(ax,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,Product_lattice,Parent_lattice,FontSize=10,Fontweight="bold"):`

Write lattice correspondence analysis to file.

### write_mohr_planes

`def write_mohr_planes(ax,Upperhalftext,Lowerhalftext,colors,markersize=8,markeredgewidth=2):`

Write Mohr circle analysis results to file.

### write_txt

`def write_txt(filename,Header,DATA):`

Write numerical data to text file.

### xyz2fractional

`def xyz2fractional(Txyz2uvw,V,frac=10,eps2=1e-2,decimals=5):`

Convert Cartesian coordinates to fractional Miller indices with reduction.

### xyz2fractional02

`def xyz2fractional02(Txyz2uvw,V):`

Simple Cartesian to fractional coordinate transformation without reduction.

### zero_normal_strains

`def zero_normal_strains(Strain, mcircles,VV,normdiri,phi_around_normdiri,Parent_xyz2hkl):`

Find planes with zero normal strain in given strain state.


---

**Total**: 76 functions

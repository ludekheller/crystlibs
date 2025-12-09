# CRYSTLIB - Comprehensive Quick Reference

**Module**: `crystlib.py`  
**Functions**: 76  
**Last Updated**: December 08, 2025

This reference provides function signatures with brief descriptions and example usage patterns.

---

## B19p_B2_lattice_correspondence

```python
def B19p_B2_lattice_correspondence(notation='Miyazaki'):
```

Generate B19'→B2 lattice correspondence matrix for NiTi.
    Returns the transformation matrix relating B19' monoclinic martensite
    to B2 cubic austenite lattice vectors. Fundamental for NiTi shape

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Get correspondence matrix
        >>> C = B19p_B2_lattice_correspondence()
# ...
```

---

## B19p_B2_lattice_correspondence_ini

```python
def B19p_B2_lattice_correspondence_ini():
```

Initialize B19'→B2 correspondence for specific transformation variant.
    Returns correspondence matrix for one of the crystallographic variants
    of the B19' martensite transformation in NiTi.

**Example**:
```python
>>> # Get correspondence for variant 1
        >>> C1 = B19p_B2_lattice_correspondence_ini(variant=1)
        >>> print("Variant 1 correspondence:")
        >>> print(C1)
# ...
```

---

## Rp_B2_lattice_correspondence

```python
def Rp_B2_lattice_correspondence():
```

Generate R-phase→B2 lattice correspondence for NiTi.
    Returns correspondence matrix for R-phase (rhombohedral) to B2
    transformation in NiTi alloys. R-phase is intermediate phase.

**Example**:
```python
>>> C_R_B2 = Rp_B2_lattice_correspondence()
        >>> print("R-phase → B2 correspondence:")
        >>> print(C_R_B2)
```

---

## an_between_vecs

```python
def an_between_vecs(v1,v2,deg=True,full2pi=False):
```

Calculate angle between two vectors.
    Computes angle using dot product formula.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # 90° angle
        >>> v1 = np.array([1, 0, 0])
# ...
```

---

## array2tuple

```python
def array2tuple(arr, decimals=2):
```

Convert a numpy array to a tuple with rounded elements.

**Example**:
```python
>>> arr = np.array([1.234, 2.567, 3.891])
        >>> result = array2tuple(arr, decimals=2)
        >>> print(result)
        (1.23, 2.57, 3.89)
# ...
```

---

## cubic2tetragonal_lattice_correspondence

```python
def cubic2tetragonal_lattice_correspondence():
```

Generate cubic→tetragonal lattice correspondence.
    Returns correspondence matrix for cubic to tetragonal phase
    transformation (e.g., FCC→FCT).

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> C = cubic2tetragonal_lattice_correspondence()
        >>> 
# ...
```

---

## cubic_lattice_vec

```python
def cubic_lattice_vec(a):
```

Generate cubic lattice vectors.
    Creates 3×3 lattice matrix for cubic crystal system with parameter a.
    All angles are 90° and all lengths are equal (a=b=c).

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # B2 austenite (NiTi)
        >>> LA = cubic_lattice_vec(3.015)
# ...
```

---

## def_gradient

```python
def def_gradient(Cd,LA, LM,StressT=np.zeros((3,3)),STA=np.zeros((3,3,3,3)),STM=np.zeros((3,3,3,3)),CId=None):
```

def_gradient - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## def_gradient_ini

```python
def def_gradient_ini(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,StressT=np.zeros((3,3)),Parent_ST=np.zeros((3,3,3,3)),Product_ST=np.zeros((3,3,3,3))):
```

def_gradient_ini - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## def_gradient_ini2

```python
def def_gradient_ini2(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,StressT=np.zeros((3,3)),Parent_ST=np.zeros((3,3,3,3)),Product_ST=np.zeros((3,3,3,3))):
```

def_gradient_ini2 - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## def_gradient_stressfree

```python
def def_gradient_stressfree(Cd,LA, LM,CId=None):
```

Calculate stress-free transformation strain (deformation gradient).
    Computes F₀ = L2 · L1⁻¹ representing lattice transformation
    without applied stress.

**Example**:
```python
>>> # B2 → B19' transformation in NiTi
        >>> L_B2 = cubic_lattice_vec(3.015)
        >>> L_B19p = monoclinic_lattice_vec(2.889, 4.120, 4.622, 96.8)
        >>> 
# ...
```

---

## def_gradient_stressfree_ini

```python
def def_gradient_stressfree_ini(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,Ci_d=None):
```

Initialize stress-free deformation gradient from lattice parameters.
    Convenience function to calculate F₀ from lattice parameters directly.

**Example**:
```python
>>> # Cubic → Tetragonal
        >>> params_cubic = (3.0, 3.0, 3.0, 90, 90, 90)
        >>> params_tetra = (3.0, 3.0, 4.0, 90, 90, 90)
        >>> 
# ...
```

---

## dir2string

```python
def dir2string(v, digits=2):
```

Format direction as string representation [u,v,w] with Miller notation.
    Converts direction vector to crystallographic notation using
    square brackets, following standard Miller index convention.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Format [110] direction
        >>> uvw = [1, 1, 0]
# ...
```

---

## equivalent_elements

```python
def equivalent_elements(element,lattice):
```

Find symmetrically equivalent elements.

---

## find_gcd

```python
def find_gcd(x, y):
```

Find greatest common divisor using recursive Euclidean algorithm.
    Fundamental mathematical operation used throughout the module for
    Miller indices reduction and fractional coordinate normalization.

**Example**:
```python
>>> # Basic GCD calculation
        >>> gcd = find_gcd(48, 18)
        >>> print(gcd)  # Output: 6
        >>> 
# ...
```

---

## flipvector

```python
def flipvector(v, Tol=1e-9):
```

Flip vector to ensure positive first non-zero component.
    Standardizes vector direction for comparison.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> v1 = np.array([-1, 2, 3])
        >>> v1_flipped = flipvector(v1)
# ...
```

---

## flipvector2negative

```python
def flipvector2negative(v, Tol=1e-9):
```

Flip vector to ensure negative first non-zero component.
    Opposite of flipvector() - ensures negative first component.

**Example**:
```python
>>> v = np.array([1, 2, 3])
        >>> v_neg = flipvector2negative(v)
        >>> print(v_neg)  # [-1, -2, -3]
```

---

## gen_twinned_lattice_points

```python
def gen_twinned_lattice_points(ParentLatticePoints,eta1,shear_angle,K1,shift=0.0,dK1=None,bvr=None,deta1=None):
```

Generate lattice points for parent and twinned regions.
    Creates atomic positions showing twin boundary.

**Example**:
```python
>>> L = cubic_lattice_vec(3.0)
        >>> K1 = np.array([1, 1, 1])  # (111) twin
        >>> 
        >>> parent_pts, twin_pts = gen_twinned_lattice_points(L, K1, 10, 10, 10)
# ...
```

---

## genallHexSys

```python
def genallHexSys():
```

Generate all hexagonal slip systems including <c+a> pyramidal.
    Comprehensive slip system generation for HCP crystals including
    basal, prismatic, and pyramidal <c+a> systems. Essential for

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Magnesium lattice
        >>> a = 3.21
# ...
```

---

## generate_hkls

```python
def generate_hkls(hklmax, syms, hkls=[]):
```

Generate unique Miller indices (hkl) considering crystal symmetry operations.
    This version rounds values to avoid floating point comparison issues.

**Example**:
```python
>>> import numpy as np
        >>> # Define cubic symmetry operations (identity and 90° rotation around z)
        >>> identity = np.eye(3)
        >>> rot_z = np.array([[0, -1, 0],
# ...
```

---

## generate_hkls01

```python
def generate_hkls01(hklmax, syms, hkls=[]):
```

Generate unique Miller indices (hkl) considering crystal symmetry operations.
    Alternative version without rounding.

**Example**:
```python
>>> import numpy as np
        >>> # Define tetragonal symmetry operations
        >>> identity = np.eye(3)
        >>> rot_90 = np.array([[0, -1, 0],
# ...
```

---

## generate_hkls02

```python
def generate_hkls02(hklmax, syms, G, hkls=[]):
```

Generate unique Miller indices (hkl) with metric tensor transformation.
    This version applies symmetry operations in the correct metric space.

**Example**:
```python
>>> import numpy as np
        >>> # Define hexagonal metric tensor
        >>> a = 3.0  # lattice parameter
        >>> c = 5.0  # c-axis parameter
# ...
```

---

## generate_lattice_faces

```python
def generate_lattice_faces(uvw2xyz,basal_dirs):
```

Generate face polygons for unit cell visualization.
    Creates vertex lists for 6 faces of unit cell for 3D rendering.

**Example**:
```python
>>> import numpy as np
        >>> from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        >>> 
        >>> L = cubic_lattice_vec(3.0)
# ...
```

---

## generate_lattice_points

```python
def generate_lattice_points(uvw2xyz,basal_dirs):
```

Generate lattice points within specified unit cell range.
    Creates array of lattice points for visualization and analysis.
    Generates n1×n2×n3 unit cells.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Generate 2×2×2 cubic lattice
        >>> L_cubic = cubic_lattice_vec(3.0)
# ...
```

---

## generate_lattice_vectors

```python
def generate_lattice_vectors(Points,uvw2xyz,S=1,Q=np.eye(3),xlim=[],ylim=[],zlim=[],fitpoints=False,shift=np.zeros(3)):
```

generate_lattice_vectors - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## generate_lattite_atom_positions

```python
def generate_lattite_atom_positions(atoms_xyz_position,uvw2xyz,S=1,Q=np.eye(3),R=np.eye(3),shift=np.zeros(3),xlim=[],ylim=[],zlim=[]):
```

generate_lattite_atom_positions - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## generate_plane_vertices

```python
def generate_plane_vertices(PlanePoints,normal,Q=np.eye(3),move=np.zeros(3)):
```

generate_plane_vertices - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## generate_product_lattice_faces

```python
def generate_product_lattice_faces(F,Parentlattices):
```

Generate faces for two lattices.
    Creates face polygons for both unit cells.

**Example**:
```python
>>> L1 = cubic_lattice_vec(3.0)
        >>> L2 = tetragonal_lattice_vec(3.0, 3.0, 4.0)
        >>> faces1, faces2 = generate_product_lattice_faces(L1, L2)
```

---

## generate_product_lattice_points

```python
def generate_product_lattice_points(F,Parentlattice_points,Q=np.eye(3)):
```

generate_product_lattice_points - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## gensystemsHex

```python
def gensystemsHex(eta1,K1,L, Lr,sm=None,eta2=None, K2=None):
```

Generate hexagonal slip systems in real space.
    Transforms slip system templates to actual Cartesian coordinates
    using the hexagonal lattice matrix. Produces slip directions and

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # HCP titanium lattice
        >>> a = 2.95
# ...
```

---

## gensystemsHexIni

```python
def gensystemsHexIni(eta1,K1,L, Lr,sm=None,eta2=None):
```

Generate initial hexagonal slip system templates.
    Creates basal, prismatic, and pyramidal slip system templates for
    hexagonal close-packed (HCP) crystal structures. Returns normalized

**Example**:
```python
>>> dirs, norms = gensystemsHexIni()
        >>> print(f"Number of slip systems: {len(dirs)}")
        >>> 
        >>> # Examine first system
# ...
```

---

## get_interface2d

```python
def get_interface2d(pointOutproj,normal,horizontalproj,verticalproj):
```

Extract interface atoms from two lattices.
    Finds atoms from both structures near interface plane.
    Used in phase boundary visualization.

**Example**:
```python
>>> # Parent and product phase atoms
        >>> L1 = cubic_lattice_vec(3.0)
        >>> L2 = tetragonal_lattice_vec(3.0, 3.0, 4.0)
        >>> atoms1 = generate_lattite_atom_positions(L1, n1=5, n2=5, n3=5)
# ...
```

---

## get_twinning_dislocation

```python
def get_twinning_dislocation(K1,eta1,eta2,L,G=None,Gr=None):
```

Calculate twinning dislocation Burgers vector.
    Determines dislocation content at twin boundary.

**Example**:
```python
>>> twin_data = niti_twinning('type1')
        >>> b_parent = np.array([1., 0., 0.])  # Parent Burgers vector
        >>> 
        >>> disloc = get_twinning_dislocation(twin_data, b_parent)
# ...
```

---

## get_twinning_plane_points

```python
def get_twinning_plane_points(K1,Pointsout,horizontal,vertical):
```

Get atoms on twinning plane.
    Identifies atoms lying on twin boundary for twinning analysis.

**Example**:
```python
>>> # Type-I twin in cubic
        >>> L = cubic_lattice_vec(3.0)
        >>> K1 = np.array([1, 1, 1])  # (111) twin plane
        >>> twin_atoms = get_twinning_plane_points(L, K1, n1=10, n2=10, n3=10)
# ...
```

---

## get_twinningdata

```python
def get_twinningdata(orim,eus,Ldir_css,twin_systems,twt,phase, tension=True):
```

Extract complete twinning data from parent and twin lattices.
    Analyzes lattice relationship to determine all twinning elements.

**Example**:
```python
>>> L_parent = monoclinic_lattice_vec(2.889, 4.120, 4.622, 96.8)
        >>> # Twin by compound operation
        >>> R_twin = rotation_from_axis_angle([1,1,1], np.radians(70.5))
        >>> L_twin = R_twin.dot(L_parent)
# ...
```

---

## get_unique_families

```python
def get_unique_families(hkls):
```

Returns unique families of Miller indices based on permutation symmetry.
    Families are considered equivalent if they are permutations of each other
    (considering absolute values).

**Example**:
```python
>>> hkls = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1)]
        >>> families = get_unique_families(hkls)
        >>> print(families)
        {(1, 0, 0): 3, (1, 1, 0): 2}
# ...
```

---

## habitplane_equation_solution

```python
def habitplane_equation_solution(Uj,Ui,Qj,n,a,tol=1e-10):
```

Solve habit plane equation for phase transformation.
    Finds habit plane normal that satisfies invariant plane strain condition:
    F·n = λn + s  (where s is shear direction)

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Martensitic transformation
        >>> F = np.diag([1.05, 0.98, 0.97])  # Deformation gradient
# ...
```

---

## hkil2hkl

```python
def hkil2hkl(hkil):
```

Convert 4-index hexagonal plane (hkil) to 3-index (hkl).
    Transforms Miller-Bravais 4-index plane notation to standard 3-index
    notation for hexagonal crystal systems.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Convert (11-20) plane
        >>> hkil = np.array([1, 1, -2, 0])
# ...
```

---

## hkl2hkil

```python
def hkl2hkil(hkl):
```

Convert 3-index hexagonal plane (hkl) to 4-index (hkil).
    Transforms standard 3-index plane notation to Miller-Bravais 4-index
    notation for hexagonal crystal systems. The redundant index i = -(h+k).

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Convert (110) plane
        >>> hkl = np.array([1, 1, 0])
# ...
```

---

## kronecker

```python
def kronecker():
```

Compute Kronecker delta δ_{ij}.
    Returns 1 if indices are equal, 0 otherwise. Fundamental in
    tensor algebra and represents the identity tensor.

**Example**:
```python
>>> # Identity matrix via Kronecker delta
        >>> import numpy as np
        >>> I = np.array([[kronecker(i,j) for j in range(3)] for i in range(3)])
        >>> print(I)
# ...
```

---

## lattice_correspondence

```python
def lattice_correspondence(LatCorr,parent_symops,product_symops):
```

Calculate lattice correspondence matrix between two crystal structures.
    Finds transformation matrix C relating two lattices: L2 = C · L1.
    Optionally optimizes to minimize lattice mismatch.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Cubic to tetragonal
        >>> L_cubic = cubic_lattice_vec(3.0)
# ...
```

---

## lattice_vec

```python
def lattice_vec(lattice_param):
```

Calculate lattice vectors (a1, a2, a3) for various crystal systems.

**Example**:
```python
>>> # Cubic crystal (e.g., Silicon)
        >>> cubic_params = {'type': 'cubic', 'a': 5.43}
        >>> a1, a2, a3 = lattice_vec(cubic_params)
        >>> print("a1:", a1)
# ...
```

---

## miller2fractional

```python
def miller2fractional(uvw,frac=10,eps2=1e-2,decimals=5):
```

Reduce Miller indices to lowest integer form.
    Takes potentially fractional Miller indices and reduces them to the
    simplest integer representation by finding rational approximations

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Reduce fractional indices
        >>> uvw = np.array([0.5, 1.0, 1.5])
# ...
```

---

## mohr_circles

```python
def mohr_circles(tensor):
```

Calculate Mohr's circles from strain or stress tensor.
    Computes principal strains/stresses and parameters for three Mohr's
    circles. Used for visualizing 3D strain state in 2D.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Uniaxial strain state
        >>> strain = np.array([[0.1,  0.0, 0.0],
# ...
```

---

## monoclinic_lattice_vec

```python
def monoclinic_lattice_vec(a,b,c,beta):
```

Generate monoclinic lattice vectors (B19' martensite).
    Creates 3×3 lattice matrix for monoclinic crystal system.
    One unique angle (β) differs from 90°, typical of B19' martensite in NiTi.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # B19' martensite (NiTi)
        >>> LM = monoclinic_lattice_vec(2.889, 4.120, 4.622, 96.8)
# ...
```

---

## niti_twinning

```python
def niti_twinning(B2_symops,B2_recsymops,B19p_recsymops,B19p_symops,Uv,Parent_uvw2xyz,Parent_hkl2xyz,Product_uvw2xyz, Product_hkl2xyz,Parent_uvw_2_Product_uvw_all, Parent_hkl_2_Product_hkl_all, Parent_uvw_2_Product_uvw_all_norm,miller='greaterthanone',Qv=None):
```

Get NiTi twinning parameters for specific variant.
    Returns crystallographic twinning elements for NiTi B19' martensite.

**Example**:
```python
>>> # Type-I twin (most common in NiTi)
        >>> twin_data = niti_twinning('type1')
        >>> print(f"Twin plane K1: {twin_data['K1']}")
        >>> print(f"Twin direction η1: {twin_data['eta1']}")
# ...
```

---

## normArrayColumns

```python
def normArrayColumns(arr):
```

Normalize each column of matrix to unit length.
    Divides each column by its Euclidean norm, creating an orthonormal
    or semi-orthonormal matrix. Commonly used in crystallography for

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Normalize lattice vectors
        >>> lattice = np.array([[3, 0, 0],
# ...
```

---

## np_kronecker

```python
def np_kronecker():
```

Compute Kronecker delta (NumPy version).
    NumPy-compatible version of kronecker().

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Create identity matrix
        >>> n = 3
# ...
```

---

## np_permut_tensor3

```python
def np_permut_tensor3():
```

Compute Levi-Civita symbol (NumPy-compatible version).
    Same as permut_tensor3 but optimized for NumPy array indexing.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Build full permutation tensor
        >>> epsilon = np.zeros((3,3,3))
# ...
```

---

## permut_tensor3

```python
def permut_tensor3():
```

Compute 3D Levi-Civita permutation symbol ε_{ijk}.
    Returns the permutation tensor (Levi-Civita symbol) for three indices.
    Used in cross products, curl operations, and tensor calculations.

**Example**:
```python
>>> # Even permutations
        >>> print(permut_tensor3(0, 1, 2))  # +1 (xyz)
        >>> print(permut_tensor3(1, 2, 0))  # +1 (yzx)
        >>> print(permut_tensor3(2, 0, 1))  # +1 (zxy)
# ...
```

---

## perpendicular_vector

```python
def perpendicular_vector(v):
```

Find arbitrary unit vector perpendicular to input 3D vector.
    Computes a normalized vector perpendicular to the input using cross product
    with a judiciously chosen auxiliary vector. Used in crystallographic

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Find perpendicular to [1,2,3]
        >>> v = np.array([1, 2, 3])
# ...
```

---

## plane2string

```python
def plane2string(v, digits=2):
```

Format plane normal as string representation (h,k,l) with Miller notation.
    Converts plane normal vector to crystallographic notation using
    parentheses, following standard Miller index convention.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Format {111} plane
        >>> hkl = [1, 1, 1]
# ...
```

---

## plane_line_intersection

```python
def plane_line_intersection(n,V0,P0,P1):
```

---

## print_correspondence

```python
def print_correspondence(Mcorr,VecA,latticeA, latticeB,planes=False,returnB=False):
```

Print lattice correspondence matrix in readable format.
    Displays correspondence matrix with phase labels for documentation
    and reporting purposes.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> C = B19p_B2_lattice_correspondence()
        >>> print_correspondence(C, 'B19 martensite', 'B2 austenite')
# ...
```

---

## read_txt

```python
def read_txt(filename,delimiter='\t',skiprows=1):
```

Read numerical data from text file.
    Loads data written by write_txt() or similar format.

**Example**:
```python
>>> # Read previously saved data
        >>> data = read_txt('orientations.txt')
        >>> print(f"Loaded {len(data)} orientations")
        >>> print(f"Shape: {data.shape}")
# ...
```

---

## reciprocal_basis

```python
def reciprocal_basis(a1,a2,a3):
```

Calculate reciprocal lattice basis vectors from real space lattice vectors.
    The reciprocal lattice vectors satisfy: a_i · b_j = 2π δ_ij
    (Note: This implementation uses the crystallographic convention without 2π factor)

**Example**:
```python
>>> import numpy as np
        >>> # Define cubic lattice
        >>> a = 5.0  # Angstroms
        >>> a1 = np.array([a, 0, 0])
# ...
```

---

## select_atomic_plane

```python
def select_atomic_plane(LatticePoints,normal,eps=1e-1,shift=0.,eps2=None):
```

Select atoms lying on or near a crystallographic plane.
    Finds atoms within tolerance distance from specified plane.

**Example**:
```python
>>> L = cubic_lattice_vec(3.0)
        >>> atoms = generate_lattite_atom_positions(L, n1=3, n2=3, n3=3)
        >>> 
        >>> # Select atoms on (001) plane at z=6
# ...
```

---

## select_atomic_region

```python
def select_atomic_region(LatticePoints,normal,side='lower',eps=1e-1,shift=0.):
```

Select atoms within spherical region.
    Finds all atoms within specified radius of center point.

**Example**:
```python
>>> atoms = generate_lattite_atom_positions(L, n1=10, n2=10, n3=10)
        >>> 
        >>> # Select atoms near point [15, 15, 15]
        >>> center = np.array([15, 15, 15])
# ...
```

---

## select_crystal_planes

```python
def select_crystal_planes(NormalsOnCircle,ShearsOnCircle,maxhkl):
```

Generate list of low-index crystallographic planes.
    Creates planes (hkl) with indices up to max_index for analysis.

**Example**:
```python
>>> L = cubic_lattice_vec(3.0)
        >>> planes = select_crystal_planes(L, max_index=2)
        >>> 
        >>> print(f"Generated {len(planes)} planes")
# ...
```

---

## select_plane

```python
def select_plane(LatticeVectors,normal,eps=1e-1,shift=0.,Q=np.eye(3)):
```

select_plane - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## strains_along_13mohrcirle

```python
def strains_along_13mohrcirle(Strain,VV,normdiri,phi_around_V2,Parent_xyz2hkl):
```

Calculate strains along maximum Mohr's circle.
    Computes normal and shear strain components for points on the
    largest Mohr's circle (ε₁-ε₃ circle).

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> strain = np.diag([0.1, 0.0, -0.05])
        >>> mohr_result = mohr_circles(strain)
# ...
```

---

## tetragonal_lattice_vec

```python
def tetragonal_lattice_vec(a,b,c):
```

Generate tetragonal lattice vectors.
    Creates 3×3 lattice matrix for tetragonal crystal system.
    Two equal parameters (a=b) and one unique (c), all angles 90°.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Tetragonal system with a=b=3.0, c=4.0
        >>> L_tet = tetragonal_lattice_vec(3.0, 3.0, 4.0)
# ...
```

---

## twin_equation_solution

```python
def twin_equation_solution(Uj,Ui,L_A,Lr_A,L_M,Lr_M, R_AM, Ci_d,Ci_p,tol=1e-10,miller='greaterthanone',printlambda=False, Qj=None,Qi=None):
```

Solve complete twin equation system.
    Finds K1, K2, η1, η2 twin elements from lattice parameters.

**Example**:
```python
>>> L_parent = cubic_lattice_vec(3.0)
        >>> # Create twinned lattice
        >>> theta = np.radians(70.5)  # Twin rotation angle
        >>> R_twin = rotation_from_axis_angle([1,1,1], theta)
# ...
```

---

## twin_equation_solution_ini

```python
def twin_equation_solution_ini(Uj,Ui,Parent_uvw2xyz,Parent_hkl2xyz,Product_uvw2xyz,Product_hkl2xyz, Parent_uvw_2_Product_uvw_rot, Parent_uvw_2_Product_uvw,Parent_hkl_2_Product_hkl,tol=1e-10,miller='greaterthanone',printlambda=False, Qj=None,Qi=None):
```

Initialize twin equation solver with default parameters.
    Sets up initial guess and parameters for twin equation solution.

**Example**:
```python
>>> params = twin_equation_solution_ini()
        >>> # Modify parameters as needed
        >>> params['max_iter'] = 1000
        >>> params['tolerance'] = 1e-8
# ...
```

---

## twinnedhabitplane

```python
def twinnedhabitplane(Ui,Uj,Qij,a1,n1,hbplanes=[],addondata={},method='bhata'):
```

Calculate twinned habit plane from lattice parameters.
    Determines twin plane (K1) and twin direction from parent and
    twin lattice geometries.

**Example**:
```python
>>> # Type-I twin in cubic
        >>> L_parent = cubic_lattice_vec(3.0)
        >>> # Twin related by mirror symmetry
        >>> C_twin = np.diag([1, 1, -1])  # Mirror in (001)
# ...
```

---

## uvtw2uvw

```python
def uvtw2uvw(uvtw):
```

Convert 4-index hexagonal direction [uvtw] to 3-index [uvw].
    Transforms Miller-Bravais 4-index notation to standard 3-index notation
    for hexagonal crystal systems. The redundant index t is removed.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Convert [11-20] direction
        >>> uvtw = np.array([1, 1, -2, 0])
# ...
```

---

## uvw2uvtw

```python
def uvw2uvtw(uvw):
```

Convert 3-index hexagonal direction [uvw] to 4-index [uvtw].
    Transforms standard 3-index notation to Miller-Bravais 4-index notation
    for hexagonal crystal systems. The redundant index t = -(u+v) is inserted.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Convert [110] direction
        >>> uvw = np.array([1, 1, 0])
# ...
```

---

## vec2string

```python
def vec2string(v, digits=2):
```

Format vector as string representation [x,y,z] with specified precision.
    Converts numpy array or list to readable string format suitable for
    display, logging, or file output in crystallographic applications.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Format with default 2 decimal places
        >>> v = [1.23456, 0.56789, 0.89123]
# ...
```

---

## vector2miller_ini

```python
def vector2miller_ini(v, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3):
```

Convert Cartesian vector to Miller indices (initial guess).
    Provides starting point for vector→Miller conversion.

**Example**:
```python
>>> L = cubic_lattice_vec(3.0)
        >>> v_cart = np.array([3, 3, 0])
        >>> uvw_guess = vector2miller_ini(v_cart, L)
        >>> print(uvw_guess)  # Initial guess for [110]
```

---

## vectors2miller

```python
def vectors2miller(V, MIN=True, Tol=1e-9,tol=1e5,text=False):
```

Convert multiple vectors to Miller indices.

---

## write_lattice_correspondence

```python
def write_lattice_correspondence(ax,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,Product_lattice,Parent_lattice,FontSize=10,Fontweight="bold"):
```

Write lattice correspondence analysis to file.
    Saves correspondence matrix and lattice mismatch information.

**Example**:
```python
>>> L_B2 = cubic_lattice_vec(3.015)
        >>> L_B19p = monoclinic_lattice_vec(2.889, 4.120, 4.622, 96.8)
        >>> C = lattice_correspondence(L_B19p, L_B2)
        >>> 
# ...
```

---

## write_mohr_planes

```python
def write_mohr_planes(ax,Upperhalftext,Lowerhalftext,colors,markersize=8,markeredgewidth=2):
```

Write Mohr circle analysis results to file.
    Saves principal strains, circles, and plane-specific strains.

**Example**:
```python
>>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> planes = [[1,0,0], [1,1,0], [1,1,1]]
        >>> L = cubic_lattice_vec(3.0)
# ...
```

---

## write_txt

```python
def write_txt(filename,Header,DATA):
```

Write numerical data to text file.
    Saves array data in formatted text file with optional header.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Save orientation data
        >>> euler_angles = np.random.rand(100, 3) * [360, 180, 360]
# ...
```

---

## xyz2fractional

```python
def xyz2fractional(Txyz2uvw,V,frac=10,eps2=1e-2,decimals=5):
```

Convert Cartesian coordinates to fractional Miller indices with reduction.
    Transforms Cartesian vector to Miller indices using transformation matrix,
    then reduces to lowest integer form. Handles fractional indices by finding

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Cubic lattice with a=3.0
        >>> a = 3.0
# ...
```

---

## xyz2fractional02

```python
def xyz2fractional02(Txyz2uvw,V):
```

Simple Cartesian to fractional coordinate transformation without reduction.
    Performs basic coordinate transformation using matrix multiplication,
    then normalizes result. Does not reduce to Miller indices or apply GCD.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Transform coordinates in cubic system
        >>> a = 3.0
# ...
```

---

## zero_normal_strains

```python
def zero_normal_strains(Strain, mcircles,VV,normdiri,phi_around_normdiri,Parent_xyz2hkl):
```

Find planes with zero normal strain in given strain state.
    Calculates planes where normal strain εₙ = nᵀ·ε·n = 0.
    These planes experience only shear.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Pure shear state
        >>> strain = np.array([[0.1,  0.05, 0.0],
# ...
```

---


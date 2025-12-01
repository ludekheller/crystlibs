# CRYSTLIB - DOCUMENTATION SUMMARY
## Organized Function Reference

**Module**: `crystlib.py` (Extended Version)  
**Total Functions**: 73  
**Version**: Extended 2024  

---

## 📊 FUNCTION OVERVIEW

### By Alphabetical Order

| Function | Purpose |
|----------|---------|
| `B19p_B2_lattice_correspondence` | Generate B19'→B2 lattice correspondence matrix for NiTi.... |
| `B19p_B2_lattice_correspondence_ini` | Initialize B19'→B2 correspondence for specific transformation variant.... |
| `Rp_B2_lattice_correspondence` | Generate R-phase→B2 lattice correspondence for NiTi.... |
| `an_between_vecs` | Calculate angle between two vectors.... |
| `array2tuple` | Convert a numpy array to a tuple with rounded elements.... |
| `cubic2tetragonal_lattice_correspondence` | Generate cubic→tetragonal lattice correspondence.... |
| `cubic_lattice_vec` | Generate cubic lattice vectors.... |
| `def_gradient` | Calculate total deformation gradient with rotation.... |
| `def_gradient_ini` | Initialize deformation gradient with Euler angle rotation.... |
| `def_gradient_ini2` | Initialize deformation gradient with axis-angle rotation.... |
| `def_gradient_stressfree` | Calculate stress-free transformation strain (deformation gradient).... |
| `def_gradient_stressfree_ini` | Initialize stress-free deformation gradient from lattice parameters.... |
| `dir2string` | Format direction as string representation [u,v,w] with Miller notation.... |
| `equivalent_elements` | Generate all symmetry-equivalent crystallographic elements.... |
| `find_gcd` | Find greatest common divisor using recursive Euclidean algorithm.... |
| `flipvector` | Flip vector to ensure positive first non-zero component.... |
| `flipvector2negative` | Flip vector to ensure negative first non-zero component.... |
| `gen_twinned_lattice_points` | Generate lattice points for parent and twinned regions.... |
| `genallHexSys` | Generate all hexagonal slip systems including <c+a> pyramidal.... |
| `generate_hkls` | Generate unique Miller indices (hkl) considering crystal symmetry operations.... |
| `generate_hkls01` | Generate unique Miller indices (hkl) considering crystal symmetry operations.... |
| `generate_hkls02` | Generate unique Miller indices (hkl) with metric tensor transformation.... |
| `generate_lattice_faces` | Generate face polygons for unit cell visualization.... |
| `generate_lattice_points` | Generate lattice points within specified unit cell range.... |
| `generate_lattice_vectors` | Generate lattice vectors for visualization.... |
| `generate_lattite_atom_positions` | Generate atomic positions in lattice with basis.... |
| `generate_plane_vertices` | Generate vertices for plotting crystallographic plane.... |
| `generate_product_lattice_faces` | Generate faces for two lattices.... |
| `generate_product_lattice_points` | Generate lattice points for two overlapping lattices.... |
| `gensystemsHex` | Generate hexagonal slip systems in real space.... |
| `gensystemsHexIni` | Generate initial hexagonal slip system templates.... |
| `get_interface2d` | Extract interface atoms from two lattices.... |
| `get_twinning_dislocation` | Calculate twinning dislocation Burgers vector.... |
| `get_twinning_plane_points` | Get atoms on twinning plane.... |
| `get_twinningdata` | Extract complete twinning data from parent and twin lattices.... |
| `get_unique_families` | Returns unique families of Miller indices based on permutation symmetry.... |
| `habitplane_equation_solution` | Solve habit plane equation for phase transformation.... |
| `hkil2hkl` | Convert 4-index hexagonal plane (hkil) to 3-index (hkl).... |
| `hkl2hkil` | Convert 3-index hexagonal plane (hkl) to 4-index (hkil).... |
| `kronecker` | Compute Kronecker delta δ_{ij}.... |
| `lattice_correspondence` | Calculate lattice correspondence matrix between two crystal structures.... |
| `lattice_vec` | Calculate lattice vectors (a1, a2, a3) for various crystal systems.... |
| `miller2fractional` | Reduce Miller indices to lowest integer form.... |
| `monoclinic_lattice_vec` | Generate monoclinic lattice vectors (B19' martensite).... |
| `niti_twinning` | Get NiTi twinning parameters for specific variant.... |
| `normArrayColumns` | Normalize each column of matrix to unit length.... |
| `np_kronecker` | Compute Kronecker delta (NumPy version).... |
| `np_permut_tensor3` | Compute Levi-Civita symbol (NumPy-compatible version).... |
| `permut_tensor3` | Compute 3D Levi-Civita permutation symbol ε_{ijk}.... |
| `perpendicular_vector` | Find arbitrary unit vector perpendicular to input 3D vector.... |
| `plane2string` | Format plane normal as string representation (h,k,l) with Miller notation.... |
| `plane_line_intersection` | Calculate intersection of plane and line in 3D.... |
| `print_correspondence` | Print lattice correspondence matrix in readable format.... |
| `read_txt` | Read numerical data from text file.... |
| `reciprocal_basis` | Calculate reciprocal lattice basis vectors from real space lattice vectors.... |
| `select_atomic_plane` | Select atoms lying on or near a crystallographic plane.... |
| `select_atomic_region` | Select atoms within spherical region.... |
| `select_crystal_planes` | Generate list of low-index crystallographic planes.... |
| `select_plane` | Select points on crystallographic plane (hkl).... |
| `symmetry_elements` | Generate symmetry operation matrices for crystal system.... |
| `tetragonal_lattice_vec` | Generate tetragonal lattice vectors.... |
| `twin_equation_solution` | Solve complete twin equation system.... |
| `twin_equation_solution_ini` | Initialize twin equation solver with default parameters.... |
| `twinnedhabitplane` | Calculate twinned habit plane from lattice parameters.... |
| `uvtw2uvw` | Convert 4-index hexagonal direction [uvtw] to 3-index [uvw].... |
| `uvw2uvtw` | Convert 3-index hexagonal direction [uvw] to 4-index [uvtw].... |
| `vec2string` | Format vector as string representation [x,y,z] with specified precision.... |
| `vector2miller_ini` | Convert Cartesian vector to Miller indices (initial guess).... |
| `vectors2miller` | Convert Cartesian vector to Miller indices with optimization.... |
| `write_lattice_correspondence` | Write lattice correspondence analysis to file.... |
| `write_txt` | Write numerical data to text file.... |
| `xyz2fractional` | Convert Cartesian coordinates to fractional Miller indices with reduction.... |
| `xyz2fractional02` | Simple Cartesian to fractional coordinate transformation without reduction.... |

---

## 📖 DETAILED FUNCTION DESCRIPTIONS

### `B19p_B2_lattice_correspondence`

**Signature**: `def B19p_B2_lattice_correspondence()`

Generate B19'→B2 lattice correspondence matrix for NiTi.
    Returns the transformation matrix relating B19' monoclinic martensite
    to B2 cubic austenite lattice vectors. Fundamental for NiTi shape
    memory alloy analysis.

---

### `B19p_B2_lattice_correspondence_ini`

**Signature**: `def B19p_B2_lattice_correspondence_ini(variant=1)`

Initialize B19'→B2 correspondence for specific transformation variant.
    Returns correspondence matrix for one of the crystallographic variants
    of the B19' martensite transformation in NiTi.

---

### `Rp_B2_lattice_correspondence`

**Signature**: `def Rp_B2_lattice_correspondence()`

Generate R-phase→B2 lattice correspondence for NiTi.
    Returns correspondence matrix for R-phase (rhombohedral) to B2
    transformation in NiTi alloys. R-phase is intermediate phase.

---

### `an_between_vecs`

**Signature**: `def an_between_vecs(v1, v2, degrees=True)`

Calculate angle between two vectors.
    Computes angle using dot product formula.

---

### `array2tuple`

**Signature**: `def array2tuple(arr, decimals=2)`

Convert a numpy array to a tuple with rounded elements.

---

### `cubic2tetragonal_lattice_correspondence`

**Signature**: `def cubic2tetragonal_lattice_correspondence()`

Generate cubic→tetragonal lattice correspondence.
    Returns correspondence matrix for cubic to tetragonal phase
    transformation (e.g., FCC→FCT).

---

### `cubic_lattice_vec`

**Signature**: `def cubic_lattice_vec(a)`

Generate cubic lattice vectors.
    Creates 3×3 lattice matrix for cubic crystal system with parameter a.
    All angles are 90° and all lengths are equal (a=b=c).

---

### `def_gradient`

**Signature**: `def def_gradient(L1, L2, R)`

Calculate total deformation gradient with rotation.
    Computes F = R · F₀ where R is rigid rotation and F₀ is lattice strain.

---

### `def_gradient_ini`

**Signature**: `def def_gradient_ini(lattice_params_1, lattice_params_2, euler_angles)`

Initialize deformation gradient with Euler angle rotation.
    Calculates F with rotation specified by Euler angles.

---

### `def_gradient_ini2`

**Signature**: `def def_gradient_ini2(L1, L2, axis, angle)`

Initialize deformation gradient with axis-angle rotation.
    Alternative initialization using axis-angle for rotation.

---

### `def_gradient_stressfree`

**Signature**: `def def_gradient_stressfree(L1, L2)`

Calculate stress-free transformation strain (deformation gradient).
    Computes F₀ = L2 · L1⁻¹ representing lattice transformation
    without applied stress.

---

### `def_gradient_stressfree_ini`

**Signature**: `def def_gradient_stressfree_ini(lattice_params_1, lattice_params_2)`

Initialize stress-free deformation gradient from lattice parameters.
    Convenience function to calculate F₀ from lattice parameters directly.

---

### `dir2string`

**Signature**: `def dir2string(v, digits=2)`

Format direction as string representation [u,v,w] with Miller notation.
    Converts direction vector to crystallographic notation using
    square brackets, following standard Miller index convention.

---

### `equivalent_elements`

**Signature**: `def equivalent_elements(element, symops)`

Generate all symmetry-equivalent crystallographic elements.
    Applies all symmetry operations to a direction or plane normal to find
    all equivalent elements. Used in texture analysis and pole figures.

---

### `find_gcd`

**Signature**: `def find_gcd(x, y)`

Find greatest common divisor using recursive Euclidean algorithm.
    Fundamental mathematical operation used throughout the module for
    Miller indices reduction and fractional coordinate normalization.
    Implements the classical Euclidean algorithm via recursion.

---

### `flipvector`

**Signature**: `def flipvector(v)`

Flip vector to ensure positive first non-zero component.
    Standardizes vector direction for comparison.

---

### `flipvector2negative`

**Signature**: `def flipvector2negative(v)`

Flip vector to ensure negative first non-zero component.
    Opposite of flipvector() - ensures negative first component.

---

### `gen_twinned_lattice_points`

**Signature**: `def gen_twinned_lattice_points(L_parent, twin_plane_normal, n1=5, n2=5, n3=5)`

Generate lattice points for parent and twinned regions.
    Creates atomic positions showing twin boundary.

---

### `genallHexSys`

**Signature**: `def genallHexSys(L)`

Generate all hexagonal slip systems including <c+a> pyramidal.
    Comprehensive slip system generation for HCP crystals including
    basal, prismatic, and pyramidal <c+a> systems. Essential for
    complete plasticity modeling of hexagonal materials.

---

### `generate_hkls`

**Signature**: `def generate_hkls(hklmax, syms, hkls=[])`

Generate unique Miller indices (hkl) considering crystal symmetry operations.
    This version rounds values to avoid floating point comparison issues.

---

### `generate_hkls01`

**Signature**: `def generate_hkls01(hklmax, syms, hkls=[])`

Generate unique Miller indices (hkl) considering crystal symmetry operations.
    Alternative version without rounding.

---

### `generate_hkls02`

**Signature**: `def generate_hkls02(hklmax, syms, G, hkls=[])`

Generate unique Miller indices (hkl) with metric tensor transformation.
    This version applies symmetry operations in the correct metric space.

---

### `generate_lattice_faces`

**Signature**: `def generate_lattice_faces(L)`

Generate face polygons for unit cell visualization.
    Creates vertex lists for 6 faces of unit cell for 3D rendering.

---

### `generate_lattice_points`

**Signature**: `def generate_lattice_points(L, n1=1, n2=1, n3=1)`

Generate lattice points within specified unit cell range.
    Creates array of lattice points for visualization and analysis.
    Generates n1×n2×n3 unit cells.

---

### `generate_lattice_vectors`

**Signature**: `def generate_lattice_vectors(L, n1=1, n2=1, n3=1, include_origin=True)`

Generate lattice vectors for visualization.
    Creates list of vectors from origin to lattice points for
    arrow-based visualization.

---

### `generate_lattite_atom_positions`

**Signature**: `def generate_lattite_atom_positions(L, basis=None, n1=1, n2=1, n3=1)`

Generate atomic positions in lattice with basis.
    Creates Cartesian coordinates of atoms including basis atoms
    within each unit cell.

---

### `generate_plane_vertices`

**Signature**: `def generate_plane_vertices(h, k, l, L, scale=1.0)`

Generate vertices for plotting crystallographic plane.
    Creates polygon vertices to visualize plane (hkl) in 3D.

---

### `generate_product_lattice_faces`

**Signature**: `def generate_product_lattice_faces(L1, L2)`

Generate faces for two lattices.
    Creates face polygons for both unit cells.

---

### `generate_product_lattice_points`

**Signature**: `def generate_product_lattice_points(L1, L2, n1=1, n2=1, n3=1)`

Generate lattice points for two overlapping lattices.
    Creates points for both lattices to visualize phase coexistence
    or transformation.

---

### `gensystemsHex`

**Signature**: `def gensystemsHex(L)`

Generate hexagonal slip systems in real space.
    Transforms slip system templates to actual Cartesian coordinates
    using the hexagonal lattice matrix. Produces slip directions and
    plane normals for deformation analysis.

---

### `gensystemsHexIni`

**Signature**: `def gensystemsHexIni()`

Generate initial hexagonal slip system templates.
    Creates basal, prismatic, and pyramidal slip system templates for
    hexagonal close-packed (HCP) crystal structures. Returns normalized
    direction and plane normal vectors.

---

### `get_interface2d`

**Signature**: `def get_interface2d(atoms1, atoms2, plane_normal, plane_point, tolerance=0.1)`

Extract interface atoms from two lattices.
    Finds atoms from both structures near interface plane.
    Used in phase boundary visualization.

---

### `get_twinning_dislocation`

**Signature**: `def get_twinning_dislocation(twin_data, b_parent)`

Calculate twinning dislocation Burgers vector.
    Determines dislocation content at twin boundary.

---

### `get_twinning_plane_points`

**Signature**: `def get_twinning_plane_points(L, twin_plane_normal, n1=5, n2=5, n3=5, tolerance=0.1)`

Get atoms on twinning plane.
    Identifies atoms lying on twin boundary for twinning analysis.

---

### `get_twinningdata`

**Signature**: `def get_twinningdata(L_parent, L_twin)`

Extract complete twinning data from parent and twin lattices.
    Analyzes lattice relationship to determine all twinning elements.

---

### `get_unique_families`

**Signature**: `def get_unique_families(hkls)`

Returns unique families of Miller indices based on permutation symmetry.
    Families are considered equivalent if they are permutations of each other
    (considering absolute values).

---

### `habitplane_equation_solution`

**Signature**: `def habitplane_equation_solution(F, s)`

Solve habit plane equation for phase transformation.
    Finds habit plane normal that satisfies invariant plane strain condition:
    F·n = λn + s  (where s is shear direction)

---

### `hkil2hkl`

**Signature**: `def hkil2hkl(hkil)`

Convert 4-index hexagonal plane (hkil) to 3-index (hkl).
    Transforms Miller-Bravais 4-index plane notation to standard 3-index
    notation for hexagonal crystal systems.

---

### `hkl2hkil`

**Signature**: `def hkl2hkil(hkl)`

Convert 3-index hexagonal plane (hkl) to 4-index (hkil).
    Transforms standard 3-index plane notation to Miller-Bravais 4-index
    notation for hexagonal crystal systems. The redundant index i = -(h+k).

---

### `kronecker`

**Signature**: `def kronecker(i, j)`

Compute Kronecker delta δ_{ij}.
    Returns 1 if indices are equal, 0 otherwise. Fundamental in
    tensor algebra and represents the identity tensor.

---

### `lattice_correspondence`

**Signature**: `def lattice_correspondence(L1, L2, optimize=False)`

Calculate lattice correspondence matrix between two crystal structures.
    Finds transformation matrix C relating two lattices: L2 = C · L1.
    Optionally optimizes to minimize lattice mismatch.

---

### `lattice_vec`

**Signature**: `def lattice_vec(lattice_param)`

Calculate lattice vectors (a1, a2, a3) for various crystal systems.

---

### `miller2fractional`

**Signature**: `def miller2fractional(uvw, frac=10, eps2=1e-2, decimals=5)`

Reduce Miller indices to lowest integer form.
    Takes potentially fractional Miller indices and reduces them to the
    simplest integer representation by finding rational approximations
    and applying GCD reduction.

---

### `monoclinic_lattice_vec`

**Signature**: `def monoclinic_lattice_vec(a, b, c, beta)`

Generate monoclinic lattice vectors (B19' martensite).
    Creates 3×3 lattice matrix for monoclinic crystal system.
    One unique angle (β) differs from 90°, typical of B19' martensite in NiTi.

---

### `niti_twinning`

**Signature**: `def niti_twinning(variant='type1')`

Get NiTi twinning parameters for specific variant.
    Returns crystallographic twinning elements for NiTi B19' martensite.

---

### `normArrayColumns`

**Signature**: `def normArrayColumns(arr)`

Normalize each column of matrix to unit length.
    Divides each column by its Euclidean norm, creating an orthonormal
    or semi-orthonormal matrix. Commonly used in crystallography for
    basis vector normalization.

---

### `np_kronecker`

**Signature**: `def np_kronecker(i, j)`

Compute Kronecker delta (NumPy version).
    NumPy-compatible version of kronecker().

---

### `np_permut_tensor3`

**Signature**: `def np_permut_tensor3(i, j, k)`

Compute Levi-Civita symbol (NumPy-compatible version).
    Same as permut_tensor3 but optimized for NumPy array indexing.

---

### `permut_tensor3`

**Signature**: `def permut_tensor3(i, j, k)`

Compute 3D Levi-Civita permutation symbol ε_{ijk}.
    Returns the permutation tensor (Levi-Civita symbol) for three indices.
    Used in cross products, curl operations, and tensor calculations.

---

### `perpendicular_vector`

**Signature**: `def perpendicular_vector(v)`

Find arbitrary unit vector perpendicular to input 3D vector.
    Computes a normalized vector perpendicular to the input using cross product
    with a judiciously chosen auxiliary vector. Used in crystallographic
    calculations requiring orthogonal basis construction.

---

### `plane2string`

**Signature**: `def plane2string(v, digits=2)`

Format plane normal as string representation (h,k,l) with Miller notation.
    Converts plane normal vector to crystallographic notation using
    parentheses, following standard Miller index convention.

---

### `plane_line_intersection`

**Signature**: `def plane_line_intersection(plane_point, plane_normal, line_point, line_direction)`

Calculate intersection of plane and line in 3D.
    Finds point where line intersects plane, if it exists.

---

### `print_correspondence`

**Signature**: `def print_correspondence(C, phase1='Phase1', phase2='Phase2')`

Print lattice correspondence matrix in readable format.
    Displays correspondence matrix with phase labels for documentation
    and reporting purposes.

---

### `read_txt`

**Signature**: `def read_txt(filename)`

Read numerical data from text file.
    Loads data written by write_txt() or similar format.

---

### `reciprocal_basis`

**Signature**: `def reciprocal_basis(a1, a2, a3)`

Calculate reciprocal lattice basis vectors from real space lattice vectors.
    The reciprocal lattice vectors satisfy: a_i · b_j = 2π δ_ij
    (Note: This implementation uses the crystallographic convention without 2π factor)

---

### `select_atomic_plane`

**Signature**: `def select_atomic_plane(atoms, plane_normal, plane_point, tolerance=0.1)`

Select atoms lying on or near a crystallographic plane.
    Finds atoms within tolerance distance from specified plane.

---

### `select_atomic_region`

**Signature**: `def select_atomic_region(atoms, center, radius)`

Select atoms within spherical region.
    Finds all atoms within specified radius of center point.

---

### `select_crystal_planes`

**Signature**: `def select_crystal_planes(lattice, max_index=3)`

Generate list of low-index crystallographic planes.
    Creates planes (hkl) with indices up to max_index for analysis.

---

### `select_plane`

**Signature**: `def select_plane(points, h, k, l, L, tolerance=0.1)`

Select points on crystallographic plane (hkl).
    Finds points lying on specified Miller plane.

---

### `symmetry_elements`

**Signature**: `def symmetry_elements(crystal_system)`

Generate symmetry operation matrices for crystal system.
    Returns list of 3×3 rotation matrices representing all symmetry
    operations for specified crystal system.

---

### `tetragonal_lattice_vec`

**Signature**: `def tetragonal_lattice_vec(a, b, c)`

Generate tetragonal lattice vectors.
    Creates 3×3 lattice matrix for tetragonal crystal system.
    Two equal parameters (a=b) and one unique (c), all angles 90°.

---

### `twin_equation_solution`

**Signature**: `def twin_equation_solution(L_parent, L_twin, params=None)`

Solve complete twin equation system.
    Finds K1, K2, η1, η2 twin elements from lattice parameters.

---

### `twin_equation_solution_ini`

**Signature**: `def twin_equation_solution_ini()`

Initialize twin equation solver with default parameters.
    Sets up initial guess and parameters for twin equation solution.

---

### `twinnedhabitplane`

**Signature**: `def twinnedhabitplane(L_parent, L_twin, correspondence_matrix)`

Calculate twinned habit plane from lattice parameters.
    Determines twin plane (K1) and twin direction from parent and
    twin lattice geometries.

---

### `uvtw2uvw`

**Signature**: `def uvtw2uvw(uvtw)`

Convert 4-index hexagonal direction [uvtw] to 3-index [uvw].
    Transforms Miller-Bravais 4-index notation to standard 3-index notation
    for hexagonal crystal systems. The redundant index t is removed.

---

### `uvw2uvtw`

**Signature**: `def uvw2uvtw(uvw)`

Convert 3-index hexagonal direction [uvw] to 4-index [uvtw].
    Transforms standard 3-index notation to Miller-Bravais 4-index notation
    for hexagonal crystal systems. The redundant index t = -(u+v) is inserted.

---

### `vec2string`

**Signature**: `def vec2string(v, digits=2)`

Format vector as string representation [x,y,z] with specified precision.
    Converts numpy array or list to readable string format suitable for
    display, logging, or file output in crystallographic applications.

---

### `vector2miller_ini`

**Signature**: `def vector2miller_ini(v, L)`

Convert Cartesian vector to Miller indices (initial guess).
    Provides starting point for vector→Miller conversion.

---

### `vectors2miller`

**Signature**: `def vectors2miller(v, L, max_index=5)`

Convert Cartesian vector to Miller indices with optimization.
    Finds best integer Miller indices representing given direction.

---

### `write_lattice_correspondence`

**Signature**: `def write_lattice_correspondence(filename, C, L1, L2, phase1='Phase1', phase2='Phase2')`

Write lattice correspondence analysis to file.
    Saves correspondence matrix and lattice mismatch information.

---

### `write_txt`

**Signature**: `def write_txt(filename, data, header='')`

Write numerical data to text file.
    Saves array data in formatted text file with optional header.

---

### `xyz2fractional`

**Signature**: `def xyz2fractional(Txyz2uvw, V, frac=10, eps2=1e-2, decimals=5)`

Convert Cartesian coordinates to fractional Miller indices with reduction.
    Transforms Cartesian vector to Miller indices using transformation matrix,
    then reduces to lowest integer form. Handles fractional indices by finding
    closest rational approximation with specified maximum denominator.

---

### `xyz2fractional02`

**Signature**: `def xyz2fractional02(Txyz2uvw, V)`

Simple Cartesian to fractional coordinate transformation without reduction.
    Performs basic coordinate transformation using matrix multiplication,
    then normalizes result. Does not reduce to Miller indices or apply GCD.

---


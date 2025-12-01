# CRYSTLIB - QUICK REFERENCE (Comprehensive)
## Fast Function Lookup

**Module**: `crystlib.py` (Extended Version)  
**Total Functions**: 73  

---

## 🔍 QUICK LOOKUP TABLE

| # | Function | Parameters | Returns |
|---|----------|------------|---------|
| 1 | `B19p_B2_lattice_correspondence` | `` | numpy.ndarray (3×3): Correspondence matr |
| 2 | `B19p_B2_lattice_correspondence_ini` | `variant=1` | numpy.ndarray (3×3): Correspondence matr |
| 3 | `Rp_B2_lattice_correspondence` | `` | numpy.ndarray (3×3): R-phase → B2 corres |
| 4 | `an_between_vecs` | `v1, v2, degrees=True` | float: Angle between vectors |
| 5 | `array2tuple` | `arr, decimals=2` | tuple - Rounded values as a tuple |
| 6 | `cubic2tetragonal_lattice_correspondence` | `` | numpy.ndarray (3×3): Correspondence matr |
| 7 | `cubic_lattice_vec` | `a` | numpy.ndarray (3×3): Lattice matrix [a1| |
| 8 | `def_gradient` | `L1, L2, R` | numpy.ndarray (3×3): Total deformation g |
| 9 | `def_gradient_ini` | `lattice_params_1, lattice_para...` | numpy.ndarray (3×3): Deformation gradien |
| 10 | `def_gradient_ini2` | `L1, L2, axis, angle` | numpy.ndarray (3×3): Deformation gradien |
| 11 | `def_gradient_stressfree` | `L1, L2` | numpy.ndarray (3×3): Stress-free deforma |
| 12 | `def_gradient_stressfree_ini` | `lattice_params_1, lattice_para...` | numpy.ndarray (3×3): Deformation gradien |
| 13 | `dir2string` | `v, digits=2` | str: Formatted string '[u.uu,v.vv,w.ww]' |
| 14 | `equivalent_elements` | `element, symops` | list: List of numpy.ndarray [3] - all un |
| 15 | `find_gcd` | `x, y` | int or float: Greatest common divisor of |
| 16 | `flipvector` | `v` | numpy.ndarray [3]: Flipped vector (if ne |
| 17 | `flipvector2negative` | `v` | numpy.ndarray [3]: Flipped vector (if ne |
| 18 | `gen_twinned_lattice_points` | `L_parent, twin_plane_normal, n...` | tuple: (parent_atoms, twin_atoms) - atom |
| 19 | `genallHexSys` | `L` | tuple: (all_directions, all_normals) |
| 20 | `generate_hkls` | `hklmax, syms, hkls=[]` | hkls: list of tuples - Unique Miller ind |
| 21 | `generate_hkls01` | `hklmax, syms, hkls=[]` | hkls: list of tuples - Unique Miller ind |
| 22 | `generate_hkls02` | `hklmax, syms, G, hkls=[]` | hkls: list of tuples - Unique Miller ind |
| 23 | `generate_lattice_faces` | `L` | list: List of 6 face vertex arrays (each |
| 24 | `generate_lattice_points` | `L, n1=1, n2=1, n3=1` | numpy.ndarray (N×3): Array of lattice po |
| 25 | `generate_lattice_vectors` | `L, n1=1, n2=1, n3=1, include_o...` | list: List of numpy arrays - lattice vec |
| 26 | `generate_lattite_atom_positions` | `L, basis=None, n1=1, n2=1, n3=...` | numpy.ndarray (N×3): Atomic positions in |
| 27 | `generate_plane_vertices` | `h, k, l, L, scale=1.0` | numpy.ndarray (N×3): Vertices of plane p |
| 28 | `generate_product_lattice_faces` | `L1, L2` | tuple: (faces1, faces2) - lists of face  |
| 29 | `generate_product_lattice_points` | `L1, L2, n1=1, n2=1, n3=1` | tuple: (points1, points2) - arrays of la |
| 30 | `gensystemsHex` | `L` | tuple: (directions, normals) in Cartesia |
| 31 | `gensystemsHexIni` | `` | tuple: (directions, normals) where each  |
| 32 | `get_interface2d` | `atoms1, atoms2, plane_normal, ...` | tuple: (indices1, indices2) - atom indic |
| 33 | `get_twinning_dislocation` | `twin_data, b_parent` | dict: { |
| 34 | `get_twinning_plane_points` | `L, twin_plane_normal, n1=5, n2...` | numpy.ndarray (N×3): Atoms on twin plane |
| 35 | `get_twinningdata` | `L_parent, L_twin` | dict: Complete twin characterization inc |
| 36 | `get_unique_families` | `hkls` | dict - Dictionary mapping representative |
| 37 | `habitplane_equation_solution` | `F, s` | dict: { |
| 38 | `hkil2hkl` | `hkil` | numpy.ndarray [3]: 3-index plane (h, k,  |
| 39 | `hkl2hkil` | `hkl` | numpy.ndarray [4]: 4-index plane (h, k,  |
| 40 | `kronecker` | `i, j` | int: 1 if i==j, 0 otherwise |
| 41 | `lattice_correspondence` | `L1, L2, optimize=False` | numpy.ndarray (3×3): Correspondence matr |
| 42 | `lattice_vec` | `lattice_param` | a1: numpy array (3,) - First lattice vec |
| 43 | `miller2fractional` | `uvw, frac=10, eps2=1e-2, decim...` | numpy.ndarray [3]: Reduced integer Mille |
| 44 | `monoclinic_lattice_vec` | `a, b, c, beta` | numpy.ndarray (3×3): Monoclinic lattice  |
| 45 | `niti_twinning` | `variant='type1'` | dict: { |
| 46 | `normArrayColumns` | `arr` | numpy.ndarray (3×3): Matrix with unit-le |
| 47 | `np_kronecker` | `i, j` | int: Kronecker delta value |
| 48 | `np_permut_tensor3` | `i, j, k` | int: Permutation symbol value |
| 49 | `permut_tensor3` | `i, j, k` | int: +1 for even permutation, -1 for odd |
| 50 | `perpendicular_vector` | `v` | numpy.ndarray [3]: Unit vector perpendic |
| 51 | `plane2string` | `v, digits=2` | str: Formatted string '(h.hh,k.kk,l.ll)' |
| 52 | `plane_line_intersection` | `plane_point, plane_normal, lin...` | numpy.ndarray [3] or None: Intersection  |
| 53 | `print_correspondence` | `C, phase1='Phase1', phase2='Ph...` | None (prints to console) |
| 54 | `read_txt` | `filename` | numpy.ndarray: Loaded data |
| 55 | `reciprocal_basis` | `a1, a2, a3` | b1: numpy array (3,) - First reciprocal  |
| 56 | `select_atomic_plane` | `atoms, plane_normal, plane_poi...` | numpy.ndarray: Indices of atoms on plane |
| 57 | `select_atomic_region` | `atoms, center, radius` | numpy.ndarray: Indices of atoms in regio |
| 58 | `select_crystal_planes` | `lattice, max_index=3` | list: List of (h,k,l) tuples |
| 59 | `select_plane` | `points, h, k, l, L, tolerance=...` | numpy.ndarray: Indices of points on plan |
| 60 | `symmetry_elements` | `crystal_system` | list: List of numpy.ndarray (3×3) symmet |
| 61 | `tetragonal_lattice_vec` | `a, b, c` | numpy.ndarray (3×3): Tetragonal lattice  |
| 62 | `twin_equation_solution` | `L_parent, L_twin, params=None` | dict: { |
| 63 | `twin_equation_solution_ini` | `` | dict: Initial parameters for twin solver |
| 64 | `twinnedhabitplane` | `L_parent, L_twin, corresponden...` | dict: { |
| 65 | `uvtw2uvw` | `uvtw` | numpy.ndarray [3]: 3-index direction [u, |
| 66 | `uvw2uvtw` | `uvw` | numpy.ndarray [4]: 4-index direction [u, |
| 67 | `vec2string` | `v, digits=2` | str: Formatted string '[x.xx,y.yy,z.zz]' |
| 68 | `vector2miller_ini` | `v, L` | numpy.ndarray [3]: Initial Miller index  |
| 69 | `vectors2miller` | `v, L, max_index=5` | numpy.ndarray [3]: Miller indices [uvw] |
| 70 | `write_lattice_correspondence` | `filename, C, L1, L2, phase1='P...` | None (writes to file) |
| 71 | `write_txt` | `filename, data, header=''` | None (writes to file) |
| 72 | `xyz2fractional` | `Txyz2uvw, V, frac=10, eps2=1e-...` | numpy.ndarray [3]: Reduced Miller indice |
| 73 | `xyz2fractional02` | `Txyz2uvw, V` | numpy.ndarray [3]: Normalized fractional |

---

## 📝 FUNCTION SUMMARIES

### `B19p_B2_lattice_correspondence`
Generate B19'→B2 lattice correspondence matrix for NiTi.

### `B19p_B2_lattice_correspondence_ini`
Initialize B19'→B2 correspondence for specific transformation variant.

### `Rp_B2_lattice_correspondence`
Generate R-phase→B2 lattice correspondence for NiTi.

### `an_between_vecs`
Calculate angle between two vectors.

### `array2tuple`
Convert a numpy array to a tuple with rounded elements.

### `cubic2tetragonal_lattice_correspondence`
Generate cubic→tetragonal lattice correspondence.

### `cubic_lattice_vec`
Generate cubic lattice vectors.

### `def_gradient`
Calculate total deformation gradient with rotation.

### `def_gradient_ini`
Initialize deformation gradient with Euler angle rotation.

### `def_gradient_ini2`
Initialize deformation gradient with axis-angle rotation.

### `def_gradient_stressfree`
Calculate stress-free transformation strain (deformation gradient).

### `def_gradient_stressfree_ini`
Initialize stress-free deformation gradient from lattice parameters.

### `dir2string`
Format direction as string representation [u,v,w] with Miller notation.

### `equivalent_elements`
Generate all symmetry-equivalent crystallographic elements.

### `find_gcd`
Find greatest common divisor using recursive Euclidean algorithm.

### `flipvector`
Flip vector to ensure positive first non-zero component.

### `flipvector2negative`
Flip vector to ensure negative first non-zero component.

### `gen_twinned_lattice_points`
Generate lattice points for parent and twinned regions.

### `genallHexSys`
Generate all hexagonal slip systems including <c+a> pyramidal.

### `generate_hkls`
Generate unique Miller indices (hkl) considering crystal symmetry operations.

### `generate_hkls01`
Generate unique Miller indices (hkl) considering crystal symmetry operations.

### `generate_hkls02`
Generate unique Miller indices (hkl) with metric tensor transformation.

### `generate_lattice_faces`
Generate face polygons for unit cell visualization.

### `generate_lattice_points`
Generate lattice points within specified unit cell range.

### `generate_lattice_vectors`
Generate lattice vectors for visualization.

### `generate_lattite_atom_positions`
Generate atomic positions in lattice with basis.

### `generate_plane_vertices`
Generate vertices for plotting crystallographic plane.

### `generate_product_lattice_faces`
Generate faces for two lattices.

### `generate_product_lattice_points`
Generate lattice points for two overlapping lattices.

### `gensystemsHex`
Generate hexagonal slip systems in real space.

### `gensystemsHexIni`
Generate initial hexagonal slip system templates.

### `get_interface2d`
Extract interface atoms from two lattices.

### `get_twinning_dislocation`
Calculate twinning dislocation Burgers vector.

### `get_twinning_plane_points`
Get atoms on twinning plane.

### `get_twinningdata`
Extract complete twinning data from parent and twin lattices.

### `get_unique_families`
Returns unique families of Miller indices based on permutation symmetry.

### `habitplane_equation_solution`
Solve habit plane equation for phase transformation.

### `hkil2hkl`
Convert 4-index hexagonal plane (hkil) to 3-index (hkl).

### `hkl2hkil`
Convert 3-index hexagonal plane (hkl) to 4-index (hkil).

### `kronecker`
Compute Kronecker delta δ_{ij}.

### `lattice_correspondence`
Calculate lattice correspondence matrix between two crystal structures.

### `lattice_vec`
Calculate lattice vectors (a1, a2, a3) for various crystal systems.

### `miller2fractional`
Reduce Miller indices to lowest integer form.

### `monoclinic_lattice_vec`
Generate monoclinic lattice vectors (B19' martensite).

### `niti_twinning`
Get NiTi twinning parameters for specific variant.

### `normArrayColumns`
Normalize each column of matrix to unit length.

### `np_kronecker`
Compute Kronecker delta (NumPy version).

### `np_permut_tensor3`
Compute Levi-Civita symbol (NumPy-compatible version).

### `permut_tensor3`
Compute 3D Levi-Civita permutation symbol ε_{ijk}.

### `perpendicular_vector`
Find arbitrary unit vector perpendicular to input 3D vector.

### `plane2string`
Format plane normal as string representation (h,k,l) with Miller notation.

### `plane_line_intersection`
Calculate intersection of plane and line in 3D.

### `print_correspondence`
Print lattice correspondence matrix in readable format.

### `read_txt`
Read numerical data from text file.

### `reciprocal_basis`
Calculate reciprocal lattice basis vectors from real space lattice vectors.

### `select_atomic_plane`
Select atoms lying on or near a crystallographic plane.

### `select_atomic_region`
Select atoms within spherical region.

### `select_crystal_planes`
Generate list of low-index crystallographic planes.

### `select_plane`
Select points on crystallographic plane (hkl).

### `symmetry_elements`
Generate symmetry operation matrices for crystal system.

### `tetragonal_lattice_vec`
Generate tetragonal lattice vectors.

### `twin_equation_solution`
Solve complete twin equation system.

### `twin_equation_solution_ini`
Initialize twin equation solver with default parameters.

### `twinnedhabitplane`
Calculate twinned habit plane from lattice parameters.

### `uvtw2uvw`
Convert 4-index hexagonal direction [uvtw] to 3-index [uvw].

### `uvw2uvtw`
Convert 3-index hexagonal direction [uvw] to 4-index [uvtw].

### `vec2string`
Format vector as string representation [x,y,z] with specified precision.

### `vector2miller_ini`
Convert Cartesian vector to Miller indices (initial guess).

### `vectors2miller`
Convert Cartesian vector to Miller indices with optimization.

### `write_lattice_correspondence`
Write lattice correspondence analysis to file.

### `write_txt`
Write numerical data to text file.

### `xyz2fractional`
Convert Cartesian coordinates to fractional Miller indices with reduction.

### `xyz2fractional02`
Simple Cartesian to fractional coordinate transformation without reduction.


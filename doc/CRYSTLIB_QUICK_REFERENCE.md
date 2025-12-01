# CRYSTLIB - QUICK REFERENCE
## Ultra-Fast Lookup

**Module**: `crystlib.py` (Extended Version)  
**Functions**: 73 | **Type**: Condensed Reference

---

## ⚡ FUNCTION LIST

**`B19p_B2_lattice_correspondence`** - Generate B19'→B2 lattice correspondence matrix for NiTi.

**`B19p_B2_lattice_correspondence_ini`** - Initialize B19'→B2 correspondence for specific transformatio

**`Rp_B2_lattice_correspondence`** - Generate R-phase→B2 lattice correspondence for NiTi.

**`an_between_vecs`** - Calculate angle between two vectors.

**`array2tuple`** - Convert a numpy array to a tuple with rounded elements.

**`cubic2tetragonal_lattice_correspondence`** - Generate cubic→tetragonal lattice correspondence.

**`cubic_lattice_vec`** - Generate cubic lattice vectors.

**`def_gradient`** - Calculate total deformation gradient with rotation.

**`def_gradient_ini`** - Initialize deformation gradient with Euler angle rotation.

**`def_gradient_ini2`** - Initialize deformation gradient with axis-angle rotation.

**`def_gradient_stressfree`** - Calculate stress-free transformation strain (deformation gra

**`def_gradient_stressfree_ini`** - Initialize stress-free deformation gradient from lattice par

**`dir2string`** - Format direction as string representation [u,v,w] with Mille

**`equivalent_elements`** - Generate all symmetry-equivalent crystallographic elements.

**`find_gcd`** - Find greatest common divisor using recursive Euclidean algor

**`flipvector`** - Flip vector to ensure positive first non-zero component.

**`flipvector2negative`** - Flip vector to ensure negative first non-zero component.

**`gen_twinned_lattice_points`** - Generate lattice points for parent and twinned regions.

**`genallHexSys`** - Generate all hexagonal slip systems including <c+a> pyramida

**`generate_hkls`** - Generate unique Miller indices (hkl) considering crystal sym

**`generate_hkls01`** - Generate unique Miller indices (hkl) considering crystal sym

**`generate_hkls02`** - Generate unique Miller indices (hkl) with metric tensor tran

**`generate_lattice_faces`** - Generate face polygons for unit cell visualization.

**`generate_lattice_points`** - Generate lattice points within specified unit cell range.

**`generate_lattice_vectors`** - Generate lattice vectors for visualization.

**`generate_lattite_atom_positions`** - Generate atomic positions in lattice with basis.

**`generate_plane_vertices`** - Generate vertices for plotting crystallographic plane.

**`generate_product_lattice_faces`** - Generate faces for two lattices.

**`generate_product_lattice_points`** - Generate lattice points for two overlapping lattices.

**`gensystemsHex`** - Generate hexagonal slip systems in real space.

**`gensystemsHexIni`** - Generate initial hexagonal slip system templates.

**`get_interface2d`** - Extract interface atoms from two lattices.

**`get_twinning_dislocation`** - Calculate twinning dislocation Burgers vector.

**`get_twinning_plane_points`** - Get atoms on twinning plane.

**`get_twinningdata`** - Extract complete twinning data from parent and twin lattices

**`get_unique_families`** - Returns unique families of Miller indices based on permutati

**`habitplane_equation_solution`** - Solve habit plane equation for phase transformation.

**`hkil2hkl`** - Convert 4-index hexagonal plane (hkil) to 3-index (hkl).

**`hkl2hkil`** - Convert 3-index hexagonal plane (hkl) to 4-index (hkil).

**`kronecker`** - Compute Kronecker delta δ_{ij}.

**`lattice_correspondence`** - Calculate lattice correspondence matrix between two crystal 

**`lattice_vec`** - Calculate lattice vectors (a1, a2, a3) for various crystal s

**`miller2fractional`** - Reduce Miller indices to lowest integer form.

**`monoclinic_lattice_vec`** - Generate monoclinic lattice vectors (B19' martensite).

**`niti_twinning`** - Get NiTi twinning parameters for specific variant.

**`normArrayColumns`** - Normalize each column of matrix to unit length.

**`np_kronecker`** - Compute Kronecker delta (NumPy version).

**`np_permut_tensor3`** - Compute Levi-Civita symbol (NumPy-compatible version).

**`permut_tensor3`** - Compute 3D Levi-Civita permutation symbol ε_{ijk}.

**`perpendicular_vector`** - Find arbitrary unit vector perpendicular to input 3D vector.

**`plane2string`** - Format plane normal as string representation (h,k,l) with Mi

**`plane_line_intersection`** - Calculate intersection of plane and line in 3D.

**`print_correspondence`** - Print lattice correspondence matrix in readable format.

**`read_txt`** - Read numerical data from text file.

**`reciprocal_basis`** - Calculate reciprocal lattice basis vectors from real space l

**`select_atomic_plane`** - Select atoms lying on or near a crystallographic plane.

**`select_atomic_region`** - Select atoms within spherical region.

**`select_crystal_planes`** - Generate list of low-index crystallographic planes.

**`select_plane`** - Select points on crystallographic plane (hkl).

**`symmetry_elements`** - Generate symmetry operation matrices for crystal system.

**`tetragonal_lattice_vec`** - Generate tetragonal lattice vectors.

**`twin_equation_solution`** - Solve complete twin equation system.

**`twin_equation_solution_ini`** - Initialize twin equation solver with default parameters.

**`twinnedhabitplane`** - Calculate twinned habit plane from lattice parameters.

**`uvtw2uvw`** - Convert 4-index hexagonal direction [uvtw] to 3-index [uvw].

**`uvw2uvtw`** - Convert 3-index hexagonal direction [uvw] to 4-index [uvtw].

**`vec2string`** - Format vector as string representation [x,y,z] with specifie

**`vector2miller_ini`** - Convert Cartesian vector to Miller indices (initial guess).

**`vectors2miller`** - Convert Cartesian vector to Miller indices with optimization

**`write_lattice_correspondence`** - Write lattice correspondence analysis to file.

**`write_txt`** - Write numerical data to text file.

**`xyz2fractional`** - Convert Cartesian coordinates to fractional Miller indices w

**`xyz2fractional02`** - Simple Cartesian to fractional coordinate transformation wit


# COMPLETE FUNCTION MAPPING REFERENCE
## From crystallography_functions.py to Library Files

This document provides the complete mapping of all 116 functions showing where each function was distributed.

---

## 📋 PART 1 FUNCTIONS (1-40)

| # | Function Name | Source | Destination | Status |
|---|---------------|--------|-------------|--------|
| 1 | `set_aspect_equal_3d` | PART1 | plotlib | ✓ Added |
| 2 | `find_gcd` | PART1 | crystlib | ✓ Added |
| 3 | `perpendicular_vector` | PART1 | crystlib | ✓ Added |
| 4 | `vec2string` | PART1 | crystlib | ✓ Added |
| 5 | `plane2string` | PART1 | crystlib | ✓ Added |
| 6 | `dir2string` | PART1 | crystlib | ✓ Added |
| 7 | `xyz2fractional` | PART1 | crystlib | ✓ Added |
| 8 | `miller2fractional` | PART1 | crystlib | ✓ Added |
| 9 | `xyz2fractional02` | PART1 | crystlib | ✓ Added |
| 10 | `normArrayColumns` | PART1 | crystlib | ✓ Added |
| 11 | `cubic_lattice_vec` | PART1 | crystlib | ✓ Added |
| 12 | `monoclinic_lattice_vec` | PART1 | crystlib | ✓ Added |
| 13 | `tetragonal_lattice_vec` | PART1 | crystlib | ✓ Added |
| 14 | `uvtw2uvw` | PART1 | crystlib | ✓ Added |
| 15 | `uvw2uvtw` | PART1 | crystlib | ✓ Added |
| 16 | `hkil2hkl` | PART1 | crystlib | ✓ Added |
| 17 | `hkl2hkil` | PART1 | crystlib | ✓ Added |
| 18 | `gensystemsHexIni` | PART1 | crystlib | ✓ Added |
| 19 | `gensystemsHex` | PART1 | crystlib | ✓ Added |
| 20 | `rotation_from_axis_angle` | PART1 | orilib | ✓ Added |
| 21 | `genallHexSys` | PART1 | crystlib | ✓ Added |
| 22 | `lattice_vec` | PART1 | crystlib | ⚠ Skipped (exists) |
| 23 | `reciprocal_basis` | PART1 | crystlib | ⚠ Skipped (exists) |
| 24 | `np_euler_matrix` | PART1 | orilib | ⚠ Skipped (exists) |
| 25 | `np_inverse_euler_matrix` | PART1 | orilib | ⚠ Skipped (exists) |
| 26 | `ol_g_rtheta_rad` | PART1 | orilib | ⚠ Skipped (exists) |
| 27 | `np_ol_g_rtheta_rad` | PART1 | orilib | ⚠ Skipped (exists) |
| 28 | `ol_rtheta_g_rad` | PART1 | orilib | ⚠ Skipped (exists) |
| 29 | `np_ol_rtheta_g_rad` | PART1 | orilib | ⚠ Skipped (exists) |
| 30 | `ol_g_R` | PART1 | orilib | ✓ Added |
| 31 | `np_ol_g_R` | PART1 | orilib | ✓ Added |
| 32 | `ol_R_g` | PART1 | orilib | ✓ Added |
| 33 | `np_ol_R_g` | PART1 | orilib | ✓ Added |
| 34 | `ol_g_R2` | PART1 | orilib | ✓ Added |
| 35 | `np_ol_g_R2` | PART1 | orilib | ✓ Added |
| 36 | `ol_R_g2` | PART1 | orilib | ✓ Added |
| 37 | `np_ol_R_g2` | PART1 | orilib | ✓ Added |
| 38 | `np_ol_R_q2` | PART1 | orilib | ✓ Added |
| 39 | `np_ol_g_q2` | PART1 | orilib | ✓ Added |
| 40 | `np_ol_q_g` | PART1 | orilib | ✓ Added |

---

## 📋 PART 2 FUNCTIONS (41-80)

| # | Function Name | Source | Destination | Status |
|---|---------------|--------|-------------|--------|
| 41 | `permut_tensor3` | PART2 | crystlib | ✓ Added |
| 42 | `np_permut_tensor3` | PART2 | crystlib | ✓ Added |
| 43 | `kronecker` | PART2 | crystlib | ✓ Added |
| 44 | `np_kronecker` | PART2 | crystlib | ✓ Added |
| 45 | `active_rotation` | PART2 | orilib | ✓ Added |
| 46 | `passive_rotation` | PART2 | orilib | ✓ Added |
| 47 | `stereoprojection_directions` | PART2 | projlib | ⚠ Skipped (exists) |
| 48 | `equalarea_directions` | PART2 | projlib | ⚠ Skipped (exists) |
| 49 | `euler_angles_reduction` | PART2 | orilib | ✓ Added |
| 50 | `symmetry_elements` | PART2 | crystlib | ✓ Added |
| 51 | `equivalent_elements` | PART2 | crystlib | ✓ Added |
| 52 | `B19p_B2_lattice_correspondence` | PART2 | crystlib | ✓ Added |
| 53 | `lattice_correspondence` | PART2 | crystlib | ✓ Added |
| 54 | `B19p_B2_lattice_correspondence_ini` | PART2 | crystlib | ✓ Added |
| 55 | `cubic2tetragonal_lattice_correspondence` | PART2 | crystlib | ✓ Added |
| 56 | `Rp_B2_lattice_correspondence` | PART2 | crystlib | ✓ Added |
| 57 | `print_correspondence` | PART2 | crystlib | ✓ Added |
| 58 | `mohr_circles` | PART2 | plotlib | ✓ Added |
| 59 | `generate_lattice_points` | PART2 | crystlib | ✓ Added |
| 60 | `plot_lattice_plane` | PART2 | plotlib | ✓ Added |
| 61 | `plot_lattice_boundaries` | PART2 | plotlib | ✓ Added |
| 62 | `generate_lattice_faces` | PART2 | crystlib | ✓ Added |
| 63 | `generate_product_lattice_points` | PART2 | crystlib | ✓ Added |
| 64 | `generate_product_lattice_faces` | PART2 | crystlib | ✓ Added |
| 65 | `plot_lattice3D` | PART2 | plotlib | ✓ Added |
| 66 | `plot_latticefaces3D` | PART2 | plotlib | ✓ Added |
| 67 | `plot_latticesfaces3D` | PART2 | plotlib | ✓ Added |
| 68 | `plot_lattice2D` | PART2 | plotlib | ✓ Added |
| 69 | `plot_lattice_2Dprojection` | PART2 | plotlib | ✓ Added |
| 70 | `zero_normal_strains` | PART2 | plotlib | ✓ Added |
| 71 | `strains_along_13mohrcirle` | PART2 | plotlib | ✓ Added |
| 72 | `select_crystal_planes` | PART2 | crystlib | ✓ Added |
| 73 | `plot_mohr_circles` | PART2 | plotlib | ✓ Added |
| 74 | `plot_planes_on_mohr_circle` | PART2 | plotlib | ✓ Added |
| 75 | `plot_planes_on_stereotriangle` | PART2 | plotlib | ✓ Added |
| 76 | `plot_planes_on_wulffnet` | PART2 | plotlib | ✓ Added |
| 77 | `plot_princip_dir_on_stereotriangle` | PART2 | plotlib | ✓ Added |
| 78 | `plot_princip_dir_on_wulffnet` | PART2 | plotlib | ✓ Added |
| 79 | `write_mohr_planes` | PART2 | plotlib | ✓ Added |
| 80 | `write_lattice_correspondence` | PART2 | crystlib | ✓ Added |

---

## 📋 PART 3 FUNCTIONS (81-116)

| # | Function Name | Source | Destination | Status |
|---|---------------|--------|-------------|--------|
| 81 | `generate_lattite_atom_positions` | PART3 | crystlib | ✓ Added |
| 82 | `generate_lattice_vectors` | PART3 | crystlib | ✓ Added |
| 83 | `plot_lattice` | PART3 | plotlib | ✓ Added |
| 84 | `plot_lattice_proj` | PART3 | plotlib | ✓ Added |
| 85 | `plot_points_proj` | PART3 | plotlib | ✓ Added |
| 86 | `select_atomic_plane` | PART3 | crystlib | ✓ Added |
| 87 | `get_interface2d` | PART3 | crystlib | ✓ Added |
| 88 | `select_plane` | PART3 | crystlib | ✓ Added |
| 89 | `generate_plane_vertices` | PART3 | crystlib | ✓ Added |
| 90 | `select_atomic_region` | PART3 | crystlib | ✓ Added |
| 91 | `plot_atomic_plane2D` | PART3 | plotlib | ✓ Added |
| 92 | `get_twinning_plane_points` | PART3 | crystlib | ✓ Added |
| 93 | `plot_atomic_plane3D` | PART3 | plotlib | ✓ Added |
| 94 | `plot_atomlattice2D` | PART3 | plotlib | ✓ Added |
| 95 | `an_between_vecs` | PART3 | crystlib | ✓ Added |
| 96 | `habitplane_equation_solution` | PART3 | crystlib | ✓ Added |
| 97 | `twinnedhabitplane` | PART3 | crystlib | ✓ Added |
| 98 | `twin_equation_solution_ini` | PART3 | crystlib | ✓ Added |
| 99 | `twin_equation_solution` | PART3 | crystlib | ✓ Added |
| 100 | `def_gradient_stressfree` | PART3 | crystlib | ✓ Added |
| 101 | `def_gradient_stressfree_ini` | PART3 | crystlib | ✓ Added |
| 102 | `def_gradient` | PART3 | crystlib | ✓ Added |
| 103 | `def_gradient_ini` | PART3 | crystlib | ✓ Added |
| 104 | `def_gradient_ini2` | PART3 | crystlib | ✓ Added |
| 105 | `niti_twinning` | PART3 | crystlib | ✓ Added |
| 106 | `get_twinningdata` | PART3 | crystlib | ✓ Added |
| 107 | `get_twinning_dislocation` | PART3 | crystlib | ✓ Added |
| 108 | `gen_twinned_lattice_points` | PART3 | crystlib | ✓ Added |
| 109 | `write_txt` | PART3 | crystlib | ✓ Added |
| 110 | `read_txt` | PART3 | crystlib | ✓ Added |
| 111 | `plane_line_intersection` | PART3 | crystlib | ✓ Added |
| 112 | `plot_cut2D` | PART3 | plotlib | ✓ Added |
| 113 | `flipvector` | PART3 | crystlib | ✓ Added |
| 114 | `flipvector2negative` | PART3 | crystlib | ✓ Added |
| 115 | `vector2miller_ini` | PART3 | crystlib | ✓ Added |
| 116 | `vectors2miller` | PART3 | crystlib | ✓ Added |

---

## 📊 SUMMARY BY DESTINATION

### CRYSTLIB (68 functions assigned, 66 added, 2 skipped)

**Utility Functions (8):**
- find_gcd, perpendicular_vector
- vec2string, plane2string, dir2string
- flipvector, flipvector2negative
- normArrayColumns

**Lattice Operations (16):**
- cubic_lattice_vec, monoclinic_lattice_vec, tetragonal_lattice_vec
- lattice_vec (skipped), reciprocal_basis (skipped)
- generate_lattice_points, generate_lattice_faces
- generate_product_lattice_points, generate_product_lattice_faces
- generate_lattite_atom_positions, generate_lattice_vectors
- generate_plane_vertices

**Miller Indices (8):**
- uvtw2uvw, uvw2uvtw, hkil2hkl, hkl2hkil
- vector2miller_ini, vectors2miller

**Coordinate Transformations (3):**
- xyz2fractional, miller2fractional, xyz2fractional02

**Tensor Operations (4):**
- permut_tensor3, np_permut_tensor3
- kronecker, np_kronecker

**Symmetry (2):**
- symmetry_elements, equivalent_elements

**Lattice Correspondence (7):**
- B19p_B2_lattice_correspondence, lattice_correspondence
- B19p_B2_lattice_correspondence_ini
- cubic2tetragonal_lattice_correspondence
- Rp_B2_lattice_correspondence
- print_correspondence, write_lattice_correspondence

**Hexagonal Systems (3):**
- gensystemsHexIni, gensystemsHex, genallHexSys

**Plane Selection (4):**
- select_atomic_plane, select_plane
- select_atomic_region, select_crystal_planes

**Twinning (7):**
- an_between_vecs, habitplane_equation_solution
- twinnedhabitplane, twin_equation_solution_ini
- twin_equation_solution, get_twinning_plane_points
- get_interface2d

**Deformation Gradients (5):**
- def_gradient_stressfree, def_gradient_stressfree_ini
- def_gradient, def_gradient_ini, def_gradient_ini2

**NiTi Specific (4):**
- niti_twinning, get_twinningdata
- get_twinning_dislocation, gen_twinned_lattice_points

**File I/O (2):**
- write_txt, read_txt

**Other (2):**
- plane_line_intersection

### ORILIB (21 functions assigned, 15 added, 6 skipped)

**Euler Angles (3, 1 new):**
- np_euler_matrix (skipped)
- np_inverse_euler_matrix (skipped)
- euler_angles_reduction (added)

**Axis-Angle (7, 1 new):**
- ol_g_rtheta_rad (skipped), np_ol_g_rtheta_rad (skipped)
- ol_rtheta_g_rad (skipped), np_ol_rtheta_g_rad (skipped)
- rotation_from_axis_angle (added)

**Rodrigues-Frank Vectors (10 new):**
- ol_g_R, np_ol_g_R, ol_R_g, np_ol_R_g
- ol_g_R2, np_ol_g_R2, ol_R_g2, np_ol_R_g2
- np_ol_R_q2, np_ol_g_q2, np_ol_q_g

**Rotation Operations (2 new):**
- active_rotation, passive_rotation

### PLOTLIB (25 functions, all added)

**Mohr Circles (5):**
- mohr_circles, plot_mohr_circles
- plot_planes_on_mohr_circle
- strains_along_13mohrcirle, zero_normal_strains
- write_mohr_planes

**3D Plotting (4):**
- set_aspect_equal_3d, plot_lattice3D
- plot_latticefaces3D, plot_latticesfaces3D

**2D Plotting (7):**
- plot_lattice2D, plot_lattice_2Dprojection
- plot_lattice, plot_lattice_proj
- plot_points_proj, plot_cut2D
- plot_lattice_plane, plot_lattice_boundaries

**Atomic Planes (3):**
- plot_atomic_plane2D, plot_atomic_plane3D
- plot_atomlattice2D

**Stereographic Plotting (4):**
- plot_planes_on_stereotriangle
- plot_planes_on_wulffnet
- plot_princip_dir_on_stereotriangle
- plot_princip_dir_on_wulffnet

### PROJLIB (2 functions assigned, 0 added, 2 skipped)

**Stereographic Projections (2, both skipped):**
- stereoprojection_directions (skipped - already exists)
- equalarea_directions (skipped - already exists)

---

## 🎯 DISTRIBUTION RATIONALE

### Why CRYSTLIB got the most functions (66 added)

CRYSTLIB is the core crystallography library handling:
- Fundamental lattice operations
- All coordinate transformations
- Miller index operations
- Symmetry operations
- Material-specific calculations (NiTi twinning)
- Deformation analysis
- File I/O for crystallographic data

### Why ORILIB got moderate additions (15 added)

ORILIB focuses specifically on orientation representations:
- Already had comprehensive quaternion operations
- Needed extended Rodrigues representations
- Needed rotation operation utilities
- Euler angle functions already present

### Why PLOTLIB got significant additions (25 added)

PLOTLIB had minimal original functions but needed:
- Complete Mohr circle visualization
- Lattice visualization in 2D and 3D
- Atomic plane plotting
- Stereographic projection plotting
- All graphical output functions

### Why PROJLIB got no additions (0 added)

PROJLIB already comprehensive:
- Had all stereographic projection functions
- The 2 functions from crystallography_functions already existed
- No gaps in functionality

---

## ✅ VERIFICATION CHECKLIST

- [x] All 116 functions accounted for
- [x] No functions lost during distribution
- [x] Existing functions preserved (10 skipped appropriately)
- [x] Complete docstrings transferred
- [x] File headers updated
- [x] Import statements checked
- [x] 106 functions successfully distributed
- [x] 10 functions appropriately skipped (already existed)
- [x] 0 functions missing

---

**Document Created**: December 2024  
**Status**: Complete and Verified  
**Total Functions**: 116 mapped  
**Successfully Distributed**: 106  
**Appropriately Skipped**: 10  


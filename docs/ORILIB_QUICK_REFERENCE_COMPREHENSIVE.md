# ORILIB - QUICK REFERENCE (Comprehensive)
## Fast Function Lookup

**Module**: `orilib.py` (Extended Version)  
**Total Functions**: 43  

---

## 🔍 QUICK LOOKUP TABLE

| # | Function | Parameters | Returns |
|---|----------|------------|---------|
| 1 | `QMatproduct` | `sym, Q` | SQ: numpy array (4, N) - Transformed qua |
| 2 | `Qlog` | `QM` | qlog: numpy array (4, N, M) - Logarithm  |
| 3 | `Qproduct` | `P, Q` | result: numpy array (4, N) - Product qua |
| 4 | `active_rotation` | `g, v` | numpy.ndarray [3]: Rotated vector v' |
| 5 | `eu2quat` | `phi1, Phi, phi2` | q: numpy array (4,) - Quaternion [w, x,  |
| 6 | `euler_angles_reduction` | `phi1, Phi, phi2, symops` | tuple: (phi1_red, Phi_red, phi2_red) - R |
| 7 | `grid_s1` | `resol, grids=6` | points: list of floats - Angles in radia |
| 8 | `hopf2quat` | `Points` | quats: list of lists - Quaternions [w, x |
| 9 | `mat_to_quat` | `R` | q: numpy array (4,) - Unit quaternion [w |
| 10 | `misori_sym_deg_quats` | `q1, q2, sym_quats` | min_angle: float - Minimum misorientatio |
| 11 | `misori_sym_deg_sample_to_crystal_fast` | `M1, M2, symops` | miso: numpy array (N,) - Minimum misorie |
| 12 | `mk_pix2xy` | `` | pix2x: list (1024,) - X-coordinate looku |
| 13 | `np_euler_matrix` | `ai, aj, ak` | g: numpy array (3, 3) - Rotation matrix |
| 14 | `np_eulers_matrices` | `data, deg=False` | g: numpy array (N, 3, 3) - Array of rota |
| 15 | `np_g2quats` | `umatsa` | Q: numpy array (4, N) - Array of quatern |
| 16 | `np_gmat2rodrigues` | `g` | rodrigues: numpy array (3,) - Rodrigues- |
| 17 | `np_inverse_euler_matrix` | `ai, aj, ak` | U: numpy array (3, 3) - Inverse rotation |
| 18 | `np_ol_R_g` | `R` | numpy.ndarray (3×3): Rotation matrix |
| 19 | `np_ol_R_g2` | `R` | numpy.ndarray (3×3): Rotation matrix |
| 20 | `np_ol_R_q2` | `R` | numpy.ndarray [4]: Unit quaternion [w, x |
| 21 | `np_ol_g_R` | `g` | numpy.ndarray [3]: Rodrigues-Frank vecto |
| 22 | `np_ol_g_R2` | `g` | numpy.ndarray [3]: Rodrigues-Frank vecto |
| 23 | `np_ol_g_q2` | `g` | numpy.ndarray [4]: Unit quaternion [w, x |
| 24 | `np_ol_g_rtheta_rad` | `g` | r: numpy array (3,) - Rotation axis (uni |
| 25 | `np_ol_q_g` | `q` | numpy.ndarray (3×3): Rotation matrix |
| 26 | `np_ol_rtheta_g_rad` | `r, theta` | g: numpy array (3, 3) - Rotation matrix |
| 27 | `np_rodrigues2gmat` | `rodrigues` | g: numpy array (3, 3) - Rotation matrix |
| 28 | `nside2npix` | `nside` | npix: int - Total number of pixels = 12  |
| 29 | `ol_R_g` | `R` | list (3×3): Rotation matrix |
| 30 | `ol_R_g2` | `R` | list (3×3): Rotation matrix |
| 31 | `ol_g_R` | `g` | list [3]: Rodrigues-Frank vector |
| 32 | `ol_g_R2` | `g` | list [3]: Rodrigues-Frank vector |
| 33 | `ol_g_rtheta_rad` | `g` | r: list (3,) - Rotation axis (unit vecto |
| 34 | `ol_rtheta_g_rad` | `r, theta` | g: list (3, 3) - Rotation matrix |
| 35 | `passive_rotation` | `g, v` | numpy.ndarray [3]: Vector components in  |
| 36 | `pix2ang_nest` | `nside, ipix, pix2x, pix2y` | theta: float - Colatitude angle in radia |
| 37 | `quat_conjugate` | `q` | q_conj: numpy array (4,) - Conjugate qua |
| 38 | `quat_misori_deg` | `q1, q2` | angle: float - Misorientation angle in d |
| 39 | `quat_mult` | `q1, q2` | q: numpy array (4,) - Product quaternion |
| 40 | `quat_multiply` | `q1, q2` | q: numpy array (4,) - Product quaternion |
| 41 | `quat_to_mat` | `q` | R: numpy array (3x3) - Rotation matrix |
| 42 | `rotation_from_axis_angle` | `axis, angle` | numpy.ndarray (3×3): Rotation matrix R |
| 43 | `simple_grid` | `resol` | quats: list of lists - Uniformly distrib |

---

## 📝 FUNCTION SUMMARIES

### `QMatproduct`
Multiply a single symmetry quaternion with multiple quaternions.

### `Qlog`
Compute the logarithm of quaternions (quaternion logarithm map).

### `Qproduct`
Compute quaternion product for arrays of quaternions.

### `active_rotation`
Perform active rotation of vector v by rotation matrix g.

### `eu2quat`
Convert Bunge Euler angles to quaternion representation.

### `euler_angles_reduction`
Reduce Euler angles to fundamental zone using symmetry operations.

### `grid_s1`
Generate uniformly distributed points on S¹ (circle).

### `hopf2quat`
Convert Hopf coordinates to quaternions.

### `mat_to_quat`
Convert a rotation matrix to a quaternion representation.

### `misori_sym_deg_quats`
Calculate minimum misorientation angle between two quaternions considering

### `misori_sym_deg_sample_to_crystal_fast`
Compute misorientation angles (deg) between orientations M1 and M2

### `mk_pix2xy`
Create lookup tables for HEALPix pixel index to (x,y) coordinate conversion.

### `np_euler_matrix`
Convert Euler angles to rotation matrix (ZXZ convention, Bunge notation).

### `np_eulers_matrices`
Convert multiple Euler angles to rotation matrices (vectorized).

### `np_g2quats`
Convert multiple rotation matrices to quaternions (vectorized).

### `np_gmat2rodrigues`
Convert rotation matrix to Rodrigues-Frank vector representation.

### `np_inverse_euler_matrix`
Convert Euler angles to inverse (transpose) rotation matrix.

### `np_ol_R_g`
Convert Rodrigues-Frank vector to rotation matrix (NumPy version).

### `np_ol_R_g2`
Alternative Rodrigues-Frank to rotation matrix (NumPy version, method 2).

### `np_ol_R_q2`
Convert Rodrigues-Frank vector to quaternion (method 2).

### `np_ol_g_R`
Convert rotation matrix to Rodrigues-Frank vector (NumPy version).

### `np_ol_g_R2`
Alternative Rodrigues-Frank conversion (NumPy version, method 2).

### `np_ol_g_q2`
Convert rotation matrix to quaternion (method 2).

### `np_ol_g_rtheta_rad`
Convert rotation matrix to axis-angle representation (NumPy optimized version).

### `np_ol_q_g`
Convert quaternion to rotation matrix.

### `np_ol_rtheta_g_rad`
Convert axis-angle representation to rotation matrix (NumPy version).

### `np_rodrigues2gmat`
Convert Rodrigues-Frank vector to rotation matrix.

### `nside2npix`
Calculate the number of pixels in a HEALPix map.

### `ol_R_g`
Convert Rodrigues-Frank vector to rotation matrix (list version).

### `ol_R_g2`
Alternative Rodrigues-Frank to rotation matrix (list version, method 2).

### `ol_g_R`
Convert rotation matrix to Rodrigues-Frank vector (list version).

### `ol_g_R2`
Alternative Rodrigues-Frank conversion (list version, method 2).

### `ol_g_rtheta_rad`
Convert rotation matrix to axis-angle representation (Rodrigues-Frank vector).

### `ol_rtheta_g_rad`
Convert axis-angle representation to rotation matrix using Rodrigues' formula.

### `passive_rotation`
Perform passive rotation (coordinate transformation).

### `pix2ang_nest`
Convert HEALPix pixel index to spherical coordinates (theta, phi).

### `quat_conjugate`
Compute the conjugate of a quaternion (inverse for unit quaternions).

### `quat_misori_deg`
Calculate the misorientation angle between two quaternions in degrees.

### `quat_mult`
Multiply two quaternions (Hamilton product).

### `quat_multiply`
Quaternion multiplication (same as quat_mult, alternative implementation).

### `quat_to_mat`
Convert a quaternion to a rotation matrix.

### `rotation_from_axis_angle`
Generate rotation matrix from axis-angle representation using Rodrigues formula.

### `simple_grid`
Generate a uniform grid of quaternions covering orientation space (SO(3)).


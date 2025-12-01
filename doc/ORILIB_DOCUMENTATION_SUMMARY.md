# ORILIB - DOCUMENTATION SUMMARY
## Organized Function Reference

**Module**: `orilib.py` (Extended Version)  
**Total Functions**: 43  
**Version**: Extended 2024  

---

## 📊 FUNCTION OVERVIEW

### By Alphabetical Order

| Function | Purpose |
|----------|---------|
| `QMatproduct` | Multiply a single symmetry quaternion with multiple quaternions.... |
| `Qlog` | Compute the logarithm of quaternions (quaternion logarithm map).... |
| `Qproduct` | Compute quaternion product for arrays of quaternions.... |
| `active_rotation` | Perform active rotation of vector v by rotation matrix g.... |
| `eu2quat` | Convert Bunge Euler angles to quaternion representation.... |
| `euler_angles_reduction` | Reduce Euler angles to fundamental zone using symmetry operations.... |
| `grid_s1` | Generate uniformly distributed points on S¹ (circle).... |
| `hopf2quat` | Convert Hopf coordinates to quaternions.... |
| `mat_to_quat` | Convert a rotation matrix to a quaternion representation.... |
| `misori_sym_deg_quats` | Calculate minimum misorientation angle between two quaternions considering... |
| `misori_sym_deg_sample_to_crystal_fast` | Compute misorientation angles (deg) between orientations M1 and M2... |
| `mk_pix2xy` | Create lookup tables for HEALPix pixel index to (x,y) coordinate conversion.... |
| `np_euler_matrix` | Convert Euler angles to rotation matrix (ZXZ convention, Bunge notation).... |
| `np_eulers_matrices` | Convert multiple Euler angles to rotation matrices (vectorized).... |
| `np_g2quats` | Convert multiple rotation matrices to quaternions (vectorized).... |
| `np_gmat2rodrigues` | Convert rotation matrix to Rodrigues-Frank vector representation.... |
| `np_inverse_euler_matrix` | Convert Euler angles to inverse (transpose) rotation matrix.... |
| `np_ol_R_g` | Convert Rodrigues-Frank vector to rotation matrix (NumPy version).... |
| `np_ol_R_g2` | Alternative Rodrigues-Frank to rotation matrix (NumPy version, method 2).... |
| `np_ol_R_q2` | Convert Rodrigues-Frank vector to quaternion (method 2).... |
| `np_ol_g_R` | Convert rotation matrix to Rodrigues-Frank vector (NumPy version).... |
| `np_ol_g_R2` | Alternative Rodrigues-Frank conversion (NumPy version, method 2).... |
| `np_ol_g_q2` | Convert rotation matrix to quaternion (method 2).... |
| `np_ol_g_rtheta_rad` | Convert rotation matrix to axis-angle representation (NumPy optimized version).... |
| `np_ol_q_g` | Convert quaternion to rotation matrix.... |
| `np_ol_rtheta_g_rad` | Convert axis-angle representation to rotation matrix (NumPy version).... |
| `np_rodrigues2gmat` | Convert Rodrigues-Frank vector to rotation matrix.... |
| `nside2npix` | Calculate the number of pixels in a HEALPix map.... |
| `ol_R_g` | Convert Rodrigues-Frank vector to rotation matrix (list version).... |
| `ol_R_g2` | Alternative Rodrigues-Frank to rotation matrix (list version, method 2).... |
| `ol_g_R` | Convert rotation matrix to Rodrigues-Frank vector (list version).... |
| `ol_g_R2` | Alternative Rodrigues-Frank conversion (list version, method 2).... |
| `ol_g_rtheta_rad` | Convert rotation matrix to axis-angle representation (Rodrigues-Frank vector).... |
| `ol_rtheta_g_rad` | Convert axis-angle representation to rotation matrix using Rodrigues' formula.... |
| `passive_rotation` | Perform passive rotation (coordinate transformation).... |
| `pix2ang_nest` | Convert HEALPix pixel index to spherical coordinates (theta, phi).... |
| `quat_conjugate` | Compute the conjugate of a quaternion (inverse for unit quaternions).... |
| `quat_misori_deg` | Calculate the misorientation angle between two quaternions in degrees.... |
| `quat_mult` | Multiply two quaternions (Hamilton product).... |
| `quat_multiply` | Quaternion multiplication (same as quat_mult, alternative implementation).... |
| `quat_to_mat` | Convert a quaternion to a rotation matrix.... |
| `rotation_from_axis_angle` | Generate rotation matrix from axis-angle representation using Rodrigues formula.... |
| `simple_grid` | Generate a uniform grid of quaternions covering orientation space (SO(3)).... |

---

## 📖 DETAILED FUNCTION DESCRIPTIONS

### `QMatproduct`

**Signature**: `def QMatproduct(sym, Q)`

Multiply a single symmetry quaternion with multiple quaternions.
    Applies the same symmetry operation to an array of orientations.

---

### `Qlog`

**Signature**: `def Qlog(QM)`

Compute the logarithm of quaternions (quaternion logarithm map).
    Maps quaternions from unit sphere S³ to tangent space at identity.

---

### `Qproduct`

**Signature**: `def Qproduct(P, Q)`

Compute quaternion product for arrays of quaternions.
    Computes P * Q for each pair of quaternions.

---

### `active_rotation`

**Signature**: `def active_rotation(g, v)`

Perform active rotation of vector v by rotation matrix g.
    Rotates the vector itself while keeping coordinate system fixed.
    v' = g · v (matrix-vector multiplication).

---

### `eu2quat`

**Signature**: `def eu2quat(phi1, Phi, phi2)`

Convert Bunge Euler angles to quaternion representation.
    Euler angles follow the ZXZ convention (Bunge notation).

---

### `euler_angles_reduction`

**Signature**: `def euler_angles_reduction(phi1, Phi, phi2, symops)`

Reduce Euler angles to fundamental zone using symmetry operations.
    Applies crystal symmetry operations to find equivalent orientation
    with Euler angles in the fundamental zone (asymmetric unit).

---

### `grid_s1`

**Signature**: `def grid_s1(resol, grids=6)`

Generate uniformly distributed points on S¹ (circle).
    Used for sampling the third Euler angle in Hopf coordinates.

---

### `hopf2quat`

**Signature**: `def hopf2quat(Points)`

Convert Hopf coordinates to quaternions.
    Hopf coordinates (θ, φ, ψ) parameterize the unit quaternion sphere S³.

---

### `mat_to_quat`

**Signature**: `def mat_to_quat(R)`

Convert a rotation matrix to a quaternion representation.
    Uses the Shepperd's method for numerical stability.

---

### `misori_sym_deg_quats`

**Signature**: `def misori_sym_deg_quats(q1, q2, sym_quats)`

Calculate minimum misorientation angle between two quaternions considering
    crystal symmetry operations. This is a Numba-optimized version.

---

### `misori_sym_deg_sample_to_crystal_fast`

**Signature**: `def misori_sym_deg_sample_to_crystal_fast(M1, M2, symops)`

Compute misorientation angles (deg) between orientations M1 and M2
    considering crystal symmetry operations using vectorized operations.
    This is the fast version optimized with Numba parallel processing.

---

### `mk_pix2xy`

**Signature**: `def mk_pix2xy()`

Create lookup tables for HEALPix pixel index to (x,y) coordinate conversion.
    These tables are used by pix2ang_nest function.

---

### `np_euler_matrix`

**Signature**: `def np_euler_matrix(ai, aj, ak)`

Convert Euler angles to rotation matrix (ZXZ convention, Bunge notation).
    Single orientation version.

---

### `np_eulers_matrices`

**Signature**: `def np_eulers_matrices(data, deg=False)`

Convert multiple Euler angles to rotation matrices (vectorized).
    Processes an array of Euler angle triplets efficiently.

---

### `np_g2quats`

**Signature**: `def np_g2quats(umatsa)`

Convert multiple rotation matrices to quaternions (vectorized).
    Handles arrays of rotation matrices efficiently.

---

### `np_gmat2rodrigues`

**Signature**: `def np_gmat2rodrigues(g)`

Convert rotation matrix to Rodrigues-Frank vector representation.
    Rodrigues vector = rotation_axis * tan(angle/2)

---

### `np_inverse_euler_matrix`

**Signature**: `def np_inverse_euler_matrix(ai, aj, ak)`

Convert Euler angles to inverse (transpose) rotation matrix.
    Equivalent to the transpose of the forward rotation matrix.

---

### `np_ol_R_g`

**Signature**: `def np_ol_R_g(R)`

Convert Rodrigues-Frank vector to rotation matrix (NumPy version).
    NumPy implementation of Rodrigues-Frank to rotation matrix conversion.

---

### `np_ol_R_g2`

**Signature**: `def np_ol_R_g2(R)`

Alternative Rodrigues-Frank to rotation matrix (NumPy version, method 2).
    NumPy implementation of alternative conversion method.

---

### `np_ol_R_q2`

**Signature**: `def np_ol_R_q2(R)`

Convert Rodrigues-Frank vector to quaternion (method 2).
    Transforms Rodrigues-Frank representation to unit quaternion [w,x,y,z].

---

### `np_ol_g_R`

**Signature**: `def np_ol_g_R(g)`

Convert rotation matrix to Rodrigues-Frank vector (NumPy version).
    NumPy implementation returning Rodrigues-Frank vector as numpy array.

---

### `np_ol_g_R2`

**Signature**: `def np_ol_g_R2(g)`

Alternative Rodrigues-Frank conversion (NumPy version, method 2).
    NumPy implementation of alternative Rodrigues-Frank calculation.

---

### `np_ol_g_q2`

**Signature**: `def np_ol_g_q2(g)`

Convert rotation matrix to quaternion (method 2).
    Extracts unit quaternion representation from rotation matrix.

---

### `np_ol_g_rtheta_rad`

**Signature**: `def np_ol_g_rtheta_rad(g)`

Convert rotation matrix to axis-angle representation (NumPy optimized version).

---

### `np_ol_q_g`

**Signature**: `def np_ol_q_g(q)`

Convert quaternion to rotation matrix.
    Transforms unit quaternion [w,x,y,z] to 3×3 rotation matrix.

---

### `np_ol_rtheta_g_rad`

**Signature**: `def np_ol_rtheta_g_rad(r, theta)`

Convert axis-angle representation to rotation matrix (NumPy version).
    Uses Rodrigues' rotation formula.

---

### `np_rodrigues2gmat`

**Signature**: `def np_rodrigues2gmat(rodrigues)`

Convert Rodrigues-Frank vector to rotation matrix.

---

### `nside2npix`

**Signature**: `def nside2npix(nside)`

Calculate the number of pixels in a HEALPix map.
    HEALPix = Hierarchical Equal Area isoLatitude Pixelization.

---

### `ol_R_g`

**Signature**: `def ol_R_g(R)`

Convert Rodrigues-Frank vector to rotation matrix (list version).
    Reconstructs rotation matrix from Rodrigues-Frank representation.

---

### `ol_R_g2`

**Signature**: `def ol_R_g2(R)`

Alternative Rodrigues-Frank to rotation matrix (list version, method 2).
    Second implementation of Rodrigues-Frank to rotation matrix conversion.

---

### `ol_g_R`

**Signature**: `def ol_g_R(g)`

Convert rotation matrix to Rodrigues-Frank vector (list version).
    Calculates Rodrigues-Frank vector R = r·tan(θ/2) from rotation matrix,
    where r is the rotation axis and θ is the rotation angle.

---

### `ol_g_R2`

**Signature**: `def ol_g_R2(g)`

Alternative Rodrigues-Frank conversion (list version, method 2).
    Second implementation of rotation matrix to Rodrigues-Frank vector.
    May use different numerical approach for stability.

---

### `ol_g_rtheta_rad`

**Signature**: `def ol_g_rtheta_rad(g)`

Convert rotation matrix to axis-angle representation (Rodrigues-Frank vector).
    Returns rotation axis and angle.

---

### `ol_rtheta_g_rad`

**Signature**: `def ol_rtheta_g_rad(r, theta)`

Convert axis-angle representation to rotation matrix using Rodrigues' formula.

---

### `passive_rotation`

**Signature**: `def passive_rotation(g, v)`

Perform passive rotation (coordinate transformation).
    Rotates the coordinate system while vector stays fixed in space.
    v' = g^T · v (transpose of rotation matrix).

---

### `pix2ang_nest`

**Signature**: `def pix2ang_nest(nside, ipix, pix2x, pix2y)`

Convert HEALPix pixel index to spherical coordinates (theta, phi).
    Uses NESTED indexing scheme.

---

### `quat_conjugate`

**Signature**: `def quat_conjugate(q)`

Compute the conjugate of a quaternion (inverse for unit quaternions).

---

### `quat_misori_deg`

**Signature**: `def quat_misori_deg(q1, q2)`

Calculate the misorientation angle between two quaternions in degrees.
    Returns the minimum rotation angle needed to go from q1 to q2.

---

### `quat_mult`

**Signature**: `def quat_mult(q1, q2)`

Multiply two quaternions (Hamilton product).
    Order matters: q1 * q2 ≠ q2 * q1 in general.

---

### `quat_multiply`

**Signature**: `def quat_multiply(q1, q2)`

Quaternion multiplication (same as quat_mult, alternative implementation).

---

### `quat_to_mat`

**Signature**: `def quat_to_mat(q)`

Convert a quaternion to a rotation matrix.

---

### `rotation_from_axis_angle`

**Signature**: `def rotation_from_axis_angle(axis, angle)`

Generate rotation matrix from axis-angle representation using Rodrigues formula.
    Creates 3×3 rotation matrix for rotation by 'angle' radians around 'axis'.
    Implements Rodrigues' rotation formula for arbitrary axis rotations.

---

### `simple_grid`

**Signature**: `def simple_grid(resol)`

Generate a uniform grid of quaternions covering orientation space (SO(3)).
    Uses HEALPix for S² sampling and uniform sampling for S¹, combined via Hopf fibration.

---


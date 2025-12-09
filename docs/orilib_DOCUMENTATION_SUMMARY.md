# ORILIB - Documentation Summary

**Module**: `orilib.py`  
**Purpose**: Orientation Analysis and Transformations  
**Total Functions**: 44  
**Last Updated**: December 08, 2025

---

## Function Overview

### QMatproduct

`def QMatproduct(sym, Q):`

Multiply a single symmetry quaternion with multiple quaternions.

### Qlog

`def Qlog(QM):`

Compute the logarithm of quaternions (quaternion logarithm map).

### Qproduct

`def Qproduct(P, Q):`

Compute quaternion product for arrays of quaternions.

### active_rotation

`def active_rotation(an, aboutaxis, deg=False):`

Perform active rotation of vector v by rotation matrix g.

### eu2quat

`def eu2quat(phi1, Phi, phi2):`

Convert Bunge Euler angles to quaternion representation.

### euler_angles_reduction

`def euler_angles_reduction(Phi1,PHI,Phi2):`

Reduce Euler angles to fundamental zone using symmetry operations.

### grid_s1

`def grid_s1(resol, grids=6):`

Generate uniformly distributed points on S¹ (circle).

### hopf2quat

`def hopf2quat(Points):`

Convert Hopf coordinates to quaternions.

### mat_to_quat

`def mat_to_quat(R):`

Convert a rotation matrix to a quaternion representation.

### misori_sym_deg_quats

`def misori_sym_deg_quats(q1, q2, sym_quats):`

Calculate minimum misorientation angle between two quaternions considering

### misori_sym_deg_sample_to_crystal_fast

`def misori_sym_deg_sample_to_crystal_fast(M1, M2, symops):`

Compute misorientation angles (deg) between orientations M1 and M2

### mk_pix2xy

`def mk_pix2xy():`

Create lookup tables for HEALPix pixel index to (x,y) coordinate conversion.

### np_euler_matrix

`def np_euler_matrix(ai, aj, ak):`

Convert Euler angles to rotation matrix (ZXZ convention, Bunge notation).

### np_eulers_matrices

`def np_eulers_matrices(data, deg=False):`

Convert multiple Euler angles to rotation matrices (vectorized).

### np_g2quats

`def np_g2quats(umatsa):`

Convert multiple rotation matrices to quaternions (vectorized).

### np_gmat2rodrigues

`def np_gmat2rodrigues(g):`

Convert rotation matrix to Rodrigues-Frank vector representation.

### np_inverse_euler_matrix

`def np_inverse_euler_matrix(ai, aj, ak):`

Convert Euler angles to inverse (transpose) rotation matrix.

### np_ol_R_g

`def np_ol_R_g (R):`

Convert Rodrigues-Frank vector to rotation matrix (NumPy version).

### np_ol_R_g2

`def np_ol_R_g2(R,epsilon, delta):`

Alternative Rodrigues-Frank to rotation matrix (NumPy version, method 2).

### np_ol_R_q2

`def np_ol_R_q2(R):`

Convert Rodrigues-Frank vector to quaternion (method 2).

### np_ol_g_R

`def np_ol_g_R(g):`

Convert rotation matrix to Rodrigues-Frank vector (NumPy version).

### np_ol_g_R2

`def np_ol_g_R2(g,epsilon, delta):`

Alternative Rodrigues-Frank conversion (NumPy version, method 2).

### np_ol_g_q2

`def np_ol_g_q2(g):`

Convert rotation matrix to quaternion (method 2).

### np_ol_g_rtheta_rad

`def np_ol_g_rtheta_rad(g):`

Convert rotation matrix to axis-angle representation (NumPy optimized version).

### np_ol_q_g

`def np_ol_q_g(q):`

Convert quaternion to rotation matrix.

### np_ol_rtheta_g_rad

`def np_ol_rtheta_g_rad(r, theta):`

Convert axis-angle representation to rotation matrix (NumPy version).

### np_rodrigues2gmat

`def np_rodrigues2gmat(rodrigues):`

Convert Rodrigues-Frank vector to rotation matrix.

### nside2npix

`def nside2npix(nside):`

Calculate the number of pixels in a HEALPix map.

### ol_R_g

`def ol_R_g (R):`

Convert Rodrigues-Frank vector to rotation matrix (list version).

### ol_R_g2

`def ol_R_g2(R):`

Alternative Rodrigues-Frank to rotation matrix (list version, method 2).

### ol_g_R

`def ol_g_R(g):`

Convert rotation matrix to Rodrigues-Frank vector (list version).

### ol_g_R2

`def ol_g_R2(g):`

Alternative Rodrigues-Frank conversion (list version, method 2).

### ol_g_rtheta_rad

`def ol_g_rtheta_rad(g):`

Convert rotation matrix to axis-angle representation (Rodrigues-Frank vector).

### ol_rtheta_g_rad

`def ol_rtheta_g_rad(r, theta):`

Convert axis-angle representation to rotation matrix using Rodrigues' formula.

### passive_rotation

`def passive_rotation(an, aboutaxis, deg=False):`

Perform passive rotation (coordinate transformation).

### pix2ang_nest

`def pix2ang_nest(nside, ipix, pix2x, pix2y):`

Convert HEALPix pixel index to spherical coordinates (theta, phi).

### quat_conjugate

`def quat_conjugate(q):`

Compute the conjugate of a quaternion (inverse for unit quaternions).

### quat_misori_deg

`def quat_misori_deg(q1, q2):`

Calculate the misorientation angle between two quaternions in degrees.

### quat_mult

`def quat_mult(q1, q2):`

Multiply two quaternions (Hamilton product).

### quat_multiply

`def quat_multiply(q1, q2):`

Quaternion multiplication (same as quat_mult, alternative implementation).

### quat_to_mat

`def quat_to_mat(q):`

Convert a quaternion to a rotation matrix.

### rotation_from_axis_angle

`def rotation_from_axis_angle(ax,an,deg=False):`

Generate rotation matrix from axis-angle representation using Rodrigues formula.

### simple_grid

`def simple_grid(resol):`

Generate a uniform grid of quaternions covering orientation space (SO(3)).

### symmetry_elements

`def symmetry_elements(lattice):`

Generate symmetry operation matrices for crystal system.


---

**Total**: 44 functions

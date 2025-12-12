# ORILIB - Complete Documentation

**Module**: `orilib.py`  
**Purpose**: Orientation Analysis and Transformations  
**Functions**: 59  
**Last Updated**: December 12, 2025

---

## Table of Contents

1. [Mat2Quat](#function-mat2quat)
2. [Mat2Quat_ini](#function-mat2quat_ini)
3. [QMatproduct](#function-qmatproduct)
4. [Qlog](#function-qlog)
5. [Qproduct](#function-qproduct)
6. [active_rotation](#function-active_rotation)
7. [disorimat](#function-disorimat)
8. [disorimat_ini](#function-disorimat_ini)
9. [disorimat_test01](#function-disorimat_test01)
10. [disorimat_test02](#function-disorimat_test02)
11. [eu2quat](#function-eu2quat)
12. [euler_angles_from_matrix](#function-euler_angles_from_matrix)
13. [euler_angles_reduction](#function-euler_angles_reduction)
14. [grid_s1](#function-grid_s1)
15. [hopf2quat](#function-hopf2quat)
16. [mat2quat02](#function-mat2quat02)
17. [mat_to_quat](#function-mat_to_quat)
18. [misori_sym_deg_quats](#function-misori_sym_deg_quats)
19. [misori_sym_deg_sample_to_crystal_fast](#function-misori_sym_deg_sample_to_crystal_fast)
20. [misorimat](#function-misorimat)
21. [misorimat_ini](#function-misorimat_ini)
22. [mk_pix2xy](#function-mk_pix2xy)
23. [np_euler_matrix](#function-np_euler_matrix)
24. [np_eulers_matrices](#function-np_eulers_matrices)
25. [np_g2quats](#function-np_g2quats)
26. [np_gmat2rodrigues](#function-np_gmat2rodrigues)
27. [np_inverse_euler_matrix](#function-np_inverse_euler_matrix)
28. [np_ol_R_g](#function-np_ol_r_g)
29. [np_ol_R_g2](#function-np_ol_r_g2)
30. [np_ol_R_q2](#function-np_ol_r_q2)
31. [np_ol_g_R](#function-np_ol_g_r)
32. [np_ol_g_R2](#function-np_ol_g_r2)
33. [np_ol_g_q2](#function-np_ol_g_q2)
34. [np_ol_g_rtheta_rad](#function-np_ol_g_rtheta_rad)
35. [np_ol_q_g](#function-np_ol_q_g)
36. [np_ol_rtheta_g_rad](#function-np_ol_rtheta_g_rad)
37. [np_rodrigues2gmat](#function-np_rodrigues2gmat)
38. [nside2npix](#function-nside2npix)
39. [ol_R_g](#function-ol_r_g)
40. [ol_R_g2](#function-ol_r_g2)
41. [ol_g_R](#function-ol_g_r)
42. [ol_g_R2](#function-ol_g_r2)
43. [ol_g_rtheta_rad](#function-ol_g_rtheta_rad)
44. [ol_rtheta_g_rad](#function-ol_rtheta_g_rad)
45. [orilistMult](#function-orilistmult)
46. [passive_rotation](#function-passive_rotation)
47. [pix2ang_nest](#function-pix2ang_nest)
48. [quat_conjugate](#function-quat_conjugate)
49. [quat_misori_deg](#function-quat_misori_deg)
50. [quat_mult](#function-quat_mult)
51. [quat_multiply](#function-quat_multiply)
52. [quat_to_mat](#function-quat_to_mat)
53. [rotation_from_axis_angle](#function-rotation_from_axis_angle)
54. [simple_grid](#function-simple_grid)
55. [symmetry_elements](#function-symmetry_elements)
56. [symmetry_reduced_oris](#function-symmetry_reduced_oris)
57. [symposMult](#function-symposmult)
58. [symposMult02](#function-symposmult02)
59. [trace_to_angle](#function-trace_to_angle)

---

## Function: Mat2Quat

**Signature**:
```python
def Mat2Quat(umatsa): #orientation matrix to quaternion
```

**Description**:

Convert rotation matrix to quaternion representation.
    General matrix to quaternion conversion. Alternative implementation
    that may use different numerical approach than mat_to_quat.

**Input**:

R: numpy array (3, 3) - Rotation matrix (orthogonal with det=1)

**Output**:

q: numpy array (4,) - Unit quaternion [w, x, y, z]

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Random rotation matrix
        >>> from scipy.spatial.transform import Rotation as Rot
        >>> R = Rot.random().as_matrix()
        >>> 
        >>> # Convert to quaternion
        >>> q = Mat2Quat(R)
        >>> print("Quaternion:", q)
        >>> print("Magnitude:", np.linalg.norm(q))  # Should be 1
        >>> 
        >>> # Round trip test
        >>> R_reconstructed = quat_to_mat(q)
        >>> print("Reconstruction match:", np.allclose(R, R_reconstructed))
        >>> # Use for orientation averaging
        >>> matrices = Rot.random(10).as_matrix()
        >>> quats = np.array([Mat2Quat(R) for R in matrices])
        >>> print("Quaternions shape:", quats.shape)  # (10, 4)
```

---

## Function: Mat2Quat_ini

**Signature**:
```python
def Mat2Quat_ini(umatsa): #orientation matrix to quaternion
```

**Description**:

Initialize matrix to quaternion conversion.
    Preliminary/initialization version of rotation matrix to quaternion
    conversion. May be an earlier implementation or setup function.

**Input**:

R: numpy array (3, 3) - Rotation matrix

**Output**:

q: numpy array (4,) - Quaternion [w, x, y, z]

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Simple 90° rotation around Z
        >>> R = np.array([[0, -1, 0],
        ...               [1,  0, 0],
        ...               [0,  0, 1]], dtype=float)
        >>> 
        >>> q = Mat2Quat_ini(R)
        >>> print("Quaternion:", q)
        >>> # Should give [0.707, 0, 0, 0.707] approximately
        >>> # Convert back to verify
        >>> R_back = quat_to_mat(q)
        >>> print("Reconstruction error:", np.max(np.abs(R - R_back)))
```

---

## Function: QMatproduct

**Signature**:
```python
def QMatproduct(sym, Q):
```

**Description**:

Multiply a single symmetry quaternion with multiple quaternions.
    Applies the same symmetry operation to an array of orientations.

**Input**:

sym: numpy array (4,) - Single symmetry quaternion [w, x, y, z]
        Q: numpy array (4, N) - Array of quaternions to transform

**Output**:

SQ: numpy array (4, N) - Transformed quaternions (sym * Q)

**Usage Example**:

```python
>>> import numpy as np
        >>> # Define a 90° rotation symmetry around Z
        >>> sym = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> 
        >>> # Multiple orientations
        >>> N = 50
        >>> Q = np.random.randn(4, N)
        >>> Q = Q / np.linalg.norm(Q, axis=0)
        >>> 
        >>> # Apply symmetry to all
        >>> SQ = QMatproduct(sym, Q)
        >>> print("Transformed shape:", SQ.shape)  # (4, 50)
```

---

## Function: Qlog

**Signature**:
```python
def Qlog(QM):
```

**Description**:

Compute the logarithm of quaternions (quaternion logarithm map).
    Maps quaternions from unit sphere S³ to tangent space at identity.

**Input**:

QM: numpy array (4, N, M) - Array of quaternions [w, x, y, z, ...]

**Output**:

qlog: numpy array (4, N, M) - Logarithm of quaternions
                                       qlog[0] = 0, qlog[1:4] = v * angle

**Usage Example**:

```python
>>> import numpy as np
        >>> # Single quaternion for 90° rotation
        >>> q = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> # Reshape for function (4, 1, 1)
        >>> QM = q.reshape(4, 1, 1)
        >>> qlog = Qlog(QM)
        >>> print("Quaternion log:")
        >>> print(qlog[:, 0, 0])
        >>> # qlog[0] should be 0, qlog[3] should be π/2
```

---

## Function: Qproduct

**Signature**:
```python
def Qproduct(P, Q):
```

**Description**:

Compute quaternion product for arrays of quaternions.
    Computes P * Q for each pair of quaternions.

**Input**:

P: numpy array (4, N) - First array of quaternions [w, x, y, z]
        Q: numpy array (4, N) - Second array of quaternions [w, x, y, z]

**Output**:

result: numpy array (4, N) - Product quaternions P * Q

**Usage Example**:

```python
>>> import numpy as np
        >>> # Multiple quaternion pairs
        >>> N = 100
        >>> P = np.random.randn(4, N)
        >>> P = P / np.linalg.norm(P, axis=0)  # Normalize
        >>> Q = np.random.randn(4, N)
        >>> Q = Q / np.linalg.norm(Q, axis=0)
        >>> 
        >>> result = Qproduct(P, Q)
        >>> print("Result shape:", result.shape)  # (4, 100)
        >>> # Verify results are normalized
        >>> norms = np.linalg.norm(result, axis=0)
        >>> print("All normalized:", np.allclose(norms, 1.0))
```

---

## Function: active_rotation

**Signature**:
```python
def active_rotation(an, aboutaxis, deg=False):
```

**Description**:

Perform active rotation of vector v by rotation matrix g.
    Rotates the vector itself while keeping coordinate system fixed.
    v' = g · v (matrix-vector multiplication).

**Input**:

g (array 3×3): Rotation matrix
        v (array [3]): Vector to rotate

**Output**:

numpy.ndarray [3]: Rotated vector v'

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # 90° rotation around Z (active)
        >>> g_z90 = np.array([[0, -1, 0],
        ...                   [1,  0, 0],
        ...                   [0,  0, 1]])
        >>> v = np.array([1, 0, 0])
        >>> v_rot = active_rotation(g_z90, v)
        >>> print(v_rot)  # [0, 1, 0]
        >>> 
        >>> # Verify: [1,0,0] rotates to [0,1,0]
        >>> assert np.allclose(v_rot, [0, 1, 0])
    Notes:
        - Active = rotate the object
        - Compare with passive_rotation (rotate coordinates)
        - Common in physics and mechanics
        - v' = g · v
    Formula:
        v'_i = g_{ij} v_j
```

---

## Function: disorimat

**Signature**:
```python
def disorimat(umatsa,symops,prnt=False,withfirst=False,eqmats=False):
```

**Description**:

Calculate disorientation matrix considering crystal symmetry.
    Computes the minimum misorientation between two orientations by
    considering all symmetrically equivalent variants. This is the
    crystallographically meaningful misorientation.

**Input**:

M1: numpy array (3, 3) - First orientation matrix
        M2: numpy array (3, 3) - Second orientation matrix
        symops: numpy array (Ns, 3, 3) - Crystal symmetry operations

**Output**:

disori: numpy array (3, 3) - Disorientation matrix (minimum misorientation)
        angle: float - Disorientation angle in degrees

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # EBSD grain orientations
        >>> grain1_euler = np.radians([120, 45, 80])
        >>> grain2_euler = np.radians([130, 50, 85])
        >>> 
        >>> M1 = np_euler_matrix(*grain1_euler)
        >>> M2 = np_euler_matrix(*grain2_euler)
        >>> 
        >>> # Cubic symmetry (get from symmetry_elements function)
        >>> symops = symmetry_elements('cubic')
        >>> 
        >>> # Calculate disorientation
        >>> disori, angle = disorimat(M1, M2, symops)
        >>> print(f"Grain boundary misorientation: {angle:.2f}°")
        >>> 
        >>> # Classify grain boundary type
        >>> if angle < 15:
        ...     print("Low-angle grain boundary")
        >>> elif angle > 15:
        ...     print("High-angle grain boundary")
        >>> # Get disorientation axis
        >>> axis, _ = np_ol_g_rtheta_rad(disori)
        >>> print(f"Rotation axis: [{axis[0]:.3f}, {axis[1]:.3f}, {axis[2]:.3f}]")
```

---

## Function: disorimat_ini

**Signature**:
```python
def disorimat_ini(umatsa,symops):
```

---

## Function: disorimat_test01

**Signature**:
```python
def disorimat_test01(umatsa,symops):
```

**Description**:

Test version 1 for disorientation matrix calculation.
    First test implementation of disorientation computation. Used for
    validating algorithms before final implementation.

**Input**:

M1: numpy array (3, 3) - First orientation matrix
        M2: numpy array (3, 3) - Second orientation matrix
        symops: numpy array (Ns, 3, 3) - Symmetry operations

**Output**:

disori: numpy array (3, 3) - Disorientation matrix
        angle: float - Disorientation angle in degrees

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Create test orientations
        >>> M1 = np_euler_matrix(0, 0, 0)
        >>> M2 = np_euler_matrix(np.pi/6, 0, 0)  # 30° rotation
        >>> 
        >>> # Identity symmetry
        >>> symops = np.array([np.eye(3)])
        >>> 
        >>> disori, angle = disorimat_test01(M1, M2, symops)
        >>> print(f"Test result angle: {angle:.2f}°")
        >>> # Compare with test02
        >>> disori2, angle2 = disorimat_test02(M1, M2, symops)
        >>> print(f"Difference: {abs(angle - angle2):.4f}°")
```

---

## Function: disorimat_test02

**Signature**:
```python
def disorimat_test02(umatsa,symops):
```

**Description**:

Test version 2 for disorientation matrix calculation.
    Experimental/testing version of disorientation calculation. The
    disorientation is the minimum misorientation considering crystal
    symmetry operations.

**Input**:

M1: numpy array (3, 3) - First orientation matrix
        M2: numpy array (3, 3) - Second orientation matrix
        symops: numpy array (Ns, 3, 3) - Symmetry operations

**Output**:

disori: numpy array (3, 3) - Disorientation matrix
        angle: float - Disorientation angle in degrees

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Test with known misorientations
        >>> M1 = np.eye(3)
        >>> M2 = np_euler_matrix(np.pi/4, 0, 0)  # 45° rotation
        >>> 
        >>> # Simple symmetry (identity only)
        >>> symops = np.array([np.eye(3)])
        >>> 
        >>> disori, angle = disorimat_test02(M1, M2, symops)
        >>> print(f"Disorientation angle: {angle:.2f}°")
        >>> print("Disorientation matrix:")
        >>> print(disori)
```

---

## Function: eu2quat

**Signature**:
```python
def eu2quat(phi1, Phi, phi2):
```

**Description**:

Convert Bunge Euler angles to quaternion representation.
    Euler angles follow the ZXZ convention (Bunge notation).

**Input**:

phi1: float - First Euler angle (rotation around Z) in radians
        Phi: float - Second Euler angle (rotation around X') in radians
        phi2: float - Third Euler angle (rotation around Z'') in radians

**Output**:

q: numpy array (4,) - Quaternion [w, x, y, z] with positive w

**Usage Example**:

```python
>>> import numpy as np
        >>> # 45-degree rotation around each axis
        >>> phi1 = np.radians(45)
        >>> Phi = np.radians(45)
        >>> phi2 = np.radians(45)
        >>> q = eu2quat(phi1, Phi, phi2)
        >>> print("Quaternion:", q)
        >>> # Identity rotation (all zeros)
        >>> q_identity = eu2quat(0, 0, 0)
        >>> print("Identity quaternion:", q_identity)
        >>> # Should give [1, 0, 0, 0]
```

---

## Function: euler_angles_from_matrix

**Signature**:
```python
def euler_angles_from_matrix(Rl,deg=False):
```

**Description**:

Extract Euler angles from a rotation matrix.
    Inverse operation of np_euler_matrix. Converts a rotation matrix back
    to Bunge Euler angles (ZXZ convention).

**Input**:

R: numpy array (3, 3) - Rotation matrix (orthogonal with det=1)

**Output**:

euler: numpy array (3,) - Euler angles [phi1, Phi, phi2] in radians

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Create rotation matrix from Euler angles
        >>> phi1, Phi, phi2 = np.radians([45, 60, 30])
        >>> R = np_euler_matrix(phi1, Phi, phi2)
        >>> 
        >>> # Extract Euler angles back
        >>> euler_recovered = euler_angles_from_matrix(R)
        >>> print("Original:", [phi1, Phi, phi2])
        >>> print("Recovered:", euler_recovered)
        >>> # Should match within numerical precision
        >>> # Handle gimbal lock cases
        >>> # When Phi = 0 or 180 degrees
        >>> R_gimbal = np_euler_matrix(0, 0, 0)
        >>> euler_gimbal = euler_angles_from_matrix(R_gimbal)
        >>> print("Gimbal case:", euler_gimbal)
        >>> # Convert to degrees for readability
        >>> euler_deg = np.degrees(euler_recovered)
        >>> print("Euler angles (degrees):", euler_deg)
```

---

## Function: euler_angles_reduction

**Signature**:
```python
def euler_angles_reduction(Phi1,PHI,Phi2):
```

**Description**:

Reduce Euler angles to fundamental zone using symmetry operations.
    Applies crystal symmetry operations to find equivalent orientation
    with Euler angles in the fundamental zone (asymmetric unit).

**Input**:

phi1, Phi, phi2 (float): Euler angles in radians (Bunge convention)
        symops (list): List of 3×3 symmetry operation matrices

**Output**:

tuple: (phi1_red, Phi_red, phi2_red) - Reduced Euler angles in radians

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Cubic symmetry (simplified - just identity for demo)
        >>> symops = [np.eye(3)]
        >>> 
        >>> # Reduce orientation
        >>> phi1, Phi, phi2 = np.radians([45, 30, 60])
        >>> phi1_r, Phi_r, phi2_r = euler_angles_reduction(phi1, Phi, phi2, symops)
        >>> print(f"Reduced: φ1={np.degrees(phi1_r):.1f}°")
    Notes:
        - Fundamental zone depends on crystal symmetry
        - Cubic: 0≤φ1≤90°, 0≤Φ≤45°, 0≤φ2≤90°
        - Reduces orientation distribution function (ODF) storage
        - Essential for texture analysis
```

---

## Function: grid_s1

**Signature**:
```python
def grid_s1(resol, grids=6):
```

**Description**:

Generate uniformly distributed points on S¹ (circle).
    Used for sampling the third Euler angle in Hopf coordinates.

**Input**:

resol: int - Resolution parameter (number of points = 2^resol * grids)
        grids: int - Grid multiplier (default: 6)

**Output**:

points: list of floats - Angles in radians from 0 to 2π

**Usage Example**:

```python
>>> points = grid_s1(resol=2, grids=6)
        >>> print("Number of points:", len(points))
        >>> # 2^2 * 6 = 24 points
        >>> print("First few points:", points[:5])
        >>> print("Last point:", points[-1])
        >>> # Higher resolution
        >>> points_fine = grid_s1(resol=4, grids=6)
        >>> print("Fine grid points:", len(points_fine))
        >>> # 2^4 * 6 = 96 points
```

---

## Function: hopf2quat

**Signature**:
```python
def hopf2quat(Points):
```

**Description**:

Convert Hopf coordinates to quaternions.
    Hopf coordinates (θ, φ, ψ) parameterize the unit quaternion sphere S³.

**Input**:

Points: list of tuples - Each tuple contains (theta, phi, psi) in radians
                                 theta ∈ [0, π], phi ∈ [0, 2π], psi ∈ [0, 2π]

**Output**:

quats: list of lists - Quaternions [w, x, y, z] for each point

**Usage Example**:

```python
>>> import numpy as np
        >>> # Single point in Hopf coordinates
        >>> points = [(np.pi/2, 0, 0)]
        >>> quats = hopf2quat(points)
        >>> print("Quaternion:", quats[0])
        >>> # Multiple points
        >>> points_multiple = [
        ...     (0, 0, 0),           # Identity
        ...     (np.pi/2, 0, 0),     # 
        ...     (np.pi/2, np.pi, 0)  # 
        ... ]
        >>> quats_multiple = hopf2quat(points_multiple)
        >>> print("Number of quaternions:", len(quats_multiple))
```

---

## Function: mat2quat02

**Signature**:
```python
def mat2quat02(matrix): #orientation matrix to quaternion
```

**Description**:

Alternative matrix to quaternion conversion (version 2).
    Different algorithm for converting rotation matrix to quaternion.
    May handle numerical edge cases differently than mat_to_quat.

**Input**:

R: numpy array (3, 3) - Rotation matrix

**Output**:

q: numpy array (4,) - Quaternion [w, x, y, z]

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Create rotation matrix
        >>> angle = np.pi/3  # 60 degrees
        >>> axis = np.array([1, 1, 1]) / np.sqrt(3)
        >>> R = np_ol_rtheta_g_rad(axis, angle)
        >>> 
        >>> # Convert to quaternion
        >>> q = mat2quat02(R)
        >>> print("Quaternion:", q)
        >>> 
        >>> # Compare with standard method
        >>> q_standard = mat_to_quat(R)
        >>> print("Standard:", q_standard)
        >>> 
        >>> # Both should give same result (within sign)
        >>> print("Match:", np.allclose(np.abs(q), np.abs(q_standard)))
```

---

## Function: mat_to_quat

**Signature**:
```python
def mat_to_quat(R):
```

**Description**:

Convert a rotation matrix to a quaternion representation.
    Uses the Shepperd's method for numerical stability.

**Input**:

R: numpy array (3x3) - Rotation matrix (orthogonal with det=1)

**Output**:

q: numpy array (4,) - Unit quaternion [w, x, y, z] where w is the scalar part

**Usage Example**:

```python
>>> import numpy as np
        >>> # 90-degree rotation around z-axis
        >>> R = np.array([[0, -1, 0],
        ...               [1,  0, 0],
        ...               [0,  0, 1]], dtype=float)
        >>> q = mat_to_quat(R)
        >>> print("Quaternion:", q)
        >>> # Should give approximately [0.707, 0, 0, 0.707]
        >>> # Identity rotation
        >>> R_identity = np.eye(3)
        >>> q_identity = mat_to_quat(R_identity)
        >>> print("Identity quaternion:", q_identity)
        >>> # Should give [1, 0, 0, 0]
```

---

## Function: misori_sym_deg_quats

**Signature**:
```python
def misori_sym_deg_quats(q1, q2, sym_quats):
```

**Description**:

Calculate minimum misorientation angle between two quaternions considering
    crystal symmetry operations. This is a Numba-optimized version.

**Input**:

q1: numpy array (4,) - First quaternion [w, x, y, z]
        q2: numpy array (4,) - Second quaternion [w, x, y, z]
        sym_quats: numpy array (Ns, 4) - Array of symmetry operation quaternions

**Output**:

min_angle: float - Minimum misorientation angle in degrees

**Usage Example**:

```python
>>> import numpy as np
        >>> # Define two orientations
        >>> q1 = np.array([1.0, 0.0, 0.0, 0.0])
        >>> q2 = np.array([0.707, 0.0, 0.0, 0.707])  # 90° around z
        >>> 
        >>> # Define cubic symmetry (simplified - just identity for demo)
        >>> sym_quats = np.array([[1.0, 0.0, 0.0, 0.0]])
        >>> 
        >>> angle = misori_sym_deg_quats(q1, q2, sym_quats)
        >>> print("Misorientation with symmetry:", angle, "degrees")
```

---

## Function: misori_sym_deg_sample_to_crystal_fast

**Signature**:
```python
def misori_sym_deg_sample_to_crystal_fast(M1, M2, symops):
```

**Description**:

Compute misorientation angles (deg) between orientations M1 and M2
    considering crystal symmetry operations using vectorized operations.
    This is the fast version optimized with Numba parallel processing.

**Input**:

M1: numpy array (N, 3, 3) - Orientation matrices (sample→crystal reference frame)
        M2: numpy array (N, 3, 3) - Orientation matrices (sample→crystal reference frame)
                                     Must have same length as M1
        symops: numpy array (Ns, 3, 3) - Crystal symmetry operation matrices

**Output**:

miso: numpy array (N,) - Minimum misorientation angle in degrees for each pair

**Usage Example**:

```python
>>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R
        >>> 
        >>> # Generate some random orientations
        >>> N = 100
        >>> M1 = R.random(N).as_matrix()
        >>> M2 = R.random(N).as_matrix()
        >>> 
        >>> # Cubic symmetry operations (just identity for example)
        >>> symops = np.array([np.eye(3)])
        >>> 
        >>> # Calculate misorientations
        >>> misorientations = misori_sym_deg_sample_to_crystal_fast(M1, M2, symops)
        >>> print("Mean misorientation:", np.mean(misorientations), "degrees")
        >>> print("Max misorientation:", np.max(misorientations), "degrees")
```

---

## Function: misorimat

**Signature**:
```python
def misorimat(umatsa):
```

**Description**:

Calculate misorientation matrix between two orientations.
    Computes the relative rotation (misorientation) between two crystal
    orientations. The misorientation matrix represents the rotation needed
    to go from orientation M1 to orientation M2.

**Input**:

M1: numpy array (3, 3) - First orientation matrix (sample→crystal)
        M2: numpy array (3, 3) - Second orientation matrix (sample→crystal)
        symops: numpy array (Ns, 3, 3) - Crystal symmetry operations (optional)

**Output**:

misori: numpy array (3, 3) - Misorientation matrix
        angle: float - Misorientation angle in degrees (if symmetry applied)

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Two grain orientations
        >>> euler1 = np.radians([30, 45, 60])
        >>> euler2 = np.radians([60, 50, 70])
        >>> M1 = np_euler_matrix(*euler1)
        >>> M2 = np_euler_matrix(*euler2)
        >>> 
        >>> # Calculate misorientation (no symmetry)
        >>> misori = misorimat(M1, M2)
        >>> 
        >>> # Get misorientation angle
        >>> trace = np.trace(misori)
        >>> angle = np.degrees(np.arccos((trace - 1) / 2))
        >>> print(f"Misorientation angle: {angle:.2f}°")
        >>> # With cubic symmetry
        >>> symops = symmetry_elements('cubic')
        >>> misori_sym, angle_sym = misorimat(M1, M2, symops)
        >>> print(f"Minimum misorientation angle: {angle_sym:.2f}°")
```

---

## Function: misorimat_ini

**Signature**:
```python
def misorimat_ini(umatsa):
```

**Description**:

Initialize misorientation matrix calculation (initialization version).
    Preliminary version of misorientation matrix computation. Sets up
    the calculation framework for determining misorientation between
    two orientations.

**Input**:

M1: numpy array (3, 3) - First orientation matrix
        M2: numpy array (3, 3) - Second orientation matrix

**Output**:

misori: numpy array (3, 3) - Misorientation matrix M2 * M1^T

**Usage Example**:

```python
>>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R
        >>> 
        >>> # Two random orientations
        >>> M1 = R.random().as_matrix()
        >>> M2 = R.random().as_matrix()
        >>> 
        >>> # Calculate misorientation
        >>> misori = misorimat_ini(M1, M2)
        >>> print("Misorientation matrix:")
        >>> print(misori)
        >>> 
        >>> # Verify it's a rotation matrix
        >>> print("Det:", np.linalg.det(misori))  # Should be 1
        >>> print("Orthogonal:", np.allclose(misori @ misori.T, np.eye(3)))
```

---

## Function: mk_pix2xy

**Signature**:
```python
def mk_pix2xy():
```

**Description**:

Create lookup tables for HEALPix pixel index to (x,y) coordinate conversion.
    These tables are used by pix2ang_nest function.

**Input**:

None

**Output**:

pix2x: list (1024,) - X-coordinate lookup table
        pix2y: list (1024,) - Y-coordinate lookup table

**Usage Example**:

```python
>>> pix2x, pix2y = mk_pix2xy()
        >>> print("Table size:", len(pix2x))  # 1024
        >>> print("First few X values:", pix2x[:5])
        >>> print("First few Y values:", pix2y[:5])
        >>> 
        >>> # Use with pix2ang_nest
        >>> nside = 8
        >>> theta, phi = pix2ang_nest(nside, ipix=100, pix2x=pix2x, pix2y=pix2y)
```

---

## Function: np_euler_matrix

**Signature**:
```python
def np_euler_matrix(ai, aj, ak):
```

**Description**:

Convert Euler angles to rotation matrix (ZXZ convention, Bunge notation).
    Single orientation version.

**Input**:

ai: float - First Euler angle (phi1) in radians
        aj: float - Second Euler angle (Phi) in radians
        ak: float - Third Euler angle (phi2) in radians

**Output**:

g: numpy array (3, 3) - Rotation matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> # 90-degree rotation around Z axis
        >>> phi1 = np.pi/2
        >>> Phi = 0
        >>> phi2 = 0
        >>> g = np_euler_matrix(phi1, Phi, phi2)
        >>> print("Rotation matrix:")
        >>> print(g)
        >>> # Identity rotation
        >>> g_identity = np_euler_matrix(0, 0, 0)
        >>> print("Identity:", g_identity)
```

---

## Function: np_eulers_matrices

**Signature**:
```python
def np_eulers_matrices(data, deg=False):
```

**Description**:

Convert multiple Euler angles to rotation matrices (vectorized).
    Processes an array of Euler angle triplets efficiently.

**Input**:

data: numpy array (N, 3) - Array of Euler angles [phi1, Phi, phi2]
        deg: bool - If True, input angles are in degrees; if False, radians (default: False)

**Output**:

g: numpy array (N, 3, 3) - Array of rotation matrices

**Usage Example**:

```python
>>> import numpy as np
        >>> # Multiple orientations in degrees
        >>> euler_angles = np.array([
        ...     [0, 0, 0],        # Identity
        ...     [90, 0, 0],       # 90° around Z
        ...     [0, 90, 0],       # 90° around X'
        ...     [45, 45, 45]      # Combined rotation
        ... ])
        >>> matrices = np_eulers_matrices(euler_angles, deg=True)
        >>> print("Shape:", matrices.shape)  # (4, 3, 3)
        >>> print("First matrix (identity):")
        >>> print(matrices[0])
        >>> # Process many orientations efficiently
        >>> N = 1000
        >>> random_eulers = np.random.rand(N, 3) * np.array([360, 180, 360])
        >>> rot_matrices = np_eulers_matrices(random_eulers, deg=True)
        >>> print(f"Processed {N} orientations")
```

---

## Function: np_g2quats

**Signature**:
```python
def np_g2quats(umatsa):
```

**Description**:

Convert multiple rotation matrices to quaternions (vectorized).
    Handles arrays of rotation matrices efficiently.

**Input**:

umatsa: numpy array (N, 3, 3) - Array of rotation matrices

**Output**:

Q: numpy array (4, N) - Array of quaternions [w, x, y, z] for each matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> # Create multiple rotation matrices
        >>> N = 100
        >>> from scipy.spatial.transform import Rotation as R
        >>> matrices = R.random(N).as_matrix()
        >>> 
        >>> # Convert all to quaternions
        >>> quats = np_g2quats(matrices)
        >>> print("Shape:", quats.shape)  # (4, 100)
        >>> print("First quaternion:", quats[:, 0])
        >>> 
        >>> # Verify quaternions are normalized
        >>> norms = np.linalg.norm(quats, axis=0)
        >>> print("All normalized:", np.allclose(norms, 1.0))
```

---

## Function: np_gmat2rodrigues

**Signature**:
```python
def np_gmat2rodrigues(g):
```

**Description**:

Convert rotation matrix to Rodrigues-Frank vector representation.
    Rodrigues vector = rotation_axis * tan(angle/2)

**Input**:

g: numpy array (3, 3) - Rotation matrix

**Output**:

rodrigues: numpy array (3,) - Rodrigues-Frank vector

**Usage Example**:

```python
>>> import numpy as np
        >>> # 90-degree rotation around Z
        >>> g = np.array([[0, -1, 0],
        ...               [1,  0, 0],
        ...               [0,  0, 1]], dtype=float)
        >>> rod = np_gmat2rodrigues(g)
        >>> print("Rodrigues vector:", rod)
        >>> # For 90° rotation: tan(45°) = 1, so rod ≈ [0, 0, 1]
        >>> # Small rotation (5 degrees around X)
        >>> angle = np.radians(5)
        >>> g_small = np.array([[1, 0, 0],
        ...                     [0, np.cos(angle), -np.sin(angle)],
        ...                     [0, np.sin(angle), np.cos(angle)]])
        >>> rod_small = np_gmat2rodrigues(g_small)
        >>> print("Small rotation Rodrigues:", rod_small)
```

---

## Function: np_inverse_euler_matrix

**Signature**:
```python
def np_inverse_euler_matrix(ai, aj, ak):
```

**Description**:

Convert Euler angles to inverse (transpose) rotation matrix.
    Equivalent to the transpose of the forward rotation matrix.

**Input**:

ai: float - First Euler angle (phi1) in radians
        aj: float - Second Euler angle (Phi) in radians
        ak: float - Third Euler angle (phi2) in radians

**Output**:

U: numpy array (3, 3) - Inverse rotation matrix (transpose of forward matrix)

**Usage Example**:

```python
>>> import numpy as np
        >>> # Forward rotation
        >>> ai, aj, ak = np.pi/4, np.pi/3, np.pi/6
        >>> g_forward = np_euler_matrix(ai, aj, ak)
        >>> g_inverse = np_inverse_euler_matrix(ai, aj, ak)
        >>> 
        >>> # Verify: forward * inverse = identity
        >>> result = g_forward @ g_inverse
        >>> print("Forward @ Inverse:")
        >>> print(result)
        >>> # Should be close to identity matrix
```

---

## Function: np_ol_R_g

**Signature**:
```python
def np_ol_R_g (R):
```

**Description**:

Convert Rodrigues-Frank vector to rotation matrix (NumPy version).
    NumPy implementation of Rodrigues-Frank to rotation matrix conversion.

**Input**:

R (numpy.ndarray [3]): Rodrigues-Frank vector

**Output**:

numpy.ndarray (3×3): Rotation matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Small rotation
        >>> R_small = np.array([0.01, 0.02, 0.03])
        >>> g = np_ol_R_g(R_small)
        >>> print(f"Rotation angle: {np.degrees(2*np.arctan(np.linalg.norm(R_small))):.2f}°")
        >>> 
        >>> # Round trip test
        >>> g_orig = rotation_from_axis_angle([1,0,0], 0.5)
        >>> R_test = np_ol_g_R(g_orig)
        >>> g_recovered = np_ol_R_g(R_test)
        >>> print(np.allclose(g_orig, g_recovered))  # True
    Notes:
        - Inverse of np_ol_g_R()
        - Handles small rotations accurately
        - Returns identity for zero vector
```

---

## Function: np_ol_R_g2

**Signature**:
```python
def np_ol_R_g2(R,epsilon, delta):
```

**Description**:

Alternative Rodrigues-Frank to rotation matrix (NumPy version, method 2).
    NumPy implementation of alternative conversion method.

**Input**:

R (numpy.ndarray [3]): Rodrigues-Frank vector

**Output**:

numpy.ndarray (3×3): Rotation matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> R = np.array([0.1, 0.2, 0.3])
        >>> g = np_ol_R_g2(R)
    Notes:
        - Alternative implementation for comparison
```

---

## Function: np_ol_R_q2

**Signature**:
```python
def np_ol_R_q2(R):
```

**Description**:

Convert Rodrigues-Frank vector to quaternion (method 2).
    Transforms Rodrigues-Frank representation to unit quaternion [w,x,y,z].

**Input**:

R (numpy.ndarray [3]): Rodrigues-Frank vector

**Output**:

numpy.ndarray [4]: Unit quaternion [w, x, y, z]

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Small rotation
        >>> R = np.array([0.01, 0, 0])
        >>> q = np_ol_R_q2(R)
        >>> print(q)  # Near [1, 0.01, 0, 0]
        >>> print(f"Norm: {np.linalg.norm(q)}")  # 1.0
    Notes:
        - Output is normalized quaternion
        - w ≥ 0 convention
        - Quaternion avoids singularities of Rodrigues-Frank
    Formula:
        θ = 2·arctan(||R||)
        q = [cos(θ/2), r·sin(θ/2)]
        where r = R/||R||
```

---

## Function: np_ol_g_R

**Signature**:
```python
def np_ol_g_R(g):
```

**Description**:

Convert rotation matrix to Rodrigues-Frank vector (NumPy version).
    NumPy implementation returning Rodrigues-Frank vector as numpy array.

**Input**:

g (numpy.ndarray 3×3): Rotation matrix

**Output**:

numpy.ndarray [3]: Rodrigues-Frank vector

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Identity rotation
        >>> g_id = np.eye(3)
        >>> R_id = np_ol_g_R(g_id)
        >>> print(R_id)  # [0, 0, 0]
        >>> 
        >>> # General rotation
        >>> g = rotation_from_axis_angle([1,1,1], np.pi/3)
        >>> R = np_ol_g_R(g)
        >>> print(f"Rodrigues vector: {R}")
    Notes:
        - Efficient numpy implementation
        - Useful in grain boundary analysis
        - Common in crystal plasticity codes
```

---

## Function: np_ol_g_R2

**Signature**:
```python
def np_ol_g_R2(g,epsilon, delta):
```

**Description**:

Alternative Rodrigues-Frank conversion (NumPy version, method 2).
    NumPy implementation of alternative Rodrigues-Frank calculation.

**Input**:

g (numpy.ndarray 3×3): Rotation matrix

**Output**:

numpy.ndarray [3]: Rodrigues-Frank vector

**Usage Example**:

```python
>>> import numpy as np
        >>> g = np.eye(3)
        >>> R = np_ol_g_R2(g)
    Notes:
        - Alternative implementation
        - For numerical stability comparison
```

---

## Function: np_ol_g_q2

**Signature**:
```python
def np_ol_g_q2(g):
```

**Description**:

Convert rotation matrix to quaternion (method 2).
    Extracts unit quaternion representation from rotation matrix.

**Input**:

g (numpy.ndarray 3×3): Rotation matrix

**Output**:

numpy.ndarray [4]: Unit quaternion [w, x, y, z]

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # 90° around Z
        >>> g = np.array([[0, -1, 0],
        ...               [1,  0, 0],
        ...               [0,  0, 1]], dtype=float)
        >>> q = np_ol_g_q2(g)
        >>> print(q)  # [0.707, 0, 0, 0.707]
    Notes:
        - Uses Rodrigues-Frank as intermediate
        - Normalized output
        - Stable for all rotation angles
    Formula:
        g → R → q
        (via np_ol_g_R and np_ol_R_q2)
```

---

## Function: np_ol_g_rtheta_rad

**Signature**:
```python
def np_ol_g_rtheta_rad(g):
```

**Description**:

Convert rotation matrix to axis-angle representation (NumPy optimized version).

**Input**:

g: numpy array (3, 3) - Rotation matrix

**Output**:

r: numpy array (3,) - Rotation axis (unit vector)
        ptheta: float - Rotation angle in radians

**Usage Example**:

```python
>>> import numpy as np
        >>> # 120-degree rotation around [1,1,1]
        >>> angle = np.radians(120)
        >>> axis = np.array([1, 1, 1]) / np.sqrt(3)
        >>> # Create rotation matrix using Rodrigues formula
        >>> K = np.array([[0, -axis[2], axis[1]],
        ...               [axis[2], 0, -axis[0]],
        ...               [-axis[1], axis[0], 0]])
        >>> g = np.eye(3) + np.sin(angle)*K + (1-np.cos(angle))*K@K
        >>> 
        >>> # Convert back to axis-angle
        >>> axis_out, angle_out = np_ol_g_rtheta_rad(g)
        >>> print("Recovered axis:", axis_out)
        >>> print("Recovered angle (deg):", np.degrees(angle_out))
```

---

## Function: np_ol_q_g

**Signature**:
```python
def np_ol_q_g(q):
```

**Description**:

Convert quaternion to rotation matrix.
    Transforms unit quaternion [w,x,y,z] to 3×3 rotation matrix.

**Input**:

q (numpy.ndarray [4]): Unit quaternion [w, x, y, z]

**Output**:

numpy.ndarray (3×3): Rotation matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Identity quaternion
        >>> q_id = np.array([1, 0, 0, 0])
        >>> g_id = np_ol_q_g(q_id)
        >>> print(np.allclose(g_id, np.eye(3)))  # True
        >>> 
        >>> # 180° around X
        >>> q_180x = np.array([0, 1, 0, 0])
        >>> g_180x = np_ol_q_g(q_180x)
        >>> print(g_180x)
        >>> # [[ 1,  0,  0],
        >>> #  [ 0, -1,  0],
        >>> #  [ 0,  0, -1]]
    Notes:
        - Input should be normalized
        - Avoids gimbal lock
        - Efficient for interpolation
    Formula:
        g₁₁ = 1 - 2(y² + z²)
        g₁₂ = 2(xy - zw)
        g₁₃ = 2(xz + yw)
        ... (full 3×3 matrix)
```

---

## Function: np_ol_rtheta_g_rad

**Signature**:
```python
def np_ol_rtheta_g_rad(r, theta):
```

**Description**:

Convert axis-angle representation to rotation matrix (NumPy version).
    Uses Rodrigues' rotation formula.

**Input**:

r: numpy array (3,) - Rotation axis (unit vector)
        theta: float - Rotation angle in radians

**Output**:

g: numpy array (3, 3) - Rotation matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> # 45-degree rotation around [1,1,1] axis
        >>> axis = np.array([1, 1, 1]) / np.sqrt(3)  # Normalize
        >>> angle = np.radians(45)
        >>> g = np_ol_rtheta_g_rad(axis, angle)
        >>> print("Rotation matrix:")
        >>> print(g)
        >>> 
        >>> # Verify it's a valid rotation matrix
        >>> det = np.linalg.det(g)
        >>> print("Determinant:", det)  # Should be 1
        >>> orthogonal = g @ g.T
        >>> print("G @ G.T (should be identity):")
        >>> print(orthogonal)
```

---

## Function: np_rodrigues2gmat

**Signature**:
```python
def np_rodrigues2gmat(rodrigues):
```

**Description**:

Convert Rodrigues-Frank vector to rotation matrix.

**Input**:

rodrigues: numpy array (3,) - Rodrigues-Frank vector (axis * tan(angle/2))

**Output**:

g: numpy array (3, 3) - Rotation matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> # Rodrigues vector for 90° around Z
        >>> rod = np.array([0, 0, 1])  # tan(45°) = 1
        >>> g = np_rodrigues2gmat(rod)
        >>> print("Rotation matrix:")
        >>> print(g)
        >>> # Should give approximately [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        >>> # Round trip test
        >>> g_original = np.array([[0.5, -0.866, 0],
        ...                        [0.866, 0.5, 0],
        ...                        [0, 0, 1]])  # 60° around Z
        >>> rod = np_gmat2rodrigues(g_original)
        >>> g_recovered = np_rodrigues2gmat(rod)
        >>> print("Original equals recovered:", np.allclose(g_original, g_recovered))
```

---

## Function: nside2npix

**Signature**:
```python
def nside2npix(nside):
```

**Description**:

Calculate the number of pixels in a HEALPix map.
    HEALPix = Hierarchical Equal Area isoLatitude Pixelization.

**Input**:

nside: int - HEALPix resolution parameter (must be power of 2)

**Output**:

npix: int - Total number of pixels = 12 * nside²

**Usage Example**:

```python
>>> # Low resolution
        >>> npix_low = nside2npix(nside=1)
        >>> print("Pixels for nside=1:", npix_low)  # 12
        >>> # Medium resolution
        >>> npix_med = nside2npix(nside=16)
        >>> print("Pixels for nside=16:", npix_med)  # 3072
        >>> # High resolution
        >>> npix_high = nside2npix(nside=128)
        >>> print("Pixels for nside=128:", npix_high)  # 196608
```

---

## Function: ol_R_g

**Signature**:
```python
def ol_R_g (R):
```

**Description**:

Convert Rodrigues-Frank vector to rotation matrix (list version).
    Reconstructs rotation matrix from Rodrigues-Frank representation.

**Input**:

R (list [3]): Rodrigues-Frank vector

**Output**:

list (3×3): Rotation matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Rodrigues vector for 90° around Z
        >>> R = [0, 0, 1]  # tan(45°) = 1
        >>> g = ol_R_g(R)
        >>> # Should give 90° rotation around Z
    Notes:
        - Inverse of ol_g_R()
        - Extracts angle: θ = 2·arctan(||R||)
        - Extracts axis: r = R/||R||
    Formula:
        θ = 2·arctan(||R||)
        r = R / ||R||
        g = ol_rtheta_g_rad(r, θ)
```

---

## Function: ol_R_g2

**Signature**:
```python
def ol_R_g2(R):
```

**Description**:

Alternative Rodrigues-Frank to rotation matrix (list version, method 2).
    Second implementation of Rodrigues-Frank to rotation matrix conversion.

**Input**:

R (list [3]): Rodrigues-Frank vector

**Output**:

list (3×3): Rotation matrix

**Usage Example**:

```python
>>> R = [0, 0, 0.5]
        >>> g = ol_R_g2(R)
    Notes:
        - Alternative implementation
        - Should match ol_R_g() results
```

---

## Function: ol_g_R

**Signature**:
```python
def ol_g_R(g):
```

**Description**:

Convert rotation matrix to Rodrigues-Frank vector (list version).
    Calculates Rodrigues-Frank vector R = r·tan(θ/2) from rotation matrix,
    where r is the rotation axis and θ is the rotation angle.

**Input**:

g (list 3×3): Rotation matrix

**Output**:

list [3]: Rodrigues-Frank vector

**Usage Example**:

```python
>>> # 90° rotation around Z
        >>> g = [[0, -1, 0],
        ...      [1,  0, 0],
        ...      [0,  0, 1]]
        >>> R = ol_g_R(g)
        >>> print(R)  # [0, 0, 1] (since tan(45°)=1)
    Notes:
        - Compact representation: 3 parameters vs 9 for matrix
        - R = axis · tan(angle/2)
        - Magnitude ||R|| = tan(θ/2)
        - Direction = rotation axis
        - Singular at θ = 180° (infinite magnitude)
    Formula:
        R = r · tan(θ/2)
        where (r, θ) from ol_g_rtheta_rad(g)
```

---

## Function: ol_g_R2

**Signature**:
```python
def ol_g_R2(g):
```

**Description**:

Alternative Rodrigues-Frank conversion (list version, method 2).
    Second implementation of rotation matrix to Rodrigues-Frank vector.
    May use different numerical approach for stability.

**Input**:

g (list 3×3): Rotation matrix

**Output**:

list [3]: Rodrigues-Frank vector

**Usage Example**:

```python
>>> g = [[0, -1, 0],
        ...      [1,  0, 0],
        ...      [0,  0, 1]]
        >>> R = ol_g_R2(g)
    Notes:
        - Alternative implementation for numerical comparison
        - Should give same result as ol_g_R()
        - Useful for validation
```

---

## Function: ol_g_rtheta_rad

**Signature**:
```python
def ol_g_rtheta_rad(g):
```

**Description**:

Convert rotation matrix to axis-angle representation (Rodrigues-Frank vector).
    Returns rotation axis and angle.

**Input**:

g: list or array (3, 3) - Rotation matrix

**Output**:

r: list (3,) - Rotation axis (unit vector)
        ptheta: float - Rotation angle in radians

**Usage Example**:

```python
>>> import numpy as np
        >>> # 90-degree rotation around Z
        >>> g = [[0, -1, 0],
        ...      [1,  0, 0],
        ...      [0,  0, 1]]
        >>> axis, angle = ol_g_rtheta_rad(g)
        >>> print("Rotation axis:", axis)
        >>> print("Rotation angle (rad):", angle)
        >>> print("Rotation angle (deg):", np.degrees(angle))
        >>> # Should give axis ≈ [0, 0, 1] and angle ≈ π/2
```

---

## Function: ol_rtheta_g_rad

**Signature**:
```python
def ol_rtheta_g_rad(r, theta):
```

**Description**:

Convert axis-angle representation to rotation matrix using Rodrigues' formula.

**Input**:

r: list or array (3,) - Rotation axis (should be unit vector)
        theta: float - Rotation angle in radians

**Output**:

g: list (3, 3) - Rotation matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> # 90-degree rotation around Z axis
        >>> axis = [0, 0, 1]
        >>> angle = np.pi/2
        >>> g = ol_rtheta_g_rad(axis, angle)
        >>> print("Rotation matrix:")
        >>> for row in g:
        ...     print(row)
        >>> # Should give approximately [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        >>> # 180-degree rotation around X axis
        >>> g_180x = ol_rtheta_g_rad([1, 0, 0], np.pi)
        >>> print("180° around X:")
        >>> for row in g_180x:
        ...     print(row)
```

---

## Function: orilistMult

**Signature**:
```python
def orilistMult(Mats,Dr):
```

**Description**:

Multiply a list of orientation matrices by symmetry operations.
    Applies symmetry operations to orientation matrices to generate all
    symmetrically equivalent orientations.

**Input**:

orilist: numpy array (N, 3, 3) - Array of orientation matrices
        symops: numpy array (Ns, 3, 3) - Array of symmetry operation matrices

**Output**:

result: numpy array (N*Ns, 3, 3) - All symmetrically equivalent orientations

**Usage Example**:

```python
>>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R
        >>> 
        >>> # Create orientation matrices
        >>> N = 10
        >>> orientations = R.random(N).as_matrix()
        >>> 
        >>> # Cubic symmetry operations (simplified - just identity for demo)
        >>> symops = np.array([np.eye(3)])
        >>> 
        >>> # Generate all equivalent orientations
        >>> equiv_oris = orilistMult(orientations, symops)
        >>> print("Original:", orientations.shape)
        >>> print("With symmetry:", equiv_oris.shape)
        >>> # For full cubic symmetry (24 operations)
        >>> # equiv_oris would have shape (N*24, 3, 3)
```

---

## Function: passive_rotation

**Signature**:
```python
def passive_rotation(an, aboutaxis, deg=False):
```

**Description**:

Perform passive rotation (coordinate transformation).
    Rotates the coordinate system while vector stays fixed in space.
    v' = g^T · v (transpose of rotation matrix).

**Input**:

g (array 3×3): Rotation matrix
        v (array [3]): Vector in old coordinates

**Output**:

numpy.ndarray [3]: Vector components in new coordinates

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Rotate coordinates 90° around Z
        >>> g_z90 = np.array([[0, -1, 0],
        ...                   [1,  0, 0],
        ...                   [0,  0, 1]])
        >>> v = np.array([0, 1, 0])
        >>> v_new_coords = passive_rotation(g_z90, v)
        >>> print(v_new_coords)  # [1, 0, 0]
    Notes:
        - Passive = rotate coordinate system
        - v' = g^T · v (transpose)
        - Common in crystallography (sample↔crystal)
        - Inverse of active rotation for same g
    Formula:
        v'_i = g_{ji} v_j = g^T_{ij} v_j
```

---

## Function: pix2ang_nest

**Signature**:
```python
def pix2ang_nest(nside, ipix, pix2x, pix2y):
```

**Description**:

Convert HEALPix pixel index to spherical coordinates (theta, phi).
    Uses NESTED indexing scheme.

**Input**:

nside: int - HEALPix resolution parameter
        ipix: int - Pixel index (0 to 12*nside² - 1)
        pix2x: list - Lookup table for x-coordinate (from mk_pix2xy)
        pix2y: list - Lookup table for y-coordinate (from mk_pix2xy)

**Output**:

theta: float - Colatitude angle in radians [0, π]
        phi: float - Azimuthal angle in radians [0, 2π]

**Usage Example**:

```python
>>> # Setup lookup tables
        >>> pix2x, pix2y = mk_pix2xy()
        >>> 
        >>> # Convert pixel 0 at nside=4
        >>> nside = 4
        >>> theta, phi = pix2ang_nest(nside, ipix=0, pix2x=pix2x, pix2y=pix2y)
        >>> print(f"Pixel 0: theta={np.degrees(theta):.2f}°, phi={np.degrees(phi):.2f}°")
        >>> 
        >>> # Convert all pixels at low resolution
        >>> nside = 2
        >>> npix = nside2npix(nside)
        >>> for i in range(npix):
        ...     theta, phi = pix2ang_nest(nside, i, pix2x, pix2y)
        ...     # Process coordinates...
```

---

## Function: quat_conjugate

**Signature**:
```python
def quat_conjugate(q):
```

**Description**:

Compute the conjugate of a quaternion (inverse for unit quaternions).

**Input**:

q: numpy array (4,) - Quaternion [w, x, y, z]

**Output**:

q_conj: numpy array (4,) - Conjugate quaternion [w, -x, -y, -z]

**Usage Example**:

```python
>>> import numpy as np
        >>> q = np.array([0.707, 0.707, 0.0, 0.0])
        >>> q_conj = quat_conjugate(q)
        >>> print("Original:", q)
        >>> print("Conjugate:", q_conj)
        >>> # Conjugate: [0.707, -0.707, 0.0, 0.0]
        >>> # Verify: q * q_conj = identity
        >>> result = quat_mult(q, q_conj)
        >>> print("q * q_conj:", result)
        >>> # Should give approximately [1, 0, 0, 0]
```

---

## Function: quat_misori_deg

**Signature**:
```python
def quat_misori_deg(q1, q2):
```

**Description**:

Calculate the misorientation angle between two quaternions in degrees.
    Returns the minimum rotation angle needed to go from q1 to q2.

**Input**:

q1: numpy array (4,) - First quaternion [w, x, y, z]
        q2: numpy array (4,) - Second quaternion [w, x, y, z]

**Output**:

angle: float - Misorientation angle in degrees (0 to 180)

**Usage Example**:

```python
>>> import numpy as np
        >>> # Same orientation
        >>> q1 = np.array([1.0, 0.0, 0.0, 0.0])
        >>> angle = quat_misori_deg(q1, q1)
        >>> print("Same orientation angle:", angle)
        >>> # Should give 0 degrees
        >>> # 90-degree misorientation
        >>> q1 = np.array([1.0, 0.0, 0.0, 0.0])
        >>> q2 = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> angle = quat_misori_deg(q1, q2)
        >>> print("Misorientation:", angle, "degrees")
        >>> # Should give approximately 90 degrees
```

---

## Function: quat_mult

**Signature**:
```python
def quat_mult(q1, q2):
```

**Description**:

Multiply two quaternions (Hamilton product).
    Order matters: q1 * q2 ≠ q2 * q1 in general.

**Input**:

q1: numpy array (4,) - First quaternion [w, x, y, z]
        q2: numpy array (4,) - Second quaternion [w, x, y, z]

**Output**:

q: numpy array (4,) - Product quaternion q1 * q2

**Usage Example**:

```python
>>> import numpy as np
        >>> # Two 90-degree rotations around z-axis = 180-degree rotation
        >>> q_90 = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> q_result = quat_mult(q_90, q_90)
        >>> print("Result quaternion:", q_result)
        >>> # Should give approximately [0, 0, 0, 1] (180-degree rotation)
        >>> # Rotation composition
        >>> q_x = np.array([np.cos(np.pi/4), np.sin(np.pi/4), 0, 0])  # 90° around x
        >>> q_z = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])  # 90° around z
        >>> q_combined = quat_mult(q_x, q_z)
        >>> print("Combined rotation:", q_combined)
```

---

## Function: quat_multiply

**Signature**:
```python
def quat_multiply(q1, q2):
```

**Description**:

Quaternion multiplication (same as quat_mult, alternative implementation).

**Input**:

q1: numpy array (4,) - First quaternion [w, x, y, z]
        q2: numpy array (4,) - Second quaternion [w, x, y, z]

**Output**:

q: numpy array (4,) - Product quaternion q1 * q2

**Usage Example**:

```python
>>> import numpy as np
        >>> q1 = np.array([1.0, 0.0, 0.0, 0.0])  # Identity
        >>> q2 = np.array([0.707, 0.707, 0.0, 0.0])  # 90° around x
        >>> result = quat_multiply(q1, q2)
        >>> print("Result:", result)
        >>> # Should give q2 since q1 is identity
```

---

## Function: quat_to_mat

**Signature**:
```python
def quat_to_mat(q):
```

**Description**:

Convert a quaternion to a rotation matrix.

**Input**:

q: numpy array (4,) - Unit quaternion [w, x, y, z]

**Output**:

R: numpy array (3x3) - Rotation matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> # Quaternion for 90-degree rotation around z-axis
        >>> q = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> R = quat_to_mat(q)
        >>> print("Rotation matrix:")
        >>> print(R)
        >>> # Should give approximately [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        >>> # Identity quaternion
        >>> q_identity = np.array([1.0, 0.0, 0.0, 0.0])
        >>> R_identity = quat_to_mat(q_identity)
        >>> print("Identity matrix:")
        >>> print(R_identity)
```

---

## Function: rotation_from_axis_angle

**Signature**:
```python
def rotation_from_axis_angle(ax,an,deg=False):
```

**Description**:

Generate rotation matrix from axis-angle representation using Rodrigues formula.
    Creates 3×3 rotation matrix for rotation by 'angle' radians around 'axis'.
    Implements Rodrigues' rotation formula for arbitrary axis rotations.

**Input**:

axis (array [3]): Rotation axis (will be normalized)
        angle (float): Rotation angle in radians

**Output**:

numpy.ndarray (3×3): Rotation matrix R

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # 90° rotation around z-axis
        >>> axis = np.array([0, 0, 1])
        >>> angle = np.pi / 2
        >>> R = rotation_from_axis_angle(axis, angle)
        >>> print(R)
        >>> # [[0, -1, 0],
        >>> #  [1,  0, 0],
        >>> #  [0,  0, 1]]
        >>> 
        >>> # Verify: rotate [1,0,0] to [0,1,0]
        >>> v = np.array([1, 0, 0])
        >>> v_rot = R.dot(v)
        >>> print(v_rot)  # [0, 1, 0]
        >>> 
        >>> # 120° rotation around [111]
        >>> axis_111 = np.array([1, 1, 1])
        >>> R_111 = rotation_from_axis_angle(axis_111, 2*np.pi/3)
    Notes:
        - Axis is automatically normalized
        - Right-hand rule: thumb along axis, fingers show rotation
        - Preserves vector lengths (orthogonal matrix)
        - det(R) = 1 (proper rotation)
        - Used in crystallographic symmetry operations
    Formula (Rodrigues):
        R = I + sin(θ)K + (1-cos(θ))K²
        where K is the skew-symmetric matrix of the axis
```

---

## Function: simple_grid

**Signature**:
```python
def simple_grid(resol):
```

**Description**:

Generate a uniform grid of quaternions covering orientation space (SO(3)).
    Uses HEALPix for S² sampling and uniform sampling for S¹, combined via Hopf fibration.

**Input**:

resol: int - Resolution parameter
                     - Number of S¹ points: 2^resol * 6
                     - HEALPix nside: 2^resol
                     - Total quaternions: 12 * 4^resol * (2^resol * 6)

**Output**:

quats: list of lists - Uniformly distributed quaternions [w, x, y, z]

**Usage Example**:

```python
>>> # Low resolution grid (fast)
        >>> quats_low = simple_grid(resol=2)
        >>> print("Low resolution quaternions:", len(quats_low))
        >>> # 12 * 4^2 * (2^2 * 6) = 12 * 16 * 24 = 4608
        >>> # Medium resolution (moderate)
        >>> quats_med = simple_grid(resol=3)
        >>> print("Medium resolution quaternions:", len(quats_med))
        >>> # 12 * 4^3 * (2^3 * 6) = 12 * 64 * 48 = 36864
        >>> # Convert to rotation matrices for use
        >>> import numpy as np
        >>> q_sample = np.array(quats_low[0])
        >>> R = quat_to_mat(q_sample)
        >>> print("Sample rotation matrix:")
        >>> print(R)
        >>> # Use for texture analysis, pole figures, etc.
        >>> # Each quaternion represents a unique crystal orientation
```

---

## Function: symmetry_elements

**Signature**:
```python
def symmetry_elements(lattice):
```

**Description**:

Generate symmetry operation matrices for crystal system.
    Returns list of 3×3 rotation matrices representing all symmetry
    operations for specified crystal system.

**Input**:

crystal_system (str): 'cubic', 'hexagonal', 'tetragonal', 
                             'orthorhombic', 'monoclinic', 'triclinic'

**Output**:

list: List of numpy.ndarray (3×3) symmetry matrices

**Usage Example**:

```python
>>> symops = symmetry_elements('cubic')
        >>> print(f"Cubic has {len(symops)} symmetry operations")
        >>> # 24 operations for cubic (point group m-3m)
        >>> 
        >>> # Verify they're proper rotations
        >>> for g in symops:
        ...     assert np.abs(np.linalg.det(g) - 1.0) < 1e-10
    Notes:
        - Cubic (Oh): 24 operations
        - Hexagonal (D6h): 12 operations (typically)
        - Tetragonal (D4h): 8 operations
        - All matrices are proper rotations (det=1)
        - Used in texture analysis and pole figures
```

---

## Function: symmetry_reduced_oris

**Signature**:
```python
def symmetry_reduced_oris(umatsa, symops):
```

**Description**:

Reduce orientations to fundamental zone using crystal symmetry.
    Applies symmetry operations to bring all orientations into the
    fundamental zone (asymmetric unit) of orientation space. This
    ensures unique representation of each orientation.

**Input**:

oris: numpy array (N, 3, 3) - Array of orientation matrices
        symops: numpy array (Ns, 3, 3) - Crystal symmetry operations

**Output**:

reduced_oris: numpy array (N, 3, 3) - Orientations in fundamental zone

**Usage Example**:

```python
>>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R
        >>> 
        >>> # Generate random orientations
        >>> N = 100
        >>> orientations = R.random(N).as_matrix()
        >>> 
        >>> # Cubic symmetry operations
        >>> symops = symmetry_elements('cubic')
        >>> 
        >>> # Reduce to fundamental zone
        >>> reduced = symmetry_reduced_oris(orientations, symops)
        >>> print("Original shape:", orientations.shape)
        >>> print("Reduced shape:", reduced.shape)
        >>> 
        >>> # All reduced orientations are now in fundamental zone
        >>> # Check Euler angle ranges after reduction
        >>> euler_reduced = np.array([euler_angles_from_matrix(R) 
        ...                           for R in reduced])
        >>> print("Phi1 range:", np.degrees(euler_reduced[:, 0]).min(), 
        ...       np.degrees(euler_reduced[:, 0]).max())
```

---

## Function: symposMult

**Signature**:
```python
def symposMult(sympos,Mats):
```

**Description**:

Multiply symmetry operations to generate composite symmetry operations.
    Computes the product of two sets of symmetry operations to generate
    all possible combinations.

**Input**:

sympos1: numpy array (N1, 3, 3) - First set of symmetry operations
        sympos2: numpy array (N2, 3, 3) - Second set of symmetry operations

**Output**:

result: numpy array (N1*N2, 3, 3) - All products of sympos1 and sympos2

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Define some symmetry operations
        >>> # 90° rotation around Z
        >>> rot_z = np.array([[0, -1, 0],
        ...                   [1,  0, 0],
        ...                   [0,  0, 1]])
        >>> 
        >>> # Mirror in XY plane
        >>> mirror_xy = np.array([[1, 0, 0],
        ...                       [0, 1, 0],
        ...                       [0, 0, -1]])
        >>> 
        >>> symops1 = np.array([np.eye(3), rot_z])
        >>> symops2 = np.array([np.eye(3), mirror_xy])
        >>> 
        >>> # Generate all combinations
        >>> combined = symposMult(symops1, symops2)
        >>> print("Combined symmetry operations:", combined.shape)
        >>> # Shape: (4, 3, 3) = 2*2 combinations
```

---

## Function: symposMult02

**Signature**:
```python
def symposMult02(sympos,Mats):
```

**Description**:

Alternative implementation of symmetry operations multiplication.
    Similar to symposMult but may use different algorithm or ordering.
    Generates all products of two symmetry operation sets.

**Input**:

symops1: numpy array (N1, 3, 3) - First set of symmetry operations
        symops2: numpy array (N2, 3, 3) - Second set of symmetry operations

**Output**:

result: numpy array (N1*N2, 3, 3) - Product symmetry operations

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Point group symmetry operations
        >>> identity = np.eye(3)
        >>> inversion = -np.eye(3)
        >>> 
        >>> group1 = np.array([identity])
        >>> group2 = np.array([identity, inversion])
        >>> 
        >>> # Combine to generate centrosymmetric group
        >>> result = symposMult02(group1, group2)
        >>> print("Result shape:", result.shape)
        >>> # Use for crystallographic point groups
        >>> # e.g., combining rotation and inversion symmetries
```

---

## Function: trace_to_angle

**Signature**:
```python
def trace_to_angle(tr, out="deg"):
```

**Description**:

Convert rotation matrix trace to rotation angle.
    Uses the trace (sum of diagonal elements) of a rotation matrix to
    calculate the rotation angle. Based on the relation:
    trace(R) = 1 + 2*cos(theta)

**Input**:

trace: float - Trace of rotation matrix (sum of diagonal elements)

**Output**:

angle: float - Rotation angle in degrees

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Create rotation matrix
        >>> theta_original = np.radians(45)
        >>> axis = np.array([0, 0, 1])
        >>> R = np_ol_rtheta_g_rad(axis, theta_original)
        >>> 
        >>> # Get trace and convert to angle
        >>> tr = np.trace(R)
        >>> angle = trace_to_angle(tr)
        >>> print(f"Original angle: {np.degrees(theta_original):.2f}°")
        >>> print(f"Recovered angle: {angle:.2f}°")
        >>> # Works for any rotation matrix
        >>> R_euler = np_euler_matrix(np.pi/3, np.pi/4, np.pi/6)
        >>> angle_euler = trace_to_angle(np.trace(R_euler))
        >>> print(f"Rotation angle from Euler: {angle_euler:.2f}°")
        >>> # Verify formula: trace = 1 + 2*cos(theta)
        >>> theta_calc = np.arccos((tr - 1) / 2)
        >>> print(f"Direct calculation: {np.degrees(theta_calc):.2f}°")
```

---


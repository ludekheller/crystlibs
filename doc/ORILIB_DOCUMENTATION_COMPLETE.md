# ORILIB - COMPLETE DOCUMENTATION
## Orientation Analysis and Transformations Library

**Module**: `orilib.py` (Extended Version)  
**Total Functions**: 43  
**Version**: Extended 2024  
**Status**: Production Ready  

---

## 📚 TABLE OF CONTENTS

1. [mat_to_quat](#mat_to_quat)
2. [quat_to_mat](#quat_to_mat)
3. [quat_mult](#quat_mult)
4. [quat_misori_deg](#quat_misori_deg)
5. [misori_sym_deg_sample_to_crystal_fast](#misori_sym_deg_sample_to_crystal_fast)
6. [quat_conjugate](#quat_conjugate)
7. [quat_multiply](#quat_multiply)
8. [misori_sym_deg_quats](#misori_sym_deg_quats)
9. [eu2quat](#eu2quat)
10. [np_euler_matrix](#np_euler_matrix)
11. [np_eulers_matrices](#np_eulers_matrices)
12. [np_inverse_euler_matrix](#np_inverse_euler_matrix)
13. [ol_g_rtheta_rad](#ol_g_rtheta_rad)
14. [np_ol_g_rtheta_rad](#np_ol_g_rtheta_rad)
15. [ol_rtheta_g_rad](#ol_rtheta_g_rad)
16. [np_ol_rtheta_g_rad](#np_ol_rtheta_g_rad)
17. [np_gmat2rodrigues](#np_gmat2rodrigues)
18. [np_rodrigues2gmat](#np_rodrigues2gmat)
19. [np_g2quats](#np_g2quats)
20. [Qlog](#qlog)
21. [Qproduct](#qproduct)
22. [QMatproduct](#qmatproduct)
23. [grid_s1](#grid_s1)
24. [hopf2quat](#hopf2quat)
25. [nside2npix](#nside2npix)
26. [pix2ang_nest](#pix2ang_nest)
27. [mk_pix2xy](#mk_pix2xy)
28. [simple_grid](#simple_grid)
29. [rotation_from_axis_angle](#rotation_from_axis_angle)
30. [ol_g_R](#ol_g_r)
31. [np_ol_g_R](#np_ol_g_r)
32. [ol_R_g](#ol_r_g)
33. [np_ol_R_g](#np_ol_r_g)
34. [ol_g_R2](#ol_g_r2)
35. [np_ol_g_R2](#np_ol_g_r2)
36. [ol_R_g2](#ol_r_g2)
37. [np_ol_R_g2](#np_ol_r_g2)
38. [np_ol_R_q2](#np_ol_r_q2)
39. [np_ol_g_q2](#np_ol_g_q2)
40. [np_ol_q_g](#np_ol_q_g)
41. [active_rotation](#active_rotation)
42. [passive_rotation](#passive_rotation)
43. [euler_angles_reduction](#euler_angles_reduction)

---

## 1. `mat_to_quat`

### Signature
```python
def mat_to_quat(R)
```

### Description
Convert a rotation matrix to a quaternion representation.
    Uses the Shepperd's method for numerical stability.

### Input Parameters
R: numpy array (3x3) - Rotation matrix (orthogonal with det=1)

### Output
q: numpy array (4,) - Unit quaternion [w, x, y, z] where w is the scalar part

### Usage Examples
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

## 2. `quat_to_mat`

### Signature
```python
def quat_to_mat(q)
```

### Description
Convert a quaternion to a rotation matrix.

### Input Parameters
q: numpy array (4,) - Unit quaternion [w, x, y, z]

### Output
R: numpy array (3x3) - Rotation matrix

### Usage Examples
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

## 3. `quat_mult`

### Signature
```python
def quat_mult(q1, q2)
```

### Description
Multiply two quaternions (Hamilton product).
    Order matters: q1 * q2 ≠ q2 * q1 in general.

### Input Parameters
q1: numpy array (4,) - First quaternion [w, x, y, z]
        q2: numpy array (4,) - Second quaternion [w, x, y, z]

### Output
q: numpy array (4,) - Product quaternion q1 * q2

### Usage Examples
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

## 4. `quat_misori_deg`

### Signature
```python
def quat_misori_deg(q1, q2)
```

### Description
Calculate the misorientation angle between two quaternions in degrees.
    Returns the minimum rotation angle needed to go from q1 to q2.

### Input Parameters
q1: numpy array (4,) - First quaternion [w, x, y, z]
        q2: numpy array (4,) - Second quaternion [w, x, y, z]

### Output
angle: float - Misorientation angle in degrees (0 to 180)

### Usage Examples
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

## 5. `misori_sym_deg_sample_to_crystal_fast`

### Signature
```python
def misori_sym_deg_sample_to_crystal_fast(M1, M2, symops)
```

### Description
Compute misorientation angles (deg) between orientations M1 and M2
    considering crystal symmetry operations using vectorized operations.
    This is the fast version optimized with Numba parallel processing.

### Input Parameters
M1: numpy array (N, 3, 3) - Orientation matrices (sample→crystal reference frame)
        M2: numpy array (N, 3, 3) - Orientation matrices (sample→crystal reference frame)
                                     Must have same length as M1
        symops: numpy array (Ns, 3, 3) - Crystal symmetry operation matrices

### Output
miso: numpy array (N,) - Minimum misorientation angle in degrees for each pair

### Usage Examples
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

## 6. `quat_conjugate`

### Signature
```python
def quat_conjugate(q)
```

### Description
Compute the conjugate of a quaternion (inverse for unit quaternions).

### Input Parameters
q: numpy array (4,) - Quaternion [w, x, y, z]

### Output
q_conj: numpy array (4,) - Conjugate quaternion [w, -x, -y, -z]

### Usage Examples
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

## 7. `quat_multiply`

### Signature
```python
def quat_multiply(q1, q2)
```

### Description
Quaternion multiplication (same as quat_mult, alternative implementation).

### Input Parameters
q1: numpy array (4,) - First quaternion [w, x, y, z]
        q2: numpy array (4,) - Second quaternion [w, x, y, z]

### Output
q: numpy array (4,) - Product quaternion q1 * q2

### Usage Examples
```python
>>> import numpy as np
        >>> q1 = np.array([1.0, 0.0, 0.0, 0.0])  # Identity
        >>> q2 = np.array([0.707, 0.707, 0.0, 0.0])  # 90° around x
        >>> result = quat_multiply(q1, q2)
        >>> print("Result:", result)
        >>> # Should give q2 since q1 is identity
```

---

## 8. `misori_sym_deg_quats`

### Signature
```python
def misori_sym_deg_quats(q1, q2, sym_quats)
```

### Description
Calculate minimum misorientation angle between two quaternions considering
    crystal symmetry operations. This is a Numba-optimized version.

### Input Parameters
q1: numpy array (4,) - First quaternion [w, x, y, z]
        q2: numpy array (4,) - Second quaternion [w, x, y, z]
        sym_quats: numpy array (Ns, 4) - Array of symmetry operation quaternions

### Output
min_angle: float - Minimum misorientation angle in degrees

### Usage Examples
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

## 9. `eu2quat`

### Signature
```python
def eu2quat(phi1, Phi, phi2)
```

### Description
Convert Bunge Euler angles to quaternion representation.
    Euler angles follow the ZXZ convention (Bunge notation).

### Input Parameters
phi1: float - First Euler angle (rotation around Z) in radians
        Phi: float - Second Euler angle (rotation around X') in radians
        phi2: float - Third Euler angle (rotation around Z'') in radians

### Output
q: numpy array (4,) - Quaternion [w, x, y, z] with positive w

### Usage Examples
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

## 10. `np_euler_matrix`

### Signature
```python
def np_euler_matrix(ai, aj, ak)
```

### Description
Convert Euler angles to rotation matrix (ZXZ convention, Bunge notation).
    Single orientation version.

### Input Parameters
ai: float - First Euler angle (phi1) in radians
        aj: float - Second Euler angle (Phi) in radians
        ak: float - Third Euler angle (phi2) in radians

### Output
g: numpy array (3, 3) - Rotation matrix

### Usage Examples
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

## 11. `np_eulers_matrices`

### Signature
```python
def np_eulers_matrices(data, deg=False)
```

### Description
Convert multiple Euler angles to rotation matrices (vectorized).
    Processes an array of Euler angle triplets efficiently.

### Input Parameters
data: numpy array (N, 3) - Array of Euler angles [phi1, Phi, phi2]
        deg: bool - If True, input angles are in degrees; if False, radians (default: False)

### Output
g: numpy array (N, 3, 3) - Array of rotation matrices

### Usage Examples
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

## 12. `np_inverse_euler_matrix`

### Signature
```python
def np_inverse_euler_matrix(ai, aj, ak)
```

### Description
Convert Euler angles to inverse (transpose) rotation matrix.
    Equivalent to the transpose of the forward rotation matrix.

### Input Parameters
ai: float - First Euler angle (phi1) in radians
        aj: float - Second Euler angle (Phi) in radians
        ak: float - Third Euler angle (phi2) in radians

### Output
U: numpy array (3, 3) - Inverse rotation matrix (transpose of forward matrix)

### Usage Examples
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

## 13. `ol_g_rtheta_rad`

### Signature
```python
def ol_g_rtheta_rad(g)
```

### Description
Convert rotation matrix to axis-angle representation (Rodrigues-Frank vector).
    Returns rotation axis and angle.

### Input Parameters
g: list or array (3, 3) - Rotation matrix

### Output
r: list (3,) - Rotation axis (unit vector)
        ptheta: float - Rotation angle in radians

### Usage Examples
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

## 14. `np_ol_g_rtheta_rad`

### Signature
```python
def np_ol_g_rtheta_rad(g)
```

### Description
Convert rotation matrix to axis-angle representation (NumPy optimized version).

### Input Parameters
g: numpy array (3, 3) - Rotation matrix

### Output
r: numpy array (3,) - Rotation axis (unit vector)
        ptheta: float - Rotation angle in radians

### Usage Examples
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

## 15. `ol_rtheta_g_rad`

### Signature
```python
def ol_rtheta_g_rad(r, theta)
```

### Description
Convert axis-angle representation to rotation matrix using Rodrigues' formula.

### Input Parameters
r: list or array (3,) - Rotation axis (should be unit vector)
        theta: float - Rotation angle in radians

### Output
g: list (3, 3) - Rotation matrix

### Usage Examples
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

## 16. `np_ol_rtheta_g_rad`

### Signature
```python
def np_ol_rtheta_g_rad(r, theta)
```

### Description
Convert axis-angle representation to rotation matrix (NumPy version).
    Uses Rodrigues' rotation formula.

### Input Parameters
r: numpy array (3,) - Rotation axis (unit vector)
        theta: float - Rotation angle in radians

### Output
g: numpy array (3, 3) - Rotation matrix

### Usage Examples
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

## 17. `np_gmat2rodrigues`

### Signature
```python
def np_gmat2rodrigues(g)
```

### Description
Convert rotation matrix to Rodrigues-Frank vector representation.
    Rodrigues vector = rotation_axis * tan(angle/2)

### Input Parameters
g: numpy array (3, 3) - Rotation matrix

### Output
rodrigues: numpy array (3,) - Rodrigues-Frank vector

### Usage Examples
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

## 18. `np_rodrigues2gmat`

### Signature
```python
def np_rodrigues2gmat(rodrigues)
```

### Description
Convert Rodrigues-Frank vector to rotation matrix.

### Input Parameters
rodrigues: numpy array (3,) - Rodrigues-Frank vector (axis * tan(angle/2))

### Output
g: numpy array (3, 3) - Rotation matrix

### Usage Examples
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

## 19. `np_g2quats`

### Signature
```python
def np_g2quats(umatsa)
```

### Description
Convert multiple rotation matrices to quaternions (vectorized).
    Handles arrays of rotation matrices efficiently.

### Input Parameters
umatsa: numpy array (N, 3, 3) - Array of rotation matrices

### Output
Q: numpy array (4, N) - Array of quaternions [w, x, y, z] for each matrix

### Usage Examples
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

## 20. `Qlog`

### Signature
```python
def Qlog(QM)
```

### Description
Compute the logarithm of quaternions (quaternion logarithm map).
    Maps quaternions from unit sphere S³ to tangent space at identity.

### Input Parameters
QM: numpy array (4, N, M) - Array of quaternions [w, x, y, z, ...]

### Output
qlog: numpy array (4, N, M) - Logarithm of quaternions
                                       qlog[0] = 0, qlog[1:4] = v * angle

### Usage Examples
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

## 21. `Qproduct`

### Signature
```python
def Qproduct(P, Q)
```

### Description
Compute quaternion product for arrays of quaternions.
    Computes P * Q for each pair of quaternions.

### Input Parameters
P: numpy array (4, N) - First array of quaternions [w, x, y, z]
        Q: numpy array (4, N) - Second array of quaternions [w, x, y, z]

### Output
result: numpy array (4, N) - Product quaternions P * Q

### Usage Examples
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

## 22. `QMatproduct`

### Signature
```python
def QMatproduct(sym, Q)
```

### Description
Multiply a single symmetry quaternion with multiple quaternions.
    Applies the same symmetry operation to an array of orientations.

### Input Parameters
sym: numpy array (4,) - Single symmetry quaternion [w, x, y, z]
        Q: numpy array (4, N) - Array of quaternions to transform

### Output
SQ: numpy array (4, N) - Transformed quaternions (sym * Q)

### Usage Examples
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

## 23. `grid_s1`

### Signature
```python
def grid_s1(resol, grids=6)
```

### Description
Generate uniformly distributed points on S¹ (circle).
    Used for sampling the third Euler angle in Hopf coordinates.

### Input Parameters
resol: int - Resolution parameter (number of points = 2^resol * grids)
        grids: int - Grid multiplier (default: 6)

### Output
points: list of floats - Angles in radians from 0 to 2π

### Usage Examples
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

## 24. `hopf2quat`

### Signature
```python
def hopf2quat(Points)
```

### Description
Convert Hopf coordinates to quaternions.
    Hopf coordinates (θ, φ, ψ) parameterize the unit quaternion sphere S³.

### Input Parameters
Points: list of tuples - Each tuple contains (theta, phi, psi) in radians
                                 theta ∈ [0, π], phi ∈ [0, 2π], psi ∈ [0, 2π]

### Output
quats: list of lists - Quaternions [w, x, y, z] for each point

### Usage Examples
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

## 25. `nside2npix`

### Signature
```python
def nside2npix(nside)
```

### Description
Calculate the number of pixels in a HEALPix map.
    HEALPix = Hierarchical Equal Area isoLatitude Pixelization.

### Input Parameters
nside: int - HEALPix resolution parameter (must be power of 2)

### Output
npix: int - Total number of pixels = 12 * nside²

### Usage Examples
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

## 26. `pix2ang_nest`

### Signature
```python
def pix2ang_nest(nside, ipix, pix2x, pix2y)
```

### Description
Convert HEALPix pixel index to spherical coordinates (theta, phi).
    Uses NESTED indexing scheme.

### Input Parameters
nside: int - HEALPix resolution parameter
        ipix: int - Pixel index (0 to 12*nside² - 1)
        pix2x: list - Lookup table for x-coordinate (from mk_pix2xy)
        pix2y: list - Lookup table for y-coordinate (from mk_pix2xy)

### Output
theta: float - Colatitude angle in radians [0, π]
        phi: float - Azimuthal angle in radians [0, 2π]

### Usage Examples
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

## 27. `mk_pix2xy`

### Signature
```python
def mk_pix2xy()
```

### Description
Create lookup tables for HEALPix pixel index to (x,y) coordinate conversion.
    These tables are used by pix2ang_nest function.

### Input Parameters
None

### Output
pix2x: list (1024,) - X-coordinate lookup table
        pix2y: list (1024,) - Y-coordinate lookup table

### Usage Examples
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

## 28. `simple_grid`

### Signature
```python
def simple_grid(resol)
```

### Description
Generate a uniform grid of quaternions covering orientation space (SO(3)).
    Uses HEALPix for S² sampling and uniform sampling for S¹, combined via Hopf fibration.

### Input Parameters
resol: int - Resolution parameter
                     - Number of S¹ points: 2^resol * 6
                     - HEALPix nside: 2^resol
                     - Total quaternions: 12 * 4^resol * (2^resol * 6)

### Output
quats: list of lists - Uniformly distributed quaternions [w, x, y, z]

### Usage Examples
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

## 29. `rotation_from_axis_angle`

### Signature
```python
def rotation_from_axis_angle(axis, angle)
```

### Description
Generate rotation matrix from axis-angle representation using Rodrigues formula.
    Creates 3×3 rotation matrix for rotation by 'angle' radians around 'axis'.
    Implements Rodrigues' rotation formula for arbitrary axis rotations.

### Input Parameters
axis (array [3]): Rotation axis (will be normalized)
        angle (float): Rotation angle in radians

### Output
numpy.ndarray (3×3): Rotation matrix R

### Usage Examples
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
```

### Notes
- Axis is automatically normalized
        - Right-hand rule: thumb along axis, fingers show rotation
        - Preserves vector lengths (orthogonal matrix)
        - det(R) = 1 (proper rotation)
        - Used in crystallographic symmetry operations
    Formula (Rodrigues):
        R = I + sin(θ)K + (1-cos(θ))K²
        where K is the skew-symmetric matrix of the axis

---

## 30. `ol_g_R`

### Signature
```python
def ol_g_R(g)
```

### Description
Convert rotation matrix to Rodrigues-Frank vector (list version).
    Calculates Rodrigues-Frank vector R = r·tan(θ/2) from rotation matrix,
    where r is the rotation axis and θ is the rotation angle.

### Input Parameters
g (list 3×3): Rotation matrix

### Output
list [3]: Rodrigues-Frank vector

### Usage Examples
```python
>>> # 90° rotation around Z
        >>> g = [[0, -1, 0],
        ...      [1,  0, 0],
        ...      [0,  0, 1]]
        >>> R = ol_g_R(g)
        >>> print(R)  # [0, 0, 1] (since tan(45°)=1)
```

### Notes
- Compact representation: 3 parameters vs 9 for matrix
        - R = axis · tan(angle/2)
        - Magnitude ||R|| = tan(θ/2)
        - Direction = rotation axis
        - Singular at θ = 180° (infinite magnitude)

### Formula
R = r · tan(θ/2)
        where (r, θ) from ol_g_rtheta_rad(g)

---

## 31. `np_ol_g_R`

### Signature
```python
def np_ol_g_R(g)
```

### Description
Convert rotation matrix to Rodrigues-Frank vector (NumPy version).
    NumPy implementation returning Rodrigues-Frank vector as numpy array.

### Input Parameters
g (numpy.ndarray 3×3): Rotation matrix

### Output
numpy.ndarray [3]: Rodrigues-Frank vector

### Usage Examples
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
```

### Notes
- Efficient numpy implementation
        - Useful in grain boundary analysis
        - Common in crystal plasticity codes

---

## 32. `ol_R_g`

### Signature
```python
def ol_R_g(R)
```

### Description
Convert Rodrigues-Frank vector to rotation matrix (list version).
    Reconstructs rotation matrix from Rodrigues-Frank representation.

### Input Parameters
R (list [3]): Rodrigues-Frank vector

### Output
list (3×3): Rotation matrix

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # Rodrigues vector for 90° around Z
        >>> R = [0, 0, 1]  # tan(45°) = 1
        >>> g = ol_R_g(R)
        >>> # Should give 90° rotation around Z
```

### Notes
- Inverse of ol_g_R()
        - Extracts angle: θ = 2·arctan(||R||)
        - Extracts axis: r = R/||R||

### Formula
θ = 2·arctan(||R||)
        r = R / ||R||
        g = ol_rtheta_g_rad(r, θ)

---

## 33. `np_ol_R_g`

### Signature
```python
def np_ol_R_g(R)
```

### Description
Convert Rodrigues-Frank vector to rotation matrix (NumPy version).
    NumPy implementation of Rodrigues-Frank to rotation matrix conversion.

### Input Parameters
R (numpy.ndarray [3]): Rodrigues-Frank vector

### Output
numpy.ndarray (3×3): Rotation matrix

### Usage Examples
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
```

### Notes
- Inverse of np_ol_g_R()
        - Handles small rotations accurately
        - Returns identity for zero vector

---

## 34. `ol_g_R2`

### Signature
```python
def ol_g_R2(g)
```

### Description
Alternative Rodrigues-Frank conversion (list version, method 2).
    Second implementation of rotation matrix to Rodrigues-Frank vector.
    May use different numerical approach for stability.

### Input Parameters
g (list 3×3): Rotation matrix

### Output
list [3]: Rodrigues-Frank vector

### Usage Examples
```python
>>> g = [[0, -1, 0],
        ...      [1,  0, 0],
        ...      [0,  0, 1]]
        >>> R = ol_g_R2(g)
```

### Notes
- Alternative implementation for numerical comparison
        - Should give same result as ol_g_R()
        - Useful for validation

---

## 35. `np_ol_g_R2`

### Signature
```python
def np_ol_g_R2(g)
```

### Description
Alternative Rodrigues-Frank conversion (NumPy version, method 2).
    NumPy implementation of alternative Rodrigues-Frank calculation.

### Input Parameters
g (numpy.ndarray 3×3): Rotation matrix

### Output
numpy.ndarray [3]: Rodrigues-Frank vector

### Usage Examples
```python
>>> import numpy as np
        >>> g = np.eye(3)
        >>> R = np_ol_g_R2(g)
```

### Notes
- Alternative implementation
        - For numerical stability comparison

---

## 36. `ol_R_g2`

### Signature
```python
def ol_R_g2(R)
```

### Description
Alternative Rodrigues-Frank to rotation matrix (list version, method 2).
    Second implementation of Rodrigues-Frank to rotation matrix conversion.

### Input Parameters
R (list [3]): Rodrigues-Frank vector

### Output
list (3×3): Rotation matrix

### Usage Examples
```python
>>> R = [0, 0, 0.5]
        >>> g = ol_R_g2(R)
```

### Notes
- Alternative implementation
        - Should match ol_R_g() results

---

## 37. `np_ol_R_g2`

### Signature
```python
def np_ol_R_g2(R)
```

### Description
Alternative Rodrigues-Frank to rotation matrix (NumPy version, method 2).
    NumPy implementation of alternative conversion method.

### Input Parameters
R (numpy.ndarray [3]): Rodrigues-Frank vector

### Output
numpy.ndarray (3×3): Rotation matrix

### Usage Examples
```python
>>> import numpy as np
        >>> R = np.array([0.1, 0.2, 0.3])
        >>> g = np_ol_R_g2(R)
```

### Notes
- Alternative implementation for comparison

---

## 38. `np_ol_R_q2`

### Signature
```python
def np_ol_R_q2(R)
```

### Description
Convert Rodrigues-Frank vector to quaternion (method 2).
    Transforms Rodrigues-Frank representation to unit quaternion [w,x,y,z].

### Input Parameters
R (numpy.ndarray [3]): Rodrigues-Frank vector

### Output
numpy.ndarray [4]: Unit quaternion [w, x, y, z]

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # Small rotation
        >>> R = np.array([0.01, 0, 0])
        >>> q = np_ol_R_q2(R)
        >>> print(q)  # Near [1, 0.01, 0, 0]
        >>> print(f"Norm: {np.linalg.norm(q)}")  # 1.0
```

### Notes
- Output is normalized quaternion
        - w ≥ 0 convention
        - Quaternion avoids singularities of Rodrigues-Frank

### Formula
θ = 2·arctan(||R||)
        q = [cos(θ/2), r·sin(θ/2)]
        where r = R/||R||

---

## 39. `np_ol_g_q2`

### Signature
```python
def np_ol_g_q2(g)
```

### Description
Convert rotation matrix to quaternion (method 2).
    Extracts unit quaternion representation from rotation matrix.

### Input Parameters
g (numpy.ndarray 3×3): Rotation matrix

### Output
numpy.ndarray [4]: Unit quaternion [w, x, y, z]

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # 90° around Z
        >>> g = np.array([[0, -1, 0],
        ...               [1,  0, 0],
        ...               [0,  0, 1]], dtype=float)
        >>> q = np_ol_g_q2(g)
        >>> print(q)  # [0.707, 0, 0, 0.707]
```

### Notes
- Uses Rodrigues-Frank as intermediate
        - Normalized output
        - Stable for all rotation angles

### Formula
g → R → q
        (via np_ol_g_R and np_ol_R_q2)

---

## 40. `np_ol_q_g`

### Signature
```python
def np_ol_q_g(q)
```

### Description
Convert quaternion to rotation matrix.
    Transforms unit quaternion [w,x,y,z] to 3×3 rotation matrix.

### Input Parameters
q (numpy.ndarray [4]): Unit quaternion [w, x, y, z]

### Output
numpy.ndarray (3×3): Rotation matrix

### Usage Examples
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
```

### Notes
- Input should be normalized
        - Avoids gimbal lock
        - Efficient for interpolation

### Formula
g₁₁ = 1 - 2(y² + z²)
        g₁₂ = 2(xy - zw)
        g₁₃ = 2(xz + yw)
        ... (full 3×3 matrix)

---

## 41. `active_rotation`

### Signature
```python
def active_rotation(g, v)
```

### Description
Perform active rotation of vector v by rotation matrix g.
    Rotates the vector itself while keeping coordinate system fixed.
    v' = g · v (matrix-vector multiplication).

### Input Parameters
g (array 3×3): Rotation matrix
        v (array [3]): Vector to rotate

### Output
numpy.ndarray [3]: Rotated vector v'

### Usage Examples
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
```

### Notes
- Active = rotate the object
        - Compare with passive_rotation (rotate coordinates)
        - Common in physics and mechanics
        - v' = g · v

### Formula
v'_i = g_{ij} v_j

---

## 42. `passive_rotation`

### Signature
```python
def passive_rotation(g, v)
```

### Description
Perform passive rotation (coordinate transformation).
    Rotates the coordinate system while vector stays fixed in space.
    v' = g^T · v (transpose of rotation matrix).

### Input Parameters
g (array 3×3): Rotation matrix
        v (array [3]): Vector in old coordinates

### Output
numpy.ndarray [3]: Vector components in new coordinates

### Usage Examples
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
```

### Notes
- Passive = rotate coordinate system
        - v' = g^T · v (transpose)
        - Common in crystallography (sample↔crystal)
        - Inverse of active rotation for same g

### Formula
v'_i = g_{ji} v_j = g^T_{ij} v_j

---

## 43. `euler_angles_reduction`

### Signature
```python
def euler_angles_reduction(phi1, Phi, phi2, symops)
```

### Description
Reduce Euler angles to fundamental zone using symmetry operations.
    Applies crystal symmetry operations to find equivalent orientation
    with Euler angles in the fundamental zone (asymmetric unit).

### Input Parameters
phi1, Phi, phi2 (float): Euler angles in radians (Bunge convention)
        symops (list): List of 3×3 symmetry operation matrices

### Output
tuple: (phi1_red, Phi_red, phi2_red) - Reduced Euler angles in radians

### Usage Examples
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
```

### Notes
- Fundamental zone depends on crystal symmetry
        - Cubic: 0≤φ1≤90°, 0≤Φ≤45°, 0≤φ2≤90°
        - Reduces orientation distribution function (ODF) storage
        - Essential for texture analysis

---


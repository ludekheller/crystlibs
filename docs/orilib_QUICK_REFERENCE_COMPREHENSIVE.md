# ORILIB - Comprehensive Quick Reference

**Module**: `orilib.py`  
**Functions**: 59  
**Last Updated**: December 12, 2025

This reference provides function signatures with brief descriptions and example usage patterns.

---

## Mat2Quat

```python
def Mat2Quat(umatsa): #orientation matrix to quaternion
```

Convert rotation matrix to quaternion representation.
    General matrix to quaternion conversion. Alternative implementation
    that may use different numerical approach than mat_to_quat.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Random rotation matrix
        >>> from scipy.spatial.transform import Rotation as Rot
# ...
```

---

## Mat2Quat_ini

```python
def Mat2Quat_ini(umatsa): #orientation matrix to quaternion
```

Initialize matrix to quaternion conversion.
    Preliminary/initialization version of rotation matrix to quaternion
    conversion. May be an earlier implementation or setup function.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Simple 90° rotation around Z
        >>> R = np.array([[0, -1, 0],
# ...
```

---

## QMatproduct

```python
def QMatproduct(sym, Q):
```

Multiply a single symmetry quaternion with multiple quaternions.
    Applies the same symmetry operation to an array of orientations.

**Example**:
```python
>>> import numpy as np
        >>> # Define a 90° rotation symmetry around Z
        >>> sym = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> 
# ...
```

---

## Qlog

```python
def Qlog(QM):
```

Compute the logarithm of quaternions (quaternion logarithm map).
    Maps quaternions from unit sphere S³ to tangent space at identity.

**Example**:
```python
>>> import numpy as np
        >>> # Single quaternion for 90° rotation
        >>> q = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> # Reshape for function (4, 1, 1)
# ...
```

---

## Qproduct

```python
def Qproduct(P, Q):
```

Compute quaternion product for arrays of quaternions.
    Computes P * Q for each pair of quaternions.

**Example**:
```python
>>> import numpy as np
        >>> # Multiple quaternion pairs
        >>> N = 100
        >>> P = np.random.randn(4, N)
# ...
```

---

## active_rotation

```python
def active_rotation(an, aboutaxis, deg=False):
```

Perform active rotation of vector v by rotation matrix g.
    Rotates the vector itself while keeping coordinate system fixed.
    v' = g · v (matrix-vector multiplication).

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # 90° rotation around Z (active)
        >>> g_z90 = np.array([[0, -1, 0],
# ...
```

---

## disorimat

```python
def disorimat(umatsa,symops,prnt=False,withfirst=False,eqmats=False):
```

Calculate disorientation matrix considering crystal symmetry.
    Computes the minimum misorientation between two orientations by
    considering all symmetrically equivalent variants. This is the

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # EBSD grain orientations
        >>> grain1_euler = np.radians([120, 45, 80])
# ...
```

---

## disorimat_ini

```python
def disorimat_ini(umatsa,symops):
```

---

## disorimat_test01

```python
def disorimat_test01(umatsa,symops):
```

Test version 1 for disorientation matrix calculation.
    First test implementation of disorientation computation. Used for
    validating algorithms before final implementation.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Create test orientations
        >>> M1 = np_euler_matrix(0, 0, 0)
# ...
```

---

## disorimat_test02

```python
def disorimat_test02(umatsa,symops):
```

Test version 2 for disorientation matrix calculation.
    Experimental/testing version of disorientation calculation. The
    disorientation is the minimum misorientation considering crystal

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Test with known misorientations
        >>> M1 = np.eye(3)
# ...
```

---

## eu2quat

```python
def eu2quat(phi1, Phi, phi2):
```

Convert Bunge Euler angles to quaternion representation.
    Euler angles follow the ZXZ convention (Bunge notation).

**Example**:
```python
>>> import numpy as np
        >>> # 45-degree rotation around each axis
        >>> phi1 = np.radians(45)
        >>> Phi = np.radians(45)
# ...
```

---

## euler_angles_from_matrix

```python
def euler_angles_from_matrix(Rl,deg=False):
```

Extract Euler angles from a rotation matrix.
    Inverse operation of np_euler_matrix. Converts a rotation matrix back
    to Bunge Euler angles (ZXZ convention).

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Create rotation matrix from Euler angles
        >>> phi1, Phi, phi2 = np.radians([45, 60, 30])
# ...
```

---

## euler_angles_reduction

```python
def euler_angles_reduction(Phi1,PHI,Phi2):
```

Reduce Euler angles to fundamental zone using symmetry operations.
    Applies crystal symmetry operations to find equivalent orientation
    with Euler angles in the fundamental zone (asymmetric unit).

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Cubic symmetry (simplified - just identity for demo)
        >>> symops = [np.eye(3)]
# ...
```

---

## grid_s1

```python
def grid_s1(resol, grids=6):
```

Generate uniformly distributed points on S¹ (circle).
    Used for sampling the third Euler angle in Hopf coordinates.

**Example**:
```python
>>> points = grid_s1(resol=2, grids=6)
        >>> print("Number of points:", len(points))
        >>> # 2^2 * 6 = 24 points
        >>> print("First few points:", points[:5])
# ...
```

---

## hopf2quat

```python
def hopf2quat(Points):
```

Convert Hopf coordinates to quaternions.
    Hopf coordinates (θ, φ, ψ) parameterize the unit quaternion sphere S³.

**Example**:
```python
>>> import numpy as np
        >>> # Single point in Hopf coordinates
        >>> points = [(np.pi/2, 0, 0)]
        >>> quats = hopf2quat(points)
# ...
```

---

## mat2quat02

```python
def mat2quat02(matrix): #orientation matrix to quaternion
```

Alternative matrix to quaternion conversion (version 2).
    Different algorithm for converting rotation matrix to quaternion.
    May handle numerical edge cases differently than mat_to_quat.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Create rotation matrix
        >>> angle = np.pi/3  # 60 degrees
# ...
```

---

## mat_to_quat

```python
def mat_to_quat(R):
```

Convert a rotation matrix to a quaternion representation.
    Uses the Shepperd's method for numerical stability.

**Example**:
```python
>>> import numpy as np
        >>> # 90-degree rotation around z-axis
        >>> R = np.array([[0, -1, 0],
        ...               [1,  0, 0],
# ...
```

---

## misori_sym_deg_quats

```python
def misori_sym_deg_quats(q1, q2, sym_quats):
```

Calculate minimum misorientation angle between two quaternions considering
    crystal symmetry operations. This is a Numba-optimized version.

**Example**:
```python
>>> import numpy as np
        >>> # Define two orientations
        >>> q1 = np.array([1.0, 0.0, 0.0, 0.0])
        >>> q2 = np.array([0.707, 0.0, 0.0, 0.707])  # 90° around z
# ...
```

---

## misori_sym_deg_sample_to_crystal_fast

```python
def misori_sym_deg_sample_to_crystal_fast(M1, M2, symops):
```

Compute misorientation angles (deg) between orientations M1 and M2
    considering crystal symmetry operations using vectorized operations.
    This is the fast version optimized with Numba parallel processing.

**Example**:
```python
>>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R
        >>> 
        >>> # Generate some random orientations
# ...
```

---

## misorimat

```python
def misorimat(umatsa):
```

Calculate misorientation matrix between two orientations.
    Computes the relative rotation (misorientation) between two crystal
    orientations. The misorientation matrix represents the rotation needed

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Two grain orientations
        >>> euler1 = np.radians([30, 45, 60])
# ...
```

---

## misorimat_ini

```python
def misorimat_ini(umatsa):
```

Initialize misorientation matrix calculation (initialization version).
    Preliminary version of misorientation matrix computation. Sets up
    the calculation framework for determining misorientation between

**Example**:
```python
>>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R
        >>> 
        >>> # Two random orientations
# ...
```

---

## mk_pix2xy

```python
def mk_pix2xy():
```

Create lookup tables for HEALPix pixel index to (x,y) coordinate conversion.
    These tables are used by pix2ang_nest function.

**Example**:
```python
>>> pix2x, pix2y = mk_pix2xy()
        >>> print("Table size:", len(pix2x))  # 1024
        >>> print("First few X values:", pix2x[:5])
        >>> print("First few Y values:", pix2y[:5])
# ...
```

---

## np_euler_matrix

```python
def np_euler_matrix(ai, aj, ak):
```

Convert Euler angles to rotation matrix (ZXZ convention, Bunge notation).
    Single orientation version.

**Example**:
```python
>>> import numpy as np
        >>> # 90-degree rotation around Z axis
        >>> phi1 = np.pi/2
        >>> Phi = 0
# ...
```

---

## np_eulers_matrices

```python
def np_eulers_matrices(data, deg=False):
```

Convert multiple Euler angles to rotation matrices (vectorized).
    Processes an array of Euler angle triplets efficiently.

**Example**:
```python
>>> import numpy as np
        >>> # Multiple orientations in degrees
        >>> euler_angles = np.array([
        ...     [0, 0, 0],        # Identity
# ...
```

---

## np_g2quats

```python
def np_g2quats(umatsa):
```

Convert multiple rotation matrices to quaternions (vectorized).
    Handles arrays of rotation matrices efficiently.

**Example**:
```python
>>> import numpy as np
        >>> # Create multiple rotation matrices
        >>> N = 100
        >>> from scipy.spatial.transform import Rotation as R
# ...
```

---

## np_gmat2rodrigues

```python
def np_gmat2rodrigues(g):
```

Convert rotation matrix to Rodrigues-Frank vector representation.
    Rodrigues vector = rotation_axis * tan(angle/2)

**Example**:
```python
>>> import numpy as np
        >>> # 90-degree rotation around Z
        >>> g = np.array([[0, -1, 0],
        ...               [1,  0, 0],
# ...
```

---

## np_inverse_euler_matrix

```python
def np_inverse_euler_matrix(ai, aj, ak):
```

Convert Euler angles to inverse (transpose) rotation matrix.
    Equivalent to the transpose of the forward rotation matrix.

**Example**:
```python
>>> import numpy as np
        >>> # Forward rotation
        >>> ai, aj, ak = np.pi/4, np.pi/3, np.pi/6
        >>> g_forward = np_euler_matrix(ai, aj, ak)
# ...
```

---

## np_ol_R_g

```python
def np_ol_R_g (R):
```

Convert Rodrigues-Frank vector to rotation matrix (NumPy version).
    NumPy implementation of Rodrigues-Frank to rotation matrix conversion.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Small rotation
        >>> R_small = np.array([0.01, 0.02, 0.03])
# ...
```

---

## np_ol_R_g2

```python
def np_ol_R_g2(R,epsilon, delta):
```

Alternative Rodrigues-Frank to rotation matrix (NumPy version, method 2).
    NumPy implementation of alternative conversion method.

**Example**:
```python
>>> import numpy as np
        >>> R = np.array([0.1, 0.2, 0.3])
        >>> g = np_ol_R_g2(R)
    Notes:
# ...
```

---

## np_ol_R_q2

```python
def np_ol_R_q2(R):
```

Convert Rodrigues-Frank vector to quaternion (method 2).
    Transforms Rodrigues-Frank representation to unit quaternion [w,x,y,z].

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Small rotation
        >>> R = np.array([0.01, 0, 0])
# ...
```

---

## np_ol_g_R

```python
def np_ol_g_R(g):
```

Convert rotation matrix to Rodrigues-Frank vector (NumPy version).
    NumPy implementation returning Rodrigues-Frank vector as numpy array.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Identity rotation
        >>> g_id = np.eye(3)
# ...
```

---

## np_ol_g_R2

```python
def np_ol_g_R2(g,epsilon, delta):
```

Alternative Rodrigues-Frank conversion (NumPy version, method 2).
    NumPy implementation of alternative Rodrigues-Frank calculation.

**Example**:
```python
>>> import numpy as np
        >>> g = np.eye(3)
        >>> R = np_ol_g_R2(g)
    Notes:
# ...
```

---

## np_ol_g_q2

```python
def np_ol_g_q2(g):
```

Convert rotation matrix to quaternion (method 2).
    Extracts unit quaternion representation from rotation matrix.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # 90° around Z
        >>> g = np.array([[0, -1, 0],
# ...
```

---

## np_ol_g_rtheta_rad

```python
def np_ol_g_rtheta_rad(g):
```

Convert rotation matrix to axis-angle representation (NumPy optimized version).

**Example**:
```python
>>> import numpy as np
        >>> # 120-degree rotation around [1,1,1]
        >>> angle = np.radians(120)
        >>> axis = np.array([1, 1, 1]) / np.sqrt(3)
# ...
```

---

## np_ol_q_g

```python
def np_ol_q_g(q):
```

Convert quaternion to rotation matrix.
    Transforms unit quaternion [w,x,y,z] to 3×3 rotation matrix.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Identity quaternion
        >>> q_id = np.array([1, 0, 0, 0])
# ...
```

---

## np_ol_rtheta_g_rad

```python
def np_ol_rtheta_g_rad(r, theta):
```

Convert axis-angle representation to rotation matrix (NumPy version).
    Uses Rodrigues' rotation formula.

**Example**:
```python
>>> import numpy as np
        >>> # 45-degree rotation around [1,1,1] axis
        >>> axis = np.array([1, 1, 1]) / np.sqrt(3)  # Normalize
        >>> angle = np.radians(45)
# ...
```

---

## np_rodrigues2gmat

```python
def np_rodrigues2gmat(rodrigues):
```

Convert Rodrigues-Frank vector to rotation matrix.

**Example**:
```python
>>> import numpy as np
        >>> # Rodrigues vector for 90° around Z
        >>> rod = np.array([0, 0, 1])  # tan(45°) = 1
        >>> g = np_rodrigues2gmat(rod)
# ...
```

---

## nside2npix

```python
def nside2npix(nside):
```

Calculate the number of pixels in a HEALPix map.
    HEALPix = Hierarchical Equal Area isoLatitude Pixelization.

**Example**:
```python
>>> # Low resolution
        >>> npix_low = nside2npix(nside=1)
        >>> print("Pixels for nside=1:", npix_low)  # 12
        >>> # Medium resolution
# ...
```

---

## ol_R_g

```python
def ol_R_g (R):
```

Convert Rodrigues-Frank vector to rotation matrix (list version).
    Reconstructs rotation matrix from Rodrigues-Frank representation.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Rodrigues vector for 90° around Z
        >>> R = [0, 0, 1]  # tan(45°) = 1
# ...
```

---

## ol_R_g2

```python
def ol_R_g2(R):
```

Alternative Rodrigues-Frank to rotation matrix (list version, method 2).
    Second implementation of Rodrigues-Frank to rotation matrix conversion.

**Example**:
```python
>>> R = [0, 0, 0.5]
        >>> g = ol_R_g2(R)
    Notes:
        - Alternative implementation
# ...
```

---

## ol_g_R

```python
def ol_g_R(g):
```

Convert rotation matrix to Rodrigues-Frank vector (list version).
    Calculates Rodrigues-Frank vector R = r·tan(θ/2) from rotation matrix,
    where r is the rotation axis and θ is the rotation angle.

**Example**:
```python
>>> # 90° rotation around Z
        >>> g = [[0, -1, 0],
        ...      [1,  0, 0],
        ...      [0,  0, 1]]
# ...
```

---

## ol_g_R2

```python
def ol_g_R2(g):
```

Alternative Rodrigues-Frank conversion (list version, method 2).
    Second implementation of rotation matrix to Rodrigues-Frank vector.
    May use different numerical approach for stability.

**Example**:
```python
>>> g = [[0, -1, 0],
        ...      [1,  0, 0],
        ...      [0,  0, 1]]
        >>> R = ol_g_R2(g)
# ...
```

---

## ol_g_rtheta_rad

```python
def ol_g_rtheta_rad(g):
```

Convert rotation matrix to axis-angle representation (Rodrigues-Frank vector).
    Returns rotation axis and angle.

**Example**:
```python
>>> import numpy as np
        >>> # 90-degree rotation around Z
        >>> g = [[0, -1, 0],
        ...      [1,  0, 0],
# ...
```

---

## ol_rtheta_g_rad

```python
def ol_rtheta_g_rad(r, theta):
```

Convert axis-angle representation to rotation matrix using Rodrigues' formula.

**Example**:
```python
>>> import numpy as np
        >>> # 90-degree rotation around Z axis
        >>> axis = [0, 0, 1]
        >>> angle = np.pi/2
# ...
```

---

## orilistMult

```python
def orilistMult(Mats,Dr):
```

Multiply a list of orientation matrices by symmetry operations.
    Applies symmetry operations to orientation matrices to generate all
    symmetrically equivalent orientations.

**Example**:
```python
>>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R
        >>> 
        >>> # Create orientation matrices
# ...
```

---

## passive_rotation

```python
def passive_rotation(an, aboutaxis, deg=False):
```

Perform passive rotation (coordinate transformation).
    Rotates the coordinate system while vector stays fixed in space.
    v' = g^T · v (transpose of rotation matrix).

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Rotate coordinates 90° around Z
        >>> g_z90 = np.array([[0, -1, 0],
# ...
```

---

## pix2ang_nest

```python
def pix2ang_nest(nside, ipix, pix2x, pix2y):
```

Convert HEALPix pixel index to spherical coordinates (theta, phi).
    Uses NESTED indexing scheme.

**Example**:
```python
>>> # Setup lookup tables
        >>> pix2x, pix2y = mk_pix2xy()
        >>> 
        >>> # Convert pixel 0 at nside=4
# ...
```

---

## quat_conjugate

```python
def quat_conjugate(q):
```

Compute the conjugate of a quaternion (inverse for unit quaternions).

**Example**:
```python
>>> import numpy as np
        >>> q = np.array([0.707, 0.707, 0.0, 0.0])
        >>> q_conj = quat_conjugate(q)
        >>> print("Original:", q)
# ...
```

---

## quat_misori_deg

```python
def quat_misori_deg(q1, q2):
```

Calculate the misorientation angle between two quaternions in degrees.
    Returns the minimum rotation angle needed to go from q1 to q2.

**Example**:
```python
>>> import numpy as np
        >>> # Same orientation
        >>> q1 = np.array([1.0, 0.0, 0.0, 0.0])
        >>> angle = quat_misori_deg(q1, q1)
# ...
```

---

## quat_mult

```python
def quat_mult(q1, q2):
```

Multiply two quaternions (Hamilton product).
    Order matters: q1 * q2 ≠ q2 * q1 in general.

**Example**:
```python
>>> import numpy as np
        >>> # Two 90-degree rotations around z-axis = 180-degree rotation
        >>> q_90 = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> q_result = quat_mult(q_90, q_90)
# ...
```

---

## quat_multiply

```python
def quat_multiply(q1, q2):
```

Quaternion multiplication (same as quat_mult, alternative implementation).

**Example**:
```python
>>> import numpy as np
        >>> q1 = np.array([1.0, 0.0, 0.0, 0.0])  # Identity
        >>> q2 = np.array([0.707, 0.707, 0.0, 0.0])  # 90° around x
        >>> result = quat_multiply(q1, q2)
# ...
```

---

## quat_to_mat

```python
def quat_to_mat(q):
```

Convert a quaternion to a rotation matrix.

**Example**:
```python
>>> import numpy as np
        >>> # Quaternion for 90-degree rotation around z-axis
        >>> q = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> R = quat_to_mat(q)
# ...
```

---

## rotation_from_axis_angle

```python
def rotation_from_axis_angle(ax,an,deg=False):
```

Generate rotation matrix from axis-angle representation using Rodrigues formula.
    Creates 3×3 rotation matrix for rotation by 'angle' radians around 'axis'.
    Implements Rodrigues' rotation formula for arbitrary axis rotations.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # 90° rotation around z-axis
        >>> axis = np.array([0, 0, 1])
# ...
```

---

## simple_grid

```python
def simple_grid(resol):
```

Generate a uniform grid of quaternions covering orientation space (SO(3)).
    Uses HEALPix for S² sampling and uniform sampling for S¹, combined via Hopf fibration.

**Example**:
```python
>>> # Low resolution grid (fast)
        >>> quats_low = simple_grid(resol=2)
        >>> print("Low resolution quaternions:", len(quats_low))
        >>> # 12 * 4^2 * (2^2 * 6) = 12 * 16 * 24 = 4608
# ...
```

---

## symmetry_elements

```python
def symmetry_elements(lattice):
```

Generate symmetry operation matrices for crystal system.
    Returns list of 3×3 rotation matrices representing all symmetry
    operations for specified crystal system.

**Example**:
```python
>>> symops = symmetry_elements('cubic')
        >>> print(f"Cubic has {len(symops)} symmetry operations")
        >>> # 24 operations for cubic (point group m-3m)
        >>> 
# ...
```

---

## symmetry_reduced_oris

```python
def symmetry_reduced_oris(umatsa, symops):
```

Reduce orientations to fundamental zone using crystal symmetry.
    Applies symmetry operations to bring all orientations into the
    fundamental zone (asymmetric unit) of orientation space. This

**Example**:
```python
>>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R
        >>> 
        >>> # Generate random orientations
# ...
```

---

## symposMult

```python
def symposMult(sympos,Mats):
```

Multiply symmetry operations to generate composite symmetry operations.
    Computes the product of two sets of symmetry operations to generate
    all possible combinations.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Define some symmetry operations
        >>> # 90° rotation around Z
# ...
```

---

## symposMult02

```python
def symposMult02(sympos,Mats):
```

Alternative implementation of symmetry operations multiplication.
    Similar to symposMult but may use different algorithm or ordering.
    Generates all products of two symmetry operation sets.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Point group symmetry operations
        >>> identity = np.eye(3)
# ...
```

---

## trace_to_angle

```python
def trace_to_angle(tr, out="deg"):
```

Convert rotation matrix trace to rotation angle.
    Uses the trace (sum of diagonal elements) of a rotation matrix to
    calculate the rotation angle. Based on the relation:

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> # Create rotation matrix
        >>> theta_original = np.radians(45)
# ...
```

---


# Orilib.py Quick Reference Guide

## Overview

The **orilib.py** module provides comprehensive utilities for crystallographic orientation analysis including quaternion operations, Euler angle conversions, Rodrigues vectors, and orientation sampling methods.

## Quick Start

```python
import numpy as np
from orilib import *

# Convert rotation matrix to quaternion
R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
q = mat_to_quat(R)

# Calculate misorientation
angle = quat_misori_deg(q1, q2)

# Generate uniform orientation grid
quats = simple_grid(resol=3)
```

## Most Commonly Used Functions

### 1. Quaternion Operations

#### Basic Conversions
```python
# Rotation matrix to quaternion
R = np.eye(3)  # 3×3 rotation matrix
q = mat_to_quat(R)  # Returns [w, x, y, z]

# Quaternion to rotation matrix
q = np.array([1, 0, 0, 0])  # [w, x, y, z]
R = quat_to_mat(q)  # Returns 3×3 matrix
```

#### Quaternion Multiplication
```python
# Multiply two quaternions (Hamilton product)
q1 = np.array([1, 0, 0, 0])
q2 = np.array([0.707, 0, 0, 0.707])
q_result = quat_mult(q1, q2)

# Alternative function
q_result = quat_multiply(q1, q2)

# Quaternion conjugate (inverse for unit quaternions)
q_inv = quat_conjugate(q)
```

#### Misorientation Calculations
```python
# Misorientation angle between two quaternions (degrees)
angle = quat_misori_deg(q1, q2)

# With crystal symmetry
sym_quats = np.array([[1,0,0,0], ...])  # Symmetry operations
min_angle = misori_sym_deg_quats(q1, q2, sym_quats)

# Fast vectorized version for many orientations
M1 = np.random.randn(100, 3, 3)  # 100 orientation matrices
M2 = np.random.randn(100, 3, 3)
symops = np.array([np.eye(3)])  # Symmetry operations
misorientations = misori_sym_deg_sample_to_crystal_fast(M1, M2, symops)
```

### 2. Euler Angle Conversions

#### Euler to Quaternion
```python
# Bunge Euler angles (ZXZ convention) to quaternion
phi1 = np.radians(45)  # First rotation around Z
Phi = np.radians(30)   # Second rotation around X'
phi2 = np.radians(60)  # Third rotation around Z''
q = eu2quat(phi1, Phi, phi2)
```

#### Euler to Rotation Matrix
```python
# Single orientation
ai, aj, ak = np.pi/4, np.pi/3, np.pi/6  # Euler angles in radians
g = np_euler_matrix(ai, aj, ak)

# Multiple orientations (vectorized)
euler_angles = np.array([
    [0, 0, 0],
    [90, 0, 0],
    [0, 90, 0]
])
matrices = np_eulers_matrices(euler_angles, deg=True)
print("Shape:", matrices.shape)  # (3, 3, 3)

# Inverse rotation matrix
g_inv = np_inverse_euler_matrix(ai, aj, ak)
```

### 3. Rodrigues-Frank Vectors

#### Matrix to Axis-Angle
```python
# Rotation matrix to axis-angle representation
g = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
axis, angle = np_ol_g_rtheta_rad(g)
print("Rotation axis:", axis)
print("Rotation angle (rad):", angle)

# Convert to Rodrigues-Frank vector
rodrigues = np_gmat2rodrigues(g)
```

#### Axis-Angle to Matrix
```python
# Axis-angle to rotation matrix (Rodrigues' formula)
axis = np.array([0, 0, 1])  # Rotation axis
angle = np.pi/2  # 90 degrees
g = np_ol_rtheta_g_rad(axis, angle)

# Rodrigues vector to matrix
rodrigues = np.array([0, 0, 1])  # axis × tan(angle/2)
g = np_rodrigues2gmat(rodrigues)
```

### 4. Orientation Sampling

#### Simple Grid (HEALPix + Hopf)
```python
# Generate uniform orientation grid
resol = 3  # Resolution parameter
quats = simple_grid(resol)
# Returns list of quaternions [w, x, y, z]
# Number of quaternions ≈ 12 × 4^resol × (2^resol × 6)

print(f"Generated {len(quats)} orientations")

# Convert to rotation matrices
matrices = np.array([quat_to_mat(np.array(q)) for q in quats])
```

#### HEALPix Sampling
```python
# Calculate number of pixels
nside = 16  # Must be power of 2
npix = nside2npix(nside)
print(f"Pixels: {npix}")  # 12 × 16² = 3072

# Create lookup tables
pix2x, pix2y = mk_pix2xy()

# Convert pixel index to spherical coordinates
nside = 8
for ipix in range(nside2npix(nside)):
    theta, phi = pix2ang_nest(nside, ipix, pix2x, pix2y)
    # theta: colatitude [0, π]
    # phi: azimuth [0, 2π]
```

#### Hopf Coordinates
```python
# Convert Hopf coordinates to quaternions
hopf_points = [
    (np.pi/2, 0, 0),
    (np.pi/2, np.pi, 0),
    (np.pi/4, np.pi/4, np.pi/4)
]
quats = hopf2quat(hopf_points)

# Generate points on S¹ (circle)
psi_points = grid_s1(resol=3)
print(f"S¹ points: {len(psi_points)}")
```

### 5. Vectorized Operations

#### Batch Quaternion Operations
```python
# Convert multiple matrices to quaternions (vectorized)
matrices = np.random.randn(100, 3, 3)  # 100 rotation matrices
Q = np_g2quats(matrices)
print("Shape:", Q.shape)  # (4, 100)

# Quaternion logarithm
QM = Q.reshape(4, 10, 10)  # Reshape for Qlog
qlog = Qlog(QM)

# Quaternion products
P = np.random.randn(4, 100)
P = P / np.linalg.norm(P, axis=0)  # Normalize
Q = np.random.randn(4, 100)
Q = Q / np.linalg.norm(Q, axis=0)
result = Qproduct(P, Q)

# Single symmetry applied to multiple quaternions
sym = np.array([0.707, 0, 0, 0.707])  # 90° around Z
Q = np.random.randn(4, 50)
SQ = QMatproduct(sym, Q)
```

## Complete Workflow Examples

### Example 1: Misorientation Analysis

```python
import numpy as np
from orilib import mat_to_quat, quat_misori_deg

# Convert orientation matrices to quaternions
ori1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # Identity
ori2 = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) # 90° around Z

q1 = mat_to_quat(ori1)
q2 = mat_to_quat(ori2)

# Calculate misorientation
angle = quat_misori_deg(q1, q2)
print(f"Misorientation: {angle:.2f}°")  # 90.00°
```

### Example 2: Euler Angle Workflow

```python
import numpy as np
from orilib import eu2quat, quat_to_mat, np_euler_matrix

# Define Euler angles (Bunge convention)
phi1 = np.radians(30)
Phi = np.radians(45)
phi2 = np.radians(60)

# Method 1: Euler → Quaternion → Matrix
q = eu2quat(phi1, Phi, phi2)
R1 = quat_to_mat(q)

# Method 2: Euler → Matrix directly
R2 = np_euler_matrix(phi1, Phi, phi2)

# Verify they're the same
print("Close match:", np.allclose(R1, R2))

# Inverse matrix
R_inv = np_inverse_euler_matrix(phi1, Phi, phi2)
print("Identity:", np.allclose(R1 @ R_inv, np.eye(3)))
```

### Example 3: Texture Analysis

```python
import numpy as np
from orilib import simple_grid, quat_to_mat, quat_misori_deg

# Generate uniform orientation grid
quats = simple_grid(resol=2)
print(f"Total orientations: {len(quats)}")

# Convert to matrices
matrices = np.array([quat_to_mat(np.array(q)) for q in quats])

# Calculate misorientations from reference
ref_quat = np.array([1, 0, 0, 0])  # Identity
misorientations = []
for q in quats:
    angle = quat_misori_deg(ref_quat, np.array(q))
    misorientations.append(angle)

# Find orientations within 10° of reference
close_oris = [i for i, angle in enumerate(misorientations) if angle < 10]
print(f"Orientations within 10°: {len(close_oris)}")
```

### Example 4: Symmetry Operations

```python
import numpy as np
from orilib import quat_mult, quat_conjugate, misori_sym_deg_quats

# Define cubic symmetry operations (simplified example)
sym_quats = np.array([
    [1, 0, 0, 0],           # Identity
    [0.707, 0.707, 0, 0],   # 90° around X
    [0.707, 0, 0.707, 0],   # 90° around Y
    [0.707, 0, 0, 0.707]    # 90° around Z
])

# Two orientations
q1 = np.array([0.924, 0.383, 0, 0])
q2 = np.array([0.924, 0, 0.383, 0])

# Calculate minimum misorientation considering symmetry
min_angle = misori_sym_deg_quats(q1, q2, sym_quats)
print(f"Minimum misorientation: {min_angle:.2f}°")
```

### Example 5: Batch Processing

```python
import numpy as np
from orilib import np_eulers_matrices, np_g2quats, misori_sym_deg_sample_to_crystal_fast

# Generate random Euler angles
N = 1000
euler1 = np.random.rand(N, 3) * np.array([360, 180, 360])
euler2 = np.random.rand(N, 3) * np.array([360, 180, 360])

# Convert to matrices (vectorized)
M1 = np_eulers_matrices(euler1, deg=True)
M2 = np_eulers_matrices(euler2, deg=True)

# Calculate misorientations with symmetry (fast!)
symops = np.array([np.eye(3)])  # Just identity for example
misorientations = misori_sym_deg_sample_to_crystal_fast(M1, M2, symops)

print(f"Processed {N} orientation pairs")
print(f"Mean misorientation: {np.mean(misorientations):.2f}°")
print(f"Max misorientation: {np.max(misorientations):.2f}°")
```

## Function Parameters Quick Reference

### Quaternion Functions

#### `mat_to_quat(R)`
- **Input**: R (3×3 numpy array) - Rotation matrix
- **Output**: q (4,) numpy array - [w, x, y, z]

#### `quat_to_mat(q)`
- **Input**: q (4,) numpy array - [w, x, y, z]
- **Output**: R (3×3 numpy array) - Rotation matrix

#### `quat_mult(q1, q2)`
- **Input**: q1, q2 (4,) arrays - Quaternions
- **Output**: q (4,) array - Product q1 × q2

#### `quat_misori_deg(q1, q2)`
- **Input**: q1, q2 (4,) arrays - Quaternions
- **Output**: float - Angle in degrees [0, 180]

### Euler Functions

#### `eu2quat(phi1, Phi, phi2)`
- **Input**: phi1, Phi, phi2 (float) - Euler angles in radians (ZXZ)
- **Output**: q (4,) array - Quaternion [w, x, y, z]

#### `np_euler_matrix(ai, aj, ak)`
- **Input**: ai, aj, ak (float) - Euler angles in radians
- **Output**: g (3×3 array) - Rotation matrix

#### `np_eulers_matrices(data, deg=False)`
- **Input**: data (N×3 array) - Euler angles, deg (bool)
- **Output**: g (N×3×3 array) - Rotation matrices

### Rodrigues Functions

#### `np_gmat2rodrigues(g)`
- **Input**: g (3×3 array) - Rotation matrix
- **Output**: rodrigues (3,) array - Rodrigues-Frank vector

#### `np_rodrigues2gmat(rodrigues)`
- **Input**: rodrigues (3,) array - Rodrigues-Frank vector
- **Output**: g (3×3 array) - Rotation matrix

#### `np_ol_g_rtheta_rad(g)`
- **Input**: g (3×3 array) - Rotation matrix
- **Output**: r (3,) array - Rotation axis, theta (float) - Angle in radians

### Orientation Sampling

#### `simple_grid(resol)`
- **Input**: resol (int) - Resolution parameter
- **Output**: list of lists - Quaternions [[w,x,y,z], ...]

#### `nside2npix(nside)`
- **Input**: nside (int) - HEALPix resolution (power of 2)
- **Output**: int - Number of pixels (12 × nside²)

#### `pix2ang_nest(nside, ipix, pix2x, pix2y)`
- **Input**: nside (int), ipix (int), pix2x/pix2y (lookups)
- **Output**: theta, phi (float) - Angles in radians

## Conversion Reference

### Quaternion ↔ Matrix
```python
# Matrix → Quaternion
q = mat_to_quat(R)

# Quaternion → Matrix
R = quat_to_mat(q)
```

### Euler ↔ Quaternion ↔ Matrix
```python
# Euler → Quaternion
q = eu2quat(phi1, Phi, phi2)

# Quaternion → Matrix
R = quat_to_mat(q)

# Euler → Matrix (direct)
R = np_euler_matrix(phi1, Phi, phi2)
```

### Matrix ↔ Rodrigues ↔ Axis-Angle
```python
# Matrix → Rodrigues
rod = np_gmat2rodrigues(R)

# Rodrigues → Matrix
R = np_rodrigues2gmat(rod)

# Matrix → Axis-Angle
axis, angle = np_ol_g_rtheta_rad(R)

# Axis-Angle → Matrix
R = np_ol_rtheta_g_rad(axis, angle)
```

## Mathematical Formulas

### Quaternion Representation
```
Quaternion: q = [w, x, y, z] = w + xi + yj + zk
Unit constraint: w² + x² + y² + z² = 1
```

### Quaternion Multiplication
```
q₁ × q₂ = [w₁w₂ - x₁x₂ - y₁y₂ - z₁z₂,
           w₁x₂ + x₁w₂ + y₁z₂ - z₁y₂,
           w₁y₂ - x₁z₂ + y₁w₂ + z₁x₂,
           w₁z₂ + x₁y₂ - y₁x₂ + z₁w₂]
```

### Misorientation Angle
```
angle = 2 × arccos(|q₁·q₂|)  # in radians
      = 2 × arccos(|w₁w₂ + x₁x₂ + y₁y₂ + z₁z₂|)
```

### Euler Angles (Bunge Convention)
```
ZXZ convention: Rotate φ₁ around Z, then Φ around X', then φ₂ around Z''
Range: φ₁ ∈ [0, 2π], Φ ∈ [0, π], φ₂ ∈ [0, 2π]
```

### Rodrigues-Frank Vector
```
r = n × tan(θ/2)
where n is rotation axis (unit vector), θ is rotation angle
```

### HEALPix
```
Total pixels: Npix = 12 × Nside²
Nside must be power of 2
```

## Common Patterns

### Pattern 1: Load EBSD Data and Calculate Misorientations
```python
import numpy as np
from orilib import np_eulers_matrices, np_g2quats, quat_misori_deg

# Load Euler angles from EBSD data
euler_angles = np.loadtxt('ebsd_data.txt')  # N×3 array

# Convert to matrices
matrices = np_eulers_matrices(euler_angles, deg=True)

# Convert to quaternions
Q = np_g2quats(matrices)  # (4, N)

# Calculate misorientation distribution
ref_quat = Q[:, 0]  # First grain as reference
misorientations = []
for i in range(1, Q.shape[1]):
    angle = quat_misori_deg(ref_quat, Q[:, i])
    misorientations.append(angle)

# Analyze
print(f"Mean: {np.mean(misorientations):.2f}°")
print(f"Std: {np.std(misorientations):.2f}°")
```

### Pattern 2: Generate Uniform Orientation Distribution Function (ODF)
```python
from orilib import simple_grid, quat_to_mat

# Generate grid
quats = simple_grid(resol=3)

# Convert to your preferred representation
matrices = [quat_to_mat(np.array(q)) for q in quats]

# Use for texture analysis, pole figure generation, etc.
print(f"Generated {len(matrices)} orientations")
```

### Pattern 3: Rotation Composition
```python
from orilib import quat_mult, quat_to_mat
import numpy as np

# Two successive rotations
q1 = np.array([0.707, 0.707, 0, 0])  # 90° around X
q2 = np.array([0.707, 0, 0.707, 0])  # 90° around Y

# Compose rotations
q_total = quat_mult(q1, q2)

# Verify with matrices
R1 = quat_to_mat(q1)
R2 = quat_to_mat(q2)
R_total = R2 @ R1  # Note: order matters
R_from_quat = quat_to_mat(q_total)

print("Match:", np.allclose(R_total, R_from_quat))
```

## Performance Tips

### 1. Use Vectorized Functions
```python
# Slow: Loop over orientations
for ori in orientations:
    q = mat_to_quat(ori)

# Fast: Vectorized conversion
Q = np_g2quats(orientations)
```

### 2. Use Numba-Optimized Functions
```python
# Functions with @njit decorator are optimized:
# - mat_to_quat
# - quat_to_mat
# - quat_mult
# - quat_misori_deg
# - misori_sym_deg_quats
# - misori_sym_deg_sample_to_crystal_fast

# These are much faster for large datasets
```

### 3. Batch Process Misorientations
```python
# Use the fast version for many pairs
misorientations = misori_sym_deg_sample_to_crystal_fast(M1, M2, symops)

# Not: loop with quat_misori_deg (much slower)
```

### 4. Pre-compute Symmetry Operations
```python
# Compute symmetry once
from crystallography_functions import cubic_symmetry
symops = cubic_symmetry()  # 48 operations

# Reuse for all calculations
for grain_pair in grain_pairs:
    misori = misori_sym_deg_sample_to_crystal_fast(..., symops)
```

## Common Pitfalls

### 1. Quaternion Convention
```python
# Orilib uses [w, x, y, z] convention
q = np.array([1, 0, 0, 0])  # Identity
#             ↑  ↑  ↑  ↑
#             w  x  y  z

# Some libraries use [x, y, z, w] - be careful!
```

### 2. Euler Angle Units
```python
# Functions expect radians
phi1_deg = 45
phi1_rad = np.radians(phi1_deg)  # Convert!
q = eu2quat(phi1_rad, Phi_rad, phi2_rad)

# For np_eulers_matrices, use deg=True flag
matrices = np_eulers_matrices(euler_degrees, deg=True)
```

### 3. Rotation Order
```python
# Quaternion multiplication is not commutative
q1 * q2 ≠ q2 * q1

# Order matters!
q_total = quat_mult(q1, q2)  # First q1, then q2
```

### 4. Matrix Shape for Vectorized Functions
```python
# np_g2quats expects (N, 3, 3)
matrices = np.array([...])  # Shape must be (N, 3, 3)
Q = np_g2quats(matrices)    # Returns (4, N)

# Not (3, 3, N) or (3, N, 3)
```

### 5. Normalization
```python
# Quaternions should be normalized
q = np.array([0.5, 0.5, 0.5, 0.5])
q_normalized = q / np.linalg.norm(q)

# Functions like quat_to_mat expect unit quaternions
```

## Integration with Other Modules

### With Projlib
```python
from orilib import simple_grid, quat_to_mat
from projlib import equalarea_directions

# Generate orientations
quats = simple_grid(resol=2)
matrices = [quat_to_mat(np.array(q)) for q in quats]

# Project one direction through all orientations
direction = np.array([0, 0, 1])
projected = np.array([M @ direction for M in matrices]).T

# Get pole figure coordinates
proj_coords = equalarea_directions(projected)
```

### With Plotlib
```python
from orilib import mat_to_quat, quat_misori_deg
from plotlib import plotter, get_colors
import matplotlib.pyplot as plt

# Calculate misorientations
misorientations = [...]  # Your calculations

# Map to colors
cmap = plt.cm.viridis
colors = get_colors(np.array(misorientations), cmap, vmin=0, vmax=90)

# Plot with plotter
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')
# Add your colored scatter plot
```

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "Invalid rotation matrix" | Check det(R)=1 and R.T @ R = I |
| Quaternion not normalized | Normalize: q = q / np.linalg.norm(q) |
| Wrong Euler angles | Check convention (ZXZ Bunge) and units (radians) |
| Slow misorientation calc | Use misori_sym_deg_sample_to_crystal_fast |
| Array shape errors | Check (N,3,3) for matrices, (4,N) for quaternions |
| Gimbal lock issues | Use quaternions instead of Euler angles |

## For More Information

- **Full documentation**: See `orilib_commented.py` for complete docstrings
- **Examples**: Check each function's usage examples in docstrings
- **Mathematical background**: See orientation representation literature
- **Performance**: Use Numba-optimized functions for large datasets

---

**Pro tip**: For EBSD/texture analysis, start with `simple_grid()` for uniform ODF, then use vectorized functions for fast processing of thousands of orientations.

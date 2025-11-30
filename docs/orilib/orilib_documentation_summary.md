# Orilib.py Documentation Summary

## Overview

The **orilib.py** module is a comprehensive orientation library providing utilities for crystallographic orientation analysis, including quaternion operations, Euler angle conversions, Rodrigues-Frank vectors, axis-angle representations, and uniform orientation sampling using HEALPix and Hopf coordinates.

## Documentation Structure

### Files
- **orilib.py**: Original source code
- **orilib_commented.py**: Fully documented version with comprehensive docstrings
- **orilib_quick_reference.md**: Quick reference guide (companion file)
- **orilib_documentation_summary.md**: This document

## Module Categories

### 1. Quaternion Utilities
Core functions for quaternion operations with Numba optimization for performance.

### 2. Euler Angles and Rotation Matrices
Conversions between Euler angles and rotation matrices using Bunge convention (ZXZ).

### 3. Rodrigues-Frank and Axis-Angle Representations
Alternative rotation representations for crystallographic analysis.

### 4. Vectorized Quaternion Operations
Batch processing functions for arrays of orientations.

### 5. Orientation Sampling and Gridding
HEALPix-based uniform sampling of orientation space (SO(3)).

## Complete Function List

### Quaternion Operations

#### `mat_to_quat(R)`
**Purpose**: Convert rotation matrix to quaternion using Shepperd's method

**Input**:
- `R`: numpy array (3×3) - Rotation matrix (orthogonal with det=1)

**Output**:
- `q`: numpy array (4,) - Unit quaternion [w, x, y, z]

**Features**:
- Numba-optimized (@njit)
- Numerically stable
- Handles all cases (trace > 0, trace < 0, etc.)

**Usage**:
```python
R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
q = mat_to_quat(R)
# Returns approximately [0.707, 0, 0, 0.707]
```

#### `quat_to_mat(q)`
**Purpose**: Convert quaternion to rotation matrix

**Input**:
- `q`: numpy array (4,) - Unit quaternion [w, x, y, z]

**Output**:
- `R`: numpy array (3×3) - Rotation matrix

**Features**:
- Numba-optimized (@njit)
- Fast conversion for rendering, transformations

**Usage**:
```python
q = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
R = quat_to_mat(q)
# 90-degree rotation around Z-axis
```

#### `quat_mult(q1, q2)` and `quat_multiply(q1, q2)`
**Purpose**: Multiply two quaternions (Hamilton product)

**Input**:
- `q1`, `q2`: numpy arrays (4,) - Quaternions [w, x, y, z]

**Output**:
- `q`: numpy array (4,) - Product quaternion q1 × q2

**Note**: Order matters (non-commutative)

**Features**:
- Both functions are Numba-optimized
- Used for rotation composition

**Usage**:
```python
q_90 = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
q_result = quat_mult(q_90, q_90)
# Two 90° rotations = 180° rotation
```

#### `quat_conjugate(q)`
**Purpose**: Compute quaternion conjugate (inverse for unit quaternions)

**Input**:
- `q`: numpy array (4,) - Quaternion [w, x, y, z]

**Output**:
- `q_conj`: numpy array (4,) - Conjugate [w, -x, -y, -z]

**Features**:
- Numba-optimized
- For unit quaternions: q × q* = [1, 0, 0, 0]

**Usage**:
```python
q = np.array([0.707, 0.707, 0, 0])
q_inv = quat_conjugate(q)
result = quat_mult(q, q_inv)  # Identity
```

#### `quat_misori_deg(q1, q2)`
**Purpose**: Calculate misorientation angle between two quaternions

**Input**:
- `q1`, `q2`: numpy arrays (4,) - Quaternions

**Output**:
- `angle`: float - Misorientation angle in degrees [0, 180]

**Formula**: angle = 2 × arccos(|q₁·q₂|)

**Features**:
- Numba-optimized
- Handles numerical precision issues
- Always returns acute angle

**Usage**:
```python
q1 = np.array([1, 0, 0, 0])  # Identity
q2 = np.array([0.707, 0, 0, 0.707])  # 90° rotation
angle = quat_misori_deg(q1, q2)  # 90.0 degrees
```

#### `misori_sym_deg_quats(q1, q2, sym_quats)`
**Purpose**: Calculate minimum misorientation considering crystal symmetry

**Input**:
- `q1`, `q2`: numpy arrays (4,) - Quaternions
- `sym_quats`: numpy array (Ns, 4) - Symmetry operation quaternions

**Output**:
- `min_angle`: float - Minimum misorientation in degrees

**Features**:
- Numba-optimized
- Tests all symmetry variants
- Essential for crystallographic analysis

**Usage**:
```python
# Cubic symmetry (simplified)
sym_quats = np.array([[1,0,0,0], [0.707,0.707,0,0], ...])
min_angle = misori_sym_deg_quats(q1, q2, sym_quats)
```

#### `misori_sym_deg_sample_to_crystal_fast(M1, M2, symops)`
**Purpose**: Fast vectorized misorientation calculation with symmetry

**Input**:
- `M1`, `M2`: numpy arrays (N, 3, 3) - Orientation matrices
- `symops`: numpy array (Ns, 3, 3) - Symmetry operation matrices

**Output**:
- `miso`: numpy array (N,) - Misorientation angles in degrees

**Features**:
- Numba parallel processing (@njit(parallel=True))
- Extremely fast for large datasets
- Uses scipy.spatial.transform internally

**Usage**:
```python
M1 = R.random(100).as_matrix()
M2 = R.random(100).as_matrix()
symops = cubic_symmetry_matrices()
misorientations = misori_sym_deg_sample_to_crystal_fast(M1, M2, symops)
# Processes 100 pairs very quickly
```

### Euler Angle Conversions

#### `eu2quat(phi1, Phi, phi2)`
**Purpose**: Convert Bunge Euler angles to quaternion

**Input**:
- `phi1`: float - First rotation around Z (radians)
- `Phi`: float - Second rotation around X' (radians)
- `phi2`: float - Third rotation around Z'' (radians)

**Output**:
- `q`: numpy array (4,) - Quaternion with positive w

**Convention**: ZXZ Bunge convention
- Range: φ₁ ∈ [0, 2π], Φ ∈ [0, π], φ₂ ∈ [0, 2π]

**Usage**:
```python
phi1 = np.radians(45)
Phi = np.radians(30)
phi2 = np.radians(60)
q = eu2quat(phi1, Phi, phi2)
```

#### `np_euler_matrix(ai, aj, ak)`
**Purpose**: Convert Euler angles to rotation matrix (single orientation)

**Input**:
- `ai`, `aj`, `ak`: floats - Euler angles in radians (ZXZ)

**Output**:
- `g`: numpy array (3×3) - Rotation matrix

**Usage**:
```python
g = np_euler_matrix(np.pi/4, np.pi/3, np.pi/6)
```

#### `np_eulers_matrices(data, deg=False)`
**Purpose**: Convert multiple Euler angles to matrices (vectorized)

**Input**:
- `data`: numpy array (N×3) - Euler angles [φ₁, Φ, φ₂]
- `deg`: bool - If True, input is in degrees (default: False)

**Output**:
- `g`: numpy array (N×3×3) - Rotation matrices

**Features**:
- Fully vectorized (no loops)
- Efficient for large EBSD datasets
- Direct array operations

**Usage**:
```python
euler_angles = np.array([[0,0,0], [90,0,0], [0,90,0]])
matrices = np_eulers_matrices(euler_angles, deg=True)
# Shape: (3, 3, 3)
```

#### `np_inverse_euler_matrix(ai, aj, ak)`
**Purpose**: Convert Euler angles to inverse rotation matrix

**Input**:
- `ai`, `aj`, `ak`: floats - Euler angles in radians

**Output**:
- `U`: numpy array (3×3) - Inverse (transpose) rotation matrix

**Note**: For rotation matrices, inverse = transpose

**Usage**:
```python
g_inv = np_inverse_euler_matrix(ai, aj, ak)
# Verify: g @ g_inv = I
```

### Rodrigues-Frank and Axis-Angle

#### `ol_g_rtheta_rad(g)` and `np_ol_g_rtheta_rad(g)`
**Purpose**: Convert rotation matrix to axis-angle representation

**Input**:
- `g`: rotation matrix (list or numpy array, 3×3)

**Output**:
- `r`: rotation axis (unit vector, 3 elements)
- `ptheta`: rotation angle in radians

**Features**:
- Handles special cases (small angles, 180° rotations)
- `np_ol_g_rtheta_rad` is NumPy-optimized version

**Usage**:
```python
g = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
axis, angle = np_ol_g_rtheta_rad(g)
# axis ≈ [0, 0, 1], angle ≈ π/2
```

#### `ol_rtheta_g_rad(r, theta)` and `np_ol_rtheta_g_rad(r, theta)`
**Purpose**: Convert axis-angle to rotation matrix (Rodrigues' formula)

**Input**:
- `r`: rotation axis (3 elements, should be unit vector)
- `theta`: rotation angle in radians

**Output**:
- `g`: rotation matrix (3×3)

**Formula**: Rodrigues' rotation formula
```
R = I + sin(θ)K + (1-cos(θ))K²
where K is the skew-symmetric matrix of the axis
```

**Usage**:
```python
axis = np.array([0, 0, 1])
angle = np.pi/2
g = np_ol_rtheta_g_rad(axis, angle)
```

#### `np_gmat2rodrigues(g)`
**Purpose**: Convert rotation matrix to Rodrigues-Frank vector

**Input**:
- `g`: numpy array (3×3) - Rotation matrix

**Output**:
- `rodrigues`: numpy array (3,) - Rodrigues vector (axis × tan(θ/2))

**Usage**:
```python
g = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
rod = np_gmat2rodrigues(g)
# For 90° rotation: tan(45°) = 1, so rod ≈ [0, 0, 1]
```

#### `np_rodrigues2gmat(rodrigues)`
**Purpose**: Convert Rodrigues-Frank vector to rotation matrix

**Input**:
- `rodrigues`: numpy array (3,) - Rodrigues vector

**Output**:
- `g`: numpy array (3×3) - Rotation matrix

**Usage**:
```python
rod = np.array([0, 0, 1])
g = np_rodrigues2gmat(rod)
```

### Vectorized Operations

#### `np_g2quats(umatsa)`
**Purpose**: Convert multiple rotation matrices to quaternions (vectorized)

**Input**:
- `umatsa`: numpy array (N×3×3) - Array of rotation matrices

**Output**:
- `Q`: numpy array (4×N) - Array of quaternions

**Features**:
- Fully vectorized (no loops)
- Handles all cases using conditional indexing
- Much faster than looping with mat_to_quat

**Usage**:
```python
matrices = np.random.randn(100, 3, 3)
Q = np_g2quats(matrices)
# Shape: (4, 100)
```

#### `Qlog(QM)`
**Purpose**: Compute logarithm of quaternions (logarithmic map)

**Input**:
- `QM`: numpy array (4, N, M) - Array of quaternions

**Output**:
- `qlog`: numpy array (4, N, M) - Logarithm of quaternions
  - qlog[0] = 0
  - qlog[1:4] = normalized vector part × angle

**Usage**:
```python
q = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
QM = q.reshape(4, 1, 1)
qlog = Qlog(QM)
```

#### `Qproduct(P, Q)`
**Purpose**: Quaternion product for arrays of quaternions

**Input**:
- `P`, `Q`: numpy arrays (4×N) - Arrays of quaternions

**Output**:
- `result`: numpy array (4×N) - Product quaternions P × Q

**Features**:
- Vectorized quaternion multiplication
- Processes N quaternion pairs simultaneously

**Usage**:
```python
P = np.random.randn(4, 100)
P = P / np.linalg.norm(P, axis=0)
Q = np.random.randn(4, 100)
Q = Q / np.linalg.norm(Q, axis=0)
result = Qproduct(P, Q)
```

#### `QMatproduct(sym, Q)`
**Purpose**: Multiply single symmetry quaternion with multiple quaternions

**Input**:
- `sym`: numpy array (4,) - Single symmetry quaternion
- `Q`: numpy array (4×N) - Array of quaternions

**Output**:
- `SQ`: numpy array (4×N) - Transformed quaternions

**Usage**:
```python
sym = np.array([0.707, 0, 0, 0.707])  # 90° around Z
Q = np.random.randn(4, 50)
SQ = QMatproduct(sym, Q)
```

### Orientation Sampling

#### `simple_grid(resol)`
**Purpose**: Generate uniform orientation grid covering SO(3)

**Input**:
- `resol`: int - Resolution parameter

**Output**:
- `quats`: list of lists - Uniformly distributed quaternions [[w,x,y,z], ...]

**Grid Size**:
- S¹ points: 2^resol × 6
- S² HEALPix nside: 2^resol
- Total quaternions: 12 × 4^resol × (2^resol × 6)

**Examples**:
- resol=1: ~576 quaternions
- resol=2: ~4,608 quaternions
- resol=3: ~36,864 quaternions
- resol=4: ~294,912 quaternions

**Features**:
- Uses HEALPix for S² (sphere) sampling
- Uniform grid on S¹ (circle) for third angle
- Combined via Hopf fibration for SO(3)

**Usage**:
```python
quats = simple_grid(resol=3)
print(f"Generated {len(quats)} orientations")
```

#### `nside2npix(nside)`
**Purpose**: Calculate number of pixels in HEALPix map

**Input**:
- `nside`: int - HEALPix resolution (must be power of 2)

**Output**:
- `npix`: int - Total pixels = 12 × nside²

**Usage**:
```python
npix = nside2npix(16)  # Returns 3072
```

#### `mk_pix2xy()`
**Purpose**: Create lookup tables for HEALPix indexing

**Output**:
- `pix2x`, `pix2y`: lists (1024 elements each) - Lookup tables

**Usage**:
```python
pix2x, pix2y = mk_pix2xy()
# Used by pix2ang_nest
```

#### `pix2ang_nest(nside, ipix, pix2x, pix2y)`
**Purpose**: Convert HEALPix pixel index to spherical coordinates

**Input**:
- `nside`: int - HEALPix resolution
- `ipix`: int - Pixel index [0, 12×nside²-1]
- `pix2x`, `pix2y`: lists - Lookup tables from mk_pix2xy()

**Output**:
- `theta`: float - Colatitude [0, π] radians
- `phi`: float - Azimuth [0, 2π] radians

**Usage**:
```python
pix2x, pix2y = mk_pix2xy()
theta, phi = pix2ang_nest(8, 100, pix2x, pix2y)
```

#### `grid_s1(resol, grids=6)`
**Purpose**: Generate uniformly distributed points on S¹ (circle)

**Input**:
- `resol`: int - Resolution parameter
- `grids`: int - Grid multiplier (default: 6)

**Output**:
- `points`: list - Angles in radians [0, 2π]

**Number of points**: 2^resol × grids

**Usage**:
```python
psi_points = grid_s1(resol=3)  # 48 points
```

#### `hopf2quat(Points)`
**Purpose**: Convert Hopf coordinates to quaternions

**Input**:
- `Points`: list of tuples - [(θ, φ, ψ), ...] in radians
  - θ ∈ [0, π]: polar angle on S²
  - φ ∈ [0, 2π]: azimuthal angle on S²
  - ψ ∈ [0, 2π]: angle on S¹

**Output**:
- `quats`: list of lists - Quaternions [[w,x,y,z], ...]

**Hopf Fibration**: Maps S³ → S² with fiber S¹
- Used to parameterize rotation space uniformly

**Usage**:
```python
points = [(np.pi/2, 0, 0), (np.pi/2, np.pi, 0)]
quats = hopf2quat(points)
```

## Mathematical Background

### Quaternion Mathematics

**Representation**:
```
q = w + xi + yj + zk = [w, x, y, z]
Unit constraint: w² + x² + y² + z² = 1
```

**Multiplication** (Hamilton product):
```
q₁ × q₂ = [w₁w₂ - v₁·v₂, w₁v₂ + w₂v₁ + v₁×v₂]
where v = [x, y, z] is vector part
```

**Rotation**:
- Rotate vector v by quaternion q: v' = q × [0,v] × q*
- Composition: q_total = q₂ × q₁ (first q₁, then q₂)

**Misorientation**:
```
Δq = q₂ × q₁*
angle = 2 × arccos(|w|)
```

### Euler Angles (Bunge Convention)

**Convention**: ZXZ (Roe convention, passive rotations)
- φ₁: Rotation around Z (sample axis)
- Φ: Rotation around X' (rotated axis)
- φ₂: Rotation around Z'' (crystal axis)

**Range**:
- φ₁ ∈ [0, 2π] or [0°, 360°]
- Φ ∈ [0, π] or [0°, 180°]
- φ₂ ∈ [0, 2π] or [0°, 360°]

**Matrix Form**:
```
R = R_Z(φ₁) × R_X(Φ) × R_Z(φ₂)
```

### Rodrigues-Frank Vector

**Definition**:
```
r = n × tan(θ/2)
```
where:
- n: unit rotation axis
- θ: rotation angle

**Properties**:
- Compact representation (3 values)
- No gimbal lock
- Singularity at θ = 180°

### Axis-Angle Representation

**Rodrigues' Formula**:
```
R = I + sin(θ)K + (1-cos(θ))K²
```
where K is skew-symmetric matrix of axis n:
```
K = [  0  -n₃  n₂]
    [ n₃   0  -n₁]
    [-n₂  n₁   0 ]
```

### HEALPix

**Hierarchical Equal Area isoLatitude Pixelization**:
- Divides sphere into equal-area pixels
- Hierarchical structure (nside = 2^k)
- Total pixels: N = 12 × nside²

**Properties**:
- Equal area (important for uniform sampling)
- Isolatitude rings (fast spherical harmonics)
- NESTED and RING indexing schemes

### Hopf Fibration

**Maps S³ → S²**:
- Base space: S² (2-sphere)
- Fiber: S¹ (circle)
- Total space: S³ (3-sphere) ≅ Unit quaternions

**Parameterization**:
```
q = [cos(θ/2)cos(ψ/2), sin(θ/2)cos(φ+ψ/2), 
     sin(θ/2)sin(φ+ψ/2), cos(θ/2)sin(ψ/2)]
```

## Integration Examples

### With Projlib
```python
from orilib import simple_grid, quat_to_mat
from projlib import equalarea_directions

# Generate orientations
quats = simple_grid(resol=2)
matrices = [quat_to_mat(np.array(q)) for q in quats]

# Project through orientations
direction = np.array([0, 0, 1])
projected = np.array([M @ direction for M in matrices]).T
coords = equalarea_directions(projected)
```

### With Plotlib
```python
from orilib import quat_misori_deg
from plotlib import plotter, get_colors
import matplotlib.pyplot as plt

# Calculate misorientations
misorientations = [...]

# Visualize
cmap = plt.cm.viridis
colors = get_colors(np.array(misorientations), cmap)
p = plotter()
p.plotProj(ProjType='equalarea')
```

### Complete Texture Analysis Pipeline
```python
import numpy as np
from orilib import np_eulers_matrices, np_g2quats, quat_misori_deg
from projlib import equalarea_directions, genprojgrid
from plotlib import plotter

# Load EBSD data
euler_angles = np.loadtxt('ebsd.txt')  # N×3 array

# Convert to matrices
matrices = np_eulers_matrices(euler_angles, deg=True)

# Project directions
direction = np.array([0, 0, 1])
projected = np.array([M @ direction for M in matrices]).T

# Create pole figure
proj_coords = equalarea_directions(projected)
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')
p.ax.scatter(proj_coords[0], proj_coords[1], s=1)
```

## Performance Considerations

### Numba Optimization
Functions with `@njit` decorator:
- `mat_to_quat`
- `quat_to_mat`
- `quat_mult`, `quat_multiply`
- `quat_conjugate`
- `quat_misori_deg`
- `misori_sym_deg_quats`
- `misori_sym_deg_sample_to_crystal_fast` (parallel)

**Benefits**:
- 10-100× speedup for large datasets
- Parallel processing for vectorized functions
- Near C-level performance

### Vectorization
Use vectorized functions instead of loops:

**Slow**:
```python
for ori in orientations:
    q = mat_to_quat(ori)
```

**Fast**:
```python
Q = np_g2quats(orientations)
```

### Memory Considerations
- Quaternions: 4 × N × 8 bytes = 32N bytes
- Matrices: 9 × N × 8 bytes = 72N bytes
- For N=100,000: ~7 MB vs ~3 MB

## Best Practices

### 1. Always Normalize Quaternions
```python
q = q / np.linalg.norm(q)
```

### 2. Use Correct Convention
- [w, x, y, z] in this library
- Some libraries use [x, y, z, w]

### 3. Check Rotation Matrix Validity
```python
assert np.allclose(np.linalg.det(R), 1)
assert np.allclose(R.T @ R, np.eye(3))
```

### 4. Use Radians for Angles
```python
# Convert degrees to radians
angle_rad = np.radians(angle_deg)
```

### 5. Batch Process When Possible
```python
# Use vectorized functions for large datasets
misorientations = misori_sym_deg_sample_to_crystal_fast(M1, M2, symops)
```

## Common Applications

### 1. EBSD Data Analysis
- Convert Euler angles to quaternions
- Calculate grain boundary misorientations
- Generate orientation distribution functions

### 2. Texture Analysis
- Generate uniform ODF grids
- Calculate pole figures
- Inverse pole figure coloring

### 3. Crystal Plasticity
- Track orientation evolution
- Calculate Schmid factors
- Model texture development

### 4. Grain Boundary Analysis
- Classify boundaries (low/high angle)
- Special boundaries (Σ3, Σ5, etc.)
- Misorientation distributions

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| Non-unit quaternion | Numerical errors | Normalize: q/‖q‖ |
| Invalid rotation matrix | Wrong construction | Check det=1, orthogonal |
| Wrong misorientation | Missing symmetry | Use misori_sym_deg_quats |
| Slow performance | Using loops | Use vectorized functions |
| Array shape error | Wrong dimension | Check (N,3,3) vs (3,3,N) |
| Gimbal lock | Euler angles | Switch to quaternions |

## Dependencies

**Required**:
- numpy
- numba (for optimization)
- scipy (for fast misorientation)

**Optional**:
- matplotlib (for visualization)
- projlib (for projections)
- plotlib (for pole figures)

## Version Information

**Python**: 2.7+ / 3.x compatible
**Numba**: Provides 10-100× speedup
**Documentation**: Complete with examples

## References

### Quaternions
- Shepperd, S. W. (1978). Quaternion from rotation matrix. *Journal of Guidance and Control*, 1(3), 223-224.

### Euler Angles
- Bunge, H. J. (2013). *Texture analysis in materials science*. Elsevier.

### HEALPix
- Górski, K. M., et al. (2005). HEALPix: A framework for high-resolution discretization and fast analysis of data distributed on the sphere. *The Astrophysical Journal*, 622(2), 759.

### Orientation Representations
- Morawiec, A. (2004). *Orientations and rotations: computations in crystallographic textures*. Springer.

## Summary

The orilib module provides:
- ✅ Complete quaternion operation suite
- ✅ Euler angle conversions (Bunge convention)
- ✅ Rodrigues-Frank vectors
- ✅ Axis-angle representations
- ✅ Numba-optimized performance
- ✅ Vectorized batch processing
- ✅ Uniform orientation sampling (HEALPix)
- ✅ Crystal symmetry support
- ✅ Comprehensive documentation

**Perfect for**: EBSD analysis, texture analysis, crystal plasticity, grain boundary characterization, orientation space sampling.

## For More Information

- **Quick start**: See orilib_quick_reference.md
- **Detailed examples**: See orilib_commented.py
- **Mathematical background**: See references above
- **Integration**: See projlib and plotlib documentation

---

**Module size**: ~43 KB, 40+ functions, all with comprehensive docstrings and usage examples.

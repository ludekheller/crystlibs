# Crystlib.py Documentation Summary

## Overview

The **crystlib.py** module provides utility functions for crystallographic calculations including Miller indices generation with crystal symmetry operations, lattice vector calculations for various crystal systems (cubic, tetragonal, hexagonal, monoclinic, triclinic), and reciprocal lattice calculations.

## Documentation Structure

### Files
- **crystlib.py**: Original source code
- **crystlib_commented.py**: Fully documented version with comprehensive docstrings
- **crystlib_quick_reference.md**: Quick reference guide (companion file)
- **crystlib_documentation_summary.md**: This document

## Module Purpose

This module serves as a fundamental building block for crystallographic analysis by providing:
1. Generation of unique Miller indices considering crystal symmetry
2. Calculation of real-space lattice vectors for any crystal system
3. Computation of reciprocal lattice vectors
4. Determination of d-spacings for diffraction analysis
5. Family grouping of crystallographically equivalent planes/directions

## Complete Function Documentation

### 1. Miller Indices Generation

#### `generate_hkls(hklmax, syms, hkls=[])`
**Purpose**: Generate unique Miller indices (hkl) considering crystal symmetry operations with rounding to avoid floating point issues

**Input**:
- `hklmax`: int - Maximum Miller index value (generates from -hklmax to +hklmax)
- `syms`: list of numpy arrays - List of 3×3 symmetry operation matrices
- `hkls`: list - Optional custom list of Miller indices to consider (default: empty list generates all)

**Output**:
- `hkls`: list of tuples - Unique Miller indices [(h₁,k₁,l₁), (h₂,k₂,l₂), ...]
- `hkls2`: dict - Dictionary mapping each unique hkl to its symmetry equivalents
  - Keys: unique hkl tuples
  - Values: lists of equivalent hkl tuples
- `fam`: dict - Dictionary of unique families with their multiplicities
  - Keys: representative hkl for each family
  - Values: number of distinct hkl in that family

**Algorithm**:
1. Generate all combinations from -hklmax to +hklmax (excluding (0,0,0))
2. Apply all symmetry operations to each hkl
3. Round values to avoid floating point comparison issues
4. Identify unique hkl by checking symmetry equivalents
5. Group into families based on permutation symmetry

**Features**:
- Handles any crystal system via symmetry operations
- Avoids duplicates through symmetry checking
- Rounds values for robust comparison
- Groups into crystallographic families

**Usage Example**:
```python
import numpy as np

# Cubic symmetry (identity + 3 90° rotations)
identity = np.eye(3)
rot_x = np.array([[1,0,0],[0,0,-1],[0,1,0]])
rot_y = np.array([[0,0,1],[0,1,0],[-1,0,0]])
rot_z = np.array([[0,-1,0],[1,0,0],[0,0,1]])
syms = [identity, rot_x, rot_y, rot_z]

# Generate hkls up to (2,2,2)
hkls, hkls2, families = generate_hkls(2, syms)

# Results
print(f"Unique planes: {len(hkls)}")
# Symmetry equivalents of (100)
print(hkls2.get((1.0, 0.0, 0.0)))
# Family information
print(families)
```

**Applications**:
- X-ray/neutron diffraction analysis
- EBSD pattern indexing
- Texture analysis
- Structure factor calculations

#### `generate_hkls01(hklmax, syms, hkls=[])`
**Purpose**: Alternative version without rounding for exact floating point comparison

**Differences from generate_hkls()**:
- No rounding of values
- Uses exact tuple comparison
- May have floating point precision issues

**When to use**:
- When exact precision is required
- When symmetry operations produce exact integer results
- For testing/verification

**Same Input/Output format as generate_hkls()**

#### `generate_hkls02(hklmax, syms, G, hkls=[])`
**Purpose**: Generate Miller indices with metric tensor transformation for non-orthogonal crystal systems

**Additional Input**:
- `G`: numpy array (3×3) - Metric tensor matrix
  - For direct space: G_ij = a_i · a_j
  - Accounts for non-orthogonal angles and different lattice parameters

**Formula**:
```
Symmetry operation in reciprocal space:
h' = G⁻¹ · S · G · h
where S is the symmetry operation matrix
```

**When to use**:
- Hexagonal/trigonal systems
- Monoclinic/triclinic systems
- Any non-orthogonal crystal system
- When proper reciprocal space symmetry is needed

**Usage Example**:
```python
# Hexagonal metric tensor
a = 3.0
c = 5.0
G = np.array([[a**2, -a**2/2, 0],
              [-a**2/2, a**2, 0],
              [0, 0, c**2]])

# Hexagonal symmetry
rot_60 = np.array([[0.5, -np.sqrt(3)/2, 0],
                   [np.sqrt(3)/2, 0.5, 0],
                   [0, 0, 1]])
syms_hex = [np.eye(3), rot_60]

hkls, hkls2, fam = generate_hkls02(2, syms_hex, G)
```

### 2. Family Grouping

#### `get_unique_families(hkls)`
**Purpose**: Group Miller indices into families based on permutation symmetry

**Input**:
- `hkls`: list or tuple - List of Miller index tuples [(h₁,k₁,l₁), ...]

**Output**:
- `dict`: Dictionary mapping representative hkl to multiplicity
  - Keys: Representative hkl tuple for each family
  - Values: Number of members in that family

**Family Definition**:
Two Miller indices belong to the same family if they are permutations of each other (considering absolute values).

**Examples**:
- (100), (010), (001) → same family {100}, multiplicity 3
- (110), (101), (011) → same family {110}, multiplicity 3
- (111) → unique family {111}, multiplicity 1
- (210), (201), (120), (102), (021), (012) → same family {210}, multiplicity 6

**Algorithm**:
1. For each hkl, compute sorted absolute values
2. Group by equality of sorted absolute values
3. Select representative (typically largest by lexicographic order)
4. Count members in each group

**Usage Example**:
```python
hkls = [(1,0,0), (0,1,0), (0,0,1), (1,1,0), (0,1,1)]
families = get_unique_families(hkls)
print(families)
# Output: {(1,0,0): 3, (1,1,0): 2}

# Negative indices
hkls2 = [(1,1,1), (-1,1,1), (1,-1,1)]
families2 = get_unique_families(hkls2)
print(families2)
# Output: {(1,1,1): 3}  # All are permutations
```

**Applications**:
- Structure factor multiplicity
- Systematic absences analysis
- Powder diffraction intensity calculations

### 3. Lattice Vector Calculation

#### `lattice_vec(lattice_param)`
**Purpose**: Calculate real-space lattice vectors (a₁, a₂, a₃) for any crystal system

**Input**:
- `lattice_param`: dict - Dictionary containing crystal system parameters

**Supported Crystal Systems**:

**1. Cubic (a = b = c, α = β = γ = 90°)**
```python
Parameters: {'type': 'cubic', 'a': float}
Vectors:
a₁ = [a, 0, 0]
a₂ = [0, a, 0]
a₃ = [0, 0, a]
```

**2. Tetragonal (a = b ≠ c, α = β = γ = 90°)**
```python
Parameters: {'type': 'tetragonal', 'a': float, 'b': float, 'c': float}
Vectors:
a₁ = [a, 0, 0]
a₂ = [0, b, 0]
a₃ = [0, 0, c]
```

**3. Trigonal/Hexagonal (a = b ≠ c, α = β = 90°, γ = 120°)**
```python
Parameters: {'type': 'trigonal', 'a': float, 'c': float}
Vectors:
a₁ = [a/2, -a√3/2, 0]
a₂ = [a/2,  a√3/2, 0]
a₃ = [0, 0, c]
Angle between a₁ and a₂: 120°
```

**4. Monoclinic (a ≠ b ≠ c, α = γ = 90°, β ≠ 90°)**
```python
Parameters: {
    'type': 'monoclinic',
    'a': float, 'b': float, 'c': float,
    'beta': float  # in radians
}
Vectors:
a₁ = [a, 0, 0]
a₂ = [0, b, 0]
a₃ = [c·cos(β), 0, c·sin(β)]
```

**5. Triclinic (a ≠ b ≠ c, α ≠ β ≠ γ)**
```python
Parameters: {
    'type': 'triclinic',
    'a': float, 'b': float, 'c': float,
    'alpha': float,  # in radians
    'beta': float,   # in radians
    'gamma': float   # in radians
}
Vectors:
a₁ = [a, 0, 0]
a₂ = [b·cos(γ), b·sin(γ), 0]
a₃ = [cₓ, cᵧ, cᵧ]
where:
  cₓ = c·cos(β)
  cᵧ = c·(cos(α) - cos(β)cos(γ))/sin(γ)
  cᵧ = √(c² - cₓ² - cᵧ²)
```

**Output**:
- `a1`, `a2`, `a3`: numpy arrays (3,) - Three lattice vectors

**Properties**:
- Vectors form a right-handed coordinate system
- Volume V = a₁ · (a₂ × a₃)
- Satisfy crystallographic conventions

**Usage Examples**:
```python
# Silicon (cubic, diamond structure)
Si_params = {'type': 'cubic', 'a': 5.43}
a1, a2, a3 = lattice_vec(Si_params)

# Titanium (hexagonal HCP)
Ti_params = {'type': 'trigonal', 'a': 2.95, 'c': 4.68}
a1, a2, a3 = lattice_vec(Ti_params)

# Check c/a ratio
c_over_a = np.linalg.norm(a3) / np.linalg.norm(a1)
print(f"c/a = {c_over_a:.3f}")  # Ideal HCP: 1.633
```

**Applications**:
- Crystal structure visualization
- Unit cell volume calculation
- Coordinate transformations
- Reciprocal lattice calculation

### 4. Reciprocal Lattice Calculation

#### `reciprocal_basis(a1, a2, a3)`
**Purpose**: Calculate reciprocal lattice basis vectors from real-space lattice vectors

**Input**:
- `a1`, `a2`, `a3`: numpy arrays (3,) - Real-space lattice vectors

**Output**:
- `b1`, `b2`, `b3`: numpy arrays (3,) - Reciprocal lattice vectors

**Mathematical Formula**:
```
b₁ = (a₂ × a₃) / V
b₂ = (a₃ × a₁) / V
b₃ = (a₁ × a₂) / V

where V = a₁ · (a₂ × a₃) is the unit cell volume
```

**Orthogonality Condition**:
```
aᵢ · bⱼ = δᵢⱼ (Kronecker delta)

Specifically:
a₁ · b₁ = 1,  a₁ · b₂ = 0,  a₁ · b₃ = 0
a₂ · b₁ = 0,  a₂ · b₂ = 1,  a₂ · b₃ = 0
a₃ · b₁ = 0,  a₃ · b₂ = 0,  a₃ · b₃ = 1
```

**Properties**:
- Reciprocal cell volume: V* = 1/V
- Verification: V × V* = 1
- For cubic: b₁ = (1/a, 0, 0), b₂ = (0, 1/a, 0), b₃ = (0, 0, 1/a)

**Convention Note**:
This implementation uses the **crystallographic convention** (without 2π).
For physics convention, multiply by 2π: b_physics = 2π × b_crystal

**Usage Example**:
```python
# Cubic crystal
a = 5.0  # Angstroms
a1 = np.array([a, 0, 0])
a2 = np.array([0, a, 0])
a3 = np.array([0, 0, a])

b1, b2, b3 = reciprocal_basis(a1, a2, a3)
print(b1)  # [0.2, 0, 0] = [1/a, 0, 0]

# Verify orthogonality
print(np.dot(a1, b1))  # 1.0
print(np.dot(a1, b2))  # 0.0

# Volume relationship
V = np.dot(a1, np.cross(a2, a3))
V_star = np.dot(b1, np.cross(b2, b3))
print(f"V × V* = {V * V_star}")  # 1.0
```

**Applications**:
- Diffraction pattern indexing
- d-spacing calculations: d_hkl = 1/|G_hkl|
- Structure factor calculations
- Ewald sphere construction
- Reciprocal space mapping

### 5. Utility Functions

#### `array2tuple(arr, decimals=2)`
**Purpose**: Convert numpy array to tuple with rounded elements

**Input**:
- `arr`: numpy array or list - Numerical values
- `decimals`: int - Decimal places for rounding (default: 2)

**Output**:
- `tuple`: Rounded values as tuple

**Features**:
- Handles both arrays and lists
- Configurable precision
- Used internally for Miller index comparison

**Usage Example**:
```python
arr = np.array([1.234567, 2.891011, 3.141592])
result = array2tuple(arr, decimals=3)
print(result)  # (1.235, 2.891, 3.142)
```

## Mathematical Background

### Crystal Systems

**7 Crystal Systems, 14 Bravais Lattices**:

1. **Cubic** (a=b=c, α=β=γ=90°)
   - Simple, Body-centered, Face-centered

2. **Tetragonal** (a=b≠c, α=β=γ=90°)
   - Simple, Body-centered

3. **Hexagonal** (a=b≠c, α=β=90°, γ=120°)
   - Simple

4. **Trigonal/Rhombohedral** (a=b=c, α=β=γ≠90°)
   - Simple

5. **Orthorhombic** (a≠b≠c, α=β=γ=90°)
   - Simple, Body-centered, Face-centered, Base-centered

6. **Monoclinic** (a≠b≠c, α=γ=90°≠β)
   - Simple, Base-centered

7. **Triclinic** (a≠b≠c, α≠β≠γ)
   - Simple

### Miller Indices

**Notation**:
- (hkl): Plane notation
- [uvw]: Direction notation
- {hkl}: Family of equivalent planes
- <uvw>: Family of equivalent directions

**Properties**:
- (hkl) plane is perpendicular to [hkl] direction (in cubic)
- d-spacing: d_hkl = 1/|h·b₁ + k·b₂ + l·b₃|
- Plane intercepts: a₁/h, a₂/k, a₃/l

### Reciprocal Lattice

**Definition**:
The reciprocal lattice is defined such that:
```
G_hkl = h·b₁ + k·b₂ + l·b₃
```
is perpendicular to the (hkl) plane in real space.

**Bragg's Law**:
```
nλ = 2d_hkl sin(θ)
```
where:
- n: order of reflection
- λ: wavelength
- d_hkl: plane spacing
- θ: Bragg angle

**d-spacing Formula**:
```
d_hkl = 1 / |G_hkl|
```

For cubic: 
```
d_hkl = a / √(h² + k² + l²)
```

For hexagonal:
```
1/d² = 4(h² + hk + k²)/(3a²) + l²/c²
```

### Metric Tensor

**Direct Space Metric Tensor**:
```
G_ij = a_i · a_j
```

For general crystal:
```
G = [a₁·a₁  a₁·a₂  a₁·a₃]
    [a₂·a₁  a₂·a₂  a₂·a₃]
    [a₃·a₁  a₃·a₂  a₃·a₃]
```

**Reciprocal Space Metric Tensor**:
```
G* = G⁻¹
```

**Vector length**:
```
|v|² = v^T · G · v (direct space)
|g|² = g^T · G* · g (reciprocal space)
```

## Integration Examples

### With Orilib (Orientation Analysis)

```python
from crystlib import lattice_vec, reciprocal_basis
from orilib import np_euler_matrix, quat_to_mat

# Define crystal structure
cubic_params = {'type': 'cubic', 'a': 3.5}
a1, a2, a3 = lattice_vec(cubic_params)

# Create basis matrix
L = np.column_stack([a1, a2, a3])

# Define crystal orientation (Euler angles)
euler = [0.5, 0.3, 0.1]  # radians
g = np_euler_matrix(*euler)

# Transform crystal direction to sample frame
direction_crystal = np.array([1, 0, 0])  # [100]
direction_sample = g @ direction_crystal

print("Crystal [100] in sample frame:", direction_sample)
```

### With Projlib (Pole Figures)

```python
from crystlib import lattice_vec, reciprocal_basis, generate_hkls
from projlib import gen_dirs_norms, equalarea_directions

# Setup crystal
cubic = {'type': 'cubic', 'a': 3.5}
a1, a2, a3 = lattice_vec(cubic)
b1, b2, b3 = reciprocal_basis(a1, a2, a3)

L = np.column_stack([a1, a2, a3])
Lr = np.column_stack([b1, b2, b3])

# Generate symmetry operations (simplified)
symops = [np.eye(3)]  # Add more for full symmetry

# Generate crystal directions
dirs, norms = gen_dirs_norms(
    L, Lr,
    uvws=[[1,0,0], [1,1,0], [1,1,1]],
    hkls=[[1,0,0], [1,1,0], [1,1,1]],
    symops=symops
)

# Project for pole figure
dir_vectors = np.array([d['vector'] for d in dirs]).T
proj_coords = equalarea_directions(dir_vectors)
```

### With Plotlib (Visualization)

```python
from crystlib import generate_hkls
from plotlib import plotter

# Generate reflections
syms = [np.eye(3)]  # Add symmetry operations
hkls, hkls2, fam = generate_hkls(3, syms)

# Setup plotter
p = plotter()
p.setAttributes(
    norms=[list(hkl) for hkl in hkls[:20]],  # First 20 planes
    ProjType='equalarea',
    sphere='half'
)

p.plotProj()
# p.plotDirsNorms()  # If method available
p.figsave(fname='pole_figure_reflections.png')
```

### Complete Diffraction Analysis

```python
from crystlib import *
import numpy as np

# Define crystal
al_params = {'type': 'cubic', 'a': 4.05}  # Aluminum
a1, a2, a3 = lattice_vec(al_params)
b1, b2, b3 = reciprocal_basis(a1, a2, a3)

# Generate reflections
cubic_syms = [np.eye(3)]  # Add full cubic symmetry
hkls, hkls2, fam = generate_hkls(3, cubic_syms)

# Calculate d-spacings
wavelength = 1.5406  # Cu Kα

print("Reflection  d-spacing (Å)  2θ (°)")
print("-" * 40)

for hkl in hkls[:15]:
    G = hkl[0]*b1 + hkl[1]*b2 + hkl[2]*b3
    d = 1.0 / np.linalg.norm(G)
    
    # Bragg angle
    sin_theta = wavelength / (2*d)
    if sin_theta <= 1.0:
        theta = np.arcsin(sin_theta)
        two_theta = 2 * np.degrees(theta)
        print(f"({int(hkl[0])}{int(hkl[1])}{int(hkl[2])})  {d:8.4f}      {two_theta:7.2f}")
```

## Best Practices

### 1. Symmetry Operations
```python
# Define complete symmetry for your crystal system
# Cubic: 48 operations
# Hexagonal: 24 operations
# Tetragonal: 16 operations
# etc.

# Verify symmetry matrices are orthogonal
for sym in syms:
    assert np.allclose(sym @ sym.T, np.eye(3))
    assert np.allclose(np.linalg.det(sym), 1)
```

### 2. Angle Units
```python
# Always use radians for lattice_vec
import numpy as np

# Convert from degrees
beta_deg = 120
beta_rad = np.radians(beta_deg)

params = {'type': 'monoclinic', 'beta': beta_rad}
```

### 3. Reciprocal Space Calculations
```python
# Calculate reciprocal basis once, reuse
a1, a2, a3 = lattice_vec(params)
b1, b2, b3 = reciprocal_basis(a1, a2, a3)

# Reuse for all d-spacing calculations
for hkl in hkls:
    G = hkl[0]*b1 + hkl[1]*b2 + hkl[2]*b3
    d = 1.0 / np.linalg.norm(G)
```

### 4. Validation
```python
# Verify reciprocal lattice
V = np.dot(a1, np.cross(a2, a3))
V_star = np.dot(b1, np.cross(b2, b3))
assert np.isclose(V * V_star, 1.0)

# Check orthogonality
for i, ai in enumerate([a1, a2, a3]):
    for j, bj in enumerate([b1, b2, b3]):
        dot = np.dot(ai, bj)
        expected = 1.0 if i == j else 0.0
        assert np.isclose(dot, expected)
```

## Common Applications

### 1. X-ray Diffraction Analysis
- Generate allowed reflections
- Calculate d-spacings
- Predict 2θ angles
- Identify systematic absences

### 2. EBSD Pattern Indexing
- Generate crystal directions
- Apply symmetry operations
- Match to experimental patterns
- Determine orientation

### 3. Texture Analysis
- Define crystal families
- Calculate pole figures
- Determine preferred orientations
- Quantify texture strength

### 4. Structure Factor Calculations
- Identify unique reflections
- Determine multiplicities
- Calculate intensities
- Analyze systematic absences

## Performance Considerations

### Memory Usage
- hklmax=3: ~100-500 reflections
- hklmax=5: ~1000-5000 reflections
- hklmax=10: ~10000-50000 reflections

### Optimization Tips
1. Use appropriate hklmax for your needs
2. Cache symmetry operations
3. Pre-compute reciprocal lattice
4. Limit to relevant families for analysis

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| KeyError in lattice_vec | Wrong 'type' key | Check spelling, use lowercase |
| Wrong angles | Degrees vs radians | Use np.radians() |
| No hkls generated | hklmax too small | Increase hklmax |
| Too many hkls | hklmax too large | Decrease hklmax |
| Wrong d-spacings | Wrong reciprocal lattice | Verify calculation |
| Symmetry errors | Non-orthogonal matrices | Check det=1, orthogonality |

## Dependencies

**Required**:
- numpy

**Optional** (for integration):
- orilib (orientation analysis)
- projlib (projections)
- plotlib (visualization)

## Version Information

**Python**: 2.7+ / 3.x compatible
**File size**: ~19 KB
**Functions**: 7 documented functions
**Lines**: ~522 lines

## References

### Crystallography
- International Tables for Crystallography, Vol. A: Space-Group Symmetry
- Giacovazzo, C. et al. (2011). *Fundamentals of Crystallography*. Oxford University Press.

### Diffraction
- Cullity, B. D., & Stock, S. R. (2014). *Elements of X-ray Diffraction*. Pearson.

### EBSD
- Schwartz, A. J., et al. (2009). *Electron Backscatter Diffraction in Materials Science*. Springer.

## Summary

The crystlib module provides:
- ✅ Miller indices generation with symmetry
- ✅ Support for all 7 crystal systems
- ✅ Lattice vector calculations
- ✅ Reciprocal lattice calculations
- ✅ Family grouping of planes/directions
- ✅ Metric tensor handling
- ✅ Comprehensive documentation

**Perfect for**: Diffraction analysis, EBSD indexing, texture analysis, crystal structure calculations, reciprocal space mapping.

## For More Information

- **Quick start**: See crystlib_quick_reference.md
- **Detailed examples**: See crystlib_commented.py
- **Integration**: See orilib, projlib, plotlib documentation
- **Crystallography**: Consult standard crystallography textbooks

---

**Module size**: ~19 KB, 7 functions for comprehensive crystallographic calculations.

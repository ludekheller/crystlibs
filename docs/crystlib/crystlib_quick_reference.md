# Crystlib.py Quick Reference Guide

## Overview

The **crystlib.py** module provides utility functions for crystallographic calculations including Miller indices generation with symmetry operations, lattice vector calculations for various crystal systems, and reciprocal basis vector calculations.

## Quick Start

```python
import numpy as np
from crystlib import *

# Generate Miller indices for cubic crystal
identity = np.eye(3)
symops = [identity]  # Add your symmetry operations
hkls, hkls2, families = generate_hkls(2, symops)

# Calculate lattice vectors for cubic crystal
cubic_params = {'type': 'cubic', 'a': 5.43}
a1, a2, a3 = lattice_vec(cubic_params)

# Calculate reciprocal lattice
b1, b2, b3 = reciprocal_basis(a1, a2, a3)
```

## Most Commonly Used Functions

### 1. Miller Indices Generation

#### `generate_hkls(hklmax, syms, hkls=[])`
**Purpose**: Generate unique Miller indices (hkl) considering crystal symmetry

**Input**:
- `hklmax`: int - Maximum Miller index value
- `syms`: list - 3×3 symmetry operation matrices
- `hkls`: list - Optional custom list (default: [])

**Output**:
- `hkls`: list of tuples - Unique Miller indices
- `hkls2`: dict - Maps each hkl to symmetry equivalents
- `fam`: dict - Families with multiplicities

**Usage**:
```python
import numpy as np

# Cubic symmetry (simplified example)
identity = np.eye(3)
rot_90_x = np.array([[1,  0,  0],
                     [0,  0, -1],
                     [0,  1,  0]])
rot_90_y = np.array([[ 0,  0,  1],
                     [ 0,  1,  0],
                     [-1,  0,  0]])
rot_90_z = np.array([[0, -1,  0],
                     [1,  0,  0],
                     [0,  0,  1]])

syms = [identity, rot_90_x, rot_90_y, rot_90_z]

# Generate up to (222)
hkls, hkls2, families = generate_hkls(2, syms)

print(f"Unique planes: {len(hkls)}")
print("First few:", hkls[:5])
print("Families:", families)

# Get symmetry equivalents of (100)
equiv_100 = hkls2.get((1.0, 0.0, 0.0))
print("Equivalents of (100):", equiv_100)
```

#### `generate_hkls01(hklmax, syms, hkls=[])`
**Purpose**: Alternative version without rounding (for exact comparison)

**Usage**: Same as `generate_hkls()` but uses exact floating point comparison

```python
# Tetragonal symmetry
syms_tetra = [identity, rot_90_z]
hkls, hkls2, fam = generate_hkls01(3, syms_tetra)
```

#### `generate_hkls02(hklmax, syms, G, hkls=[])`
**Purpose**: Generate hkls with metric tensor transformation

**Additional Input**:
- `G`: numpy array (3×3) - Metric tensor matrix

**Usage**:
```python
# Hexagonal metric tensor
a = 3.0
c = 5.0
G = np.array([[a**2,     -a**2/2,  0],
              [-a**2/2,   a**2,     0],
              [0,         0,        c**2]])

# Hexagonal symmetry (60° rotation)
rot_60_z = np.array([[0.5, -np.sqrt(3)/2, 0],
                     [np.sqrt(3)/2, 0.5, 0],
                     [0, 0, 1]])

syms_hex = [identity, rot_60_z]
hkls, hkls2, fam = generate_hkls02(2, syms_hex, G)
```

### 2. Lattice Vector Calculations

#### `lattice_vec(lattice_param)`
**Purpose**: Calculate real space lattice vectors for various crystal systems

**Input**: `lattice_param` - Dictionary with crystal system parameters

**Crystal Systems Supported**:

**Cubic**:
```python
cubic_params = {
    'type': 'cubic',
    'a': 5.43  # Angstroms
}
a1, a2, a3 = lattice_vec(cubic_params)
# Returns: orthogonal vectors of length a
```

**Tetragonal**:
```python
tetra_params = {
    'type': 'tetragonal',
    'a': 3.0,
    'b': 3.0,
    'c': 5.0
}
a1, a2, a3 = lattice_vec(tetra_params)
```

**Trigonal/Hexagonal**:
```python
hex_params = {
    'type': 'trigonal',
    'a': 3.0,
    'c': 5.0
}
a1, a2, a3 = lattice_vec(hex_params)
# Returns 120° angle between a1 and a2
```

**Monoclinic**:
```python
mono_params = {
    'type': 'monoclinic',
    'a': 5.0,
    'b': 6.0,
    'c': 7.0,
    'beta': np.radians(120)  # β angle in radians
}
a1, a2, a3 = lattice_vec(mono_params)
```

**Triclinic**:
```python
triclinic_params = {
    'type': 'triclinic',
    'a': 5.0,
    'b': 6.0,
    'c': 7.0,
    'alpha': np.radians(90),   # α angle
    'beta':  np.radians(95),   # β angle
    'gamma': np.radians(100)   # γ angle
}
a1, a2, a3 = lattice_vec(triclinic_params)
```

**Output**: Three lattice vectors as numpy arrays (3,)

### 3. Reciprocal Lattice

#### `reciprocal_basis(a1, a2, a3)`
**Purpose**: Calculate reciprocal lattice vectors

**Input**:
- `a1`, `a2`, `a3`: numpy arrays (3,) - Real space lattice vectors

**Output**:
- `b1`, `b2`, `b3`: numpy arrays (3,) - Reciprocal lattice vectors

**Formula**: 
```
b1 = (a2 × a3) / (a1 · (a2 × a3))
b2 = (a3 × a1) / (a2 · (a3 × a1))
b3 = (a1 × a2) / (a3 · (a1 × a2))
```

**Orthogonality**: a_i · b_j = δ_ij (Kronecker delta)

**Usage**:
```python
# Cubic example
a = 5.0
a1 = np.array([a, 0, 0])
a2 = np.array([0, a, 0])
a3 = np.array([0, 0, a])

b1, b2, b3 = reciprocal_basis(a1, a2, a3)
print("b1:", b1)  # [1/a, 0, 0]

# Verify orthogonality
print("a1 · b1 =", np.dot(a1, b1))  # 1.0
print("a1 · b2 =", np.dot(a1, b2))  # 0.0

# Hexagonal example
hex_params = {'type': 'trigonal', 'a': 3.0, 'c': 5.0}
a1_hex, a2_hex, a3_hex = lattice_vec(hex_params)
b1_hex, b2_hex, b3_hex = reciprocal_basis(a1_hex, a2_hex, a3_hex)

print("c* length:", np.linalg.norm(b3_hex))
```

### 4. Utility Functions

#### `array2tuple(arr, decimals=2)`
**Purpose**: Convert array to rounded tuple

**Input**:
- `arr`: numpy array or list
- `decimals`: int - Decimal places (default: 2)

**Output**: tuple - Rounded values

**Usage**:
```python
arr = np.array([1.234, 2.567, 3.891])
result = array2tuple(arr, decimals=2)
print(result)  # (1.23, 2.57, 3.89)
```

#### `get_unique_families(hkls)`
**Purpose**: Group Miller indices into families by permutation symmetry

**Input**: hkls - List of Miller index tuples

**Output**: dict - {representative_hkl: multiplicity}

**Usage**:
```python
hkls = [(1,0,0), (0,1,0), (0,0,1), (1,1,0), (0,1,1)]
families = get_unique_families(hkls)
print(families)
# {(1,0,0): 3, (1,1,0): 2}
# (100), (010), (001) → same family
# (110), (011) → same family
```

## Complete Workflow Examples

### Example 1: Cubic Crystal Analysis

```python
import numpy as np
from crystlib import *

# Define silicon (cubic, a = 5.43 Å)
lattice_params = {'type': 'cubic', 'a': 5.43}

# Get real space lattice
a1, a2, a3 = lattice_vec(lattice_params)
print("Lattice vectors:")
print("a1:", a1)
print("a2:", a2)
print("a3:", a3)

# Calculate reciprocal lattice
b1, b2, b3 = reciprocal_basis(a1, a2, a3)
print("\nReciprocal lattice:")
print("b1:", b1)
print("b1 magnitude:", np.linalg.norm(b1))

# Volume of unit cell
V = np.dot(a1, np.cross(a2, a3))
print(f"\nUnit cell volume: {V:.3f} Ų")

# Volume of reciprocal cell
V_recip = np.dot(b1, np.cross(b2, b3))
print(f"Reciprocal cell volume: {V_recip:.6f} ų⁻³")
print(f"V × V* = {V * V_recip:.6f}")  # Should be 1.0
```

### Example 2: Generate Diffraction Planes with Symmetry

```python
import numpy as np
from crystlib import *

# Define cubic symmetry operations (48 operations total)
# Simplified example with 4 operations
identity = np.eye(3)

# 90° rotations
rot_90_x = np.array([[1,  0,  0],
                     [0,  0, -1],
                     [0,  1,  0]])

rot_90_y = np.array([[ 0,  0,  1],
                     [ 0,  1,  0],
                     [-1,  0,  0]])

rot_90_z = np.array([[0, -1,  0],
                     [1,  0,  0],
                     [0,  0,  1]])

# 180° rotations
rot_180_x = rot_90_x @ rot_90_x
rot_180_y = rot_90_y @ rot_90_y
rot_180_z = rot_90_z @ rot_90_z

# Combine
cubic_syms = [identity, rot_90_x, rot_90_y, rot_90_z,
              rot_180_x, rot_180_y, rot_180_z]

# Generate reflections
hkls, hkls2, families = generate_hkls(3, cubic_syms)

print(f"Total unique reflections up to (333): {len(hkls)}")
print("\nFirst 10 reflections:")
for i, hkl in enumerate(hkls[:10]):
    print(f"{i+1}. {hkl}")

print("\nFamilies:")
for family, mult in sorted(families.items())[:5]:
    print(f"{family}: multiplicity {mult}")

# Check specific reflection
if (1.0, 0.0, 0.0) in hkls2:
    equiv = hkls2[(1.0, 0.0, 0.0)]
    print(f"\n(100) has {len(equiv)} equivalent reflections:")
    print(equiv[:6])  # Show first 6
```

### Example 3: Hexagonal Crystal (e.g., Titanium)

```python
import numpy as np
from crystlib import *

# Titanium (hexagonal, a=2.95 Å, c=4.68 Å)
hex_params = {
    'type': 'trigonal',
    'a': 2.95,
    'c': 4.68
}

# Get lattice vectors
a1, a2, a3 = lattice_vec(hex_params)

print("Real space lattice vectors:")
print(f"a1: {a1}")
print(f"a2: {a2}")
print(f"a3 (c-axis): {a3}")

# Check angles
angle_a1_a2 = np.arccos(np.dot(a1, a2) / 
                        (np.linalg.norm(a1) * np.linalg.norm(a2)))
print(f"\nAngle between a1 and a2: {np.degrees(angle_a1_a2):.1f}°")
# Should be 120°

# c/a ratio
c_over_a = np.linalg.norm(a3) / np.linalg.norm(a1)
print(f"c/a ratio: {c_over_a:.3f}")

# Reciprocal lattice
b1, b2, b3 = reciprocal_basis(a1, a2, a3)
print("\nReciprocal lattice:")
print(f"b1: {b1}")
print(f"b2: {b2}")
print(f"b3 (c*): {b3}")

# d-spacing for (002) reflection
hkl = np.array([0, 0, 2])
d_002 = 1.0 / np.linalg.norm(hkl[0]*b1 + hkl[1]*b2 + hkl[2]*b3)
print(f"\nd-spacing for (002): {d_002:.3f} Å")
```

### Example 4: Diffraction Analysis

```python
import numpy as np
from crystlib import *

# Setup crystal
cubic_params = {'type': 'cubic', 'a': 3.5}
a1, a2, a3 = lattice_vec(cubic_params)
b1, b2, b3 = reciprocal_basis(a1, a2, a3)

# Calculate d-spacings for common reflections
reflections = [(1,0,0), (1,1,0), (1,1,1), (2,0,0), (2,1,1), (2,2,0)]

print("Reflection  d-spacing (Å)  |G| (ų⁻¹)")
print("-" * 40)

for hkl in reflections:
    # Reciprocal lattice vector
    G = hkl[0]*b1 + hkl[1]*b2 + hkl[2]*b3
    G_mag = np.linalg.norm(G)
    
    # d-spacing
    d = 1.0 / G_mag if G_mag > 0 else np.inf
    
    print(f"({hkl[0]}{hkl[1]}{hkl[2]})       {d:8.4f}      {G_mag:6.4f}")

# X-ray wavelength (Cu Kα = 1.5406 Å)
wavelength = 1.5406

print(f"\nBragg angles (λ = {wavelength} Å, Cu Kα):")
print("Reflection  2θ (degrees)")
print("-" * 30)

for hkl in reflections:
    G = hkl[0]*b1 + hkl[1]*b2 + hkl[2]*b3
    d = 1.0 / np.linalg.norm(G)
    
    # Bragg's law: λ = 2d sin(θ)
    sin_theta = wavelength / (2 * d)
    if sin_theta <= 1.0:
        theta = np.arcsin(sin_theta)
        two_theta = 2 * np.degrees(theta)
        print(f"({hkl[0]}{hkl[1]}{hkl[2]})       {two_theta:7.2f}")
```

### Example 5: Miller Index Families

```python
import numpy as np
from crystlib import *

# Generate all equivalents for a cubic crystal
identity = np.eye(3)

# Add more cubic symmetry operations
# (in practice, use all 48 operations)
syms = [identity]  # Simplified

hkls, hkls2, families = generate_hkls(2, syms)

print("Miller Index Families:")
print("=" * 50)

for family_rep, multiplicity in sorted(families.items()):
    # Find all members of this family
    members = []
    for key, equiv_list in hkls2.items():
        for equiv in equiv_list:
            if get_unique_families([equiv]) == {family_rep: 1}:
                if key not in members:
                    members.append(key)
    
    print(f"\nFamily {{{family_rep[0]}{family_rep[1]}{family_rep[2]}}}")
    print(f"Multiplicity: {multiplicity}")
    print(f"Members: {members[:6]}...")  # Show first 6
```

## Function Parameters Quick Reference

### generate_hkls()
```python
generate_hkls(
    hklmax,      # int: max Miller index
    syms,        # list: symmetry matrices (3×3)
    hkls=[]      # list: optional custom hkls
)
→ (hkls, hkls2, fam)
```

### lattice_vec()
```python
lattice_vec(
    lattice_param  # dict: crystal parameters
)
→ (a1, a2, a3)  # numpy arrays (3,)
```

### reciprocal_basis()
```python
reciprocal_basis(
    a1, a2, a3   # numpy arrays (3,): lattice vectors
)
→ (b1, b2, b3)  # numpy arrays (3,): reciprocal vectors
```

## Crystal System Parameters

| System | Parameters Required | Example |
|--------|-------------------|---------|
| Cubic | a | {'type': 'cubic', 'a': 5.0} |
| Tetragonal | a, b, c | {'type': 'tetragonal', 'a': 3.0, 'b': 3.0, 'c': 5.0} |
| Trigonal/Hex | a, c | {'type': 'trigonal', 'a': 3.0, 'c': 5.0} |
| Monoclinic | a, b, c, beta | {'type': 'monoclinic', 'a': 5.0, 'b': 6.0, 'c': 7.0, 'beta': 1.57} |
| Triclinic | a, b, c, α, β, γ | {'type': 'triclinic', 'a': 5.0, 'b': 6.0, 'c': 7.0, 'alpha': 1.57, 'beta': 1.57, 'gamma': 1.57} |

**Note**: Angles must be in radians. Use `np.radians(degrees)` to convert.

## Common Crystal Structures

### Cubic Systems (a only)
```python
# Silicon (diamond cubic)
Si = {'type': 'cubic', 'a': 5.43}

# Aluminum (FCC)
Al = {'type': 'cubic', 'a': 4.05}

# Iron (BCC, room temp)
Fe = {'type': 'cubic', 'a': 2.87}

# Copper (FCC)
Cu = {'type': 'cubic', 'a': 3.61}
```

### Hexagonal Systems (a, c)
```python
# Titanium (HCP)
Ti = {'type': 'trigonal', 'a': 2.95, 'c': 4.68}

# Zinc (HCP)
Zn = {'type': 'trigonal', 'a': 2.66, 'c': 4.95}

# Magnesium (HCP)
Mg = {'type': 'trigonal', 'a': 3.21, 'c': 5.21}
```

## Mathematical Formulas

### Reciprocal Lattice Vectors
```
b₁ = (a₂ × a₃) / V
b₂ = (a₃ × a₁) / V
b₃ = (a₁ × a₂) / V

where V = a₁ · (a₂ × a₃) is the unit cell volume
```

### Orthogonality Condition
```
aᵢ · bⱼ = δᵢⱼ

where δᵢⱼ is the Kronecker delta (1 if i=j, 0 otherwise)
```

### d-spacing Calculation
```
d_hkl = 1 / |G_hkl|

where G_hkl = h·b₁ + k·b₂ + l·b₃
```

### Bragg's Law
```
nλ = 2d sin(θ)

where:
- n = order of diffraction (usually 1)
- λ = X-ray wavelength
- d = d-spacing
- θ = Bragg angle
```

### Unit Cell Volume
```
V = a₁ · (a₂ × a₃) = |det([a₁ a₂ a₃])|
```

### Reciprocal Cell Volume
```
V* = b₁ · (b₂ × b₃) = 1/V
```

## Common Pitfalls

### 1. Angle Units
```python
# Wrong - using degrees
params = {'type': 'monoclinic', 'beta': 120}

# Correct - convert to radians
import numpy as np
params = {'type': 'monoclinic', 'beta': np.radians(120)}
```

### 2. Symmetry Matrix Format
```python
# Ensure symmetry matrices are numpy arrays
syms = [np.eye(3), np.array([[0,-1,0],[1,0,0],[0,0,1]])]

# Not Python lists
syms = [[1,0,0], [0,1,0], [0,0,1]]  # Wrong!
```

### 3. Miller Index Format
```python
# Functions return tuples, not lists
hkls, hkls2, fam = generate_hkls(2, syms)
print(type(hkls[0]))  # tuple

# Access like: hkls2[(1.0, 0.0, 0.0)]
# Not: hkls2[[1, 0, 0]]
```

### 4. Reciprocal Lattice Convention
```python
# This implementation uses crystallographic convention
# (without 2π factor)
# For physics convention, multiply by 2π:
b1_physics = 2 * np.pi * b1
```

## Integration with Other Modules

### With Projlib
```python
from crystlib import lattice_vec, reciprocal_basis
from projlib import gen_dirs_norms

# Setup crystal
cubic = {'type': 'cubic', 'a': 3.5}
a1, a2, a3 = lattice_vec(cubic)
b1, b2, b3 = reciprocal_basis(a1, a2, a3)

# Create basis matrices
L = np.column_stack([a1, a2, a3])
Lr = np.column_stack([b1, b2, b3])

# Generate directions and normals
dirs, norms = gen_dirs_norms(
    L, Lr,
    uvws=[[1,0,0], [1,1,0]],
    hkls=[[1,0,0], [1,1,1]]
)
```

### With Plotlib
```python
from crystlib import generate_hkls
from plotlib import plotter

# Generate reflections
hkls, hkls2, fam = generate_hkls(2, syms)

# Plot pole figure with reflections
p = plotter()
p.setAttributes(norms=[list(hkl) for hkl in hkls[:10]])
p.plotProj(ProjType='equalarea')
```

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "KeyError" in lattice_vec | Check 'type' key spelling and case |
| Angles in wrong units | Convert degrees to radians: `np.radians(deg)` |
| No symmetry equivalents | Verify symmetry matrices are orthogonal |
| Wrong d-spacings | Check reciprocal lattice calculation |
| Empty hkls list | Ensure hklmax > 0 and syms is not empty |

## Performance Tips

1. **Pre-compute symmetry operations**: Generate once, reuse many times
2. **Limit hklmax**: Start small (2-3), increase as needed
3. **Cache reciprocal lattice**: Calculate once per crystal
4. **Use appropriate function**: `generate_hkls()` vs `generate_hkls02()` based on needs

## For More Information

- **Full documentation**: See `crystlib_commented.py`
- **Examples**: Check function docstrings
- **Integration**: See projlib and plotlib documentation
- **Crystallography**: Refer to standard crystallography textbooks

---

**Pro tip**: For EBSD analysis, combine crystlib with orilib and projlib to generate full orientation relationships and pole figures for any crystal system!

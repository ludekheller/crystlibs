# GetPhases.py Quick Reference Guide

## Overview

The **getphases.py** module provides comprehensive utilities for managing crystallographic phases and their transformations, with specific focus on shape memory alloys like NiTi. It handles phase definitions, orientation relationships, deformation gradients, twin systems, and habit plane variants.

## Quick Start

```python
from getphases import getPhases
import numpy as np

# Initialize phase manager
gp = getPhases()

# Load phases from CIF files
gp.fromCif('austenite.cif', 'A')
gp.fromCif('martensite.cif', 'M')

# Or define from lattice parameters
gp.fromLatticeParams('A', a=3.015)  # Cubic austenite
gp.fromLatticeParams('M', a=2.889, b=4.120, c=4.622,
                      beta=np.radians(96.8))  # Monoclinic martensite

# Get orientation relationship
gp.getOR(name='NiTi', OR='NiTi')

# Calculate deformation gradients
gp.getDefGrad(name='NiTi')

# Get habit plane variants
gp.getHBVs(name='NiTi')
```

## Most Commonly Used Methods

### 1. Phase Definition

#### `fromCif(file, phasename)`
**Purpose**: Load crystal structure from CIF file

**Usage**:
```python
gp = getPhases()

# Load B2 austenite
gp.fromCif('B2_NiTi.cif', 'A')

# Load B19' martensite
gp.fromCif('B19prime_NiTi.cif', 'M')

# Access phase information
print("Austenite lattice parameter:", gp.phases['A']['a'])
print("Number of symmetry operations:", len(gp.phases['A']['symops']))
```

#### `fromLatticeParams(phasename, a, b=None, c=None, alpha=π/2, beta=π/2, gamma=π/2)`
**Purpose**: Build lattice matrices from parameters

**Usage**:
```python
# Cubic B2 austenite (a = 3.015 Ã…)
gp.fromLatticeParams('A', a=3.015)

# Monoclinic B19' martensite
gp.fromLatticeParams('M', 
    a=2.889,  # Ã…
    b=4.120,  # Ã…
    c=4.622,  # Ã…
    beta=np.radians(96.8)  # β angle
)

# Hexagonal (for comparison)
gp.fromLatticeParams('Ti', 
    a=2.95, 
    c=4.68,
    gamma=np.radians(120)
)

# Access lattice matrices
L_austenite = gp.phases['A']['L']      # Direct lattice
Lr_martensite = gp.phases['M']['Lr']   # Reciprocal lattice
```

### 2. Orientation Relationships

#### `getOR(name=None, OR='NiTi')`
**Purpose**: Get orientation relationship matrices between phases

**Usage**:
```python
# Get NiTi orientation relationship
gp.getOR(name='NiTi', OR='NiTi')

# Access correspondence matrices
Cd = gp.OR['NiTi']['Cd']    # Direction correspondence (M→A)
CId = gp.OR['NiTi']['CId']  # Inverse direction correspondence (A→M)
Cp = gp.OR['NiTi']['Cp']    # Plane correspondence (M→A)
CIp = gp.OR['NiTi']['CIp']  # Inverse plane correspondence (A→M)

# Example: Map martensite [100] to austenite frame
mart_dir = np.array([1, 0, 0])
aust_dir = Cd @ mart_dir
print("Martensite [100] corresponds to austenite:", aust_dir)
```

### 3. Deformation Gradients

#### `getDefGrad(name=None, OR='NiTi')`
**Purpose**: Compute deformation gradients for stress-free transformation

**Usage**:
```python
# Calculate deformation gradients
gp.getDefGrad(name='NiTi')

# Access deformation gradient for variant 1
F1 = gp.defGrad['NiTi'][0]['F']
print("Deformation gradient shape:", F1.shape)  # (3, 3)

# Calculate transformation strain
epsilon = (F1.T @ F1 - np.eye(3)) / 2
print("Transformation strain:\n", epsilon)

# All 12 correspondence variants
print("Number of variants:", len(gp.defGrad['NiTi']))
for i, variant in enumerate(gp.defGrad['NiTi']):
    F = variant['F']
    print(f"Variant {i}: det(F) = {np.linalg.det(F):.6f}")
```

#### `getTrStrain(phase=None, name='NiTi', oris=None)`
**Purpose**: Calculate transformation strain for orientations

**Usage**:
```python
from orilib import np_eulers_matrices

# Random orientations
euler = np.random.rand(100, 3) * [360, 180, 360]
oris = np_eulers_matrices(euler, deg=True)

# Calculate transformation strains
strains = gp.getTrStrain(phase='A', name='NiTi', oris=oris)

print("Transformation strain shape:", strains.shape)  # (100, 6)
# 6 components: [ε11, ε22, ε33, ε23, ε13, ε12]

# Analyze strain distribution
import matplotlib.pyplot as plt
plt.hist(strains[:, 0], bins=30)  # ε11 distribution
plt.xlabel('ε11')
plt.ylabel('Frequency')
plt.show()
```

### 4. Habit Plane Variants

#### `getHBVs(name='NiTi')`
**Purpose**: Calculate habit plane variants

**Usage**:
```python
# Get habit plane variants
gp.getHBVs(name='NiTi')

# Access Type I twins
print("Type I twin systems:", len(gp.twinSys['NiTi']['I']))

# Access habit plane normals
for i, hbv in enumerate(gp.vardict['NiTi']['HBV']):
    n = hbv['n_hpv']  # Habit plane normal
    m = hbv['m_hpv']  # Shear direction
    print(f"HBV {i}: n = {n}, m = {m}")

# Plot habit plane variants
gp.plotHBVII011(name='NiTi', ProjType='equalarea')
```

#### `getHPVII011(name='NiTi')`
**Purpose**: Get {011} Type II twin habit plane variants

**Usage**:
```python
# Calculate {011} Type II twin HPVs
gp.getHPVII011(name='NiTi')

# Access HPV data
hpvs = gp.vardict['NiTi']['HPVII011']
print(f"Number of {011} Type II HPVs:", len(hpvs))

# Check habit plane orientations
for hpv in hpvs:
    print(f"Habit plane: {hpv['hkl_hpp']}")
    print(f"Shear direction: {hpv['uvw_hpv']}")
    print(f"Transformation strain: {hpv['lbd']:.4f}")
```

### 5. Miller Indices Generation

#### `generate_hkls(hklmax, phase, hkls=[])`
**Purpose**: Generate unique Miller indices with symmetry

**Usage**:
```python
# Generate {hkl} up to h,k,l = 2
hkls_data = gp.generate_hkls(hklmax=2, phase='A')

# Unpack results
unique_hkls = hkls_data[0]        # List of unique (h,k,l) tuples
equiv_dict = hkls_data[1]         # Dictionary of equivalent indices
families = hkls_data[2]           # Family information
all_hkls = hkls_data[3]           # All symmetry-equivalent indices

print("Unique planes:", unique_hkls[:5])
print("Total unique:", len(unique_hkls))

# Get family multiplicity
for family, count in families.items():
    print(f"{{{family[0]}{family[1]}{family[2]}}}: {count} variants")
```

#### `generate_uvws(uvwmax, phase, uvws=[])`
**Purpose**: Generate unique direction indices with symmetry

**Usage**:
```python
# Generate <uvw> up to u,v,w = 2
uvws_data = gp.generate_uvws(uvwmax=2, phase='M')

unique_uvws = uvws_data[0]
print("Unique directions:", unique_uvws[:5])

# Example for slip system analysis
slip_dirs = gp.generate_uvws(uvwmax=1, phase='M')
print("<111> and <100> type directions:", slip_dirs[0])
```

### 6. Equivalent Directions/Planes with Projections

#### `getEqHKL(phasename, hkls, R2Proj=np.eye(3), hemisphere='upper')`
**Purpose**: Get equivalent HKL planes with stereographic projections

**Usage**:
```python
# Get equivalent {100} planes with projections
hkls_to_check = [[1,0,0], [0,1,0], [0,0,1]]
eq_hkls = gp.getEqHKL('A', hkls_to_check)

for eq in eq_hkls:
    print(f"Plane: {eq['hkl']}")
    print(f"Normal vector: {eq['vector']}")
    print(f"Equal-area projection: {eq['equalarea'][:2, 0]}")
    print(f"Stereographic projection: {eq['stereo'][:2, 0]}")
    print("---")

# Use in pole figure
from plotlib import plotter
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

for eq in eq_hkls:
    proj = eq['equalarea']
    p.ax.plot(proj[0], proj[1], 'ro', markersize=10)
    p.ax.text(proj[0, 0], proj[1, 0], eq['label'])
```

#### `getEqUVW(phasename, uvws, R2Proj=np.eye(3), hemisphere='upper')`
**Purpose**: Get equivalent UVW directions with projections

**Usage**:
```python
# Get equivalent <111> directions
uvws_to_check = [[1,1,1], [-1,1,1], [1,-1,1], [1,1,-1]]
eq_uvws = gp.getEqUVW('A', uvws_to_check)

print(f"Total equivalent directions: {len(eq_uvws)}")

# Plot in inverse pole figure
from projlib import stereotriangle
fig, ax = stereotriangle(equalarea=True)

for eq in eq_uvws:
    proj = eq['equalarea']
    ax.plot(proj[0], proj[1], 'bs', markersize=8)
```

### 7. Slip System Generation

#### `genSlipSys(phasename, SlipHKL, SlipUVW, mag=1)`
**Purpose**: Generate slip system from plane and direction

**Usage**:
```python
# Define slip system: {110}<111>
slip_plane = [1, 1, 0]
slip_dir = [1, 1, 1]

slip_sys = gp.genSlipSys('A', slip_plane, slip_dir)

print("Slip plane normal:", slip_sys['n'])
print("Slip direction:", slip_sys['m'])
print("Schmid tensor:\n", slip_sys['P'])

# Calculate Schmid factor
stress = np.array([[100, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0]])  # Uniaxial tension in x
schmid_factor = np.tensordot(slip_sys['P'], stress)
print(f"Schmid factor: {schmid_factor:.4f}")
```

#### `genEqSlipSys(phasename, SlipHKLs, SlipUVWs, getIndependentOnly=True, tol=1e-10, mag=1)`
**Purpose**: Generate all equivalent slip systems

**Usage**:
```python
# Define slip system family
SlipHKLs = [[1,1,0], [-1,1,0], [1,-1,0]]  # {110} planes
SlipUVWs = [[1,1,1]]                       # <111> directions

# Generate all equivalent systems
slip_systems = gp.genEqSlipSys('A', SlipHKLs, SlipUVWs, 
                                getIndependentOnly=True)

print(f"Number of independent slip systems: {len(slip_systems)}")

# For FCC: {111}<110> has 12 systems
SlipHKLs_fcc = [[1,1,1]]
SlipUVWs_fcc = [[1,1,0]]
slip_fcc = gp.genEqSlipSys('A', SlipHKLs_fcc, SlipUVWs_fcc)
print(f"FCC {{111}}<110> systems: {len(slip_fcc)}")  # Should be 12
```

### 8. Visualization

#### `plotHBVII011(name='NiTi', ProjType='stereo')`
**Purpose**: Plot {011} Type II twin habit plane variants

**Usage**:
```python
import matplotlib.pyplot as plt

# Plot on stereographic projection
gp.plotHBVII011(name='NiTi', ProjType='stereo')
plt.title('Type II {011} Habit Plane Variants - Stereographic')
plt.show()

# Plot on equal-area projection
gp.plotHBVII011(name='NiTi', ProjType='equalarea')
plt.title('Type II {011} Habit Plane Variants - Equal-Area')
plt.show()
```

#### `printHBVII011(name='NiTi')`
**Purpose**: Print habit plane variant information

**Usage**:
```python
# Print detailed HPV information
gp.printHBVII011(name='NiTi')

# Output shows:
# - Habit plane Miller indices
# - Shear direction
# - Transformation strain magnitude
# - Correspondence variant indices
```

### 9. Correspondence Variants

#### `getCVdict(name='NiTi')`
**Purpose**: Get correspondence variant dictionary

**Usage**:
```python
# Get correspondence variants
gp.getCVdict(name='NiTi')

# Access CV data
cv_dict = gp.CV['NiTi']
print(f"Number of correspondence variants: {len(cv_dict)}")

# Each CV contains lattice correspondence information
for cv_idx, cv in cv_dict.items():
    print(f"CV {cv_idx}:")
    print(f"  Cd matrix:\n{cv['Cd']}")
    print(f"  Determinant: {np.linalg.det(cv['Cd']):.6f}")
```

### 10. Elastic Properties

#### `getEl(name='NiTi')`
**Purpose**: Get elastic constants for phases

**Usage**:
```python
# Load elastic constants
gp.getEl(name='NiTi')

# Access austenite stiffness matrix (Voigt notation)
C_aust = gp.el['NiTi']['A']['C']
print("Austenite stiffness matrix (GPa):\n", C_aust / 1e3)

# Access martensite compliance
S_mart = gp.el['NiTi']['M']['S']
print("Martensite compliance:\n", S_mart)

# Tensor forms (for more complex calculations)
C_tensor_aust = gp.el['NiTi']['A']['CT']  # 3×3×3×3 tensor
S_tensor_mart = gp.el['NiTi']['M']['ST']  # 3×3×3×3 tensor

# Calculate Young's modulus in direction [100]
direction = np.array([1, 0, 0])
E_100 = 1.0 / (direction @ S_mart[:3, :3] @ direction)
print(f"Young's modulus E[100]: {E_100/1e3:.1f} GPa")
```

## Complete Workflow Examples

### Example 1: NiTi Shape Memory Alloy Analysis

```python
from getphases import getPhases
import numpy as np
import matplotlib.pyplot as plt

# Initialize
gp = getPhases()

# Define NiTi phases
gp.fromLatticeParams('A', a=3.015)  # B2 austenite
gp.fromLatticeParams('M', a=2.889, b=4.120, c=4.622,
                      beta=np.radians(96.8))  # B19' martensite

# Get orientation relationship
gp.getOR(name='NiTi', OR='NiTi')

# Calculate deformation gradients
gp.getDefGrad(name='NiTi')

# Get elastic properties
gp.getEl(name='NiTi')

# Calculate habit plane variants
gp.getHBVs(name='NiTi')
gp.getHPVII011(name='NiTi')

# Analyze transformation
print(f"Number of correspondence variants: {len(gp.defGrad['NiTi'])}")
print(f"Number of Type I twins: {len(gp.twinSys['NiTi']['I'])}")
print(f"Number of Type II twins: {len(gp.twinSys['NiTi']['II'])}")

# Visualize
gp.plotHBVII011(name='NiTi', ProjType='equalarea')
plt.title('NiTi {011} Type II Habit Plane Variants')
plt.show()
```

### Example 2: Transformation Strain Analysis

```python
from orilib import np_eulers_matrices
import matplotlib.pyplot as plt

# Setup phases
gp = getPhases()
gp.fromLatticeParams('A', a=3.015)
gp.fromLatticeParams('M', a=2.889, b=4.120, c=4.622, beta=np.radians(96.8))
gp.getOR(name='NiTi')
gp.getDefGrad(name='NiTi')

# Generate random austenite orientations
N = 1000
euler = np.random.rand(N, 3) * [360, 180, 360]
oris = np_eulers_matrices(euler, deg=True)

# Calculate transformation strains
strains = gp.getTrStrain(phase='A', name='NiTi', oris=oris)

# Analyze strain components
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
labels = ['ε₁₁', 'ε₂₂', 'ε₃₃', 'ε₂₃', 'ε₁₃', 'ε₁₂']

for i, (ax, label) in enumerate(zip(axes.flat, labels)):
    ax.hist(strains[:, i], bins=50, alpha=0.7, edgecolor='black')
    ax.set_xlabel(label)
    ax.set_ylabel('Frequency')
    ax.set_title(f'Distribution of {label}')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

print(f"Mean ε11: {np.mean(strains[:, 0]):.4f}")
print(f"Std ε11: {np.std(strains[:, 0]):.4f}")
```

### Example 3: Slip System Analysis

```python
# Setup
gp = getPhases()
gp.fromLatticeParams('A', a=3.015)  # Cubic

# Define slip systems for B2 (CsCl structure)
# Typical: {110}<001> and {100}<001>
SlipHKLs = [[1,1,0], [1,-1,0], [0,1,1], [0,1,-1], [1,0,1], [-1,0,1]]
SlipUVWs = [[0,0,1]]

# Generate all equivalent slip systems
slip_systems = gp.genEqSlipSys('A', SlipHKLs, SlipUVWs)

print(f"Total slip systems: {len(slip_systems)}")

# Calculate Schmid factors for uniaxial tension
tension_axis = np.array([1, 0, 0])  # Tension along [100]
stress = np.outer(tension_axis, tension_axis)  # σ = n ⊗ n

schmid_factors = []
for ss in slip_systems:
    P = ss['P']  # Schmid tensor
    schmid = np.tensordot(P, stress)
    schmid_factors.append(abs(schmid))

# Find maximum Schmid factor
max_schmid = max(schmid_factors)
max_idx = schmid_factors.index(max_schmid)

print(f"Maximum Schmid factor: {max_schmid:.4f}")
print(f"Active slip system:")
print(f"  Plane: {slip_systems[max_idx]['hkl']}")
print(f"  Direction: {slip_systems[max_idx]['uvw']}")

# Plot Schmid factor distribution
plt.figure(figsize=(8, 6))
plt.hist(schmid_factors, bins=20, edgecolor='black')
plt.xlabel('Schmid Factor')
plt.ylabel('Number of Systems')
plt.title('Schmid Factor Distribution for [100] Tension')
plt.axvline(max_schmid, color='r', linestyle='--', label=f'Max = {max_schmid:.3f}')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

### Example 4: Pole Figure Generation

```python
from plotlib import plotter

# Setup
gp = getPhases()
gp.fromLatticeParams('M', a=2.889, b=4.120, c=4.622, beta=np.radians(96.8))

# Generate {100} and {010} poles
hkls = [[1,0,0], [0,1,0], [0,0,1]]
eq_hkls = gp.getEqHKL('M', hkls)

# Create pole figure
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

# Plot each family in different color
colors = ['red', 'blue', 'green']
markers = ['o', 's', '^']

for i, hkl in enumerate(hkls):
    # Filter equivalents for this family
    family_poles = [eq for eq in eq_hkls if list(eq['hklf']) == hkl]
    
    for eq in family_poles:
        proj = eq['equalarea']
        p.ax.plot(proj[0], proj[1], markers[i], 
                 color=colors[i], markersize=10, 
                 label=f'{{{hkl[0]}{hkl[1]}{hkl[2]}}}' if eq == family_poles[0] else '')

p.ax.legend()
p.ax.set_title('Martensite {100}, {010}, {001} Pole Figure')
plt.show()
```

## Function Parameters Quick Reference

| Method | Key Parameters | Returns |
|--------|---------------|---------|
| `fromCif` | file, phasename | None (populates phases dict) |
| `fromLatticeParams` | phasename, a, b, c, α, β, γ | None (populates phases dict) |
| `getOR` | name, OR | None (populates OR dict) |
| `getDefGrad` | name, OR | None (populates defGrad dict) |
| `getTrStrain` | phase, name, oris | strain array (N, 6) |
| `getHBVs` | name | None (populates vardict) |
| `generate_hkls` | hklmax, phase, hkls | tuple (unique, equiv, families, all) |
| `generate_uvws` | uvwmax, phase, uvws | tuple (unique, equiv, families, all) |
| `getEqHKL` | phasename, hkls, R2Proj, hemisphere | list of dicts |
| `getEqUVW` | phasename, uvws, R2Proj, hemisphere | list of dicts |
| `genSlipSys` | phasename, SlipHKL, SlipUVW, mag | dict with n, m, P |
| `genEqSlipSys` | phasename, SlipHKLs, SlipUVWs, flags | list of slip systems |
| `plotHBVII011` | name, ProjType | None (creates plot) |

## Mathematical Formulas

### Deformation Gradient
```
F = L_product @ L_parent^(-1)

where:
  L_product: Product phase lattice matrix
  L_parent: Parent phase lattice matrix
```

### Transformation Strain
```
ε = (F^T @ F - I) / 2

where:
  F: Deformation gradient
  I: Identity matrix
```

### Schmid Tensor
```
P = m ⊗ n

where:
  m: Slip direction (unit vector)
  n: Slip plane normal (unit vector)
  ⊗: Outer product
```

### Schmid Factor
```
μ = P : σ = Σᵢⱼ Pᵢⱼ σᵢⱼ

where:
  σ: Applied stress tensor
  :: Double contraction
```

## Common Pitfalls

1. **Angle Units**: Always use radians for angles (α, β, γ)
   ```python
   # ✗ Wrong
   gp.fromLatticeParams('M', a=2.889, b=4.120, c=4.622, beta=96.8)
   
   # ✓ Correct
   gp.fromLatticeParams('M', a=2.889, b=4.120, c=4.622, beta=np.radians(96.8))
   ```

2. **Phase Names**: Use consistent phase naming
   ```python
   # Define clear names
   gp.setAttributes(austenite='Austenite', martensite='Martensite')
   # Or stick with defaults 'A' and 'M'
   ```

3. **OR Calculation Order**: Get OR before calculating deformation gradients
   ```python
   # ✓ Correct order
   gp.getOR(name='NiTi')
   gp.getDefGrad(name='NiTi')  # Requires OR
   ```

4. **Matrix Dimensions**: Check lattice matrix column format
   ```python
   # Lattice vectors are COLUMNS
   L = gp.phases['A']['L']
   a1 = L[:, 0]  # First lattice vector
   a2 = L[:, 1]  # Second lattice vector
   a3 = L[:, 2]  # Third lattice vector
   ```

## Integration with Other Modules

### With orilib
```python
from orilib import np_eulers_matrices

# Generate orientations
euler = np.random.rand(100, 3) * [360, 180, 360]
oris = np_eulers_matrices(euler, deg=True)

# Calculate transformation strains
strains = gp.getTrStrain(phase='A', name='NiTi', oris=oris)
```

### With projlib
```python
from projlib import equalarea_directions, stereotriangle

# Get equivalent directions
eq_uvws = gp.getEqUVW('A', [[1,1,1]])

# Project and plot
fig, ax = stereotriangle(equalarea=True)
for eq in eq_uvws:
    proj = eq['equalarea']
    ax.plot(proj[0], proj[1], 'ro')
```

### With plotlib
```python
from plotlib import plotter

# Create pole figure
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

# Add phase-specific poles
eq_hkls = gp.getEqHKL('M', [[1,0,0]])
for eq in eq_hkls:
    proj = eq['equalarea']
    p.ax.plot(proj[0], proj[1], 'bs', markersize=10)
```

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| "KeyError: 'A'" | Phase not defined | Call `fromCif()` or `fromLatticeParams()` first |
| "KeyError: 'NiTi' in OR" | OR not calculated | Call `getOR(name='NiTi')` |
| Wrong β angle | Degrees used instead of radians | Use `np.radians()` for angles |
| Empty defGrad dict | OR not set | Call `getOR()` before `getDefGrad()` |
| Asymmetric lattice matrix | Wrong parameter order | Check a, b, c, α, β, γ assignments |

## Performance Tips

1. **Pre-calculate symmetry operations** once and reuse
2. **Cache generated Miller indices** if used repeatedly
3. **Use vectorized orientation arrays** for batch strain calculations
4. **Store phase definitions** in configuration files

## Best Practices

1. Always define phases before orientation relationships
2. Use descriptive phase names ('Austenite', 'Martensite') for clarity
3. Document lattice parameters with units (Ã… for lengths, radians for angles)
4. Validate lattice matrices (check orthogonality for cubic, angles for monoclinic)
5. Store elastic constants with references to source papers
6. Use consistent coordinate systems throughout analysis

---

*This quick reference covers the most commonly used features of getphases.py. For complete documentation, see getphases_documentation_summary.md*

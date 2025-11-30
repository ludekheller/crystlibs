# GetPhases.py Documentation Summary

## Module Overview

**getphases.py** is a comprehensive Python module for managing crystallographic phases and their transformations, with specific emphasis on shape memory alloys (particularly NiTi). The module provides tools for phase definition, orientation relationship calculation, deformation gradient computation, twin system analysis, and habit plane variant determination.

### Purpose

The module serves as the central framework for:
- Defining crystal structures from CIF files or lattice parameters
- Computing lattice correspondence matrices between phases
- Calculating stress-free transformation deformation gradients
- Analyzing twin systems (Type I and Type II)
- Determining habit plane variants
- Generating crystallographic indices with symmetry operations
- Computing transformation strains for orientation distributions
- Visualizing phase transformations on pole figures

### Key Applications

1. **Shape Memory Alloys**: NiTi, CuAlNi, FeMnSi analysis
2. **Phase Transformation Crystallography**: Martensitic transformations
3. **Twin Analysis**: Deformation and transformation twins
4. **EBSD Data Analysis**: Multiphase material characterization
5. **Texture Analysis**: Transformed microstructures
6. **Habit Plane Prediction**: Validation against experimental observations

---

## Class: getPhases

The main class providing phase management functionality.

---

## Initialization Methods

### `__init__()`

**Purpose**: Initialize the phase management system

**Input**:
- None

**Output**:
- None (creates instance with empty dictionaries)

**Attributes Initialized**:
```python
self.phases = {}       # Phase definitions
self.OR = {}          # Orientation relationships
self.vardict = {}     # Variant information
self.austenite = 'A'  # Austenite phase identifier
self.martensite = 'M' # Martensite phase identifier
self.defGrad = {}     # Deformation gradients
self.CV = {}          # Correspondence variants
self.twinSys = {}     # Twin systems
self.el = {}          # Elastic properties
```

**Usage Example**:
```python
from getphases import getPhases

# Create instance
gp = getPhases()

# Verify initialization
print("Phases:", gp.phases)
print("Austenite label:", gp.austenite)
print("Martensite label:", gp.martensite)
```

**Implementation Details**:
- Automatically calls `getOR()` during initialization
- Creates empty dictionaries for all phase-related data
- Sets default phase identifiers ('A' for austenite, 'M' for martensite)

---

### `setAttributes(**kwargs)`

**Purpose**: Dynamically set arbitrary attributes

**Input**:
- `**kwargs`: keyword arguments - Any attribute-value pairs

**Output**:
- None (updates instance attributes)

**Usage Examples**:
```python
# Set custom phase names
gp.setAttributes(austenite='Austenite', martensite='Martensite')

# Set analysis parameters
gp.setAttributes(
    temperature=300,  # K
    pressure=1.0,     # atm
    analysis_name='NiTi_Study_001'
)

# Set multiple attributes
gp.setAttributes(
    verbose=True,
    output_path='/results',
    plot_format='png'
)

# Verify
print("Temperature:", gp.temperature)
```

---

## Phase Definition Methods

### `fromCif(file, phasename)`

**Purpose**: Load crystal structure from Crystallographic Information File

**Input**:
- `file`: str - Path to CIF file
- `phasename`: str - Name for this phase (e.g., 'A', 'M', 'austenite')

**Output**:
- None (populates `self.phases[phasename]` dictionary)

**Phase Dictionary Structure**:
```python
self.phases[phasename] = {
    'cif': Crystal object,           # from crystals library
    'ciffile': str,                  # Path to CIF file
    'symops': list,                  # Direct symmetry operations (3×3)
    'recsymops': list,               # Reciprocal symmetry operations (3×3)
    'a': float,                      # Lattice parameter (Å)
    'b': float,                      # Lattice parameter (Å)
    'c': float,                      # Lattice parameter (Å)
    'alpha': float,                  # Angle (radians)
    'beta': float,                   # Angle (radians)
    'gamma': float,                  # Angle (radians)
    'L': ndarray (3,3),             # Direct lattice matrix
    'Lr': ndarray (3,3),            # Reciprocal lattice matrix
    'G': ndarray (3,3),             # Direct metric tensor
    'Gr': ndarray (3,3),            # Reciprocal metric tensor
    'Gi': ndarray (3,3),            # Inverse direct metric
    'Gir': ndarray (3,3)            # Inverse reciprocal metric
}
```

**Usage Examples**:
```python
# Load B2 austenite from CIF
gp.fromCif('B2_NiTi.cif', 'A')

# Load B19' martensite
gp.fromCif('B19prime_NiTi.cif', 'M')

# Access loaded data
print("Space group:", gp.phases['A']['cif'].spacegroup)
print("Lattice a:", gp.phases['A']['a'], "Å")
print("Number of symmetry ops:", len(gp.phases['A']['symops']))

# Access lattice matrices
L_austenite = gp.phases['A']['L']
print("Austenite lattice matrix:\n", L_austenite)
```

**Notes**:
- Requires `crystals` library for CIF parsing
- Automatically converts angles from degrees to radians
- Automatically calls `fromLatticeParams()` to build matrices
- Extracts both direct and reciprocal space symmetry operations

---

### `fromLatticeParams(phasename, a, b=None, c=None, alpha=π/2, beta=π/2, gamma=π/2)`

**Purpose**: Construct lattice matrices from lattice parameters

**Input**:
- `phasename`: str - Name for this phase
- `a`: float - Lattice parameter a (Å)
- `b`: float - Lattice parameter b (Å), default: None (uses a for cubic)
- `c`: float - Lattice parameter c (Å), default: None (uses a for cubic)
- `alpha`: float - Angle between b and c (radians), default: π/2
- `beta`: float - Angle between a and c (radians), default: π/2
- `gamma`: float - Angle between a and b (radians), default: π/2

**Output**:
- None (populates/updates `self.phases[phasename]` dictionary)

**Mathematical Formulas**:

Direct lattice matrix construction:
```
L[:,0] = [a, 0, 0]
L[:,1] = [b·cos(γ), b·sin(γ), 0]
L[:,2] = [cx, cy, cz]

where:
  cx = c·cos(β)
  cy = c·(cos(α) - cos(β)·cos(γ))/sin(γ)
  cz = sqrt(c² - cx² - cy²)
```

Reciprocal lattice matrix:
```
Lr = L^(-T) = (L^(-1))^T
```

Metric tensors:
```
G = L^T @ L              (direct space)
Gr = Lr^T @ Lr           (reciprocal space)
```

**Usage Examples**:
```python
import numpy as np

# Cubic B2 austenite (a = 3.015 Å)
gp.fromLatticeParams('A', a=3.015)

# Monoclinic B19' martensite
gp.fromLatticeParams('M', 
    a=2.889,              # Å
    b=4.120,              # Å
    c=4.622,              # Å
    beta=np.radians(96.8) # Convert degrees to radians
)

# Hexagonal (for comparison with Ti)
gp.fromLatticeParams('Ti_alpha', 
    a=2.95, 
    c=4.68,
    gamma=np.radians(120)  # 120° for hexagonal
)

# Orthorhombic
gp.fromLatticeParams('Ortho', 
    a=3.0, 
    b=4.0, 
    c=5.0
    # All angles default to 90°
)

# Triclinic (most general)
gp.fromLatticeParams('Triclinic',
    a=3.0, b=4.0, c=5.0,
    alpha=np.radians(85),
    beta=np.radians(95),
    gamma=np.radians(100)
)

# Access constructed matrices
L = gp.phases['M']['L']
Lr = gp.phases['M']['Lr']
G = gp.phases['M']['G']

# Verify reciprocal lattice
print("L @ Lr.T =\n", L @ Lr.T)  # Should be identity

# Calculate unit cell volume
V = np.abs(np.linalg.det(L))
print(f"Unit cell volume: {V:.4f} ų")

# Check angle
a_vec = L[:, 0]
c_vec = L[:, 2]
beta_calc = np.arccos(np.dot(a_vec, c_vec) / 
                      (np.linalg.norm(a_vec) * np.linalg.norm(c_vec)))
print(f"Beta angle: {np.degrees(beta_calc):.2f}°")
```

**Crystal System Support**:

| System | Parameters | Example |
|--------|-----------|---------|
| Cubic | a | B2 NiTi austenite |
| Tetragonal | a, c | - |
| Orthorhombic | a, b, c | - |
| Hexagonal | a, c, γ=120° | Ti (α-phase) |
| Monoclinic | a, b, c, β | B19' NiTi martensite |
| Triclinic | a, b, c, α, β, γ | General case |

**Notes**:
- All angles MUST be in radians
- Lattice matrix L has lattice vectors as COLUMNS
- Default angles are π/2 (90°) for all
- If b or c not specified, uses a (for cubic/tetragonal)
- Automatically computes inverse metrics for crystallographic calculations

---

## Orientation Relationship Methods

### `getOR(name=None, OR='NiTi')`

**Purpose**: Get orientation relationship matrices between phases

**Input**:
- `name`: str - Name for storing this OR (default: None, uses OR value)
- `OR`: str - Type of orientation relationship (default: 'NiTi')

**Output**:
- None (populates `self.OR[name]` dictionary)

**OR Dictionary Structure**:
```python
self.OR[name] = {
    'Cd': ndarray (3,3),   # Direction correspondence matrix (M→A)
    'CId': ndarray (3,3),  # Inverse direction correspondence (A→M)
    'Cp': ndarray (3,3),   # Plane correspondence matrix (M→A)
    'CIp': ndarray (3,3)   # Inverse plane correspondence (A→M)
}
```

**Mathematical Relationships**:
```
Martensite direction → Austenite direction:
  u_A = Cd @ u_M

Austenite direction → Martensite direction:
  u_M = CId @ u_A = Cd^(-1) @ u_A

Plane normals (similar):
  n_A = Cp @ n_M
  n_M = CIp @ n_A
```

**Usage Examples**:
```python
# Get NiTi orientation relationship
gp.getOR(name='NiTi', OR='NiTi')

# Access correspondence matrices
Cd = gp.OR['NiTi']['Cd']
CId = gp.OR['NiTi']['CId']
Cp = gp.OR['NiTi']['Cp']
CIp = gp.OR['NiTi']['CIp']

# Example: Map martensite [100]_M to austenite frame
u_M = np.array([1, 0, 0])
u_A = Cd @ u_M
print(f"Martensite [100] → Austenite {u_A}")

# Example: Map austenite (110)_A to martensite frame
n_A = np.array([1, 1, 0])
n_M = CIp @ n_A
print(f"Austenite (110) → Martensite {n_M}")

# Verify inverse relationship
u_M_recovered = CId @ u_A
print("Recovered:", u_M_recovered)
print("Match:", np.allclose(u_M, u_M_recovered))

# Example with multiple directions
directions_M = np.array([[1,0,0], [0,1,0], [0,0,1]]).T
directions_A = Cd @ directions_M
print("Martensite basis in austenite frame:\n", directions_A)
```

**Current Implementations**:
- 'NiTi': B19' ↔ B2 lattice correspondence (Waitz notation)

**Notes**:
- Must call this before `getDefGrad()`
- Correspondence matrices relate crystallographic directions/planes
- Uses literature-based correspondence for NiTi
- Can be extended for other material systems

---

## Deformation Gradient Methods

### `getDefGrad(name=None, OR='NiTi')`

**Purpose**: Compute deformation gradients for stress-free phase transformation

**Input**:
- `name`: str - Name for this set of deformation gradients
- `OR`: str - Orientation relationship to use (default: 'NiTi')

**Output**:
- None (populates `self.defGrad[name]` list)

**Deformation Gradient List Structure**:
```python
self.defGrad[name] = [
    {
        'F': ndarray (3,3),      # Deformation gradient matrix
        'CVid': int,             # Correspondence variant ID
        'Cd': ndarray (3,3),     # Direction correspondence for this CV
        'Cp': ndarray (3,3),     # Plane correspondence for this CV
        # ... additional variant information
    },
    # ... 12 total correspondence variants for NiTi
]
```

**Mathematical Formula**:
```
F = L_product @ Cd @ L_parent^(-1)

where:
  L_product: Product phase lattice matrix
  L_parent: Parent phase lattice matrix
  Cd: Direction correspondence matrix for this variant
```

**Properties**:
- Volume change: det(F) ≈ 1.01 for NiTi (slight expansion)
- Transformation strain: ε = (F^T @ F - I) / 2
- 12 correspondence variants for NiTi (from 3 lattice correspondences × 4 B2 domains)

**Usage Examples**:
```python
# Calculate deformation gradients
gp.getDefGrad(name='NiTi', OR='NiTi')

# Access first variant
F1 = gp.defGrad['NiTi'][0]['F']
print("Deformation gradient F1:\n", F1)
print("Determinant:", np.linalg.det(F1))

# Calculate transformation strain for variant 1
epsilon = (F1.T @ F1 - np.eye(3)) / 2
print("Transformation strain ε:\n", epsilon)

# Get principal strains
eigenvalues, eigenvectors = np.linalg.eig(epsilon)
print("Principal strains:", np.sort(eigenvalues))

# Analyze all variants
print(f"Total correspondence variants: {len(gp.defGrad['NiTi'])}")

for i, variant in enumerate(gp.defGrad['NiTi']):
    F = variant['F']
    det_F = np.linalg.det(F)
    print(f"Variant {i}: det(F) = {det_F:.6f}, CV ID = {variant['CVid']}")

# Calculate transformation strain distribution
strains = []
for variant in gp.defGrad['NiTi']:
    F = variant['F']
    eps = (F.T @ F - np.eye(3)) / 2
    strains.append(eps)

# Stack for analysis
strains = np.array(strains)
print("Strain shape:", strains.shape)  # (12, 3, 3)
```

**Notes**:
- Requires prior call to `getOR(name, OR)`
- Generates all correspondence variants
- Each variant corresponds to a specific lattice correspondence
- Deformation gradients are stress-free (no applied stress)

---

### `getTrStrain(phase=None, name='NiTi', oris=None)`

**Purpose**: Calculate transformation strain for orientation distribution

**Input**:
- `phase`: str - Parent phase identifier (default: None, uses self.austenite)
- `name`: str - Name of OR/deformation gradient set to use
- `oris`: ndarray (N, 3, 3) - Orientation matrices (parent→sample frame)

**Output**:
- `strains`: ndarray (N, 6) - Transformation strains in Voigt notation
  - Components: [ε₁₁, ε₂₂, ε₃₃, ε₂₃, ε₁₃, ε₁₂]

**Mathematical Formula**:
```
For each orientation g and each correspondence variant F_CV:
  F_total = g @ F_CV @ g^T
  ε_sample = (F_total^T @ F_total - I) / 2
  
Convert ε_sample (3×3 tensor) to Voigt notation (6 components)
```

**Usage Examples**:
```python
from orilib import np_eulers_matrices
import matplotlib.pyplot as plt

# Generate random orientations
N = 1000
euler_angles = np.random.rand(N, 3) * [360, 180, 360]
oris = np_eulers_matrices(euler_angles, deg=True)

# Calculate transformation strains
strains = gp.getTrStrain(phase='A', name='NiTi', oris=oris)

print("Strains shape:", strains.shape)  # (1000, 6)
print("Strain components: ε11, ε22, ε33, ε23, ε13, ε12")

# Analyze ε11 distribution
plt.figure(figsize=(10, 6))
plt.hist(strains[:, 0], bins=50, alpha=0.7, edgecolor='black')
plt.xlabel('ε₁₁')
plt.ylabel('Frequency')
plt.title('Distribution of Normal Strain ε₁₁')
plt.grid(True, alpha=0.3)
plt.show()

# Statistics for all components
labels = ['ε₁₁', 'ε₂₂', 'ε₃₃', 'ε₂₃', 'ε₁₃', 'ε₁₂']
for i, label in enumerate(labels):
    mean = np.mean(strains[:, i])
    std = np.std(strains[:, i])
    print(f"{label}: mean = {mean:.4f}, std = {std:.4f}")

# Find orientation with maximum ε11
max_idx = np.argmax(strains[:, 0])
print(f"Maximum ε11: {strains[max_idx, 0]:.4f}")
print(f"Orientation (Euler): {euler_angles[max_idx]}")

# Multiplot of all components
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
for i, (ax, label) in enumerate(zip(axes.flat, labels)):
    ax.hist(strains[:, i], bins=30, alpha=0.7, edgecolor='black')
    ax.set_xlabel(label)
    ax.set_ylabel('Frequency')
    ax.set_title(f'Distribution of {label}')
    ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

**Notes**:
- Calculates strains for ALL correspondence variants
- Returns strain in sample reference frame
- Voigt notation ordering: [11, 22, 33, 23, 13, 12]
- Useful for texture-transformation strain coupling analysis

---

## Twin System Methods

### `getTwinSys(name=None, OR='NiTi')`

**Purpose**: Calculate twin systems (Type I and Type II)

**Input**:
- `name`: str - Name for this twin system set
- `OR`: str - Orientation relationship to use

**Output**:
- None (populates `self.twinSys[name]` dictionary)

**Twin System Dictionary Structure**:
```python
self.twinSys[name] = {
    'I': list,   # Type I twin systems
    'II': list   # Type II twin systems
}
```

**Twin Types**:
- **Type I**: Shear plane is rational in both parent and product (K₁ = K₂)
- **Type II**: Shear direction is rational in both phases (η₁ = η₂)

**Usage Examples**:
```python
# Calculate twin systems
gp.getTwinSys(name='NiTi', OR='NiTi')

# Access Type I twins
type_I_twins = gp.twinSys['NiTi']['I']
print(f"Number of Type I twin systems: {len(type_I_twins)}")

# Access Type II twins
type_II_twins = gp.twinSys['NiTi']['II']
print(f"Number of Type II twin systems: {len(type_II_twins)}")

# Analyze a Type I twin
if type_I_twins:
    twin = type_I_twins[0]
    print("Type I Twin System 0:")
    for key in twin.keys():
        print(f"  {key}: {twin[key]}")

# Typical output includes:
# - Twin plane (K1, K2)
# - Shear direction (η1, η2)
# - Transformation matrices
# - Correspondence variant indices
```

---

## Habit Plane Variant Methods

### `getHBVs(name='NiTi')`

**Purpose**: Calculate habit plane variants

**Input**:
- `name`: str - Name for this HPV set

**Output**:
- None (populates `self.vardict[name]['HBV']` list)

**Usage Examples**:
```python
# Get habit plane variants
gp.getHBVs(name='NiTi')

# Access HPV data
hbvs = gp.vardict['NiTi']['HBV']
print(f"Number of habit plane variants: {len(hbvs)}")

# Analyze first HPV
hbv = hbvs[0]
print("Habit plane normal:", hbv['n_hpv'])
print("Shear direction:", hbv['m_hpv'])
print("Transformation strain magnitude:", hbv.get('lambda', 'N/A'))
```

---

### `getHPVII011(name='NiTi')`

**Purpose**: Get {011} Type II twin habit plane variants

**Input**:
- `name`: str - Name for this HPV set

**Output**:
- None (populates `self.vardict[name]['HPVII011']` list)

**Usage Examples**:
```python
# Calculate {011} Type II twin HPVs
gp.getHPVII011(name='NiTi')

# Access data
hpvs_011 = gp.vardict['NiTi']['HPVII011']
print(f"Number of {{011}} Type II HPVs: {len(hpvs_011)}")

# Detailed analysis
for i, hpv in enumerate(hpvs_011):
    print(f"\nHPV {i}:")
    print(f"  Habit plane (hkl): {hpv.get('hkl_hpp', 'N/A')}")
    print(f"  Shear direction (uvw): {hpv.get('uvw_hpv', 'N/A')}")
    print(f"  Lambda: {hpv.get('lbd', 'N/A')}")
```

---

## Continued in Part 2...

*Due to length, the documentation summary continues in a second file with remaining methods including Miller indices generation, visualization, slip systems, and advanced features.*

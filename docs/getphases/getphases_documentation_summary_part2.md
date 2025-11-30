# GetPhases.py Documentation Summary - Part 2

## Miller Indices Generation Methods

### `generate_hkls(hklmax, phase, hkls=[])`

**Purpose**: Generate unique Miller indices up to maximum value considering symmetry

**Input**:
- `hklmax`: int - Maximum Miller index value (e.g., 2 gives indices from -2 to +2)
- `phase`: str - Phase name to get symmetry operations from
- `hkls`: list - Existing HKL list to extend (default: empty list)

**Output**:
- tuple containing:
  - `[0]`: list of tuples - Unique (h,k,l) representatives
  - `[1]`: dict - Maps each unique hkl to list of symmetry equivalents
  - `[2]`: dict - Family information with multiplicities
  - `[3]`: ndarray - All symmetry-equivalent indices as array

**Usage Examples**:
```python
# Generate {hkl} up to index 2
hkls_data = gp.generate_hkls(hklmax=2, phase='A')

# Unpack results
unique_hkls = hkls_data[0]
equiv_dict = hkls_data[1]
families = hkls_data[2]
all_hkls = hkls_data[3]

print(f"Unique planes: {len(unique_hkls)}")
print("First 5 unique:", unique_hkls[:5])

# Check family multiplicities
for family, multiplicity in families.items():
    print(f"{{family}}: {multiplicity} variants")

# Example output for cubic:
# {100}: 6 variants  (±100, ±010, ±001)
# {110}: 12 variants
# {111}: 8 variants

# Get all equivalents of a specific plane
hkl_100 = (1.0, 0.0, 0.0)
if hkl_100 in equiv_dict:
    equivalents = equiv_dict[hkl_100]
    print(f"Equivalents of (100): {equivalents}")

# Use for diffraction analysis
for hkl in unique_hkls[:10]:
    h, k, l = hkl
    # Calculate d-spacing, structure factor, etc.
    Lr = gp.phases['A']['Lr']
    G_hkl = Lr @ np.array([h, k, l])
    d_spacing = 1.0 / np.linalg.norm(G_hkl)
    print(f"({h:.0f}{k:.0f}{l:.0f}): d = {d_spacing:.3f} Å")
```

**Notes**:
- Uses reciprocal space symmetry operations
- Returns both unique representatives and all equivalents
- Families are grouped by permutation symmetry
- Useful for systematic diffraction analysis

---

### `generate_uvws(uvwmax, phase, uvws=[])`

**Purpose**: Generate unique direction indices [uvw] with symmetry

**Input**:
- `uvwmax`: int - Maximum direction index
- `phase`: str - Phase name
- `uvws`: list - Existing list to extend

**Output**:
- tuple (same structure as `generate_hkls`)

**Usage Examples**:
```python
# Generate <uvw> up to index 2
uvws_data = gp.generate_uvws(uvwmax=2, phase='M')

unique_uvws = uvws_data[0]
print(f"Unique directions: {len(unique_uvws)}")

# For slip system analysis
slip_dirs = gp.generate_uvws(uvwmax=1, phase='A')
print("Low-index directions:", slip_dirs[0])

# Typical for BCC: <111>, <100>
# Typical for FCC: <110>, <100>

# Check direction families
families_uvw = uvws_data[2]
for family, count in families_uvw.items():
    print(f"<{family[0]}{family[1]}{family[2]}>: {count} directions")
```

**Notes**:
- Uses direct space symmetry operations (not reciprocal)
- Essential for slip system and twin analysis
- Provides foundation for texture analysis

---

## Projection and Visualization Methods

### `getEqHKL(phasename, hkls, R2Proj=np.eye(3), hemisphere='upper', **kwargs)`

**Purpose**: Get equivalent HKL planes with stereographic and equal-area projections

**Input**:
- `phasename`: str - Phase name
- `hkls`: list - List of [h,k,l] Miller indices
- `R2Proj`: ndarray (3,3) - Rotation to projection reference frame (default: identity)
- `hemisphere`: str - 'upper', 'lower', or 'both' (default: 'upper')
- `**kwargs`: Additional arguments

**Output**:
- list of dictionaries, each containing:
  - `'vector'`: ndarray (3,) - Plane normal in projection frame
  - `'hkl'`: list - Miller indices [h,k,l]
  - `'hklf'`: original Miller indices
  - `'label'`: str - Formatted label "(hkl)"
  - `'equalarea'`: ndarray (3, 1) - Equal-area projection coordinates
  - `'stereo'`: ndarray (3, 1) - Stereographic projection coordinates
  - `'equalarea plane'`: Plane trace for equal-area projection
  - `'stereo plane'`: Plane trace for stereographic projection

**Usage Examples**:
```python
from plotlib import plotter

# Get {100} family with projections
hkls_100 = [[1,0,0], [0,1,0], [0,0,1]]
eq_hkls = gp.getEqHKL('M', hkls_100)

print(f"Total equivalent planes: {len(eq_hkls)}")

# Create pole figure
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

# Plot each plane
for eq in eq_hkls:
    proj = eq['equalarea']
    p.ax.plot(proj[0], proj[1], 'ro', markersize=10)
    p.ax.text(proj[0, 0] + 0.05, proj[1, 0], eq['label'])

p.ax.set_title('Martensite {100} Pole Figure')
plt.show()

# Access plane traces
for eq in eq_hkls[:3]:
    trace = eq['equalarea plane']
    p.ax.plot(trace[0], trace[1], 'b-', linewidth=1, alpha=0.5)

# With custom rotation
R_custom = np.array([[0, -1, 0],
                     [1,  0, 0],
                     [0,  0, 1]])  # 90° around z
eq_hkls_rotated = gp.getEqHKL('M', hkls_100, R2Proj=R_custom)
```

**Notes**:
- Automatically applies crystal symmetry operations
- Provides both projection coordinates and plane traces
- Essential for pole figure generation
- Coordinate system can be rotated via R2Proj

---

### `getEqUVW(phasename, uvws, R2Proj=np.eye(3), hemisphere='upper', **kwargs)`

**Purpose**: Get equivalent UVW directions with projections

**Input**:
- `phasename`: str - Phase name
- `uvws`: list - List of [u,v,w] direction indices
- `R2Proj`: ndarray (3,3) - Rotation matrix
- `hemisphere`: str - Projection hemisphere

**Output**:
- list of dictionaries (similar structure to `getEqHKL`)

**Usage Examples**:
```python
from projlib import stereotriangle

# Get <111> family
uvws_111 = [[1,1,1], [-1,1,1], [1,-1,1], [1,1,-1]]
eq_uvws = gp.getEqUVW('A', uvws_111)

print(f"Equivalent <111> directions: {len(eq_uvws)}")

# Plot in inverse pole figure
fig, ax = stereotriangle(equalarea=True)

for eq in eq_uvws:
    proj = eq['equalarea']
    ax.plot(proj[0], proj[1], 'bs', markersize=8)
    ax.text(proj[0, 0], proj[1, 0], eq['label'], fontsize=8)

ax.set_title('Austenite <111> Directions')
plt.show()

# Fiber texture analysis
# Get directions along [001] fiber
fiber_dirs = [[0,0,1]]
eq_fiber = gp.getEqUVW('A', fiber_dirs)

for eq in eq_fiber:
    print(f"Direction: {eq['uvw']}")
    print(f"Equal-area coords: {eq['equalarea'][:2, 0]}")
```

---

### `plotHBVII011(name='NiTi', ProjType='stereo')`

**Purpose**: Plot {011} Type II twin habit plane variants on pole figure

**Input**:
- `name`: str - Name of HPV set to plot
- `ProjType`: str - 'stereo' or 'equalarea' projection type

**Output**:
- None (creates matplotlib figure)

**Usage Examples**:
```python
import matplotlib.pyplot as plt

# Calculate HPVs first
gp.getHPVII011(name='NiTi')

# Plot on stereographic projection
gp.plotHBVII011(name='NiTi', ProjType='stereo')
plt.title('Type II {011} Habit Plane Variants - Stereographic')
plt.savefig('hpv_stereo.png', dpi=300)
plt.show()

# Plot on equal-area projection
gp.plotHBVII011(name='NiTi', ProjType='equalarea')
plt.title('Type II {011} Habit Plane Variants - Equal-Area')
plt.savefig('hpv_equalarea.png', dpi=300)
plt.show()

# Customize further
fig, ax = plt.subplots(figsize=(8, 8))
gp.plotHBVII011(name='NiTi', ProjType='equalarea')
ax.set_title('NiTi Habit Planes', fontsize=16)
ax.legend()
plt.tight_layout()
plt.show()
```

**Notes**:
- Automatically creates pole figure with proper projection
- Shows habit plane normals
- Useful for comparison with experimental EBSD data
- Can be customized with matplotlib commands

---

## Slip System Methods

### `genSlipSys(phasename, SlipHKL, SlipUVW, mag=1)`

**Purpose**: Generate slip system from plane and direction

**Input**:
- `phasename`: str - Phase name
- `SlipHKL`: list [h,k,l] - Slip plane Miller indices
- `SlipUVW`: list [u,v,w] - Slip direction indices
- `mag`: float - Magnitude scaling (default: 1)

**Output**:
- dict containing:
  - `'n'`: ndarray (3,) - Slip plane normal (unit vector)
  - `'m'`: ndarray (3,) - Slip direction (unit vector)
  - `'P'`: ndarray (3,3) - Schmid tensor (m ⊗ n)
  - `'hkl'`: list - Slip plane indices
  - `'uvw'`: list - Slip direction indices

**Mathematical Formulas**:
```
Slip plane normal: n = (Lr @ hkl) / |Lr @ hkl|
Slip direction: m = (L @ uvw) / |L @ uvw|
Schmid tensor: P = m ⊗ n
Schmid factor: μ = P : σ
```

**Usage Examples**:
```python
# Define {110}<111> slip system (BCC)
slip_plane = [1, 1, 0]
slip_dir = [1, 1, 1]

slip_sys = gp.genSlipSys('A', slip_plane, slip_dir)

print("Slip plane normal:", slip_sys['n'])
print("Slip direction:", slip_sys['m'])
print("Schmid tensor:\n", slip_sys['P'])

# Verify orthogonality (for most slip systems, not required)
dot_product = np.dot(slip_sys['n'], slip_sys['m'])
print(f"n·m = {dot_product:.6f}")  # Should be ~0 for most systems

# Calculate Schmid factor for uniaxial tension
# Stress tensor for tension along [100]
stress = np.array([[100, 0, 0],
                   [0,   0, 0],
                   [0,   0, 0]])  # MPa

P = slip_sys['P']
schmid_factor = np.tensordot(P, stress)
print(f"Schmid factor: {schmid_factor:.4f}")

# Maximum Schmid factor is 0.5 for uniaxial tension
# Occurs when plane normal and slip direction are both 45° to stress axis

# Example with different loading
# Pure shear in x-y plane
shear_stress = np.array([[0,  50, 0],
                         [50, 0,  0],
                         [0,  0,  0]])

schmid_shear = np.tensordot(P, shear_stress)
print(f"Schmid factor (shear): {schmid_shear:.4f}")
```

**Notes**:
- Plane normal uses reciprocal lattice (Lr)
- Direction uses direct lattice (L)
- Schmid tensor is outer product m ⊗ n
- Used for slip activity prediction in crystal plasticity

---

### `genEqSlipSys(phasename, SlipHKLs, SlipUVWs, getIndependentOnly=True, tol=1e-10, mag=1)`

**Purpose**: Generate all equivalent slip systems considering crystal symmetry

**Input**:
- `phasename`: str - Phase name
- `SlipHKLs`: list of lists - Multiple slip plane families [[h1,k1,l1], ...]
- `SlipUVWs`: list of lists - Multiple slip direction families [[u1,v1,w1], ...]
- `getIndependentOnly`: bool - Return only independent systems (default: True)
- `tol`: float - Tolerance for identifying duplicates (default: 1e-10)
- `mag`: float - Magnitude scaling

**Output**:
- list of slip system dictionaries (same format as `genSlipSys`)

**Usage Examples**:
```python
# BCC {110}<111> slip (48 combinations, 12 independent)
SlipHKLs_bcc = [[1,1,0], [-1,1,0], [1,0,1], 
                [-1,0,1], [0,1,1], [0,-1,1]]
SlipUVWs_bcc = [[1,1,1]]

slip_bcc = gp.genEqSlipSys('A', SlipHKLs_bcc, SlipUVWs_bcc,
                            getIndependentOnly=True)

print(f"BCC {{110}}<111> systems: {len(slip_bcc)}")  # Should be 12

# FCC {111}<110> slip (12 systems)
SlipHKLs_fcc = [[1,1,1]]
SlipUVWs_fcc = [[1,1,0], [1,-1,0], [0,1,1]]

slip_fcc = gp.genEqSlipSys('A', SlipHKLs_fcc, SlipUVWs_fcc,
                            getIndependentOnly=True)

print(f"FCC {{111}}<110> systems: {len(slip_fcc)}")  # Should be 12

# Calculate Schmid factors for all systems
tension = np.array([1, 0, 0])  # [100] direction
stress_tensor = np.outer(tension, tension) * 100  # 100 MPa

schmid_factors = []
for i, ss in enumerate(slip_bcc):
    schmid = np.tensordot(ss['P'], stress_tensor)
    schmid_factors.append(abs(schmid))
    print(f"System {i}: ({ss['hkl'][0]}{ss['hkl'][1]}{ss['hkl'][2]})"
          f"[{ss['uvw'][0]}{ss['uvw'][1]}{ss['uvw'][2]}], "
          f"μ = {abs(schmid):.4f}")

# Find most favorably oriented system
max_schmid = max(schmid_factors)
max_idx = schmid_factors.index(max_schmid)
active_system = slip_bcc[max_idx]

print(f"\nMost active system:")
print(f"  Plane: ({active_system['hkl'][0]}{active_system['hkl'][1]}"
      f"{active_system['hkl'][2]})")
print(f"  Direction: [{active_system['uvw'][0]}{active_system['uvw'][1]}"
      f"{active_system['uvw'][2]}]")
print(f"  Schmid factor: {max_schmid:.4f}")

# Analyze distribution
import matplotlib.pyplot as plt
plt.figure(figsize=(8, 6))
plt.bar(range(len(schmid_factors)), schmid_factors)
plt.xlabel('Slip System')
plt.ylabel('Schmid Factor')
plt.title(f'Schmid Factors for [100] Tension')
plt.axhline(y=0.5, color='r', linestyle='--', label='Maximum (0.5)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

**Notes**:
- Automatically applies crystal symmetry to generate all variants
- `getIndependentOnly=True` removes duplicates (same Schmid tensor)
- Essential for crystal plasticity finite element modeling
- Number of systems depends on crystal structure and slip family

---

## Utility and Helper Methods

### `getUnique(a)`

**Purpose**: Get unique rows from array

**Input**:
- `a`: ndarray - Input array

**Output**:
- ndarray - Array with unique rows only

**Usage Example**:
```python
# Array with duplicate rows
data = np.array([[1, 0, 0],
                 [0, 1, 0],
                 [1, 0, 0],  # Duplicate
                 [0, 0, 1]])

unique_data = gp.getUnique(data)
print("Original shape:", data.shape)      # (4, 3)
print("Unique shape:", unique_data.shape) # (3, 3)
```

---

### `getCVdict(name='NiTi')`

**Purpose**: Get correspondence variant dictionary

**Input**:
- `name`: str - Name for CV set

**Output**:
- None (populates `self.CV[name]` dictionary)

**Usage Example**:
```python
gp.getCVdict(name='NiTi')

# Access correspondence variants
cv_dict = gp.CV['NiTi']
print(f"Number of CVs: {len(cv_dict)}")

for cv_idx, cv_data in cv_dict.items():
    print(f"CV {cv_idx}:")
    print(f"  Cd matrix:\n{cv_data['Cd']}")
```

---

### `printHBVII011(name='NiTi')`

**Purpose**: Print {011} Type II habit plane variant information

**Input**:
- `name`: str - HPV set name

**Output**:
- None (prints to console)

**Usage Example**:
```python
gp.printHBVII011(name='NiTi')

# Typical output:
# HPV 0: (011)[100], λ=0.045
# HPV 1: (0-11)[100], λ=0.045
# ...
```

---

### `printCV(name='NiTi')`

**Purpose**: Print correspondence variant information

**Input**:
- `name`: str - CV set name

**Output**:
- None (prints to console)

**Usage Example**:
```python
gp.printCV(name='NiTi')

# Shows correspondence matrices for each CV
```

---

### `setPlotAttributes(defs={}, name='NiTi')`

**Purpose**: Set plotting attributes for visualizations

**Input**:
- `defs`: dict - Dictionary of plotting attributes
- `name`: str - Name to associate attributes with

**Output**:
- None (updates internal plot settings)

**Usage Example**:
```python
# Customize plot appearance
plot_settings = {
    'marker': 'o',
    'markersize': 10,
    'color': 'red',
    'linestyle': '-',
    'linewidth': 2
}

gp.setPlotAttributes(defs=plot_settings, name='NiTi')
```

---

## Complete Application Example

### NiTi Shape Memory Alloy Full Analysis

```python
#!/usr/bin/env python3
"""
Complete NiTi shape memory alloy crystallographic analysis
"""
from getphases import getPhases
from orilib import np_eulers_matrices
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 1. Initialize and Define Phases
# ============================================================
gp = getPhases()

# Define B2 austenite (cubic, a = 3.015 Å)
gp.fromLatticeParams('A', a=3.015)

# Define B19' martensite (monoclinic)
gp.fromLatticeParams('M', 
    a=2.889,  # Å
    b=4.120,  # Å
    c=4.622,  # Å
    beta=np.radians(96.8)
)

print("Phases defined:")
print(f"  Austenite: cubic, a = {gp.phases['A']['a']:.3f} Å")
print(f"  Martensite: monoclinic, a={gp.phases['M']['a']:.3f}, "
      f"b={gp.phases['M']['b']:.3f}, c={gp.phases['M']['c']:.3f} Å, "
      f"β={np.degrees(gp.phases['M']['beta']):.1f}°")

# ============================================================
# 2. Orientation Relationship and Deformation Gradients
# ============================================================
gp.getOR(name='NiTi', OR='NiTi')
gp.getDefGrad(name='NiTi')

print(f"\nOrientation relationship established")
print(f"Correspondence variants: {len(gp.defGrad['NiTi'])}")

# Analyze deformation gradients
for i, cv in enumerate(gp.defGrad['NiTi']):
    F = cv['F']
    det_F = np.linalg.det(F)
    print(f"  CV{i}: det(F) = {det_F:.6f}")

# ============================================================
# 3. Elastic Properties
# ============================================================
gp.getEl(name='NiTi')

C_aust = gp.el['NiTi']['A']['C']
C_mart = gp.el['NiTi']['M']['C']

print("\nElastic stiffness (GPa):")
print(f"  Austenite C11: {C_aust[0,0]/1e3:.1f}")
print(f"  Martensite C11: {C_mart[0,0]/1e3:.1f}")

# ============================================================
# 4. Twin Systems
# ============================================================
gp.getTwinSys(name='NiTi')

print(f"\nTwin systems:")
print(f"  Type I: {len(gp.twinSys['NiTi']['I'])}")
print(f"  Type II: {len(gp.twinSys['NiTi']['II'])}")

# ============================================================
# 5. Habit Plane Variants
# ============================================================
gp.getHBVs(name='NiTi')
gp.getHPVII011(name='NiTi')

hpvs = gp.vardict['NiTi']['HPVII011']
print(f"\n{{011}} Type II HPVs: {len(hpvs)}")

# ============================================================
# 6. Transformation Strain Analysis
# ============================================================
# Generate 1000 random austenite orientations
N = 1000
euler_angles = np.random.rand(N, 3) * [360, 180, 360]
oris = np_eulers_matrices(euler_angles, deg=True)

# Calculate transformation strains
strains = gp.getTrStrain(phase='A', name='NiTi', oris=oris)

print(f"\nTransformation strain analysis ({N} orientations):")
labels = ['ε₁₁', 'ε₂₂', 'ε₃₃', 'ε₂₃', 'ε₁₃', 'ε₁₂']
for i, label in enumerate(labels):
    mean = np.mean(strains[:, i])
    std = np.std(strains[:, i])
    print(f"  {label}: mean={mean:.4f}, std={std:.4f}")

# ============================================================
# 7. Visualization
# ============================================================
# Plot habit plane variants
gp.plotHBVII011(name='NiTi', ProjType='equalarea')
plt.title('NiTi {011} Type II Habit Plane Variants')
plt.savefig('niti_hpv.png', dpi=300, bbox_inches='tight')
plt.show()

# Plot strain distribution
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
for i, (ax, label) in enumerate(zip(axes.flat, labels)):
    ax.hist(strains[:, i], bins=50, alpha=0.7, edgecolor='black')
    ax.set_xlabel(label, fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title(f'Distribution of {label}', fontsize=14)
    ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('niti_strain_dist.png', dpi=300, bbox_inches='tight')
plt.show()

print("\n✓ Analysis complete!")
print("  - Figures saved: niti_hpv.png, niti_strain_dist.png")
```

---

## Dependencies

Required Python packages:
- `numpy`: Array operations and linear algebra
- `crystals`: CIF file parsing
- `orilib`: Orientation operations (quaternions, Euler angles)
- `projlib`: Stereographic projections
- `crystlib`: Crystallographic utilities
- `matplotlib`: Visualization
- `orix` (optional): Advanced crystallographic operations

---

## Best Practices

1. **Always use radians** for lattice angles
2. **Define phases before** calculating ORs
3. **Calculate OR before** deformation gradients
4. **Verify lattice matrices** (check det(L @ Lr.T) = I)
5. **Store configurations** for reproducibility
6. **Document parameter sources** (literature references)
7. **Validate against** experimental data when possible

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| KeyError for phase | Define phase with `fromCif()` or `fromLatticeParams()` |
| Wrong β angle | Use `np.radians()` to convert degrees |
| Empty defGrad | Call `getOR()` before `getDefGrad()` |
| Asymmetric lattice | Check parameter order: a, b, c, α, β, γ |
| Mismatch in calculations | Verify column format of lattice matrix |

---

*This comprehensive documentation covers all methods in getphases.py. For quick reference, see getphases_quick_reference.md*

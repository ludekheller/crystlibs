# Crystlib.py Documentation - Available Files Summary

## ✅ Documentation Already Created!

I have already created comprehensive documentation for **crystlib.py** in this session. Here's what's available:

---

## 📦 Available Files

### 1. **crystlib_quick_reference.md** (17 KB, ~900 lines)

**Purpose**: Quick reference guide for crystallographic calculations

**Contents**:
- ✅ Quick start examples
- ✅ Most commonly used functions (7 functions)
- ✅ Complete workflow examples (5 detailed examples)
- ✅ Crystal system parameters reference
- ✅ Function parameters quick reference
- ✅ Common crystal structures (Si, Al, Ti, Fe, Cu, Zn, Mg)
- ✅ Mathematical formulas (reciprocal lattice, d-spacing, Bragg's law)
- ✅ Common pitfalls and solutions
- ✅ Integration examples (with projlib, plotlib, orilib)
- ✅ Troubleshooting guide
- ✅ Performance tips

**Major Sections**:

1. **Overview** - Module purpose and capabilities
2. **Quick Start** - Immediate usage examples
3. **Most Commonly Used Functions**:
   - Miller Indices Generation (3 functions)
   - Lattice Vector Calculations (1 function, 5 crystal systems)
   - Reciprocal Lattice (1 function)
   - Utility Functions (2 functions)
4. **Complete Workflow Examples**:
   - Cubic Crystal Analysis
   - Generate Diffraction Planes with Symmetry
   - Hexagonal Crystal (Titanium)
   - Diffraction Analysis
   - Miller Index Families
5. **Function Parameters Quick Reference**
6. **Crystal System Parameters** (all 7 systems)
7. **Common Crystal Structures** (examples for Si, Al, Fe, Cu, Ti, Zn, Mg)
8. **Mathematical Formulas** (reciprocal lattice, d-spacing, Bragg's law, volume)
9. **Common Pitfalls** (angle units, symmetry matrices, Miller index format)
10. **Integration with Other Modules** (projlib, plotlib)
11. **Troubleshooting Table**
12. **Performance Tips**

---

### 2. **crystlib_documentation_summary.md** (19 KB, ~1,000 lines)

**Purpose**: Comprehensive function-by-function documentation

**Contents**:
- ✅ Overview and module purpose
- ✅ Documentation structure
- ✅ Complete function documentation (7 functions)
- ✅ Detailed function specifications:
  - Purpose and algorithm
  - Input parameters (types, formats, constraints)
  - Output values (types, meanings)
  - Features and properties
  - Usage examples (3-5 per function)
  - Applications
- ✅ Mathematical background:
  - Crystal systems (7 systems, 14 Bravais lattices)
  - Miller indices notation
  - Reciprocal lattice theory
  - Metric tensor operations
  - d-spacing formulas
  - Bragg's law
- ✅ Integration examples (with orilib, projlib, plotlib)
- ✅ Best practices (symmetry ops, angle units, validation)
- ✅ Common applications (XRD, EBSD, texture, structure factor)
- ✅ Performance considerations
- ✅ Troubleshooting guide
- ✅ Dependencies and references

**Functions Documented**:

#### 1. Miller Indices Generation Functions:

**`generate_hkls(hklmax, syms, hkls=[])`**
- Generate unique Miller indices with symmetry
- Handles any crystal system via symmetry operations
- Returns unique hkls, equivalents dictionary, families
- Includes rounding for robust floating point comparison
- **Usage**: X-ray/neutron diffraction, EBSD indexing

**`generate_hkls01(hklmax, syms, hkls=[])`**
- Alternative version without rounding
- Exact floating point comparison
- Same output format as generate_hkls()
- **Usage**: When exact precision required

**`generate_hkls02(hklmax, syms, G, hkls=[])`**
- Generates hkls with metric tensor transformation
- For non-orthogonal crystal systems
- Applies proper reciprocal space symmetry
- **Usage**: Hexagonal, monoclinic, triclinic systems

#### 2. Family Grouping:

**`get_unique_families(hkls)`**
- Groups Miller indices into families
- Based on permutation symmetry
- Returns representative hkl and multiplicity
- **Usage**: Structure factor calculations, powder diffraction

#### 3. Lattice Calculations:

**`lattice_vec(lattice_param)`**
- Calculates real-space lattice vectors
- Supports all 7 crystal systems:
  - Cubic (a = b = c, α = β = γ = 90°)
  - Tetragonal (a = b ≠ c, α = β = γ = 90°)
  - Hexagonal/Trigonal (a = b ≠ c, α = β = 90°, γ = 120°)
  - Orthorhombic (a ≠ b ≠ c, α = β = γ = 90°)
  - Monoclinic (a ≠ b ≠ c, α = γ = 90° ≠ β)
  - Triclinic (a ≠ b ≠ c, α ≠ β ≠ γ)
- Returns three lattice vectors as numpy arrays
- **Usage**: Crystal structure visualization, unit cell calculations

**`reciprocal_basis(a1, a2, a3)`**
- Calculates reciprocal lattice vectors
- Satisfies orthogonality: aᵢ · bⱼ = δᵢⱼ
- Uses crystallographic convention (without 2π)
- Returns b1, b2, b3
- **Usage**: Diffraction analysis, d-spacing calculations

#### 4. Utility Functions:

**`array2tuple(arr, decimals=2)`**
- Converts array to rounded tuple
- Configurable precision
- Used internally for Miller index comparison

---

### 3. **crystlib_commented.py** (19 KB, 522 lines)

**Purpose**: Fully documented source code

**Contents**:
- ✅ Module-level docstring with overview
- ✅ Complete docstrings for all 7 functions
- ✅ Input/output descriptions
- ✅ Usage examples for each function
- ✅ Mathematical context
- ✅ Crystal system implementations

---

## 📊 Documentation Statistics

| Aspect | Details |
|--------|---------|
| **Functions Documented** | 7 main functions |
| **Crystal Systems Supported** | 7 systems (all Bravais lattices) |
| **Total Documentation** | ~1,900 lines across 3 files |
| **File Sizes** | Quick Ref: 17 KB, Summary: 19 KB, Code: 19 KB |
| **Usage Examples** | 3-5 per function (~30 total examples) |
| **Workflow Examples** | 5 complete workflows |
| **Mathematical Formulas** | 10+ formulas with explanations |
| **Crystal Structure Examples** | 7+ common materials |

---

## 🎯 What's Covered

### Core Functionality:
✅ **Miller indices generation** - With crystal symmetry  
✅ **All 7 crystal systems** - Cubic to triclinic  
✅ **Lattice vector calculations** - For any crystal  
✅ **Reciprocal lattice** - Crystallographic convention  
✅ **Family grouping** - Permutation symmetry  
✅ **Metric tensor operations** - For non-orthogonal systems  

### Crystal Systems:
✅ **Cubic** - Si, Al, Fe, Cu (examples provided)  
✅ **Tetragonal** - Parameter specifications  
✅ **Hexagonal/Trigonal** - Ti, Zn, Mg (examples provided)  
✅ **Orthorhombic** - General case  
✅ **Monoclinic** - With β angle  
✅ **Triclinic** - Full generality  

### Mathematical Background:
✅ Reciprocal lattice vectors  
✅ Orthogonality conditions (aᵢ · bⱼ = δᵢⱼ)  
✅ d-spacing formulas  
✅ Bragg's law (nλ = 2d sin θ)  
✅ Unit cell volume calculations  
✅ Metric tensor operations  

### Practical Applications:
✅ X-ray/neutron diffraction analysis  
✅ EBSD pattern indexing  
✅ Texture analysis  
✅ Structure factor calculations  
✅ Reciprocal space mapping  

---

## 📖 Quick Navigation

### For Beginners:
👉 **Start with**: `crystlib_quick_reference.md`
- Read "Quick Start" section
- Try "Common Crystal Structures" examples
- Follow "Complete Workflow Examples"

### For Crystal System Setup:
👉 **Use**: "Crystal System Parameters" table
- Find your crystal system (cubic, hexagonal, etc.)
- Copy parameter dictionary format
- Plug in your lattice parameters

### For Diffraction Analysis:
👉 **Reference**: "Diffraction Analysis" workflow
- Generate Miller indices
- Calculate d-spacings
- Compute Bragg angles
- Analyze reflections

### For Advanced Users:
👉 **Reference**: `crystlib_documentation_summary.md`
- Mathematical background section
- Metric tensor operations
- Integration patterns
- Full function specifications

---

## 🚀 Key Features Documented

### Crystal System Support:

**Cubic Systems** (a only):
```python
Si = {'type': 'cubic', 'a': 5.43}  # Silicon
Al = {'type': 'cubic', 'a': 4.05}  # Aluminum
Fe = {'type': 'cubic', 'a': 2.87}  # Iron (BCC)
Cu = {'type': 'cubic', 'a': 3.61}  # Copper (FCC)
```

**Hexagonal Systems** (a, c):
```python
Ti = {'type': 'trigonal', 'a': 2.95, 'c': 4.68}  # Titanium
Zn = {'type': 'trigonal', 'a': 2.66, 'c': 4.95}  # Zinc
Mg = {'type': 'trigonal', 'a': 3.21, 'c': 5.21}  # Magnesium
```

**Monoclinic** (a, b, c, β):
```python
params = {
    'type': 'monoclinic',
    'a': 5.0, 'b': 6.0, 'c': 7.0,
    'beta': np.radians(120)  # Must be in radians!
}
```

**Triclinic** (a, b, c, α, β, γ):
```python
params = {
    'type': 'triclinic',
    'a': 5.0, 'b': 6.0, 'c': 7.0,
    'alpha': np.radians(90),
    'beta': np.radians(95),
    'gamma': np.radians(100)
}
```

### Miller Indices Examples:

```python
# Generate reflections for cubic crystal
identity = np.eye(3)
# Add cubic symmetry operations (48 total)
syms = [identity, ...]  # Add more operations

hkls, hkls2, families = generate_hkls(3, syms)

# Results:
# hkls: [(1.0, 0.0, 0.0), (1.0, 1.0, 0.0), ...]
# hkls2: {(1.0, 0.0, 0.0): [(1,0,0), (0,1,0), (0,0,1), ...]}
# families: {(1, 0, 0): 6, (1, 1, 0): 12, (1, 1, 1): 8}
```

---

## 📥 Download Files

**Quick Reference**:
[crystlib_quick_reference.md](computer:///mnt/user-data/outputs/crystlib_quick_reference.md)

**Documentation Summary**:
[crystlib_documentation_summary.md](computer:///mnt/user-data/outputs/crystlib_documentation_summary.md)

**Commented Source Code**:
[crystlib_commented.py](computer:///mnt/user-data/outputs/crystlib_commented.py)

**All Files**:
[View All Files](computer:///mnt/user-data/outputs/)

---

## 🔍 Sample Content Preview

### From Quick Reference:

```python
# Example 1: Cubic Crystal Analysis
import numpy as np
from crystlib import *

# Define silicon (cubic, a = 5.43 Å)
lattice_params = {'type': 'cubic', 'a': 5.43}

# Get real space lattice
a1, a2, a3 = lattice_vec(lattice_params)
print("Lattice vectors:")
print("a1:", a1)  # [5.43, 0, 0]
print("a2:", a2)  # [0, 5.43, 0]
print("a3:", a3)  # [0, 0, 5.43]

# Calculate reciprocal lattice
b1, b2, b3 = reciprocal_basis(a1, a2, a3)
print("\nReciprocal lattice:")
print("b1:", b1)  # [1/5.43, 0, 0]
print("b1 magnitude:", np.linalg.norm(b1))  # 0.184

# Volume of unit cell
V = np.dot(a1, np.cross(a2, a3))
print(f"\nUnit cell volume: {V:.3f} Ų")  # 160.1 Ų
```

### From Documentation Summary:

#### Function: `lattice_vec(lattice_param)`

**Purpose**: Calculate real-space lattice vectors for any crystal system

**Supported Crystal Systems**:

1. **Cubic** (a = b = c, α = β = γ = 90°)
   - Parameters: `{'type': 'cubic', 'a': float}`
   - Examples: Si, Al, Fe, Cu

2. **Hexagonal** (a = b ≠ c, α = β = 90°, γ = 120°)
   - Parameters: `{'type': 'trigonal', 'a': float, 'c': float}`
   - Examples: Ti, Zn, Mg
   - c/a ratio: ~1.633 for ideal HCP

3. **Monoclinic** (a ≠ b ≠ c, α = γ = 90° ≠ β)
   - Parameters: `{'type': 'monoclinic', 'a': float, 'b': float, 'c': float, 'beta': float}`
   - **Critical**: beta must be in radians!

**Output**:
- `a1`, `a2`, `a3`: numpy arrays (3,) - Three lattice vectors

**Mathematical Formulas**:

For **Hexagonal**:
```
a₁ = [a/2, -a√3/2, 0]
a₂ = [a/2,  a√3/2, 0]
a₃ = [0, 0, c]
Angle between a₁ and a₂ = 120°
```

For **Cubic**:
```
a₁ = [a, 0, 0]
a₂ = [0, a, 0]
a₃ = [0, 0, a]
```

---

## 💡 Pro Tips

1. **Angle Units**: Always use radians! Convert with `np.radians(degrees)`
2. **Symmetry Operations**: Ensure matrices are orthogonal (det = 1)
3. **Miller Indices**: Functions return tuples like `(1.0, 0.0, 0.0)`, not lists
4. **Reciprocal Convention**: This uses crystallographic convention (no 2π)
5. **Validate Results**: Check orthogonality (aᵢ · bⱼ = δᵢⱼ) and V × V* = 1
6. **Cache Calculations**: Compute reciprocal lattice once, reuse for all d-spacings
7. **Start Simple**: Use cubic examples first, then progress to complex systems

---

## 🎯 Common Use Cases Covered

### 1. X-ray Diffraction Analysis
```python
from crystlib import *
import numpy as np

# Setup crystal
cubic = {'type': 'cubic', 'a': 5.43}
a1, a2, a3 = lattice_vec(cubic)
b1, b2, b3 = reciprocal_basis(a1, a2, a3)

# Generate reflections
syms = [np.eye(3)]  # Add full cubic symmetry
hkls, hkls2, fam = generate_hkls(3, syms)

# Calculate d-spacings
wavelength = 1.5406  # Cu Kα
for hkl in hkls[:10]:
    G = hkl[0]*b1 + hkl[1]*b2 + hkl[2]*b3
    d = 1.0 / np.linalg.norm(G)
    
    # Bragg angle
    sin_theta = wavelength / (2*d)
    if sin_theta <= 1.0:
        two_theta = 2 * np.degrees(np.arcsin(sin_theta))
        print(f"({int(hkl[0])}{int(hkl[1])}{int(hkl[2])})  d={d:.3f}Å  2θ={two_theta:.2f}°")
```

### 2. EBSD Pattern Indexing
```python
from crystlib import lattice_vec, reciprocal_basis
from projlib import gen_dirs_norms

# Setup crystal structure
cubic = {'type': 'cubic', 'a': 3.5}
a1, a2, a3 = lattice_vec(cubic)
b1, b2, b3 = reciprocal_basis(a1, a2, a3)

# Create basis matrices
L = np.column_stack([a1, a2, a3])
Lr = np.column_stack([b1, b2, b3])

# Generate directions for pattern simulation
dirs, norms = gen_dirs_norms(
    L, Lr,
    uvws=[[1,0,0], [1,1,0], [1,1,1]],
    hkls=[[1,0,0], [1,1,0], [1,1,1]]
)
```

### 3. Hexagonal Crystal Setup
```python
from crystlib import lattice_vec, reciprocal_basis
import numpy as np

# Titanium (hexagonal)
Ti_params = {'type': 'trigonal', 'a': 2.95, 'c': 4.68}
a1, a2, a3 = lattice_vec(Ti_params)

# Check c/a ratio
c_over_a = np.linalg.norm(a3) / np.linalg.norm(a1)
print(f"c/a = {c_over_a:.3f}")  # ~1.586 (ideal: 1.633)

# Verify 120° angle
angle = np.arccos(np.dot(a1, a2) / (np.linalg.norm(a1)**2))
print(f"Angle: {np.degrees(angle):.1f}°")  # 120°

# Reciprocal lattice
b1, b2, b3 = reciprocal_basis(a1, a2, a3)
print(f"c* length: {np.linalg.norm(b3):.4f} ų⁻¹")
```

---

## 🔧 Integration Examples

### With Orilib (Full Workflow):
```python
from crystlib import lattice_vec, reciprocal_basis
from orilib import np_euler_matrix

# Define crystal
cubic = {'type': 'cubic', 'a': 3.5}
a1, a2, a3 = lattice_vec(cubic)
L = np.column_stack([a1, a2, a3])

# Define orientation
euler = [0.5, 0.3, 0.1]  # radians
g = np_euler_matrix(*euler)

# Transform direction
dir_crystal = np.array([1, 0, 0])
dir_sample = g @ dir_crystal
print("Crystal [100] in sample frame:", dir_sample)
```

### With Projlib (Pole Figure):
```python
from crystlib import lattice_vec, reciprocal_basis, generate_hkls
from projlib import equalarea_directions

cubic = {'type': 'cubic', 'a': 3.5}
a1, a2, a3 = lattice_vec(cubic)
b1, b2, b3 = reciprocal_basis(a1, a2, a3)

# Generate plane normals
syms = [np.eye(3)]
hkls, hkls2, fam = generate_hkls(2, syms)

# Convert to vectors (plane normals in reciprocal space)
normals = np.array([h*b1 + k*b2 + l*b3 
                    for h, k, l in hkls[:20]]).T

# Project for pole figure
proj = equalarea_directions(normals)
```

---

## ✅ Summary

**You have complete, comprehensive documentation for crystlib.py including:**

✅ Quick reference guide (17 KB)  
✅ Documentation summary (19 KB)  
✅ Commented source code (19 KB)  
✅ 7 functions documented  
✅ 30+ usage examples  
✅ 5 complete workflows  
✅ All 7 crystal systems  
✅ Mathematical formulas  
✅ Integration patterns  
✅ Troubleshooting guides  
✅ Common crystal structures  

**Everything is ready to use for:**
- X-ray/neutron diffraction
- EBSD analysis
- Crystal structure calculations
- Reciprocal space mapping
- Texture analysis
- Research and publication

---

*All files are available in `/mnt/user-data/outputs/`*

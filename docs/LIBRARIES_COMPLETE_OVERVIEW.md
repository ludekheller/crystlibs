# CRYSTALLOGRAPHY LIBRARIES - COMPLETE OVERVIEW
## Extended Toolkit for Materials Science Research

**Version**: Extended 2024  
**Total Functions**: 192  
**Total Size**: 356 KB  
**Language**: Python 2/3 Compatible  

---

## 📚 LIBRARY MODULES

### 1. crystlib.py
**Crystal Structure and Lattice Operations**

- **Functions**: 73
- **File Size**: 116 KB
- **Primary Focus**: Core crystallographic calculations and lattice operations

**Key Capabilities**:
- Lattice vector generation for all 7 crystal systems
- Miller index operations and conversions
- Coordinate transformations (Cartesian ↔ fractional)
- Crystal symmetry operations
- Lattice correspondence matrices
- Twinning and deformation analysis
- NiTi shape memory alloy specific functions
- Atomic position generation
- Hexagonal slip systems

**Use Cases**:
- Phase transformation analysis
- Twinning studies
- Deformation mechanics
- Lattice correspondence calculations
- Crystallographic indexing

---

### 2. orilib.py
**Orientation Analysis and Transformations**

- **Functions**: 43
- **File Size**: 60 KB
- **Primary Focus**: Rotation representations and orientation mathematics

**Key Capabilities**:
- Quaternion operations (multiplication, conjugation, misorientation)
- Euler angle conversions (Bunge convention, ZXZ)
- Rodrigues-Frank vector representations
- Rotation matrices and axis-angle conversions
- Orientation sampling (HEALPix, Hopf coordinates)
- Misorientation calculations with symmetry
- Active and passive rotations
- Orientation space gridding

**Use Cases**:
- Texture analysis
- EBSD data processing
- Misorientation calculations
- Orientation distribution functions
- Grain boundary characterization

---

### 3. plotlib.py
**Crystallographic Plotting and Visualization**

- **Functions**: 30
- **File Size**: 68 KB
- **Primary Focus**: Publication-quality crystallographic visualization

**Key Capabilities**:
- Stereographic projections (Schmidt/Wulff nets)
- Pole figures and inverse pole figures
- Mohr circles for 3D strain analysis
- 3D lattice visualization
- 2D lattice projections
- Atomic plane plotting
- Orientation density plots
- Custom colormaps and color scales
- Interactive data annotations

**Use Cases**:
- Publication figures
- Texture visualization
- Strain analysis visualization
- Crystal orientation mapping
- Lattice structure presentation

---

### 4. projlib.py
**Stereographic Projections and Coordinate Transformations**

- **Functions**: 46
- **File Size**: 112 KB
- **Primary Focus**: Stereographic projection operations

**Key Capabilities**:
- Equal-area (Schmidt) projections
- Equal-angle (Wulff) projections
- Spherical to Cartesian conversions
- Stereographic grid generation
- Pole figure calculations
- Inverse pole figure coloring
- Standard stereographic triangles
- Direction and plane projections

**Use Cases**:
- Pole figure generation
- IPF color mapping
- Crystallographic orientation display
- Texture component identification
- Orientation relationship analysis

---

## 🎯 INTEGRATED WORKFLOW

These four libraries work together seamlessly:

```python
import crystlib
import orilib
import plotlib
import projlib

# 1. Define crystal structure (crystlib)
lattice = crystlib.cubic_lattice_vec(a=3.6)
reciprocal = crystlib.reciprocal_basis(lattice)

# 2. Generate orientations (orilib)
euler_angles = [45, 30, 60]  # degrees
rotation_matrix = orilib.np_euler_matrix(*np.radians(euler_angles))

# 3. Project directions (projlib)
directions = [[1,0,0], [1,1,0], [1,1,1]]
x, y = projlib.equalarea_directions(directions)

# 4. Visualize (plotlib)
p = plotlib.plotter()
p.plotProj(ProjType='equalarea', sphere='half')
p.setAttributes(dirs=directions)
p.plotDirsNorms()
p.figsave(fname='pole_figure.png')
```

---

## 📊 FUNCTION DISTRIBUTION

### By Category

| Category | crystlib | orilib | plotlib | projlib | Total |
|----------|----------|--------|---------|---------|-------|
| **Lattice Operations** | 25 | 0 | 0 | 0 | 25 |
| **Orientations** | 0 | 43 | 0 | 0 | 43 |
| **Plotting** | 0 | 0 | 30 | 0 | 30 |
| **Projections** | 0 | 0 | 0 | 46 | 46 |
| **Utilities** | 15 | 0 | 0 | 0 | 15 |
| **Twinning** | 12 | 0 | 0 | 0 | 12 |
| **Deformation** | 5 | 0 | 0 | 0 | 5 |
| **Miller Indices** | 10 | 0 | 0 | 0 | 10 |
| **Symmetry** | 6 | 0 | 0 | 0 | 6 |
| **TOTAL** | **73** | **43** | **30** | **46** | **192** |

---

## 🚀 QUICK START EXAMPLES

### Example 1: Basic Lattice Operations
```python
import crystlib
import numpy as np

# Create cubic lattice
a = 3.6  # Angstroms
lattice = crystlib.cubic_lattice_vec(a)

# Get reciprocal lattice
reciprocal = crystlib.reciprocal_basis(lattice)

# Convert Miller indices to Cartesian
hkl = [1, 1, 1]
cart = crystlib.miller2fractional(hkl, lattice)

print(f"Lattice vectors:\n{lattice}")
print(f"[111] direction: {cart}")
```

### Example 2: Orientation Analysis
```python
import orilib
import numpy as np

# Create rotation from Euler angles
phi1, Phi, phi2 = np.radians([45, 30, 60])
g = orilib.np_euler_matrix(phi1, Phi, phi2)

# Convert to axis-angle
axis, angle = orilib.np_ol_g_rtheta_rad(g)

# Convert to quaternion
q = orilib.mat_to_quat(g)

print(f"Rotation axis: {axis}")
print(f"Rotation angle: {np.degrees(angle):.2f}°")
print(f"Quaternion: {q}")
```

### Example 3: Pole Figure Creation
```python
import plotlib
import numpy as np

# Create plotter
p = plotlib.plotter()

# Setup equal-area projection
p.plotProj(ProjType='equalarea', sphere='half', figsize=(6,6))

# Plot crystal directions
directions = [[1,0,0], [0,1,0], [0,0,1], [1,1,1]]
p.setAttributes(dirs=directions)
p.plotDirsNorms()

# Save figure
p.figsave(fname='pole_figure.png', imformats=['png', 'pdf'])
```

### Example 4: Mohr Circle Analysis
```python
import plotlib
import numpy as np

# Define strain tensor
strain = np.array([[0.02, 0.01, 0.00],
                   [0.01, 0.03, 0.00],
                   [0.00, 0.00, 0.01]])

# Calculate Mohr circles
mohr_data = plotlib.mohr_circles(strain)

# Plot Mohr circles
plotlib.plot_mohr_circles(mohr_data)
```

### Example 5: NiTi Twinning Analysis
```python
import crystlib
import numpy as np

# B19' lattice parameters
a, b, c = 2.89, 4.12, 4.62  # Angstroms
beta = 96.8  # degrees

lattice_B19p = crystlib.monoclinic_lattice_vec(a, b, c, beta)

# Calculate twinning data
twin_data = crystlib.niti_twinning(lattice_B19p)

# Get habit plane
habit_plane = crystlib.habitplane_equation_solution(twin_data)

print(f"Twinning plane: {habit_plane}")
```

---

## 📦 INSTALLATION & SETUP

### Requirements
```
numpy>=1.16.0
matplotlib>=3.0.0
scipy>=1.3.0
numba>=0.45.0 (for orilib optimization)
```

### Installation
```bash
# Place the four .py files in your Python path
# Or add to your project directory

# Import in your scripts
import crystlib
import orilib
import plotlib
import projlib
```

---

## 📖 DOCUMENTATION STRUCTURE

For each library, we provide:

1. **DOCUMENTATION_COMPLETE.md**
   - Comprehensive reference with all functions
   - Detailed descriptions and mathematical background
   - Complete input/output specifications
   - Multiple usage examples per function
   - Implementation notes and caveats

2. **DOCUMENTATION_SUMMARY.md**
   - Organized overview by category
   - Function groups and relationships
   - Workflow examples
   - Best practices guide

3. **QUICK_REFERENCE_COMPREHENSIVE.md**
   - Alphabetical function listing
   - One-line descriptions
   - Quick syntax reference
   - Parameter summary
   - Return value summary

4. **QUICK_REFERENCE.md**
   - Ultra-condensed reference
   - Essential functions only
   - Minimal examples
   - Quick lookup table

---

## 🔧 COMMON WORKFLOWS

### Workflow 1: EBSD Data Analysis
```
1. Import orientation data
2. Convert to rotation matrices (orilib)
3. Calculate misorientations (orilib)
4. Create pole figures (plotlib + projlib)
5. Generate IPF maps (projlib)
```

### Workflow 2: Phase Transformation
```
1. Define parent and product lattices (crystlib)
2. Calculate correspondence matrix (crystlib)
3. Determine variant orientations (orilib)
4. Plot orientation relationships (plotlib)
5. Analyze habit planes (crystlib)
```

### Workflow 3: Texture Analysis
```
1. Load orientation data
2. Apply crystal symmetry (orilib)
3. Calculate ODF (orilib + projlib)
4. Create pole figures (plotlib)
5. Identify texture components (projlib)
```

### Workflow 4: Deformation Analysis
```
1. Define deformation gradient (crystlib)
2. Calculate strain tensor (crystlib)
3. Find principal strains (crystlib)
4. Create Mohr circles (plotlib)
5. Plot strain distribution (plotlib)
```

---

## 📊 PERFORMANCE NOTES

### Optimized Functions
- **orilib**: Numba-accelerated quaternion operations
- **crystlib**: Vectorized NumPy operations for lattice calculations
- **projlib**: Fast stereographic projection algorithms

### Memory Considerations
- Large orientation datasets: Use streaming for >100k orientations
- High-resolution plots: Adjust `nump` parameter
- Batch processing: Process in chunks for memory efficiency

### Speed Tips
- Use NumPy versions (`np_*`) instead of list versions (`ol_*`)
- Pre-compute symmetry operations
- Cache frequently used projections
- Use lower resolution for preview, high for final output

---

## 🐛 TROUBLESHOOTING

### Common Issues

**Import Errors**
```python
# Solution: Ensure all dependencies installed
pip install numpy scipy matplotlib numba
```

**Projection Errors**
```python
# Issue: Invalid directions for projection
# Solution: Normalize directions first
dirs = dirs / np.linalg.norm(dirs, axis=1, keepdims=True)
```

**Figure Not Displaying**
```python
# Solution: Add plt.show()
import matplotlib.pyplot as plt
plt.show()
```

**Quaternion Normalization**
```python
# Solution: Normalize quaternions
q = q / np.linalg.norm(q)
```

---

## 📚 RELATED DOCUMENTATION

- [CRYSTLIB_DOCUMENTATION_COMPLETE.md](computer:///mnt/user-data/outputs/CRYSTLIB_DOCUMENTATION_COMPLETE.md)
- [ORILIB_DOCUMENTATION_COMPLETE.md](computer:///mnt/user-data/outputs/ORILIB_DOCUMENTATION_COMPLETE.md)
- [PLOTLIB_DOCUMENTATION_COMPLETE.md](computer:///mnt/user-data/outputs/PLOTLIB_DOCUMENTATION_COMPLETE.md)
- [PROJLIB_DOCUMENTATION_COMPLETE.md](computer:///mnt/user-data/outputs/PROJLIB_DOCUMENTATION_COMPLETE.md)

---

**Last Updated**: December 2024  
**Status**: Production Ready  
**Total Documentation**: 356 KB of code + comprehensive documentation suite

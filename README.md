# Crystallographic Analysis Toolkit

A comprehensive Python library for crystallographic calculations, orientation analysis, stereographic projections, and visualization in materials science research.

## 📚 Overview

This toolkit consists of four main modules designed for crystallographic analysis, particularly focused on:
- Crystal structure operations and lattice transformations
- Orientation representations and calculations
- Stereographic projections and pole figures
- Publication-quality crystallographic visualizations

**Key Applications:**
- Shape memory alloys (NiTi, Cu-based)
- Martensitic phase transformations
- EBSD data analysis
- Crystallographic texture analysis
- Deformation mechanics and twinning

---

## 📦 Modules

### 1. crystlib.py - Crystal Structure & Lattice Operations

Comprehensive functions for crystallographic calculations including lattice vectors, Miller indices, deformation analysis, and twinning operations.

**Key Features:**
- Lattice vector generation (cubic, monoclinic, tetragonal, hexagonal)
- Miller index operations and conversions
- Lattice correspondence matrices (B19'↔B2, cubic↔tetragonal)
- Deformation gradient and stress-free strain calculations
- Twinning analysis for shape memory alloys
- Mohr circle analysis for 3D strains

**Documentation:**
- 📘 [Complete Documentation](./docs/crystlib_DOCUMENTATION_COMPLETE.md)
- 📄 [Documentation Summary](./docs/crystlib_DOCUMENTATION_SUMMARY.md)
- 📋 [Quick Reference](./docs/crystlib_QUICK_REFERENCE.md)
- 📖 [Comprehensive Quick Reference](./docs/crystlib_QUICK_REFERENCE_COMPREHENSIVE.md)

---

### 2. orilib.py - Orientation Analysis & Transformations

Tools for crystallographic orientation representations including Euler angles, quaternions, rotation matrices, and Rodrigues-Frank vectors.

**Key Features:**
- Euler angle ↔ rotation matrix conversions
- Quaternion operations and misorientation calculations
- Rodrigues-Frank vector representations
- Axis-angle conversions
- Active and passive rotations
- HEALPix orientation sampling
- Crystal symmetry operations
- **Misorientation and disorientation analysis**
- **Symmetry-reduced orientations**

**Documentation:****Documentation:**
- 📘 [Complete Documentation](./docs/orilib_DOCUMENTATION_COMPLETE.md)
- 📄 [Documentation Summary](./docs/orilib_DOCUMENTATION_SUMMARY.md)
- 📋 [Quick Reference](./docs/orilib_QUICK_REFERENCE.md)
- 📖 [Comprehensive Quick Reference](./docs/orilib_QUICK_REFERENCE_COMPREHENSIVE.md)

---

### 3. plotlib.py - Crystallographic Visualization

Plotting functions for creating publication-quality crystallographic figures including pole figures, lattice structures, and Mohr circles.

**Key Features:**
- 3D lattice structure visualization
- 2D lattice projections
- Mohr circle plotting for strain analysis
- Atomic plane visualization
- Stereographic projection plotting
- Pole figure generation with density contours
- Custom colormaps and color scales

**Documentation:**
- 📘 [Complete Documentation](./docs/plotlib_DOCUMENTATION_COMPLETE.md)
- 📄 [Documentation Summary](./docs/plotlib_DOCUMENTATION_SUMMARY.md)
- 📋 [Quick Reference](./docs/plotlib_QUICK_REFERENCE.md)
- 📖 [Comprehensive Quick Reference](./docs/plotlib_QUICK_REFERENCE_COMPREHENSIVE.md)

---

### 4. projlib.py - Stereographic Projections

Comprehensive stereographic projection utilities for crystallographic data including Schmidt nets, Wulff nets, and inverse pole figures.

**Key Features:**
- Equal-area (Schmidt) and equal-angle (Wulff) projections
- Stereographic triangle for cubic systems
- Inverse pole figure (IPF) coloring
- Pole figure density calculations
- Orientation density plots with contours
- Spherical KDE for texture analysis
- HEALPix grid generation

**Documentation:**
- 📘 [Complete Documentation](./docs/projlib_DOCUMENTATION_COMPLETE.md)
- 📄 [Documentation Summary](./docs/projlib_DOCUMENTATION_SUMMARY.md)
- 📋 [Quick Reference](./docs/projlib_QUICK_REFERENCE.md)
- 📖 [Comprehensive Quick Reference](./docs/projlib_QUICK_REFERENCE_COMPREHENSIVE.md)

---

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone <repository-url>
cd crystallographic-toolkit

# Install dependencies
pip install numpy scipy matplotlib numba
```

### Basic Usage

```python
import numpy as np
import crystlib
import orilib
import plotlib
import projlib

# Example 1: Crystal structure analysis
lattice_b2 = crystlib.cubic_lattice_vec(3.015)  # NiTi B2 phase
lattice_b19p = crystlib.monoclinic_lattice_vec(2.89, 4.12, 4.62, 96.8)  # NiTi B19'
Cd = crystlib.B19p_B2_lattice_correspondence(lattice_b2, lattice_b19p)

# Example 2: Orientation analysis
euler_angles = np.array([45, 30, 60]) * np.pi/180  # radians
rotation_matrix = orilib.np_euler_matrix(*euler_angles)
axis, angle = orilib.np_ol_g_rtheta_rad(rotation_matrix)

# Example 3: Stereographic projection
directions = [[1,0,0], [0,1,0], [0,0,1]]
x, y = projlib.equalarea_directions(directions, np.eye(3))

# Example 4: Visualization
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Plot lattice structure
plotlib.set_aspect_equal_3d(ax)
plt.show()
```

---

## 📖 Documentation Guide

Each module has four levels of documentation:

### 1. Complete Documentation (`*_DOCUMENTATION_COMPLETE.md`)
- Comprehensive reference for all functions
- Detailed mathematical formulas and algorithms
- Extensive usage examples for each function
- Input/output specifications
- Implementation notes and warnings

### 2. Documentation Summary (`*_DOCUMENTATION_SUMMARY.md`)
- Overview of all functions grouped by category
- Brief descriptions and primary use cases
- Quick parameter reference
- Cross-references between related functions

### 3. Quick Reference (`*_QUICK_REFERENCE.md`)
- Compact function listing with signatures
- One-line descriptions
- Common usage patterns
- Typical parameter values
- Perfect for quick lookups

### 4. Comprehensive Quick Reference (`*_QUICK_REFERENCE_COMPREHENSIVE.md`)
- Expanded quick reference with more details
- Example code snippets for each function
- Parameter ranges and constraints
- Related functions and alternatives

---

## 🔧 Dependencies

### Required
- **numpy** (≥1.19.0) - Numerical computing
- **scipy** (≥1.5.0) - Scientific computing (matrix operations)
- **matplotlib** (≥3.3.0) - Plotting and visualization

### Optional
- **numba** (≥0.51.0) - JIT compilation for performance (highly recommended)
- **spherical_kde** - Spherical kernel density estimation
- **sympy** - Symbolic mathematics
- **wand** - Image manipulation (for PDF cropping)

### Installation

```bash
# Minimum installation
pip install numpy scipy matplotlib

# Recommended installation (with performance optimizations)
pip install numpy scipy matplotlib numba

# Full installation (all optional dependencies)
pip install numpy scipy matplotlib numba spherical-kde sympy Wand
```

---

## 📂 Project Structure

```
crystallographic-toolkit/
│
├── README.md                           # This file
├── crystlib.py                         # Crystal structure module
├── orilib.py                           # Orientation analysis module
├── plotlib.py                          # Visualization module
├── projlib.py                          # Stereographic projection module
│
└── docs/                               # Documentation directory
    ├── crystlib_DOCUMENTATION_COMPLETE.md
    ├── crystlib_DOCUMENTATION_SUMMARY.md
    ├── crystlib_QUICK_REFERENCE.md
    ├── crystlib_QUICK_REFERENCE_COMPREHENSIVE.md
    ├── orilib_DOCUMENTATION_COMPLETE.md
    ├── orilib_DOCUMENTATION_SUMMARY.md
    ├── orilib_QUICK_REFERENCE.md
    ├── orilib_QUICK_REFERENCE_COMPREHENSIVE.md
    ├── plotlib_DOCUMENTATION_COMPLETE.md
    ├── plotlib_DOCUMENTATION_SUMMARY.md
    ├── plotlib_QUICK_REFERENCE.md
    ├── plotlib_QUICK_REFERENCE_COMPREHENSIVE.md
    ├── projlib_DOCUMENTATION_COMPLETE.md
    ├── projlib_DOCUMENTATION_SUMMARY.md
    ├── projlib_QUICK_REFERENCE.md
    └── projlib_QUICK_REFERENCE_COMPREHENSIVE.md
```

---

## 🎯 Common Use Cases

### NiTi Shape Memory Alloy Analysis
```python
# Define lattice parameters
a_b2 = 3.015  # B2 cubic
a_b19p, b_b19p, c_b19p, beta = 2.89, 4.12, 4.62, 96.8  # B19' monoclinic

# Create lattices
lattice_b2 = crystlib.cubic_lattice_vec(a_b2)
lattice_b19p = crystlib.monoclinic_lattice_vec(a_b19p, b_b19p, c_b19p, beta)

# Calculate correspondence
Cd = crystlib.B19p_B2_lattice_correspondence(lattice_b2, lattice_b19p)

# Deformation analysis
F, U, Q = crystlib.def_gradient_stressfree(Cd, lattice_b2, lattice_b19p)

# Twinning analysis
twin_data = crystlib.niti_twinning(lattice_b19p)
```

### EBSD Texture Analysis
```python
# Load orientation data (Euler angles)
euler_angles = np.loadtxt('ebsd_data.txt')  # Nx3 array

# Convert to rotation matrices
rotation_matrices = orilib.np_eulers_matrices(euler_angles, deg=True)

# Calculate misorientations
from scipy.spatial.transform import Rotation as R
misorientations = []
for i in range(len(rotation_matrices)-1):
    M1 = rotation_matrices[i]
    M2 = rotation_matrices[i+1]
    # Calculate misorientation...
```

### Pole Figure Generation
```python
import matplotlib.pyplot as plt

# Create pole figure
fig, ax = plt.subplots(figsize=(8, 8))

# Setup equal-area projection
from projlib import schmidtnet_half
fig, ax = schmidtnet_half(ax=ax)

# Plot crystal directions
directions = [[1,0,0], [0,1,0], [0,0,1]]
x, y = projlib.equalarea_directions(directions, np.eye(3))
ax.scatter(x, y, s=100, c='red', marker='o')

plt.title('{100} Pole Figure')
plt.show()
```

---

## 📊 Module Statistics

| Module | Functions | Lines of Code | Size | Primary Focus |
|--------|-----------|---------------|------|---------------|
| crystlib.py | 71 | 7,089 | 275 KB | Crystal structures, lattice operations |
| orilib.py | 60 | 2,500+ | 80 KB | Orientation transformations |
| plotlib.py | 26 | 2,379 | 87 KB | Visualization and plotting |
| projlib.py | 45 | 3,285 | 118 KB | Stereographic projections |
| **Total** | **202** | **15,400+** | **624 KB** | Complete toolkit |

---

## 🔬 Research Applications

This toolkit has been used for:
- **Shape Memory Alloys**: NiTi, Cu-based alloys, phase transformation analysis
- **Martensitic Transformations**: Habit plane calculations, orientation relationships
- **EBSD Analysis**: Texture characterization, grain boundary analysis
- **Crystal Plasticity**: Slip system analysis, deformation modeling
- **Twinning Studies**: Twin variant selection, habit plane determination

---

## 📝 Citation

If you use this toolkit in your research, please cite:

```bibtex
@software{crystallographic_toolkit,
  title = {Crystallographic Analysis Toolkit},
  author = {Heller, L.},
  year = {2019-2024},
  note = {Python toolkit for crystallographic analysis and visualization}
}
```

---

## 👥 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes with tests
4. Submit a pull request

For major changes, please open an issue first to discuss proposed changes.

---

## 📄 License

This project is provided for academic and research use. Please contact the author for licensing information.

---

## 📧 Contact

**Author**: L. Heller  
**Institution**: Materials Science Research  
**Created**: 2019-2024

For questions, bug reports, or feature requests, please open an issue on the repository.

---

## 🔗 Related Resources

- **Crystallography Textbooks**: 
  - Kocks, U. F., Tomé, C. N., & Wenk, H. R. (2000). *Texture and anisotropy*
  - Bunge, H. J. (2013). *Texture analysis in materials science*

- **Software Tools**:
  - MTEX - MATLAB toolbox for texture analysis
  - DREAM.3D - Microstructure analysis platform
  - AIMSGB - Grain boundary analysis

---

## ⚡ Performance Notes

- Functions with `np_` prefix use NumPy vectorization for better performance
- Functions with `@njit` decorator use Numba JIT compilation
- For large datasets (>10,000 orientations), use vectorized functions
- HEALPix sampling is optimized for uniform orientation distributions

---

## 🆕 Version History

**Current Version**: 2024.12 (December 2024)
- Complete redistribution of functions across modules
- Comprehensive documentation for all 186 functions
- Added 4-level documentation structure
- Module-specific organization (crystal, orientation, plot, projection)

**Previous Versions**:
- 2019-2023: Initial development and crystallography_functions.py
- Various improvements and additions for NiTi analysis

---

## 🙏 Acknowledgments

This toolkit was developed for materials science research with focus on:
- Shape memory alloy characterization
- Phase transformation analysis
- Crystallographic texture studies

Special thanks to the materials science community for feedback and contributions.

---

**Last Updated**: December 8, 2024

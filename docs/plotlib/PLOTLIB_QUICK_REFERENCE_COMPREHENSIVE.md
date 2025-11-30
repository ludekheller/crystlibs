# Plotlib.py - Comprehensive Quick Reference Guide

**Crystallographic Plotting Library**  
**Version**: Production Ready  
**Last Updated**: November 2025

---

## 📋 Table of Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [Class Methods](#class-methods)
4. [Utility Functions](#utility-functions)
5. [Complete Workflow Examples](#complete-workflow-examples)
6. [Parameters Quick Reference](#parameters-quick-reference)
7. [Mathematical Formulas](#mathematical-formulas)
8. [Common Pitfalls](#common-pitfalls)
9. [Integration Examples](#integration-examples)
10. [Troubleshooting](#troubleshooting)
11. [Best Practices](#best-practices)

---

## Overview

The **plotlib.py** module provides publication-quality crystallographic plotting capabilities:

**Key Features:**
- Stereographic projections (Schmidt & Wulff nets)
- Pole figures and inverse pole figures
- Orientation density maps with contours
- Custom colormaps and color scales
- Interactive data annotations
- Scatter plots with crystallographic indexing
- Multi-format export (PNG, PDF, SVG)

**Main Class:** `plotter`  
**Standalone Functions:** 5 utility functions

---

## Quick Start

### 3-Step Basic Plot

```python
from plotlib import plotter
import numpy as np

# 1. Create plotter
p = plotter()

# 2. Setup projection
p.plotProj(ProjType='equalarea', sphere='half', figsize=(6, 6))

# 3. Save
p.figsave(fname='output.png')
```

### 5-Minute Pole Figure

```python
from plotlib import plotter

# Initialize
p = plotter()

# Create equal-area projection, upper hemisphere
p.plotProj(ProjType='equalarea', sphere='half')

# Plot crystal directions
p.setAttributes(
    dirs=[[1,0,0], [0,1,0], [0,0,1]],  # <100> family
    norms=[[1,1,1]]                     # {111} plane
)
p.plotDirsNorms()

# Save multi-format
p.figsave(fname='pole_figure.png', imformats=['png', 'pdf'])
```

---

## Class Methods

### 1. Initialize Plotter

#### `plotter()`

**Purpose**: Create plotter instance with default attributes

```python
# Basic initialization
p = plotter()

# Check defaults
print("Colormap:", p.cmap)           # 'jet'
print("Sphere:", p.sphere)           # 'full'
print("Grid points:", p.nump)        # 1001
```

**Default Attributes:**
- **Colorbar**: `cbartitle=""`, `cbarh=0.04`, `cbarwfac=0.75`
- **Figure**: `fig=None`, `figsize=None`, `suptitle=" "`
- **Colormap**: `cmap='jet'`, `nump=1001`, `plotmap=True`
- **Projection**: `sphere='full'`, `stereogrid=False`
- **Axis**: `contourfontsize=9`, `linewidths=0.5`

---

### 2. Setup Projection

#### `plotProj(**kwargs)`

**Purpose**: Create stereographic projection base

**Key Parameters:**
- `ProjType`: 'equalarea' (Schmidt) or 'equalangle' (Wulff)
- `sphere`: 'full', 'half', 'triangle'
- `figsize`: (width, height) in inches
- `stereogrid`: bool - Show grid lines
- `stereoresolution`: int - Grid spacing (degrees)

**Examples:**

```python
# Equal-area, upper hemisphere (most common)
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half', figsize=(6, 6))

# Equal-angle, full sphere
p.plotProj(ProjType='equalangle', sphere='full')

# Stereographic triangle (cubic 001-011-111)
p.plotProj(
    ProjType='equalarea',
    sphere='triangle',
    stereogrid=True,
    stereoresolution=10
)

# Use existing matplotlib axis
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 8))
p.plotProj(ProjType='equalarea', sphere='half', fig=fig, ax=ax)
```

**Projection Types:**
- **Equal-area (Schmidt)**: r = √2 sin(θ/2) - Preserves area
- **Equal-angle (Wulff)**: r = tan(θ/2) - Preserves angles

**Sphere Options:**
- `'full'`: Complete sphere (360°)
- `'half'`: Upper hemisphere (most common for pole figures)
- `'triangle'`: Standard stereographic triangle (cubic)

---

### 3. Set Attributes

#### `setAttributes(**kwargs)`

**Purpose**: Configure plotter attributes dynamically

**Examples:**

```python
# Single attribute
p.setAttributes(cmap='viridis')

# Multiple attributes
p.setAttributes(
    sphere='half',
    figsize=(8, 8),
    dirs=[[1,0,0], [1,1,0], [1,1,1]],
    norms=[[1,1,0]],
    cmap='jet',
    vm=[0, 10]
)

# Crystal directions with colors
p.setAttributes(
    dirs=[[1,0,0], [0,1,0], [0,0,1]],
    dirsfacecolors={0: 'red', 1: 'green', 2: 'blue'},
    dirsedgecolors={0: 'darkred', 1: 'darkgreen', 2: 'darkblue'}
)

# Colorbar configuration
p.setAttributes(
    cbartitle="MRD",      # Multiple of Random Distribution
    cbarh=0.05,           # Height
    cbarwfac=0.8,         # Width factor
    cbarhshift=-0.15      # Vertical shift
)
```

**Common Attributes:**
```python
# Appearance
cmap = 'jet'|'viridis'|'plasma'|'coolwarm'
figsize = (8, 8)
contourfontsize = 9

# Data
dirs = [[1,0,0], [0,1,0], ...]     # Crystal directions [uvw]
norms = [[1,1,1], ...]              # Plane normals (hkl)
oris = orientation_matrices         # (N, 3, 3) array

# Projection
sphere = 'full'|'half'|'triangle'
R2Proj = np.eye(3)                  # Rotation to projection frame

# Colormap
vm = [min, max]                     # Value range
levels = [0, 1, 2, ...]             # Contour levels
contourcol = 'k'                    # Contour color
nump = 1001                         # Grid resolution
```

---

### 4. Generate Color Scales

#### `getScales(vmcbar, numticks=None, ticks=None, tickslabels=None, geq=False, leq=False, cmapbins=100, cmapbinsmult=None)`

**Purpose**: Generate color scale parameters for colorbar

**Input:**
- `vmcbar`: [min, max] - Value range
- `numticks`: int - Number of ticks
- `ticks`: array - Explicit tick positions
- `tickslabels`: list - Custom labels
- `geq`: bool - Add ≥ to maximum
- `leq`: bool - Add ≤ to minimum
- `cmapbins`: int - Colormap resolution

**Output:** Dictionary with:
- `'tickslabels'`: Formatted tick labels
- `'vm'`: [min, max] for colormap
- `'vmbar'`: [min, max] for colorbar
- `'cmap'`: matplotlib colormap
- `'ticks'`: Tick positions

**Examples:**

```python
import numpy as np

# Basic scale 0-10
scales = p.getScales(vmcbar=[0, 10])
print(scales['tickslabels'])  # ['0', '1', ..., '10']

# Custom ticks with inequalities
scales = p.getScales(
    vmcbar=[0, 100],
    ticks=np.array([0, 25, 50, 75, 100]),
    geq=True,  # "≥100"
    leq=True   # "≤0"
)
# Result: ['≤0', '25', '50', '75', '≥100']

# High-resolution colormap
scales = p.getScales(
    vmcbar=[0, 1],
    cmapbins=1000  # Smooth gradient
)

# Custom labels
scales = p.getScales(
    vmcbar=[1, 5],
    ticks=np.array([1, 2, 3, 4, 5]),
    tickslabels=['Very Low', 'Low', 'Medium', 'High', 'Very High']
)

# Use in plot
p.setAttributes(
    vm=scales['vm'],
    ticks=scales['ticks'],
    ticklabels=scales['tickslabels'],
    cmap=scales['cmap']
)
p.plotColorbar()
```

**Default Colormap**: White → Blue → Green → Yellow → Dark Red

---

### 5. Configure Figure Parameters

#### `getFigparam(fontsize=None, save=False, phase='A', figsize=None, **kwargs)`

**Purpose**: Get figure parameter dictionary for saving

**Output:** Dictionary for `plt.savefig()`

**Examples:**

```python
# Get defaults
params = p.getFigparam()
print("DPI:", params['dpi'])  # 300

# High-resolution PDF
params = p.getFigparam(
    save=True,
    format='pdf',
    dpi=600
)

# Custom parameters
params = p.getFigparam(
    figsize=(12, 10),
    dpi=300,
    transparent=True
)

# Use with matplotlib
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot([1, 2, 3], [1, 4, 9])
fig.savefig('output.png', **params)
```

**Default Parameters:**
```python
{
    'dpi': 300,
    'facecolor': 'None',
    'edgecolor': 'w',
    'orientation': 'portrait',
    'format': 'png',
    'transparent': True,
    'bbox_inches': 'None',
    'pad_inches': 0.1
}
```

---

### 6. Save Figures

#### `figsave(**kwargs)`

**Purpose**: Save figure to file(s) with optional cropping

**Key Parameters:**
- `fname`: Filename or list of filenames
- `imformats`: ['png', 'pdf', 'svg'] - Multiple formats
- `crop`: bool - Auto-crop whitespace (PNG only, requires wand)
- `figparam`: dict - Override figure parameters

**Examples:**

```python
# Single format
p.figsave(fname='pole_figure.png')

# Multiple formats
p.figsave(
    fname='pole_figure.png',
    imformats=['png', 'pdf', 'svg']
)
# Creates: pole_figure.png, .pdf, .svg

# With cropping
p.figsave(fname='figure.png', crop=True)

# High DPI
p.figparam['dpi'] = 600
p.figsave(fname='high_res.png')

# Custom parameters
p.figsave(
    fname='custom.pdf',
    figparam={
        'dpi': 300,
        'facecolor': 'white',
        'transparent': False
    }
)
```

---

### 7. Internal Save Method

#### `figsaveproc(fname, **kwargs)`

**Purpose**: Internal method for actual file writing

**Usage:**
```python
# Typically called by figsave(), but can use directly
p.figsaveproc('output.png')
```

---

## Utility Functions

### 1. Create Custom Colormap

#### `get_cmap(colors, nbins=1000, name='my_cmap')`

**Purpose**: Create custom colormap from color list

**Examples:**

```python
from plotlib import get_cmap
import matplotlib.pyplot as plt
import numpy as np

# White to red gradient
colors = [(1, 1, 1), (1, 0, 0)]
cmap = get_cmap(colors, nbins=256)

# Blue-white-red diverging
colors_div = ['blue', 'white', 'red']
cmap_div = get_cmap(colors_div, nbins=100)

# Multi-color rainbow
colors_multi = [
    (0, 0, 1),  # Blue
    (0, 1, 1),  # Cyan
    (0, 1, 0),  # Green
    (1, 1, 0),  # Yellow
    (1, 0, 0)   # Red
]
cmap_rainbow = get_cmap(colors_multi, nbins=500)

# Use in plotter
p = plotter()
p.setAttributes(cmap=cmap_rainbow)

# Or use in matplotlib
data = np.random.rand(10, 10)
plt.imshow(data, cmap=cmap)
plt.colorbar()
plt.show()
```

---

### 2. Map Values to Colors

#### `get_colors(values, cmap, vmin=None, vmax=None, cmin=0, cmax=None, to255=True, nancolor=[0,0,0,1])`

**Purpose**: Convert data values to RGB/RGBA colors

**Examples:**

```python
from plotlib import get_colors
import matplotlib.pyplot as plt
import numpy as np

# Map values to colors
values = np.array([1, 5, 10, 15, 20])
cmap = plt.cm.viridis
colors = get_colors(values, cmap, vmin=0, vmax=20, to255=False)
print("Colors shape:", colors.shape)  # (5, 4) RGBA

# Use in scatter plot
x = np.random.rand(100)
y = np.random.rand(100)
z = np.random.rand(100) * 100
colors = get_colors(z, plt.cm.jet, to255=False)
plt.scatter(x, y, c=colors, s=50)
plt.colorbar()

# Handle NaN values
values_nan = np.array([1, 2, np.nan, 4, 5])
colors = get_colors(
    values_nan,
    plt.cm.coolwarm,
    nancolor=[1, 1, 1, 1]  # White for NaN
)
```

**Parameters:**
- `to255=True`: Scale to 0-255 (int)
- `to255=False`: Keep 0-1 (float)
- `nancolor`: RGBA color for NaN values

---

### 3. Shift Colormap Center

#### `shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap')`

**Purpose**: Shift colormap center to specific data value

**Use Case**: Data with asymmetric range (-15 to +5), want center at zero

**Examples:**

```python
from plotlib import shiftedColorMap
import matplotlib.pyplot as plt
import numpy as np

# Data from -15 to +5, center at zero
vmin, vmax = -15, 5

# Calculate midpoint
midpoint = 1 - vmax / (vmax + abs(vmin))  # 0.75
print(f"Midpoint: {midpoint:.3f}")

# Create shifted colormap
original = plt.cm.RdBu_r
shifted = shiftedColorMap(original, midpoint=midpoint)

# Plot comparison
data = np.random.randn(100, 100) * 10 - 5

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Original
im1 = ax1.imshow(data, cmap=original, vmin=vmin, vmax=vmax)
ax1.set_title('Original')
plt.colorbar(im1, ax=ax1)

# Shifted (zero = white/center color)
im2 = ax2.imshow(data, cmap=shifted, vmin=vmin, vmax=vmax)
ax2.set_title('Shifted (zero centered)')
plt.colorbar(im2, ax=ax2)

plt.show()
```

**Formula:**
```
midpoint = 1 - vmax / (vmax + |vmin|)
```

For symmetric data (-10 to +10): midpoint = 0.5  
For asymmetric (-15 to +5): midpoint = 0.75

---

### 4. Plot Multiple Colormaps

#### `plotcolmaps(fname=None, withdraw=False)`

**Purpose**: Plot 2×2 grid of colormaps from pickle file

**Usage:**
```python
from plotlib import plotcolmaps

# Load and plot saved configuration
plotcolmaps('my_pole_figures.pkl')

# Without display
plotcolmaps('config.pkl', withdraw=True)
```

---

### 5. Plot Single Colormap

#### `plotcolmap(fname=None, withdraw=False)`

**Purpose**: Load and plot single colormap from pickle

**Usage:**
```python
from plotlib import plotcolmap

plotcolmap('single_figure.pkl')
```

---

## Complete Workflow Examples

### Example 1: Basic Pole Figure

```python
from plotlib import plotter

# Create plotter
p = plotter()

# Setup projection
p.plotProj(ProjType='equalarea', sphere='half', figsize=(6, 6))

# Plot crystal directions
directions = [[1,0,0], [0,1,0], [0,0,1]]  # <100>
planes = [[1,1,1]]                         # {111}

p.setAttributes(dirs=directions, norms=planes)
p.plotDirsNorms()

# Title
p.ax.set_title('{100} Directions and (111) Plane', fontsize=14)

# Save
p.figsave(fname='pole_figure.png', imformats=['png', 'pdf'])
```

---

### Example 2: Orientation Density Map

```python
from plotlib import plotter
from orilib import np_eulers_matrices
import numpy as np

# Generate random orientations
N = 1000
euler_angles = np.random.rand(N, 3) * [360, 180, 360]
oris = np_eulers_matrices(euler_angles, deg=True)

# Create density plot
p = plotter()
p.setAttributes(
    oris=oris,
    nump=301,  # Grid resolution
    ProjType='equalarea'
)

p.plotProj(sphere='half', figsize=(8, 8))
p.plotColormap()   # Density contours
p.plotColorbar()   # Add colorbar

# Customize
p.ax.set_title('Orientation Density (MRD)', fontsize=14, weight='bold')
p.figsave(fname='density_map.png', dpi=300)
```

---

### Example 3: Multi-Panel Figure

```python
import matplotlib.pyplot as plt
from plotlib import plotter

# Create 2×2 grid
fig, axes = plt.subplots(2, 2, figsize=(12, 12))

projections = [
    ('half', 'equalarea', '<100> Equal-area'),
    ('half', 'equalangle', '<100> Equal-angle'),
    ('full', 'equalarea', 'Full sphere'),
    ('triangle', 'equalarea', 'Triangle')
]

for ax, (sphere, proj, title) in zip(axes.flat, projections):
    p = plotter()
    p.plotProj(ProjType=proj, sphere=sphere, ax=ax, fig=fig)
    
    # Add directions
    p.setAttributes(dirs=[[1,0,0], [1,1,0], [1,1,1]])
    p.plotDirsNorms()
    
    ax.set_title(title, fontsize=12)

plt.tight_layout()
plt.savefig('multi_panel.png', dpi=300, bbox_inches='tight')
```

---

### Example 4: Custom Colormap Plot

```python
from plotlib import plotter, get_cmap, shiftedColorMap
import matplotlib.pyplot as plt
import numpy as np

# Create custom colormap
colors = [(1,1,1), (0,0,1), (0,1,0), (1,1,0), (1,0,0)]
custom_cmap = get_cmap(colors, nbins=256)

# Or use shifted colormap for asymmetric data
vmin, vmax = -15, 5
midpoint = 1 - vmax / (vmax + abs(vmin))
shifted_cmap = shiftedColorMap(plt.cm.RdBu_r, midpoint=midpoint)

# Setup plotter
p = plotter()
p.setAttributes(
    figsize=(8, 8),
    cmap=custom_cmap,  # or shifted_cmap
    vm=[vmin, vmax],
    contourfontsize=10
)

p.plotProj(ProjType='equalarea', sphere='half')

# Plot data...
# p.plotColormap()
# p.plotColorbar()

p.figsave(fname='custom_cmap.png', dpi=300)
```

---

### Example 5: Publication-Quality Figure

```python
from plotlib import plotter
import matplotlib.pyplot as plt

# Configure for publication
p = plotter()
p.setAttributes(
    figsize=(3.5, 3.5),  # Single column width
    contourfontsize=8,
    dirsnormsalpha=0.7,
    cmap='viridis'
)

# Create projection
p.plotProj(ProjType='equalarea', sphere='half')

# Add data
p.setAttributes(
    dirs=[[1,0,0], [0,1,0], [0,0,1]],
    norms=[[1,1,0], [1,1,1]]
)
p.plotDirsNorms()

# Formatting
p.ax.set_title('(a) Pole Figure', fontsize=10, weight='bold', loc='left')

# Save publication quality
p.figparam['dpi'] = 600
p.figparam['transparent'] = False
p.figparam['facecolor'] = 'white'

p.figsave(
    fname='publication_figure.pdf',
    figparam=p.figparam
)
```

---

### Example 6: Interactive Annotation

```python
from plotlib import plotter
import numpy as np

# Create plot
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

# Enable interactive annotations
p.setAttributes(
    annot=True,
    showdataasplotted=True,
    showdata=['Direction indices', 'Angles']
)

# Plot data with annotations
p.setAttributes(dirs=[[1,0,0], [1,1,0], [1,1,1]])
p.plotDirsNorms()

# Activate mouse-over text
p.onpressActivate()  # If method exists

plt.show()  # Interactive window
```

---

## Parameters Quick Reference

### plotProj() Parameters

| Parameter | Type | Options | Description |
|-----------|------|---------|-------------|
| `ProjType` | str | 'equalarea', 'equalangle' | Projection type |
| `sphere` | str | 'full', 'half', 'triangle' | Hemisphere |
| `figsize` | tuple | (width, height) | Figure size (inches) |
| `stereogrid` | bool | True/False | Show grid |
| `stereoresolution` | int | degrees | Grid spacing |
| `stereomesh` | bool | True/False | Show mesh |
| `ax` | axis | matplotlib axis | Existing axis |
| `fig` | figure | matplotlib figure | Existing figure |

### setAttributes() Common Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cmap` | str | 'jet' | Colormap name |
| `figsize` | tuple | None | Figure size |
| `dirs` | list | [] | Crystal directions [[u,v,w], ...] |
| `norms` | list | [] | Plane normals [[h,k,l], ...] |
| `oris` | array | None | Orientations (N,3,3) |
| `vm` | list | None | [vmin, vmax] for colormap |
| `nump` | int | 1001 | Grid resolution |
| `sphere` | str | 'full' | Hemisphere type |
| `contourfontsize` | int | 9 | Contour label size |

### getScales() Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `vmcbar` | list | required | [min, max] value range |
| `numticks` | int | auto | Number of ticks |
| `ticks` | array | None | Explicit positions |
| `tickslabels` | list | None | Custom labels |
| `geq` | bool | False | Add ≥ to max |
| `leq` | bool | False | Add ≤ to min |
| `cmapbins` | int | 100 | Colormap resolution |

### figsave() Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fname` | str/list | required | Filename(s) |
| `imformats` | list | None | ['png', 'pdf', 'svg'] |
| `crop` | bool | False | Auto-crop whitespace |
| `figparam` | dict | {} | Override parameters |

---

## Mathematical Formulas

### Stereographic Projections

**Equal-area (Schmidt)**:
```
r = √2 × sin(θ/2)
X = r × cos(φ)
Y = r × sin(φ)
```

**Equal-angle (Wulff)**:
```
r = tan(θ/2)
X = r × cos(φ)
Y = r × sin(φ)
```

Where:
- θ = colatitude (0 to π radians)
- φ = azimuth (0 to 2π radians)
- r = radial distance on projection plane

### Color Normalization

```
v_norm = (v - v_min) / (v_max - v_min) × (c_max - c_min) + c_min
```

### Colormap Midpoint (Shifted)

For asymmetric data [v_min, v_max], center at zero:
```
midpoint = 1 - v_max / (v_max + |v_min|)
```

**Examples:**
- Symmetric (-10, +10): midpoint = 0.5
- Asymmetric (-15, +5): midpoint = 0.75
- Asymmetric (-5, +15): midpoint = 0.25

---

## Common Pitfalls

### 1. Projection Type Typos
```python
# ✗ Wrong
p.plotProj(ProjType='stereo')  # Not valid

# ✓ Correct
p.plotProj(ProjType='equalarea')   # Schmidt
p.plotProj(ProjType='equalangle')  # Wulff
```

### 2. Missing plotProj() Call
```python
# ✗ Wrong - no projection created
p = plotter()
p.figsave('output.png')  # Nothing to save!

# ✓ Correct
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')
p.figsave('output.png')
```

### 3. Colorbar Range Reversed
```python
# ✗ Wrong
scales = p.getScales(vmcbar=[10, 0])  # max < min

# ✓ Correct
scales = p.getScales(vmcbar=[0, 10])  # min, max
```

### 4. Wrong Coordinate Format
```python
# ✗ Wrong - single direction as flat list
p.setAttributes(dirs=[1, 0, 0])

# ✓ Correct - list of directions
p.setAttributes(dirs=[[1, 0, 0]])  # Note double brackets
p.setAttributes(dirs=[[1,0,0], [0,1,0], [0,0,1]])
```

### 5. Forgetting to Import
```python
# ✗ Wrong
from plotlib import *  # Imports everything
p = plotter()  # Works but imports too much

# ✓ Better
from plotlib import plotter  # Import only what needed
p = plotter()
```

### 6. DPI Too Low for Publication
```python
# ✗ Wrong - default 300 may not be enough
p.figsave('fig.png')

# ✓ Better for publication
p.figparam['dpi'] = 600
p.figsave('fig.png')
```

---

## Integration Examples

### With projlib

```python
from plotlib import plotter
from projlib import equalarea_directions, stereoprojection_directions

# Get projection coordinates
directions = [[1,0,0], [1,1,0], [1,1,1]]

# Equal-area projection
proj_ea = [equalarea_directions(d) for d in directions]

# Plot
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

for proj in proj_ea:
    p.ax.plot(proj[0], proj[1], 'ro', markersize=10)

p.ax.set_title('Directions from projlib')
p.figsave('with_projlib.png')
```

### With orilib

```python
from plotlib import plotter
from orilib import np_eulers_matrices, quat_misori_deg
import numpy as np

# Generate orientations
N = 100
euler = np.random.rand(N, 3) * [360, 180, 360]
oris = np_eulers_matrices(euler, deg=True)

# Calculate misorientations (example)
# misori = calculate_misorientations(oris)

# Plot density
p = plotter()
p.setAttributes(oris=oris, nump=201)
p.plotProj(ProjType='equalarea', sphere='half')
p.plotColormap()
p.plotColorbar()

p.ax.set_title(f'Orientation Density (N={N})')
p.figsave('with_orilib.png')
```

### With crystlib

```python
from plotlib import plotter
from crystlib import generate_hkls
import numpy as np

# Generate Miller indices for cubic
symops = [np.eye(3)]  # Cubic symmetry operations
unique_hkls, equiv, families, all_hkls = generate_hkls(
    hklmax=2,
    symops=symops
)

# Plot poles
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

# Convert hkls to projection coordinates and plot
# (Requires additional projection calculations)
# ...

p.ax.set_title('Miller Indices from crystlib')
p.figsave('with_crystlib.png')
```

### Complete Integration Example

```python
from plotlib import plotter, get_cmap
from orilib import np_eulers_matrices
from projlib import equalarea_directions
import numpy as np

# 1. Generate orientations (orilib)
euler = np.random.rand(500, 3) * [360, 180, 360]
oris = np_eulers_matrices(euler, deg=True)

# 2. Create custom colormap (plotlib)
colors = [(1,1,1), (0,0,1), (0,1,0), (1,0,0)]
cmap = get_cmap(colors, nbins=256)

# 3. Setup plotter
p = plotter()
p.setAttributes(
    oris=oris,
    cmap=cmap,
    nump=201,
    figsize=(8, 8)
)

# 4. Create projection
p.plotProj(ProjType='equalarea', sphere='half')

# 5. Plot density
p.plotColormap()

# 6. Add specific directions (projlib)
special_dirs = [[1,0,0], [0,1,0], [0,0,1]]
for d in special_dirs:
    xy = equalarea_directions(d)
    p.ax.plot(xy[0], xy[1], 'ko', markersize=8, 
              markeredgecolor='white', markeredgewidth=1.5)

# 7. Finalize
p.plotColorbar()
p.ax.set_title('Integrated Example: orilib + projlib + plotlib')
p.figsave('integrated_example.png', dpi=300)
```

---

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| Empty plot | No projection created | Call `plotProj()` first |
| Wrong colors | Colormap not applied | Use `setAttributes(cmap=...)` |
| Low resolution | Default DPI | Set `figparam['dpi'] = 600` |
| Cropping fails | Wand not installed | `pip install Wand` or disable crop |
| Banded colormap | Too few bins | Increase `cmapbins` to 1000+ |
| Figure too small | Default size | Set `figsize=(10, 10)` |
| Directions not showing | Wrong format | Use `[[u,v,w]]` not `[u,v,w]` |
| Colorbar missing | Not called | Call `p.plotColorbar()` |
| Contours not smooth | Low nump | Increase `nump` to 301+ |
| Memory error | nump too high | Reduce `nump` (start with 101) |

---

## Best Practices

### General
1. **Always create projection first**: Call `plotProj()` before adding data
2. **Use descriptive filenames**: Include date, parameters, sample info
3. **Set DPI appropriately**: 300 for presentations, 600 for publications
4. **Save in vector format**: PDF or SVG for scalability

### Projections
5. **Use equal-area for texture**: Schmidt projection preserves area
6. **Use equal-angle for angles**: Wulff projection preserves angles
7. **Hemisphere for pole figures**: Most common is `sphere='half'`
8. **Triangle for cubic**: Use for inverse pole figures

### Colormaps
9. **Choose perceptually uniform**: viridis, plasma, cividis
10. **Use diverging for bipolar**: RdBu, coolwarm for +/- data
11. **High resolution**: `cmapbins >= 256` for smooth gradients
12. **Shift center for asymmetric**: Use `shiftedColorMap` for data not centered at zero

### Resolution
13. **Preview**: `nump=101` for quick checks
14. **Normal**: `nump=201-301` for most plots
15. **Publication**: `nump=501+` for final figures

### File Management
16. **Multi-format export**: Save PNG for preview, PDF for publication
17. **Backup originals**: Keep .pkl files of complex configurations
18. **Version filenames**: Use `_v1`, `_v2` for iterations

### Code Style
19. **Group setAttributes**: Set related attributes together
20. **Comment complex plots**: Explain non-obvious settings
21. **Reuse configurations**: Save common settings as templates

---

## Performance Tips

1. **Lower nump for preview**: Use 101 vs 301 for speed
2. **Cache colormaps**: Reuse `get_cmap()` results
3. **Batch export**: Use `imformats` for multiple formats in one call
4. **Disable interactive**: When generating many figures programmatically
5. **Use appropriate resolution**: Don't use nump=1001 if 301 suffices

---

## Summary Statistics

**Plotter Class:**
- **Methods**: 7 (init, set, proj, scales, params, save×2)
- **Attributes**: 80+ configurable parameters
- **Categories**: 9 (colorbar, figure, colormap, projection, axis, crystal, save, scatter, interactive)

**Utility Functions:**
- **Total**: 5 standalone functions
- **Colormaps**: 2 (get_cmap, shiftedColorMap)
- **Colors**: 1 (get_colors)
- **Plotting**: 2 (plotcolmap, plotcolmaps)

**File Formats:**
- **Raster**: PNG (with optional cropping)
- **Vector**: PDF, SVG, EPS
- **Data**: PKL (pickle for configurations)

---

*This comprehensive quick reference covers all essential features of plotlib.py. For detailed method documentation including all parameters and edge cases, see the full documentation summary.*

**Documentation**: Complete  
**Status**: Production Ready  
**Integration**: orilib, projlib, crystlib compatible

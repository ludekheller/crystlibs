# Plotlib.py - Comprehensive Documentation Summary (Part 1 of 2)

**Crystallographic Plotting Library - Detailed Method Documentation**  
**Version**: Production Ready  
**Last Updated**: November 2025

---

## 📑 Table of Contents - Part 1

1. [Module Overview](#module-overview)
2. [Class: plotter](#class-plotter)
3. [Method 1: \_\_init\_\_](#method-1-__init__)
4. [Method 2: setAttributes](#method-2-setattributes)
5. [Method 3: plotProj](#method-3-plotproj)
6. [Method 4: getScales](#method-4-getscales)

**See Part 2 for**: getFigparam, figsave, figsaveproc, and all utility functions

---

## Module Overview

### Purpose

The **plotlib.py** module provides a comprehensive framework for creating publication-quality crystallographic plots. It specializes in stereographic projections, pole figures, inverse pole figures, and orientation density maps commonly used in materials science and crystallography.

### Key Features

- **Stereographic Projections**: Equal-area (Schmidt) and equal-angle (Wulff) nets
- **Flexible Hemispheres**: Full sphere, upper hemisphere, or stereographic triangle
- **Density Mapping**: Orientation distribution functions with contours
- **Custom Colormaps**: Create and shift colormaps for optimal data visualization
- **Multi-format Export**: PNG, PDF, SVG with automatic cropping
- **Interactive Annotations**: Mouse-over data display
- **Scatter Plots**: Crystallographic indexing and labeling

### Applications

- **Texture Analysis**: Pole figures and ODF plots
- **EBSD Data Visualization**: Inverse pole figures and orientation maps
- **Phase Transformation Studies**: Habit plane variant visualization
- **Crystal Plasticity**: Slip system analysis
- **Materials Characterization**: Anisotropy visualization

### Dependencies

**Required:**
- numpy: Numerical operations
- matplotlib: Base plotting functionality
- projlib: Stereographic projection calculations

**Optional:**
- scipy: Interpolation for smooth contours
- wand: Image cropping (PNG only)
- spherical_kde: Spherical kernel density estimation
- crystallography_functions: Crystallographic utilities

### Module Structure

```
plotlib.py
│
├── class plotter
│   ├── __init__()           # Initialize with defaults
│   ├── setAttributes()      # Dynamic attribute setting
│   ├── plotProj()           # Create stereographic projection
│   ├── getScales()          # Generate color scales
│   ├── getFigparam()        # Figure parameters for saving
│   ├── figsave()            # Save figure(s)
│   └── figsaveproc()        # Internal save method
│
└── Utility Functions
    ├── get_cmap()           # Create custom colormap
    ├── get_colors()         # Map values to colors
    ├── shiftedColorMap()    # Shift colormap center
    ├── plotcolmap()         # Plot single saved config
    └── plotcolmaps()        # Plot multiple saved configs
```

---

## Class: plotter

### Class Description

The `plotter` class is the central interface for all crystallographic plotting operations. It manages figure creation, projection setup, data plotting, colormap configuration, and file export.

**Design Philosophy:**
- Attribute-based configuration (not method chaining)
- Flexible customization through `setAttributes()`
- Separation of concerns (projection, data, styling, export)
- Integration with matplotlib for advanced customization

### Attribute Categories

The plotter maintains 9 categories of attributes:

1. **Colorbar Attributes**
   - `cbartitle`: Colorbar title text
   - `vmbar`: Value range for colorbar
   - `cbarh`: Colorbar height
   - `cbarwfac`: Width factor
   - `cbarhshift`: Vertical shift

2. **Figure Attributes**
   - `fig`: Matplotlib figure object
   - `suptitle`: Figure super title
   - `figsize`: Figure dimensions (width, height)
   - `annotationtext`: Additional annotations
   - `datadeviders`: Data separators for multi-panel

3. **Colormap Attributes**
   - `colmapdata`: Data for colormap
   - `oris`: Orientation matrices (N, 3, 3)
   - `plotmap`: Boolean to plot colormap
   - `contourcol`: Contour line color
   - `nump`: Grid resolution (e.g., 1001, 301, 101)
   - `oris2`: Secondary orientation data

4. **Projection Attributes**
   - `sphere`: 'full', 'half', 'triangle'
   - `R2Proj`: Rotation matrix to projection frame
   - `stereogrid`: Show stereographic grid
   - `stereoresolution`: Grid spacing (degrees)
   - `stereomesh`: Display mesh

5. **Axis Attributes**
   - `ax`: Matplotlib axis object
   - `ticks`: Tick positions
   - `levels`: Contour levels
   - `ticklabels`: Tick label strings
   - `cmap`: Colormap name or object
   - `vm`: Value min/max [vmin, vmax]
   - `contourfontsize`: Font size for contour labels
   - `dirsnormsalpha`: Transparency for dirs/norms
   - `contourlabel`: Show contour labels
   - `linewidths`: Contour line width

6. **Crystal Attributes**
   - `symops`: Direct space symmetry operations
   - `recsymops`: Reciprocal space symmetry operations
   - `dirs`: Crystal directions [[u,v,w], ...]
   - `norms`: Plane normals [[h,k,l], ...]
   - `dhkl`: d-spacing for planes

7. **Save Attributes**
   - `SAVE`: Boolean save flag
   - `fname`: Filename or list of filenames
   - `crop`: Auto-crop whitespace
   - `imformats`: List of output formats
   - `figparam`: Dictionary of figure parameters

8. **Scatter Attributes**
   - `scatterplot`: Boolean to enable scatter
   - `scatterdata`: Data array for scatter
   - `scattercolscale`: Color scale for points
   - `scattercolscalevm`: Value range for colors
   - `scattersizescale`: Size scale for points
   - `scatteridxs`: Indices to plot
   - Additional styling options

9. **Interactive Attributes**
   - `annot`: Enable annotations
   - `onmovetext`: Text display on mouse move
   - `showdataasplotted`: Show data as plotted
   - `showdatanames`: Names for data display
   - `showdata`: Data values to show

---

## Method 1: \_\_init\_\_

### Description

Initialize the plotter instance with default attributes for all plotting categories. Sets up empty dictionaries and default values for 80+ configurable parameters.

**When to use:** Every time you start a new plot.

### Input

**Parameters:** None

### Output

**Returns:** None  
**Side Effects:** Creates plotter instance with initialized attributes

### Initialized Attributes

**Colorbar Defaults:**
```python
cbartitle = ""           # Empty title
vmbar = None             # Auto-determine range
cbarh = 0.04             # Height ratio
cbarwfac = 0.75          # Width factor
cbarhshift = -0.15       # Vertical position shift
```

**Figure Defaults:**
```python
fig = None               # No figure yet
suptitle = " "           # Empty super title
figsize = None           # Auto-determined
annotationtext = ""      # No annotations
```

**Colormap Defaults:**
```python
colmapdata = None        # No data yet
oris = None              # No orientations
plotmap = True           # Plot colormap by default
contourcol = 'k'         # Black contour lines
nump = 1001              # High resolution grid
```

**Projection Defaults:**
```python
sphere = 'full'          # Complete sphere
R2Proj = np.eye(3)       # Identity rotation
stereogrid = False       # No grid
stereoresolution = None  # Auto-determined
stereomesh = False       # No mesh
```

**Axis Defaults:**
```python
ax = None                # No axis yet
ticks = None             # Auto ticks
levels = None            # Auto levels
ticklabels = None        # Auto labels
cmap = 'jet'             # Default colormap
vm = None                # Auto value range
contourfontsize = 9      # Contour label size
dirsnormsalpha = 0.5     # Semi-transparent
contourlabel = True      # Show labels
linewidths = 0.5         # Thin contour lines
```

**Crystal Defaults:**
```python
symops = [np.eye(3)]     # Identity symmetry
recsymops = [np.eye(3)]  # Identity reciprocal
dirs = []                # No directions
norms = []               # No normals
```

**Save Defaults:**
```python
SAVE = False             # Don't auto-save
fname = r".\EA.png"      # Default filename
crop = False             # No cropping
imformats = None         # Single format
figparam = {             # Figure parameters
    'dpi': 300,
    'facecolor': 'none',
    'edgecolor': 'none',
    'orientation': 'portrait',
    'transparent': 'False',
    'bbox_inches': 'None',
    'pad_inches': 0.1
}
```

**Scatter Defaults:**
```python
scatterplot = False      # No scatter plot
scatterdata = None       # No data
# ... many more scatter options
```

**Interactive Defaults:**
```python
annot = False            # No annotations
onmovetext = None        # No mouse-over text
showdataasplotted = True # Show as plotted
```

### Usage Examples

#### Example 1: Basic Initialization

```python
from plotlib import plotter

# Create plotter
p = plotter()

# Check default values
print("Default colormap:", p.cmap)           # 'jet'
print("Default sphere:", p.sphere)           # 'full'
print("Default resolution:", p.nump)         # 1001
print("Grid enabled:", p.stereogrid)         # False
print("Contour font size:", p.contourfontsize)  # 9
```

#### Example 2: Initialization with Immediate Configuration

```python
from plotlib import plotter

# Initialize
p = plotter()

# Immediately configure
p.setAttributes(
    cmap='viridis',
    sphere='half',
    nump=301,
    figsize=(8, 8),
    contourfontsize=10
)

# Now ready to plot
p.plotProj(ProjType='equalarea')
```

#### Example 3: Multiple Plotter Instances

```python
from plotlib import plotter

# Create multiple independent plotters
p1 = plotter()  # For pole figure
p2 = plotter()  # For inverse pole figure
p3 = plotter()  # For density map

# Configure each independently
p1.setAttributes(sphere='half', cmap='jet')
p2.setAttributes(sphere='triangle', cmap='viridis')
p3.setAttributes(sphere='half', cmap='plasma', nump=501)

# Each has independent settings
print(p1.cmap, p2.cmap, p3.cmap)  # jet viridis plasma
```

### Notes

- All attributes can be changed later using `setAttributes()`
- Default `nump=1001` provides high resolution but slower computation
- Default `cmap='jet'` is common but not perceptually uniform
- `sphere='full'` shows both hemispheres (less common for pole figures)
- Multiple plotter instances are completely independent

---

## Method 2: setAttributes

### Description

Set or update any plotter attribute dynamically. This is the primary method for configuring the plotter after initialization. Accepts any keyword argument and updates the instance's `__dict__`.

**When to use:** 
- After initialization to configure plotting
- To change settings between plots
- To set multiple related attributes together

### Input

**Parameters:**
- `**kwargs`: Any attribute name and value pairs

**Common Keyword Arguments:**
```python
# Appearance
cmap='jet'|'viridis'|'plasma'|'coolwarm'
figsize=(width, height)
contourfontsize=int

# Data
dirs=[[u,v,w], ...]          # List of directions
norms=[[h,k,l], ...]         # List of plane normals
oris=array(N, 3, 3)          # Orientation matrices

# Projection
sphere='full'|'half'|'triangle'
ProjType='equalarea'|'equalangle'
R2Proj=array(3, 3)           # Rotation matrix

# Colormap
vm=[vmin, vmax]              # Value range
levels=[...]                 # Contour levels
contourcol='k'|'w'|'r'       # Contour color
nump=101|201|301|501|1001    # Grid resolution

# Colorbar
cbartitle='string'
cbarh=float                  # Height
cbarwfac=float              # Width factor

# Save
fname='filename.png'
crop=True|False
imformats=['png', 'pdf', 'svg']
```

### Output

**Returns:** None  
**Side Effects:** Updates instance attributes

### Usage Examples

#### Example 1: Set Single Attribute

```python
from plotlib import plotter

p = plotter()

# Change colormap
p.setAttributes(cmap='viridis')

# Change grid resolution
p.setAttributes(nump=301)  # Medium resolution

# Change figure size
p.setAttributes(figsize=(10, 10))
```

#### Example 2: Set Multiple Related Attributes

```python
from plotlib import plotter
import numpy as np

p = plotter()

# Configure projection and appearance together
p.setAttributes(
    sphere='half',
    ProjType='equalarea',
    figsize=(8, 8),
    stereogrid=True,
    stereoresolution=10
)

# Configure crystal data
p.setAttributes(
    dirs=[[1,0,0], [0,1,0], [0,0,1]],
    norms=[[1,1,1], [1,1,0]],
    symops=[np.eye(3)]  # Cubic symmetry
)

# Configure colormap
p.setAttributes(
    cmap='jet',
    vm=[0, 10],
    levels=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    contourcol='k',
    contourfontsize=9
)
```

#### Example 3: Configure for Publication

```python
from plotlib import plotter

p = plotter()

# Publication-quality settings
p.setAttributes(
    figsize=(3.5, 3.5),      # Single column width
    cmap='viridis',           # Perceptually uniform
    contourfontsize=8,        # Small labels
    dirsnormsalpha=0.7,       # Semi-transparent markers
    nump=501                  # High resolution
)

# Update figure parameters for high DPI
p.figparam['dpi'] = 600
p.figparam['transparent'] = False
p.figparam['facecolor'] = 'white'
```

#### Example 4: Set Crystal Directions with Colors

```python
from plotlib import plotter

p = plotter()

# Set directions
directions = [[1,0,0], [0,1,0], [0,0,1]]

p.setAttributes(
    dirs=directions,
    dirsfacecolors={0: 'red', 1: 'green', 2: 'blue'},
    dirsedgecolors={0: 'darkred', 1: 'darkgreen', 2: 'darkblue'}
)

# Each direction gets a specific color
# Index 0 ([1,0,0]) → red
# Index 1 ([0,1,0]) → green  
# Index 2 ([0,0,1]) → blue
```

#### Example 5: Dynamic Configuration in Loop

```python
from plotlib import plotter
import numpy as np

# Create multiple plots with different settings
cmaps = ['jet', 'viridis', 'plasma', 'coolwarm']

for i, cmap_name in enumerate(cmaps):
    p = plotter()
    
    # Dynamic configuration
    p.setAttributes(
        cmap=cmap_name,
        sphere='half',
        figsize=(6, 6),
        nump=201
    )
    
    p.plotProj(ProjType='equalarea')
    # ... add data ...
    p.figsave(fname=f'plot_{i}_{cmap_name}.png')
```

### Notes

- Can set any attribute, even ones not in `__init__`
- Changes take effect immediately
- Multiple calls accumulate (don't reset other attributes)
- Can override attributes set in previous calls
- Dictionaries for colors use index as key (0-based)

---

## Method 3: plotProj

### Description

Create a stereographic projection (pole figure) with appropriate grid. This method sets up the figure, axes, and projection type (equal-area Schmidt net or equal-angle Wulff net). It handles full sphere, hemisphere, or standard stereographic triangle projections.

**When to use:** 
- First step after initializing plotter
- Before plotting any data
- To create the base projection framework

### Input

**Parameters:**
- `**kwargs`: Keyword arguments for projection setup

**Key Keyword Arguments:**
```python
ProjType : str
    'equalarea' - Schmidt net (preserves area)
    'equalangle' - Wulff net (preserves angles)

sphere : str
    'full' - Complete sphere (two circles)
    'half' - Upper hemisphere only
    'triangle' - Stereographic triangle (cubic)

figsize : tuple
    (width, height) in inches, e.g., (6, 6), (8, 8)

stereogrid : bool
    True - Show stereographic grid lines
    False - No grid (default)

stereoresolution : int
    Grid line spacing in degrees (e.g., 10, 15, 30)

stereomesh : bool
    True - Display mesh
    False - No mesh (default)

ax : matplotlib.axes.Axes
    Existing axis to plot on (optional)

fig : matplotlib.figure.Figure
    Existing figure (optional)
```

### Output

**Returns:** None  
**Side Effects:** 
- Creates or modifies `self.fig` (matplotlib figure)
- Creates or modifies `self.ax` (matplotlib axis)
- Sets `self.equalarea` boolean flag
- Sets label offset parameters (`dx1`, `dy1`, `dx2`, `dy2`)

### Mathematical Background

**Equal-area (Schmidt) Projection:**
```
r = √2 × sin(θ/2)
X = r × cos(φ)
Y = r × sin(φ)
```
- Preserves area on the sphere
- Best for density distributions and texture analysis
- Circles on sphere → ellipses on projection
- Used for pole density plots

**Equal-angle (Wulff) Projection:**
```
r = tan(θ/2)
X = r × cos(φ)
Y = r × sin(φ)
```
- Preserves angles between great circles
- Best for angular relationships
- Circles on sphere → circles on projection
- Used for crystallographic geometry

Where:
- θ = colatitude (0 at north pole, π at south pole)
- φ = azimuthal angle (0 to 2π)
- r = radial distance on projection plane

### Usage Examples

#### Example 1: Basic Equal-area Projection

```python
from plotlib import plotter
import matplotlib.pyplot as plt

# Most common: equal-area, upper hemisphere
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half', figsize=(6, 6))

# Add title
p.ax.set_title('Equal-area Projection (Schmidt Net)')

plt.show()
```

#### Example 2: Equal-angle Projection

```python
from plotlib import plotter

# Equal-angle (Wulff) projection
p = plotter()
p.plotProj(ProjType='equalangle', sphere='half', figsize=(6, 6))

# Better for angular measurements
p.ax.set_title('Equal-angle Projection (Wulff Net)')

p.figsave('wulff_net.png')
```

#### Example 3: Full Sphere Projection

```python
from plotlib import plotter

# Both hemispheres
p = plotter()
p.plotProj(ProjType='equalarea', sphere='full', figsize=(8, 6))

# Shows upper and lower hemispheres as two circles
p.ax.set_title('Full Sphere Projection')

p.figsave('full_sphere.png')
```

#### Example 4: Stereographic Triangle

```python
from plotlib import plotter

# Standard stereographic triangle for cubic crystals
p = plotter()
p.plotProj(
    ProjType='equalarea',
    sphere='triangle',
    figsize=(6, 6),
    stereogrid=True,          # Show grid
    stereoresolution=10       # 10° spacing
)

# Triangle bounded by 001, 011, 111 for cubic
p.ax.set_title('Stereographic Triangle (Cubic)')

p.figsave('triangle.png')
```

#### Example 5: Using Existing Figure

```python
import matplotlib.pyplot as plt
from plotlib import plotter

# Create figure with specific properties
fig, ax = plt.subplots(figsize=(8, 8), facecolor='white')

# Use existing axis
p = plotter()
p.plotProj(
    ProjType='equalarea',
    sphere='half',
    fig=fig,
    ax=ax
)

# Can add more to the same figure
ax.set_title('Using Existing Figure')

plt.savefig('existing_fig.png', dpi=300)
```

#### Example 6: Multi-panel with Different Projections

```python
import matplotlib.pyplot as plt
from plotlib import plotter

# Create 1×2 subplot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Equal-area on left
p1 = plotter()
p1.plotProj(ProjType='equalarea', sphere='half', fig=fig, ax=ax1)
ax1.set_title('Equal-area (Schmidt)')

# Equal-angle on right
p2 = plotter()
p2.plotProj(ProjType='equalangle', sphere='half', fig=fig, ax=ax2)
ax2.set_title('Equal-angle (Wulff)')

plt.tight_layout()
plt.savefig('comparison.png', dpi=300)
```

### Notes

- **Order matters**: Call `plotProj()` before plotting data
- **Equal-area most common**: Used for ~90% of pole figures
- **Triangle for inverse pole figures**: Shows fundamental zone for cubic
- **Grid optional**: Usually omitted for cleaner appearance
- **Resolution trade-off**: Higher resolution = slower but smoother contours

**Default Label Offsets (set automatically):**
- Full sphere equal-area: `dy1=-0.12`, `dx1=-0.08`, `dy2=0.02`, `dx2=-0.05`
- Half sphere: `dy1=-0.14`, `dx1=-0.08`, `dy2=0.02`, `dx2=-0.05`
- Triangle equal-area: `dy1=-0.04`, `dx1=0.0`, `dy2=0.01`, `dx2=0.01`

---

## Method 4: getScales

### Description

Generate color scale parameters including ticks, labels, and colormap. Creates a custom gradient colormap from white → blue → green → yellow → dark red with customizable tick positions and labels. Essential for creating colorbars with proper scales.

**When to use:**
- Before plotting colormap/colorbar
- To customize colorbar appearance
- To add inequality symbols (≥, ≤)
- To create custom tick labels

### Input

**Parameters:**

```python
vmcbar : list or array
    [min, max] - Value range for colorbar
    Example: [0, 10], [-5, 5], [0.5, 2.5]

numticks : int, optional
    Number of tick marks
    Default: None (auto-computed as max - min + 1)

ticks : numpy array, optional
    Explicit tick positions
    Example: np.array([0, 2.5, 5, 7.5, 10])
    Default: None (uses linspace)

tickslabels : list of str, optional
    Custom tick labels
    Example: ['Low', 'Medium', 'High']
    Default: None (uses numeric values)

geq : bool, optional
    Add ≥ symbol to maximum tick
    Default: False

leq : bool, optional
    Add ≤ symbol to minimum tick
    Default: False

cmapbins : int, optional
    Number of colormap bins (higher = smoother)
    Default: 100
    Recommended: 256-1000 for publication

cmapbinsmult : int, optional
    Multiplier for bins based on number of ticks
    If set, cmapbins = (numticks - 1) * cmapbinsmult
    Default: None
```

### Output

**Returns:** Dictionary with keys:

```python
{
    'tickslabels': list of str
        Formatted tick labels (with ≥, ≤ if requested)
    
    'vm': list [min, max]
        Colormap value range
    
    'vmbar': list [min, max]
        Colorbar value range (same as vmcbar)
    
    'cmap': matplotlib.colors.LinearSegmentedColormap
        Generated colormap object
    
    'ticks': numpy array
        Tick positions
}
```

### Colormap Gradient

The default colormap uses 5 key colors:
1. White (1, 1, 1)
2. Blue (0, 0, 1)
3. Green (0, 1, 0)
4. Yellow (1, 1, 0)
5. Dark Red (0.545, 0, 0)

This creates a smooth gradient: White → Blue → Green → Yellow → Dark Red

### Usage Examples

#### Example 1: Basic Scale

```python
from plotlib import plotter
import numpy as np

p = plotter()

# Basic scale from 0 to 10
scales = p.getScales(vmcbar=[0, 10])

print("Tick labels:", scales['tickslabels'])
# Output: ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

print("Colormap range:", scales['vm'])
# Output: [0, 10]

print("Colormap type:", type(scales['cmap']))
# Output: <class 'matplotlib.colors.LinearSegmentedColormap'>

# Use in plot
p.setAttributes(
    vm=scales['vm'],
    ticks=scales['ticks'],
    ticklabels=scales['tickslabels'],
    cmap=scales['cmap']
)
```

#### Example 2: Custom Ticks with Inequalities

```python
from plotlib import plotter
import numpy as np

p = plotter()

# MRD scale 0-100 with custom ticks
scales = p.getScales(
    vmcbar=[0, 100],
    ticks=np.array([0, 25, 50, 75, 100]),
    geq=True,  # Add ≥ to maximum
    leq=True   # Add ≤ to minimum
)

print("Labels:", scales['tickslabels'])
# Output: ['≤0', '25', '50', '75', '≥100']

# This indicates:
# - Values ≤ 0 get minimum color
# - Values ≥ 100 get maximum color
# - Values between interpolate
```

#### Example 3: High-Resolution Smooth Colormap

```python
from plotlib import plotter

p = plotter()

# Very smooth gradient for publication
scales = p.getScales(
    vmcbar=[0, 1],
    numticks=11,      # 0.0, 0.1, 0.2, ..., 1.0
    cmapbins=1000     # High resolution
)

print("Number of bins:", 1000)
# Smooth gradients with no visible banding
```

#### Example 4: Custom Tick Labels

```python
from plotlib import plotter
import numpy as np

p = plotter()

# Descriptive labels instead of numbers
scales = p.getScales(
    vmcbar=[1, 5],
    ticks=np.array([1, 2, 3, 4, 5]),
    tickslabels=['Very Low', 'Low', 'Medium', 'High', 'Very High']
)

print("Custom labels:", scales['tickslabels'])
# Output: ['Very Low', 'Low', 'Medium', 'High', 'Very High']
```

#### Example 5: Fractional Values

```python
from plotlib import plotter
import numpy as np

p = plotter()

# Non-integer ticks
scales = p.getScales(
    vmcbar=[0.0, 2.5],
    ticks=np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
)

print("Fractional labels:", scales['tickslabels'])
# Output: ['0.0', '0.5', '1.0', '1.5', '2.0', '2.5']

# For cleaner labels with no decimals:
scales2 = p.getScales(
    vmcbar=[0, 10],
    ticks=np.array([0, 2, 4, 6, 8, 10])
)
print("Integer labels:", scales2['tickslabels'])
# Output: ['0', '2', '4', '6', '8', '10']
```

#### Example 6: Automatic Bins Based on Ticks

```python
from plotlib import plotter
import numpy as np

p = plotter()

# Bins automatically scaled to number of ticks
scales = p.getScales(
    vmcbar=[0, 10],
    ticks=np.array([0, 2, 4, 6, 8, 10]),  # 6 ticks
    cmapbinsmult=20  # 20 bins per interval
)

# Effective bins = (6 - 1) * 20 = 100
# Creates 100 bins for 5 intervals
```

### Formula for Tick Label Formatting

```python
# For each tick value:
if tick == round(tick):
    label = str(round(tick))  # Integer: "5"
else:
    label = str(tick)          # Float: "5.5"

# Then apply inequalities:
if geq:
    labels[-1] = '≥' + labels[-1]  # Maximum
if leq:
    labels[0] = '≤' + labels[0]     # Minimum
```

### Notes

- **Always sorted**: `vmcbar` is automatically sorted to [min, max]
- **Inequality symbols**: Use `geq`/`leq` when data extends beyond range
- **Bin resolution**: 100 bins = visible banding, 256+ = smooth gradient
- **Default auto-ticks**: Creates `max - min + 1` ticks if not specified
- **Colormap reusable**: Can store and reuse `scales['cmap']` object

**Common Patterns:**
- Texture data (MRD): `vmcbar=[0, max_MRD]`, `geq=True`
- Strain data: `vmcbar=[min_strain, max_strain]`, both `geq` and `leq`
- Normalized data: `vmcbar=[0, 1]`, `cmapbins=256`

---

*This completes Part 1 of the plotlib documentation. Part 1 covers the class overview and the first 4 methods (__init__, setAttributes, plotProj, getScales). See Part 2 for the remaining methods (getFigparam, figsave, figsaveproc) and all utility functions.*

**Documentation Status**: Part 1 Complete  
**Methods Covered**: 4 of 7 class methods  
**Next**: Part 2 - Remaining methods and utilities

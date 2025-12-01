# Plotlib.py Documentation Summary

## Overview

The `plotlib.py` module provides a comprehensive framework for creating publication-quality crystallographic plots including stereographic projections, pole figures, and orientation density maps.

## Documentation Structure

### Original Files
- **plotlib.py**: Original source code (uncommented)
- **plotlib_commented.py**: Fully documented version with comprehensive docstrings

### Documentation Files
- **plotter_class_documentation.pdf**: Complete LaTeX-generated reference manual
- **plotter_class_documentation.tex**: LaTeX source for the manual
- **plotlib_quick_reference.md**: Quick reference guide (this companion file)
- **plotlib_documentation_summary.md**: This document

## Main Components

### 1. The `plotter` Class

The central class for all plotting operations, providing:
- Stereographic projection setup (Schmidt/Wulff nets)
- Pole figure creation with density contours
- Color mapping and colorbar generation
- Crystal direction and plane visualization
- Interactive data annotations
- Multi-format figure export

### 2. Standalone Utility Functions

Helper functions for color manipulation and plotting:
- `get_cmap()`: Create custom colormaps
- `get_colors()`: Map data values to colors
- `shiftedColorMap()`: Shift colormap center for asymmetric data
- `plotcolmaps()`, `plotcolmap()`: Load and display saved configurations

## Documented Methods

### Core Methods

#### `__init__()`
**Purpose**: Initialize plotter with default attributes

**Features**:
- Sets up all attribute categories (colorbar, figure, projection, etc.)
- Provides sensible defaults for all parameters
- Ready to use immediately or customize via `setAttributes()`

**Attribute Categories**:
1. **Colorbar**: Color scale display and formatting
2. **Figure**: Layout, titles, annotations
3. **Colormap**: Data-to-color mapping, contours
4. **Projection**: Stereographic type and settings
5. **Axis**: Plot styling, ticks, labels
6. **Crystal**: Symmetry operations, directions
7. **Save**: File output configuration
8. **Scatter**: Scatter plot data and styling
9. **Interactive**: Mouse-over annotations

#### `setAttributes(**kwargs)`
**Purpose**: Dynamically set or update any plotter attribute

**Input**:
- `**kwargs`: Any attribute name-value pairs

**Features**:
- Flexible configuration without knowing internal structure
- Can set single or multiple attributes
- Updates bracket notation based on projection type
- Validates attribute updates

**Usage**:
```python
# Single attribute
p.setAttributes(cmap='viridis')

# Multiple attributes
p.setAttributes(
    sphere='half',
    ProjType='equalarea',
    figsize=(8, 8),
    dirs=[[1,0,0], [1,1,0]]
)
```

#### `plotProj(**kwargs)`
**Purpose**: Create stereographic projection with appropriate grid

**Input Parameters**:
- `ProjType`: 'equalarea' (Schmidt) or 'equalangle' (Wulff)
- `sphere`: 'full', 'half', or 'triangle'
- `figsize`: Figure size (width, height) in inches
- `stereogrid`: Display stereographic grid lines
- `stereoresolution`: Grid line resolution
- `ax`, `fig`: Existing matplotlib objects (optional)

**Output**:
- Creates/modifies `self.fig` and `self.ax`
- Sets up projection-specific parameters

**Projection Types**:
1. **Equal-area (Schmidt)**: Preserves area, ideal for texture analysis
2. **Equal-angle (Wulff)**: Preserves angles, ideal for crystallography
3. **Triangle**: Standard stereographic triangle for cubic systems

#### `getScales(vmcbar, ...)`
**Purpose**: Generate color scale parameters for colorbars

**Input Parameters**:
- `vmcbar`: Value range [min, max]
- `numticks`: Number of tick marks
- `ticks`: Explicit tick positions
- `tickslabels`: Custom tick labels
- `geq`, `leq`: Add inequality symbols (≥, ≤)
- `cmapbins`: Colormap resolution (default: 100)

**Output Dictionary**:
- `'tickslabels'`: Formatted tick labels
- `'vm'`: Colormap value range [min, max]
- `'vmbar'`: Colorbar value range [min, max]
- `'cmap'`: Generated colormap object
- `'ticks'`: Tick positions array

**Color Gradient**:
- White → Blue → Green → Yellow → Dark Red
- Customizable number of bins for smooth transitions

#### `getFigparam(...)`
**Purpose**: Get figure parameters for high-quality output

**Input Parameters**:
- `fontsize`: Text font size
- `save`: Use save-optimized size
- `phase`: Phase identifier (affects default size)
- `figsize`: Custom figure size
- `**kwargs`: Override any default parameter

**Output Dictionary**:
- `'dpi'`: Resolution (default: 300)
- `'facecolor'`: Background color
- `'edgecolor'`: Edge color
- `'orientation'`: 'portrait' or 'landscape'
- `'format'`: Image format
- `'transparent'`: Transparent background
- `'bbox_inches'`: Bounding box setting
- `'pad_inches'`: Padding around figure

#### `figsave(**kwargs)`
**Purpose**: Save figure to file(s) with optional cropping

**Input Parameters**:
- `fname`: Filename (string or list)
- `imformats`: List of formats ['png', 'pdf', 'svg']
- `crop`: Auto-crop whitespace (requires wand)
- `figparam`: Figure parameters dictionary

**Features**:
- Multi-format export in single call
- Automatic whitespace cropping
- High-resolution output support
- Vector and raster format support

**Supported Formats**:
- Raster: PNG, JPG, TIFF
- Vector: PDF, SVG, EPS
- Other: PS, SVGZ

#### `figsaveproc(fname, **kwargs)`
**Purpose**: Internal method for actual file writing

**Called by**: `figsave()`

**Features**:
- Handles format detection from filename
- Applies tight_layout if requested
- Uses matplotlib's savefig with configured parameters

## Standalone Utility Functions

### `get_cmap(colors, nbins=1000, name='my_cmap')`
**Purpose**: Create custom colormap from color list

**Input**:
- `colors`: List of RGB tuples or color names
- `nbins`: Number of discrete bins (default: 1000)
- `name`: Colormap name (default: 'my_cmap')

**Output**:
- LinearSegmentedColormap object

**Usage**:
```python
colors = [(1, 1, 1), (1, 0, 0), (0, 0, 0)]  # White→Red→Black
cmap = get_cmap(colors, nbins=256)
```

### `get_colors(values, cmap, vmin=None, vmax=None, ...)`
**Purpose**: Map data values to colors using colormap

**Input**:
- `values`: Data values array
- `cmap`: Matplotlib colormap
- `vmin`, `vmax`: Value range (auto if None)
- `cmin`, `cmax`: Colormap index range
- `to255`: Scale to 0-255 range
- `nancolor`: RGBA color for NaN values

**Output**:
- Colors array (RGBA or RGB)

**Usage**:
```python
values = np.array([1, 5, 10, 15, 20])
colors = get_colors(values, plt.cm.viridis, vmin=0, vmax=20)
```

### `shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, ...)`
**Purpose**: Shift colormap center to specific data value

**Input**:
- `cmap`: Original colormap
- `start`: Offset from lowest point [0.0, midpoint]
- `midpoint`: New center position [0.0, 1.0]
- `stop`: Offset from highest point [midpoint, 1.0]
- `name`: New colormap name

**Output**:
- Shifted LinearSegmentedColormap

**Key Formula**:
```
midpoint = 1 - v_max/(v_max + |v_min|)
```

**Usage**:
```python
# Data: -15 to +5, want zero = white
midpoint = 1 - 5/(5 + 15)  # 0.75
shifted_cmap = shiftedColorMap(plt.cm.RdBu_r, midpoint=0.75)
```

### `plotcolmaps(fname=None, withdraw=False)`
**Purpose**: Plot multiple colormaps from saved pickle file

**Input**:
- `fname`: Path to pickle file
- `withdraw`: Withdraw plot from display

**Features**:
- Creates 2x2 subplot figure
- Loads plot configurations from pickle
- Displays 4 pole figures simultaneously

### `plotcolmap(fname=None, withdraw=False)`
**Purpose**: Plot single colormap from saved pickle file

**Input**:
- `fname`: Path to pickle file
- `withdraw`: Withdraw plot from display

**Features**:
- Loads and executes plot configuration
- Restores figure size from saved data

## Complete Usage Examples

### Example 1: Basic Pole Figure
```python
import numpy as np
from plotlib import plotter

# Create and configure
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')
p.setAttributes(
    dirs=[[1,0,0], [0,1,0], [0,0,1]],
    norms=[[1,1,1]]
)

# Save
p.figsave(fname='pole_figure.png')
```

### Example 2: Density Contour Plot
```python
import numpy as np
from plotlib import plotter

orientations = np.random.randn(1000, 3, 3)
p = plotter()
p.setAttributes(oris=orientations, ProjType='equalarea')
p.plotProj()
p.figsave(fname='density_plot.png', imformats=['png', 'pdf'])
```

### Example 3: Multi-Panel Figure
```python
import matplotlib.pyplot as plt
from plotlib import plotter

fig, axes = plt.subplots(2, 2, figsize=(12, 12))

for i, ax in enumerate(axes.flat):
    p = plotter()
    p.plotProj(ProjType='equalarea', sphere='half', ax=ax, fig=fig)
    dirs = [[i+1, 0, 0], [0, i+1, 0], [0, 0, i+1]]
    p.setAttributes(dirs=dirs)

plt.tight_layout()
plt.savefig('multi_panel.png', dpi=300)
```

### Example 4: High-Quality Publication Figure
```python
from plotlib import plotter

p = plotter()
p.setAttributes(figsize=(8, 8))
p.plotProj(ProjType='equalarea', sphere='half')

# High-resolution settings
p.figparam['dpi'] = 600
p.setAttributes(tight_layout=True)

# Save multiple formats
p.figsave(
    fname='publication.png',
    imformats=['png', 'pdf', 'svg'],
    crop=True
)
```

## Best Practices

### Figure Quality
1. **DPI Settings**:
   - Preview: 150 DPI
   - Screen: 300 DPI  
   - Print: 600 DPI
   - Default: 300 DPI

2. **File Formats**:
   - Vector (PDF, SVG): Scalable, best for simple plots
   - Raster (PNG): Good for complex plots, photos
   - Multi-format: Export all at once

3. **Layout**:
   - Use `tight_layout=True` for multi-panel figures
   - Enable `crop=True` to remove whitespace
   - Set appropriate `figsize` before plotting

### Color Scales
1. **Colormap Selection**:
   - Continuous data: viridis, plasma, inferno
   - Diverging data: RdBu, RdYlGn, coolwarm
   - Sequential data: YlOrRd, PuBu, Greens
   - Custom: Use `get_cmap()` for brand colors

2. **Resolution**:
   - Quick preview: cmapbins=100
   - Standard: cmapbins=256
   - Smooth gradient: cmapbins=1000+

3. **Asymmetric Data**:
   - Calculate midpoint: `1 - vmax/(vmax + |vmin|)`
   - Use `shiftedColorMap()` to center at zero
   - Good for strain, stress, deformation data

### Projections
1. **Equal-area (Schmidt)**:
   - Use for: Texture analysis, density plots
   - Property: Preserves area
   - Best for: Quantitative analysis

2. **Equal-angle (Wulff)**:
   - Use for: Angular relationships
   - Property: Preserves angles
   - Best for: Crystallographic angles

3. **Triangle**:
   - Use for: Cubic systems
   - Shows: Standard stereographic triangle
   - Best for: Inverse pole figures

### Performance
1. **Memory**:
   - Use appropriate grid resolution
   - Lower nump for preview
   - Higher nump for final output

2. **Speed**:
   - Cache computed projections
   - Batch export formats
   - Use simpler colormaps for iteration

3. **Storage**:
   - PNG for complex plots (smaller file size)
   - PDF for simple plots (scalable)
   - SVG for web/interactive use

## Mathematical Background

### Stereographic Projection
Maps points from sphere to plane:
```
(X, Y) = (x/(1-z), y/(1-z))  # From north pole
```

### Equal-Area Projection (Schmidt)
Preserves area:
```
(X, Y) = (x/√(1+z), y/√(1+z))
```

### Color Normalization
Maps data values to colormap indices:
```
v_n = (v - v_min)/(v_max - v_min) × (c_max - c_min) + c_min
```

### Colormap Shifting
For asymmetric ranges:
```
m = 1 - v_max/(v_max + |v_min|)
```

Example: Data from -15 to +5
```
m = 1 - 5/(5 + 15) = 0.75
```

## Attribute Reference

### Commonly Used Attributes

**Figure Control**:
- `figsize`: (width, height) in inches
- `suptitle`: Figure title
- `tight_layout`: Boolean for tight layout

**Projection**:
- `ProjType`: 'equalarea' or 'equalangle'
- `sphere`: 'full', 'half', 'triangle'
- `R2Proj`: 3×3 rotation matrix

**Colors**:
- `cmap`: Colormap name or object
- `vm`: Value range [min, max]
- `levels`: Contour levels
- `ticks`: Colorbar tick positions

**Data**:
- `oris`: Orientation matrices (N, 3, 3)
- `dirs`: Crystal directions [[u,v,w], ...]
- `norms`: Plane normals [[h,k,l], ...]

**Save**:
- `fname`: Filename string
- `crop`: Boolean for auto-crop
- `imformats`: List of format strings
- `figparam`: Dictionary of save parameters

## Troubleshooting

### Common Issues

**Figure not displaying**:
```python
# Solution: Call plt.show()
import matplotlib.pyplot as plt
p.plotProj()
plt.show()
```

**Low resolution output**:
```python
# Solution: Increase DPI
p.figparam['dpi'] = 600
p.figsave(fname='high_res.png')
```

**Colormap appears banded**:
```python
# Solution: Increase bins
scales = p.getScales(vmcbar=[0,10], cmapbins=1000)
```

**Cropping not working**:
```python
# Solution: Install wand library
# pip install Wand
p.figsave(fname='figure.png', crop=True)
```

**Multiple formats fail**:
```python
# Solution: Use list notation
p.figsave(fname='figure.png', imformats=['png', 'pdf', 'svg'])
# Not: imformats='png,pdf,svg'
```

## Integration Examples

### With Matplotlib
```python
import matplotlib.pyplot as plt
from plotlib import plotter

fig, ax = plt.subplots(figsize=(8, 8))
p = plotter()
p.plotProj(fig=fig, ax=ax)
ax.set_title('Custom Title')
plt.savefig('integrated.png')
```

### With NumPy Arrays
```python
import numpy as np
from plotlib import plotter

# Generate data
orientations = np.random.randn(1000, 3, 3)
density = np.random.rand(1000)

# Plot
p = plotter()
p.setAttributes(oris=orientations)
p.plotProj()
```

### With Custom Colormaps
```python
from plotlib import plotter, get_cmap
import matplotlib.colors as mcolors

# Create custom colormap
colors = [mcolors.to_rgb('white'),
          mcolors.to_rgb('blue'),
          mcolors.to_rgb('red')]
cmap = get_cmap(colors, nbins=256)

# Use in plotter
p = plotter()
p.setAttributes(cmap=cmap)
```

## File Organization

The plotlib module contains:
1. **Imports** (lines 1-50): Required libraries
2. **plotter class** (~2100 lines): Main plotting class
3. **Utility functions** (~100 lines): Helper functions
4. **Documentation** (added): Comprehensive docstrings

## Related Documentation

- **plotter_class_documentation.pdf**: Complete reference manual with LaTeX formatting
- **plotter_class_documentation.tex**: LaTeX source
- **plotlib_quick_reference.md**: Quick reference guide
- **plotlib_commented.py**: Fully documented source code

## References

### Crystallographic Projections
- Kocks, U. F., Tomé, C. N., & Wenk, H. R. (2000). *Texture and anisotropy*. Cambridge University Press.
- Bunge, H. J. (2013). *Texture analysis in materials science*. Elsevier.

### Color Theory
- Moreland, K. (2016). Why we use bad color maps and what you can do about it. *Electronic Imaging*, 2016(16), 1-6.

### Scientific Visualization
- Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9(3), 90-95.

## Summary

The plotlib module provides:
- ✅ Comprehensive plotter class for crystallographic plots
- ✅ Multiple projection types (Schmidt, Wulff, triangle)
- ✅ Flexible color mapping and scaling
- ✅ High-quality multi-format export
- ✅ Utility functions for color manipulation
- ✅ Full documentation with examples
- ✅ Integration with matplotlib and numpy
- ✅ Publication-ready output

## For More Information

- **Quick start**: See plotlib_quick_reference.md
- **Detailed reference**: See plotter_class_documentation.pdf
- **Source code**: See plotlib_commented.py
- **Examples**: Check docstrings in each method

The module is production-ready for crystallographic texture analysis, pole figure generation, and orientation visualization!

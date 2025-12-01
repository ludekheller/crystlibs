# Plotlib.py Quick Reference Guide

## Overview

The `plotter` class is the main interface for creating publication-quality crystallographic plots including stereographic projections, pole figures, and orientation density maps.

## Quick Start

```python
from plotlib import plotter
import numpy as np

# Create plotter instance
p = plotter()

# Setup projection
p.plotProj(ProjType='equalarea', sphere='half')

# Plot and save
p.figsave(fname='output.png')
```

## Most Commonly Used Methods

### 1. Create Plotter Instance

```python
# Initialize plotter with defaults
p = plotter()

# All attributes can be modified using setAttributes()
p.setAttributes(cmap='viridis', figsize=(8, 8))
```

### 2. Setup Projection

```python
# Equal-area projection (Schmidt net), upper hemisphere
p.plotProj(ProjType='equalarea', sphere='half', figsize=(6, 6))

# Equal-angle projection (Wulff net), full sphere
p.plotProj(ProjType='equalangle', sphere='full')

# Stereographic triangle (for cubic systems)
p.plotProj(ProjType='equalarea', sphere='triangle', stereogrid=True)
```

### 3. Set Attributes

```python
# Set single attribute
p.setAttributes(cmap='jet')

# Set multiple attributes at once
p.setAttributes(
    dirs=[[1,0,0], [1,1,0], [1,1,1]],
    norms=[[1,0,0], [1,1,0]],
    cmap='viridis',
    figsize=(10, 10)
)
```

### 4. Generate Color Scales

```python
# Basic scale from 0 to 10
scales = p.getScales(vmcbar=[0, 10])

# Custom ticks with inequality symbols
scales = p.getScales(
    vmcbar=[0, 100],
    ticks=np.array([0, 25, 50, 75, 100]),
    geq=True,  # Add ≥ to maximum
    leq=True   # Add ≤ to minimum
)

# Access scale components
ticklabels = scales['tickslabels']
cmap = scales['cmap']
ticks = scales['ticks']
vm = scales['vm']  # [min, max] for colormap
```

### 5. Save Figures

```python
# Save single format
p.figsave(fname='pole_figure.png')

# Save multiple formats simultaneously
p.figsave(
    fname='pole_figure.png',
    imformats=['png', 'pdf', 'svg']
)
# Creates: pole_figure.png, pole_figure.pdf, pole_figure.svg

# Save with cropping
p.figsave(fname='figure.png', crop=True)

# High DPI
p.figparam['dpi'] = 600
p.figsave(fname='high_res.png')
```

## Complete Workflow Examples

### Example 1: Basic Pole Figure

```python
import numpy as np
from plotlib import plotter

# Create plotter
p = plotter()

# Setup equal-area projection, upper hemisphere
p.plotProj(ProjType='equalarea', sphere='half')

# Define crystal directions to plot
p.setAttributes(
    dirs=[[1,0,0], [0,1,0], [0,0,1]],  # <100> directions
    norms=[[1,1,1]]                      # {111} plane
)

# Plot directions and normals (if plotDirsNorms method exists)
# p.plotDirsNorms()

# Save
p.figsave(fname='basic_pole_figure.png')
```

### Example 2: Density Contour Plot

```python
import numpy as np
from plotlib import plotter

# Generate random orientation data
N = 1000
orientations = np.random.randn(N, 3, 3)

# Create plotter with orientations
p = plotter()
p.setAttributes(
    oris=orientations,
    ProjType='equalarea',
    sphere='half'
)

# Setup projection
p.plotProj()

# Plot colormap (if plotColormap method exists)
# p.plotColormap()

# Add colorbar
# p.plotColorbar()

# Save
p.figsave(fname='density_plot.png', imformats=['png', 'pdf'])
```

### Example 3: Multi-Panel Figure

```python
import matplotlib.pyplot as plt
import numpy as np
from plotlib import plotter

# Create 2x2 subplot
fig, axes = plt.subplots(2, 2, figsize=(12, 12))

# Create different pole figures in each subplot
for i, ax in enumerate(axes.flat):
    p = plotter()
    p.plotProj(
        ProjType='equalarea',
        sphere='half',
        ax=ax,
        fig=fig
    )
    
    # Different directions for each subplot
    dirs = [[i+1, 0, 0], [0, i+1, 0], [0, 0, i+1]]
    p.setAttributes(dirs=dirs)
    # p.plotDirsNorms()

plt.tight_layout()
plt.savefig('multi_panel.png', dpi=300)
```

### Example 4: Custom Colormap and Scale

```python
from plotlib import plotter, get_cmap
import numpy as np

# Create custom colormap
colors = [(1, 1, 1), (0, 0, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0)]
custom_cmap = get_cmap(colors, nbins=256)

# Create plotter
p = plotter()
p.setAttributes(cmap=custom_cmap)

# Generate color scale
scales = p.getScales(
    vmcbar=[0, 50],
    numticks=6,
    cmapbins=500
)

# Use in plotting
p.plotProj(ProjType='equalarea')
# Apply colormap with scales
```

## Attribute Categories

### Colorbar Attributes
```python
p.setAttributes(
    cbartitle="Density (MRD)",  # Colorbar title
    vmbar=[0, 10],               # Value range
    cbarh=0.04,                  # Height
    cbarwfac=0.75,               # Width factor
    cbarhshift=-0.15             # Horizontal shift
)
```

### Figure Attributes
```python
p.setAttributes(
    fig=None,                    # Figure object (created if None)
    suptitle="Pole Figure",      # Figure title
    figsize=(8, 8),              # Figure size (width, height)
    annotationtext="Sample A"    # Annotation text
)
```

### Projection Attributes
```python
p.setAttributes(
    sphere='half',               # 'full', 'half', 'triangle'
    R2Proj=np.eye(3),           # Rotation to projection frame
    stereogrid=True,             # Show grid
    stereoresolution=10          # Grid resolution
)
```

### Axis Attributes
```python
p.setAttributes(
    cmap='jet',                  # Colormap name
    vm=[0, 10],                  # Value range
    levels=np.linspace(0,10,11), # Contour levels
    ticks=[0, 2.5, 5, 7.5, 10], # Colorbar ticks
    ticklabels=['0','2.5','5','7.5','10']
)
```

### Crystal Attributes
```python
p.setAttributes(
    dirs=[[1,0,0], [1,1,0]],    # Crystal directions [uvw]
    norms=[[1,0,0], [1,1,1]],   # Plane normals (hkl)
    symops=[np.eye(3)],          # Symmetry operations
    dhkl=False                   # Show d-spacing
)
```

### Save Attributes
```python
p.setAttributes(
    fname='output.png',          # Filename
    crop=False,                  # Auto-crop whitespace
    imformats=['png','pdf'],     # Multiple formats
    figparam={'dpi': 300}        # Figure parameters
)
```

## Standalone Utility Functions

### 1. Custom Colormap Creation

```python
from plotlib import get_cmap

# Create colormap from color list
colors = [(1, 1, 1), (1, 0, 0), (0, 0, 0)]  # White→Red→Black
cmap = get_cmap(colors, nbins=256)

# Use in plots
import matplotlib.pyplot as plt
data = np.random.rand(10, 10)
plt.imshow(data, cmap=cmap)
plt.colorbar()
```

### 2. Map Values to Colors

```python
from plotlib import get_colors
import matplotlib.pyplot as plt

# Map values to colors
values = np.array([1, 5, 10, 15, 20])
cmap = plt.cm.viridis
colors = get_colors(values, cmap, vmin=0, vmax=20, to255=False)

# Use with scatter plot
x = np.random.rand(100)
y = np.random.rand(100)
z = np.random.rand(100) * 100
colors = get_colors(z, plt.cm.jet, to255=False)
plt.scatter(x, y, c=colors)
```

### 3. Shift Colormap Center

```python
from plotlib import shiftedColorMap
import matplotlib.pyplot as plt
import numpy as np

# Data ranging from -15 to +5
data = np.random.randn(100, 100) * 10 - 5
vmin, vmax = -15, 5

# Calculate midpoint to center at zero
midpoint = 1 - vmax / (vmax + abs(vmin))  # 0.75

# Create shifted colormap
original_cmap = plt.cm.RdBu_r
shifted_cmap = shiftedColorMap(original_cmap, midpoint=midpoint)

# Plot with shifted colormap (zero = white)
plt.imshow(data, cmap=shifted_cmap, vmin=vmin, vmax=vmax)
plt.colorbar()
plt.title('Shifted Colormap (zero = white)')
```

## Method Parameters Quick Reference

### plotProj()
```python
p.plotProj(
    ProjType='equalarea',     # 'equalarea' or 'equalangle'
    sphere='half',            # 'full', 'half', 'triangle'
    figsize=(6, 6),          # Figure size in inches
    stereogrid=False,         # Show grid
    stereoresolution=None,    # Grid resolution
    ax=None,                  # Existing axis (optional)
    fig=None                  # Existing figure (optional)
)
```

### getScales()
```python
scales = p.getScales(
    vmcbar=[0, 10],          # Value range [min, max]
    numticks=None,            # Number of ticks (auto if None)
    ticks=None,               # Explicit tick positions
    tickslabels=None,         # Custom tick labels
    geq=False,                # Add ≥ to max tick
    leq=False,                # Add ≤ to min tick
    cmapbins=100,             # Colormap resolution
    cmapbinsmult=None         # Bins multiplier
)
```

### getFigparam()
```python
params = p.getFigparam(
    fontsize=None,            # Font size
    save=False,               # Use save-optimized size
    phase='A',                # Phase identifier
    figsize=None,             # Custom size
    dpi=300,                  # Resolution
    format='png'              # Image format
)
```

### figsave()
```python
p.figsave(
    fname='output.png',       # Filename
    imformats=None,           # ['png','pdf','svg']
    crop=False,               # Auto-crop whitespace
    figparam=None             # Figure parameters dict
)
```

## Projection Types

### Equal-Area (Schmidt Net)
- **Use for**: Texture analysis, pole figures, density plots
- **Property**: Preserves area (equal solid angles → equal areas)
- **Best for**: Quantitative texture analysis
- **Projection**: `ProjType='equalarea'`

### Equal-Angle (Wulff Net)
- **Use for**: Angular relationships, crystallography
- **Property**: Preserves angles (conformal projection)
- **Best for**: Studying angular relationships
- **Projection**: `ProjType='equalangle'`

## Sphere Options

```python
# Full sphere projection
sphere='full'          # Complete stereographic projection

# Upper hemisphere only
sphere='half'          # Most common for pole figures

# Standard stereographic triangle
sphere='triangle'      # For cubic crystal systems
```

## Color Scale Tips

### Creating Effective Color Scales

```python
# For continuous data (texture intensity)
p.setAttributes(cmap='jet')  # or 'viridis', 'plasma'

# For diverging data (strain, stress)
p.setAttributes(cmap='RdBu_r')  # Red-Blue diverging

# For sequential data (density)
p.setAttributes(cmap='YlOrRd')  # Yellow-Orange-Red

# Custom gradient
colors = [(1,1,1), (0,0,1), (0,1,0), (1,1,0), (1,0,0)]
custom_cmap = get_cmap(colors, nbins=256)
p.setAttributes(cmap=custom_cmap)
```

### Handling Asymmetric Data

```python
# Data range: -15 to +5 (want zero = white)
vmin, vmax = -15, 5
midpoint = 1 - vmax / (vmax + abs(vmin))  # 0.75

shifted_cmap = shiftedColorMap(plt.cm.RdBu_r, midpoint=midpoint)
p.setAttributes(cmap=shifted_cmap, vm=[vmin, vmax])
```

## Common Workflows

### Workflow 1: Quick Pole Figure

```python
from plotlib import plotter

p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')
p.setAttributes(dirs=[[1,0,0], [1,1,0], [1,1,1]])
# p.plotDirsNorms()
p.figsave(fname='quick_pf.png')
```

### Workflow 2: High-Quality Publication Figure

```python
from plotlib import plotter

p = plotter()
p.setAttributes(figsize=(8, 8))
p.plotProj(ProjType='equalarea', sphere='half')

# High resolution settings
p.figparam['dpi'] = 600
p.setAttributes(
    cmap='viridis',
    vm=[0, 10],
    tight_layout=True
)

# Save in multiple formats
p.figsave(
    fname='publication_figure.png',
    imformats=['png', 'pdf', 'svg'],
    crop=True
)
```

### Workflow 3: Custom Styled Figure

```python
from plotlib import plotter, get_cmap

# Create custom colormap
colors = [(1,1,1), (0,0,1), (0,1,0), (1,1,0), (0.5,0,0)]
custom_cmap = get_cmap(colors, nbins=500)

p = plotter()
p.setAttributes(
    figsize=(10, 10),
    cmap=custom_cmap,
    suptitle='Custom Pole Figure'
)

# Generate custom scales
scales = p.getScales(
    vmcbar=[0, 20],
    ticks=np.array([0, 5, 10, 15, 20]),
    geq=True,
    cmapbins=1000
)

p.plotProj(ProjType='equalarea', sphere='half')
p.figsave(fname='custom_style.png')
```

## Common Pitfalls and Solutions

### 1. Figure Size Issues

```python
# Problem: Figure too small
# Solution: Set figsize before plotProj
p.setAttributes(figsize=(10, 10))
p.plotProj()
```

### 2. Low Resolution Output

```python
# Problem: Blurry saved figure
# Solution: Increase DPI
p.figparam['dpi'] = 600
p.figsave(fname='high_res.png')
```

### 3. Colormap Appears Banded

```python
# Problem: Discrete color bands visible
# Solution: Increase colormap bins
scales = p.getScales(vmcbar=[0,10], cmapbins=1000)
```

### 4. Can't Crop PNG

```python
# Problem: crop=True not working
# Solution: Install wand library
# pip install Wand
p.figsave(fname='figure.png', crop=True)
```

### 5. Whitespace Around Figure

```python
# Problem: Too much whitespace
# Solution: Use tight_layout
p.setAttributes(tight_layout=True)
p.figsave(fname='tight.png')
```

## Performance Tips

1. **Use appropriate DPI**:
   - Preview: 150 DPI
   - Screen: 300 DPI
   - Print: 600 DPI

2. **Colormap resolution**:
   - Quick preview: cmapbins=100
   - Standard: cmapbins=256
   - Smooth gradient: cmapbins=1000

3. **Figure formats**:
   - Raster: PNG (good for photos/complex plots)
   - Vector: PDF, SVG (scalable, good for simple plots)

4. **Batch export**:
   ```python
   # Export all formats at once
   p.figsave(fname='figure.png', imformats=['png','pdf','svg'])
   ```

## Integration with Matplotlib

```python
import matplotlib.pyplot as plt
from plotlib import plotter

# Create matplotlib figure first
fig, ax = plt.subplots(figsize=(8, 8))

# Use with plotter
p = plotter()
p.plotProj(
    ProjType='equalarea',
    sphere='half',
    fig=fig,
    ax=ax
)

# Continue with matplotlib commands
ax.set_title('My Custom Title', fontsize=14)
plt.savefig('integrated.png', dpi=300)
```

## Mathematical Formulas

### Color Scale Normalization
```
v_normalized = (v - v_min)/(v_max - v_min) × (c_max - c_min) + c_min
```

### Colormap Midpoint (for asymmetric data)
```
midpoint = 1 - v_max/(v_max + |v_min|)
```

Example: For data ranging from -15 to +5:
```
midpoint = 1 - 5/(5 + 15) = 1 - 5/20 = 0.75
```

## Best Practices

### For Publication Figures
1. Use DPI ≥ 300 (600 for print)
2. Save as vector format (PDF, SVG) when possible
3. Use perceptually uniform colormaps (viridis, plasma)
4. Enable tight_layout=True
5. Use crop=True to remove whitespace
6. Export multiple formats simultaneously

### For Colormaps
1. **Continuous data**: Use sequential colormaps (viridis, YlOrRd)
2. **Diverging data**: Use diverging colormaps (RdBu, RdYlGn)
3. **Categorical data**: Use qualitative colormaps (Set1, tab10)
4. **Asymmetric ranges**: Use shiftedColorMap()
5. **High resolution**: Use cmapbins ≥ 256

### For Projections
1. **Texture analysis**: Use equal-area (Schmidt net)
2. **Angular relationships**: Use equal-angle (Wulff net)
3. **Cubic systems**: Consider stereographic triangle
4. **Hemisphere vs full**: Use hemisphere for most pole figures

## Quick Troubleshooting

| Problem | Solution |
|---------|----------|
| Figure not showing | Call `plt.show()` or check `withdraw` parameter |
| Low resolution | Increase `figparam['dpi']` |
| Banded colormap | Increase `cmapbins` in getScales() |
| Wrong projection | Check `ProjType='equalarea'` or `'equalangle'` |
| Cropping fails | Install wand: `pip install Wand` |
| File not saved | Check `fname` path is valid |
| Colors look wrong | Verify `vm` range matches data |
| Multiple formats fail | Check `imformats` list syntax |

## For More Information

- **Full documentation**: See `plotlib_commented.py` for complete docstrings
- **Class documentation**: See `plotter_class_documentation.pdf` for detailed reference
- **Examples**: Check the Usage Examples in each function's docstring
- **Integration**: See examples for matplotlib, numpy, scipy integration

---

**Quick tip**: Start with default settings, then customize. The plotter class is designed to work with minimal configuration while allowing extensive customization when needed.

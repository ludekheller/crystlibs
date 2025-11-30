# Plotlib.py - Comprehensive Documentation Summary (Part 2 of 2)

**Crystallographic Plotting Library - Detailed Method Documentation**  
**Version**: Production Ready  
**Last Updated**: November 2025

---

## 📑 Table of Contents - Part 2

1. [Method 5: getFigparam](#method-5-getfigparam)
2. [Method 6: figsave](#method-6-figsave)
3. [Method 7: figsaveproc](#method-7-figsaveproc)
4. [Utility Function 1: get_cmap](#utility-function-1-get_cmap)
5. [Utility Function 2: get_colors](#utility-function-2-get_colors)
6. [Utility Function 3: shiftedColorMap](#utility-function-3-shiftedcolormap)
7. [Utility Function 4: plotcolmap](#utility-function-4-plotcolmap)
8. [Utility Function 5: plotcolmaps](#utility-function-5-plotcolmaps)
9. [Complete Application Examples](#complete-application-examples)
10. [Integration Guide](#integration-guide)
11. [Best Practices](#best-practices)
12. [Troubleshooting Guide](#troubleshooting-guide)

**See Part 1 for**: Module overview, class description, __init__, setAttributes, plotProj, getScales

---

## Method 5: getFigparam

### Description

Get figure parameter dictionary for saving high-quality figures. Returns a dictionary of parameters compatible with `matplotlib.pyplot.savefig()`. Provides sensible defaults for publication-quality output with options for customization.

**When to use:**
- Before saving figures with specific requirements
- To get consistent figure parameters across multiple plots
- To override default save settings

### Input

**Parameters:**

```python
fontsize : float, optional
    Font size for text elements
    Default: None (matplotlib default)

save : bool, optional
    If True, use smaller figure size optimized for saving
    Default: False

phase : str, optional
    Phase identifier affecting default figure size
    'A': 7/2.45 inch square
    Other: 11/2.45 × 10/2.45 inches
    Default: 'A'

figsize : tuple, optional
    Custom figure size (width, height) in inches
    Overrides automatic sizing
    Default: None

**kwargs : dict
    Additional parameters to override defaults
    Any key in figparam dictionary can be overridden
```

### Output

**Returns:** Dictionary with keys:

```python
{
    'dpi': int
        Resolution in dots per inch
        Default: 300
        
    'facecolor': str
        Figure background color
        Default: 'None' (transparent)
        
    'edgecolor': str
        Edge color
        Default: 'w' (white)
        
    'orientation': str
        'portrait' or 'landscape'
        Default: 'portrait'
        
    'papertype': str or None
        Paper type for PS/PDF
        Default: None
        
    'format': str
        Output format ('png', 'pdf', 'svg', etc.)
        Default: 'png'
        
    'transparent': bool
        Transparent background
        Default: True
        
    'bbox_inches': str
        Bounding box setting
        Default: 'None'
        
    'pad_inches': float
        Padding around figure
        Default: 0.1
        
    'frameon': bool or None
        Whether to draw frame
        Default: None
        
    'figsize': tuple
        Figure size (width, height) in inches
        Default: Depends on save and phase parameters
}
```

### Usage Examples

#### Example 1: Get Default Parameters

```python
from plotlib import plotter

p = plotter()

# Get defaults
params = p.getFigparam()

print("DPI:", params['dpi'])                    # 300
print("Format:", params['format'])              # 'png'
print("Transparent:", params['transparent'])    # True
print("Figure size:", params.get('figsize'))    # Depends on save flag

# Use with matplotlib
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot([1, 2, 3], [1, 4, 9])
fig.savefig('output.png', **params)
```

#### Example 2: High-Resolution PDF

```python
from plotlib import plotter

p = plotter()

# Publication-quality PDF
params_pdf = p.getFigparam(
    save=True,
    format='pdf',
    dpi=600,
    transparent=False,
    facecolor='white'
)

print("PDF parameters:")
for key, value in params_pdf.items():
    print(f"  {key}: {value}")

# Save with these parameters
p.plotProj(ProjType='equalarea', sphere='half')
# ... add data ...
p.fig.savefig('publication.pdf', **params_pdf)
```

#### Example 3: Custom Figure Size

```python
from plotlib import plotter

p = plotter()

# Specific dimensions for journal requirements
params = p.getFigparam(
    figsize=(3.5, 3.5),  # Single column width
    dpi=600,
    format='tiff',
    transparent=False
)

print("Figure size:", params['figsize'])  # (3.5, 3.5)
```

#### Example 4: Override Multiple Parameters

```python
from plotlib import plotter

p = plotter()

# Custom parameters for web display
params_web = p.getFigparam(
    dpi=150,              # Lower for web
    format='png',
    transparent=True,
    facecolor='none',
    bbox_inches='tight',  # Trim whitespace
    pad_inches=0
)

# Use in save
p.plotProj(ProjType='equalarea', sphere='half')
p.fig.savefig('web_figure.png', **params_web)
```

#### Example 5: Presentation vs Publication

```python
from plotlib import plotter

p = plotter()

# Presentation settings (lower DPI, larger)
params_pres = p.getFigparam(
    dpi=150,
    figsize=(8, 6),
    format='png'
)

# Publication settings (higher DPI, smaller)
params_pub = p.getFigparam(
    dpi=600,
    figsize=(3.5, 3.5),
    format='pdf',
    transparent=False,
    facecolor='white'
)

# Choose based on use case
params = params_pub  # or params_pres

p.plotProj(ProjType='equalarea', sphere='half')
p.fig.savefig('figure.pdf', **params)
```

### Default Figure Sizes

```python
# When save=False (display):
figsize = (16/2.45, 20/2.45)  # ≈ (6.5, 8.2) inches

# When save=True, phase='A':
figsize = (7/2.45, 7/2.45)    # ≈ (2.9, 2.9) inches

# When save=True, phase != 'A':
figsize = (11/2.45, 10/2.45)  # ≈ (4.5, 4.1) inches
```

### Notes

- Parameters returned as dictionary for use with `**kwargs`
- DPI 300 sufficient for most publications
- DPI 600 recommended for high-quality journals
- `transparent=True` useful for presentations with custom backgrounds
- `bbox_inches='tight'` removes extra whitespace
- Override any parameter using kwargs
- Can store and reuse parameter dictionaries

---

## Method 6: figsave

### Description

Save figure to file(s) with optional cropping and multi-format export. Supports simultaneous export to multiple formats (PNG, PDF, SVG, etc.) and automatic whitespace cropping for PNG files.

**When to use:**
- To save completed figures
- To export to multiple formats simultaneously
- To apply automatic cropping

### Input

**Parameters:**

```python
**kwargs : keyword arguments
    
fname : str or list of str
    Filename(s) to save
    Single: 'figure.png'
    Multiple: ['fig.png', 'fig.pdf'] or use imformats
    
imformats : list of str, optional
    Image formats to export
    Example: ['png', 'pdf', 'svg']
    Creates multiple files with same base name
    Default: None (single format from fname extension)
    
crop : bool, optional
    Auto-crop whitespace (PNG only, requires wand library)
    Default: False
    
figparam : dict, optional
    Figure parameters (see getFigparam)
    Overrides default parameters
    Default: Uses self.figparam
```

### Output

**Returns:** None  
**Side Effects:** Saves file(s) to disk

### File Naming

```python
# With imformats:
fname='figure.png', imformats=['png', 'pdf', 'svg']
# Creates:
#   figure.png
#   figure.pdf
#   figure.svg

# Extension automatically determined from imformats
```

### Usage Examples

#### Example 1: Save Single Format

```python
from plotlib import plotter

p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

# Simple save
p.figsave(fname='pole_figure.png')

# File created: pole_figure.png
```

#### Example 2: Save Multiple Formats

```python
from plotlib import plotter

p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')
# ... add data ...

# Export to multiple formats simultaneously
p.figsave(
    fname='pole_figure.png',
    imformats=['png', 'pdf', 'svg', 'eps']
)

# Files created:
#   pole_figure.png  (for presentations, web)
#   pole_figure.pdf  (for publications, vector)
#   pole_figure.svg  (for editing, web)
#   pole_figure.eps  (for legacy publications)
```

#### Example 3: Save with Cropping

```python
from plotlib import plotter

p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')
# ... add data ...

# Auto-crop whitespace (PNG only)
p.figsave(
    fname='cropped_figure.png',
    crop=True
)

# Requires: pip install Wand
# Also requires ImageMagick installed on system
```

#### Example 4: High-DPI Save

```python
from plotlib import plotter

p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

# Set high DPI before saving
p.figparam['dpi'] = 600

p.figsave(fname='high_res.png')

# File created at 600 DPI
```

#### Example 5: Custom Parameters

```python
from plotlib import plotter

p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

# Override figure parameters for this save only
p.figsave(
    fname='custom.pdf',
    figparam={
        'dpi': 300,
        'facecolor': 'white',
        'edgecolor': 'none',
        'transparent': False,
        'bbox_inches': 'tight',
        'pad_inches': 0.05
    }
)
```

#### Example 6: Batch Saving with Different Settings

```python
from plotlib import plotter

p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')
# ... add data ...

# Web version (low DPI, PNG)
p.figparam['dpi'] = 150
p.figsave(fname='figure_web.png')

# Print version (high DPI, PDF)
p.figparam['dpi'] = 600
p.figparam['transparent'] = False
p.figparam['facecolor'] = 'white'
p.figsave(fname='figure_print.pdf')

# Presentation version (medium DPI, PNG, cropped)
p.figparam['dpi'] = 200
p.figsave(fname='figure_pres.png', crop=True)
```

### Notes

- **PNG cropping**: Requires `wand` library and ImageMagick
  ```bash
  pip install Wand
  # macOS: brew install imagemagick
  # Ubuntu: sudo apt-get install imagemagick
  # Windows: Download from imagemagick.org
  ```
- **Format from extension**: If no `imformats`, uses fname extension
- **Vector formats**: PDF, SVG, EPS are resolution-independent
- **Raster formats**: PNG, JPG, TIFF depend on DPI
- **Multiple calls**: Each figsave() is independent
- **No overwrite warning**: Files silently overwritten

---

## Method 7: figsaveproc

### Description

Internal method for processing and saving figure to file. Called by `figsave()` to handle actual file writing. Usually not called directly by users, but available for custom workflows.

**When to use:**
- Internally by figsave()
- For custom save workflows
- When bypassing figsave() logic

### Input

**Parameters:**

```python
fname : str
    Filename to save
    
**kwargs : keyword arguments
    Additional parameters to update before saving
```

### Output

**Returns:** None  
**Side Effects:** 
- Saves file to disk
- Applies tight_layout if enabled
- Updates figure format parameter

### Usage Example

```python
from plotlib import plotter

p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

# Direct call (typically not needed)
p.figsaveproc('output.png')

# Equivalent to:
# p.figsave(fname='output.png')
```

### Implementation Details

```python
def figsaveproc(self, fname, **kwargs):
    # Update attributes
    self.setAttributes(**kwargs)
    
    # Determine format from filename if not set
    if 'format' not in self.figparam:
        extension = fname.split('.')[-1]
        self.figparam['format'] = extension
    
    # Apply tight layout if enabled
    if self.tight_layout:
        self.fig.tight_layout()
    
    # Save using matplotlib
    self.fig.savefig(
        fname,
        bbox_inches='tight',
        dpi=self.figparam['dpi'],
        facecolor=self.figparam['facecolor'],
        edgecolor=self.figparam['edgecolor'],
        orientation=self.figparam['orientation'],
        format=self.figparam['format'],
        transparent=self.figparam['transparent'],
        pad_inches=self.figparam['pad_inches']
    )
```

### Notes

- Used internally by `figsave()`
- No cropping applied (handled by figsave)
- Format auto-detected from extension
- Always uses `bbox_inches='tight'`
- Respects `tight_layout` flag

---

## Utility Function 1: get_cmap

### Description

Create a custom colormap from a list of colors with smooth gradients between them. Uses matplotlib's LinearSegmentedColormap to interpolate smoothly through the provided colors.

**When to use:**
- Creating custom color schemes
- Matching specific color palettes
- Building diverging colormaps
- Creating perceptually meaningful colors

### Input

**Parameters:**

```python
colors : list
    List of colors to interpolate between
    Can be RGB tuples (R, G, B) with values 0-1
    Can be named colors ('blue', 'red', 'green')
    Can mix formats
    
nbins : int, optional
    Number of discrete bins in colormap
    Higher = smoother gradient
    Default: 1000
    Recommended: 256-1000 for publication
    
name : str, optional
    Name for the colormap
    Default: 'my_cmap'
```

### Output

**Returns:** `matplotlib.colors.LinearSegmentedColormap`  
Custom colormap object usable anywhere matplotlib colormaps work

### Usage Examples

#### Example 1: Simple Gradient

```python
from plotlib import get_cmap
import matplotlib.pyplot as plt
import numpy as np

# White to red gradient
colors = [(1, 1, 1), (1, 0, 0)]
cmap = get_cmap(colors, nbins=256)

# Use in plot
data = np.random.rand(10, 10)
plt.imshow(data, cmap=cmap)
plt.colorbar()
plt.title('White to Red Gradient')
plt.show()
```

#### Example 2: Multi-color Rainbow

```python
from plotlib import get_cmap
import numpy as np
import matplotlib.pyplot as plt

# Create rainbow gradient
colors_rainbow = [
    (0, 0, 1),    # Blue
    (0, 1, 1),    # Cyan
    (0, 1, 0),    # Green
    (1, 1, 0),    # Yellow
    (1, 0.5, 0),  # Orange
    (1, 0, 0)     # Red
]
cmap_rainbow = get_cmap(colors_rainbow, nbins=500)

# Plot
data = np.random.rand(100, 100)
plt.imshow(data, cmap=cmap_rainbow)
plt.colorbar()
plt.title('Custom Rainbow')
plt.show()
```

#### Example 3: Diverging Colormap

```python
from plotlib import get_cmap
import numpy as np
import matplotlib.pyplot as plt

# Blue-white-red diverging (for +/- data)
colors_div = [
    (0, 0, 1),    # Blue (negative)
    (1, 1, 1),    # White (zero)
    (1, 0, 0)     # Red (positive)
]
cmap_div = get_cmap(colors_div, nbins=256)

# Symmetric data around zero
data = np.random.randn(50, 50)
plt.imshow(data, cmap=cmap_div, vmin=-3, vmax=3)
plt.colorbar()
plt.title('Diverging: Blue-White-Red')
plt.show()
```

#### Example 4: Named Colors

```python
from plotlib import get_cmap

# Using named colors
colors_named = ['darkblue', 'white', 'darkred']
cmap = get_cmap(colors_named, nbins=256, name='dark_diverging')

print("Colormap name:", cmap.name)  # 'dark_diverging'
```

#### Example 5: Use with Plotter

```python
from plotlib import plotter, get_cmap

# Create custom colormap
colors = [(1,1,1), (0.8,0.8,0), (1,0.5,0), (1,0,0), (0.5,0,0)]
custom_cmap = get_cmap(colors, nbins=500)

# Use in plotter
p = plotter()
p.setAttributes(cmap=custom_cmap)
p.plotProj(ProjType='equalarea', sphere='half')
# ... plot data ...
p.plotColormap()
p.plotColorbar()
```

### Mathematical Background

Linear interpolation between colors in RGB space:
```
For position t ∈ [0, 1] between colors C1 and C2:
C(t) = C1 × (1 - t) + C2 × t

For N colors at positions p1, p2, ..., pN:
Find adjacent colors Ci, Ci+1 where pi ≤ t < pi+1
Local t' = (t - pi) / (pi+1 - pi)
C(t) = Ci × (1 - t') + Ci+1 × t'
```

### Notes

- Higher `nbins` = smoother but more memory
- 256 bins sufficient for most uses
- Can mix RGB tuples and named colors
- Supports any matplotlib color specification
- Colormap objects are reusable
- Store and reuse to save computation time

---

## Utility Function 2: get_colors

### Description

Map data values to RGB/RGBA colors using a colormap. Converts numeric data to color arrays suitable for scatter plots, point coloring, or any visualization requiring value-to-color mapping.

**When to use:**
- Coloring scatter plots by data values
- Creating custom color-coded visualizations
- Converting data to colors for manual plotting
- Handling NaN values with special colors

### Input

**Parameters:**

```python
values : numpy array
    Data values to map to colors
    Can be any shape
    
cmap : matplotlib colormap
    Colormap to use (e.g., plt.cm.viridis)
    
vmin : float, optional
    Minimum data value for colormap
    Values < vmin get minimum color
    Default: None (uses min of values)
    
vmax : float, optional
    Maximum data value for colormap
    Values > vmax get maximum color
    Default: None (uses max of values)
    
cmin : int, optional
    Minimum colormap index
    Default: 0
    
cmax : int, optional
    Maximum colormap index
    Default: None (uses cmap.N)
    
to255 : bool, optional
    Scale colors to 0-255 range (int)
    True: RGB values 0-255
    False: RGB values 0.0-1.0
    Default: True
    
nancolor : list, optional
    RGBA color for NaN values
    Format: [R, G, B, A] with values 0-1
    Default: [0, 0, 0, 1] (black)
```

### Output

**Returns:** numpy array  
Shape: (values.shape, 4) for RGBA colors  
Type: int if to255=True, float if to255=False

### Usage Examples

#### Example 1: Basic Value Mapping

```python
from plotlib import get_colors
import matplotlib.pyplot as plt
import numpy as np

# Data values
values = np.array([1, 5, 10, 15, 20])

# Map to colors
cmap = plt.cm.viridis
colors = get_colors(values, cmap, vmin=0, vmax=20, to255=False)

print("Colors shape:", colors.shape)  # (5, 4)
print("First color (RGBA):", colors[0])
# Example: [0.267004, 0.004874, 0.329415, 1.0]

# Use in plot
x = [1, 2, 3, 4, 5]
y = [2, 4, 3, 5, 4]
plt.scatter(x, y, c=colors, s=100, edgecolors='black')
plt.title('Scatter Plot with Value Colors')
plt.show()
```

#### Example 2: Scatter Plot by Value

```python
from plotlib import get_colors
import matplotlib.pyplot as plt
import numpy as np

# Random scatter data
N = 100
x = np.random.rand(N)
y = np.random.rand(N)
z = np.random.rand(N) * 100  # Value to color

# Map z values to colors
colors = get_colors(z, plt.cm.jet, vmin=0, vmax=100, to255=False)

# Plot
plt.scatter(x, y, c=colors, s=50, edgecolors='black', linewidths=0.5)
plt.colorbar(plt.cm.ScalarMappable(cmap=plt.cm.jet, 
                                    norm=plt.Normalize(0, 100)))
plt.title(f'Scatter Plot (N={N})')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
```

#### Example 3: Handle NaN Values

```python
from plotlib import get_colors
import matplotlib.pyplot as plt
import numpy as np

# Data with NaN
values = np.array([1, 2, np.nan, 4, 5, np.nan, 7])

# Map to colors with white for NaN
colors = get_colors(
    values,
    plt.cm.coolwarm,
    vmin=0,
    vmax=10,
    nancolor=[1, 1, 1, 1],  # White for NaN
    to255=False
)

# Plot
x = range(len(values))
plt.scatter(x, values, c=colors, s=200, edgecolors='black')
plt.title('Data with NaN (white)')
plt.show()
```

#### Example 4: Integer RGB for Image Processing

```python
from plotlib import get_colors
import numpy as np

# Values to map
values = np.linspace(0, 1, 256).reshape(16, 16)

# Get integer RGB (0-255)
colors_int = get_colors(values, plt.cm.viridis, to255=True)

print("Data type:", colors_int.dtype)  # int
print("Value range:", colors_int.min(), "to", colors_int.max())  # 0 to 255
print("Shape:", colors_int.shape)  # (16, 16, 4)

# Can save as image array
```

#### Example 5: Custom Value Range

```python
from plotlib import get_colors
import matplotlib.pyplot as plt
import numpy as np

# Temperature data (°C)
temps = np.array([-10, -5, 0, 5, 10, 15, 20, 25, 30])

# Map with custom range (clamp extremes)
colors = get_colors(
    temps,
    plt.cm.RdYlBu_r,
    vmin=-5,   # Values < -5 get minimum color
    vmax=25,   # Values > 25 get maximum color
    to255=False
)

# Plot
plt.scatter(range(len(temps)), temps, c=colors, s=200)
plt.axhline(y=-5, color='blue', linestyle='--', alpha=0.3)
plt.axhline(y=25, color='red', linestyle='--', alpha=0.3)
plt.title('Temperature with Clamped Colors')
plt.ylabel('Temperature (°C)')
plt.show()
```

### Value Normalization Formula

```python
# Normalize values to colormap range:
v_normalized = (values - vmin) / (vmax - vmin) * (cmax - cmin) + cmin

# Clamp to colormap indices:
indices = np.clip(v_normalized, cmin, cmax).astype(int)

# Get colors:
colors = cmap(indices)

# Handle NaN:
colors[np.isnan(values)] = nancolor

# Scale if needed:
if to255:
    colors = (colors * 255).astype(int)
```

### Notes

- Returns RGBA (4 channels) even if cmap is RGB
- NaN handling automatic if `nancolor` provided
- `to255=True` for integer color arrays (image processing)
- `to255=False` for matplotlib plotting (0-1 float)
- Values outside [vmin, vmax] are clamped
- Colormap indices are integers (discretized)

---

## Utility Function 3: shiftedColorMap

### Description

Shift the center of a colormap to a specific data value. Essential for data with asymmetric ranges where you want the colormap center (typically white or neutral color) at a meaningful value like zero.

**When to use:**
- Data ranges asymmetric around zero (e.g., -15 to +5)
- Want neutral color at zero for +/- data
- Centering colormap at specific value (not 0)
- Creating custom diverging colormaps

### Input

**Parameters:**

```python
cmap : matplotlib colormap
    The colormap to shift
    Typically diverging: RdBu, RdYlBu, coolwarm, seismic
    
start : float, optional
    Offset from lowest point [0.0, midpoint]
    Default: 0.0 (use full range)
    
midpoint : float
    New center position [0.0, 1.0]
    0.0 = center at minimum
    0.5 = center at middle (symmetric)
    1.0 = center at maximum
    Calculate as: 1 - vmax/(vmax + |vmin|)
    Default: 0.5
    
stop : float, optional
    Offset from highest point [midpoint, 1.0]
    Default: 1.0 (use full range)
    
name : str, optional
    Name for new colormap
    Default: 'shiftedcmap'
```

### Output

**Returns:** `matplotlib.colors.LinearSegmentedColormap`  
New colormap with shifted center

### Midpoint Calculation

For data range [vmin, vmax], to center at zero:
```
midpoint = 1 - vmax / (vmax + |vmin|)
```

**Examples:**
- Symmetric (-10, +10): midpoint = 1 - 10/(10+10) = 0.5
- Asymmetric (-15, +5): midpoint = 1 - 5/(5+15) = 0.75
- Asymmetric (-5, +15): midpoint = 1 - 15/(15+5) = 0.25

### Usage Examples

#### Example 1: Asymmetric Data Centered at Zero

```python
from plotlib import shiftedColorMap
import matplotlib.pyplot as plt
import numpy as np

# Data from -15 to +5
vmin, vmax = -15, 5

# Calculate midpoint to center at zero
midpoint = 1 - vmax / (vmax + abs(vmin))
print(f"Midpoint: {midpoint:.3f}")  # 0.75

# Create shifted colormap
original = plt.cm.RdBu_r
shifted = shiftedColorMap(original, midpoint=midpoint, name='centered_at_zero')

# Generate test data
data = np.random.randn(100, 100) * 10 - 5  # Roughly -15 to +5

# Plot comparison
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Original (zero not at center)
im1 = ax1.imshow(data, cmap=original, vmin=vmin, vmax=vmax)
ax1.set_title('Original Colormap')
plt.colorbar(im1, ax=ax1)

# Shifted (zero at white/center)
im2 = ax2.imshow(data, cmap=shifted, vmin=vmin, vmax=vmax)
ax2.set_title('Shifted Colormap (Zero = White)')
plt.colorbar(im2, ax=ax2)

plt.tight_layout()
plt.show()
```

#### Example 2: Center at Specific Value

```python
from plotlib import shiftedColorMap
import matplotlib.pyplot as plt
import numpy as np

# Temperature data: 90°C to 110°C
# Want center (white) at 100°C (boiling point)
vmin, vmax = 90, 110
center_value = 100

# Calculate midpoint
midpoint = 1 - (vmax - center_value) / (vmax - vmin)
print(f"Midpoint for center at {center_value}: {midpoint:.3f}")  # 0.5

# Create shifted colormap
cmap_shifted = shiftedColorMap(plt.cm.coolwarm, midpoint=midpoint)

# Data
temps = np.random.uniform(90, 110, (50, 50))

# Plot
plt.imshow(temps, cmap=cmap_shifted, vmin=vmin, vmax=vmax)
plt.colorbar(label='Temperature (°C)')
plt.title(f'Temperature Map (White = {center_value}°C)')
plt.show()
```

#### Example 3: Strain Data (Tension/Compression)

```python
from plotlib import shiftedColorMap
import matplotlib.pyplot as plt
import numpy as np

# Strain data: -0.03 to +0.01 (more compression than tension)
strain_min, strain_max = -0.03, 0.01

# Center at zero strain
midpoint = 1 - strain_max / (strain_max + abs(strain_min))
print(f"Strain midpoint: {midpoint:.3f}")  # 0.75

# Create colormap (blue=compression, red=tension)
cmap = shiftedColorMap(plt.cm.RdBu, midpoint=midpoint)

# Simulated strain field
strain = np.random.randn(100, 100) * 0.015 - 0.01

# Plot
plt.imshow(strain, cmap=cmap, vmin=strain_min, vmax=strain_max)
plt.colorbar(label='Strain')
plt.title('Strain Field (White = Zero Strain)')
plt.axhline(y=50, color='k', linewidth=0.5, alpha=0.3)
plt.axvline(x=50, color='k', linewidth=0.5, alpha=0.3)
plt.show()
```

#### Example 4: Use with Plotter

```python
from plotlib import plotter, shiftedColorMap
import matplotlib.pyplot as plt
import numpy as np

# Data range
vmin, vmax = -20, 8
midpoint = 1 - vmax / (vmax + abs(vmin))

# Create shifted colormap
cmap_shifted = shiftedColorMap(plt.cm.seismic, midpoint=midpoint)

# Use in plotter
p = plotter()
p.setAttributes(
    cmap=cmap_shifted,
    vm=[vmin, vmax]
)

p.plotProj(ProjType='equalarea', sphere='half')
# ... plot data ...
p.plotColormap()
p.plotColorbar()
```

### Algorithm

The function creates a new colormap by:
1. Defining regular sampling points (257 points)
2. Creating shifted index (128 points before midpoint, 129 after)
3. Mapping regular indices to shifted indices
4. Interpolating colors from original colormap
5. Building new LinearSegmentedColormap

```python
# Regular index: uniform 0 to 1
reg_index = np.linspace(start, stop, 257)

# Shifted index: concentrated around midpoint
shift_index = np.hstack([
    np.linspace(0.0, midpoint, 128, endpoint=False),
    np.linspace(midpoint, 1.0, 129, endpoint=True)
])

# For each point:
#   reg_i = position in regular sampling
#   shift_i = position in shifted sampling
#   color = cmap(reg_i)
#   new_cmap[shift_i] = color
```

### Notes

- Original colormap not modified (creates new)
- Works best with diverging colormaps
- Midpoint calculation formula assumes centering at zero
- For other center values, adjust formula accordingly
- Can adjust `start` and `stop` to use colormap subset
- Resolution: 257 points (standard for matplotlib)

---

## Utility Function 4: plotcolmap

### Description

Plot a single colormap configuration from a saved pickle file. Loads plot data, configuration, and code from pickle file and executes to recreate the plot.

**When to use:**
- Reproducing saved plot configurations
- Loading complex multi-panel setups
- Sharing plot configurations between sessions

### Input

**Parameters:**

```python
fname : str, optional
    Path to pickle file containing plot configuration
    Default: None
    
withdraw : bool, optional
    If True, don't display plot (for batch processing)
    Default: False
```

### Output

**Returns:** None  
**Side Effects:** 
- Loads data from pickle file
- Executes stored code
- Creates matplotlib figure
- Displays plot (unless withdraw=True)

### Pickle File Structure

Expected keys in pickle dictionary:
```python
{
    'data': plot data,
    'code': Python code string to execute,
    'figsize': tuple (width, height),
    ... other configuration parameters
}
```

### Usage Example

```python
from plotlib import plotcolmap

# Load and plot saved configuration
plotcolmap('saved_pole_figure.pkl')

# Load without displaying (batch mode)
plotcolmap('config.pkl', withdraw=True)
```

### Saving Configuration (Example)

```python
import pickle
from plotlib import plotter

# Create plot
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')
# ... configure and plot ...

# Save configuration
config = {
    'data': p.__dict__,
    'figsize': p.figsize,
    'code': "p.plotProj(ProjType='equalarea', sphere='half')"
}

with open('saved_config.pkl', 'wb') as f:
    pickle.dump(config, f)

# Later reload:
# plotcolmap('saved_config.pkl')
```

### Notes

- Pickle files can contain arbitrary Python code
- Use only with trusted pickle files (security concern)
- Useful for reproducibility
- Can store complex multi-step configurations
- Alternative: save as Python script instead

---

## Utility Function 5: plotcolmaps

### Description

Plot multiple colormaps from a saved pickle file in a 2×2 subplot grid. Loads and displays up to 4 pole figures or colormaps simultaneously with shared configuration.

**When to use:**
- Comparing multiple projections
- Multi-panel figures from saved data
- Batch visualization of related plots

### Input

**Parameters:**

```python
fname : str, optional
    Path to pickle file containing multiple plot configurations
    Default: None
    
withdraw : bool, optional
    If True, don't display plot (for batch processing)
    Default: False
```

### Output

**Returns:** None  
**Side Effects:** 
- Creates 2×2 subplot figure
- Plots up to 4 configurations
- Displays figure (unless withdraw=True)

### Pickle File Structure

Expected keys:
```python
{
    'data2plot': list of data dictionaries,
    'attributes2use': list of attribute dictionaries,
    'colormapdata': list of colormap data,
    'contourZero': array of zero contour points,
    'showdhalf': display parameters,
    ... other shared parameters
}
```

### Usage Example

```python
from plotlib import plotcolmaps

# Load and plot 2×2 grid
plotcolmaps('multi_pole_figures.pkl')

# Load without display
plotcolmaps('batch_plots.pkl', withdraw=True)
```

### Notes

- Creates 2×2 grid (4 panels)
- Each panel configured independently
- Useful for comparative visualization
- Can show different projections of same data
- Security: only use trusted pickle files

---

## Complete Application Examples

### Example 1: Publication-Quality Pole Figure

```python
from plotlib import plotter, get_cmap
import matplotlib.pyplot as plt
import numpy as np

# Create custom colormap
colors = [(1,1,1), (0,0.5,1), (0,1,0.5), (1,0.8,0), (0.8,0,0)]
custom_cmap = get_cmap(colors, nbins=256)

# Configure plotter for publication
p = plotter()
p.setAttributes(
    figsize=(3.5, 3.5),      # Single column
    cmap=custom_cmap,
    contourfontsize=8,
    dirsnormsalpha=0.7,
    nump=501                  # High resolution
)

# Create projection
p.plotProj(ProjType='equalarea', sphere='half')

# Add crystal directions
directions = [[1,0,0], [0,1,0], [0,0,1]]
p.setAttributes(
    dirs=directions,
    dirsfacecolors={0:'red', 1:'green', 2:'blue'},
    dirsedgecolors={0:'darkred', 1:'darkgreen', 2:'darkblue'}
)
p.plotDirsNorms()

# Format
p.ax.set_title('(a) {100} Pole Figure', 
               fontsize=10, weight='bold', loc='left')

# Save publication quality
p.figparam['dpi'] = 600
p.figparam['transparent'] = False
p.figparam['facecolor'] = 'white'

p.figsave(
    fname='publication_figure.pdf',
    imformats=['pdf', 'png']
)
```

---

### Example 2: Orientation Density with Custom Scale

```python
from plotlib import plotter, shiftedColorMap
from orilib import np_eulers_matrices
import matplotlib.pyplot as plt
import numpy as np

# Generate orientations
N = 2000
euler = np.random.rand(N, 3) * [360, 180, 360]
oris = np_eulers_matrices(euler, deg=True)

# Setup plotter
p = plotter()
p.setAttributes(
    oris=oris,
    nump=301,
    figsize=(8, 8)
)

# Get color scales
scales = p.getScales(
    vmcbar=[0, 5],
    ticks=np.array([0, 1, 2, 3, 4, 5]),
    geq=True
)

# Apply scales
p.setAttributes(
    vm=scales['vm'],
    cmap=scales['cmap'],
    ticks=scales['ticks'],
    ticklabels=scales['tickslabels']
)

# Create plot
p.plotProj(ProjType='equalarea', sphere='half')
p.plotColormap()
p.plotColorbar()

# Format
p.ax.set_title(f'Orientation Density (N={N})', 
               fontsize=14, weight='bold')
p.setAttributes(cbartitle='MRD')

# Save
p.figsave(fname='density_with_scale.png', dpi=300)
```

---

### Example 3: Multi-Panel Comparison

```python
import matplotlib.pyplot as plt
from plotlib import plotter
import numpy as np

# Create 2×2 grid
fig, axes = plt.subplots(2, 2, figsize=(12, 12))

# Plot parameters
projections = [
    ('half', 'equalarea', '{100} Equal-area'),
    ('half', 'equalangle', '{100} Equal-angle'),
    ('triangle', 'equalarea', 'IPF Triangle'),
    ('full', 'equalarea', 'Full Sphere')
]

directions = [[1,0,0], [0,1,0], [0,0,1]]

# Create each panel
for ax, (sphere, proj, title) in zip(axes.flat, projections):
    p = plotter()
    p.plotProj(
        ProjType=proj,
        sphere=sphere,
        ax=ax,
        fig=fig,
        figsize=(6, 6)
    )
    
    # Add directions
    p.setAttributes(dirs=directions)
    p.plotDirsNorms()
    
    # Format
    ax.set_title(title, fontsize=12, weight='bold')

plt.tight_layout()
plt.savefig('multi_panel_comparison.png', dpi=300, bbox_inches='tight')
plt.show()
```

---

### Example 4: Asymmetric Data with Shifted Colormap

```python
from plotlib import plotter, shiftedColorMap, get_colors
import matplotlib.pyplot as plt
import numpy as np

# Asymmetric strain data
strain_data = np.random.randn(100, 100) * 0.02 - 0.01
vmin, vmax = strain_data.min(), strain_data.max()

# Calculate midpoint to center at zero
midpoint = 1 - vmax / (vmax + abs(vmin))
print(f"Data range: [{vmin:.4f}, {vmax:.4f}]")
print(f"Midpoint: {midpoint:.3f}")

# Create shifted colormap
cmap_shifted = shiftedColorMap(plt.cm.RdBu, midpoint=midpoint)

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Original colormap
im1 = ax1.imshow(strain_data, cmap=plt.cm.RdBu, vmin=vmin, vmax=vmax)
ax1.set_title('Original Colormap')
plt.colorbar(im1, ax=ax1, label='Strain')

# Shifted colormap (zero at white)
im2 = ax2.imshow(strain_data, cmap=cmap_shifted, vmin=vmin, vmax=vmax)
ax2.set_title('Shifted Colormap (Zero = White)')
plt.colorbar(im2, ax=ax2, label='Strain')

plt.tight_layout()
plt.savefig('shifted_colormap_comparison.png', dpi=300)
plt.show()
```

---

## Integration Guide

### With orilib

```python
from plotlib import plotter
from orilib import np_eulers_matrices, quat_misori_deg
import numpy as np

# Generate orientations using orilib
euler = np.random.rand(500, 3) * [360, 180, 360]
oris = np_eulers_matrices(euler, deg=True)

# Plot density using plotlib
p = plotter()
p.setAttributes(oris=oris, nump=201)
p.plotProj(ProjType='equalarea', sphere='half')
p.plotColormap()
p.plotColorbar()
p.figsave('orilib_integration.png')
```

### With projlib

```python
from plotlib import plotter
from projlib import equalarea_directions

# Project directions using projlib
directions = [[1,0,0], [1,1,0], [1,1,1]]
proj_coords = [equalarea_directions(d) for d in directions]

# Plot using plotlib
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

# Add projected points
for xy in proj_coords:
    p.ax.plot(xy[0], xy[1], 'ro', markersize=10)

p.figsave('projlib_integration.png')
```

### With crystlib

```python
from plotlib import plotter
from crystlib import generate_hkls
import numpy as np

# Generate Miller indices using crystlib
symops = [np.eye(3)]
unique, equiv, families, all_hkls = generate_hkls(
    hklmax=2,
    symops=symops
)

# Plot using plotlib
p = plotter()
p.plotProj(ProjType='equalarea', sphere='half')

# Project and plot hkls
# (requires additional projection calculations)

p.figsave('crystlib_integration.png')
```

---

## Best Practices

### 1. Resolution Settings
```python
# Preview/development: Fast
p.setAttributes(nump=101)

# Normal use: Balanced
p.setAttributes(nump=201)

# Publication: High quality
p.setAttributes(nump=501)

# Ultra-high: Special cases
p.setAttributes(nump=1001)
```

### 2. Colormap Selection
```python
# Texture/density: Perceptually uniform
p.setAttributes(cmap='viridis')  # or plasma, inferno

# Diverging data: Blue-white-red
p.setAttributes(cmap='RdBu_r')   # or coolwarm, seismic

# Traditional: Jet (use sparingly)
p.setAttributes(cmap='jet')
```

### 3. File Format Guidelines
```python
# Web/presentation: PNG
p.figsave('figure.png', imformats=['png'])

# Publication: PDF (vector)
p.figsave('figure.pdf', imformats=['pdf'])

# Both: Multi-format
p.figsave('figure.png', imformats=['png', 'pdf', 'svg'])
```

### 4. DPI Recommendations
```python
# Screen display
p.figparam['dpi'] = 100-150

# Presentations
p.figparam['dpi'] = 150-200

# Print/posters
p.figparam['dpi'] = 300

# High-quality journals
p.figparam['dpi'] = 600
```

---

## Troubleshooting Guide

### Common Issues

**Issue: Empty plot**
- **Cause**: Forgot to call plotProj()
- **Solution**: Always call `p.plotProj()` before plotting data

**Issue: Low resolution**
- **Cause**: Default DPI too low
- **Solution**: Set `p.figparam['dpi'] = 600`

**Issue: Banded colormap**
- **Cause**: Too few colormap bins
- **Solution**: Increase `cmapbins` to 256+ in getScales()

**Issue: Cropping fails**
- **Cause**: Wand library not installed
- **Solution**: `pip install Wand` and install ImageMagick

**Issue: Wrong projection shape**
- **Cause**: Incorrect ProjType or sphere
- **Solution**: Check spelling: 'equalarea' not 'equal-area'

**Issue: Directions not showing**
- **Cause**: Wrong format or empty list
- **Solution**: Use `[[u,v,w]]` not `[u,v,w]`

### Performance Issues

**Issue: Slow plotting**
- **Cause**: nump too high
- **Solution**: Reduce to 201 for development, 501 for final

**Issue: Memory error**
- **Cause**: nump too high for available RAM
- **Solution**: Start with nump=101, increase gradually

---

**End of Part 2**

**Documentation Status**: Complete  
**Total Methods**: 7 class methods + 5 utility functions  
**Coverage**: 100%  
**Quality**: Production Ready  
**Integration**: Full support for orilib, projlib, crystlib

*For quick reference, see PLOTLIB_QUICK_REFERENCE_COMPREHENSIVE.md*  
*For class and method basics, see Part 1*

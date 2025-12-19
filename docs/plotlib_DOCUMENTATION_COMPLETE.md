# PLOTLIB - Complete Documentation

**Module**: `plotlib.py`  
**Purpose**: Crystallographic Visualization  
**Total Functions**: 57  
**Last Updated**: December 19, 2025

---

## Table of Contents

1. [dataAnnot](#function-dataannot)
2. [dataShow](#function-datashow)
3. [figsave](#function-figsave)
4. [figsaveproc](#function-figsaveproc)
5. [format_annot](#function-format_annot)
6. [format_coord](#function-format_coord)
7. [format_coord_test](#function-format_coord_test)
8. [genPoris](#function-genporis)
9. [generateSphericalHistSampleData](#function-generatesphericalhistsampledata)
10. [generateSphericalKDESampleData](#function-generatesphericalkdesampledata)
11. [getColormap](#function-getcolormap)
12. [getFigparam](#function-getfigparam)
13. [getScales](#function-getscales)
14. [get_cmap](#function-get_cmap)
15. [get_colors](#function-get_colors)
16. [onclicactivate](#function-onclicactivate)
17. [onclick](#function-onclick)
18. [onclick2](#function-onclick2)
19. [onclick3](#function-onclick3)
20. [onmove](#function-onmove)
21. [onpress](#function-onpress)
22. [onpressActivate](#function-onpressactivate)
23. [plotColorbar](#function-plotcolorbar)
24. [plotColormap](#function-plotcolormap)
25. [plotColormaps](#function-plotcolormaps)
26. [plotDirsNorms](#function-plotdirsnorms)
27. [plotHist](#function-plothist)
28. [plotProj](#function-plotproj)
29. [plotScatter](#function-plotscatter)
30. [plotScatterAsHist](#function-plotscatterashist)
31. [plot_atomic_plane2D](#function-plot_atomic_plane2d)
32. [plot_atomic_plane3D](#function-plot_atomic_plane3d)
33. [plot_atomlattice2D](#function-plot_atomlattice2d)
34. [plot_cut2D](#function-plot_cut2d)
35. [plot_lattice](#function-plot_lattice)
36. [plot_lattice2D](#function-plot_lattice2d)
37. [plot_lattice3D](#function-plot_lattice3d)
38. [plot_lattice_2Dprojection](#function-plot_lattice_2dprojection)
39. [plot_lattice_boundaries](#function-plot_lattice_boundaries)
40. [plot_lattice_plane](#function-plot_lattice_plane)
41. [plot_lattice_proj](#function-plot_lattice_proj)
42. [plot_latticefaces3D](#function-plot_latticefaces3d)
43. [plot_latticesfaces3D](#function-plot_latticesfaces3d)
44. [plot_mohr_circles](#function-plot_mohr_circles)
45. [plot_planes_on_mohr_circle](#function-plot_planes_on_mohr_circle)
46. [plot_planes_on_stereotriangle](#function-plot_planes_on_stereotriangle)
47. [plot_planes_on_wulffnet](#function-plot_planes_on_wulffnet)
48. [plot_points_proj](#function-plot_points_proj)
49. [plot_princip_dir_on_stereotriangle](#function-plot_princip_dir_on_stereotriangle)
50. [plot_princip_dir_on_wulffnet](#function-plot_princip_dir_on_wulffnet)
51. [plotcolmap](#function-plotcolmap)
52. [plotcolmaps](#function-plotcolmaps)
53. [processScatterData](#function-processscatterdata)
54. [scatterDataAnnot](#function-scatterdataannot)
55. [setAttributes](#function-setattributes)
56. [set_aspect_equal_3d](#function-set_aspect_equal_3d)
57. [shiftedColorMap](#function-shiftedcolormap)

---

## Function: dataAnnot

**Signature**:
```python
def dataAnnot(self,**kwargs):
```

**Description**:

Annotate data points with crystallographic information.

---

## Function: dataShow

**Signature**:
```python
def dataShow(self,**kwargs):
```

**Description**:

Display data information in plot or console.

---

## Function: figsave

**Signature**:
```python
def figsave(self, **kwargs):
```

**Description**:

Save figure to file(s) with optional cropping.
        Supports multiple image formats and automatic cropping of whitespace.
        Can save to multiple formats simultaneously.

**Input**:

**kwargs: keyword arguments
                fname: str or list - Filename(s) to save
                imformats: list of str - Image formats ['png', 'pdf', 'svg', etc.]
                crop: bool - Auto-crop whitespace (requires wand library)
                figparam: dict - Figure parameters (see getFigparam)

**Output**:

None (saves file(s) to disk)

**Usage Example**:

```python
>>> p = plotter()
            >>> p.plotProj(ProjType='equalarea')
            >>> 
            >>> # Save single format
            >>> p.figsave(fname='pole_figure.png')
            >>> 
            >>> # Save multiple formats
            >>> p.figsave(
            ...     fname='pole_figure.png',
            ...     imformats=['png', 'pdf', 'svg']
            ... )
            >>> # Creates: pole_figure.png, pole_figure.pdf, pole_figure.svg
            >>> 
            >>> # Save with cropping
            >>> p.figsave(
            ...     fname='cropped_figure.png',
            ...     crop=True
            ... )
            >>> 
            >>> # High DPI save
            >>> p.figparam['dpi'] = 600
            >>> p.figsave(fname='high_res.png')
```

---

## Function: figsaveproc

**Signature**:
```python
def figsaveproc(self, fname, **kwargs):
```

**Description**:

Process and save figure to file.
        Internal method called by figsave() to handle actual file writing.

**Input**:

fname: str - Filename to save
            **kwargs: Additional parameters to update

**Output**:

None (saves file to disk)

**Usage Example**:

```python
>>> # Typically called internally by figsave()
            >>> # But can be used directly:
            >>> p = plotter()
            >>> p.plotProj()
            >>> p.figsaveproc('output.png')
```

---

## Function: format_annot

**Signature**:
```python
def format_annot(self,x, y,**kwargs):
```

**Description**:

Format annotation text for data points.

---

## Function: format_coord

**Signature**:
```python
def format_coord(self,x, y,**kwargs):
```

**Description**:

Format coordinates for display in plot toolbar.

---

## Function: format_coord_test

**Signature**:
```python
def format_coord_test(self,x, y,**kwargs):
```

**Description**:

Test version of coordinate formatting for debugging.

---

## Function: genPoris

**Signature**:
```python
def genPoris(self,**kwargs):
```

**Description**:

Generate pole figure orientation data for given crystal direction.
            Creates dense sampling of orientations where specified crystal direction
            aligns with projection direction, useful for inverse pole figures.

**Input**:

cr_dir: array (3,) - Crystal direction [u, v, w]
                symops: list - Crystal symmetry operations
                resol: int - HEALPix resolution (default: 4)

**Output**:

oris: array (N, 3, 3) - Orientation matrices

**Usage Example**:

```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> 
                >>> # Generate orientations for [001] pole figure
                >>> p = plotter()
                >>> cubic_symops = [np.eye(3)]  # Simplified
                >>> oris_001 = p.genPoris([0, 0, 1], cubic_symops, resol=3)
```

---

## Function: generateSphericalHistSampleData

**Signature**:
```python
def generateSphericalHistSampleData(self,**kwargs):
```

**Description**:

Generate spherical histogram from orientation data.
            Bins orientations on sphere using equal-area binning scheme,
            useful for discrete texture representation.

**Input**:

oris: array (N, 3, 3) - Orientation matrices
                cr_dir: array (3,) - Crystal direction for projection
                bins: int - Number of bins (default: 128)
                histnorm: bool - Normalize histogram (default: True)
                symops: list - Symmetry operations

**Output**:

histogram: array - Binned density values

**Usage Example**:

```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Generate histogram
                >>> orientations = R.random(1000).as_matrix()
                >>> 
                >>> p = plotter()
                >>> hist = p.generateSphericalHistSampleData(
                ...     orientations,
                ...     cr_dir=[0, 0, 1],
                ...     bins=64
                ... )
```

---

## Function: generateSphericalKDESampleData

**Signature**:
```python
def generateSphericalKDESampleData(self,**kwargs):
```

**Description**:

Generate spherical kernel density estimate from orientation data.
            Computes density using Von Mises-Fisher distributions on the sphere,
            accounting for crystal symmetry and bandwidth optimization.

**Input**:

oris: array (N, 3, 3) - Orientation matrices
                cr_dir: array (3,) - Crystal direction for projection
                sdbandwidth: float - KDE bandwidth (default: 0.15)
                sdweights: array - Sample weights (optional)
                symops: list - Symmetry operations

**Output**:

density: array - Density values on sphere

**Usage Example**:

```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Generate KDE from orientations
                >>> orientations = R.random(500).as_matrix()
                >>> 
                >>> p = plotter()
                >>> p.setAttributes(oris=orientations)
                >>> density = p.generateSphericalKDESampleData(
                ...     orientations,
                ...     cr_dir=[0, 0, 1],
                ...     sdbandwidth=0.15
                ... )
```

---

## Function: getColormap

**Signature**:
```python
def getColormap(self,**kwargs):
```

**Description**:

Generate orientation density colormap from sample data.
            Computes orientation density using spherical KDE or histograms,
            interpolates to grid, and prepares colormap data.

**Input**:

oris: array (N, 3, 3) - Orientation matrices
                nump: int - Grid resolution (default: 1001)
                colmapFromSampleData: bool - Generate from data (default: False)
                sdbandwidth: float - KDE bandwidth (default: 0.15)
                sdweights: array - Sample weights (optional)

**Output**:

colmapdata: array - Density values on grid

**Usage Example**:

```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Generate colormap from orientations
                >>> orientations = R.random(1000).as_matrix()
                >>> 
                >>> p = plotter()
                >>> p.setAttributes(oris=orientations)
                >>> colormap_data = p.getColormap(nump=501)
```

---

## Function: getFigparam

**Signature**:
```python
def getFigparam(self, fontsize=None, save=False, phase='A,', figsize=None, **kwargs):
```

**Description**:

Get figure parameter dictionary for saving high-quality figures.

**Input**:

fontsize: float - Font size for text (optional)
            save: bool - If True, use smaller figure size (default: False)
            phase: str - Phase identifier ('A' or other) affects default size
            figsize: tuple - Custom figure size (width, height) in inches
            **kwargs: Additional parameters to override defaults

**Output**:

dict - Figure parameters for plt.savefig()
                'dpi': int - Resolution (default: 300)
                'facecolor': str - Background color
                'edgecolor': str - Edge color
                'orientation': str - 'portrait' or 'landscape'
                'format': str - Image format
                'transparent': bool - Transparent background
                'bbox_inches': str - Bounding box setting
                'pad_inches': float - Padding around figure

**Usage Example**:

```python
>>> p = plotter()
            >>> 
            >>> # Get default parameters
            >>> params = p.getFigparam()
            >>> print("DPI:", params['dpi'])
            >>> 
            >>> # High-resolution PDF
            >>> params_pdf = p.getFigparam(
            ...     save=True,
            ...     format='pdf',
            ...     dpi=600
            ... )
            >>> 
            >>> # Use with savefig
            >>> fig, ax = plt.subplots()
            >>> ax.plot([1, 2, 3], [1, 4, 9])
            >>> fig.savefig('output.png', **params)
```

---

## Function: getScales

**Signature**:
```python
def getScales(self, vmcbar, numticks=None, ticks=None, tickslabels=None, geq=False, leq=False, cmapbins=100, cmapbinsmult=None):
```

**Description**:

Generate color scale parameters including ticks, labels, and colormap.
        Creates a gradient colormap from white→blue→green→yellow→dark red
        with customizable tick positions and labels.

**Input**:

vmcbar: list/array [min, max] - Value range for colorbar
            numticks: int - Number of tick marks (optional, auto-computed if None)
            ticks: array - Explicit tick positions (optional)
            tickslabels: list of str - Custom tick labels (optional)
            geq: bool - Add ≥ symbol to maximum tick (default: False)
            leq: bool - Add ≤ symbol to minimum tick (default: False)
            cmapbins: int - Number of colormap bins (default: 100)
            cmapbinsmult: int - Multiplier for bins based on ticks (optional)

**Output**:

dict - Dictionary containing:
                'tickslabels': list of str - Formatted tick labels
                'vm': list [min, max] - Colormap value range
                'vmbar': list [min, max] - Colorbar value range
                'cmap': matplotlib colormap - Generated colormap
                'ticks': array - Tick positions

**Usage Example**:

```python
>>> p = plotter()
            >>> 
            >>> # Basic scale from 0 to 10
            >>> scales = p.getScales(vmcbar=[0, 10])
            >>> print("Tick labels:", scales['tickslabels'])
            >>> print("Colormap range:", scales['vm'])
            >>> 
            >>> # Custom ticks with inequality symbols
            >>> scales2 = p.getScales(
            ...     vmcbar=[0, 100],
            ...     ticks=np.array([0, 25, 50, 75, 100]),
            ...     geq=True,  # Add ≥ to maximum
            ...     leq=True   # Add ≤ to minimum
            ... )
            >>> print("Labels:", scales2['tickslabels'])
            >>> # Output: ['≤0', '25', '50', '75', '≥100']
            >>> # High-resolution colormap
            >>> scales3 = p.getScales(
            ...     vmcbar=[0, 1],
            ...     cmapbins=1000  # Smooth gradient
            ... )
```

---

## Function: get_cmap

**Signature**:
```python
def get_cmap(colors, nbins=1000, name='my_cmap'):
```

**Description**:

Create a custom colormap from a list of colors with smooth gradients.

**Input**:

colors: list of tuples/colors - Colors to interpolate between
                Can be RGB tuples (R, G, B) or named colors
        nbins: int - Number of discrete bins in colormap (default: 1000)
        name: str - Name for the colormap (default: 'my_cmap')

**Output**:

cmap: matplotlib.colors.LinearSegmentedColormap - Custom colormap

**Usage Example**:

```python
>>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> 
        >>> # Create colormap from white to red to black
        >>> colors = [(1, 1, 1), (1, 0, 0), (0, 0, 0)]
        >>> cmap = get_cmap(colors, nbins=256)
        >>> 
        >>> # Use in plot
        >>> data = np.random.rand(10, 10)
        >>> plt.imshow(data, cmap=cmap)
        >>> plt.colorbar()
        >>> plt.show()
        >>> # Blue-white-red diverging colormap
        >>> colors_div = ['blue', 'white', 'red']
        >>> cmap_div = get_cmap(colors_div, nbins=100)
        >>> 
        >>> # Multi-color gradient
        >>> colors_multi = [
        ...     (0, 0, 1),      # Blue
        ...     (0, 1, 1),      # Cyan
        ...     (0, 1, 0),      # Green
        ...     (1, 1, 0),      # Yellow
        ...     (1, 0, 0)       # Red
        ... ]
        >>> cmap_rainbow = get_cmap(colors_multi, nbins=500)
```

---

## Function: get_colors

**Signature**:
```python
def get_colors(values, cmap, vmin=None, vmax=None, cmin=0, cmax=None, to255=True, nancolor=[0, 0, 0, 1]):
```

**Description**:

Map data values to colors using a colormap.

**Input**:

values: numpy array - Data values to map to colors
        cmap: matplotlib colormap - Colormap to use
        vmin: float - Minimum data value (default: min of values)
        vmax: float - Maximum data value (default: max of values)
        cmin: int - Minimum colormap index (default: 0)
        cmax: int - Maximum colormap index (default: cmap.N)
        to255: bool - Scale colors to 0-255 range (default: True)
        nancolor: list - RGBA color for NaN values (default: black)

**Output**:

Colors: numpy array - Array of colors (RGBA or RGB depending on cmap)
                Shape: (values.shape, 4) if RGBA

**Usage Example**:

```python
>>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> 
        >>> # Map values to colors
        >>> values = np.array([1, 5, 10, 15, 20])
        >>> cmap = plt.cm.viridis
        >>> colors = get_colors(values, cmap, vmin=0, vmax=20, to255=False)
        >>> print("Colors shape:", colors.shape)
        >>> print("First color:", colors[0])
        >>> 
        >>> # Use with scatter plot
        >>> x = np.random.rand(100)
        >>> y = np.random.rand(100)
        >>> z = np.random.rand(100) * 100
        >>> colors = get_colors(z, plt.cm.jet, to255=False)
        >>> plt.scatter(x, y, c=colors)
        >>> plt.show()
        >>> # Handle NaN values
        >>> values_with_nan = np.array([1, 2, np.nan, 4, 5])
        >>> colors_nan = get_colors(values_with_nan, plt.cm.coolwarm, 
        ...                          nancolor=[1, 1, 1, 1])  # White for NaN
```

---

## Function: onclicactivate

**Signature**:
```python
def onclicactivate(self,**kwargs):
```

**Description**:

Activate mouse click event handling for the plot.

---

## Function: onclick

**Signature**:
```python
def onclick(self,event):
```

**Description**:

Mouse click event handler for data selection.

---

## Function: onclick2

**Signature**:
```python
def onclick2(self,event):
```

**Description**:

Alternative click handler for different interaction mode.

---

## Function: onclick3

**Signature**:
```python
def onclick3(self,event):
```

**Description**:

Third click handler variant for specialized interactions.

---

## Function: onmove

**Signature**:
```python
def onmove(self,event):
```

**Description**:

Mouse move event handler for interactive data display.

---

## Function: onpress

**Signature**:
```python
def onpress(self,event):
```

**Description**:

Keyboard press event handler for plot interaction.

---

## Function: onpressActivate

**Signature**:
```python
def onpressActivate(self):
```

**Description**:

Activate keyboard event handling for the plot.

---

## Function: plotColorbar

**Signature**:
```python
def plotColorbar(self,**kwargs):
```

**Description**:

Add colorbar to the plot with custom ticks and labels.
            Creates a horizontal or vertical colorbar with formatted tick labels,
            including support for inequality symbols and custom positioning.

**Input**:

cbartitle: str - Colorbar title (default: "")
                vmbar: list [min, max] - Value range (default: None, auto-detect)
                cbarh: float - Colorbar height (default: 0.04)
                cbarwfac: float - Width factor (default: 0.75)
                cbarhshift: float - Horizontal shift (default: -0.15)
                ticks: array - Tick positions (optional)
                ticklabels: list - Tick labels (optional)

**Output**:

None (adds colorbar to self.fig)

**Usage Example**:

```python
>>> from plotlib import plotter
                >>> import numpy as np
                >>> 
                >>> # Basic colorbar
                >>> p = plotter()
                >>> p.plotProj(ProjType='equalarea')
                >>> orientations = np.random.randn(500, 3, 3)
                >>> p.setAttributes(oris=orientations)
                >>> p.plotColormap()
                >>> p.plotColorbar(cbartitle="MUD")
```

---

## Function: plotColormap

**Signature**:
```python
def plotColormap(self,**kwargs):
```

**Description**:

Plot orientation density colormap with contours on stereographic projection.
            Creates density maps from orientation data using spherical KDE or histograms,
            with contour lines showing texture intensity.

**Input**:

nump: int - Number of points for grid (default: 1001)
                contourcol: str - Contour line color (default: 'k')
                plotmap: bool - Whether to plot colormap (default: True)
                colmapdata: array - Pre-computed colormap data (optional)

**Output**:

None (plots colormap on self.ax)

**Usage Example**:

```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> 
                >>> # Create plotter with orientation data
                >>> p = plotter()
                >>> orientations = np.random.randn(1000, 3, 3)  # Random orientations
                >>> 
                >>> # Setup projection
                >>> p.plotProj(ProjType='equalarea', sphere='half')
                >>> 
                >>> # Set orientation data
                >>> p.setAttributes(oris=orientations)
                >>> 
                >>> # Plot colormap with contours
                >>> p.plotColormap(nump=501, contourcol='black')
                >>> 
                >>> # High-resolution colormap
                >>> p2 = plotter()
                >>> p2.setAttributes(oris=orientations, ProjType='equalarea')
                >>> p2.plotProj()
                >>> p2.plotColormap(nump=1001)
```

---

## Function: plotColormaps

**Signature**:
```python
def plotColormaps(self,sel='all',**kwargs):
```

**Description**:

Plot multiple colormaps in multi-panel figure (typically 2x2).
            Creates comparative visualization of different pole figures,
            useful for showing texture evolution or multi-phase analysis.

**Input**:

Multiple attributes set via setAttributes for each panel

**Output**:

None (creates multi-panel figure)

**Usage Example**:

```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Create 2x2 pole figure comparison
                >>> import matplotlib.pyplot as plt
                >>> fig, axes = plt.subplots(2, 2, figsize=(12, 12))
                >>> 
                >>> # Four different texture states
                >>> oris1 = R.random(500).as_matrix()
                >>> oris2 = R.random(500).as_matrix()
```

---

## Function: plotDirsNorms

**Signature**:
```python
def plotDirsNorms(self,**kwargs):
```

**Description**:

Plot crystallographic directions and plane normals on stereographic projection.
            Displays Miller indices with appropriate notation, supports family notation,
            and handles multiple phases with correspondence relationships.

**Input**:

dirs: list - Crystal directions [[u,v,w], ...] (optional)
                norms: list - Plane normals [[h,k,l], ...] (optional)
                printasfamily: bool - Use family notation {hkl} <uvw> (default: True)
                dhkl: bool - Use decimal Miller indices (default: False)
                Cd: array (3,3) - Direction transformation matrix (default: eye(3))
                Cp: array (3,3) - Plane transformation matrix (default: eye(3))

**Output**:

None (plots directions and normals on self.ax)

**Usage Example**:

```python
>>> from plotlib import plotter
                >>> import numpy as np
                >>> 
                >>> # Basic pole figure with cubic directions
                >>> p = plotter()
                >>> p.plotProj(ProjType='equalarea', sphere='half')
                >>> 
                >>> # Plot low-index directions
                >>> dirs = [[1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,1,1]]
                >>> p.setAttributes(dirs=dirs)
                >>> p.plotDirsNorms()
                >>> 
                >>> # Plot with plane normals
                >>> norms = [[1,1,1], [1,1,0], [1,0,0]]
                >>> p.setAttributes(dirs=dirs, norms=norms)
                >>> p.plotDirsNorms()
```

---

## Function: plotHist

**Signature**:
```python
def plotHist():
```

**Description**:

Plot histogram of orientation data on stereographic projection.
            Creates binned representation of texture with discrete color levels,
            alternative to continuous colormap visualization.

**Input**:

oris: array (N, 3, 3) - Orientation matrices
                bins: int - Number of histogram bins (default: 128)
                histnorm: bool - Normalize histogram (default: True)
                histscale: list [min, max] - Value range (optional)

**Output**:

None (plots histogram on self.ax)

**Usage Example**:

```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Basic histogram
                >>> orientations = R.random(1000).as_matrix()
                >>> 
                >>> p = plotter()
                >>> p.plotProj(ProjType='equalarea', sphere='half')
                >>> p.setAttributes(oris=orientations)
                >>> p.plotHist(bins=64)
```

---

## Function: plotProj

**Signature**:
```python
def plotProj(self, **kwargs):
```

**Description**:

Create a stereographic projection (pole figure) with appropriate grid.
        Sets up the figure, axes, and projection type (equal-area Schmidt net
        or equal-angle Wulff net). Handles full sphere, hemisphere, or standard
        stereographic triangle projections.

**Input**:

**kwargs: keyword arguments
                ProjType: str - 'equalarea' (Schmidt) or 'equalangle' (Wulff)
                sphere: str - 'full', 'half', or 'triangle'
                figsize: tuple - Figure size (width, height) in inches
                stereogrid: bool - Display stereographic grid
                stereoresolution: int - Grid resolution
                stereomesh: bool - Display mesh
                ax: matplotlib axis - Existing axis to plot on (optional)
                fig: matplotlib figure - Existing figure (optional)

**Output**:

None (creates/modifies self.fig and self.ax)

**Usage Example**:

```python
>>> import matplotlib.pyplot as plt
            >>> 
            >>> # Equal-area projection, upper hemisphere
            >>> p = plotter()
            >>> p.plotProj(ProjType='equalarea', sphere='half', figsize=(6, 6))
            >>> plt.show()
            >>> # Equal-angle projection, full sphere
            >>> p2 = plotter()
            >>> p2.plotProj(ProjType='equalangle', sphere='full')
            >>> plt.show()
            >>> # Standard stereographic triangle (for cubic)
            >>> p3 = plotter()
            >>> p3.plotProj(
            ...     ProjType='equalarea',
            ...     sphere='triangle',
            ...     stereogrid=True,
            ...     stereoresolution=10
            ... )
            >>> plt.show()
```

---

## Function: plotScatter

**Signature**:
```python
def plotScatter(self,**kwargs):
```

**Description**:

Create scatter plot of crystallographic orientations or data points.
            Plots data points with color and size scaling based on additional parameters.
            Useful for showing orientation distributions with properties.

**Input**:

scatterdata: array (N, 3, 3) - Orientation matrices to plot
                scattercolscale: array (N,) - Values for color scaling (optional)
                scattersizescale: array (N,) - Values for size scaling (optional)
                scatterfacecolor: str/array - Point colors (default: 'k')
                scatteredgecolor: str - Edge color (default: 'w')
                scatterpointsize: float - Base point size (default: 70)

**Output**:

None (plots scatter points on self.ax)

**Usage Example**:

```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Generate random orientations
                >>> N = 500
                >>> orientations = R.random(N).as_matrix()
                >>> 
                >>> # Create scatter plot
                >>> p = plotter()
                >>> p.plotProj(ProjType='equalarea', sphere='half')
                >>> p.setAttributes(scatterdata=orientations)
                >>> p.processScatterData()
                >>> p.plotScatter()
```

---

## Function: plotScatterAsHist

**Signature**:
```python
def plotScatterAsHist(self,**kwargs):
```

**Description**:

Display scatter plot data as histogram representation.
            Converts scattered orientation points to binned histogram display,
            useful for comparing discrete and continuous representations.

**Input**:

scatterdata: array (N, 3, 3) - Orientation matrices
                bins: int - Number of bins (default: 128)
                histnorm: bool - Normalize (default: True)

**Output**:

None (plots histogram on self.ax)

**Usage Example**:

```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Scatter as histogram
                >>> orientations = R.random(500).as_matrix()
                >>> 
                >>> p = plotter()
                >>> p.plotProj(ProjType='equalarea', sphere='half')
                >>> p.setAttributes(scatterdata=orientations)
                >>> p.plotScatterAsHist(bins=64)
```

---

## Function: plot_atomic_plane2D

**Signature**:
```python
def plot_atomic_plane2D(LatticePoints,normal,vertical,ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],plot=True,salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5],linewidths=[1,1,1],markersizes=[200,200,200], Q=np.eye(3),xlim=[],ylim=[],out=False,zorder=1):
```

**Description**:

plot_atomic_plane2D - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: plot_atomic_plane3D

**Signature**:
```python
def plot_atomic_plane3D(LatticePoints,ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5],Q=np.eye(3)):
```

**Description**:

plot_atomic_plane3D - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: plot_atomlattice2D

**Signature**:
```python
def plot_atomlattice2D(atoms_xyz_position,uvw2xyz,normal,vertical,S=1,R=np.eye(3),ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5]):
```

**Description**:

plot_atomlattice2D - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: plot_cut2D

**Signature**:
```python
def plot_cut2D(ax,Lattice_points,normal,horizontal,vertical,col,alpha=1):
```

**Description**:

Plot 2D cross-section of 3D points.
    Selects points near plane and plots 2D projection.

**Input**:

points (array N×3): 3D points
        cut_plane_normal (array [3]): Cutting plane normal
        cut_plane_point (array [3]): Point on plane
        **kwargs: Plot parameters

**Output**:

fig, ax: Figure and 2D axis

**Usage Example**:

```python
>>> atoms = generate_lattite_atom_positions(L, n1=10, n2=10, n3=10)
        >>> 
        >>> # Cut at z=15
        >>> fig, ax = plot_cut2D(atoms, [0,0,1], [0,0,15])
        >>> ax.set_title('Cross-section at z=15')
        >>> plt.show()
    Notes:
        - Shows 2D slice of 3D structure
        - Useful for interface analysis
        - Can reveal internal structure
```

---

## Function: plot_lattice

**Signature**:
```python
def plot_lattice(Points,LatticeVectors,ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5],Q=np.eye(3),shift=np.zeros(3),atoms=True,linewidth=1,move=np.zeros(3),normal=np.array([0,0,0]),halfspace='upper',s=200,plot=True):
```

**Description**:

plot_lattice - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: plot_lattice2D

**Signature**:
```python
def plot_lattice2D(ax,VV,description,Parentlattice_points,Parent_lattice,\ Productlattice_points,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,linewidth=2,xlim=None):
```

**Description**:

Plot 2D projection of lattice.
    Projects 3D lattice onto 2D plane for simplified visualization.

**Input**:

ax (Axes): 2D matplotlib axis
        L (array 3×3): Lattice matrix
        n1, n2 (int): Cell range
        projection_axis (int): Axis to project onto (0=YZ, 1=XZ, 2=XY)
        **kwargs: Plot parameters

**Output**:

None (modifies axis)

**Usage Example**:

```python
>>> fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        >>> 
        >>> L = lattice_vec(2.95, 2.95, 4.68, 90, 90, 120)
        >>> 
        >>> # Three projections
        >>> plot_lattice2D(axes[0], L, projection_axis=2)  # XY
        >>> axes[0].set_title('XY projection')
        >>> 
        >>> plot_lattice2D(axes[1], L, projection_axis=1)  # XZ
        >>> axes[1].set_title('XZ projection')
        >>> 
        >>> plot_lattice2D(axes[2], L, projection_axis=0)  # YZ
        >>> axes[2].set_title('YZ projection')
    Notes:
        - Shows 2D pattern
        - Useful for symmetry visualization
        - Faster than 3D rendering
```

---

## Function: plot_lattice3D

**Signature**:
```python
def plot_lattice3D(ax,VV,description,Parentlattice_points,Productlattice_points,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,linewidth=2):
```

**Description**:

Complete 3D lattice visualization.
    Plots lattice points and/or edges in single function call.

**Input**:

ax (Axes3D): 3D axis
        L (array 3×3): Lattice matrix
        n1, n2, n3 (int): Cell range
        show_points (bool): Plot lattice points
        show_edges (bool): Plot cell edges
        **kwargs: Plot styling parameters

**Output**:

None (modifies axis)

**Usage Example**:

```python
>>> fig = plt.figure(figsize=(10, 10))
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L = lattice_vec(2.95, 2.95, 4.68, 90, 90, 120)  # HCP
        >>> plot_lattice3D(ax, L, n1=2, n2=2, n3=2, 
        ...                show_points=True, show_edges=True,
        ...                color='blue', s=50)
        >>> 
        >>> ax.set_xlabel('X')
        >>> ax.set_ylabel('Y')
        >>> ax.set_zlabel('Z')
        >>> set_aspect_equal_3d(ax)
        >>> plt.show()
    Notes:
        - Convenience function
        - Combines points and edges
        - Customizable appearance
```

---

## Function: plot_lattice_2Dprojection

**Signature**:
```python
def plot_lattice_2Dprojection(ax,VV,description,Parentlattice_points,Parent_lattice,\ Productlattice_points,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,normals,verticals, linewidth=2,xlim=None):
```

**Description**:

Project lattice along specific crystallographic direction.
    Creates 2D projection viewed along given direction vector.

**Input**:

L (array 3×3): Lattice matrix
        direction (array [3]): Viewing direction
        **kwargs: Plot parameters

**Output**:

None (creates new figure)

**Usage Example**:

```python
>>> L = cubic_lattice_vec(3.0)
        >>> 
        >>> # View along [111]
        >>> plot_lattice_2Dprojection(L, direction=[1,1,1])
        >>> plt.title('View along [111]')
        >>> plt.show()
        >>> 
        >>> # View along [110]
        >>> plot_lattice_2Dprojection(L, direction=[1,1,0])
        >>> plt.title('View along [110]')
    Notes:
        - Arbitrary viewing direction
        - Creates orthogonal projection
        - Shows atomic arrangements
```

---

## Function: plot_lattice_boundaries

**Signature**:
```python
def plot_lattice_boundaries(axl,LatticePointsNew,allPoints=None,polygon=False,tol=1e-1,**kwargs):
```

**Description**:

Plot unit cell boundaries in 3D.
    Draws edges of unit cells to visualize lattice structure.

**Input**:

ax (Axes3D): Matplotlib 3D axis
        L (array 3×3): Lattice matrix
        n1, n2, n3 (int): Number of cells in each direction
        **kwargs: Line plot parameters

**Output**:

None (modifies axis)

**Usage Example**:

```python
>>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L = cubic_lattice_vec(3.0)
        >>> plot_lattice_boundaries(ax, L, n1=2, n2=2, n3=2, color='black')
        >>> 
        >>> set_aspect_equal_3d(ax)
        >>> plt.show()
    Notes:
        - Draws wireframe of unit cells
        - Helps visualize lattice structure
        - Combine with plot_lattice_points for complete view
```

---

## Function: plot_lattice_plane

**Signature**:
```python
def plot_lattice_plane(axl,PlanePoints,**kwargs):
```

**Description**:

Plot a crystallographic plane in 3D lattice.
    Draws plane (hkl) intersecting lattice unit cell using matplotlib 3D.

**Input**:

ax (Axes3D): Matplotlib 3D axis
        L (array 3×3): Lattice matrix
        h, k, l (int): Miller indices of plane
        **kwargs: Additional matplotlib plot parameters (color, alpha, etc.)

**Output**:

None (modifies axis)

**Usage Example**:

```python
>>> import matplotlib.pyplot as plt
        >>> from mpl_toolkits.mplot3d import Axes3D
        >>> import numpy as np
        >>> 
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L = cubic_lattice_vec(3.0)
        >>> plot_lattice_plane(ax, L, 1, 1, 1, color='blue', alpha=0.3)
        >>> plot_lattice_plane(ax, L, 1, 0, 0, color='red', alpha=0.3)
        >>> 
        >>> set_aspect_equal_3d(ax)
        >>> plt.show()
    Notes:
        - Plane intersects unit cell
        - Uses Miller indices for specification
        - Transparency recommended (alpha < 1)
        - Multiple planes can be overlaid
```

---

## Function: plot_lattice_proj

**Signature**:
```python
def plot_lattice_proj(LatticeVectors,normalproj,verticalproj, ax=None, linewidth=2,color='b',eps=1e-1,Q=np.eye(3),Qprojr=np.eye(2),shift=np.zeros(3),shiftproj=0,move=np.zeros(3),normal=np.array([0,0,0]),shifthalfspace=np.zeros(3), halfspace='upper',shiftplot=np.array([0,0]),out=False):
```

**Description**:

plot_lattice_proj - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: plot_latticefaces3D

**Signature**:
```python
def plot_latticefaces3D(ax,Parentlattices,linewidth=2,alpha=0.15,edgecolor='r',linestyle='-',facecolor=(1, 0, 0, 0.15)):
```

**Description**:

plot_latticefaces3D - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: plot_latticesfaces3D

**Signature**:
```python
def plot_latticesfaces3D(ax,VV,description,Parentlattices,Productlattices,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,linewidth=2,alpha=0.15,xlim=[-2,2]):
```

**Description**:

Plot two lattices with transparent faces.
    Visualizes two crystal structures simultaneously for comparison.

**Input**:

ax (Axes3D): 3D axis
        L1, L2 (array 3×3): Lattice matrices
        alpha1, alpha2 (float): Transparency values
        **kwargs: Additional styling

**Output**:

None (modifies axis)

**Usage Example**:

```python
>>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L_parent = cubic_lattice_vec(3.0)
        >>> L_product = tetragonal_lattice_vec(3.0, 3.0, 4.0)
        >>> 
        >>> plot_latticesfaces3D(ax, L_parent, L_product,
        ...                      alpha1=0.2, alpha2=0.2,
        ...                      facecolor1='blue', facecolor2='red')
        >>> 
        >>> ax.set_title('Parent vs Product Phase')
        >>> set_aspect_equal_3d(ax)
    Notes:
        - Overlays two structures
        - Different colors recommended
        - Shows transformation relationship
```

---

## Function: plot_mohr_circles

**Signature**:
```python
def plot_mohr_circles(mcircles,VV,DD,xyz2uvw,scale,xticks=None,yticks=None,ax=None,Parent_lattice='B2'):
```

**Description**:

Plot all three Mohr's circles.
    Visualizes complete 3D strain state via Mohr's circles.

**Input**:

mohr_result (dict): Output from mohr_circles()
        **kwargs: Plot styling parameters

**Output**:

None (creates matplotlib figure)

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> strain = np.array([[0.10,  0.02, 0.01],
        ...                    [0.02, -0.03, 0.00],
        ...                    [0.01,  0.00, -0.05]])
        >>> 
        >>> result = mohr_circles(strain)
        >>> plot_mohr_circles(result)
        >>> plt.title('Mohr Circles - 3D Strain State')
        >>> plt.show()
    Notes:
        - Three circles shown
        - Largest circle envelope
        - Principal strains marked
        - Standard in mechanics
```

---

## Function: plot_planes_on_mohr_circle

**Signature**:
```python
def plot_planes_on_mohr_circle(ax,scale,SelNormalsOnCircle,SelShearsOnCircle,Parent_xyz2hkl, colors,text=False,Parent_lattice='B2'):
```

**Description**:

Plot specific crystallographic planes on Mohr's circle.
    Shows where particular crystal planes plot on Mohr diagram.

**Input**:

mohr_result (dict): From mohr_circles()
        plane_normals (list): List of plane normal vectors
        lattice (array 3×3): Lattice matrix
        **kwargs: Plot parameters

**Output**:

None (adds to current figure)

**Usage Example**:

```python
>>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> 
        >>> # Select planes
        >>> planes = [[1,0,0], [0,1,0], [1,1,0], [1,1,1]]
        >>> 
        >>> plot_mohr_circles(result)
        >>> plot_planes_on_mohr_circle(result, planes, lattice)
        >>> plt.show()
    Notes:
        - Shows strain state on specific planes
        - Useful in transformation analysis
        - Identifies critical planes
```

---

## Function: plot_planes_on_stereotriangle

**Signature**:
```python
def plot_planes_on_stereotriangle(ax,SelNormalsOnCircle,SelShearsOnCircle,Parent_xyz2hkl,colors):
```

**Description**:

Plot planes on stereographic triangle.
    Projects plane normals onto standard stereographic triangle
    for cubic systems.

**Input**:

planes (list): List of (h,k,l) tuples
        **kwargs: Plot parameters

**Output**:

None (creates figure with stereographic triangle)

**Usage Example**:

```python
>>> planes = [(1,0,0), (1,1,0), (1,1,1), (2,1,0)]
        >>> plot_planes_on_stereotriangle(planes)
        >>> plt.title('Low-Index Planes')
        >>> plt.show()
    Notes:
        - Standard triangle for cubic
        - Shows plane distribution
        - Used in texture analysis
```

---

## Function: plot_planes_on_wulffnet

**Signature**:
```python
def plot_planes_on_wulffnet(ax,SelNormalsOnCircle,SelShearsOnCircle,Parent_xyz2hkl,colors):
```

**Description**:

Plot plane normals on Wulff net (equal-angle projection).
    Projects planes onto full Wulff stereographic net.

**Input**:

planes (list): Plane Miller indices
        lattice (array 3×3): Lattice matrix
        **kwargs: Plot parameters

**Output**:

None (creates Wulff net figure)

**Usage Example**:

```python
>>> L = cubic_lattice_vec(3.0)
        >>> planes = [(1,0,0), (0,1,0), (0,0,1), (1,1,1)]
        >>> plot_planes_on_wulffnet(planes, L)
        >>> plt.title('Planes on Wulff Net')
    Notes:
        - Equal-angle projection
        - Preserves angular relationships
        - Full sphere projection
```

---

## Function: plot_points_proj

**Signature**:
```python
def plot_points_proj(Points,normalproj,verticalproj, ax=None, marker="o",markersize=10, color='b',Q=np.eye(3),Qprojr=np.eye(2),shift=np.zeros(3),move=np.zeros(3)):
```

**Description**:

plot_points_proj - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: plot_princip_dir_on_stereotriangle

**Signature**:
```python
def plot_princip_dir_on_stereotriangle(ax,VV,description,markersize=10,markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5):
```

**Description**:

Plot principal strain/stress directions on stereographic triangle.
    Projects principal directions onto standard triangle.

**Input**:

principal_directions (array 3×3): Principal direction matrix
        **kwargs: Plot styling

**Output**:

None (adds to current figure or creates new)

**Usage Example**:

```python
>>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> 
        >>> plot_planes_on_stereotriangle([(1,0,0), (1,1,0), (1,1,1)])
        >>> plot_princip_dir_on_stereotriangle(result['directions'], 
        ...                                      marker='*', s=200, c='red')
        >>> plt.title('Principal Directions')
    Notes:
        - Shows orientation of principal axes
        - Overlays on crystal planes
        - Important for anisotropy analysis
```

---

## Function: plot_princip_dir_on_wulffnet

**Signature**:
```python
def plot_princip_dir_on_wulffnet(ax,VV,description,markersize=10,markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5):
```

**Description**:

Plot principal directions on Wulff net.
    Projects principal strain/stress axes onto Wulff stereographic net.

**Input**:

principal_directions (array 3×3): Principal directions
        **kwargs: Plot parameters

**Output**:

None (creates or modifies figure)

**Usage Example**:

```python
>>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> plot_princip_dir_on_wulffnet(result['directions'])
        >>> plt.title('Principal Strain Axes')
    Notes:
        - Full sphere projection
        - Shows 3D orientation clearly
```

---

## Function: plotcolmap

**Signature**:
```python
def plotcolmap(fname=None, withdraw=False):
```

**Description**:

Plot a single colormap from a saved pickle file.
    Loads and executes plot configuration from pickle file.

**Input**:

fname: str - Path to pickle file (optional)
        withdraw: bool - Withdraw plot from display (default: False)

**Output**:

None (creates matplotlib figure)

**Usage Example**:

```python
>>> plotcolmap('single_pole_figure.pkl')
        >>> # Load without displaying
        >>> plotcolmap('config.pkl', withdraw=True)
```

---

## Function: plotcolmaps

**Signature**:
```python
def plotcolmaps(fname=None, withdraw=False):
```

**Description**:

Plot multiple colormaps from a saved pickle file.
    Creates a 2x2 subplot figure with multiple pole figures/colormaps.
    Loads data from pickle file containing plot configurations.

**Input**:

fname: str - Path to pickle file containing plot data (optional)
        withdraw: bool - Withdraw plot from display (default: False)

**Output**:

None (creates matplotlib figure)

**Usage Example**:

```python
>>> # Assuming you have a saved plot configuration
        >>> plotcolmaps('my_pole_figures.pkl')
        >>> # Creates 2x2 subplot with 4 pole figures
        >>> # With custom display
        >>> plotcolmaps('strain_maps.pkl', withdraw=True)
```

---

## Function: processScatterData

**Signature**:
```python
def processScatterData(self,**kwargs):
```

**Description**:

Process scatter plot data by transforming orientations to projection coordinates.
            Converts orientation matrices to stereographic projection coordinates,
            applies crystal symmetry, and prepares data for scatter plotting.

**Input**:

scatterdata: array (N, 3, 3) - Orientation matrices
                scatteroris: array (N, 3, 3) - Alternative orientation data (optional)
                symops: list of arrays - Symmetry operations (default: [eye(3)])
                R2Proj: array (3, 3) - Projection rotation matrix (default: eye(3))

**Output**:

None (sets self.scatterproj with processed coordinates)

**Usage Example**:

```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Process orientation data
                >>> orientations = R.random(100).as_matrix()
                >>> 
                >>> p = plotter()
                >>> p.setAttributes(scatterdata=orientations)
                >>> p.processScatterData()
                >>> print("Projection coords:", p.scatterproj.shape)
```

---

## Function: scatterDataAnnot

**Signature**:
```python
def scatterDataAnnot(self,**kwargs):
```

**Description**:

Annotate scatter plot data points with information.

---

## Function: setAttributes

**Signature**:
```python
def setAttributes(self, **kwargs):
```

**Description**:

Set or update plotter attributes dynamically.
        This method allows flexible configuration of any plotter attribute
        without needing to know the internal structure.

**Input**:

**kwargs: keyword arguments - Any attribute name and value pairs

**Output**:

None (updates instance attributes)

**Usage Example**:

```python
>>> p = plotter()
            >>> 
            >>> # Set single attribute
            >>> p.setAttributes(cmap='viridis')
            >>> 
            >>> # Set multiple attributes
            >>> p.setAttributes(
            ...     sphere='half',
            ...     ProjType='equalarea',
            ...     figsize=(8, 8),
            ...     dirs=[[1,0,0], [1,1,0]],
            ...     cmap='jet'
            ... )
            >>> 
            >>> # Verify changes
            >>> print("Sphere:", p.sphere)
            >>> print("Colormap:", p.cmap)
```

---

## Function: set_aspect_equal_3d

**Signature**:
```python
def set_aspect_equal_3d(ax):
```

**Description**:

Fix equal aspect ratio bug for 3D matplotlib plots.
    Adjusts axis limits to ensure equal scaling in all three dimensions,
    preventing distortion in 3D visualizations of crystal structures.
    Essential for accurate representation of lattice geometries.

**Input**:

ax (matplotlib.axes._subplots.Axes3DSubplot): 3D axis object from mpl_toolkits.mplot3d

**Output**:

None (modifies axis object in place)

**Usage Example**:

```python
>>> from mpl_toolkits.mplot3d import Axes3D
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> 
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> # Plot crystal lattice points
        >>> points = np.random.rand(100, 3)
        >>> ax.scatter(points[:,0], points[:,1], points[:,2])
        >>> set_aspect_equal_3d(ax)
        >>> plt.show()
    Notes:
        - Must be called AFTER all plotting operations
        - Prevents elongation or compression artifacts
        - Critical for crystallographic accuracy
        - Works by finding maximum range and applying to all axes
```

---

## Function: shiftedColorMap

**Signature**:
```python
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
```

**Description**:

Shift the center of a colormap to a specific data value.
    Useful for data with asymmetric ranges (e.g., -15 to +5) where you want
    the colormap center (typically white or neutral color) at zero.

**Input**:

cmap: matplotlib colormap - The colormap to shift
        start: float - Offset from lowest point [0.0, midpoint] (default: 0.0)
        midpoint: float - New center position [0.0, 1.0] (default: 0.5)
                         Calculate as: 1 - vmax/(vmax + abs(vmin))
        stop: float - Offset from highest point [midpoint, 1.0] (default: 1.0)
        name: str - Name for new colormap (default: 'shiftedcmap')

**Output**:

newcmap: matplotlib.colors.LinearSegmentedColormap - Shifted colormap

**Usage Example**:

```python
>>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> 
        >>> # Data ranging from -15 to +5
        >>> data = np.random.randn(100, 100) * 10 - 5
        >>> vmin, vmax = -15, 5
        >>> 
        >>> # Calculate midpoint to center at zero
        >>> midpoint = 1 - vmax / (vmax + abs(vmin))
        >>> print(f"Midpoint: {midpoint:.3f}")  # 0.75
        >>> 
        >>> # Create shifted colormap
        >>> original_cmap = plt.cm.RdBu_r
        >>> shifted_cmap = shiftedColorMap(original_cmap, midpoint=midpoint)
        >>> 
        >>> # Plot with original vs shifted
        >>> fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        >>> im1 = ax1.imshow(data, cmap=original_cmap, vmin=vmin, vmax=vmax)
        >>> ax1.set_title('Original colormap')
        >>> plt.colorbar(im1, ax=ax1)
        >>> 
        >>> im2 = ax2.imshow(data, cmap=shifted_cmap, vmin=vmin, vmax=vmax)
        >>> ax2.set_title('Shifted colormap (zero = white)')
        >>> plt.colorbar(im2, ax=ax2)
        >>> plt.show()
        >>> # For symmetric data around different center
        >>> # E.g., data from 90 to 110, want center at 100
        >>> vmin, vmax = 90, 110
        >>> midpoint = 1 - (vmax - 100) / (vmax - vmin)  # 0.5
        >>> cmap_temp = shiftedColorMap(plt.cm.coolwarm, midpoint=midpoint)
```

---


# PLOTLIB - Comprehensive Quick Reference

**Module**: `plotlib.py`  
**Functions**: 57  
**Last Updated**: December 19, 2025

This reference provides function signatures with brief descriptions and example usage patterns.

---

## dataAnnot

```python
def dataAnnot(self,**kwargs):
```

Annotate data points with crystallographic information.

---

## dataShow

```python
def dataShow(self,**kwargs):
```

Display data information in plot or console.

---

## figsave

```python
def figsave(self, **kwargs):
```

Save figure to file(s) with optional cropping.
        Supports multiple image formats and automatic cropping of whitespace.
        Can save to multiple formats simultaneously.

**Example**:
```python
>>> p = plotter()
            >>> p.plotProj(ProjType='equalarea')
            >>> 
            >>> # Save single format
# ...
```

---

## figsaveproc

```python
def figsaveproc(self, fname, **kwargs):
```

Process and save figure to file.
        Internal method called by figsave() to handle actual file writing.

**Example**:
```python
>>> # Typically called internally by figsave()
            >>> # But can be used directly:
            >>> p = plotter()
            >>> p.plotProj()
# ...
```

---

## format_annot

```python
def format_annot(self,x, y,**kwargs):
```

Format annotation text for data points.

---

## format_coord

```python
def format_coord(self,x, y,**kwargs):
```

Format coordinates for display in plot toolbar.

---

## format_coord_test

```python
def format_coord_test(self,x, y,**kwargs):
```

Test version of coordinate formatting for debugging.

---

## genPoris

```python
def genPoris(self,**kwargs):
```

Generate pole figure orientation data for given crystal direction.
            Creates dense sampling of orientations where specified crystal direction
            aligns with projection direction, useful for inverse pole figures.

**Example**:
```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> 
                >>> # Generate orientations for [001] pole figure
# ...
```

---

## generateSphericalHistSampleData

```python
def generateSphericalHistSampleData(self,**kwargs):
```

Generate spherical histogram from orientation data.
            Bins orientations on sphere using equal-area binning scheme,
            useful for discrete texture representation.

**Example**:
```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
# ...
```

---

## generateSphericalKDESampleData

```python
def generateSphericalKDESampleData(self,**kwargs):
```

Generate spherical kernel density estimate from orientation data.
            Computes density using Von Mises-Fisher distributions on the sphere,
            accounting for crystal symmetry and bandwidth optimization.

**Example**:
```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
# ...
```

---

## getColormap

```python
def getColormap(self,**kwargs):
```

Generate orientation density colormap from sample data.
            Computes orientation density using spherical KDE or histograms,
            interpolates to grid, and prepares colormap data.

**Example**:
```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
# ...
```

---

## getFigparam

```python
def getFigparam(self, fontsize=None, save=False, phase='A,', figsize=None, **kwargs):
```

Get figure parameter dictionary for saving high-quality figures.

**Example**:
```python
>>> p = plotter()
            >>> 
            >>> # Get default parameters
            >>> params = p.getFigparam()
# ...
```

---

## getScales

```python
def getScales(self, vmcbar, numticks=None, ticks=None, tickslabels=None, geq=False, leq=False, cmapbins=100, cmapbinsmult=None):
```

Generate color scale parameters including ticks, labels, and colormap.
        Creates a gradient colormap from white→blue→green→yellow→dark red
        with customizable tick positions and labels.

**Example**:
```python
>>> p = plotter()
            >>> 
            >>> # Basic scale from 0 to 10
            >>> scales = p.getScales(vmcbar=[0, 10])
# ...
```

---

## get_cmap

```python
def get_cmap(colors, nbins=1000, name='my_cmap'):
```

Create a custom colormap from a list of colors with smooth gradients.

**Example**:
```python
>>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> 
        >>> # Create colormap from white to red to black
# ...
```

---

## get_colors

```python
def get_colors(values, cmap, vmin=None, vmax=None, cmin=0, cmax=None, to255=True, nancolor=[0, 0, 0, 1]):
```

Map data values to colors using a colormap.

**Example**:
```python
>>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> 
        >>> # Map values to colors
# ...
```

---

## onclicactivate

```python
def onclicactivate(self,**kwargs):
```

Activate mouse click event handling for the plot.

---

## onclick

```python
def onclick(self,event):
```

Mouse click event handler for data selection.

---

## onclick2

```python
def onclick2(self,event):
```

Alternative click handler for different interaction mode.

---

## onclick3

```python
def onclick3(self,event):
```

Third click handler variant for specialized interactions.

---

## onmove

```python
def onmove(self,event):
```

Mouse move event handler for interactive data display.

---

## onpress

```python
def onpress(self,event):
```

Keyboard press event handler for plot interaction.

---

## onpressActivate

```python
def onpressActivate(self):
```

Activate keyboard event handling for the plot.

---

## plotColorbar

```python
def plotColorbar(self,**kwargs):
```

Add colorbar to the plot with custom ticks and labels.
            Creates a horizontal or vertical colorbar with formatted tick labels,
            including support for inequality symbols and custom positioning.

**Example**:
```python
>>> from plotlib import plotter
                >>> import numpy as np
                >>> 
                >>> # Basic colorbar
# ...
```

---

## plotColormap

```python
def plotColormap(self,**kwargs):
```

Plot orientation density colormap with contours on stereographic projection.
            Creates density maps from orientation data using spherical KDE or histograms,
            with contour lines showing texture intensity.

**Example**:
```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> 
                >>> # Create plotter with orientation data
# ...
```

---

## plotColormaps

```python
def plotColormaps(self,sel='all',**kwargs):
```

Plot multiple colormaps in multi-panel figure (typically 2x2).
            Creates comparative visualization of different pole figures,
            useful for showing texture evolution or multi-phase analysis.

**Example**:
```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
# ...
```

---

## plotDirsNorms

```python
def plotDirsNorms(self,**kwargs):
```

Plot crystallographic directions and plane normals on stereographic projection.
            Displays Miller indices with appropriate notation, supports family notation,
            and handles multiple phases with correspondence relationships.

**Example**:
```python
>>> from plotlib import plotter
                >>> import numpy as np
                >>> 
                >>> # Basic pole figure with cubic directions
# ...
```

---

## plotHist

```python
def plotHist():
```

Plot histogram of orientation data on stereographic projection.
            Creates binned representation of texture with discrete color levels,
            alternative to continuous colormap visualization.

**Example**:
```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
# ...
```

---

## plotProj

```python
def plotProj(self, **kwargs):
```

Create a stereographic projection (pole figure) with appropriate grid.
        Sets up the figure, axes, and projection type (equal-area Schmidt net
        or equal-angle Wulff net). Handles full sphere, hemisphere, or standard

**Example**:
```python
>>> import matplotlib.pyplot as plt
            >>> 
            >>> # Equal-area projection, upper hemisphere
            >>> p = plotter()
# ...
```

---

## plotScatter

```python
def plotScatter(self,**kwargs):
```

Create scatter plot of crystallographic orientations or data points.
            Plots data points with color and size scaling based on additional parameters.
            Useful for showing orientation distributions with properties.

**Example**:
```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
# ...
```

---

## plotScatterAsHist

```python
def plotScatterAsHist(self,**kwargs):
```

Display scatter plot data as histogram representation.
            Converts scattered orientation points to binned histogram display,
            useful for comparing discrete and continuous representations.

**Example**:
```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
# ...
```

---

## plot_atomic_plane2D

```python
def plot_atomic_plane2D(LatticePoints,normal,vertical,ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],plot=True,salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5],linewidths=[1,1,1],markersizes=[200,200,200], Q=np.eye(3),xlim=[],ylim=[],out=False,zorder=1):
```

plot_atomic_plane2D - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## plot_atomic_plane3D

```python
def plot_atomic_plane3D(LatticePoints,ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5],Q=np.eye(3)):
```

plot_atomic_plane3D - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## plot_atomlattice2D

```python
def plot_atomlattice2D(atoms_xyz_position,uvw2xyz,normal,vertical,S=1,R=np.eye(3),ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5]):
```

plot_atomlattice2D - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## plot_cut2D

```python
def plot_cut2D(ax,Lattice_points,normal,horizontal,vertical,col,alpha=1):
```

Plot 2D cross-section of 3D points.
    Selects points near plane and plots 2D projection.

**Example**:
```python
>>> atoms = generate_lattite_atom_positions(L, n1=10, n2=10, n3=10)
        >>> 
        >>> # Cut at z=15
        >>> fig, ax = plot_cut2D(atoms, [0,0,1], [0,0,15])
# ...
```

---

## plot_lattice

```python
def plot_lattice(Points,LatticeVectors,ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5],Q=np.eye(3),shift=np.zeros(3),atoms=True,linewidth=1,move=np.zeros(3),normal=np.array([0,0,0]),halfspace='upper',s=200,plot=True):
```

plot_lattice - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## plot_lattice2D

```python
def plot_lattice2D(ax,VV,description,Parentlattice_points,Parent_lattice,\ Productlattice_points,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,linewidth=2,xlim=None):
```

Plot 2D projection of lattice.
    Projects 3D lattice onto 2D plane for simplified visualization.

**Example**:
```python
>>> fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        >>> 
        >>> L = lattice_vec(2.95, 2.95, 4.68, 90, 90, 120)
        >>> 
# ...
```

---

## plot_lattice3D

```python
def plot_lattice3D(ax,VV,description,Parentlattice_points,Productlattice_points,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,linewidth=2):
```

Complete 3D lattice visualization.
    Plots lattice points and/or edges in single function call.

**Example**:
```python
>>> fig = plt.figure(figsize=(10, 10))
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L = lattice_vec(2.95, 2.95, 4.68, 90, 90, 120)  # HCP
# ...
```

---

## plot_lattice_2Dprojection

```python
def plot_lattice_2Dprojection(ax,VV,description,Parentlattice_points,Parent_lattice,\ Productlattice_points,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,normals,verticals, linewidth=2,xlim=None):
```

Project lattice along specific crystallographic direction.
    Creates 2D projection viewed along given direction vector.

**Example**:
```python
>>> L = cubic_lattice_vec(3.0)
        >>> 
        >>> # View along [111]
        >>> plot_lattice_2Dprojection(L, direction=[1,1,1])
# ...
```

---

## plot_lattice_boundaries

```python
def plot_lattice_boundaries(axl,LatticePointsNew,allPoints=None,polygon=False,tol=1e-1,**kwargs):
```

Plot unit cell boundaries in 3D.
    Draws edges of unit cells to visualize lattice structure.

**Example**:
```python
>>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L = cubic_lattice_vec(3.0)
# ...
```

---

## plot_lattice_plane

```python
def plot_lattice_plane(axl,PlanePoints,**kwargs):
```

Plot a crystallographic plane in 3D lattice.
    Draws plane (hkl) intersecting lattice unit cell using matplotlib 3D.

**Example**:
```python
>>> import matplotlib.pyplot as plt
        >>> from mpl_toolkits.mplot3d import Axes3D
        >>> import numpy as np
        >>> 
# ...
```

---

## plot_lattice_proj

```python
def plot_lattice_proj(LatticeVectors,normalproj,verticalproj, ax=None, linewidth=2,color='b',eps=1e-1,Q=np.eye(3),Qprojr=np.eye(2),shift=np.zeros(3),shiftproj=0,move=np.zeros(3),normal=np.array([0,0,0]),shifthalfspace=np.zeros(3), halfspace='upper',shiftplot=np.array([0,0]),out=False):
```

plot_lattice_proj - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## plot_latticefaces3D

```python
def plot_latticefaces3D(ax,Parentlattices,linewidth=2,alpha=0.15,edgecolor='r',linestyle='-',facecolor=(1, 0, 0, 0.15)):
```

plot_latticefaces3D - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## plot_latticesfaces3D

```python
def plot_latticesfaces3D(ax,VV,description,Parentlattices,Productlattices,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,linewidth=2,alpha=0.15,xlim=[-2,2]):
```

Plot two lattices with transparent faces.
    Visualizes two crystal structures simultaneously for comparison.

**Example**:
```python
>>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L_parent = cubic_lattice_vec(3.0)
# ...
```

---

## plot_mohr_circles

```python
def plot_mohr_circles(mcircles,VV,DD,xyz2uvw,scale,xticks=None,yticks=None,ax=None,Parent_lattice='B2'):
```

Plot all three Mohr's circles.
    Visualizes complete 3D strain state via Mohr's circles.

**Example**:
```python
>>> import numpy as np
        >>> 
        >>> strain = np.array([[0.10,  0.02, 0.01],
        ...                    [0.02, -0.03, 0.00],
# ...
```

---

## plot_planes_on_mohr_circle

```python
def plot_planes_on_mohr_circle(ax,scale,SelNormalsOnCircle,SelShearsOnCircle,Parent_xyz2hkl, colors,text=False,Parent_lattice='B2'):
```

Plot specific crystallographic planes on Mohr's circle.
    Shows where particular crystal planes plot on Mohr diagram.

**Example**:
```python
>>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> 
        >>> # Select planes
# ...
```

---

## plot_planes_on_stereotriangle

```python
def plot_planes_on_stereotriangle(ax,SelNormalsOnCircle,SelShearsOnCircle,Parent_xyz2hkl,colors):
```

Plot planes on stereographic triangle.
    Projects plane normals onto standard stereographic triangle
    for cubic systems.

**Example**:
```python
>>> planes = [(1,0,0), (1,1,0), (1,1,1), (2,1,0)]
        >>> plot_planes_on_stereotriangle(planes)
        >>> plt.title('Low-Index Planes')
        >>> plt.show()
# ...
```

---

## plot_planes_on_wulffnet

```python
def plot_planes_on_wulffnet(ax,SelNormalsOnCircle,SelShearsOnCircle,Parent_xyz2hkl,colors):
```

Plot plane normals on Wulff net (equal-angle projection).
    Projects planes onto full Wulff stereographic net.

**Example**:
```python
>>> L = cubic_lattice_vec(3.0)
        >>> planes = [(1,0,0), (0,1,0), (0,0,1), (1,1,1)]
        >>> plot_planes_on_wulffnet(planes, L)
        >>> plt.title('Planes on Wulff Net')
# ...
```

---

## plot_points_proj

```python
def plot_points_proj(Points,normalproj,verticalproj, ax=None, marker="o",markersize=10, color='b',Q=np.eye(3),Qprojr=np.eye(2),shift=np.zeros(3),move=np.zeros(3)):
```

plot_points_proj - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## plot_princip_dir_on_stereotriangle

```python
def plot_princip_dir_on_stereotriangle(ax,VV,description,markersize=10,markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5):
```

Plot principal strain/stress directions on stereographic triangle.
    Projects principal directions onto standard triangle.

**Example**:
```python
>>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> 
        >>> plot_planes_on_stereotriangle([(1,0,0), (1,1,0), (1,1,1)])
# ...
```

---

## plot_princip_dir_on_wulffnet

```python
def plot_princip_dir_on_wulffnet(ax,VV,description,markersize=10,markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5):
```

Plot principal directions on Wulff net.
    Projects principal strain/stress axes onto Wulff stereographic net.

**Example**:
```python
>>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> plot_princip_dir_on_wulffnet(result['directions'])
        >>> plt.title('Principal Strain Axes')
# ...
```

---

## plotcolmap

```python
def plotcolmap(fname=None, withdraw=False):
```

Plot a single colormap from a saved pickle file.
    Loads and executes plot configuration from pickle file.

**Example**:
```python
>>> plotcolmap('single_pole_figure.pkl')
        >>> # Load without displaying
        >>> plotcolmap('config.pkl', withdraw=True)
```

---

## plotcolmaps

```python
def plotcolmaps(fname=None, withdraw=False):
```

Plot multiple colormaps from a saved pickle file.
    Creates a 2x2 subplot figure with multiple pole figures/colormaps.
    Loads data from pickle file containing plot configurations.

**Example**:
```python
>>> # Assuming you have a saved plot configuration
        >>> plotcolmaps('my_pole_figures.pkl')
        >>> # Creates 2x2 subplot with 4 pole figures
        >>> # With custom display
# ...
```

---

## processScatterData

```python
def processScatterData(self,**kwargs):
```

Process scatter plot data by transforming orientations to projection coordinates.
            Converts orientation matrices to stereographic projection coordinates,
            applies crystal symmetry, and prepares data for scatter plotting.

**Example**:
```python
>>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
# ...
```

---

## scatterDataAnnot

```python
def scatterDataAnnot(self,**kwargs):
```

Annotate scatter plot data points with information.

---

## setAttributes

```python
def setAttributes(self, **kwargs):
```

Set or update plotter attributes dynamically.
        This method allows flexible configuration of any plotter attribute
        without needing to know the internal structure.

**Example**:
```python
>>> p = plotter()
            >>> 
            >>> # Set single attribute
            >>> p.setAttributes(cmap='viridis')
# ...
```

---

## set_aspect_equal_3d

```python
def set_aspect_equal_3d(ax):
```

Fix equal aspect ratio bug for 3D matplotlib plots.
    Adjusts axis limits to ensure equal scaling in all three dimensions,
    preventing distortion in 3D visualizations of crystal structures.

**Example**:
```python
>>> from mpl_toolkits.mplot3d import Axes3D
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> 
# ...
```

---

## shiftedColorMap

```python
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
```

Shift the center of a colormap to a specific data value.
    Useful for data with asymmetric ranges (e.g., -15 to +5) where you want
    the colormap center (typically white or neutral color) at zero.

**Example**:
```python
>>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> 
        >>> # Data ranging from -15 to +5
# ...
```

---


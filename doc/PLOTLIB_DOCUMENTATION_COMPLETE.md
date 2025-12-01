# PLOTLIB - COMPLETE DOCUMENTATION
## Crystallographic Plotting and Visualization Library

**Module**: `plotlib.py` (Extended Version)  
**Total Functions**: 30  
**Version**: Extended 2024  
**Status**: Production Ready  

---

## 📚 TABLE OF CONTENTS

1. [get_cmap](#get_cmap)
2. [get_colors](#get_colors)
3. [plotcolmaps](#plotcolmaps)
4. [plotcolmap](#plotcolmap)
5. [shiftedColorMap](#shiftedcolormap)
6. [mohr_circles](#mohr_circles)
7. [plot_mohr_circles](#plot_mohr_circles)
8. [plot_planes_on_mohr_circle](#plot_planes_on_mohr_circle)
9. [write_mohr_planes](#write_mohr_planes)
10. [strains_along_13mohrcirle](#strains_along_13mohrcirle)
11. [zero_normal_strains](#zero_normal_strains)
12. [plot_lattice3D](#plot_lattice3d)
13. [plot_latticefaces3D](#plot_latticefaces3d)
14. [plot_latticesfaces3D](#plot_latticesfaces3d)
15. [set_aspect_equal_3d](#set_aspect_equal_3d)
16. [plot_lattice2D](#plot_lattice2d)
17. [plot_lattice_2Dprojection](#plot_lattice_2dprojection)
18. [plot_lattice](#plot_lattice)
19. [plot_lattice_proj](#plot_lattice_proj)
20. [plot_points_proj](#plot_points_proj)
21. [plot_lattice_plane](#plot_lattice_plane)
22. [plot_lattice_boundaries](#plot_lattice_boundaries)
23. [plot_atomic_plane2D](#plot_atomic_plane2d)
24. [plot_atomic_plane3D](#plot_atomic_plane3d)
25. [plot_atomlattice2D](#plot_atomlattice2d)
26. [plot_cut2D](#plot_cut2d)
27. [plot_planes_on_stereotriangle](#plot_planes_on_stereotriangle)
28. [plot_planes_on_wulffnet](#plot_planes_on_wulffnet)
29. [plot_princip_dir_on_stereotriangle](#plot_princip_dir_on_stereotriangle)
30. [plot_princip_dir_on_wulffnet](#plot_princip_dir_on_wulffnet)

---

## 1. `get_cmap`

### Signature
```python
def get_cmap(colors, nbins=1000, name='my_cmap')
```

### Description
Create a custom colormap from a list of colors with smooth gradients.

### Input Parameters
colors: list of tuples/colors - Colors to interpolate between
                Can be RGB tuples (R, G, B) or named colors
        nbins: int - Number of discrete bins in colormap (default: 1000)
        name: str - Name for the colormap (default: 'my_cmap')

### Output
cmap: matplotlib.colors.LinearSegmentedColormap - Custom colormap

### Usage Examples
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

## 2. `get_colors`

### Signature
```python
def get_colors(values, cmap, vmin=None, vmax=None, cmin=0, cmax=None, 
               to255=True, nancolor=[0, 0, 0, 1])
```

### Description
Map data values to colors using a colormap.

### Input Parameters
values: numpy array - Data values to map to colors
        cmap: matplotlib colormap - Colormap to use
        vmin: float - Minimum data value (default: min of values)
        vmax: float - Maximum data value (default: max of values)
        cmin: int - Minimum colormap index (default: 0)
        cmax: int - Maximum colormap index (default: cmap.N)
        to255: bool - Scale colors to 0-255 range (default: True)
        nancolor: list - RGBA color for NaN values (default: black)

### Output
Colors: numpy array - Array of colors (RGBA or RGB depending on cmap)
                Shape: (values.shape, 4) if RGBA

### Usage Examples
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

## 3. `plotcolmaps`

### Signature
```python
def plotcolmaps(fname=None, withdraw=False)
```

### Description
Plot multiple colormaps from a saved pickle file.
    Creates a 2x2 subplot figure with multiple pole figures/colormaps.
    Loads data from pickle file containing plot configurations.

### Input Parameters
fname: str - Path to pickle file containing plot data (optional)
        withdraw: bool - Withdraw plot from display (default: False)

### Output
None (creates matplotlib figure)

### Usage Examples
```python
>>> # Assuming you have a saved plot configuration
        >>> plotcolmaps('my_pole_figures.pkl')
        >>> # Creates 2x2 subplot with 4 pole figures
        >>> # With custom display
        >>> plotcolmaps('strain_maps.pkl', withdraw=True)
```

---

## 4. `plotcolmap`

### Signature
```python
def plotcolmap(fname=None, withdraw=False)
```

### Description
Plot a single colormap from a saved pickle file.
    Loads and executes plot configuration from pickle file.

### Input Parameters
fname: str - Path to pickle file (optional)
        withdraw: bool - Withdraw plot from display (default: False)

### Output
None (creates matplotlib figure)

### Usage Examples
```python
>>> plotcolmap('single_pole_figure.pkl')
        >>> # Load without displaying
        >>> plotcolmap('config.pkl', withdraw=True)
```

---

## 5. `shiftedColorMap`

### Signature
```python
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap')
```

### Description
Shift the center of a colormap to a specific data value.
    Useful for data with asymmetric ranges (e.g., -15 to +5) where you want
    the colormap center (typically white or neutral color) at zero.

### Input Parameters
cmap: matplotlib colormap - The colormap to shift
        start: float - Offset from lowest point [0.0, midpoint] (default: 0.0)
        midpoint: float - New center position [0.0, 1.0] (default: 0.5)
                         Calculate as: 1 - vmax/(vmax + abs(vmin))
        stop: float - Offset from highest point [midpoint, 1.0] (default: 1.0)
        name: str - Name for new colormap (default: 'shiftedcmap')

### Output
newcmap: matplotlib.colors.LinearSegmentedColormap - Shifted colormap

### Usage Examples
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

## 6. `mohr_circles`

### Signature
```python
def mohr_circles(strain_tensor)
```

### Description
Calculate Mohr's circles from strain or stress tensor.
    Computes principal strains/stresses and parameters for three Mohr's
    circles. Used for visualizing 3D strain state in 2D.

### Input Parameters
strain_tensor (array 3×3): Symmetric strain or stress tensor

### Output
dict: {
            'principal': array [3] - principal values (sorted)
            'directions': array 3×3 - principal directions
            'circles': list of 3 tuples - (center, radius) for each circle
        }

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # Uniaxial strain state
        >>> strain = np.array([[0.1,  0.0, 0.0],
        ...                    [0.0, -0.03, 0.0],
        ...                    [0.0,  0.0, -0.03]])
        >>> 
        >>> result = mohr_circles(strain)
        >>> print("Principal strains:", result['principal'])
        >>> print("Circle 1 (max):", result['circles'][0])
        >>> 
        >>> # Visualize
        >>> plot_mohr_circles(result)
```

### Notes
- Three circles for 3D state
        - Largest circle: (ε₁ - ε₃)/2
        - Used in failure analysis
        - Shows all possible strain states on planes

### Formula
Circle i: center = (σᵢ + σⱼ)/2, radius = |σᵢ - σⱼ|/2
        Three circles: 1-2, 2-3, 1-3 planes

---

## 7. `plot_mohr_circles`

### Signature
```python
def plot_mohr_circles(mohr_result, **kwargs)
```

### Description
Plot all three Mohr's circles.
    Visualizes complete 3D strain state via Mohr's circles.

### Input Parameters
mohr_result (dict): Output from mohr_circles()
        **kwargs: Plot styling parameters

### Output
None (creates matplotlib figure)

### Usage Examples
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
```

### Notes
- Three circles shown
        - Largest circle envelope
        - Principal strains marked
        - Standard in mechanics

---

## 8. `plot_planes_on_mohr_circle`

### Signature
```python
def plot_planes_on_mohr_circle(mohr_result, plane_normals, lattice, **kwargs)
```

### Description
Plot specific crystallographic planes on Mohr's circle.
    Shows where particular crystal planes plot on Mohr diagram.

### Input Parameters
mohr_result (dict): From mohr_circles()
        plane_normals (list): List of plane normal vectors
        lattice (array 3×3): Lattice matrix
        **kwargs: Plot parameters

### Output
None (adds to current figure)

### Usage Examples
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
```

### Notes
- Shows strain state on specific planes
        - Useful in transformation analysis
        - Identifies critical planes

---

## 9. `write_mohr_planes`

### Signature
```python
def write_mohr_planes(filename, mohr_result, planes, lattice)
```

### Description
Write Mohr circle analysis results to file.
    Saves principal strains, circles, and plane-specific strains.

### Input Parameters
filename (str): Output file path
        mohr_result (dict): From mohr_circles()
        planes (list): Plane normals analyzed
        lattice (array 3×3): Lattice matrix

### Output
None (writes to file)

### Usage Examples
```python
>>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> planes = [[1,0,0], [1,1,0], [1,1,1]]
        >>> L = cubic_lattice_vec(3.0)
        >>> 
        >>> write_mohr_planes('mohr_analysis.txt', result, planes, L)
        >>> # Creates text file with complete analysis
```

### Notes
- Formatted text output
        - Includes principal values
        - Lists strain on each plane
        - Suitable for reports

---

## 10. `strains_along_13mohrcirle`

### Signature
```python
def strains_along_13mohrcirle(mohr_result, n_points=100)
```

### Description
Calculate strains along maximum Mohr's circle.
    Computes normal and shear strain components for points on the
    largest Mohr's circle (ε₁-ε₃ circle).

### Input Parameters
mohr_result (dict): Output from mohr_circles()
        n_points (int): Number of points to sample (default: 100)

### Output
dict: {
            'normal': array - normal strain values
            'shear': array - shear strain values
            'angles': array - rotation angles (radians)
        }

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> strain = np.diag([0.1, 0.0, -0.05])
        >>> mohr_result = mohr_circles(strain)
        >>> 
        >>> points = strains_along_13mohrcirle(mohr_result, n_points=50)
        >>> 
        >>> # Plot
        >>> plt.plot(points['normal'], points['shear'])
        >>> plt.xlabel('Normal strain')
        >>> plt.ylabel('Shear strain')
        >>> plt.axis('equal')
        >>> plt.title('Maximum Mohr Circle')
```

### Notes
- Maximum circle between ε₁ and ε₃
        - Shows all possible strain states
        - Used in failure analysis

### Formula
εₙ = center + radius·cos(2θ)
        γ/2 = radius·sin(2θ)

---

## 11. `zero_normal_strains`

### Signature
```python
def zero_normal_strains(strain_tensor)
```

### Description
Find planes with zero normal strain in given strain state.
    Calculates planes where normal strain εₙ = nᵀ·ε·n = 0.
    These planes experience only shear.

### Input Parameters
strain_tensor (array 3×3): Strain tensor

### Output
list: List of plane normals with zero normal strain

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # Pure shear state
        >>> strain = np.array([[0.1,  0.05, 0.0],
        ...                    [0.05, -0.1, 0.0],
        ...                    [0.0,  0.0,  0.0]])
        >>> 
        >>> planes = zero_normal_strains(strain)
        >>> print(f"Found {len(planes)} planes with zero normal strain")
        >>> for p in planes:
        ...     print(f"Plane: {p}")
```

### Notes
- Important in transformation theory
        - Related to habit planes
        - Used in invariant plane strain analysis

### Formula
εₙ = nᵀ·ε·n = 0
        Solve eigenvalue problem

---

## 12. `plot_lattice3D`

### Signature
```python
def plot_lattice3D(ax, L, n1=1, n2=1, n3=1, show_points=True, show_edges=True, **kwargs)
```

### Description
Complete 3D lattice visualization.
    Plots lattice points and/or edges in single function call.

### Input Parameters
ax (Axes3D): 3D axis
        L (array 3×3): Lattice matrix
        n1, n2, n3 (int): Cell range
        show_points (bool): Plot lattice points
        show_edges (bool): Plot cell edges
        **kwargs: Plot styling parameters

### Output
None (modifies axis)

### Usage Examples
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
```

### Notes
- Convenience function
        - Combines points and edges
        - Customizable appearance

---

## 13. `plot_latticefaces3D`

### Signature
```python
def plot_latticefaces3D(ax, L, alpha=0.3, **kwargs)
```

### Description
Plot lattice unit cell faces with transparency.
    Visualizes unit cell as transparent parallelepiped.

### Input Parameters
ax (Axes3D): 3D axis
        L (array 3×3): Lattice matrix
        alpha (float): Transparency (0-1)
        **kwargs: Poly3DCollection parameters

### Output
None (modifies axis)

### Usage Examples
```python
>>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L = monoclinic_lattice_vec(2.9, 4.1, 4.6, 97)
        >>> plot_latticefaces3D(ax, L, alpha=0.2, facecolor='cyan', edgecolor='black')
        >>> 
        >>> set_aspect_equal_3d(ax)
        >>> plt.show()
```

### Notes
- Shows 3D cell shape clearly
        - Transparency recommended
        - Can overlay multiple cells

---

## 14. `plot_latticesfaces3D`

### Signature
```python
def plot_latticesfaces3D(ax, L1, L2, alpha1=0.2, alpha2=0.2, **kwargs)
```

### Description
Plot two lattices with transparent faces.
    Visualizes two crystal structures simultaneously for comparison.

### Input Parameters
ax (Axes3D): 3D axis
        L1, L2 (array 3×3): Lattice matrices
        alpha1, alpha2 (float): Transparency values
        **kwargs: Additional styling

### Output
None (modifies axis)

### Usage Examples
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
```

### Notes
- Overlays two structures
        - Different colors recommended
        - Shows transformation relationship

---

## 15. `set_aspect_equal_3d`

### Signature
```python
def set_aspect_equal_3d(ax)
```

### Description
Fix equal aspect ratio bug for 3D matplotlib plots.
    Adjusts axis limits to ensure equal scaling in all three dimensions,
    preventing distortion in 3D visualizations of crystal structures.
    Essential for accurate representation of lattice geometries.

### Input Parameters
ax (matplotlib.axes._subplots.Axes3DSubplot): 3D axis object from mpl_toolkits.mplot3d

### Output
None (modifies axis object in place)

### Usage Examples
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
```

### Notes
- Must be called AFTER all plotting operations
        - Prevents elongation or compression artifacts
        - Critical for crystallographic accuracy
        - Works by finding maximum range and applying to all axes

---

## 16. `plot_lattice2D`

### Signature
```python
def plot_lattice2D(ax, L, n1=2, n2=2, projection_axis=2, **kwargs)
```

### Description
Plot 2D projection of lattice.
    Projects 3D lattice onto 2D plane for simplified visualization.

### Input Parameters
ax (Axes): 2D matplotlib axis
        L (array 3×3): Lattice matrix
        n1, n2 (int): Cell range
        projection_axis (int): Axis to project onto (0=YZ, 1=XZ, 2=XY)
        **kwargs: Plot parameters

### Output
None (modifies axis)

### Usage Examples
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
```

### Notes
- Shows 2D pattern
        - Useful for symmetry visualization
        - Faster than 3D rendering

---

## 17. `plot_lattice_2Dprojection`

### Signature
```python
def plot_lattice_2Dprojection(L, direction=[0,0,1], **kwargs)
```

### Description
Project lattice along specific crystallographic direction.
    Creates 2D projection viewed along given direction vector.

### Input Parameters
L (array 3×3): Lattice matrix
        direction (array [3]): Viewing direction
        **kwargs: Plot parameters

### Output
None (creates new figure)

### Usage Examples
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
```

### Notes
- Arbitrary viewing direction
        - Creates orthogonal projection
        - Shows atomic arrangements

---

## 18. `plot_lattice`

### Signature
```python
def plot_lattice(L, n1=1, n2=1, n3=1, basis=None, ax=None, **kwargs)
```

### Description
Complete lattice plotting function with atoms and unit cells.
    High-level function combining atomic positions, cell boundaries,
    and styling in single call.

### Input Parameters
L (array 3×3): Lattice matrix
        n1, n2, n3 (int): Cell range
        basis (list): Atomic basis
        ax (Axes3D): 3D axis (creates new if None)
        **kwargs: Styling parameters

### Output
fig, ax: Matplotlib figure and axis objects

### Usage Examples
```python
>>> # Plot BCC structure
        >>> L = cubic_lattice_vec(2.87)  # Fe
        >>> basis = [[0,0,0], [0.5,0.5,0.5]]
        >>> fig, ax = plot_lattice(L, n1=2, n2=2, n3=2, basis=basis,
        ...                        atom_color='blue', atom_size=50)
        >>> ax.set_title('BCC Iron')
        >>> plt.show()
```

### Notes
- All-in-one plotting function
        - Combines multiple visualization elements
        - Customizable appearance

---

## 19. `plot_lattice_proj`

### Signature
```python
def plot_lattice_proj(L, direction, n1=2, n2=2, basis=None, **kwargs)
```

### Description
Plot lattice projection along specific direction.
    Creates 2D projection of 3D lattice viewed along given direction.

### Input Parameters
L (array 3×3): Lattice matrix
        direction (array [3]): Viewing direction
        n1, n2 (int): Cell range
        basis (list): Atomic basis
        **kwargs: Plot parameters

### Output
fig, ax: Figure and axis

### Usage Examples
```python
>>> L = cubic_lattice_vec(3.0)
        >>> # View along [111]
        >>> fig, ax = plot_lattice_proj(L, [1,1,1], n1=3, n2=3)
        >>> ax.set_title('Cubic lattice - [111] projection')
        >>> plt.show()
```

### Notes
- Shows 2D atomic arrangement
        - Useful for structure visualization
        - Reveals symmetry patterns

---

## 20. `plot_points_proj`

### Signature
```python
def plot_points_proj(points, direction, ax=None, **kwargs)
```

### Description
Project and plot arbitrary 3D points.
    General function to project any set of 3D points onto 2D plane.

### Input Parameters
points (array N×3): 3D point coordinates
        direction (array [3]): Projection direction
        ax (Axes): 2D axis (creates if None)
        **kwargs: Scatter plot parameters

### Output
fig, ax: Figure and axis

### Usage Examples
```python
>>> points = np.random.rand(100, 3)
        >>> fig, ax = plot_points_proj(points, [0,0,1])
        >>> ax.set_title('Random points - XY projection')
```

### Notes
- Generic projection function
        - Works with any point set
        - Not specific to lattices

---

## 21. `plot_lattice_plane`

### Signature
```python
def plot_lattice_plane(ax, L, h, k, l, **kwargs)
```

### Description
Plot a crystallographic plane in 3D lattice.
    Draws plane (hkl) intersecting lattice unit cell using matplotlib 3D.

### Input Parameters
ax (Axes3D): Matplotlib 3D axis
        L (array 3×3): Lattice matrix
        h, k, l (int): Miller indices of plane
        **kwargs: Additional matplotlib plot parameters (color, alpha, etc.)

### Output
None (modifies axis)

### Usage Examples
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
```

### Notes
- Plane intersects unit cell
        - Uses Miller indices for specification
        - Transparency recommended (alpha < 1)
        - Multiple planes can be overlaid

---

## 22. `plot_lattice_boundaries`

### Signature
```python
def plot_lattice_boundaries(ax, L, n1=1, n2=1, n3=1, **kwargs)
```

### Description
Plot unit cell boundaries in 3D.
    Draws edges of unit cells to visualize lattice structure.

### Input Parameters
ax (Axes3D): Matplotlib 3D axis
        L (array 3×3): Lattice matrix
        n1, n2, n3 (int): Number of cells in each direction
        **kwargs: Line plot parameters

### Output
None (modifies axis)

### Usage Examples
```python
>>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L = cubic_lattice_vec(3.0)
        >>> plot_lattice_boundaries(ax, L, n1=2, n2=2, n3=2, color='black')
        >>> 
        >>> set_aspect_equal_3d(ax)
        >>> plt.show()
```

### Notes
- Draws wireframe of unit cells
        - Helps visualize lattice structure
        - Combine with plot_lattice_points for complete view

---

## 23. `plot_atomic_plane2D`

### Signature
```python
def plot_atomic_plane2D(atoms_on_plane, plane_normal, **kwargs)
```

### Description
Plot 2D view of atomic plane.
    Creates 2D plot of atoms lying on crystallographic plane.

### Input Parameters
atoms_on_plane (array N×3): Coordinates of atoms on plane
        plane_normal (array [3]): Plane normal (for orientation)
        **kwargs: Scatter plot parameters

### Output
fig, ax: Figure and 2D axis

### Usage Examples
```python
>>> L = cubic_lattice_vec(3.0)
        >>> atoms = generate_lattite_atom_positions(L, n1=5, n2=5, n3=5)
        >>> indices = select_plane(atoms, 1, 1, 1, L)
        >>> plane_atoms = atoms[indices]
        >>> 
        >>> fig, ax = plot_atomic_plane2D(plane_atoms, [1,1,1], s=100, c='blue')
        >>> ax.set_title('(111) Atomic Plane')
        >>> plt.show()
```

### Notes
- 2D projection of 3D plane
        - Shows atomic arrangement
        - Useful for interface analysis

---

## 24. `plot_atomic_plane3D`

### Signature
```python
def plot_atomic_plane3D(atoms_on_plane, plane_normal, L, ax=None, **kwargs)
```

### Description
Plot atomic plane in 3D context.
    Visualizes atoms on plane with lattice context in 3D.

### Input Parameters
atoms_on_plane (array N×3): Atoms on plane
        plane_normal (array [3]): Plane normal
        L (array 3×3): Lattice matrix (for context)
        ax (Axes3D): 3D axis (creates if None)
        **kwargs: Scatter parameters

### Output
fig, ax: Figure and 3D axis

### Usage Examples
```python
>>> L = cubic_lattice_vec(3.0)
        >>> atoms_all = generate_lattite_atom_positions(L, n1=5, n2=5, n3=5)
        >>> indices = select_plane(atoms_all, 1, 0, 0, L)
        >>> plane_atoms = atoms_all[indices]
        >>> 
        >>> fig, ax = plot_atomic_plane3D(plane_atoms, [1,0,0], L, c='red', s=100)
        >>> # Also plot all atoms faintly
        >>> ax.scatter(atoms_all[:,0], atoms_all[:,1], atoms_all[:,2], 
        ...            c='gray', s=10, alpha=0.3)
        >>> set_aspect_equal_3d(ax)
```

### Notes
- Shows plane in 3D lattice context
        - Highlights selected atoms
        - Can overlay multiple planes

---

## 25. `plot_atomlattice2D`

### Signature
```python
def plot_atomlattice2D(L, basis, n1, n2, direction=[0,0,1], **kwargs)
```

### Description
Plot 2D atomic lattice projection.
    Complete 2D visualization of atomic structure with basis.

### Input Parameters
L (array 3×3): Lattice matrix
        basis (list): Atomic basis
        n1, n2 (int): Cell range
        direction (array [3]): Viewing direction
        **kwargs: Plot styling

### Output
fig, ax: Figure and axis

### Usage Examples
```python
>>> # FCC structure
        >>> L = cubic_lattice_vec(3.615)  # Al
        >>> basis = [[0,0,0], [0.5,0.5,0], [0.5,0,0.5], [0,0.5,0.5]]
        >>> 
        >>> fig, ax = plot_atomlattice2D(L, basis, n1=4, n2=4, 
        ...                               direction=[1,1,1], s=100)
        >>> ax.set_title('FCC [111] projection')
        >>> plt.show()
```

### Notes
- Shows 2D atomic arrangement
        - Includes all basis atoms
        - Useful for structure determination

---

## 26. `plot_cut2D`

### Signature
```python
def plot_cut2D(points, cut_plane_normal, cut_plane_point, **kwargs)
```

### Description
Plot 2D cross-section of 3D points.
    Selects points near plane and plots 2D projection.

### Input Parameters
points (array N×3): 3D points
        cut_plane_normal (array [3]): Cutting plane normal
        cut_plane_point (array [3]): Point on plane
        **kwargs: Plot parameters

### Output
fig, ax: Figure and 2D axis

### Usage Examples
```python
>>> atoms = generate_lattite_atom_positions(L, n1=10, n2=10, n3=10)
        >>> 
        >>> # Cut at z=15
        >>> fig, ax = plot_cut2D(atoms, [0,0,1], [0,0,15])
        >>> ax.set_title('Cross-section at z=15')
        >>> plt.show()
```

### Notes
- Shows 2D slice of 3D structure
        - Useful for interface analysis
        - Can reveal internal structure

---

## 27. `plot_planes_on_stereotriangle`

### Signature
```python
def plot_planes_on_stereotriangle(planes, **kwargs)
```

### Description
Plot planes on stereographic triangle.
    Projects plane normals onto standard stereographic triangle
    for cubic systems.

### Input Parameters
planes (list): List of (h,k,l) tuples
        **kwargs: Plot parameters

### Output
None (creates figure with stereographic triangle)

### Usage Examples
```python
>>> planes = [(1,0,0), (1,1,0), (1,1,1), (2,1,0)]
        >>> plot_planes_on_stereotriangle(planes)
        >>> plt.title('Low-Index Planes')
        >>> plt.show()
```

### Notes
- Standard triangle for cubic
        - Shows plane distribution
        - Used in texture analysis

---

## 28. `plot_planes_on_wulffnet`

### Signature
```python
def plot_planes_on_wulffnet(planes, lattice, **kwargs)
```

### Description
Plot plane normals on Wulff net (equal-angle projection).
    Projects planes onto full Wulff stereographic net.

### Input Parameters
planes (list): Plane Miller indices
        lattice (array 3×3): Lattice matrix
        **kwargs: Plot parameters

### Output
None (creates Wulff net figure)

### Usage Examples
```python
>>> L = cubic_lattice_vec(3.0)
        >>> planes = [(1,0,0), (0,1,0), (0,0,1), (1,1,1)]
        >>> plot_planes_on_wulffnet(planes, L)
        >>> plt.title('Planes on Wulff Net')
```

### Notes
- Equal-angle projection
        - Preserves angular relationships
        - Full sphere projection

---

## 29. `plot_princip_dir_on_stereotriangle`

### Signature
```python
def plot_princip_dir_on_stereotriangle(principal_directions, **kwargs)
```

### Description
Plot principal strain/stress directions on stereographic triangle.
    Projects principal directions onto standard triangle.

### Input Parameters
principal_directions (array 3×3): Principal direction matrix
        **kwargs: Plot styling

### Output
None (adds to current figure or creates new)

### Usage Examples
```python
>>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> 
        >>> plot_planes_on_stereotriangle([(1,0,0), (1,1,0), (1,1,1)])
        >>> plot_princip_dir_on_stereotriangle(result['directions'], 
        ...                                      marker='*', s=200, c='red')
        >>> plt.title('Principal Directions')
```

### Notes
- Shows orientation of principal axes
        - Overlays on crystal planes
        - Important for anisotropy analysis

---

## 30. `plot_princip_dir_on_wulffnet`

### Signature
```python
def plot_princip_dir_on_wulffnet(principal_directions, **kwargs)
```

### Description
Plot principal directions on Wulff net.
    Projects principal strain/stress axes onto Wulff stereographic net.

### Input Parameters
principal_directions (array 3×3): Principal directions
        **kwargs: Plot parameters

### Output
None (creates or modifies figure)

### Usage Examples
```python
>>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> plot_princip_dir_on_wulffnet(result['directions'])
        >>> plt.title('Principal Strain Axes')
```

### Notes
- Full sphere projection
        - Shows 3D orientation clearly

---


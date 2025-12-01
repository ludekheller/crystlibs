# PLOTLIB - DOCUMENTATION SUMMARY
## Organized Function Reference

**Module**: `plotlib.py` (Extended Version)  
**Total Functions**: 30  
**Version**: Extended 2024  

---

## 📊 FUNCTION OVERVIEW

### By Alphabetical Order

| Function | Purpose |
|----------|---------|
| `get_cmap` | Create a custom colormap from a list of colors with smooth gradients.... |
| `get_colors` | Map data values to colors using a colormap.... |
| `mohr_circles` | Calculate Mohr's circles from strain or stress tensor.... |
| `plot_atomic_plane2D` | Plot 2D view of atomic plane.... |
| `plot_atomic_plane3D` | Plot atomic plane in 3D context.... |
| `plot_atomlattice2D` | Plot 2D atomic lattice projection.... |
| `plot_cut2D` | Plot 2D cross-section of 3D points.... |
| `plot_lattice` | Complete lattice plotting function with atoms and unit cells.... |
| `plot_lattice2D` | Plot 2D projection of lattice.... |
| `plot_lattice3D` | Complete 3D lattice visualization.... |
| `plot_lattice_2Dprojection` | Project lattice along specific crystallographic direction.... |
| `plot_lattice_boundaries` | Plot unit cell boundaries in 3D.... |
| `plot_lattice_plane` | Plot a crystallographic plane in 3D lattice.... |
| `plot_lattice_proj` | Plot lattice projection along specific direction.... |
| `plot_latticefaces3D` | Plot lattice unit cell faces with transparency.... |
| `plot_latticesfaces3D` | Plot two lattices with transparent faces.... |
| `plot_mohr_circles` | Plot all three Mohr's circles.... |
| `plot_planes_on_mohr_circle` | Plot specific crystallographic planes on Mohr's circle.... |
| `plot_planes_on_stereotriangle` | Plot planes on stereographic triangle.... |
| `plot_planes_on_wulffnet` | Plot plane normals on Wulff net (equal-angle projection).... |
| `plot_points_proj` | Project and plot arbitrary 3D points.... |
| `plot_princip_dir_on_stereotriangle` | Plot principal strain/stress directions on stereographic triangle.... |
| `plot_princip_dir_on_wulffnet` | Plot principal directions on Wulff net.... |
| `plotcolmap` | Plot a single colormap from a saved pickle file.... |
| `plotcolmaps` | Plot multiple colormaps from a saved pickle file.... |
| `set_aspect_equal_3d` | Fix equal aspect ratio bug for 3D matplotlib plots.... |
| `shiftedColorMap` | Shift the center of a colormap to a specific data value.... |
| `strains_along_13mohrcirle` | Calculate strains along maximum Mohr's circle.... |
| `write_mohr_planes` | Write Mohr circle analysis results to file.... |
| `zero_normal_strains` | Find planes with zero normal strain in given strain state.... |

---

## 📖 DETAILED FUNCTION DESCRIPTIONS

### `get_cmap`

**Signature**: `def get_cmap(colors, nbins=1000, name='my_cmap')`

Create a custom colormap from a list of colors with smooth gradients.

---

### `get_colors`

**Signature**: `def get_colors(values, cmap, vmin=None, vmax=None, cmin=0, cmax=None, 
               to255=True, nancolor=[0, 0, 0, 1])`

Map data values to colors using a colormap.

---

### `mohr_circles`

**Signature**: `def mohr_circles(strain_tensor)`

Calculate Mohr's circles from strain or stress tensor.
    Computes principal strains/stresses and parameters for three Mohr's
    circles. Used for visualizing 3D strain state in 2D.

---

### `plot_atomic_plane2D`

**Signature**: `def plot_atomic_plane2D(atoms_on_plane, plane_normal, **kwargs)`

Plot 2D view of atomic plane.
    Creates 2D plot of atoms lying on crystallographic plane.

---

### `plot_atomic_plane3D`

**Signature**: `def plot_atomic_plane3D(atoms_on_plane, plane_normal, L, ax=None, **kwargs)`

Plot atomic plane in 3D context.
    Visualizes atoms on plane with lattice context in 3D.

---

### `plot_atomlattice2D`

**Signature**: `def plot_atomlattice2D(L, basis, n1, n2, direction=[0,0,1], **kwargs)`

Plot 2D atomic lattice projection.
    Complete 2D visualization of atomic structure with basis.

---

### `plot_cut2D`

**Signature**: `def plot_cut2D(points, cut_plane_normal, cut_plane_point, **kwargs)`

Plot 2D cross-section of 3D points.
    Selects points near plane and plots 2D projection.

---

### `plot_lattice`

**Signature**: `def plot_lattice(L, n1=1, n2=1, n3=1, basis=None, ax=None, **kwargs)`

Complete lattice plotting function with atoms and unit cells.
    High-level function combining atomic positions, cell boundaries,
    and styling in single call.

---

### `plot_lattice2D`

**Signature**: `def plot_lattice2D(ax, L, n1=2, n2=2, projection_axis=2, **kwargs)`

Plot 2D projection of lattice.
    Projects 3D lattice onto 2D plane for simplified visualization.

---

### `plot_lattice3D`

**Signature**: `def plot_lattice3D(ax, L, n1=1, n2=1, n3=1, show_points=True, show_edges=True, **kwargs)`

Complete 3D lattice visualization.
    Plots lattice points and/or edges in single function call.

---

### `plot_lattice_2Dprojection`

**Signature**: `def plot_lattice_2Dprojection(L, direction=[0,0,1], **kwargs)`

Project lattice along specific crystallographic direction.
    Creates 2D projection viewed along given direction vector.

---

### `plot_lattice_boundaries`

**Signature**: `def plot_lattice_boundaries(ax, L, n1=1, n2=1, n3=1, **kwargs)`

Plot unit cell boundaries in 3D.
    Draws edges of unit cells to visualize lattice structure.

---

### `plot_lattice_plane`

**Signature**: `def plot_lattice_plane(ax, L, h, k, l, **kwargs)`

Plot a crystallographic plane in 3D lattice.
    Draws plane (hkl) intersecting lattice unit cell using matplotlib 3D.

---

### `plot_lattice_proj`

**Signature**: `def plot_lattice_proj(L, direction, n1=2, n2=2, basis=None, **kwargs)`

Plot lattice projection along specific direction.
    Creates 2D projection of 3D lattice viewed along given direction.

---

### `plot_latticefaces3D`

**Signature**: `def plot_latticefaces3D(ax, L, alpha=0.3, **kwargs)`

Plot lattice unit cell faces with transparency.
    Visualizes unit cell as transparent parallelepiped.

---

### `plot_latticesfaces3D`

**Signature**: `def plot_latticesfaces3D(ax, L1, L2, alpha1=0.2, alpha2=0.2, **kwargs)`

Plot two lattices with transparent faces.
    Visualizes two crystal structures simultaneously for comparison.

---

### `plot_mohr_circles`

**Signature**: `def plot_mohr_circles(mohr_result, **kwargs)`

Plot all three Mohr's circles.
    Visualizes complete 3D strain state via Mohr's circles.

---

### `plot_planes_on_mohr_circle`

**Signature**: `def plot_planes_on_mohr_circle(mohr_result, plane_normals, lattice, **kwargs)`

Plot specific crystallographic planes on Mohr's circle.
    Shows where particular crystal planes plot on Mohr diagram.

---

### `plot_planes_on_stereotriangle`

**Signature**: `def plot_planes_on_stereotriangle(planes, **kwargs)`

Plot planes on stereographic triangle.
    Projects plane normals onto standard stereographic triangle
    for cubic systems.

---

### `plot_planes_on_wulffnet`

**Signature**: `def plot_planes_on_wulffnet(planes, lattice, **kwargs)`

Plot plane normals on Wulff net (equal-angle projection).
    Projects planes onto full Wulff stereographic net.

---

### `plot_points_proj`

**Signature**: `def plot_points_proj(points, direction, ax=None, **kwargs)`

Project and plot arbitrary 3D points.
    General function to project any set of 3D points onto 2D plane.

---

### `plot_princip_dir_on_stereotriangle`

**Signature**: `def plot_princip_dir_on_stereotriangle(principal_directions, **kwargs)`

Plot principal strain/stress directions on stereographic triangle.
    Projects principal directions onto standard triangle.

---

### `plot_princip_dir_on_wulffnet`

**Signature**: `def plot_princip_dir_on_wulffnet(principal_directions, **kwargs)`

Plot principal directions on Wulff net.
    Projects principal strain/stress axes onto Wulff stereographic net.

---

### `plotcolmap`

**Signature**: `def plotcolmap(fname=None, withdraw=False)`

Plot a single colormap from a saved pickle file.
    Loads and executes plot configuration from pickle file.

---

### `plotcolmaps`

**Signature**: `def plotcolmaps(fname=None, withdraw=False)`

Plot multiple colormaps from a saved pickle file.
    Creates a 2x2 subplot figure with multiple pole figures/colormaps.
    Loads data from pickle file containing plot configurations.

---

### `set_aspect_equal_3d`

**Signature**: `def set_aspect_equal_3d(ax)`

Fix equal aspect ratio bug for 3D matplotlib plots.
    Adjusts axis limits to ensure equal scaling in all three dimensions,
    preventing distortion in 3D visualizations of crystal structures.
    Essential for accurate representation of lattice geometries.

---

### `shiftedColorMap`

**Signature**: `def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap')`

Shift the center of a colormap to a specific data value.
    Useful for data with asymmetric ranges (e.g., -15 to +5) where you want
    the colormap center (typically white or neutral color) at zero.

---

### `strains_along_13mohrcirle`

**Signature**: `def strains_along_13mohrcirle(mohr_result, n_points=100)`

Calculate strains along maximum Mohr's circle.
    Computes normal and shear strain components for points on the
    largest Mohr's circle (ε₁-ε₃ circle).

---

### `write_mohr_planes`

**Signature**: `def write_mohr_planes(filename, mohr_result, planes, lattice)`

Write Mohr circle analysis results to file.
    Saves principal strains, circles, and plane-specific strains.

---

### `zero_normal_strains`

**Signature**: `def zero_normal_strains(strain_tensor)`

Find planes with zero normal strain in given strain state.
    Calculates planes where normal strain εₙ = nᵀ·ε·n = 0.
    These planes experience only shear.

---


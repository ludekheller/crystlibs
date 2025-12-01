# PLOTLIB - QUICK REFERENCE (Comprehensive)
## Fast Function Lookup

**Module**: `plotlib.py` (Extended Version)  
**Total Functions**: 30  

---

## 🔍 QUICK LOOKUP TABLE

| # | Function | Parameters | Returns |
|---|----------|------------|---------|
| 1 | `get_cmap` | `colors, nbins=1000, name='my_c...` | cmap: matplotlib.colors.LinearSegmentedC |
| 2 | `get_colors` | `values, cmap, vmin=None, vmax=...` | Colors: numpy array - Array of colors (R |
| 3 | `mohr_circles` | `strain_tensor` | dict: { |
| 4 | `plot_atomic_plane2D` | `atoms_on_plane, plane_normal, ...` | fig, ax: Figure and 2D axis |
| 5 | `plot_atomic_plane3D` | `atoms_on_plane, plane_normal, ...` | fig, ax: Figure and 3D axis |
| 6 | `plot_atomlattice2D` | `L, basis, n1, n2, direction=[0...` | fig, ax: Figure and axis |
| 7 | `plot_cut2D` | `points, cut_plane_normal, cut_...` | fig, ax: Figure and 2D axis |
| 8 | `plot_lattice` | `L, n1=1, n2=1, n3=1, basis=Non...` | fig, ax: Matplotlib figure and axis obje |
| 9 | `plot_lattice2D` | `ax, L, n1=2, n2=2, projection_...` | None (modifies axis) |
| 10 | `plot_lattice3D` | `ax, L, n1=1, n2=1, n3=1, show_...` | None (modifies axis) |
| 11 | `plot_lattice_2Dprojection` | `L, direction=[0,0,1], **kwargs` | None (creates new figure) |
| 12 | `plot_lattice_boundaries` | `ax, L, n1=1, n2=1, n3=1, **kwa...` | None (modifies axis) |
| 13 | `plot_lattice_plane` | `ax, L, h, k, l, **kwargs` | None (modifies axis) |
| 14 | `plot_lattice_proj` | `L, direction, n1=2, n2=2, basi...` | fig, ax: Figure and axis |
| 15 | `plot_latticefaces3D` | `ax, L, alpha=0.3, **kwargs` | None (modifies axis) |
| 16 | `plot_latticesfaces3D` | `ax, L1, L2, alpha1=0.2, alpha2...` | None (modifies axis) |
| 17 | `plot_mohr_circles` | `mohr_result, **kwargs` | None (creates matplotlib figure) |
| 18 | `plot_planes_on_mohr_circle` | `mohr_result, plane_normals, la...` | None (adds to current figure) |
| 19 | `plot_planes_on_stereotriangle` | `planes, **kwargs` | None (creates figure with stereographic  |
| 20 | `plot_planes_on_wulffnet` | `planes, lattice, **kwargs` | None (creates Wulff net figure) |
| 21 | `plot_points_proj` | `points, direction, ax=None, **...` | fig, ax: Figure and axis |
| 22 | `plot_princip_dir_on_stereotriangle` | `principal_directions, **kwargs` | None (adds to current figure or creates  |
| 23 | `plot_princip_dir_on_wulffnet` | `principal_directions, **kwargs` | None (creates or modifies figure) |
| 24 | `plotcolmap` | `fname=None, withdraw=False` | None (creates matplotlib figure) |
| 25 | `plotcolmaps` | `fname=None, withdraw=False` | None (creates matplotlib figure) |
| 26 | `set_aspect_equal_3d` | `ax` | None (modifies axis object in place) |
| 27 | `shiftedColorMap` | `cmap, start=0, midpoint=0.5, s...` | newcmap: matplotlib.colors.LinearSegment |
| 28 | `strains_along_13mohrcirle` | `mohr_result, n_points=100` | dict: { |
| 29 | `write_mohr_planes` | `filename, mohr_result, planes,...` | None (writes to file) |
| 30 | `zero_normal_strains` | `strain_tensor` | list: List of plane normals with zero no |

---

## 📝 FUNCTION SUMMARIES

### `get_cmap`
Create a custom colormap from a list of colors with smooth gradients.

### `get_colors`
Map data values to colors using a colormap.

### `mohr_circles`
Calculate Mohr's circles from strain or stress tensor.

### `plot_atomic_plane2D`
Plot 2D view of atomic plane.

### `plot_atomic_plane3D`
Plot atomic plane in 3D context.

### `plot_atomlattice2D`
Plot 2D atomic lattice projection.

### `plot_cut2D`
Plot 2D cross-section of 3D points.

### `plot_lattice`
Complete lattice plotting function with atoms and unit cells.

### `plot_lattice2D`
Plot 2D projection of lattice.

### `plot_lattice3D`
Complete 3D lattice visualization.

### `plot_lattice_2Dprojection`
Project lattice along specific crystallographic direction.

### `plot_lattice_boundaries`
Plot unit cell boundaries in 3D.

### `plot_lattice_plane`
Plot a crystallographic plane in 3D lattice.

### `plot_lattice_proj`
Plot lattice projection along specific direction.

### `plot_latticefaces3D`
Plot lattice unit cell faces with transparency.

### `plot_latticesfaces3D`
Plot two lattices with transparent faces.

### `plot_mohr_circles`
Plot all three Mohr's circles.

### `plot_planes_on_mohr_circle`
Plot specific crystallographic planes on Mohr's circle.

### `plot_planes_on_stereotriangle`
Plot planes on stereographic triangle.

### `plot_planes_on_wulffnet`
Plot plane normals on Wulff net (equal-angle projection).

### `plot_points_proj`
Project and plot arbitrary 3D points.

### `plot_princip_dir_on_stereotriangle`
Plot principal strain/stress directions on stereographic triangle.

### `plot_princip_dir_on_wulffnet`
Plot principal directions on Wulff net.

### `plotcolmap`
Plot a single colormap from a saved pickle file.

### `plotcolmaps`
Plot multiple colormaps from a saved pickle file.

### `set_aspect_equal_3d`
Fix equal aspect ratio bug for 3D matplotlib plots.

### `shiftedColorMap`
Shift the center of a colormap to a specific data value.

### `strains_along_13mohrcirle`
Calculate strains along maximum Mohr's circle.

### `write_mohr_planes`
Write Mohr circle analysis results to file.

### `zero_normal_strains`
Find planes with zero normal strain in given strain state.


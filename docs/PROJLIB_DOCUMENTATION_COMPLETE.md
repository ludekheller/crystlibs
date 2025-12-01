# PROJLIB - COMPLETE DOCUMENTATION
## Stereographic Projections and Coordinate Transformations Library

**Module**: `projlib.py` (Extended Version)  
**Total Functions**: 46  
**Version**: Extended 2024  
**Status**: Production Ready  

---

## 📚 TABLE OF CONTENTS

1. [genoritri](#genoritri)
2. [genori](#genori)
3. [xyz2spher](#xyz2spher)
4. [spher2xyz](#spher2xyz)
5. [genprojgrid](#genprojgrid)
6. [gen_dirs_norms](#gen_dirs_norms)
7. [fullcirc_hist](#fullcirc_hist)
8. [rp2xyz](#rp2xyz)
9. [stereoprojection_directions](#stereoprojection_directions)
10. [equalarea_planes](#equalarea_planes)
11. [stereo2xyz](#stereo2xyz)
12. [equalarea2xyz](#equalarea2xyz)
13. [equalarea_directions](#equalarea_directions)
14. [wulffnet](#wulffnet)
15. [wulffnet_half](#wulffnet_half)
16. [wulffnet_quarter](#wulffnet_quarter)
17. [schmidtnet](#schmidtnet)
18. [wulffnet_regular_grid](#wulffnet_regular_grid)
19. [schmidtnet_half](#schmidtnet_half)
20. [wulffnet_regular_grid](#wulffnet_regular_grid)
21. [schmidt_regular_grid](#schmidt_regular_grid)
22. [pf_cmap02](#pf_cmap02)
23. [pf](#pf)
24. [pf_cmap](#pf_cmap)
25. [pf_cmap_cscale](#pf_cmap_cscale)
26. [ipf](#ipf)
27. [stereotriangle](#stereotriangle)
28. [colored_stereotriangle](#colored_stereotriangle)
29. [filled_colored_stereotriangle](#filled_colored_stereotriangle)
30. [colors4stereotriangle](#colors4stereotriangle)
31. [stereotriangle_colors_from_d_IPF](#stereotriangle_colors_from_d_ipf)
32. [stereotriangle_colors_from_eumats_dir](#stereotriangle_colors_from_eumats_dir)
33. [stereotriangle_colors](#stereotriangle_colors)
34. [equivalent_elements](#equivalent_elements)
35. [stereoprojection_intotriangle_ini](#stereoprojection_intotriangle_ini)
36. [stereoprojection_intotriangle_fast](#stereoprojection_intotriangle_fast)
37. [stereoprojection_intotriangle](#stereoprojection_intotriangle)
38. [equalarea_intotriangle_fast](#equalarea_intotriangle_fast)
39. [equalarea_intotriangle](#equalarea_intotriangle)
40. [stereoprojection_planes](#stereoprojection_planes)
41. [iszero](#iszero)
42. [gcd](#gcd)
43. [gcdarr](#gcdarr)
44. [vector2miller](#vector2miller)
45. [vector2millerround](#vector2millerround)
46. [vectors2miller](#vectors2miller)

---

## 1. `genoritri`

### Signature
```python
def genoritri(resolution=1.0, mesh="spherified_cube_edge")
```

### Description
Generate orientation grid for stereographic triangle using orix library.
    Creates a uniform grid of crystallographic directions suitable for plotting
    in the standard stereographic triangle. Uses diffsims beam direction generator
    with orix quaternion/vector operations.

### Input Parameters
resolution: float - Angular resolution in degrees (default: 1.0)
                           Smaller values create finer grids but more points
        mesh: str - Mesh generation method (default: "spherified_cube_edge")
                   Options:
                   - "spherified_cube_edge": Uniform angular spacing
                   - "spherified_cube_corner": Corner-focused spacing
                   - "cubochoric": Volume-preserving grid

### Output
trioris: numpy array (3, N) - Direction vectors in Cartesian coordinates [x, y, z]
                                       Each column is a unit vector

### Usage Examples
```python
>>> # Fine grid for smooth contours
        >>> oris_fine = genoritri(resolution=0.5)
        >>> print("Grid shape:", oris_fine.shape)
        >>> # Typical output: (3, ~50000) for resolution=0.5
        >>> # Coarse grid for quick visualization
        >>> oris_coarse = genoritri(resolution=2.0)
        >>> print("Number of orientations:", oris_coarse.shape[1])
        >>> # Typical output: ~3000 points
        >>> # Use in pole figure
        >>> from plotlib import plotter
        >>> p = plotter()
        >>> p.plotProj(ProjType='equalarea', sphere='triangle')
        >>> proj_coords = equalarea_directions(oris_fine)
        >>> p.ax.scatter(proj_coords[0], proj_coords[1], s=1, c='red')
        >>> # Different mesh types
        >>> oris_corner = genoritri(resolution=1.0, mesh="spherified_cube_corner")
        >>> oris_cubo = genoritri(resolution=1.0, mesh="cubochoric")
```

---

## 2. `genori`

### Signature
```python
def genori(dangle=1.0,hemi='both', tol=1e-2, rot=np.eye(3), half='no')
```

### Description
Generate a set of orientations by rotating around two perpendicular axes.
    Creates a grid of directions by rotating around Z-axis (φ₁: 0-360°) and 
    Y-axis (φ₂: 0-180°), then optionally applies rotation matrix and filters
    by hemisphere. Useful for generating uniform orientation distributions.

### Input Parameters
dangle: float - Angular resolution in degrees (default: 1.0)
                       Number of orientations = (360/dangle) × (180/dangle)
        hemi: str - Hemisphere selection in real space (default: 'both')
                   Options:
                   - 'both': All orientations (full sphere)
                   - 'upper': z > -tol (upper hemisphere)
                   - 'lower': z < tol (lower hemisphere)
        tol: float - Tolerance for hemisphere filtering (default: 1e-2)
        rot: numpy array (3, 3) - Rotation matrix to apply to all orientations
                                   (default: identity matrix)
        half: str - Additional filter in projection space (default: 'no')
                   Options:
                   - 'no': No additional filtering
                   - 'upper': Project and keep only upper projection hemisphere
                   - 'lower': Project and keep only lower projection hemisphere

### Output
oris: numpy array (3, N) - Direction vectors [x, y, z]
                                    Each column is a unit vector
                                    N depends on dangle and filtering options

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # Generate orientations with 5° resolution
        >>> oris = genori(dangle=5.0)
        >>> print("Shape:", oris.shape)
        >>> # Output: (3, 12960) = (3, 72×180)
        >>> print("Number of orientations:", oris.shape[1])
        >>> # Upper hemisphere only (e.g., for pole figures)
        >>> oris_upper = genori(dangle=2.0, hemi='upper')
        >>> print("All z > -tol:", np.all(oris_upper[2,:] > -1e-2))
        >>> # Output: True
        >>> # Apply 90° rotation around Z-axis
        >>> R_z90 = np.array([[0, -1, 0],
        ...                   [1,  0, 0],
        ...                   [0,  0, 1]])
        >>> oris_rot = genori(dangle=10.0, rot=R_z90)
        >>> 
        >>> # Fine grid for texture analysis
        >>> oris_texture = genori(dangle=1.0, hemi='upper', half='upper')
        >>> print("Fine grid size:", oris_texture.shape[1])
        >>> # Coarse grid for quick preview
        >>> oris_preview = genori(dangle=10.0, hemi='upper')
        >>> print("Preview grid size:", oris_preview.shape[1])
        >>> # Use with pole figure plotting
        >>> import matplotlib.pyplot as plt
        >>> proj = equalarea_directions(oris_upper)
        >>> plt.figure(figsize=(6,6))
        >>> plt.scatter(proj[0], proj[1], s=1, alpha=0.5)
        >>> plt.axis('equal')
        >>> plt.show()
```

---

## 3. `xyz2spher`

### Signature
```python
def xyz2spher(xyz,deg=False)
```

### Description
Convert Cartesian coordinates to spherical coordinates (polar and azimuthal angles).
    Transforms 3D Cartesian vectors to spherical coordinate system using standard
    physics convention: θ (polar/colatitude angle from +Z axis), φ (azimuthal angle in XY plane).

### Input Parameters
xyz: numpy array (N, 3) - Cartesian coordinates [x, y, z]
                                   Each row is one point
        deg: bool - If True, return angles in degrees; if False, radians (default: False)

### Output
polar_angle: numpy array (N,) - Polar angle θ (colatitude from z-axis)
                                        Range: [0, π] radians or [0, 180°]
        azimuth_angle: numpy array (N,) - Azimuthal angle φ (longitude in xy-plane)
                                          Range: [-π, π] radians or [-180°, 180°]

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # Convert single point
        >>> xyz = np.array([[1, 0, 0]])  # Point on +X axis
        >>> theta, phi = xyz2spher(xyz, deg=True)
        >>> print(f"θ={theta[0]:.1f}°, φ={phi[0]:.1f}°")
        >>> # Output: θ=90.0°, φ=0.0° (equator, 0° longitude)
        >>> # Convert multiple points
        >>> xyz_points = np.array([
        ...     [1, 0, 0],      # +X axis
        ...     [0, 1, 0],      # +Y axis
        ...     [0, 0, 1],      # +Z axis (north pole)
        ...     [1, 1, 1]       # Diagonal
        ... ])
        >>> theta, phi = xyz2spher(xyz_points, deg=True)
        >>> for i, (t, p) in enumerate(zip(theta, phi)):
        ...     print(f"Point {i}: θ={t:.1f}°, φ={p:.1f}°")
        >>> # Output:
        >>> # Point 0: θ=90.0°, φ=0.0°
        >>> # Point 1: θ=90.0°, φ=90.0°
        >>> # Point 2: θ=0.0°, φ=0.0°
        >>> # Point 3: θ=54.7°, φ=45.0°
        >>> # Verify round-trip conversion
        >>> xyz_back = spher2xyz(theta, phi, deg=True)
        >>> print("Close match:", np.allclose(xyz_points, xyz_back))
        >>> # Output: True
        >>> # Use for texture analysis - convert orientations to angles
        >>> oris = genori(dangle=5.0, hemi='upper')
        >>> theta_all, phi_all = xyz2spher(oris.T, deg=True)
        >>> print("Polar angle range:", theta_all.min(), "to", theta_all.max())
        >>> print("Azimuthal angle range:", phi_all.min(), "to", phi_all.max())
```

---

## 4. `spher2xyz`

### Signature
```python
def spher2xyz(polar_angle,azimuth_angle,deg=False)
```

### Description
Convert spherical coordinates to Cartesian coordinates.
    Inverse transformation of xyz2spher. Converts spherical coordinates (θ, φ) to
    Cartesian (x, y, z) using standard physics convention.

### Input Parameters
polar_angle: float or array - Polar angle θ (colatitude from z-axis)
                                      Expected range: [0, π] or [0, 180°]
        azimuth_angle: float or array - Azimuthal angle φ (longitude in xy-plane)
                                        Expected range: [-π, π] or [-180°, 180°]
        deg: bool - If True, input angles are in degrees; if False, radians (default: False)

### Output
xyz: numpy array (N, 3) - Cartesian coordinates [x, y, z]
                                   Each row is one point (unit vector if inputs were from unit sphere)

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # Convert single point (θ=90°, φ=45°)
        >>> theta = 90.0  # On equator
        >>> phi = 45.0    # 45° from +X axis in XY plane
        >>> xyz = spher2xyz(theta, phi, deg=True)
        >>> print("Cartesian:", xyz[0])
        >>> # Output: [0.707, 0.707, 0.0] (≈[√2/2, √2/2, 0])
        >>> # Convert multiple points
        >>> theta_arr = np.array([0, 90, 180])  # North pole, equator, south pole
        >>> phi_arr = np.array([0, 0, 0])
        >>> xyz_arr = spher2xyz(theta_arr, phi_arr, deg=True)
        >>> print("Shape:", xyz_arr.shape)  # (3, 3)
        >>> print("North pole:", xyz_arr[0])  # [0, 0, 1]
        >>> print("Equator:", xyz_arr[1])     # [1, 0, 0]
        >>> print("South pole:", xyz_arr[2])  # [0, 0, -1]
        >>> # Create points on equator circle
        >>> theta_eq = np.ones(8) * 90  # All at θ=90° (equator)
        >>> phi_eq = np.linspace(0, 315, 8)  # Around z-axis
        >>> xyz_eq = spher2xyz(theta_eq, phi_eq, deg=True)
        >>> print("Points on equator:", xyz_eq.shape)  # (8, 3)
        >>> # Verify round-trip with xyz2spher
        >>> theta_back, phi_back = xyz2spher(xyz_eq, deg=True)
        >>> print("Round-trip close:", np.allclose(theta_eq, theta_back))
        >>> # Output: True
        >>> # Generate uniform sphere sampling
        >>> n_theta = 18  # 10° resolution in polar
        >>> n_phi = 36    # 10° resolution in azimuthal
        >>> theta_grid = np.linspace(0, 180, n_theta)
        >>> phi_grid = np.linspace(-180, 180, n_phi)
        >>> theta_mesh, phi_mesh = np.meshgrid(theta_grid, phi_grid)
        >>> xyz_sphere = spher2xyz(theta_mesh.flatten(), phi_mesh.flatten(), deg=True)
        >>> print("Sphere points:", xyz_sphere.shape)
```

---

## 5. `genprojgrid`

### Signature
```python
def genprojgrid(oris,gdata=None,nump=1001,proj='equalarea',method2='linear',gdout=False,poris=None,minmax='notfull')
```

### Description
Generate a regular 2D grid for projected orientations with optional density data.
    Creates a square grid covering the projection space and interpolates orientation
    data onto this grid. Essential for creating smooth contour plots, density maps,
    and pole figures.

### Input Parameters
oris: numpy array (3, N) - Orientation vectors in Cartesian coordinates [x, y, z]
        gdata: numpy array (N,) - Data values for each orientation (default: None)
                                   If None, returns only grid structure without data
        nump: int - Number of grid points in each direction (default: 1001)
                   Higher values give smoother contours but slower computation
        proj: str - Projection type (default: 'equalarea')
                   Options: 'equalarea' (Schmidt net), 'stereo' (Wulff net)
        method2: str - Interpolation method for griddata (default: 'linear')
                      Options: 'linear', 'nearest', 'cubic'
        gdout: bool - Return gridded data (default: False)
                     Set True when gdata is provided
        poris: numpy array (2, N) - Pre-computed projected orientations (optional)
                                     If None, computed from oris using proj type
        minmax: str - Grid extent mode (default: 'notfull')
                     'full': Use maximum theoretical projection extent
                     'notfull': Use actual data extent (tighter bounds)

### Output
If gdout=False or gdata=None:
            grid_x: numpy array (nump, nump) - X-coordinates of grid points
            grid_y: numpy array (nump, nump) - Y-coordinates of grid points
            nummask: numpy array (nump, nump) - Numerical mask (1=valid region, 0=outside)
            mask: numpy array (nump, nump) - Boolean mask (True=masked/invalid, False=valid)
        If gdout=True and gdata provided:
            grid_x: numpy array (nump, nump) - X-coordinates
            grid_y: numpy array (nump, nump) - Y-coordinates
            grid_z: masked array (nump, nump) - Interpolated data values on grid
            nummask: numpy array (nump, nump) - Numerical mask
            mask: numpy array (nump, nump) - Boolean mask

### Usage Examples
```python
>>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> 
        >>> # Generate orientation grid
        >>> oris = genori(dangle=5.0, hemi='upper')
        >>> print("Number of orientations:", oris.shape[1])
        >>> 
        >>> # Create projection grid structure only (no data)
        >>> grid_x, grid_y, nummask, mask = genprojgrid(
        ...     oris,
        ...     nump=501,
        ...     proj='equalarea'
        ... )
        >>> print("Grid shape:", grid_x.shape)
        >>> 
        >>> # With density data - create random texture
        >>> density = np.random.rand(oris.shape[1]) * 100
        >>> grid_x, grid_y, grid_z, nummask, mask = genprojgrid(
        ...     oris,
        ...     gdata=density,
        ...     nump=301,
        ...     proj='equalarea',
        ...     method2='linear',
        ...     gdout=True
        ... )
        >>> 
        >>> # Plot contour map
        >>> plt.figure(figsize=(8, 8))
        >>> levels = np.linspace(0, 100, 21)
        >>> plt.contourf(grid_x, grid_y, grid_z, levels=levels, cmap='jet')
        >>> plt.colorbar(label='Density')
        >>> plt.axis('equal')
        >>> plt.title('Pole Figure Density Map')
        >>> plt.show()
        >>> 
        >>> # High-resolution grid for publication
        >>> grid_x_hr, grid_y_hr, grid_z_hr, _, _ = genprojgrid(
        ...     oris,
        ...     gdata=density,
        ...     nump=2001,  # Very fine grid
        ...     proj='equalarea',
        ...     method2='cubic',  # Smooth interpolation
        ...     minmax='full'     # Full projection extent
        ... )
        >>> 
        >>> # Stereographic projection instead
        >>> grid_x_st, grid_y_st, grid_z_st, _, _ = genprojgrid(
        ...     oris,
        ...     gdata=density,
        ...     proj='stereo',
        ...     gdout=True
        ... )
        >>> 
        >>> # Use with pre-computed projections (faster)
        >>> poris_precomp = equalarea_directions(oris)
        >>> grid_x, grid_y, grid_z, _, _ = genprojgrid(
        ...     oris,
        ...     gdata=density,
        ...     poris=poris_precomp,  # Skip projection step
        ...     gdout=True
        ... )
```

---

## 6. `gen_dirs_norms`

### Signature
```python
def gen_dirs_norms(L, Lr, uvws,hkls, R2Proj=np.eye(3),symops=None,recsymops=None,hemisphere = "upper", **kwargs)
```

### Description
Generate all symmetry-equivalent crystal directions and plane normals with projections.
    Creates comprehensive lists of crystallographic directions [uvw] and normals (hkl) including
    all symmetry equivalents. For each direction/normal, computes:
    - 3D vector in Cartesian coordinates
    - Equal-area projection coordinates
    - Stereographic projection coordinates
    - Plane traces (for normals)

### Input Parameters
L: numpy array (3, 3) - Direct lattice basis matrix
                               Columns are lattice vectors a, b, c
                               [uvw] direction = L @ [u, v, w]
        Lr: numpy array (3, 3) - Reciprocal lattice basis matrix
                                (hkl) normal = Lr @ [h, k, l]
        uvws: list of lists - Crystal directions [[u1,v1,w1], [u2,v2,w2], ...]
                             Examples: [[1,0,0], [1,1,0], [1,1,1]]
        hkls: list of lists - Plane normals (Miller indices) [[h1,k1,l1], [h2,k2,l2], ...]
                             Examples: [[1,0,0], [1,1,0], [1,1,1]]
        R2Proj: numpy array (3, 3) - Additional rotation to projection frame (default: identity)
                                      Applied before projection
        symops: list of numpy arrays - Direct space symmetry operations (default: None)
                                       Each element is a 3×3 symmetry matrix
                                       If None, only input directions are used
        recsymops: list of numpy arrays - Reciprocal space symmetry operations (default: None)
                                          For cubic: same as symops
                                          For non-cubic: different symmetry
        hemisphere: str - Hemisphere for plane traces (default: "upper")
                         Options: "upper", "lower", "both", "triangle"

### Output
dirs: list of dicts - Crystal directions, each with keys:
            'vector': numpy array (3,) - Direction vector in Cartesian coordinates
            'uvw': list [u,v,w] - Miller indices
            'uvwf': list - Family representative (original input)
            'label': str - Text label for plotting (e.g., "[110]")
            'equalarea': numpy array (2,) - Equal-area projection coordinates
            'stereo': numpy array (2,) - Stereographic projection coordinates
            'textshift': list [dx, dy] - Label position offset for plotting
        normals: list of dicts - Plane normals, each with keys:
            'vector': numpy array (3,) - Normal vector in Cartesian coordinates
            'hkl': list [h,k,l] - Miller indices
            'hklf': list - Family representative (original input)
            'label': str - Text label for plotting (e.g., "(110)")
            'equalarea': numpy array (2,) - Normal direction projection
            'stereo': numpy array (2,) - Normal direction projection
            'equalarea plane': numpy array (2, N) - Plane trace coordinates (equal-area)
            'stereo plane': numpy array (2, N) - Plane trace coordinates (stereographic)
            'textshift': list [dx, dy] - Label position offset

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # Cubic crystal (e.g., Fe, Al, Cu)
        >>> a = 3.5  # Lattice parameter in Angstrom
        >>> L = a * np.eye(3)  # Direct lattice basis
        >>> Lr = (2*np.pi/a) * np.eye(3)  # Reciprocal lattice basis
        >>> 
        >>> # Define important crystallographic directions and planes
        >>> uvws = [[1,0,0], [1,1,0], [1,1,1]]  # <100>, <110>, <111> families
        >>> hkls = [[1,0,0], [1,1,0], [1,1,1]]  # {100}, {110}, {111} families
        >>> 
        >>> # Generate without symmetry (only input directions)
        >>> dirs, normals = gen_dirs_norms(L, Lr, uvws, hkls)
        >>> print(f"Directions: {len(dirs)}")  # Output: 3
        >>> print(f"Normals: {len(normals)}")   # Output: 3
        >>> 
        >>> # With cubic symmetry (48 symmetry operations)
        >>> from orilib import symmetry_elements
        >>> symops = symmetry_elements('cubic')
        >>> dirs_sym, normals_sym = gen_dirs_norms(
        ...     L, Lr, uvws, hkls,
        ...     symops=symops,
        ...     recsymops=symops  # For cubic: direct = reciprocal symmetry
        ... )
        >>> print(f"With symmetry - Directions: {len(dirs_sym)}")
        >>> # Output: 26 (6×<100> + 12×<110> + 8×<111>)
        >>> print(f"With symmetry - Normals: {len(normals_sym)}")
        >>> # Output: 26 (6×{100} + 12×{110} + 8×{111})
        >>> 
        >>> # Access projection data
        >>> dir0 = dirs[0]
        >>> print(f"Direction: {dir0['uvw']}")
        >>> print(f"Vector: {dir0['vector']}")
        >>> print(f"Equal-area projection: {dir0['equalarea']}")
        >>> print(f"Label: {dir0['label']}")
        >>> 
        >>> # Plot on pole figure
        >>> import matplotlib.pyplot as plt
        >>> fig, ax = plt.subplots(figsize=(8, 8))
        >>> 
        >>> # Draw projection circle
        >>> circle = plt.Circle((0, 0), np.sqrt(2), fill=False, color='k')
        >>> ax.add_patch(circle)
        >>> 
        >>> # Plot plane traces
        >>> for normal in normals_sym:
        ...     plane_trace = normal['equalarea plane']
        ...     ax.plot(plane_trace[0,:], plane_trace[1,:], 'k-', linewidth=0.5)
        >>> 
        >>> # Plot direction points
        >>> for dir in dirs_sym:
        ...     proj = dir['equalarea']
        ...     ax.plot(proj[0], proj[1], 'ro', markersize=8)
        ...     ax.text(proj[0]+0.05, proj[1]+0.05, dir['label'], fontsize=8)
        >>> 
        >>> ax.set_aspect('equal')
        >>> ax.set_xlim(-1.5, 1.5)
        >>> ax.set_ylim(-1.5, 1.5)
        >>> plt.title('Pole Figure with Cubic Symmetry')
        >>> plt.show()
        >>> 
        >>> # Hexagonal crystal example
        >>> a_hex = 3.2
        >>> c_hex = 5.2
        >>> L_hex = np.array([[a_hex, -a_hex/2, 0],
        ...                   [0, a_hex*np.sqrt(3)/2, 0],
        ...                   [0, 0, c_hex]])
        >>> Lr_hex = 2*np.pi * np.linalg.inv(L_hex).T
        >>> 
        >>> uvws_hex = [[1,0,0], [1,1,0], [0,0,1]]
        >>> hkls_hex = [[1,0,0], [1,1,0], [0,0,1]]
        >>> 
        >>> symops_hex = symmetry_elements('hexagonal')
        >>> dirs_hex, normals_hex = gen_dirs_norms(
        ...     L_hex, Lr_hex, uvws_hex, hkls_hex,
        ...     symops=symops_hex,
        ...     recsymops=symops_hex
        ... )
        >>> print(f"Hexagonal directions: {len(dirs_hex)}")
```

---

## 7. `fullcirc_hist`

### Signature
```python
def fullcirc_hist(Mats, Dr=[0,0,1], symops=None, equalarea=False, scale='sqrt', nlevels=10, lvls=None,bins=128, ax=None, title=None, ret=False, 
    """
    Create a pole figure with density contours from orientation matrices.
    
    Generates inverse pole figures showing the distribution of a specific crystal
    direction in the sample reference frame. Supports stereographic and equal-area
    projections with optional kernel density estimation.
    
    Input:
        Mats: numpy array (N, 3, 3) - Orientation matrices (crystal→sample)
        Dr: list/array (3,) - Reference direction in sample frame (default: [0,0,1])
        symops: numpy array (Ns, 3, 3) - Crystal symmetry operations (default: None)
        equalarea: bool - Use equal-area if True, stereographic if False
        bins: int - Histogram bins for density (default: 128)
        kernel: bool - Use spherical KDE instead of histogram (default: False)
        bandwidth: float - KDE bandwidth (default: None, auto)
        **kwargs: Additional arguments
    
    Output:
        fig, ax: matplotlib figure and axis objects
    
    Usage Example:
        >>> from orilib import np_eulers_matrices
        >>> euler = np.random.rand(1000, 3) * [360, 180, 360]
        >>> Mats = np_eulers_matrices(euler, deg=True)
        >>> fig, ax = fullcirc_hist(Mats, Dr=[0,0,1], equalarea=True)
    """
                  kernel=False,  weights=None,Lim=None,interp=True,interpn=1000, smooth=False, vmin=None, 
                  bandwidth=None, vmax=None,colorbar=True,ticks=None, R2Proj=None,contour=True, mrd=False, **kwargs)
```

### Description
No docstring available.

---

## 8. `rp2xyz`

### Signature
```python
def rp2xyz(r,p)
```

### Description
Convert polar coordinates (r, phi) to 3D Cartesian coordinates.

### Input Parameters
r: float/array - Polar angle in degrees (0-180)
        p: float/array - Azimuthal angle in degrees (0-360)

### Output
x, y, z: tuple - Cartesian coordinates

---

## 9. `stereoprojection_directions`

### Signature
```python
def stereoprojection_directions(dirs)
```

### Description
Project 3D direction vectors onto 2D stereographic projection plane (Wulff net).
    Uses stereographic projection from south pole: r = tan(θ/2), where θ is the
    angle from the north pole. Points on lower hemisphere project outside unit circle.

### Input Parameters
dirs: numpy array (3, N) or (3,) - Direction vectors [x, y, z]
                                           Will be automatically normalized
                                           Single vector (3,) will be reshaped to (3, 1)

### Output
proj_dirs: numpy array (3, N) - Projected coordinates [X, Y, 0]
                                         X, Y are 2D projection coordinates
                                         Third row is zeros (kept for compatibility)

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # Project single direction
        >>> dir_vec = np.array([1, 0, 0])  # +X axis
        >>> proj = stereoprojection_directions(dir_vec)
        >>> print("Projected:", proj[:2, 0])  # [X, Y]
        >>> 
        >>> # Project multiple directions
        >>> dirs = np.array([
        ...     [1, 0, 0],      # +X axis (on equator)
        ...     [0, 1, 0],      # +Y axis (on equator)
        ...     [0, 0, 1],      # +Z axis (north pole)
        ...     [0, 0, -1],     # -Z axis (south pole)
        ...     [1, 1, 1]       # Upper hemisphere
        ... ]).T
        >>> proj_all = stereoprojection_directions(dirs)
        >>> print("Shape:", proj_all.shape)  # (3, 5)
        >>> print("North pole projects to:", proj_all[:2, 2])  # [0, 0]
        >>> 
        >>> # Visualize projection
        >>> import matplotlib.pyplot as plt
        >>> fig, ax = plt.subplots(figsize=(8, 8))
        >>> 
        >>> # Draw unit circle
        >>> circle = plt.Circle((0, 0), 1, fill=False, color='k')
        >>> ax.add_patch(circle)
        >>> 
        >>> # Project grid of directions
        >>> oris = genori(dangle=10.0, hemi='upper')
        >>> proj_grid = stereoprojection_directions(oris)
        >>> ax.scatter(proj_grid[0], proj_grid[1], s=1, alpha=0.5)
        >>> 
        >>> ax.set_aspect('equal')
        >>> ax.set_xlim(-1.2, 1.2)
        >>> ax.set_ylim(-1.2, 1.2)
        >>> plt.title('Stereographic Projection (Wulff Net)')
        >>> plt.show()
        >>> 
        >>> # Compare with equal-area
        >>> proj_stereo = stereoprojection_directions(dirs)
        >>> proj_ea = equalarea_directions(dirs)
        >>> print("Stereographic radii:", np.linalg.norm(proj_stereo[:2], axis=0))
        >>> print("Equal-area radii:", np.linalg.norm(proj_ea[:2], axis=0))
```

---

## 10. `equalarea_planes`

### Signature
```python
def equalarea_planes(normals,arclength=360.,iniangle=0.,hemisphere="both")
```

### Description
Project plane traces onto equal-area (Schmidt) net.

### Input Parameters
normals: numpy array (3, N) - Plane normal vectors
        arclength: float - Arc length in degrees (default: 360)
        iniangle: float - Starting angle in degrees (default: 0)
        hemisphere: str - 'both', 'upper', or 'lower' (default: 'both')

### Output
traces: list of arrays - Trace points for each plane

---

## 11. `stereo2xyz`

### Signature
```python
def stereo2xyz(projdir)
```

### Description
Convert 2D stereographic projection to 3D unit vectors.
    Inverse of stereoprojection_directions().

### Input Parameters
projdir: numpy array (2, N) - Projection coordinates [X, Y]

### Output
dirs: numpy array (3, N) - Unit direction vectors

---

## 12. `equalarea2xyz`

### Signature
```python
def equalarea2xyz(projdir)
```

### Description
Convert 2D equal-area projection to 3D unit vectors.
    Inverse of equalarea_directions().

### Input Parameters
projdir: numpy array (2, N) - Projection coordinates

### Output
dirs: numpy array (3, N) - Unit direction vectors

---

## 13. `equalarea_directions`

### Signature
```python
def equalarea_directions(dirs)
```

### Description
Project 3D direction vectors onto 2D equal-area projection plane (Schmidt net).
    Uses Lambert azimuthal equal-area projection: r = √2 × sin(θ/2), where θ is the
    angle from the north pole. Preserves area, making it ideal for texture analysis.
    Maximum radius is √2 for hemisphere projection.

### Input Parameters
dirs: numpy array (3, N) or (3,) - Direction vectors [x, y, z]
                                           Will be automatically normalized
                                           Single vector (3,) will be reshaped to (3, 1)

### Output
proj_dirs: numpy array (3, N) - Projected coordinates [X, Y, 0]
                                         X, Y are 2D equal-area projection coordinates
                                         Third row is zeros (kept for compatibility)
                                         Range: [-√2, √2] for full sphere

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # Project single direction
        >>> dir_vec = np.array([1, 0, 0])  # +X axis
        >>> proj = equalarea_directions(dir_vec)
        >>> print("Projected:", proj[:2, 0])
        >>> print("Radius:", np.linalg.norm(proj[:2, 0]))  # Should be 1.0 (on equator)
        >>> 
        >>> # Project multiple directions
        >>> dirs = np.array([
        ...     [1, 0, 0],      # +X axis (equator)
        ...     [0, 1, 0],      # +Y axis (equator)
        ...     [0, 0, 1],      # +Z axis (north pole)
        ...     [1, 1, 1]       # Upper hemisphere
        ... ]).T
        >>> proj_all = equalarea_directions(dirs)
        >>> print("Shape:", proj_all.shape)  # (3, 4)
        >>> print("North pole:", proj_all[:2, 2])  # [0, 0]
        >>> 
        >>> # Create pole figure
        >>> import matplotlib.pyplot as plt
        >>> oris = genori(dangle=5.0, hemi='upper')
        >>> proj = equalarea_directions(oris)
        >>> 
        >>> fig, ax = plt.subplots(figsize=(8, 8))
        >>> circle = plt.Circle((0, 0), np.sqrt(2), fill=False, color='k', linewidth=2)
        >>> ax.add_patch(circle)
        >>> ax.scatter(proj[0], proj[1], s=1, alpha=0.5, c='blue')
        >>> ax.set_aspect('equal')
        >>> ax.set_xlim(-1.5, 1.5)
        >>> ax.set_ylim(-1.5, 1.5)
        >>> plt.title('Equal-Area Projection (Schmidt Net)')
        >>> plt.show()
```

---

## 14. `wulffnet`

### Signature
```python
def wulffnet(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.))
```

### Description
Draw stereographic (Wulff) net - full circle.

### Input Parameters
ax: matplotlib axis - Existing axis (default: None)
        basedirs: bool - Plot base directions (default: False)
        facecolor: tuple - Background color RGB

### Output
fig, ax: matplotlib figure and axis

---

## 15. `wulffnet_half`

### Signature
```python
def wulffnet_half(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.))
```

### Description
Draw stereographic (Wulff) net - upper hemisphere.

### Input Parameters
ax: matplotlib axis - Existing axis (default: None)
        basedirs: bool - Plot base directions (default: False)
        facecolor: tuple - Background color RGB

### Output
fig, ax: matplotlib figure and axis

---

## 16. `wulffnet_quarter`

### Signature
```python
def wulffnet_quarter(ax=None,basedirs=False)
```

### Description
Draw quarter stereographic (Wulff) net.

### Input Parameters
ax: matplotlib axis - Existing axis (default: None)
        basedirs: bool - Plot base directions (default: False)

### Output
fig, ax: matplotlib figure and axis

---

## 17. `schmidtnet`

### Signature
```python
def schmidtnet(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.))
```

### Description
Draw equal-area (Schmidt) net - full circle.

### Input Parameters
ax: matplotlib axis - Existing axis (default: None)
        basedirs: bool - Plot base directions (default: False)
        facecolor: tuple - Background color RGB

### Output
fig, ax: matplotlib figure and axis

---

## 18. `wulffnet_regular_grid`

### Signature
```python
def wulffnet_regular_grid(ax,dangle,dirout=False, plot=True)
```

### Description
Create regular angular grid on Wulff net.

### Input Parameters
ax: matplotlib axis
        dangle: float - Angular spacing in degrees
        dirout: bool - Return directions if True (default: False)
        plot: bool - Plot grid if True (default: True)

### Output
If dirout=True: numpy array (3, N) - Grid directions

---

## 19. `schmidtnet_half`

### Signature
```python
def schmidtnet_half(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.))
```

### Description
Draw equal-area (Schmidt) net - upper hemisphere.

### Input Parameters
ax: matplotlib axis - Existing axis (default: None)
        basedirs: bool - Plot base directions (default: False)
        facecolor: tuple - Background color RGB

### Output
fig, ax: matplotlib figure and axis

---

## 20. `wulffnet_regular_grid`

### Signature
```python
def wulffnet_regular_grid(ax,dangle)
```

### Description
Create regular angular grid on Wulff net.

### Input Parameters
ax: matplotlib axis
        dangle: float - Angular spacing in degrees
        dirout: bool - Return directions if True (default: False)
        plot: bool - Plot grid if True (default: True)

### Output
If dirout=True: numpy array (3, N) - Grid directions

---

## 21. `schmidt_regular_grid`

### Signature
```python
def schmidt_regular_grid(ax,Na=72,Nr=20,plot=True)
```

### Description
Create regular grid on Schmidt net.

### Input Parameters
ax: matplotlib axis
        Na: int - Azimuthal divisions (default: 72)
        Nr: int - Radial divisions (default: 20)
        plot: bool - Plot grid (default: True)

### Output
grid: tuple - Grid parameters

---

## 22. `pf_cmap02`

### Signature
```python
def pf_cmap02(GridX,GridY,GridR,GridPhi,AreaRatio,Intensity,NoCont=10,GridSize=1000,cmap='jet',method='cubic')
```

### Description
Create contour pole figure with colormap (internal function).

### Input Parameters
GridX, GridY, GridR, GridPhi, AreaRatio, Intensity: Grid and density data
        NoCont, GridSize, cmap, method: Plotting parameters

### Output
Contour plot object

---

## 23. `pf`

### Signature
```python
def pf(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20,syms=True,s=50,facecolor=(210./255.,235./255.,255./255.),plot=True)
```

### Description
Generate pole figure from Euler angles.

### Input Parameters
gPhi1, gPHI, gPhi2: arrays - Bunge Euler angles (degrees)
        Dc: array (3,) - Crystal direction [uvw]
        lattice: str - Crystal system
        **kwargs: Additional parameters

### Output
fig, ax: matplotlib figure and axis

---

## 24. `pf_cmap`

### Signature
```python
def pf_cmap(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20,syms=True,s=50,NoCont=10,GridSize=1000,cmap='jet',method='cubic')
```

### Description
Generate pole figure with density colormap.

### Input Parameters
gPhi1, gPHI, gPhi2: arrays - Bunge Euler angles (degrees)
        Dc: array (3,) - Crystal direction [uvw]
        lattice: str - Crystal system
        **kwargs: Additional parameters

### Output
fig, ax: matplotlib figure and axis

---

## 25. `pf_cmap_cscale`

### Signature
```python
def pf_cmap_cscale(fig,ax2,cmin,cmax,cmap)
```

### Description
Add colorbar to pole figure.

### Input Parameters
fig, ax2, cmin, cmax, cmap: Figure and colorbar parameters

### Output
None (modifies figure)

---

## 26. `ipf`

### Signature
```python
def ipf(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20,syms=True)
```

### Description
Generate inverse pole figure from Euler angles.

### Input Parameters
gPhi1, gPHI, gPhi2: arrays - Bunge Euler angles (degrees)
        Dc: array (3,) - Sample direction [xyz]
        lattice: str - Crystal system
        **kwargs: Additional parameters

### Output
fig, ax: matplotlib figure and axis

---

## 27. `stereotriangle`

### Signature
```python
def stereotriangle(ax=None,basedirs=False,equalarea=False,grid=False,resolution=None,gridmarkersize=None,gridmarkercol=None,gridzorder=None,mesh=False)
```

### Description
Draw standard stereographic triangle for cubic system.

### Input Parameters
ax: matplotlib axis - Existing axis (default: None)
        basedirs: bool - Label corners (default: False)
        equalarea: bool - Use equal-area (default: False)
        grid: bool - Draw grid (default: False)
        **kwargs: Additional parameters

### Output
fig, ax: matplotlib figure and axis

---

## 28. `colored_stereotriangle`

### Signature
```python
def colored_stereotriangle(basedirs=False,resolution = 1, markersize=1)
```

### Description
Draw stereographic triangle with IPF coloring.

### Input Parameters
basedirs: bool - Label corners (default: False)
        resolution: float - Grid resolution (default: 1)
        markersize: float - Point size (default: 1)

### Output
fig, ax: matplotlib figure and axis

---

## 29. `filled_colored_stereotriangle`

### Signature
```python
def filled_colored_stereotriangle(basedirs=False,resolution = 1, markersize=1,ax=None,**kwargs)
```

### Description
Draw filled triangle with smooth IPF coloring.

### Input Parameters
basedirs, resolution, markersize, ax, **kwargs

### Output
fig, ax: matplotlib figure and axis

---

## 30. `colors4stereotriangle`

### Signature
```python
def colors4stereotriangle(resolution = 1)
```

### Description
Generate RGB colors for IPF coloring.

### Input Parameters
resolution: float - Angular resolution (default: 1)

### Output
proj, RGB: Projection coordinates and RGB colors

---

## 31. `stereotriangle_colors_from_d_IPF`

### Signature
```python
def stereotriangle_colors_from_d_IPF(d_IPF)
```

### Description
Get IPF colors for crystal directions.

### Input Parameters
d_IPF: numpy array (3, N) - Directions in crystal frame

### Output
RGB: numpy array (N, 3) - RGB colors

---

## 32. `stereotriangle_colors_from_eumats_dir`

### Signature
```python
def stereotriangle_colors_from_eumats_dir(eumats,d=[1,0,0])
```

### Description
Get IPF colors from orientation matrices.

### Input Parameters
eumats: numpy array (N, 3, 3) - Orientation matrices
        d: list/array (3,) - Sample direction (default: [1,0,0])

### Output
RGB: numpy array (N, 3) - RGB colors

---

## 33. `stereotriangle_colors`

### Signature
```python
def stereotriangle_colors(proj_Ds)
```

### Description
Get IPF colors from projection coordinates.

### Input Parameters
proj_Ds: numpy array (2, N) - Projection coordinates

### Output
RGB: numpy array (N, 3) - RGB colors

---

## 34. `equivalent_elements`

### Signature
```python
def equivalent_elements(element,lattice)
```

### Description
Find symmetrically equivalent elements.

### Input Parameters
element: numpy array (3,) - Direction/normal vector
        lattice: str - Crystal system

### Output
equivalents: list of arrays - Equivalent vectors

---

## 35. `stereoprojection_intotriangle_ini`

### Signature
```python
def stereoprojection_intotriangle_ini(dirs,eps=1.0e-5)
```

### Description
Map directions into standard triangle (initial version).

### Input Parameters
dirs: numpy array (3, N) - Direction vectors
        eps: float - Tolerance (default: 1e-5)

### Output
proj: numpy array (2, N) - Projection coordinates

---

## 36. `stereoprojection_intotriangle_fast`

### Signature
```python
def stereoprojection_intotriangle_fast(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False,Rin=None,symops=None)
```

### Description
Fast mapping of directions into standard triangle.

### Input Parameters
dirs: numpy array (3, N) - Direction vectors
        eps, geteqdirs, geteqmats, Rin, symops: Optional parameters

### Output
proj: numpy array (2, N) - Projection coordinates

---

## 37. `stereoprojection_intotriangle`

### Signature
```python
def stereoprojection_intotriangle(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False,Rin=None,symops=None)
```

### Description
Map directions into standard stereographic triangle.

### Input Parameters
dirs: numpy array (3, N) - Direction vectors
        eps, geteqdirs, geteqmats, Rin, symops: Optional parameters

### Output
proj: numpy array (2, N) - Projection coordinates

---

## 38. `equalarea_intotriangle_fast`

### Signature
```python
def equalarea_intotriangle_fast(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False,Rin=None,symops=None)
```

### Description
Fast equal-area projection into triangle.

### Input Parameters
dirs: numpy array (3, N) - Direction vectors
        eps, geteqdirs, geteqmats, Rin, symops: Optional parameters

### Output
proj: numpy array (2, N) - Projection coordinates

---

## 39. `equalarea_intotriangle`

### Signature
```python
def equalarea_intotriangle(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False)
```

### Description
Equal-area projection into standard triangle.

### Input Parameters
dirs: numpy array (3, N) - Direction vectors
        eps, geteqdirs, geteqmats: Optional parameters

### Output
proj: numpy array (2, N) - Projection coordinates

---

## 40. `stereoprojection_planes`

### Signature
```python
def stereoprojection_planes(normals,arclength=360.,iniangle=0.,hemisphere='both',getpoints=False,R=np.eye(3))
```

### Description
Project plane traces onto stereographic projection.

### Input Parameters
normals: numpy array (3, N) - Plane normals
        arclength, iniangle, hemisphere, getpoints, R: Optional parameters

### Output
traces: list of arrays or plots

---

## 41. `iszero`

### Signature
```python
def iszero(a)
```

### Description
Check if value is approximately zero (|a| < 1e-9).

---

## 42. `gcd`

### Signature
```python
def gcd(a,b)
```

### Description
Compute greatest common divisor of two numbers.

---

## 43. `gcdarr`

### Signature
```python
def gcdarr(arr)
```

### Description
Compute GCD of all array elements.

---

## 44. `vector2miller`

### Signature
```python
def vector2miller(arr, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3)
```

### Description
Convert vector to Miller indices using GCD reduction.

### Input Parameters
arr: array (3,) - Vector [x, y, z]
        MIN, Tol, tol, text, decimals: Optional parameters

### Output
miller: array (3,) - Miller indices [h, k, l]

---

## 45. `vector2millerround`

### Signature
```python
def vector2millerround(v, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3)
```

### Description
Convert vector to Miller indices with rounding.

### Input Parameters
v: array (3,) - Vector [x, y, z]
        MIN, Tol, tol, text, decimals: Optional parameters

### Output
miller: array (3,) - Miller indices

---

## 46. `vectors2miller`

### Signature
```python
def vectors2miller(V, MIN=True, Tol=1e-9,tol=1e5,text=False)
```

### Description
Convert multiple vectors to Miller indices.

### Input Parameters
V: numpy array (3, N) - Vectors
        MIN, Tol, tol, text: Optional parameters

### Output
VM: numpy array (3, N) - Miller indices

---


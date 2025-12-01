# PROJLIB - QUICK REFERENCE (Comprehensive)
## Fast Function Lookup

**Module**: `projlib.py` (Extended Version)  
**Total Functions**: 46  

---

## 🔍 QUICK LOOKUP TABLE

| # | Function | Parameters | Returns |
|---|----------|------------|---------|
| 1 | `colored_stereotriangle` | `basedirs=False,resolution = 1,...` | fig, ax: matplotlib figure and axis |
| 2 | `colors4stereotriangle` | `resolution = 1` | proj, RGB: Projection coordinates and RG |
| 3 | `equalarea2xyz` | `projdir` | dirs: numpy array (3, N) - Unit directio |
| 4 | `equalarea_directions` | `dirs` | proj_dirs: numpy array (3, N) - Projecte |
| 5 | `equalarea_intotriangle` | `dirs,eps=1.0e-5,geteqdirs=Fals...` | proj: numpy array (2, N) - Projection co |
| 6 | `equalarea_intotriangle_fast` | `dirs,eps=1.0e-5,geteqdirs=Fals...` | proj: numpy array (2, N) - Projection co |
| 7 | `equalarea_planes` | `normals,arclength=360.,iniangl...` | traces: list of arrays - Trace points fo |
| 8 | `equivalent_elements` | `element,lattice` | equivalents: list of arrays - Equivalent |
| 9 | `filled_colored_stereotriangle` | `basedirs=False,resolution = 1,...` | fig, ax: matplotlib figure and axis |
| 10 | `fullcirc_hist` | `Mats, Dr=[0,0,1], symops=None,...` | See docs |
| 11 | `gcd` | `a,b` | See docs |
| 12 | `gcdarr` | `arr` | See docs |
| 13 | `gen_dirs_norms` | `L, Lr, uvws,hkls, R2Proj=np.ey...` | dirs: list of dicts - Crystal directions |
| 14 | `genori` | `dangle=1.0,hemi='both', tol=1e...` | oris: numpy array (3, N) - Direction vec |
| 15 | `genoritri` | `resolution=1.0, mesh="spherifi...` | trioris: numpy array (3, N) - Direction  |
| 16 | `genprojgrid` | `oris,gdata=None,nump=1001,proj...` | If gdout=False or gdata=None: |
| 17 | `ipf` | `gPhi1,gPHI,gPhi2,Dc,lattice,Na...` | fig, ax: matplotlib figure and axis |
| 18 | `iszero` | `a` | See docs |
| 19 | `pf` | `gPhi1,gPHI,gPhi2,Dc,lattice,Na...` | fig, ax: matplotlib figure and axis |
| 20 | `pf_cmap` | `gPhi1,gPHI,gPhi2,Dc,lattice,Na...` | fig, ax: matplotlib figure and axis |
| 21 | `pf_cmap02` | `GridX,GridY,GridR,GridPhi,Area...` | Contour plot object |
| 22 | `pf_cmap_cscale` | `fig,ax2,cmin,cmax,cmap` | None (modifies figure) |
| 23 | `rp2xyz` | `r,p` | x, y, z: tuple - Cartesian coordinates |
| 24 | `schmidt_regular_grid` | `ax,Na=72,Nr=20,plot=True` | grid: tuple - Grid parameters |
| 25 | `schmidtnet` | `ax=None,basedirs=False,facecol...` | fig, ax: matplotlib figure and axis |
| 26 | `schmidtnet_half` | `ax=None,basedirs=False,facecol...` | fig, ax: matplotlib figure and axis |
| 27 | `spher2xyz` | `polar_angle,azimuth_angle,deg=...` | xyz: numpy array (N, 3) - Cartesian coor |
| 28 | `stereo2xyz` | `projdir` | dirs: numpy array (3, N) - Unit directio |
| 29 | `stereoprojection_directions` | `dirs` | proj_dirs: numpy array (3, N) - Projecte |
| 30 | `stereoprojection_intotriangle` | `dirs,eps=1.0e-5,geteqdirs=Fals...` | proj: numpy array (2, N) - Projection co |
| 31 | `stereoprojection_intotriangle_fast` | `dirs,eps=1.0e-5,geteqdirs=Fals...` | proj: numpy array (2, N) - Projection co |
| 32 | `stereoprojection_intotriangle_ini` | `dirs,eps=1.0e-5` | proj: numpy array (2, N) - Projection co |
| 33 | `stereoprojection_planes` | `normals,arclength=360.,iniangl...` | traces: list of arrays or plots |
| 34 | `stereotriangle` | `ax=None,basedirs=False,equalar...` | fig, ax: matplotlib figure and axis |
| 35 | `stereotriangle_colors` | `proj_Ds` | RGB: numpy array (N, 3) - RGB colors |
| 36 | `stereotriangle_colors_from_d_IPF` | `d_IPF` | RGB: numpy array (N, 3) - RGB colors |
| 37 | `stereotriangle_colors_from_eumats_dir` | `eumats,d=[1,0,0]` | RGB: numpy array (N, 3) - RGB colors |
| 38 | `vector2miller` | `arr, MIN=True, Tol=1e-9,tol=1e...` | miller: array (3,) - Miller indices [h,  |
| 39 | `vector2millerround` | `v, MIN=True, Tol=1e-9,tol=1e5,...` | miller: array (3,) - Miller indices |
| 40 | `vectors2miller` | `V, MIN=True, Tol=1e-9,tol=1e5,...` | VM: numpy array (3, N) - Miller indices |
| 41 | `wulffnet` | `ax=None,basedirs=False,facecol...` | fig, ax: matplotlib figure and axis |
| 42 | `wulffnet_half` | `ax=None,basedirs=False,facecol...` | fig, ax: matplotlib figure and axis |
| 43 | `wulffnet_quarter` | `ax=None,basedirs=False` | fig, ax: matplotlib figure and axis |
| 44 | `wulffnet_regular_grid` | `ax,dangle,dirout=False, plot=T...` | If dirout=True: numpy array (3, N) - Gri |
| 45 | `wulffnet_regular_grid` | `ax,dangle` | If dirout=True: numpy array (3, N) - Gri |
| 46 | `xyz2spher` | `xyz,deg=False` | polar_angle: numpy array (N,) - Polar an |

---

## 📝 FUNCTION SUMMARIES

### `colored_stereotriangle`
Draw stereographic triangle with IPF coloring.

### `colors4stereotriangle`
Generate RGB colors for IPF coloring.

### `equalarea2xyz`
Convert 2D equal-area projection to 3D unit vectors.

### `equalarea_directions`
Project 3D direction vectors onto 2D equal-area projection plane (Schmidt net).

### `equalarea_intotriangle`
Equal-area projection into standard triangle.

### `equalarea_intotriangle_fast`
Fast equal-area projection into triangle.

### `equalarea_planes`
Project plane traces onto equal-area (Schmidt) net.

### `equivalent_elements`
Find symmetrically equivalent elements.

### `filled_colored_stereotriangle`
Draw filled triangle with smooth IPF coloring.

### `fullcirc_hist`
No docstring available.

### `gcd`
Compute greatest common divisor of two numbers.

### `gcdarr`
Compute GCD of all array elements.

### `gen_dirs_norms`
Generate all symmetry-equivalent crystal directions and plane normals with projections.

### `genori`
Generate a set of orientations by rotating around two perpendicular axes.

### `genoritri`
Generate orientation grid for stereographic triangle using orix library.

### `genprojgrid`
Generate a regular 2D grid for projected orientations with optional density data.

### `ipf`
Generate inverse pole figure from Euler angles.

### `iszero`
Check if value is approximately zero (|a| < 1e-9).

### `pf`
Generate pole figure from Euler angles.

### `pf_cmap`
Generate pole figure with density colormap.

### `pf_cmap02`
Create contour pole figure with colormap (internal function).

### `pf_cmap_cscale`
Add colorbar to pole figure.

### `rp2xyz`
Convert polar coordinates (r, phi) to 3D Cartesian coordinates.

### `schmidt_regular_grid`
Create regular grid on Schmidt net.

### `schmidtnet`
Draw equal-area (Schmidt) net - full circle.

### `schmidtnet_half`
Draw equal-area (Schmidt) net - upper hemisphere.

### `spher2xyz`
Convert spherical coordinates to Cartesian coordinates.

### `stereo2xyz`
Convert 2D stereographic projection to 3D unit vectors.

### `stereoprojection_directions`
Project 3D direction vectors onto 2D stereographic projection plane (Wulff net).

### `stereoprojection_intotriangle`
Map directions into standard stereographic triangle.

### `stereoprojection_intotriangle_fast`
Fast mapping of directions into standard triangle.

### `stereoprojection_intotriangle_ini`
Map directions into standard triangle (initial version).

### `stereoprojection_planes`
Project plane traces onto stereographic projection.

### `stereotriangle`
Draw standard stereographic triangle for cubic system.

### `stereotriangle_colors`
Get IPF colors from projection coordinates.

### `stereotriangle_colors_from_d_IPF`
Get IPF colors for crystal directions.

### `stereotriangle_colors_from_eumats_dir`
Get IPF colors from orientation matrices.

### `vector2miller`
Convert vector to Miller indices using GCD reduction.

### `vector2millerround`
Convert vector to Miller indices with rounding.

### `vectors2miller`
Convert multiple vectors to Miller indices.

### `wulffnet`
Draw stereographic (Wulff) net - full circle.

### `wulffnet_half`
Draw stereographic (Wulff) net - upper hemisphere.

### `wulffnet_quarter`
Draw quarter stereographic (Wulff) net.

### `wulffnet_regular_grid`
Create regular angular grid on Wulff net.

### `wulffnet_regular_grid`
Create regular angular grid on Wulff net.

### `xyz2spher`
Convert Cartesian coordinates to spherical coordinates (polar and azimuthal angles).


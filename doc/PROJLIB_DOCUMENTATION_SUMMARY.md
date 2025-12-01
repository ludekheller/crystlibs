# PROJLIB - DOCUMENTATION SUMMARY
## Organized Function Reference

**Module**: `projlib.py` (Extended Version)  
**Total Functions**: 46  
**Version**: Extended 2024  

---

## 📊 FUNCTION OVERVIEW

### By Alphabetical Order

| Function | Purpose |
|----------|---------|
| `colored_stereotriangle` | Draw stereographic triangle with IPF coloring.... |
| `colors4stereotriangle` | Generate RGB colors for IPF coloring.... |
| `equalarea2xyz` | Convert 2D equal-area projection to 3D unit vectors.... |
| `equalarea_directions` | Project 3D direction vectors onto 2D equal-area projection plane (Schmidt net).... |
| `equalarea_intotriangle` | Equal-area projection into standard triangle.... |
| `equalarea_intotriangle_fast` | Fast equal-area projection into triangle.... |
| `equalarea_planes` | Project plane traces onto equal-area (Schmidt) net.... |
| `equivalent_elements` | Find symmetrically equivalent elements.... |
| `filled_colored_stereotriangle` | Draw filled triangle with smooth IPF coloring.... |
| `fullcirc_hist` | No docstring available.... |
| `gcd` | Compute greatest common divisor of two numbers.... |
| `gcdarr` | Compute GCD of all array elements.... |
| `gen_dirs_norms` | Generate all symmetry-equivalent crystal directions and plane normals with proje... |
| `genori` | Generate a set of orientations by rotating around two perpendicular axes.... |
| `genoritri` | Generate orientation grid for stereographic triangle using orix library.... |
| `genprojgrid` | Generate a regular 2D grid for projected orientations with optional density data... |
| `ipf` | Generate inverse pole figure from Euler angles.... |
| `iszero` | Check if value is approximately zero (|a| < 1e-9).... |
| `pf` | Generate pole figure from Euler angles.... |
| `pf_cmap` | Generate pole figure with density colormap.... |
| `pf_cmap02` | Create contour pole figure with colormap (internal function).... |
| `pf_cmap_cscale` | Add colorbar to pole figure.... |
| `rp2xyz` | Convert polar coordinates (r, phi) to 3D Cartesian coordinates.... |
| `schmidt_regular_grid` | Create regular grid on Schmidt net.... |
| `schmidtnet` | Draw equal-area (Schmidt) net - full circle.... |
| `schmidtnet_half` | Draw equal-area (Schmidt) net - upper hemisphere.... |
| `spher2xyz` | Convert spherical coordinates to Cartesian coordinates.... |
| `stereo2xyz` | Convert 2D stereographic projection to 3D unit vectors.... |
| `stereoprojection_directions` | Project 3D direction vectors onto 2D stereographic projection plane (Wulff net).... |
| `stereoprojection_intotriangle` | Map directions into standard stereographic triangle.... |
| `stereoprojection_intotriangle_fast` | Fast mapping of directions into standard triangle.... |
| `stereoprojection_intotriangle_ini` | Map directions into standard triangle (initial version).... |
| `stereoprojection_planes` | Project plane traces onto stereographic projection.... |
| `stereotriangle` | Draw standard stereographic triangle for cubic system.... |
| `stereotriangle_colors` | Get IPF colors from projection coordinates.... |
| `stereotriangle_colors_from_d_IPF` | Get IPF colors for crystal directions.... |
| `stereotriangle_colors_from_eumats_dir` | Get IPF colors from orientation matrices.... |
| `vector2miller` | Convert vector to Miller indices using GCD reduction.... |
| `vector2millerround` | Convert vector to Miller indices with rounding.... |
| `vectors2miller` | Convert multiple vectors to Miller indices.... |
| `wulffnet` | Draw stereographic (Wulff) net - full circle.... |
| `wulffnet_half` | Draw stereographic (Wulff) net - upper hemisphere.... |
| `wulffnet_quarter` | Draw quarter stereographic (Wulff) net.... |
| `wulffnet_regular_grid` | Create regular angular grid on Wulff net.... |
| `wulffnet_regular_grid` | Create regular angular grid on Wulff net.... |
| `xyz2spher` | Convert Cartesian coordinates to spherical coordinates (polar and azimuthal angl... |

---

## 📖 DETAILED FUNCTION DESCRIPTIONS

### `colored_stereotriangle`

**Signature**: `def colored_stereotriangle(basedirs=False,resolution = 1, markersize=1)`

Draw stereographic triangle with IPF coloring.

---

### `colors4stereotriangle`

**Signature**: `def colors4stereotriangle(resolution = 1)`

Generate RGB colors for IPF coloring.

---

### `equalarea2xyz`

**Signature**: `def equalarea2xyz(projdir)`

Convert 2D equal-area projection to 3D unit vectors.
    Inverse of equalarea_directions().

---

### `equalarea_directions`

**Signature**: `def equalarea_directions(dirs)`

Project 3D direction vectors onto 2D equal-area projection plane (Schmidt net).
    Uses Lambert azimuthal equal-area projection: r = √2 × sin(θ/2), where θ is the
    angle from the north pole. Preserves area, making it ideal for texture analysis.
    Maximum radius is √2 for hemisphere projection.

---

### `equalarea_intotriangle`

**Signature**: `def equalarea_intotriangle(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False)`

Equal-area projection into standard triangle.

---

### `equalarea_intotriangle_fast`

**Signature**: `def equalarea_intotriangle_fast(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False,Rin=None,symops=None)`

Fast equal-area projection into triangle.

---

### `equalarea_planes`

**Signature**: `def equalarea_planes(normals,arclength=360.,iniangle=0.,hemisphere="both")`

Project plane traces onto equal-area (Schmidt) net.

---

### `equivalent_elements`

**Signature**: `def equivalent_elements(element,lattice)`

Find symmetrically equivalent elements.

---

### `filled_colored_stereotriangle`

**Signature**: `def filled_colored_stereotriangle(basedirs=False,resolution = 1, markersize=1,ax=None,**kwargs)`

Draw filled triangle with smooth IPF coloring.

---

### `fullcirc_hist`

**Signature**: `def fullcirc_hist(Mats, Dr=[0,0,1], symops=None, equalarea=False, scale='sqrt', nlevels=10, lvls=None,bins=128, ax=None, title=None, ret=False, 
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
                  bandwidth=None, vmax=None,colorbar=True,ticks=None, R2Proj=None,contour=True, mrd=False, **kwargs)`

No docstring available.

---

### `gcd`

**Signature**: `def gcd(a,b)`

Compute greatest common divisor of two numbers.

---

### `gcdarr`

**Signature**: `def gcdarr(arr)`

Compute GCD of all array elements.

---

### `gen_dirs_norms`

**Signature**: `def gen_dirs_norms(L, Lr, uvws,hkls, R2Proj=np.eye(3),symops=None,recsymops=None,hemisphere = "upper", **kwargs)`

Generate all symmetry-equivalent crystal directions and plane normals with projections.
    Creates comprehensive lists of crystallographic directions [uvw] and normals (hkl) including
    all symmetry equivalents. For each direction/normal, computes:
    - 3D vector in Cartesian coordinates
    - Equal-area projection coordinates
    - Stereographic projection coordinates
    - Plane traces (for normals)

---

### `genori`

**Signature**: `def genori(dangle=1.0,hemi='both', tol=1e-2, rot=np.eye(3), half='no')`

Generate a set of orientations by rotating around two perpendicular axes.
    Creates a grid of directions by rotating around Z-axis (φ₁: 0-360°) and 
    Y-axis (φ₂: 0-180°), then optionally applies rotation matrix and filters
    by hemisphere. Useful for generating uniform orientation distributions.

---

### `genoritri`

**Signature**: `def genoritri(resolution=1.0, mesh="spherified_cube_edge")`

Generate orientation grid for stereographic triangle using orix library.
    Creates a uniform grid of crystallographic directions suitable for plotting
    in the standard stereographic triangle. Uses diffsims beam direction generator
    with orix quaternion/vector operations.

---

### `genprojgrid`

**Signature**: `def genprojgrid(oris,gdata=None,nump=1001,proj='equalarea',method2='linear',gdout=False,poris=None,minmax='notfull')`

Generate a regular 2D grid for projected orientations with optional density data.
    Creates a square grid covering the projection space and interpolates orientation
    data onto this grid. Essential for creating smooth contour plots, density maps,
    and pole figures.

---

### `ipf`

**Signature**: `def ipf(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20,syms=True)`

Generate inverse pole figure from Euler angles.

---

### `iszero`

**Signature**: `def iszero(a)`

Check if value is approximately zero (|a| < 1e-9).

---

### `pf`

**Signature**: `def pf(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20,syms=True,s=50,facecolor=(210./255.,235./255.,255./255.),plot=True)`

Generate pole figure from Euler angles.

---

### `pf_cmap`

**Signature**: `def pf_cmap(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20,syms=True,s=50,NoCont=10,GridSize=1000,cmap='jet',method='cubic')`

Generate pole figure with density colormap.

---

### `pf_cmap02`

**Signature**: `def pf_cmap02(GridX,GridY,GridR,GridPhi,AreaRatio,Intensity,NoCont=10,GridSize=1000,cmap='jet',method='cubic')`

Create contour pole figure with colormap (internal function).

---

### `pf_cmap_cscale`

**Signature**: `def pf_cmap_cscale(fig,ax2,cmin,cmax,cmap)`

Add colorbar to pole figure.

---

### `rp2xyz`

**Signature**: `def rp2xyz(r,p)`

Convert polar coordinates (r, phi) to 3D Cartesian coordinates.

---

### `schmidt_regular_grid`

**Signature**: `def schmidt_regular_grid(ax,Na=72,Nr=20,plot=True)`

Create regular grid on Schmidt net.

---

### `schmidtnet`

**Signature**: `def schmidtnet(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.))`

Draw equal-area (Schmidt) net - full circle.

---

### `schmidtnet_half`

**Signature**: `def schmidtnet_half(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.))`

Draw equal-area (Schmidt) net - upper hemisphere.

---

### `spher2xyz`

**Signature**: `def spher2xyz(polar_angle,azimuth_angle,deg=False)`

Convert spherical coordinates to Cartesian coordinates.
    Inverse transformation of xyz2spher. Converts spherical coordinates (θ, φ) to
    Cartesian (x, y, z) using standard physics convention.

---

### `stereo2xyz`

**Signature**: `def stereo2xyz(projdir)`

Convert 2D stereographic projection to 3D unit vectors.
    Inverse of stereoprojection_directions().

---

### `stereoprojection_directions`

**Signature**: `def stereoprojection_directions(dirs)`

Project 3D direction vectors onto 2D stereographic projection plane (Wulff net).
    Uses stereographic projection from south pole: r = tan(θ/2), where θ is the
    angle from the north pole. Points on lower hemisphere project outside unit circle.

---

### `stereoprojection_intotriangle`

**Signature**: `def stereoprojection_intotriangle(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False,Rin=None,symops=None)`

Map directions into standard stereographic triangle.

---

### `stereoprojection_intotriangle_fast`

**Signature**: `def stereoprojection_intotriangle_fast(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False,Rin=None,symops=None)`

Fast mapping of directions into standard triangle.

---

### `stereoprojection_intotriangle_ini`

**Signature**: `def stereoprojection_intotriangle_ini(dirs,eps=1.0e-5)`

Map directions into standard triangle (initial version).

---

### `stereoprojection_planes`

**Signature**: `def stereoprojection_planes(normals,arclength=360.,iniangle=0.,hemisphere='both',getpoints=False,R=np.eye(3))`

Project plane traces onto stereographic projection.

---

### `stereotriangle`

**Signature**: `def stereotriangle(ax=None,basedirs=False,equalarea=False,grid=False,resolution=None,gridmarkersize=None,gridmarkercol=None,gridzorder=None,mesh=False)`

Draw standard stereographic triangle for cubic system.

---

### `stereotriangle_colors`

**Signature**: `def stereotriangle_colors(proj_Ds)`

Get IPF colors from projection coordinates.

---

### `stereotriangle_colors_from_d_IPF`

**Signature**: `def stereotriangle_colors_from_d_IPF(d_IPF)`

Get IPF colors for crystal directions.

---

### `stereotriangle_colors_from_eumats_dir`

**Signature**: `def stereotriangle_colors_from_eumats_dir(eumats,d=[1,0,0])`

Get IPF colors from orientation matrices.

---

### `vector2miller`

**Signature**: `def vector2miller(arr, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3)`

Convert vector to Miller indices using GCD reduction.

---

### `vector2millerround`

**Signature**: `def vector2millerround(v, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3)`

Convert vector to Miller indices with rounding.

---

### `vectors2miller`

**Signature**: `def vectors2miller(V, MIN=True, Tol=1e-9,tol=1e5,text=False)`

Convert multiple vectors to Miller indices.

---

### `wulffnet`

**Signature**: `def wulffnet(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.))`

Draw stereographic (Wulff) net - full circle.

---

### `wulffnet_half`

**Signature**: `def wulffnet_half(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.))`

Draw stereographic (Wulff) net - upper hemisphere.

---

### `wulffnet_quarter`

**Signature**: `def wulffnet_quarter(ax=None,basedirs=False)`

Draw quarter stereographic (Wulff) net.

---

### `wulffnet_regular_grid`

**Signature**: `def wulffnet_regular_grid(ax,dangle,dirout=False, plot=True)`

Create regular angular grid on Wulff net.

---

### `wulffnet_regular_grid`

**Signature**: `def wulffnet_regular_grid(ax,dangle)`

Create regular angular grid on Wulff net.

---

### `xyz2spher`

**Signature**: `def xyz2spher(xyz,deg=False)`

Convert Cartesian coordinates to spherical coordinates (polar and azimuthal angles).
    Transforms 3D Cartesian vectors to spherical coordinate system using standard
    physics convention: θ (polar/colatitude angle from +Z axis), φ (azimuthal angle in XY plane).

---


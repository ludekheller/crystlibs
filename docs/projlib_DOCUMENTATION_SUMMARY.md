# PROJLIB - Documentation Summary

**Module**: `projlib.py`  
**Purpose**: Stereographic Projections  
**Total Functions**: 45  
**Last Updated**: December 08, 2025

---

## Function Overview

### colored_stereotriangle

`def colored_stereotriangle(basedirs=False,resolution = 1, markersize=1):`

Draw stereographic triangle with IPF coloring.

### colors4stereotriangle

`def colors4stereotriangle(resolution = 1):`

Generate RGB colors for IPF coloring.

### equalarea2xyz

`def equalarea2xyz(projdir):`

Convert 2D equal-area projection to 3D unit vectors.

### equalarea_directions

`def equalarea_directions(dirs):`

Project 3D direction vectors onto 2D equal-area projection plane (Schmidt net).

### equalarea_intotriangle

`def equalarea_intotriangle(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False):`

Equal-area projection into standard triangle.

### equalarea_intotriangle_fast

`def equalarea_intotriangle_fast(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False,Rin=None,symops=None):`

Fast equal-area projection into triangle.

### equalarea_planes

`def equalarea_planes(normals,arclength=360.,iniangle=0.,hemisphere="both"):`

Project plane traces onto equal-area (Schmidt) net.

### equivalent_elements

`def equivalent_elements(element,lattice):`

Find symmetrically equivalent elements.

### filled_colored_stereotriangle

`def filled_colored_stereotriangle(basedirs=False,resolution = 1, markersize=1,ax=None,**kwargs):`

Draw filled triangle with smooth IPF coloring.

### fullcirc_hist

`def fullcirc_hist(Mats, Dr=[0,0,1], symops=None, equalarea=False, scale='sqrt', nlevels=10, lvls=None,bins=128, ax=None, title=None, ret=False, kernel=False,  weights=None,Lim=None,interp=True,interpn=1000, smooth=False, vmin=None, bandwidth=None, vmax=None,colorbar=True,ticks=None, R2Proj=None,contour=True, mrd=False, **kwargs):`

Create a pole figure with density contours from orientation matrices.

### gcd

`def gcd(a,b):`

Compute greatest common divisor of two numbers.

### gcdarr

`def gcdarr(arr):`

Compute GCD of all array elements.

### gen_dirs_norms

`def gen_dirs_norms(L, Lr, uvws,hkls, R2Proj=np.eye(3),symops=None,recsymops=None,hemisphere = "upper", **kwargs):`

Generate all symmetry-equivalent crystal directions and plane normals with projections.

### genori

`def genori(dangle=1.0,hemi='both', tol=1e-2, rot=np.eye(3), half='no'):`

Generate a set of orientations by rotating around two perpendicular axes.

### genoritri

`def genoritri(resolution=1.0, mesh="spherified_cube_edge"):`

Generate orientation grid for stereographic triangle using orix library.

### genprojgrid

`def genprojgrid(oris,gdata=None,nump=1001,proj='equalarea',method2='linear',gdout=False,poris=None,minmax='notfull'):`

Generate a regular 2D grid for projected orientations with optional density data.

### ipf

`def ipf(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20,syms=True):`

Generate inverse pole figure from Euler angles.

### iszero

`def iszero(a):`

Check if value is approximately zero (|a| < 1e-9).

### pf

`def pf(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20,syms=True,s=50,facecolor=(210./255.,235./255.,255./255.),plot=True):`

Generate pole figure from Euler angles.

### pf_cmap

`def pf_cmap(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20,syms=True,s=50,NoCont=10,GridSize=1000,cmap='jet',method='cubic'):`

Generate pole figure with density colormap.

### pf_cmap02

`def pf_cmap02(GridX,GridY,GridR,GridPhi,AreaRatio,Intensity,NoCont=10,GridSize=1000,cmap='jet',method='cubic'):`

Create contour pole figure with colormap (internal function).

### pf_cmap_cscale

`def pf_cmap_cscale(fig,ax2,cmin,cmax,cmap):`

Add colorbar to pole figure.

### rp2xyz

`def rp2xyz(r,p):`

Convert polar coordinates (r, phi) to 3D Cartesian coordinates.

### schmidt_regular_grid

`def schmidt_regular_grid(ax,Na=72,Nr=20,plot=True):`

Create regular grid on Schmidt net.

### schmidtnet

`def schmidtnet(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.)):`

Draw equal-area (Schmidt) net - full circle.

### schmidtnet_half

`def schmidtnet_half(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.)):`

Draw equal-area (Schmidt) net - upper hemisphere.

### spher2xyz

`def spher2xyz(polar_angle,azimuth_angle,deg=False):`

Convert spherical coordinates to Cartesian coordinates.

### stereo2xyz

`def stereo2xyz(projdir):`

Convert 2D stereographic projection to 3D unit vectors.

### stereoprojection_directions

`def stereoprojection_directions(dirs):`

Project 3D direction vectors onto 2D stereographic projection plane (Wulff net).

### stereoprojection_intotriangle

`def stereoprojection_intotriangle(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False,Rin=None,symops=None):`

Map directions into standard stereographic triangle.

### stereoprojection_intotriangle_fast

`def stereoprojection_intotriangle_fast(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False,Rin=None,symops=None):`

Fast mapping of directions into standard triangle.

### stereoprojection_intotriangle_ini

`def stereoprojection_intotriangle_ini(dirs,eps=1.0e-5):`

Map directions into standard triangle (initial version).

### stereoprojection_planes

`def stereoprojection_planes(normals,arclength=360.,iniangle=0.,hemisphere='both',getpoints=False,R=np.eye(3)):`

Project plane traces onto stereographic projection.

### stereotriangle

`def stereotriangle(ax=None,basedirs=False,equalarea=False,grid=False,resolution=None,gridmarkersize=None,gridmarkercol=None,gridzorder=None,mesh=False):`

Draw standard stereographic triangle for cubic system.

### stereotriangle_colors

`def stereotriangle_colors(proj_Ds):`

Get IPF colors from projection coordinates.

### stereotriangle_colors_from_d_IPF

`def stereotriangle_colors_from_d_IPF(d_IPF):`

Get IPF colors for crystal directions.

### stereotriangle_colors_from_eumats_dir

`def stereotriangle_colors_from_eumats_dir(eumats,d=[1,0,0]):`

Get IPF colors from orientation matrices.

### vector2miller

`def vector2miller(arr, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3):`

Convert vector to Miller indices using GCD reduction.

### vector2millerround

`def vector2millerround(v, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3):`

Convert vector to Miller indices with rounding.

### vectors2miller

`def vectors2miller(V, MIN=True, Tol=1e-9,tol=1e5,text=False):`

Convert multiple vectors to Miller indices.

### wulffnet

`def wulffnet(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.)):`

Draw stereographic (Wulff) net - full circle.

### wulffnet_half

`def wulffnet_half(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.)):`

Draw stereographic (Wulff) net - upper hemisphere.

### wulffnet_quarter

`def wulffnet_quarter(ax=None,basedirs=False):`

Draw quarter stereographic (Wulff) net.

### wulffnet_regular_grid

`def wulffnet_regular_grid(ax,dangle):`

Create regular angular grid on Wulff net.

### xyz2spher

`def xyz2spher(xyz,deg=False):`

Convert Cartesian coordinates to spherical coordinates (polar and azimuthal angles).


---

**Total**: 45 functions

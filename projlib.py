#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Crystallographic Projection Library Module

This module provides comprehensive utilities for crystallographic stereographic projections including:
- Stereographic projections (equal-area Schmidt net, equal-angle Wulff net)
- Direction and plane projections onto stereographic nets
- Orientation generation and grid creation
- Coordinate transformations (Cartesian ↔ spherical ↔ projection)
- Pole figure and inverse pole figure generation
- Symmetry operations and equivalent direction finding
- Stereographic triangle projections for cubic systems
- Miller indices operations and vector conversions
- Histogram and density plotting for orientation data
- Color mapping for inverse pole figures

Key Projection Formulas:
- Equal-area (Schmidt): r = √2 × sin(θ/2), where θ is angle from north pole
- Stereographic (Wulff): r = tan(θ/2)
- Inverse projections: θ = 2×arctan(r) for Wulff, θ = 2×arcsin(r/√2) for Schmidt

Module Functions by Category:

Orientation Generation:
- genori(): Generate orientation grid by rotation around two axes
- genoritri(): Generate cubic triangle orientation grid using orix

Coordinate Transformations:
- xyz2spher(): Cartesian → spherical (θ, φ)
- spher2xyz(): Spherical → Cartesian
- stereoprojection_directions(): 3D → 2D stereographic projection
- equalarea_directions(): 3D → 2D equal-area projection
- stereo2xyz(): 2D stereographic → 3D
- equalarea2xyz(): 2D equal-area → 3D

Grid and Plot Generation:
- genprojgrid(): Generate regular 2D grid for density plotting
- fullcirc_hist(): Create pole figure with density contours
- wulffnet(), schmidtnet(): Draw projection net circles
- stereotriangle(): Draw stereographic triangle

Crystallography:
- gen_dirs_norms(): Generate crystal directions/normals with symmetry
- equalarea_planes(), stereoprojection_planes(): Project plane traces
- stereoprojection_intotriangle(): Map directions to standard triangle
- vector2miller(), vectors2miller(): Convert to Miller indices

Color Mapping:
- colors4stereotriangle(): Generate IPF colors for directions
- colored_stereotriangle(): Draw colored inverse pole figure triangle
- ipf(): Create full inverse pole figure with data

Created on Wed Sep 11 10:13:57 2019
@author: lheller
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from orilib import *
from matplotlib.patches import Wedge
from functools import reduce
try:
    import random
except:
    pass
try:
    import scipy
    from scipy.interpolate import griddata

except:
    pass
try:
    from orix.quaternion import Rotation
except:
    pass
try:
    from orix.vector import Vector3d
except:
    pass
try:
    from orix.projections import StereographicProjection
except:
    pass
try:
    from diffsims.generators.rotation_list_generators import get_beam_directions_grid
except:
    pass
try:
    from matplotlib.colors import ListedColormap
except:
    pass
from spherical_kde import SphericalKDE

try:
    from wand.image import Image
except:
    pass

def genoritri(resolution=1.0, mesh="spherified_cube_edge"):
    """
    Generate orientation grid for stereographic triangle using orix library.
    
    Creates a uniform grid of crystallographic directions suitable for plotting
    in the standard stereographic triangle. Uses diffsims beam direction generator
    with orix quaternion/vector operations.
    
    Input:
        resolution: float - Angular resolution in degrees (default: 1.0)
                           Smaller values create finer grids but more points
        mesh: str - Mesh generation method (default: "spherified_cube_edge")
                   Options:
                   - "spherified_cube_edge": Uniform angular spacing
                   - "spherified_cube_corner": Corner-focused spacing
                   - "cubochoric": Volume-preserving grid
    
    Output:
        trioris: numpy array (3, N) - Direction vectors in Cartesian coordinates [x, y, z]
                                       Each column is a unit vector
    
    Usage Example:
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
    """
    grid_cub = get_beam_directions_grid("cubic", resolution, mesh=mesh)
    trioris = Rotation.from_euler(np.deg2rad(grid_cub))*Vector3d.zvector()
    return trioris.data.T

#Set of all orientations in space defined by rotation around 2 perpendicular axis and an angle resolution

def genori(dangle=1.0,hemi='both', tol=1e-2, rot=np.eye(3), half='no'):
    """
    Generate a set of orientations by rotating around two perpendicular axes.
    
    Creates a grid of directions by rotating around Z-axis (φ₁: 0-360°) and 
    Y-axis (φ₂: 0-180°), then optionally applies rotation matrix and filters
    by hemisphere. Useful for generating uniform orientation distributions.
    
    Input:
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
    
    Output:
        oris: numpy array (3, N) - Direction vectors [x, y, z]
                                    Each column is a unit vector
                                    N depends on dangle and filtering options
    
    Usage Example:
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
    """
    #dangle=1.
    Phi1=np.linspace(0.,360.-dangle,int(360./dangle))
    Phi2=np.linspace(0.,180.-dangle,int(180./dangle))
    GridX=[];
    GridY=[];
    Dc=np.array([1,0,0]);
    oris=[]
    for phi1 in Phi1:
        for phi2 in Phi2:        
            RotZ = active_rotation(phi1, 'z', deg=True) 
            RotY = active_rotation(phi2, 'y', deg=True)
            oris.append(RotY.dot(RotZ.dot(Dc)))
    oris=rot.dot(np.asarray(oris).T)
    if half=='upper':
        poris=equalarea_directions(oris)
        oris=oris[:,poris[1,:]>-tol]
    elif half=='lower':
        poris=equalarea_directions(oris)
        oris=oris[:,poris[1,:]<tol]
    if hemi=='upper':
        return oris[:,oris[2,:]>-tol]
    elif hemi=='lower':
        return oris[:,oris[2,:]<tol]
    else:    
        return oris

#generate regular grid as masked array npxnp covering the whole range of projected oris

def xyz2spher(xyz,deg=False):
    """
    Convert Cartesian coordinates to spherical coordinates (polar and azimuthal angles).
    
    Transforms 3D Cartesian vectors to spherical coordinate system using standard
    physics convention: θ (polar/colatitude angle from +Z axis), φ (azimuthal angle in XY plane).
    
    Input:
        xyz: numpy array (N, 3) - Cartesian coordinates [x, y, z]
                                   Each row is one point
        deg: bool - If True, return angles in degrees; if False, radians (default: False)
    
    Output:
        polar_angle: numpy array (N,) - Polar angle θ (colatitude from z-axis)
                                        Range: [0, π] radians or [0, 180°]
        azimuth_angle: numpy array (N,) - Azimuthal angle φ (longitude in xy-plane)
                                          Range: [-π, π] radians or [-180°, 180°]
    
    Usage Example:
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
    """
    polar_angle=np.arctan2(np.sqrt(xyz[:,0]**2+xyz[:,1]**2),xyz[:,2])
    azimuth_angle=np.arctan2(xyz[:,1],xyz[:,0])        
    if deg:
        polar_angle*=180/np.pi
        azimuth_angle*=180/np.pi    
    return polar_angle,azimuth_angle 

def spher2xyz(polar_angle,azimuth_angle,deg=False):
    """
    Convert spherical coordinates to Cartesian coordinates.
    
    Inverse transformation of xyz2spher. Converts spherical coordinates (θ, φ) to
    Cartesian (x, y, z) using standard physics convention.
    
    Input:
        polar_angle: float or array - Polar angle θ (colatitude from z-axis)
                                      Expected range: [0, π] or [0, 180°]
        azimuth_angle: float or array - Azimuthal angle φ (longitude in xy-plane)
                                        Expected range: [-π, π] or [-180°, 180°]
        deg: bool - If True, input angles are in degrees; if False, radians (default: False)
    
    Output:
        xyz: numpy array (N, 3) - Cartesian coordinates [x, y, z]
                                   Each row is one point (unit vector if inputs were from unit sphere)
    
    Usage Example:
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
    """
    if deg:
        #polar_angle*=np.pi/180
        #azimuth_angle*=np.pi/180   
        z=np.cos(polar_angle*np.pi/180)
        xy=np.sin(polar_angle*np.pi/180)
        y=np.sin(azimuth_angle*np.pi/180)*xy
        x=np.cos(azimuth_angle*np.pi/180)*xy
    else:
        z=np.cos(polar_angle)
        xy=np.sin(polar_angle)
        y=np.sin(azimuth_angle)*xy
        x=np.cos(azimuth_angle)*xy

    return np.vstack((x,y,z)).T


def genprojgrid(oris,gdata=None,nump=1001,proj='equalarea',method2='linear',gdout=False,poris=None,minmax='notfull'):
    """
    Generate a regular 2D grid for projected orientations with optional density data.
    
    Creates a square grid covering the projection space and interpolates orientation
    data onto this grid. Essential for creating smooth contour plots, density maps,
    and pole figures.
    
    Input:
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
    
    Output:
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
    
    Usage Example:
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
    """
    #project


    if poris is None:
        if proj=='stereo':
            poris=stereoprojection_directions(oris)
        elif proj=='equalarea':
            poris=equalarea_directions(oris)
        
    
    if minmax=='full':
        if proj=='stereo':
            minv=-1
            maxv=1
            minv2=-1
            maxv2=1        
        elif proj=='equalarea':
            minv=-np.sqrt(2)
            maxv=np.sqrt(2)
            minv2=-np.sqrt(2)
            maxv2=np.sqrt(2)
    else:
        minv=min(poris[0,:])
        maxv=max(poris[0,:])
        minv2=min(poris[1,:])
        maxv2=max(poris[1,:])

    
    #make regular grid
    grid_x, grid_y = np.mgrid[minv:maxv:nump*1j, minv2:maxv2:nump*1j]
    #Grided data
    gdout=True
    if gdata is None:
        gdata=poris[1,:].flatten()*0+1
        gdout=False
    gdata=griddata((poris[0,:].flatten(),poris[1,:].flatten()), gdata, (grid_x, grid_y), method=method2)
    #numerical mask of the hemisphere
    nummask = np.nan_to_num(gdata*0+1,nan=0.0)
    #boolean mask of the hemisphere
    mask = np.nan_to_num(gdata*0,nan=1.0).astype(bool)
    mask2 = np.isnan(gdata)
    #print(mask2)
    #masked gridded data
    gdata=np.nan_to_num(gdata,nan=0.0)
    grid_z =np.ma.array(gdata,mask=mask,fill_value=0)
    #grid_z1 = gdata*nummask
    if gdout:
        return grid_x,grid_y, grid_z, nummask, mask
    else:
        return grid_x,grid_y, nummask, mask



def gen_dirs_norms(L, Lr, uvws,hkls, R2Proj=np.eye(3),symops=None,recsymops=None,hemisphere = "upper", **kwargs):
    """
    Generate all symmetry-equivalent crystal directions and plane normals with projections.
    
    Creates comprehensive lists of crystallographic directions [uvw] and normals (hkl) including
    all symmetry equivalents. For each direction/normal, computes:
    - 3D vector in Cartesian coordinates
    - Equal-area projection coordinates
    - Stereographic projection coordinates
    - Plane traces (for normals)
    
    Input:
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
    
    Output:
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
    
    Usage Example:
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
    """
    #generates all symmetry equivalent directions and normals and their projections
    normals=[]
    
    for hkl in hkls:
            nv=R2Proj.dot(Lr.dot(hkl))
            nv/=np.linalg.norm(nv)
            isin=False
            for n in normals:
                if list(n['hkl']) == list(hkl):#np.linalg.norm(n['vector']-nvs)<1e-10:
                    isin=True
                    break
            if not isin:
                
                normals.append({'vector':nv,'hkl':hkl,'hklf':hkl,'label':str(hkl).replace('[','(').replace(']',')').replace(' ',''),
                                'equalarea':equalarea_directions(nv),'stereo':stereoprojection_directions(nv),
                               'equalarea plane':equalarea_planes(nv,hemisphere=hemisphere),
                                'stereo plane':stereoprojection_planes(nv,hemisphere=hemisphere),'textshift':[0,0]})
                if not recsymops is None:
                    for rs in recsymops:
                        
                        hklsym=np.round(rs.dot(hkl))
                        isin=False
                        for n in normals:
                            if True:#n['hklf']==hkl:
                                if list(hklsym) == list(n['hkl']):#np.linalg.norm(n['vector']-nvs)<1e-10:
                                    isin=True
                                    break
                        if not isin:
                            #print(list(hklsym))
                            #print(hkl)
                            #print('------------------------------')
                            nvs=R2Proj.dot(rs.dot(Lr.dot(hkl)))
                            nvs/=np.linalg.norm(nvs)
                            normals.append({'vector':nvs,'hkl':list(hklsym),'hklf':hkl,'label':str(hklsym).replace('[','(').replace(']',')').replace(' ',''),
                                            'equalarea':equalarea_directions(nvs),'stereo':stereoprojection_directions(nvs),
                                           'equalarea plane':equalarea_planes(nvs,hemisphere=hemisphere),
                                           'stereo plane':stereoprojection_planes(nv,hemisphere=hemisphere),'textshift':[0,0]})
            #print('=============================================')  
    dirs=[]
    for uvw in uvws:
            dv=R2Proj.dot(L.dot(uvw))
            dv/=np.linalg.norm(dv)
            isin=False
            for d in dirs:
                if list(d['uvw']) == list(uvw):#np.linalg.norm(n['vector']-nvs)<1e-10:
                    isin=True
                    break
            if not isin:
                dirs.append({'vector':dv,'uvw':uvw,'uvwf':uvw,'label':str(uvw),
                                'equalarea':equalarea_directions(dv),'stereo':stereoprojection_directions(dv),'textshift':[0,0]})
                if not symops is None:
                    for rs in symops:
                        #dvs=rs.dot(dv)
                        uvwsym=np.round(rs.dot(uvw))
                        isin=False
                        for d in dirs:
                            if True:#d['uvwf']==uvw:
                                if list(uvwsym) == list(d['uvw']):#np.linalg.norm(d['vector']-dvs)<1e-10:
                                    isin=True
                                    break
                        if not isin:
                            dvs=R2Proj.dot(rs.dot(L.dot(uvw)))
                            dvs/=np.linalg.norm(dvs)
                            dirs.append({'vector':dvs,'uvw':list(uvwsym),'uvwf':uvw,'label':str(uvwsym).replace('[','(').replace(']',')').replace(' ',''),
                                            'equalarea':equalarea_directions(dvs),'stereo':stereoprojection_directions(dvs),'textshift':[0,0]})
                        
    return dirs,normals



def fullcirc_hist(Mats, Dr=[0,0,1], symops=None, equalarea=False, scale='sqrt', nlevels=10, lvls=None,bins=128, ax=None, title=None, ret=False, kernel=False,  weights=None,Lim=None,interp=True,interpn=1000, smooth=False, vmin=None,
                  bandwidth=None, vmax=None,colorbar=True,ticks=None, R2Proj=None,contour=True, mrd=False, **kwargs):
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
    #generating inverse poles of Dr from orientation matrices Mats
    Dr=np.array(Dr)/np.linalg.norm(Dr)
    if symops is None:
        Mr = np.reshape(Mats, (Mats.shape[0]*Mats.shape[1],Mats.shape[2]))
        data = Mr.dot(Dr)
        data = np.reshape(data,(int(data.shape[0]/3),3)).T
    else:
        Drs = set([tuple(v) for v in np.dot(symops, Dr)])
        Drs=np.asarray(list(Drs))
        data=np.tensordot(Mats, Drs.T, axes=[[-1], [-2]]).transpose([0, 2, 1])
        data = data.reshape(-1,3).T
        
    #matrix Rmat time aray of matrices Mats[N,3,3]: (to be verified!)
    #Mr=np.reshape(Mats, (Mats.shape[0]*Mats.shape[1],Mats.shape[2])).T
    #Mr=Rmat.T.dot(Mr)
    #Mats=Mr.T.reshape(-1,3,3)
    
    
    
    #print(data.shape)
    #Changing the coordinate system to align projection x,y axes as we want
    if R2Proj is not None:
        data = R2Proj.dot(data)
    dataini=data.copy() 
    #inverting vectors projecting to lower hemisphere
    data[:,data[2,:]<0]=-1*data[:,data[2,:]<0]
    data = data[:,data[2,:]>=0]
    #plotting projection circles and setting corresponding radius of the circe (lim)
    #calculating projection of data (data->spsel)
    if equalarea:
        proj='equalarea'
        lim=np.sqrt(2)
        spsel = equalarea_directions(data)
        if Lim=='half':
            fig,ax = schmidtnet_half(ax=ax,basedirs=False)     
        elif Lim=='tri':
            fig,ax = stereotriangle(ax=ax,basedirs=False,equalarea=True)   
        else:
            fig,ax = schmidtnet(ax=ax,basedirs=False)
    else:
        proj='stereo'
        lim=1.
        spsel = stereoprojection_directions(data)
        if Lim=='half':
            fig,ax = wulffnet_half(ax=ax,basedirs=False)       
        elif Lim=='tri':
            fig,ax = stereotriangle(ax=ax,basedirs=False,equalarea=False)
        else:
            fig,ax = wulffnet(ax=ax,basedirs=False)
    if title !='':
        ax.title.set_text(title)       
        
    #calculating weighted histogram of the projected data withing the limit of the circle 
    hist, yedges, xedges = np.histogram2d(spsel[1,:], spsel[0,:], bins=bins,range=[[-lim, lim], [-lim, lim]], weights= weights)
    #histraw=hist.copy()
    
#    if scale=='sqrt':
#        hist = hist**0.5
#    elif scale=='log':
#        hist = np.log(hist)
        
    #nlevels=nlevels
    #lvls = np.linspace(0, np.max(hist), nlevels)
    
    #generating X, Y grid aligned with bins
    X, Y = np.meshgrid((xedges[:-1] + xedges[1:])/2.,
                       (yedges[:-1] + yedges[1:])/2.)
    
    #circle = X**2 + Y**2 >= 2.0
    #histraw[circle] = np.nan  
    #c=ax.contourf(histraw, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], **kwargs)
    if kernel:
        #spherical coordinates theta (elevation), phi of inverse poles
        xy = dataini[0,:]**2 + dataini[1,:]**2
        theta = np.arctan2(np.sqrt(xy), dataini[2,:])#elevation
        phi=np.arctan2(dataini[1,:], dataini[0,:])

        if bandwidth is None:
            bandwidth=0.15
        #calculation of the logarithm of kernel density function of inverse poles using corresponding spherical coordinates theta (elevation), phi
        #I do not know why SphericalKDE provides logarithm
        logkde = SphericalKDE(phi, theta, bandwidth=bandwidth,weights=weights)
        
        #generation of orientation whithin whole space where we will evaluate logkde
        dangle=5.0
        oris=genori(dangle=dangle,hemi='both', tol=1e-2, rot=np.eye(3), half='no')
        xy = oris[0,:]**2 + oris[1,:]**2
        thtg = np.arctan2(np.sqrt(xy), oris[2,:])#elevation==polar angle
        phig=np.arctan2(oris[1,:], oris[0,:])#azimuth angle
        #evaluation of kde for all orientations
        H=np.exp(logkde(phig, thtg))    
        
        if mrd:
            #Normalization of H for MRD
            #normalizing data to multiples of random distribution
            #hist/(sum(hist)*pxarea)*numpixels*pxarea  
            #integral Hds=1
            Norm=(H*np.sin(thtg)*(dangle*np.pi/180)**2/4/np.pi).sum()
            H/=Norm
            #H=H/H.sum()*H.shape[0]
            #check normalization
            #print((H*np.sin(thtg)*(dangle*np.pi/180)**2/4/np.pi).sum())
            
        #interpolation of kde over regular square grid on the projection circle - it is masked or nan outside the circle
        X,Y, hist, nummask, mask=genprojgrid(oris[:,oris[2,:]>=0],gdata=H[oris[2,:]>=0],nump=1001,proj=proj,method2='linear',gdout=True)
        #if mrd:
        #    hist=hist/np.nansum(hist)*np.where(~np.isnan(hist.data))[0].shape[0]
        bins=1001
        #if True:
        #hist=hist.data
        #k = scipy.stats.gaussian_kde([spsel[0,:], spsel[1,:]],weights= weights)
        #hist = k(np.vstack([X.flatten(), Y.flatten()]))
        #hist = hist.reshape(X.shape)
        #print(hist)
    #print(hist)
    #refine histogram by interpolation
    if interp and not kernel:
        #interpolation of histogram data to increase spatial resolution thus making nicer. It is faster and smoother than increasing number of bins 
        xedges=np.linspace(((xedges[:-1] + xedges[1:])/2.)[0],((xedges[:-1] + xedges[1:])/2.)[-1],interpn)
        yedges=np.linspace(((yedges[:-1] + yedges[1:])/2.)[0],((yedges[:-1] + yedges[1:])/2.)[-1],interpn)
        grid_x, grid_y = np.meshgrid(xedges,yedges)    
    
        hist = scipy.interpolate.griddata((X.flatten(), Y.flatten()), hist.flatten(), (grid_x, grid_y), method='linear')
        bins=interpn
        X=grid_x
        Y=grid_y

    #print(hist)
    #hist[np.where((hist<0) & (np.abs(hist)<1e-10))]=0
    #histogram scaling - should not be used if MRD required
    hist[np.where(hist<0)]=0
    if scale=='sqrt':
        hist = hist**0.5
    elif scale=='log':
        hist = np.log(hist)
    #masking data on the grid points outside the projection circle
    if equalarea:
        circle = X**2 + Y**2 >= 2.0
        xlim=[-np.sqrt(2)*1.05,np.sqrt(2)*1.05]
        ylim=[-np.sqrt(2)*1.05,np.sqrt(2)*1.05]
    else:
        circle = X**2 + Y**2 >= 1
        xlim=[-1.05,1.05]
        ylim=[-1.05,1.05]
    hist[circle] = np.nan  
    
    #normalizing data to multiples of random distribution
    #hist/(sum(hist)*pxarea)*numpixels*pxarea
    if mrd and not kernel:
        hist=hist/np.nansum(hist)*np.where(~np.isnan(hist.data))[0].shape[0]

    #calculating levels for isolines and contourse    
    nlevels=nlevels
    if vmin is None:
        vmin=0
    if vmax is None:
        vmax=np.max(hist[~np.isnan(hist)])
    if lvls is None:
        lvls = np.linspace(vmin, vmax, nlevels)
    else:
        nlevels=lvls.shape[0]
    kwargs = {}
    kwargs['levels'] = lvls[1:]
    
    #cutting projection area to half circle of stereotriangle
    if Lim is not None:
        if Lim=='half':
            ylim[0]=-1e-5
            cut=Y<-1e-5
            hist[cut] = np.nan 
        if Lim=='tri':
            oristri=genoritri(resolution=0.1,mesh="spherified_cube_edge")
            grid_x,grid_y, nummask, mask=genprojgrid(oristri,nump=bins,proj=proj,method2='linear',gdout=False)
            if kernel:
                hist[nummask==0]=np.nan
            else:
                hist[nummask.T==0]=np.nan
            
            
    #print(hist)
    #Only visualization below
    if kernel:
        if smooth:
            sc=ax.pcolor(X, Y,hist,vmin=vmin,vmax=vmax)
        else:
            #print('ok')
            #sc=ax.contourf(histraw, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], **kwargs)
            sc=ax.contourf(hist.T, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], **kwargs)
        if contour:
            CS=ax.contour(X, Y,hist,levels=lvls,colors='k')
            ax.clabel(CS, fontsize=9, inline=1,colors='k')

    else:
        if smooth:
            sc=ax.pcolor(X, Y,hist,vmin=vmin,vmax=vmax)
        else:
            sc=ax.contourf(hist, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], **kwargs)        
        if contour:
            CS=ax.contour(X, Y,hist,levels=nlevels,colors='k')
            ax.clabel(CS, fontsize=9, inline=1,colors='k')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if colorbar:
        pos = ax.get_position()
        fac=0.75
        cbarh=0.04
        cbar_ax = fig.add_axes([pos.width*(1-fac)/2+pos.x0, pos.y0+cbarh-0.1, pos.width*fac,cbarh])
        #cbar_ax =AX[1]
        cbar = fig.colorbar(sc, cax=cbar_ax, orientation='horizontal')   
        if mrd:
            cbar.ax.set_xlabel("MRD")
        else:
            cbar.ax.set_xlabel("Probability")
        #cbar.ax.set_xlim([0,30])
        xticks=(cbar.ax.get_xticks()[0:-1]+cbar.ax.get_xticks()[1:])/2
        #cbar.ax.set_xticks(np.round(xticks*100)/100)
        if ticks is not None:
            cbar.ax.set_xticks(ticks)

    if ret:
        return hist, xedges, yedges, fig, ax
#convert projected points into xyz

def rp2xyz(r,p):
    """
    Convert polar coordinates (r, phi) to 3D Cartesian coordinates.
    
    Input:
        r: float/array - Polar angle in degrees (0-180)
        p: float/array - Azimuthal angle in degrees (0-360)
    
    Output:
        x, y, z: tuple - Cartesian coordinates
    """
    npatan2d = lambda x,y: 180.*np.arctan2(x,y)/np.pi
    z = npcosd(r)
    xy = np.sqrt(1.-z**2)
    return xy*npsind(p),xy*npcosd(p),z


def stereoprojection_directions(dirs):

    """
    Project 3D direction vectors onto 2D stereographic projection plane (Wulff net).
    
    Uses stereographic projection from south pole: r = tan(θ/2), where θ is the
    angle from the north pole. Points on lower hemisphere project outside unit circle.
    
    Input:
        dirs: numpy array (3, N) or (3,) - Direction vectors [x, y, z]
                                           Will be automatically normalized
                                           Single vector (3,) will be reshaped to (3, 1)
    
    Output:
        proj_dirs: numpy array (3, N) - Projected coordinates [X, Y, 0]
                                         X, Y are 2D projection coordinates
                                         Third row is zeros (kept for compatibility)
    
    Usage Example:
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
    """
    #dirs = [x1,x2,...,xn;y1,y2,...,yn;z1,z2,...,zn];
    #example: dirs = [0,1,2,3;1,2,3,0;0,3,2,1]
    
    #normalize and project
    if len(dirs.shape)==1:
        dirs = np.expand_dims(dirs,axis=1)

    dirs = dirs.astype(float)
    #normalizing dirs
    dirs /= np.sqrt((dirs ** 2).sum(0))
    #check
    
    proj_dirs = np.vstack((1./(np.sign(dirs[2,:])*dirs[2,:]+1)*dirs[0,:], 
               1./(np.sign(dirs[2,:])*dirs[2,:]+1)*dirs[1,:],
                np.zeros(dirs[2,:].shape)))
    #check
#    alpha = np.arccos(np.sign(dirs[2,:])*dirs[2,:])
#    eps=1e-6;
#    idxs = np.where(alpha<eps)
#    alpha[idxs]=1.
#    beta = np.arccos(dirs[0,:]/np.sin(alpha))
#    rx=np.tan(alpha/2)*np.cos(beta)
#    ry=np.tan(alpha/2)*np.sin(beta)
#    rx[idxs]=0.
#    ry[idxs]=0.
#    abs(proj_dirs[0,:]-rx)<eps
#    abs(proj_dirs[1,:]-ry)<eps
    return proj_dirs

# def stereoprojection_intotriangle(dirs):

    """
    Map directions into standard stereographic triangle.
    
    Input:
        dirs: numpy array (3, N) - Direction vectors
        eps, geteqdirs, geteqmats, Rin, symops: Optional parameters
    
    Output:
        proj: numpy array (2, N) - Projection coordinates
    """
#     eps=1.0e-2
#     normals = np.array([-1,0,1]);
#     arclength = 40.#-np.arccos(np.sqrt(2)/np.sqrt(3))*180/np.pi;
#     proj_normals, points = stereoprojection_planes(normals,arclength=arclength,iniangle=90)
#     proj_tans = np.arctan(proj_normals[1,:]/proj_normals[0,:])

#     if len(dirs.shape)==1:
#         dirs = np.expand_dims(dirs,axis=1)

#     proj_dirs = np.zeros(dirs.shape)
#     inc=-1
#     for co,diri in enumerate(dirs.T):
#         #print('Direction {} from {}'.format(co+1,dirs.shape[1]))
# #        print('===================================================')
# #        print(diri)
# #        print('===================================================')
#         inc+=1
#         el=equivalent_elements(diri,'cubic')
#         #print(el)
#         for eli in el:
#             proj_eli = stereoprojection_directions(eli)
# #            print(eli)
# #            print(np.arctan(proj_eli[1,0]/proj_eli[0,0])-np.arccos(1./np.sqrt(3.)))
# #            print(np.arctan(proj_eli[1,0]/proj_eli[0,0]))#-np.pi/4)
# #            print((np.arccos(abs(eli[2])/np.sqrt(eli.dot(eli)))-np.arccos(1./np.sqrt(3.))))
# #            if (eli>=-eps).all() and (np.arccos(abs(eli[2])/np.sqrt(eli.dot(eli)))-np.arccos(1./np.sqrt(3.)))<eps:
# #            if (proj_eli[:,0]>=-eps).all() and (np.arccos(abs(eli[2])/np.sqrt(eli.dot(eli)))-np.arccos(1./np.sqrt(3.)))<eps:
#             if ((proj_eli[:,0])>=-eps).all():                #proj_eli = stereoprojection_directions(eli)
#                 atan=np.arctan2(proj_eli[1,0],proj_eli[0,0])
#                 if (atan-np.pi/4)<eps:
#                     idx=np.where(abs(proj_tans-atan)==min(abs(proj_tans-atan)))[0][0]
# #                    print(proj_eli[:,0].dot(proj_eli[:,0])) 
# #                    print(proj_normals[:,idx].dot(proj_normals[:,idx])) 
# #                    print((proj_eli[:,0].dot(proj_eli[:,0])-proj_normals[:,idx].dot(proj_normals[:,idx])))
# #                    print((proj_eli[:,0].dot(proj_eli[:,0])-proj_normals[:,idx].dot(proj_normals[:,idx]))<eps)
#                     if (proj_eli[:,0].dot(proj_eli[:,0])-proj_normals[:,idx].dot(proj_normals[:,idx]))<eps:
#                         proj_dirs[:,inc]=proj_eli[:,0]
# #                        print('OK')
# #                        print(proj_dirs[:,inc])
#                         break
#     return proj_dirs

#def coodinate_from_equalarea_proj(projdirs):
    #dirs = [x1,x2,...,xn;y1,y2,...,yn];
    #example: dirs = [0,1,2,3;1,2,3,0]
#    if len(projdirs.shape)==1:
#        projdirs = np.expand_dims(projdirs,axis=1)
#    z=-2+np.sqrt(8-projdirs[0,:]**2-projdirs[1,:]**2)

#    z=-2+np.sqrt()
    
    
    

def equalarea_planes(normals,arclength=360.,iniangle=0.,hemisphere="both"):
    """
    Project plane traces onto equal-area (Schmidt) net.
    
    Input:
        normals: numpy array (3, N) - Plane normal vectors
        arclength: float - Arc length in degrees (default: 360)
        iniangle: float - Starting angle in degrees (default: 0)
        hemisphere: str - 'both', 'upper', or 'lower' (default: 'both')
    
    Output:
        traces: list of arrays - Trace points for each plane
    """
    #%normals = [x1,x2,...,xn;y1,y2,...,yn;z1,z2,...,zn];
    #%varargin{1} arclength in deg
    #normals = np.transpose(np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1]]))
    #
    if len(normals.shape)==1:
        normals = np.expand_dims(normals,axis=1)

    normals = normals.astype(float)
    normals /= np.sqrt((normals ** 2).sum(0))

    proj_normals = equalarea_directions(normals)

    idxs = np.where(abs(normals[0,:])+abs(normals[1,:])==0)[0]
    
    inplanedirs = np.vstack((-normals[1,:],normals[0,:],np.zeros(normals[0,:].shape)));
    inplanedirs[:,idxs] = np.vstack((np.zeros(normals[0,idxs].shape), -normals[2,idxs],normals[1,idxs]));
    
    inplanedirs /= np.sqrt((inplanedirs ** 2).sum(0))

    thirdaxis=np.cross(normals,inplanedirs,axisa=0,axisb=0,axisc=0)
#    thirdaxis = np.vstack((normals[1,:]*inplanedirs[2,:]-normals[2,:]*inplanedirs[1,:],
#                           -1*(normals[0,:]*inplanedirs[2,:]-normals[2,:]*inplanedirs[0,:]),
#                           normals[0,:]*inplanedirs[1,:]-normals[1,:]*inplanedirs[0,:]));
    t=np.linspace(iniangle,iniangle+arclength,180*2+1)*np.pi/180;
    basicarc = np.vstack((np.cos(t),np.sin(t),np.zeros(t.shape)));
    
    proj_planes=[];
    Zdir=[]
    #print(hemisphere)
    for i in range(0,normals.shape[1]):
        Rot2Global = np.transpose(np.vstack((inplanedirs[:,i],thirdaxis[:,i],normals[:,i])));
        Ccp = np.matmul(Rot2Global,basicarc);
        Zdir=Ccp[2]
        if hemisphere == "both":            
            Ds = equalarea_directions(Ccp)
        elif hemisphere == "triangle":  
            #idxs = np.where(Ccp[2,:]>=0)[0]
            Ds=equalarea_intotriangle(Ccp)#[:,idxs])
        else:
            if hemisphere == "upper":
                idxs = np.where(Ccp[2,:]>=0)[0]
                Ds = equalarea_directions(Ccp[:,idxs])
            elif hemisphere == "lower":
                idxs = np.where(Ccp[2,:]<=0)[0]
                Ds = equalarea_directions(Ccp[:,idxs])
        proj_planes.append(Ds)
    
    if len(proj_planes)==1:
        return proj_planes[0]
    else:          
        return proj_planes

def stereo2xyz(projdir):
    """
    Convert 2D stereographic projection to 3D unit vectors.
    
    Inverse of stereoprojection_directions().
    
    Input:
        projdir: numpy array (2, N) - Projection coordinates [X, Y]
    
    Output:
        dirs: numpy array (3, N) - Unit direction vectors
    """
    r,p = 2.*np.arctan(np.sqrt(projdir[0]**2+projdir[1]**2)),np.arctan2(projdir[1],projdir[0])
    z = np.cos(r)
    xy = np.sqrt(1.-z**2)
    return xy*np.cos(p),xy*np.sin(p),z


def equalarea2xyz(projdir):
    """
    Convert 2D equal-area projection to 3D unit vectors.
    
    Inverse of equalarea_directions().
    
    Input:
        projdir: numpy array (2, N) - Projection coordinates
    
    Output:
        dirs: numpy array (3, N) - Unit direction vectors
    """
    if projdir[0]==0:
        an=np.pi/2
    else:
        an=np.arctan(projdir[1]/projdir[0])
    dirsxy1=np.cos(an)
    dirsxy2=np.sin(an)
    if projdir[0]==0:
        alpha=np.arcsin(projdir[1]/dirsxy2/2)*2
    else:
        alpha=np.arcsin(projdir[0]/dirsxy1/2)*2
    z=np.cos(alpha)
    if projdir[0]==0:
        x=0
        y=np.sin(alpha)
    elif projdir[1]==0:
        y=0
        x=np.sin(alpha)
    else:
        x=dirsxy1*np.sin(alpha)
        y=dirsxy2*np.sin(alpha)
    print(alpha*180/np.pi)
    return x,y,z

def equalarea_arr2xyz(projdir):
    """
    Convert 2D equal-area (Schmidt) projection coordinates to 3D unit vectors.
    
    Performs inverse stereographic projection from 2D projected points back to 
    3D unit vectors on the sphere. This is the inverse operation of 
    equalarea_directions().
    
    Mathematical Background:
    For equal-area (Lambert azimuthal) projection:
    - r = 2*sin(θ/2) where θ is the angle from projection pole
    - φ is the azimuthal angle in the projection plane
    The 3D unit vector is: (sin(θ)cos(φ), sin(θ)sin(φ), cos(θ))
    
    Input:
        projdir: numpy array (2, N) - Projection coordinates [x, y]
                 First row: x-coordinates in projection plane
                 Second row: y-coordinates in projection plane
                 Points should be within unit circle for upper hemisphere
    
    Output:
        dirs: numpy array (3, N) - Unit direction vectors [x, y, z]
              First row: x-components
              Second row: y-components  
              Third row: z-components
              All vectors have unit length: x² + y² + z² = 1
    
    Usage Example:
        >>> import numpy as np
        >>> from projlib import equalarea2xyz, equalarea_directions
        >>> 
        >>> # Example 1: Convert single projected point back to 3D
        >>> # Project [1,0,0] direction to 2D, then back to 3D
        >>> original_dir = np.array([[1], [0], [0]])
        >>> proj_coords = equalarea_directions(original_dir.T)  # Project to 2D
        >>> recovered_dir = equalarea2xyz(proj_coords)  # Back to 3D
        >>> 
        >>> print("Original:", original_dir.T)
        >>> print("Projected:", proj_coords)
        >>> print("Recovered:", recovered_dir)
        >>> # Should match original direction
        
        >>> # Example 2: Convert multiple projected points
        >>> # Create a grid of projected points
        >>> n_points = 100
        >>> theta = np.linspace(0, 2*np.pi, n_points)
        >>> r = 0.5  # Radius in projection plane
        >>> proj_x = r * np.cos(theta)
        >>> proj_y = r * np.sin(theta)
        >>> proj_coords = np.vstack([proj_x, proj_y])  # Shape: (2, 100)
        >>> 
        >>> # Convert to 3D unit vectors
        >>> dirs_3d = equalarea2xyz(proj_coords)  # Shape: (3, 100)
        >>> 
        >>> print(f"Input shape: {proj_coords.shape}")
        >>> print(f"Output shape: {dirs_3d.shape}")
        >>> # Verify unit length
        >>> lengths = np.sqrt(dirs_3d[0]**2 + dirs_3d[1]**2 + dirs_3d[2]**2)
        >>> print(f"All unit vectors: {np.allclose(lengths, 1.0)}")
        
        >>> # Example 3: Round-trip test for multiple directions
        >>> # Start with 3D directions
        >>> original_dirs = np.array([
        ...     [1, 0, 0],
        ...     [0, 1, 0],
        ...     [0, 0, 1],
        ...     [1, 1, 1],
        ...     [1, 1, 0]
        ... ]).T  # Shape: (3, 5)
        >>> 
        >>> # Normalize to unit vectors
        >>> original_dirs = original_dirs / np.linalg.norm(original_dirs, axis=0)
        >>> 
        >>> # Project to 2D
        >>> projected = equalarea_directions(original_dirs.T)  # Shape: (5, 2)
        >>> 
        >>> # Convert back to 3D
        >>> recovered = equalarea2xyz(projected.T)  # Shape: (3, 5)
        >>> 
        >>> # Check if we recovered original directions
        >>> print("Original vs Recovered:")
        >>> for i in range(5):
        ...     orig = original_dirs[:, i]
        ...     rec = recovered[:, i]
        ...     dot_product = np.dot(orig, rec)
        ...     print(f"  Direction {i}: dot product = {dot_product:.6f}")
        ...     # Should be close to 1.0 (parallel vectors)
        
        >>> # Example 4: Convert pole figure data
        >>> # Simulate scattered pole figure measurements
        >>> n_measurements = 500
        >>> # Random points in projection plane (within unit circle)
        >>> angle = np.random.rand(n_measurements) * 2 * np.pi
        >>> radius = np.sqrt(np.random.rand(n_measurements))  # Uniform on disk
        >>> proj_x = radius * np.cos(angle)
        >>> proj_y = radius * np.sin(angle)
        >>> proj_coords = np.vstack([proj_x, proj_y])
        >>> 
        >>> # Convert to 3D crystal directions
        >>> crystal_dirs = equalarea2xyz(proj_coords)
        >>> 
        >>> # Now can use these for texture analysis
        >>> print(f"Converted {n_measurements} pole figure points to 3D")
        >>> print(f"Z-components range: [{crystal_dirs[2].min():.3f}, "
        ...       f"{crystal_dirs[2].max():.3f}]")
        >>> # Z should be in [0, 1] for upper hemisphere
        
        >>> # Example 5: Interactive pole figure point selection
        >>> # User clicks on pole figure, get 3D direction
        >>> import matplotlib.pyplot as plt
        >>> from projlib import schmidtnet_half
        >>> 
        >>> # Setup pole figure
        >>> fig, ax = plt.subplots(figsize=(8, 8))
        >>> fig, ax = schmidtnet_half(ax=ax)
        >>> 
        >>> # Simulate click at position (0.3, 0.4)
        >>> click_x, click_y = 0.3, 0.4
        >>> click_coords = np.array([[click_x], [click_y]])
        >>> 
        >>> # Get corresponding 3D direction
        >>> direction_3d = equalarea2xyz(click_coords)
        >>> 
        >>> print(f"Clicked at ({click_x}, {click_y})")
        >>> print(f"Corresponds to direction: [{direction_3d[0,0]:.3f}, "
        ...       f"{direction_3d[1,0]:.3f}, {direction_3d[2,0]:.3f}]")
        >>> 
        >>> # Plot the point
        >>> ax.plot(click_x, click_y, 'ro', markersize=10)
        >>> ax.text(click_x, click_y+0.05, 
        ...         f'[{direction_3d[0,0]:.2f}, {direction_3d[1,0]:.2f}, '
        ...         f'{direction_3d[2,0]:.2f}]', ha='center')
    """
    # Handle input as 2D array (2, N)
    if projdir.ndim == 1:
        # Single point - reshape to (2, 1)
        projdir = projdir.reshape(2, 1)
    
    # Extract x and y coordinates
    x_proj = projdir[0, :]  # Shape: (N,)
    y_proj = projdir[1, :]  # Shape: (N,)
    
    # Calculate radial distance in projection plane
    r = np.sqrt(x_proj**2 + y_proj**2)
    
    # Calculate polar angle (angle from north pole)
    # For equal-area projection: r = 2*sin(theta/2)
    # Therefore: theta = 2*arcsin(r/2)
    theta = 2 * np.arcsin(r / 2)
    
    # Handle points at origin separately
    at_origin = (r < 1e-10)
    
    # Calculate azimuthal angle
    phi = np.arctan2(y_proj, x_proj)
    
    # Convert to 3D Cartesian coordinates
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    
    # Handle origin points (map to [0, 0, 1])
    x[at_origin] = 0
    y[at_origin] = 0
    z[at_origin] = 1
    
    # Stack into 3xN array
    dirs = np.vstack([x, y, z])
    
    return dirs

def equalarea_directions(dirs):

    """
    Project 3D direction vectors onto 2D equal-area projection plane (Schmidt net).
    
    Uses Lambert azimuthal equal-area projection: r = √2 × sin(θ/2), where θ is the
    angle from the north pole. Preserves area, making it ideal for texture analysis.
    Maximum radius is √2 for hemisphere projection.
    
    Input:
        dirs: numpy array (3, N) or (3,) - Direction vectors [x, y, z]
                                           Will be automatically normalized
                                           Single vector (3,) will be reshaped to (3, 1)
    
    Output:
        proj_dirs: numpy array (3, N) - Projected coordinates [X, Y, 0]
                                         X, Y are 2D equal-area projection coordinates
                                         Third row is zeros (kept for compatibility)
                                         Range: [-√2, √2] for full sphere
    
    Usage Example:
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
    """
    #dirs = [x1,x2,...,xn;y1,y2,...,yn;z1,z2,...,zn];
    #example: dirs = np.array([[0,1,2,3],[1,2,3,0],[0,3,2,1]])
    
    #normalize and project
       
    if len(dirs.shape)==1:
        dirs = np.expand_dims(dirs,axis=1)

    dirs = dirs.astype(float)
    #normalizing dirs
    dirs /= np.sqrt((dirs ** 2).sum(0))
    dirsxy = dirs[0:2,:];
    #print(dirsxy)
    eps=1e-6
    normdirsxy = np.sqrt((dirsxy ** 2).sum(0))
    idxs=np.where(normdirsxy<eps)
    normdirsxy[idxs]=1.
    dirsxy /= normdirsxy
    
    alpha=np.arccos(np.sign(dirs[2,:])*dirs[2,:]);
    proj_dirs = np.vstack((np.sin(alpha/2)*2.*dirsxy[0,:], 
               np.sin(alpha/2)*2.*dirsxy[1,:],
                np.zeros(dirs[2,:].shape)))
    proj_dirs[:,idxs]=0.
    #check
#    if False:
#        phi2=90*np.pi/180.
#        phi1=0*np.pi/180.
#        RotX = active_rotation(phi1, 'x') 
#        RotY = active_rotation(phi2, 'y')
#        Dc=[0.,0.,1.];
#        Ds = np.matmul(RotY,RotX).dot(Dc)
#        proj_Ds = equalarea_directions(Ds)  
#        #np.pi*proj_Ds[0][0]**2-(-np.cos(phi2)+np.cos(0))*2*np.pi
#        phi22=90*np.pi/180.
#        phi21=85*np.pi/180.
#        phi1=0*np.pi/180.
#        RotX = active_rotation(phi1, 'x') 
#        RotY = active_rotation(phi21, 'y')
#        RotY2 = active_rotation(phi22, 'y')
#        Dc=[0.,0.,1.];
#        Ds = np.matmul(RotY,RotX).dot(Dc)
#        proj_Ds = equalarea_directions(Ds)  
#        Ds = np.matmul(RotY2,RotX).dot(Dc)
#        proj_Ds2 = equalarea_directions(Ds)  
#        #np.pi*(proj_Ds2[0][0]**2-proj_Ds[0][0]**2)-(-np.cos(phi22)+np.cos(phi21))*2*np.pi
#        phi2=45*np.pi/180.
#        dphi2=5*np.pi/180.
#        dphi1=5.*np.pi/180.
#        RotY = active_rotation(phi2-dphi2, 'y')
#        RotY2 = active_rotation(phi2+dphi2, 'y')
#        Dc=[0.,0.,1.];
#        Ds = RotY.dot(Dc)
#        proj_Ds = equalarea_directions(Ds)  
#        Ds = RotY2.dot(Dc)
#        proj_Ds2 = equalarea_directions(Ds)  
#        #1./2.*dphi1*(proj_Ds2[0][0]**2-proj_Ds[0][0]**2)-(+np.cos(phi2-dphi2)-np.cos(phi2+dphi2))*dphi1
        
        
    return proj_dirs
# def equalarea_planes(normals,arclength=360.,iniangle=0.,hemisphere="both"):

    """
    Project plane traces onto equal-area (Schmidt) net.
    
    Input:
        normals: numpy array (3, N) - Plane normal vectors
        arclength: float - Arc length in degrees (default: 360)
        iniangle: float - Starting angle in degrees (default: 0)
        hemisphere: str - 'both', 'upper', or 'lower' (default: 'both')
    
    Output:
        traces: list of arrays - Trace points for each plane
    """
#     #%normals = [x1,x2,...,xn;y1,y2,...,yn;z1,z2,...,zn];
#     #%varargin{1} arclength in deg
#     #normals = np.transpose(np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1]]))
#     #
#     if len(normals.shape)==1:
#         normals = np.expand_dims(normals,axis=1)

#     normals = normals.astype(float)
#     normals /= np.sqrt((normals ** 2).sum(0))

#     proj_normals = equalarea_directions(normals)

#     idxs = np.where(abs(normals[0,:])+abs(normals[1,:])==0)[0]
    
#     inplanedirs = np.vstack((-normals[1,:],normals[0,:],np.zeros(normals[0,:].shape)));
#     inplanedirs[:,idxs] = np.vstack((np.zeros(normals[0,idxs].shape), -normals[2,idxs],normals[1,idxs]));
    
#     inplanedirs /= np.sqrt((inplanedirs ** 2).sum(0))

#     thirdaxis=np.cross(normals,inplanedirs,axisa=0,axisb=0,axisc=0)
# #    thirdaxis = np.vstack((normals[1,:]*inplanedirs[2,:]-normals[2,:]*inplanedirs[1,:],
# #                           -1*(normals[0,:]*inplanedirs[2,:]-normals[2,:]*inplanedirs[0,:]),
# #                           normals[0,:]*inplanedirs[1,:]-normals[1,:]*inplanedirs[0,:]));
#     t=np.linspace(iniangle,iniangle+arclength,180*2+1)*np.pi/180;
#     basicarc = np.vstack((np.cos(t),np.sin(t),np.zeros(t.shape)));
    
#     proj_planes=[];
#     Zdir=[]
#     for i in range(0,normals.shape[1]):
#         Rot2Global = np.transpose(np.vstack((inplanedirs[:,i],thirdaxis[:,i],normals[:,i])));
#         Ccp = np.matmul(Rot2Global,basicarc);
#         Zdir=Ccp[2]
#         if hemisphere == "both":            
#             Ds = equalarea_directions(Ccp)
#         else:
#             if hemisphere == "upper":
#                 idxs = np.where(Ccp[2,:]>=0)[0]
#                 Ds = equalarea_directions(Ccp[:,idxs])
#             elif hemisphere == "lower":
#                 idxs = np.where(Ccp[2,:]<=0)[0]
#                 Ds = equalarea_directions(Ccp[:,idxs])
#         proj_planes.append(Ds)
#     if len(proj_planes)==1:
#         return proj_planes[0]
#     else:          
#         return proj_planes


# def stereoprojection_planes(normals,arclength=360.,iniangle=0.,hemisphere="both"):
#     #%normals = [x1,x2,...,xn;y1,y2,...,yn;z1,z2,...,zn];
#     #%varargin{1} arclength in deg
#     #normals = np.transpose(np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1]]))
#     #
#     if len(normals.shape)==1:
#         normals = np.expand_dims(normals,axis=1)

#     normals = normals.astype(float)
#     normals /= np.sqrt((normals ** 2).sum(0))

#     proj_normals = stereoprojection_directions(normals)

#     idxs = np.where(abs(normals[0,:])+abs(normals[1,:])==0)[0]
    
#     inplanedirs = np.vstack((-normals[1,:],normals[0,:],np.zeros(normals[0,:].shape)));
#     inplanedirs[:,idxs] = np.vstack((np.zeros(normals[0,idxs].shape), -normals[2,idxs],normals[1,idxs]));
    
#     inplanedirs /= np.sqrt((inplanedirs ** 2).sum(0))

#     thirdaxis=np.cross(normals,inplanedirs,axisa=0,axisb=0,axisc=0)
# #    thirdaxis = np.vstack((normals[1,:]*inplanedirs[2,:]-normals[2,:]*inplanedirs[1,:],
# #                           -1*(normals[0,:]*inplanedirs[2,:]-normals[2,:]*inplanedirs[0,:]),
# #                           normals[0,:]*inplanedirs[1,:]-normals[1,:]*inplanedirs[0,:]));
#     t=np.linspace(iniangle,iniangle+arclength,180*2+1)*np.pi/180;
#     basicarc = np.vstack((np.cos(t),np.sin(t),np.zeros(t.shape)));
    
#     proj_planes=[];
#     Zdir=[]
#     points=[]
#     for i in range(0,normals.shape[1]):
#         Rot2Global = np.transpose(np.vstack((inplanedirs[:,i],thirdaxis[:,i],normals[:,i])));
#         Ccp = np.matmul(Rot2Global,basicarc);
#         points.append(Ccp)
#         Zdir=Ccp[2]
#         if hemisphere == "both":            
#             Ds = stereoprojection_directions(Ccp)
#         else:
#             if hemisphere == "upper":
#                 idxs = np.where(Ccp[2,:]>=0)[0]
#                 Ds = stereoprojection_directions(Ccp[:,idxs])
#             elif hemisphere == "lower":
#                 idxs = np.where(Ccp[2,:]<=0)[0]
#                 Ds = stereoprojection_directions(Ccp[:,idxs])
#         proj_planes.append(Ds)
#     if len(proj_planes)==1:
#         return proj_planes[0],points[0]
#     else:          
#         return proj_planes, points


def wulffnet(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.)):
    """
    Draw stereographic (Wulff) net - full circle.
    
    Input:
        ax: matplotlib axis - Existing axis (default: None)
        basedirs: bool - Plot base directions (default: False)
        facecolor: tuple - Background color RGB
    
    Output:
        fig, ax: matplotlib figure and axis
    """
    if ax==None:
        fig, ax = plt.subplots()
    else:
        fig=ax.get_figure()
    if basedirs:
        basicdirections = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1]]);
        basicdirections = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[1,1,0],[-1,1,0],[1,-1,0],[-1,-1,0],[1,1,1],[-1,1,1],[1,-1,1],[-1,-1,1],[0,1,1],[1,0,1]]);
        #basicdirections = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];
        basicdirectionstext = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1]]);
        basicdirectionstext = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[1,1,0],[-1,1,0],[1,-1,0],[-1,-1,0],[1,1,1],[-1,1,1],[1,-1,1],[-1,-1,1],[0,1,1],[1,0,1]]);
        #basicdirectionstext = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];



    #longitude lines
    #fig, ax = plt.subplots()
    ax.tick_params(
        axis='both',
        which='both',
        bottom=False,
        top=False,
        left=False,
        labelbottom=False,
        labelleft=False)
    ax.plot(0, 0, 'k+')
    circ = plt.Circle((0, 0), 1.0, facecolor=facecolor, edgecolor='black')
    ax.add_patch(circ)

    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
  # equal aspect ratio
    ax.axis('off')  # remove the box
    ##plt.show()
    
    t=np.linspace(0,180,180*2+1)*np.pi/180;
    xc = np.sin(t);
    yc = np.cos(t);
    AltitudeAngle = np.linspace(0,180,37)*np.pi/180;
    for an in AltitudeAngle:
        RotY = passive_rotation(an, 'y')
        Ccp = np.matmul(RotY,np.vstack((xc,yc,np.zeros(yc.shape))))

        proj_dirs = stereoprojection_directions(Ccp)
        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')

    #Latitude Lines
    LatitudeAngle = np.linspace(-90,90,37)*np.pi/180; 
    #t=np.linspace(0,180,360*2+1)*np.pi/180;
    zc = np.sin(t);
    xc = np.cos(t);
    for an in LatitudeAngle:#[0.]:#LatitudeAngle:
        Rmeridian = np.cos(an);
        px = Rmeridian*xc;
        py = np.sin(an)*np.ones(t.shape);
        pz = Rmeridian*zc;

        proj_dirs = stereoprojection_directions(np.vstack((px,py,pz)))

        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
    
    if basedirs:
        an=45.*np.pi/180.;
        an=0.*np.pi/180.;
        Rotz = active_rotation(an, 'z')
        an=np.arccos(1/np.sqrt(3));
        an=0.*np.pi/180.;
        Rotx = active_rotation(an, 'x')
        dirs = np.matmul(np.matmul(Rotx,Rotz),np.transpose(basicdirections))
        proj_dirs = stereoprojection_directions(dirs)
        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color='b',marker='o',linestyle='')
        for diri,proj_diri in zip(basicdirectionstext,np.transpose(proj_dirs)):
            ax.text(0.03+proj_diri[0],0.03+proj_diri[1],str(diri))

    ax.set_xlim((-1.05,1.05))
    ax.set_ylim((-1.05,1.05))

    return fig,ax


def wulffnet_half(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.)):
    """
    Draw stereographic (Wulff) net - upper hemisphere.
    
    Input:
        ax: matplotlib axis - Existing axis (default: None)
        basedirs: bool - Plot base directions (default: False)
        facecolor: tuple - Background color RGB
    
    Output:
        fig, ax: matplotlib figure and axis
    """
    if ax==None:
        fig, ax = plt.subplots()
    else:
        fig=ax.get_figure()
    if basedirs:
        basicdirections = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1],[1,0,2]]);
        #basicdirections = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];
        basicdirectionstext = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1],[1,0,2]]);
        #basicdirectionstext = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];



    #longitude lines
    #fig, ax = plt.subplots()
    ax.tick_params(
        axis='both',
        which='both',
        bottom=False,
        top=False,
        left=False,
        labelbottom=False,
        labelleft=False)
    ax.plot(0, 0, 'k+')
    
    
    
    
    w1 = Wedge((0,0), 1.0, 0, 180, fc=facecolor, edgecolor='black')
    ax.add_artist(w1)
    
            
    
    
#    circ = plt.Circle((0, 0), 1.0, facecolor=(210./255.,235./255.,255./255.), edgecolor='black')
#    ax.add_patch(circ)

    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
  # equal aspect ratio
    ax.axis('off')  # remove the box
   ##plt.show()
    
    t=np.linspace(0,90,180*2+1)*np.pi/180;
    xc = np.sin(t);
    yc = np.cos(t);
    AltitudeAngle = np.linspace(0,180,37)*np.pi/180;
    for an in AltitudeAngle:
        RotY = passive_rotation(an, 'y')
        Ccp = np.matmul(RotY,np.vstack((xc,yc,np.zeros(yc.shape))))

        proj_dirs = stereoprojection_directions(Ccp)
        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')

    #Latitude Lines
    LatitudeAngle = np.linspace(0,90,37)*np.pi/90; 
    #t=np.linspace(0,180,360*2+1)*np.pi/180;
    zc = np.sin(t);
    xc = np.cos(t);
    for an in LatitudeAngle:#[0.]:#LatitudeAngle:
        Rmeridian = np.cos(an);
        px = Rmeridian*xc;
        py = np.sin(an)*np.ones(t.shape);
        pz = Rmeridian*zc;

        proj_dirs = stereoprojection_directions(np.vstack((px,py,pz)))

        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
    
    if basedirs:
        an=45.*np.pi/180.;
        an=0.*np.pi/180.;
        Rotz = active_rotation(an, 'z')
        an=np.arccos(1/np.sqrt(3));
        an=0.*np.pi/180.;
        Rotx = active_rotation(an, 'x')
        dirs = np.matmul(np.matmul(Rotx,Rotz),np.transpose(basicdirections))
        proj_dirs = stereoprojection_directions(dirs)
        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color='b',marker='o',linestyle='')
        for diri,proj_diri in zip(basicdirectionstext,np.transpose(proj_dirs)):
            ax.text(0.03+proj_diri[0],0.03+proj_diri[1],str(diri))

    #ax.set_xlim((-1.05,1.05))
    #ax.set_ylim((-1.05,1.05))
    dd=0.05
    RR=1
    ax.set_xlim([-RR-dd,RR+dd])
    ax.set_ylim([0-dd,RR+dd])    


    return fig,ax

def wulffnet_quarter(ax=None,basedirs=False):
    """
    Draw quarter stereographic (Wulff) net.
    
    Input:
        ax: matplotlib axis - Existing axis (default: None)
        basedirs: bool - Plot base directions (default: False)
    
    Output:
        fig, ax: matplotlib figure and axis
    """
    if ax==None:
        fig, ax = plt.subplots()
    else:
        fig=ax.get_figure()
    if basedirs:
        basicdirections = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1],[1,0,2]]);
        #basicdirections = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];
        basicdirectionstext = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1],[1,0,2]]);
        #basicdirectionstext = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];



    #longitude lines
    #fig, ax = plt.subplots()
    ax.tick_params(
        axis='both',
        which='both',
        bottom=False,
        top=False,
        left=False,
        labelbottom=False,
        labelleft=False)
    ax.plot(0, 0, 'k+')
    
    
    
    
    w1 = Wedge((0,0), 1.0, 0, 90, fc=(210./255.,235./255.,255./255.), edgecolor='black')
    ax.add_artist(w1)
    
            
    
    
#    circ = plt.Circle((0, 0), 1.0, facecolor=(210./255.,235./255.,255./255.), edgecolor='black')
#    ax.add_patch(circ)

    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
  # equal aspect ratio
    ax.axis('off')  # remove the box
    #plt.show()
    
    t=np.linspace(0,90,180*2+1)*np.pi/180;
    xc = np.sin(t);
    yc = np.cos(t);
    AltitudeAngle = np.linspace(0,90,19)*np.pi/180;
    for an in AltitudeAngle:
        RotY = passive_rotation(an, 'y')
        Ccp = np.matmul(RotY,np.vstack((xc,yc,np.zeros(yc.shape))))

        proj_dirs = stereoprojection_directions(Ccp)
        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')

    #Latitude Lines
    LatitudeAngle = np.linspace(0,45,19)*np.pi/90; 
    t=np.linspace(0,90,180*2+1)*np.pi/180;
    zc = np.sin(t);
    xc = np.cos(t);
    for an in LatitudeAngle:#[0.]:#LatitudeAngle:
        Rmeridian = np.cos(an);
        px = Rmeridian*xc;
        py = np.sin(an)*np.ones(t.shape);
        pz = Rmeridian*zc;

        proj_dirs = stereoprojection_directions(np.vstack((px,py,pz)))

        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
    
    if basedirs:
        an=45.*np.pi/180.;
        an=0.*np.pi/180.;
        Rotz = active_rotation(an, 'z')
        an=np.arccos(1/np.sqrt(3));
        an=0.*np.pi/180.;
        Rotx = active_rotation(an, 'x')
        dirs = np.matmul(np.matmul(Rotx,Rotz),np.transpose(basicdirections))
        proj_dirs = stereoprojection_directions(dirs)
        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color='b',marker='o',linestyle='')
        for diri,proj_diri in zip(basicdirectionstext,np.transpose(proj_dirs)):
            ax.text(0.03+proj_diri[0],0.03+proj_diri[1],str(diri))


    return fig,ax


def schmidtnet(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.)):
    """
    Draw equal-area (Schmidt) net - full circle.
    
    Input:
        ax: matplotlib axis - Existing axis (default: None)
        basedirs: bool - Plot base directions (default: False)
        facecolor: tuple - Background color RGB
    
    Output:
        fig, ax: matplotlib figure and axis
    """
    if ax==None:
        fig, ax = plt.subplots()
    else:
        fig=ax.get_figure()

    if basedirs:
        basicdirections = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1],[1,0,2]]);
        #basicdirections = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];
        basicdirectionstext = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1],[1,0,2]]);
        #basicdirectionstext = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];



    #longitude lines
    #fig, ax = plt.subplots()
    ax.tick_params(
        axis='both',
        which='both',
        bottom=False,
        top=False,
        left=False,
        labelbottom=False,
        labelleft=False)
    ax.plot(0, 0, 'k+')
    equaarea_factor = 2./np.sqrt(2)
    circ = plt.Circle((0, 0), equaarea_factor*1.0, facecolor=facecolor, edgecolor='black')
    ax.add_patch(circ)

    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
    ax.axis('off')  # remove the box
    #plt.show()
    
    t=np.linspace(0,180,180*2+1)*np.pi/180;
    xc = np.sin(t);
    yc = np.cos(t);
    AltitudeAngle = np.linspace(0,180,37)*np.pi/180;
    for an in AltitudeAngle:
        RotY = passive_rotation(an, 'y')
        Ccp = np.matmul(RotY,np.vstack((xc,yc,np.zeros(yc.shape))))

        proj_dirs = equalarea_directions(Ccp)
        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')

    #Latitude Lines
    LatitudeAngle = np.linspace(-90,90,37)*np.pi/180; 
    #t=np.linspace(0,180,360*2+1)*np.pi/180;
    zc = np.sin(t);
    xc = np.cos(t);
    for an in LatitudeAngle:#[0.]:#LatitudeAngle:
        Rmeridian = np.cos(an);
        px = Rmeridian*xc;
        py = np.sin(an)*np.ones(t.shape);
        pz = Rmeridian*zc;

        proj_dirs = equalarea_directions(np.vstack((px,py,pz)))

        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
    
    if basedirs:
        an=45.*np.pi/180.;
        an=0.*np.pi/180.;
        Rotz = active_rotation(an, 'z')
        an=np.arccos(1/np.sqrt(3));
        an=0.*np.pi/180.;
        Rotx = active_rotation(an, 'x')
        dirs = np.matmul(np.matmul(Rotx,Rotz),np.transpose(basicdirections))
        proj_dirs = equalarea_directions(dirs)
        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color='b',marker='o',linestyle='')
        for diri,proj_diri in zip(basicdirectionstext,np.transpose(proj_dirs)):
            ax.text(0.03+proj_diri[0],0.03+proj_diri[1],str(diri))

    return fig,ax



def wulffnet_regular_grid(ax,dangle):
    """
    Create regular angular grid on Wulff net.
    
    Input:
        ax: matplotlib axis
        dangle: float - Angular spacing in degrees
        dirout: bool - Return directions if True (default: False)
        plot: bool - Plot grid if True (default: True)
    
    Output:
        If dirout=True: numpy array (3, N) - Grid directions
    """
    #dphi=10.deg
    #dtheta=10.deg
    #dangle = 10.
    
    Phi1=np.linspace(0.,360.-dangle,int(360./dangle))
    Phi2=np.linspace(0.,180.-dangle,int(180./dangle))
    GridX=[];
    GridY=[];
    Dc=[1.,0.,0.];
    for phi1 in Phi1:
        for phi2 in Phi2:        
            RotZ = active_rotation(phi1, 'z', deg=True) 
            RotY = active_rotation(phi2, 'y', deg=True)
            Ds = np.matmul(RotY,RotZ).dot(Dc)
            proj_Ds = stereoprojection_directions(Ds)
            GridX.append(proj_Ds[0,0])
            GridY.append(proj_Ds[1,0])
                
    

    #fig,ax = wulffnet()
    ax.plot(GridX,GridY,'.',color='r',markersize=1)
    
    return GridX,GridY

def schmidtnet_half(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.)):
    """
    Draw equal-area (Schmidt) net - upper hemisphere.
    
    Input:
        ax: matplotlib axis - Existing axis (default: None)
        basedirs: bool - Plot base directions (default: False)
        facecolor: tuple - Background color RGB
    
    Output:
        fig, ax: matplotlib figure and axis
    """
    if ax==None:
        fig, ax = plt.subplots()
    else:
        fig=ax.get_figure()
    if basedirs:
        basicdirections = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1],[1,0,2]]);
        #basicdirections = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];
        basicdirectionstext = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1],[1,0,2]]);
        #basicdirectionstext = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];



    #longitude lines
    #fig, ax = plt.subplots()
    ax.tick_params(
        axis='both',
        which='both',
        bottom=False,
        top=False,
        left=False,
        labelbottom=False,
        labelleft=False)
    ax.plot(0, 0, 'k+')
    equaarea_factor = 2./np.sqrt(2)

    w1 = Wedge((0,0), equaarea_factor*1.0, 0, 180, fc=facecolor, edgecolor='black')
    ax.add_artist(w1)

    #circ = plt.Circle((0, 0), equaarea_factor*1.0, facecolor=facecolor, edgecolor='black')
    #ax.add_patch(circ)

    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
    ax.axis('off')  # remove the box
    #plt.show()
    
    t=np.linspace(0,90,180*2+1)*np.pi/180;
    xc = np.sin(t);
    yc = np.cos(t);
    AltitudeAngle = np.linspace(0,180,37)*np.pi/180;
    for an in AltitudeAngle:
        RotY = passive_rotation(an, 'y')
        Ccp = np.matmul(RotY,np.vstack((xc,yc,np.zeros(yc.shape))))

        proj_dirs = equalarea_directions(Ccp)
        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')

    #Latitude Lines
    LatitudeAngle = np.linspace(0,180,37)*np.pi/180; 
    #t=np.linspace(0,180,360*2+1)*np.pi/180;
    zc = np.sin(t);
    xc = np.cos(t);
    for an in LatitudeAngle:#[0.]:#LatitudeAngle:
        Rmeridian = np.cos(an);
        px = Rmeridian*xc;
        py = np.sin(an)*np.ones(t.shape);
        pz = Rmeridian*zc;

        proj_dirs = equalarea_directions(np.vstack((px,py,pz)))

        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
    
    if basedirs:
        an=45.*np.pi/180.;
        an=0.*np.pi/180.;
        Rotz = active_rotation(an, 'z')
        an=np.arccos(1/np.sqrt(3));
        an=0.*np.pi/180.;
        Rotx = active_rotation(an, 'x')
        dirs = np.matmul(np.matmul(Rotx,Rotz),np.transpose(basicdirections))
        proj_dirs = equalarea_directions(dirs)
        ax.plot(proj_dirs[0,:],proj_dirs[1,:],color='b',marker='o',linestyle='')
        for diri,proj_diri in zip(basicdirectionstext,np.transpose(proj_dirs)):
            ax.text(0.03+proj_diri[0],0.03+proj_diri[1],str(diri))

    return fig,ax

def schmidt_regular_area_grid(ax,Na=72,Nr=20,plot=True):
    """
    Generate a regular grid on Schmidt (equal-area) stereographic projection.
    
    Creates a polar grid with radial and angular divisions that maintains approximately
    equal area for each grid cell. The grid density increases toward the center to
    compensate for the stereographic projection distortion.
    
    Input:
        ax: matplotlib axis - Axis object for plotting
        Na: int - Number of angular divisions at outermost radius (default: 72)
                  Controls angular resolution; should be divisible by 8
        Nr: int - Number of radial divisions (default: 20)
                  Controls radial resolution from center to edge
        plot: bool - If True, plot grid points on axis (default: True)
    
    Output:
        GridX: list - X-coordinates of grid points in projection plane
        GridY: list - Y-coordinates of grid points in projection plane
        GridR: list - Radial distances of grid points from center
        GridPhi: list - Angular positions of grid points in degrees [0, 360)
        AreaRatio: list - Area ratio of each grid cell relative to total hemisphere
                          Sum of all AreaRatio values equals 1.0
    
    Usage Example:
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> from projlib import schmidt_regular_grid, schmidtnet_half
        >>> 
        >>> # Create equal-area projection with regular grid
        >>> fig, ax = plt.subplots(figsize=(8, 8))
        >>> fig, ax = schmidtnet_half(ax=ax)
        >>> 
        >>> # Generate and plot grid
        >>> GridX, GridY, GridR, GridPhi, AreaRatio = schmidt_regular_grid(
        ...     ax, Na=72, Nr=20, plot=True
        ... )
        >>> 
        >>> print(f"Total grid points: {len(GridX)}")
        >>> print(f"Sum of area ratios: {sum(AreaRatio):.4f}")  # Should be ~1.0
        >>> plt.title('Schmidt Net Regular Grid')
        >>> plt.show()
        
        >>> # High-resolution grid for density calculations
        >>> fig, ax = plt.subplots(figsize=(10, 10))
        >>> fig, ax = schmidtnet_half(ax=ax)
        >>> GridX, GridY, GridR, GridPhi, AreaRatio = schmidt_regular_grid(
        ...     ax, Na=144, Nr=40, plot=True
        ... )
        >>> print(f"High-res grid points: {len(GridX)}")
        
        >>> # Use grid for texture analysis
        >>> # Calculate pole density at each grid point
        >>> fig, ax = plt.subplots()
        >>> fig, ax = schmidtnet_half(ax=ax)
        >>> GridX, GridY, GridR, GridPhi, AreaRatio = schmidt_regular_grid(
        ...     ax, Na=72, Nr=20, plot=False
        ... )
        >>> 
        >>> # Example: compute density from orientation data
        >>> # density[i] = number of poles near (GridX[i], GridY[i])
        >>> # weighted by AreaRatio[i]
        >>> densities = np.random.rand(len(GridX)) * 100  # Example data
        >>> scatter = ax.scatter(GridX, GridY, c=densities, s=20, 
        ...                      cmap='jet', vmin=0, vmax=100)
        >>> plt.colorbar(scatter, ax=ax, label='Density (MUD)')
        >>> plt.title('Pole Figure Density Map')
        
        >>> # Coarse grid for quick visualization
        >>> fig, ax = plt.subplots(figsize=(6, 6))
        >>> fig, ax = schmidtnet_half(ax=ax)
        >>> GridX, GridY, GridR, GridPhi, AreaRatio = schmidt_regular_grid(
        ...     ax, Na=36, Nr=10, plot=True
        ... )
        >>> # Mark grid cells with their area ratios
        >>> for i, (x, y, area) in enumerate(zip(GridX, GridY, AreaRatio)):
        ...     if i % 50 == 0:  # Label every 50th point
        ...         ax.text(x, y, f'{area:.4f}', fontsize=6, ha='center')
        >>> plt.title('Grid with Area Ratios')
        
        >>> # Extract grid points in specific angular sector
        >>> GridX, GridY, GridR, GridPhi, AreaRatio = schmidt_regular_grid(
        ...     ax, Na=72, Nr=20, plot=False
        ... )
        >>> # Select points in 0-90 degree sector
        >>> mask = (np.array(GridPhi) >= 0) & (np.array(GridPhi) <= 90)
        >>> sector_points = sum(mask)
        >>> sector_area = sum([a for a, m in zip(AreaRatio, mask) if m])
        >>> print(f"Points in first quadrant: {sector_points}")
        >>> print(f"Area fraction: {sector_area:.4f}")  # Should be ~0.25
    """
    dphi1=360/Na
    phi1=np.linspace(0,360-dphi1,Na)
    R=equalarea_directions(np.array([1,0,0]))[0,0]
    TotalArea=np.pi*R**2
    dr=R/(Nr+0.5)
    r=np.linspace(0,R-dr/2,Nr+1)
#    GridX=[0.];
#    GridY=[0.];
#    Weight= np.pi*(r[1]/2)**2/TotalArea
    AreaRatio=[]
    for ri in r:
        Nari=int(Na*(ri/r[-1]))
        
        if Nari<8:
            Nari=8
        else:
            Nari=Nari-Nari%8    
        #print(Nari)
        #Nari=8
        dphi1=360./(Nari)
        phi1=np.linspace(0,360-dphi1,Nari)
        phi1=np.linspace(dphi1/2,360-dphi1/2,Nari)
        #phi1=phi1[0:-1]
        #print(2*np.pi*((ri+dr/2)**2/2-(ri-dr/2)**2/2))
        #tot=0
        for phi1i in phi1:        
            if ri==0:
                GridX=[0.];
                GridPhi=[0.]
                GridY=[0.];
                GridR=[r[1]-dr/2]
                GridR=[ri]
                AreaRatio= [np.pi*(r[1]-dr/2.)**2/TotalArea]
            else:
                GridX.append(ri*np.cos(phi1i*np.pi/180.))
                GridY.append(ri*np.sin(phi1i*np.pi/180.))
                AreaRatio.append(dphi1*np.pi/180.*((ri+dr/2)**2/2.-(ri-dr/2)**2/2.)/TotalArea)
                GridR.append(ri)
                phi=np.arctan2(GridY[-1],GridX[-1])*180./np.pi
                GridPhi.append(phi)
                #tot+=dphi1*np.pi/180.*((ri+dr/2)**2/2-(ri-dr/2)**2/2)
        #print(tot)
    #sum(AreaRatio)
    if plot:
        ax.plot(GridX,GridY,'.',color='r',markersize=1)
    
    return GridX,GridY,GridR,GridPhi,AreaRatio, equalarea_arr2xyz(np.vstack((GridX,GridY)))


def schmidt_regular_grid(ax,Na=72,Nr=20,plot=True):
    """
    Create regular grid on Schmidt net.
    
    Input:
        ax: matplotlib axis
        Na: int - Azimuthal divisions (default: 72)
        Nr: int - Radial divisions (default: 20)
        plot: bool - Plot grid (default: True)
    
    Output:
        grid: tuple - Grid parameters
    """
    dphi1=360/Na
    phi1=np.linspace(0,360-dphi1,int(Na))
    R=equalarea_directions(np.array([1,0,0]))[0,0]
    TotalArea=np.pi*R**2
    dr=R/(Nr+0.5)
    r=np.linspace(0,R-dr/2,int(Nr+1))
#    GridX=[0.];
#    GridY=[0.];
#    Weight= np.pi*(r[1]/2)**2/TotalArea
    AreaRatio=[]
    for ri in r:
        Nari=int(Na*(ri/r[-1]))
        
        if Nari<8:
            Nari=8
        else:
            Nari=Nari-Nari%8    
        #print(Nari)
        #Nari=8
        dphi1=360./(Nari)
        phi1=np.linspace(0,360-dphi1,Nari)
        phi1=np.linspace(dphi1/2,360-dphi1/2,Nari)
        #phi1=phi1[0:-1]
        #print(2*np.pi*((ri+dr/2)**2/2-(ri-dr/2)**2/2))
        #tot=0
        for phi1i in phi1:        
            if ri==0:
                GridX=[0.];
                GridPhi=[0.]
                GridY=[0.];
                GridR=[r[1]-dr/2]
                GridR=[ri]
                AreaRatio= [np.pi*(r[1]-dr/2.)**2/TotalArea]
            else:
                GridX.append(ri*np.cos(phi1i*np.pi/180.))
                GridY.append(ri*np.sin(phi1i*np.pi/180.))
                AreaRatio.append(dphi1*np.pi/180.*((ri+dr/2)**2/2.-(ri-dr/2)**2/2.)/TotalArea)
                GridR.append(ri)
                phi=np.arctan2(GridY[-1],GridX[-1])*180./np.pi
                GridPhi.append(phi)
                #tot+=dphi1*np.pi/180.*((ri+dr/2)**2/2-(ri-dr/2)**2/2)
        #print(tot)
    #sum(AreaRatio)
    if plot:
        ax.plot(GridX,GridY,'.',color='r',markersize=1)
    
    return GridX,GridY,GridR,GridPhi,AreaRatio


def pf_cmap02(GridX,GridY,GridR,GridPhi,AreaRatio,Intensity,NoCont=10,GridSize=1000,cmap='jet',method='cubic'):
    
    
    """
    Create contour pole figure with colormap (internal function).
    
    Input:
        GridX, GridY, GridR, GridPhi, AreaRatio, Intensity: Grid and density data
        NoCont, GridSize, cmap, method: Plotting parameters
    
    Output:
        Contour plot object
    """
#    r = [np.sqrt(x**2 + y**2) for x,y in zip(GridX,GridY)]
#    theta = [np.arctan2(y,x) for x,y in zip(GridX,GridY)]
#    theta = theta-min(theta)
#    equaarea_factor = 2./np.sqrt(2)
#    
#    ri = np.linspace(0,equaarea_factor,100)
#    thetai = np.linspace(0,2.0*np.pi,100)
#    thetai,ri = np.mgrid[0:2*np.pi:GridSize, 0:equaarea_factor:GridSize]
    
    equaarea_factor = 2./np.sqrt(2)
    xi,yi = np.mgrid[-equaarea_factor:equaarea_factor:GridSize, -equaarea_factor:equaarea_factor:GridSize]
    Ii = scipy.interpolate.griddata((GridX,GridY),Intensity,(xi,yi),method=method)
    
    
    #Ii=np.nan_to_num(Ii)
    fig = plt.figure()
    #fig,ax = schmidtnet(basedirs=False,facecolor='white')
    ax = fig.add_subplot(111)
    
    palette = plt.cm.jet
    palette.set_bad ('w',1.0) # Bad values (i.e., masked, set to grey 0.8
    A = np.ma.array ( Ii, mask=np.isnan(Ii))
    
    CS=plt.contourf(xi,yi,A,NoCont,cmap=cmap)
    theta = np.linspace(0,2.0*np.pi,100)
    xc=equaarea_factor*np.cos(theta)
    yc=equaarea_factor*np.sin(theta)
    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
  # equal aspect ratio

    plt.plot(xc,yc, color='black',linewidth=2)
    
    
    
    
    
    
    ax.axis('off')  # remove the box

#    cb = plt.colorbar(p2, cax=ax)
    cb=plt.colorbar(CS)
    ax.set_xlim([-equaarea_factor,1*equaarea_factor])
    ax.set_ylim([-equaarea_factor,1*equaarea_factor])
    #plt.show()
    cb.remove()
    cmin=np.nanmin(Ii)
    cmax=np.nanmax(Ii)
    
    ax2= fig.add_axes(ax.get_position())
    ax2.set_position([0.85, 0.1,0.03, 0.8])
    #ax2.axis('off')  # remove the box
    
    norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
    cmap = mpl.cm.ScalarMappable(
          norm = norm, 
          cmap = cmap)
    cmap.set_array([])
    cb=fig.colorbar(cmap,cax=ax2)
    cb.ax.set_ylabel('Multiple of random distribution')

    textvar1=ax.text(0., equaarea_factor, "Y", size=15, rotation=0.,
             ha="center", va="center",
             bbox=dict(boxstyle="round",
                       ec=(0., 0, 0),
                       fc=(1., 1, 1),
                       ))
    textvar2=ax.text(equaarea_factor, 0.,"X", size=15, rotation=0.,
             ha="center", va="center",
             bbox=dict(boxstyle="round",
                       ec=(0., 0, 0),
                       fc=(1., 1, 1),
                       ))

    plt.tight_layout()
    #plt.show()
    return fig,ax,ax2,cb,xi,yi,A,CS,cmap
    



def pf(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20,syms=True,s=50,facecolor=(210./255.,235./255.,255./255.),plot=True):
    """
    Generate pole figure from Euler angles.
    
    Input:
        gPhi1, gPHI, gPhi2: arrays - Bunge Euler angles (degrees)
        Dc: array (3,) - Crystal direction [uvw]
        lattice: str - Crystal system
        **kwargs: Additional parameters
    
    Output:
        fig, ax: matplotlib figure and axis
    """
    #fig,ax = schmidtnet()
    GridX,GridY,GridR,GridPhi,AreaRatio=schmidt_regular_grid([],Na=Na,Nr=Nr,plot=False)
    
    Intensity = np.array(GridX)*0
    inc=0
    #Dc=[1,0,0]
    Symsall = symmetry_elements(lattice)
    if syms:
        Syms=Symsall
    
    for Phi1,PHI,Phi2 in zip(gPhi1,gPHI,gPhi2):
        inc+=1
        print(str(inc)+'/'+str(len(gPhi1)))
        #R = np.array(np_euler_matrix(Phi1, PHI,Phi2))
        #Phi1,PHI,Phi2 = euler_angles_reduction(Phi1,PHI,Phi2)
        U = np_inverse_euler_matrix(Phi1, PHI,Phi2)
        if not syms:
            Syms = [Symsall[random.randint(0,len(Symsall)-1)]]
        for Sym in Syms:
            #Ri=np.matmul(Sym,R)
            #Phi1,PHI,Phi2=euler_angles_from_matrix(Ri)
            #Phi1,PHI,Phi2 = euler_angles_reduction(Phi1,PHI,Phi2)
            Ds = np.array(U).dot(Sym.dot(Dc))
            proj_Ds = equalarea_directions(Ds)
            phi = np.arctan2(proj_Ds[1],proj_Ds[0])[0]*180./np.pi
            r=np.sqrt(proj_Ds[:,0].dot(proj_Ds[:,0]))
    
            dr=np.array(GridR)-r
            idxmin=np.where(abs(dr)==min(abs(dr)))
            ri=GridR[idxmin[0][0]]
    
            idxr=np.where(np.array(GridR)==ri)[0]
            dtan = np.array(GridPhi)[idxr]-phi
            idxtan = np.where(abs(dtan)==min(abs(dtan)))
            idxmin3=idxr[idxtan[0]][0]
            
    #        dx = np.array(GridX)-proj_Ds[0];
    #        dy = np.array(GridY)-proj_Ds[1];
    #        dr = dx**2+dy**2
    #        idxmin = np.where(dr==min(dr))[0][0]
            #Intensity[idxmin3]+=(1./AreaRatio[idxmin3])
            Intensity[idxmin3]+=(1./(len(Syms)*len(gPhi1)*AreaRatio[idxmin3]))
#            Intensity[idxmin3]+=1./len(gPhi1)
    #        GridX[idxmin]
    #        GridY[idxmin]
    
    
#    Intensity2 = np.array(GridX)*0
#    for Int,Int2,AR in zip()   
    #fig, ax = plt.subplots()
    if plot:
        fig,ax = schmidtnet(basedirs=False,facecolor=facecolor)
        GridX,GridY,GridR,GridPhi,AreaRatio=schmidt_regular_grid(ax,Na=Na,Nr=Nr,plot=True)
    
        plt.scatter(GridX,GridY, c=Intensity, s=s, edgecolor='',zorder=10,cmap='jet')
        cb=plt.colorbar()
        #plt.show()
        return fig,ax,cb,Intensity
    else:
        return GridX,GridY,GridR,GridPhi,AreaRatio,Intensity





def pf_cmap(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20,syms=True,s=50,NoCont=10,GridSize=1000,cmap='jet',method='cubic'):
    
    """
    Generate pole figure with density colormap.
    
    Input:
        gPhi1, gPHI, gPhi2: arrays - Bunge Euler angles (degrees)
        Dc: array (3,) - Crystal direction [uvw]
        lattice: str - Crystal system
        **kwargs: Additional parameters
    
    Output:
        fig, ax: matplotlib figure and axis
    """
    GridX,GridY,GridR,GridPhi,AreaRatio,Intensity=pf(gPhi1,gPHI,gPhi2,Dc,lattice,Na=Na,Nr=Nr,syms=syms,s=s,facecolor="None",plot=False)
    
    r = [np.sqrt(x**2 + y**2) for x,y in zip(GridX,GridY)]
    theta = [np.arctan2(y,x) for x,y in zip(GridX,GridY)]
    theta = theta-min(theta)
    equaarea_factor = 2./np.sqrt(2)
    
    ri = np.linspace(0,equaarea_factor,100)
    thetai = np.linspace(0,2.0*np.pi,100)
    thetai,ri = np.mgrid[0:2*np.pi:GridSize, 0:equaarea_factor:GridSize]
    # grid the data.
    Ii = scipy.interpolate.griddata((np.concatenate((theta-2*np.pi,theta,theta+2*np.pi), axis=0),np.concatenate((r,r,r), axis=0)),
                                    np.concatenate((Intensity,Intensity,Intensity), axis=0),(thetai,ri),method=method)
    
    
    Ii=np.nan_to_num(Ii)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="polar")
    
        
    CS=plt.contourf(thetai,ri,Ii,NoCont,cmap=cmap)
    plt.plot(np.linspace(0,2.0*np.pi,1000),np.linspace(0,2.0*np.pi,1000)*0+equaarea_factor, color='black',linewidth=2)
    ax.axis('off')  # remove the box

    #cb = plt.colorbar(p2, cax=ax)
    cb=plt.colorbar(CS)
    ax.set_xlim([0,2*np.pi])
    ax.set_ylim([0,1*equaarea_factor])
    #plt.show()
    cb.remove()
    cmin=np.nanmin(Ii)
    cmax=np.nanmax(Ii)
    
    ax2= fig.add_axes(ax.get_position())
    ax2.set_position([0.85, 0.1,0.03, 0.8])
    #ax2.axis('off')  # remove the box
    
    norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
    cmap = mpl.cm.ScalarMappable(
          norm = norm, 
          cmap = cmap)
    cmap.set_array([])
    cb=fig.colorbar(cmap,cax=ax2)
    cb.ax.set_ylabel('Multiple of random distribution')

    textvar1=ax.text(0., 2./np.sqrt(2), "X", size=15, rotation=0.,
             ha="center", va="center",
             bbox=dict(boxstyle="round",
                       ec=(0., 0, 0),
                       fc=(1., 1, 1),
                       ))
    textvar2=ax.text(np.pi/2, 2./np.sqrt(2), "Y", size=15, rotation=0.,
             ha="center", va="center",
             bbox=dict(boxstyle="round",
                       ec=(0., 0, 0),
                       fc=(1., 1, 1),
                       ))

    plt.tight_layout()
    #plt.show()
    return fig,ax,ax2,cb,thetai,ri,Ii,CS,cmap
    


def pf_cmap_cscale(fig,ax2,cmin,cmax,cmap):
    """
    Add colorbar to pole figure.
    
    Input:
        fig, ax2, cmin, cmax, cmap: Figure and colorbar parameters
    
    Output:
        None (modifies figure)
    """
    #cmin=0
    #cmax=5
    plt.clim(cmin,cmax)
    
    cmap = mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax),
          cmap = cmap)
    cmap.set_array([])
    cb=fig.colorbar(cmap,cax=ax2)
    cb.ax.set_ylabel('Multiple of random distribution')


    

def ipf(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20,syms=True):
    """
    Generate inverse pole figure from Euler angles.
    
    Input:
        gPhi1, gPHI, gPhi2: arrays - Bunge Euler angles (degrees)
        Dc: array (3,) - Sample direction [xyz]
        lattice: str - Crystal system
        **kwargs: Additional parameters
    
    Output:
        fig, ax: matplotlib figure and axis
    """
    #fig,ax = schmidtnet()
    GridX,GridY,GridR,GridPhi,AreaRatio=schmidt_regular_grid([],Na=Na,Nr=Nr,plot=False)
    
    Intensity = np.array(GridX)*0
    inc=0
    #Dc=[1,0,0]
    Symsall = symmetry_elements(lattice)
    if syms:
        Syms=Symsall
    
    for Phi1,PHI,Phi2 in zip(gPhi1,gPHI,gPhi2):
        inc+=1
        print(str(inc)+'/'+str(len(gPhi1)))
        #R = np.array(np_euler_matrix(Phi1, PHI,Phi2))
        #Phi1,PHI,Phi2 = euler_angles_reduction(Phi1,PHI,Phi2)
        U = np_euler_matrix(Phi1, PHI,Phi2)
        if not syms:
            Syms = [Symsall[random.randint(0,len(Symsall)-1)]]
        for Sym in Syms:
            #Ri=np.matmul(Sym,R)
            #Phi1,PHI,Phi2=euler_angles_from_matrix(Ri)
            #Phi1,PHI,Phi2 = euler_angles_reduction(Phi1,PHI,Phi2)
            Ds = np.array(U).dot(Sym.dot(Dc))
            proj_Ds = equalarea_directions(Ds)
            phi = np.arctan2(proj_Ds[1],proj_Ds[0])[0]*180./np.pi
            r=np.sqrt(proj_Ds[:,0].dot(proj_Ds[:,0]))
    
            dr=np.array(GridR)-r
            idxmin=np.where(abs(dr)==min(abs(dr)))
            ri=GridR[idxmin[0][0]]
    
            idxr=np.where(np.array(GridR)==ri)[0]
            dtan = np.array(GridPhi)[idxr]-phi
            idxtan = np.where(abs(dtan)==min(abs(dtan)))
            idxmin3=idxr[idxtan[0]][0]
            
    #        dx = np.array(GridX)-proj_Ds[0];
    #        dy = np.array(GridY)-proj_Ds[1];
    #        dr = dx**2+dy**2
    #        idxmin = np.where(dr==min(dr))[0][0]
            #Intensity[idxmin3]+=(1./AreaRatio[idxmin3])
            Intensity[idxmin3]+=(1./(len(Syms)*len(gPhi1)*AreaRatio[idxmin3]))
#            Intensity[idxmin3]+=1./len(gPhi1)
    #        GridX[idxmin]
    #        GridY[idxmin]
    
    
#    Intensity2 = np.array(GridX)*0
#    for Int,Int2,AR in zip()   
    #fig, ax = plt.subplots()
    fig,ax = schmidtnet()
    GridX,GridY,GridR,GridPhi,AreaRatio=schmidt_regular_grid(ax,Na=Na,Nr=Nr,plot=True)

    plt.scatter(GridX,GridY, c=Intensity, s=50, edgecolor='',zorder=10,cmap='jet')
    cb=plt.colorbar()
    #plt.show()
    return fig,ax,cb,Intensity
        

def stereotriangle(ax=None,basedirs=False,equalarea=False,grid=False,resolution=None,gridmarkersize=None,gridmarkercol=None,gridzorder=None,mesh=False):
    """
    Draw standard stereographic triangle for cubic system.
    
    Input:
        ax: matplotlib axis - Existing axis (default: None)
        basedirs: bool - Label corners (default: False)
        equalarea: bool - Use equal-area (default: False)
        grid: bool - Draw grid (default: False)
        **kwargs: Additional parameters
    
    Output:
        fig, ax: matplotlib figure and axis
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig=ax.get_figure()
        
    ax.tick_params(
        axis='both',
        which='both',
        bottom=False,
        top=False,
        left=False,
        labelbottom=False,
        labelleft=False)
    
    if resolution is None:
        resolution=5
    if gridmarkersize is None:
        gridmarkersize=5
    if gridmarkercol is None:
        gridmarkercol='k'
    if gridzorder is None:
        gridzorder=50000
    
    ax.plot(0, 0, 'k+')
    #circ = plt.Circle((0, 0), 1.0, facecolor=(210./255.,235./255.,255./255.), edgecolor='black')
    #ax.add_patch(circ)

    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
    ax.axis('off')  # remove the box
    #plt.show()


    #print(max(proj_normals[1,:])/max(proj_normals[0,:]))
    #testing = generate grid of crystallographic direction = integer Miller indexes
    if False:
        an1=np.arctan(1/np.array(list(range(1,8))))
        for ii,a1 in enumerate(an1):
            vertical=[]
            an2=np.arctan(np.linspace(0,1,ii+2)*np.tan(a1))
            an2max=an2[-1]
            an22=np.arctan(1/np.linspace(1,15-ii,15-ii))[::-1]
            an22=np.append([0],an22)
            an22=an22[an22<an2max]
            an22=np.append(an22,an2max)
            for a2 in an22:                
                vertical.append(np.array([np.tan(a1),np.tan(a2),1]))
        
            vertical=np.array(vertical).T
            if equalarea:
                proj_dirs=equalarea_directions(vertical)   
            else:
                proj_dirs=stereoprojection_directions(vertical)   
            ax.plot(proj_dirs[0,:],proj_dirs[1,:],color='r',markeredgecolor='None',marker="o",markersize=10,zorder=50000)

    if mesh:
        dan=resolution
        dan2=resolution
        an1=[dan*(ii+1)/180*np.pi for ii in range(int(45/dan))]
        for ii,a1 in enumerate(an1):
            vertical=[]
            an2=[dan2*(jj)/180*np.pi for jj in range(int((2+ii)*dan/dan2))]
            for a2 in an2:
                vertical.append(np.array([np.tan(a1),np.tan(a2),1]))
            vertical=np.array(vertical).T
            if equalarea:
                proj_dirs=equalarea_directions(vertical)   
            else:
                proj_dirs=stereoprojection_directions(vertical)   
            ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--',zorder=gridzorder)
        
        for jj in range(int(45/dan)+1):
            horizontal=[]
            for ii in range(int(45/dan)-jj+1):
                horizontal.append(np.array([np.tan(dan*(ii+jj)/180*np.pi),np.tan(dan*(jj)/180*np.pi),1]))
            horizontal=np.array(horizontal).T
            if equalarea:
                proj_dirs=equalarea_directions(horizontal)   
            else:
                proj_dirs=stereoprojection_directions(horizontal)   
            ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--',zorder=gridzorder)
        
     
    normals = np.array([1,-1,0]);
#    normals = np.array([1,1,1]);
    if equalarea:
        proj_111=equalarea_directions(np.array([1,1,1]))
    else:
        proj_111=stereoprojection_directions(np.array([1,1,1]))
    arclength = 55#-np.arccos(np.sqrt(2)/np.sqrt(3))*180/np.pi;
    
    if equalarea:
        proj_normals = equalarea_planes(normals,arclength=arclength,iniangle=35)
    else:
        proj_normals = stereoprojection_planes(normals,arclength=arclength,iniangle=35)
        
    ax.plot(proj_normals[0,:], proj_normals[1,:], 'k')
        
    
    
    normals = np.array([-1,0,1]);
    arclength = 35#-np.arccos(np.sqrt(2)/np.sqrt(3))*180/np.pi;
    if equalarea:
        proj_normals = equalarea_planes(normals,arclength=arclength,iniangle=90)
    else:
        proj_normals = stereoprojection_planes(normals,arclength=arclength,iniangle=90)
    ax.plot(proj_normals[0,:], proj_normals[1,:], 'k')
    
    if equalarea:
        R = equalarea_directions(np.array([1,0,1]))[0]
    else:
        R = stereoprojection_directions(np.array([1,0,1]))[0]
    
    ax.plot([0,R[0]],[0,0], 'k')
#    dirs = np.column_stack([[0,0,1],[1,1,0],[1,0,0]]);

    dirs = np.column_stack([[0,0,1],[1,1,1],[1,0,1]]);
    
    
    if equalarea:
        proj_dirs = equalarea_directions(dirs)
    else:
        proj_dirs = stereoprojection_directions(dirs)
    
    ax.plot(proj_dirs[0,:], proj_dirs[1,:], 'ko',zorder=50)
    if basedirs:
        #ax.plot(proj_dirs[0,:], proj_dirs[1,:], 'ko',zorder=50)
        for diri,proj_diri in zip(np.transpose(dirs),np.transpose(proj_dirs)):
            if (diri==[0,0,1]).all():
                ax.text(-0.1+proj_diri[0],0.01+proj_diri[1],str(diri))
            else:
                ax.text(0.01+proj_diri[0],0.01+proj_diri[1],str(diri))
    
    if grid:
        if resolution is None:
            resolution=5
        if gridmarkersize is None:
            gridmarkersize=5
        if gridmarkercol is None:
            gridmarkercol='k'
        if gridzorder is None:
            gridzorder=50000

        grid_cub = get_beam_directions_grid("cubic", resolution, mesh="spherified_cube_edge")
        grid_stereo = Rotation.from_euler(np.deg2rad(grid_cub))*Vector3d.zvector()
    
        if equalarea:
            proj_Ds = equalarea_directions(grid_stereo.data.T)
        else:
            proj_Ds = stereoprojection_directions(grid_stereo.data.T)
        
        # #Colors=stereotriangle_colors(proj_Ds)
        # #fig,ax=stereotriangle(ax=None,basedirs=basedirs)
        ax.scatter(proj_Ds[0,:],proj_Ds[1,:],c=gridmarkercol,s=gridmarkersize,zorder=gridzorder)#,'.',color='r',markersize=1)
        ax.scatter([0],[0],c=gridmarkercol,s=gridmarkersize,zorder=gridzorder)#,'.',color='r',markersize=1)
    
    dirs= np.column_stack([[1,1,1],[1,0,1]]);
    if equalarea:
        proj_dirs = equalarea_directions(dirs)
    else:
        proj_dirs = stereoprojection_directions(dirs)
    
    dd=0.02
    #ax.set_xlim([0-dd,proj_dirs[0,1]+dd])
    #ax.set_ylim([0-dd,proj_dirs[1,0]+dd])    
    
    return fig,ax


def colored_stereotriangle(basedirs=False,resolution = 1, markersize=1):
    """
    Draw stereographic triangle with IPF coloring.
    
    Input:
        basedirs: bool - Label corners (default: False)
        resolution: float - Grid resolution (default: 1)
        markersize: float - Point size (default: 1)
    
    Output:
        fig, ax: matplotlib figure and axis
    """
#resolution = 1 
    grid_cub = get_beam_directions_grid("cubic", resolution, mesh="spherified_cube_edge")
    grid_stereo = Rotation.from_euler(np.deg2rad(grid_cub))*Vector3d.zvector()

    proj_Ds = stereoprojection_directions(grid_stereo.data.T)
    Colors=stereotriangle_colors(proj_Ds)
    fig,ax=stereotriangle(ax=None,basedirs=basedirs)
    ax.scatter(proj_Ds[0,:],proj_Ds[1,:],c=Colors,s=markersize)#,'.',color='r',markersize=1)
    #plt.show()
    return fig,ax

def filled_colored_stereotriangle(basedirs=False,resolution = 1, markersize=1,ax=None,**kwargs):
    """
    Draw filled triangle with smooth IPF coloring.
    
    Input:
        basedirs, resolution, markersize, ax, **kwargs
    
    Output:
        fig, ax: matplotlib figure and axis
    """
#resolution = 1 
    #import matplotlib.tri as tri
    
    grid_cub = get_beam_directions_grid("cubic", resolution, mesh="spherified_cube_edge")
    grid_stereo = Rotation.from_euler(np.deg2rad(grid_cub))*Vector3d.zvector()

    proj_Ds = stereoprojection_directions(grid_stereo.data.T)
    Colors=stereotriangle_colors(proj_Ds)
    cmap=ListedColormap(Colors)
    fig,ax=stereotriangle(ax=ax,basedirs=basedirs)
    #ax.tricontourf(proj_Ds[0,:],proj_Ds[1,:], range(0,len(Z)),vmin=0,vmax=len(Z),cmap=ListedColormap(Colors))
    ax.tripcolor(proj_Ds[0,:],proj_Ds[1,:], list(range(0,Colors.shape[0])),cmap=cmap, shading='gouraud',**kwargs)
    #ax.scatter(proj_Ds[0,:],proj_Ds[1,:],c=Colors,s=markersize)#,'.',color='r',markersize=1)  
    
    
    return fig,ax

def colors4stereotriangle(resolution = 1):
    """
    Generate RGB colors for IPF coloring.
    
    Input:
        resolution: float - Angular resolution (default: 1)
    
    Output:
        proj, RGB: Projection coordinates and RGB colors
    """
#resolution = 1 
    grid_cub = get_beam_directions_grid("cubic", resolution, mesh="spherified_cube_edge")
    grid_stereo = Rotation.from_euler(np.deg2rad(grid_cub))*Vector3d.zvector()

    proj_Ds = stereoprojection_directions(grid_stereo.data.T)
    Colors=stereotriangle_colors(proj_Ds)
    return proj_Ds,Colors


def stereotriangle_colors_from_d_IPF(d_IPF):
    """
    Get IPF colors for crystal directions.
    
    Input:
        d_IPF: numpy array (3, N) - Directions in crystal frame
    
    Output:
        RGB: numpy array (N, 3) - RGB colors
    """
    #d - list of ipf directions shape=(3,N)
    d_IPF = np.abs(d_IPF)
    d_IPF = np.sort(d_IPF, axis=0)
    d_IPF = d_IPF[[1, 0, 2],:]
    
    proj_Ds = stereoprojection_directions(d_IPF)
    Colors=stereotriangle_colors(proj_Ds)    
    
    return Colors

def stereotriangle_colors_from_eumats_dir(eumats,d=[1,0,0]):
    """
    Get IPF colors from orientation matrices.
    
    Input:
        eumats: numpy array (N, 3, 3) - Orientation matrices
        d: list/array (3,) - Sample direction (default: [1,0,0])
    
    Output:
        RGB: numpy array (N, 3) - RGB colors
    """
    #d - list of sample vector
    #eumats shape=(N,3,3)
    #d=[0,1,0]
    d_IPF = orilistMult(eumats,d)
    d_IPF = np.abs(d_IPF)
    d_IPF = np.sort(d_IPF, axis=0)
    d_IPF = d_IPF[[1, 0, 2],:]
    
    proj_Ds = stereoprojection_directions(d_IPF)
    Colors=stereotriangle_colors(proj_Ds)    
    
    return Colors


def stereotriangle_colors(proj_Ds):
    """
    Get IPF colors from projection coordinates.
    
    Input:
        proj_Ds: numpy array (2, N) - Projection coordinates
    
    Output:
        RGB: numpy array (N, 3) - RGB colors
    """
    #adapted from https://mathematica.stackexchange.com/questions/47492/how-to-create-an-inverse-pole-figure-color-map
    #Red point
    Rp=stereoprojection_directions(np.array([0,0,1]))
    XR=Rp[0]
    YR=Rp[1]
    #Green point
    Gp=stereoprojection_directions(np.array([1,0,1]))
    XG=Gp[0]
    YG=Gp[1]
    #Blue point
    Bp=stereoprojection_directions(np.array([1,1,1]))
    XB=Bp[0]
    YB=Bp[1]
    
    #Point O to be colored
    try:
        
        if len(proj_Ds.shape)==2:
            XO=proj_Ds[0,:]
            YO=proj_Ds[1,:]
        else:
            XO=proj_Ds[0]
            YO=proj_Ds[1]
    except:
            XO=proj_Ds[0]
            YO=proj_Ds[1]
            
    #XO=XB
    #YO=YB
    
    #Intersection Gx of GO with RB
    K1GO=(XG-XO)*(YR-YB)-(YG-YO)*(XR-XB)
    #O==G
    OeG=np.where(K1GO==0)
    K1GO[OeG]=1
    XGx=((XG*YO-YG*XO)*(XR-XB)-(XG-XO)*(XR*YB-YR*XB))/K1GO
    YGx=((XG*YO-YG*XO)*(YR-YB)-(YG-YO)*(XR*YB-YR*XB))/K1GO
    XGx[OeG]=XR
    YGx[OeG]=YR
    #Intersection Bx of BO with RG
    K1BO=(XB-XO)*(YR-YG)-(YB-YO)*(XR-XG)
    #O==B
    OeB=np.where(K1BO==0)
    K1BO[OeB]=1
    XBx=((XB*YO-YB*XO)*(XR-XG)-(XB-XO)*(XR*YG-YR*XG))/K1BO
    YBx=((XB*YO-YB*XO)*(YR-YG)-(YB-YO)*(XR*YG-YR*XG))/K1BO
    XBx[OeB]=XR
    YBx[OeB]=YR
    
    
    #Intersection Rx of RO with arc [101]-[111]
    #O==R
    OeR=np.where(XO==0)
    XO2=XO.copy()
    XO2[OeR]=1.0
    K1=(YO/XO2)**2
    XRx=(np.sqrt(K1+2)-1)/(K1+1)
    YRx=YO/XO2*XRx
    XRx[OeR]=XB
    YRx[OeR]=YB
    
    #ratios |ORx|//|RRx| |OGx|//|GGx| |OBx|//|BBx|
    RED=np.sqrt((XO-XRx)**2+(YO-YRx)**2)/np.sqrt((XR-XRx)**2+(YR-YRx)**2)
    GREEN=np.sqrt((XO-XGx)**2+(YO-YGx)**2)/np.sqrt((XG-XGx)**2+(YG-YGx)**2)
    BLUE=np.sqrt((XO-XBx)**2+(YO-YBx)**2)/np.sqrt((XB-XBx)**2+(YB-YBx)**2)
    
    
    
    Colors=np.vstack((RED,GREEN,BLUE)).T
    #Colors=Colors/np.max(Colors,axis=1)
    
    Colors=np.vstack((RED/np.max(Colors,axis=1),GREEN/np.max(Colors,axis=1),BLUE/np.max(Colors,axis=1))).T
    
    return Colors



def equivalent_elements(element,lattice):
    """
    Find symmetrically equivalent elements.
    
    Input:
        element: numpy array (3,) - Direction/normal vector
        lattice: str - Crystal system
    
    Output:
        equivalents: list of arrays - Equivalent vectors
    """
    eq_elements=[]
    R=symmetry_elements(lattice)
    eq_elements.append(R[0].dot(element))
    for u in R[1:]:
        el=u.dot(element)
        isin=False
        for els in eq_elements:
            if (els==el).all():
                isin=True
                break
        if not isin:
            eq_elements.append(el)
        #eq_elements.append(el)
         
    return eq_elements


def stereoprojection_intotriangle_ini(dirs,eps=1.0e-5):
    """
    Map directions into standard triangle (initial version).
    
    Input:
        dirs: numpy array (3, N) - Direction vectors
        eps: float - Tolerance (default: 1e-5)
    
    Output:
        proj: numpy array (2, N) - Projection coordinates
    """
    normals = np.array([-1,0,1]);
    arclength = 40.#-np.arccos(np.sqrt(2)/np.sqrt(3))*180/np.pi;
    proj_normals = stereoprojection_planes(normals,arclength=arclength,iniangle=90)
    proj_tans = np.arctan(proj_normals[1,:]/proj_normals[0,:])
    
    if len(dirs.shape)==1:
        dirs = np.expand_dims(dirs,axis=1)

    proj_dirs = np.zeros(dirs.shape)
    inc=-1
    for diri in dirs.T:
#        print('===================================================')
#        print(diri)
#        print('===================================================')
        inc+=1
        el=equivalent_elements(diri,'cubic')
        #print(el)
        for eli in el:
            proj_eli = stereoprojection_directions(eli)
#            print(eli)
#            print(np.arctan(proj_eli[1,0]/proj_eli[0,0])-np.arccos(1./np.sqrt(3.)))
#            print(np.arctan(proj_eli[1,0]/proj_eli[0,0]))#-np.pi/4)
#            print((np.arccos(abs(eli[2])/np.sqrt(eli.dot(eli)))-np.arccos(1./np.sqrt(3.))))
#            if (eli>=-eps).all() and (np.arccos(abs(eli[2])/np.sqrt(eli.dot(eli)))-np.arccos(1./np.sqrt(3.)))<eps:
#            if (proj_eli[:,0]>=-eps).all() and (np.arccos(abs(eli[2])/np.sqrt(eli.dot(eli)))-np.arccos(1./np.sqrt(3.)))<eps:
            if ((proj_eli[:,0])>=-eps).all():                #proj_eli = stereoprojection_directions(eli)
                atan=np.arctan2(proj_eli[1,0],proj_eli[0,0])
                if (atan-np.pi/4)<eps:
                    idx=np.where(abs(proj_tans-atan)==min(abs(proj_tans-atan)))[0][0]
#                    print(proj_eli[:,0].dot(proj_eli[:,0])) 
#                    print(proj_normals[:,idx].dot(proj_normals[:,idx])) 
#                    print((proj_eli[:,0].dot(proj_eli[:,0])-proj_normals[:,idx].dot(proj_normals[:,idx])))
#                    print((proj_eli[:,0].dot(proj_eli[:,0])-proj_normals[:,idx].dot(proj_normals[:,idx]))<eps)
                    if (proj_eli[:,0].dot(proj_eli[:,0])-proj_normals[:,idx].dot(proj_normals[:,idx]))<eps:
                        proj_dirs[:,inc]=proj_eli[:,0]
#                        print('OK')
#                        print(proj_dirs[:,inc])
                        #break
    return proj_dirs

def stereoprojection_intotriangle_fast(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False,Rin=None,symops=None):
    """
    Fast mapping of directions into standard triangle.
    
    Input:
        dirs: numpy array (3, N) - Direction vectors
        eps, geteqdirs, geteqmats, Rin, symops: Optional parameters
    
    Output:
        proj: numpy array (2, N) - Projection coordinates
    """
    etamax=np.arctan2(1,1)*180./np.pi
    if len(dirs.shape)==1:
        dirs = np.expand_dims(dirs,axis=1)
    if symops is None:
        symops = symmetry_elements('cubic')
    eqmats=np.array([np.eye(3) for ii in range(dirs.shape[1])])
    eqdirs = np.zeros(dirs.shape)
    Rout=np.zeros(eqmats.shape)
    RTout=np.zeros(eqmats.shape)
    
    for sym in symops:
	#could be faster if  if we remove directions found in each itteration
        #sym=np.eye(3)
        datas=sym.dot(dirs)
        idxs = np.where((datas[0,:]>=0 ) & (datas[1,:]>=0) & (datas[2,:]>=0))[0]
        eta=np.arctan2(datas[0,idxs],np.abs(datas[2,idxs]))*180./np.pi
        chi=np.arctan2(datas[1,idxs],datas[0,idxs])*180./np.pi
        idxs=idxs[np.where((eta<=etamax) & (chi<=etamax))[0]]
        eqmats[idxs,:,:]=sym
        eqdirs[:,idxs]=datas[:,idxs]
        #dirs=np.delete(dirs,idxs,1)
        if Rin is not None:
            for idxsi in idxs:
                Rout[idxsi,:,:]=sym.dot(Rin[idxsi])
                RTout[idxsi,:,:]=Rout[idxsi,:,:].T
        
        if dirs.shape[1]==0:
            break    #print('test')
    proj_dirs=stereoprojection_directions(eqdirs)
    out={}
    if Rin is not None:
        out['Rout']=Rout  
        out['RTout']=RTout  
    if geteqdirs:
        out['eqdirs']=eqdirs
    if geteqmats:
        out['eqmats']=eqmats
    if len(out)>0:
        return proj_dirs, out
    else:
        return proj_dirs

def stereoprojection_intotriangle(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False,Rin=None,symops=None):
    """
    Map directions into standard stereographic triangle.
    
    Input:
        dirs: numpy array (3, N) - Direction vectors
        eps, geteqdirs, geteqmats, Rin, symops: Optional parameters
    
    Output:
        proj: numpy array (2, N) - Projection coordinates
    """
    etamax=np.arctan2(1,1)*180./np.pi
    if len(dirs.shape)==1:
        dirs = np.expand_dims(dirs,axis=1)
    if symops is None:
        symops = symmetry_elements('cubic')
    proj_dirs = np.zeros(dirs.shape)
    eqdirs = np.zeros(dirs.shape)
    eqmats=[]
    inc=-1
    Rout=[]
    RTout=[]
    idx=-1
    for diri in dirs.T:
        idx+=1
        inc+=1
        br=False
        for sym in symops:
            Ds = sym.dot(diri)
            
            if Ds[0]>=0 and Ds[1]>=0 and Ds[2]>=0:
                eta=np.arctan2(Ds[0],np.abs(Ds[2]))*180./np.pi
                chi=np.arctan2(Ds[1],Ds[0])*180./np.pi#np.arcsin(Ds[2])*180./np.pi
                if eta<=etamax and chi<=etamax:# and np.abs(chi<=etamax)<=1.0e-5: #and eta<=np.pi/2 and chi<=etamax:# and chi>=chimax:   
                    break
                    br=True
        #if not br:
        #    print("sdddddddddddddddddd")
        proj_dirs[:,inc] = stereoprojection_directions(Ds)[:,0]
        eqmats.append(sym)
        if Rin is not None:
            Rout.append(sym.dot(Rin[idx]))
            #if np.linalg.det(Rout[-1])<0:
                #Rout[-1]=-1*Rout[-1]
            RTout.append(Rout[-1].T)
        eqdirs[:,inc]=Ds
    #print('test')
    out={}
    if Rin is not None:
        out['Rout']=Rout  
        out['RTout']=RTout  
    if geteqdirs:
        out['eqdirs']=eqdirs
    if geteqmats:
        out['eqmats']=eqmats
    if len(out)>0:
        return proj_dirs, out
    else:
        return proj_dirs

def equalarea_intotriangle_fast(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False,Rin=None,symops=None):
    """
    Fast equal-area projection into triangle.
    
    Input:
        dirs: numpy array (3, N) - Direction vectors
        eps, geteqdirs, geteqmats, Rin, symops: Optional parameters
    
    Output:
        proj: numpy array (2, N) - Projection coordinates
    """
    etamax=np.arctan2(1,1)*180./np.pi
    if len(dirs.shape)==1:
        dirs = np.expand_dims(dirs,axis=1)
    if symops is None:
        symops = symmetry_elements('cubic')
    eqmats=np.array([np.eye(3) for ii in range(dirs.shape[1])])
    eqdirs = np.zeros(dirs.shape)
    Rout=np.zeros(eqmats.shape)
    RTout=np.zeros(eqmats.shape)
    
    for sym in symops:
        #sym=np.eye(3)
        datas=sym.dot(dirs)
        idxs = np.where((datas[0,:]>=0 ) & (datas[1,:]>=0) & (datas[2,:]>=0))[0]
        eta=np.arctan2(datas[0,idxs],np.abs(datas[2,idxs]))*180./np.pi
        chi=np.arctan2(datas[1,idxs],datas[0,idxs])*180./np.pi
        idxs=idxs[np.where((eta<=etamax) & (chi<=etamax))[0]]
        eqmats[idxs,:,:]=sym
        eqdirs[:,idxs]=datas[:,idxs]
        dirs=np.delete(dirs,idxs,1)
        if Rin is not None:
            for idxsi in idxs:
                Rout[idxsi,:,:]=sym.dot(Rin[idxsi])
                Rout[idxsi,:,:]=Rout[idxsi,:,:].T
        
        if dirs.shape[1]==0:
            break    #print('test')
    proj_dirs=equalarea_directions(eqdirs)
    out={}
    if Rin is not None:
        out['Rout']=Rout  
        out['RTout']=RTout  
    if geteqdirs:
        out['eqdirs']=eqdirs
    if geteqmats:
        out['eqmats']=eqmats
    if len(out)>0:
        return proj_dirs, out
    else:
        return proj_dirs


def equalarea_intotriangle(dirs,eps=1.0e-5,geteqdirs=False,geteqmats=False):
    """
    Equal-area projection into standard triangle.
    
    Input:
        dirs: numpy array (3, N) - Direction vectors
        eps, geteqdirs, geteqmats: Optional parameters
    
    Output:
        proj: numpy array (2, N) - Projection coordinates
    """
    etamax=np.arctan2(1,1)*180./np.pi
    if len(dirs.shape)==1:
        dirs = np.expand_dims(dirs,axis=1)

    symops = symmetry_elements('cubic')
    proj_dirs = np.zeros(dirs.shape)
    eqdirs = np.zeros(dirs.shape)
    eqmats=[]
    inc=-1
    
    for diri in dirs.T:
        inc+=1
        br=False
        for sym in symops:
            Ds = sym.dot(diri)
            
            if Ds[0]>=0 and Ds[1]>=0 and Ds[2]>=0:
                eta=np.arctan2(Ds[0],np.abs(Ds[2]))*180./np.pi
                chi=np.arctan2(Ds[1],Ds[0])*180./np.pi#np.arcsin(Ds[2])*180./np.pi
                if eta<=etamax and chi<=etamax:# and np.abs(chi<=etamax)<=1.0e-5: #and eta<=np.pi/2 and chi<=etamax:# and chi>=chimax:   
                    break
                    br=True
        #if not br:
        #    print("sdddddddddddddddddd")
        proj_dirs[:,inc] = equalarea_directions(Ds)[:,0]
        eqmats.append(sym)
        eqdirs[:,inc]=Ds
#    if geteqdirs and not geteqmats:
#        return proj_dirs,eqdirs
#    elif geteqmats and not geteqdirs:
#        return proj_dirs,eqmats
#    elif eqmats and geteqdirs:
#        return proj_dirs,eqdirs,eqmats
#    else:
#        return proj_dirs


    out={}
    if geteqdirs:
        out['eqdirs']=eqdirs
    if geteqmats:
        out['eqmats']=eqmats
    if len(out)>0:
        return proj_dirs, out
    else:
        return proj_dirs




def stereoprojection_planes(normals,arclength=360.,iniangle=0.,hemisphere='both',getpoints=False,R=np.eye(3)):
    """
    Project plane traces onto stereographic projection.
    
    Input:
        normals: numpy array (3, N) - Plane normals
        arclength, iniangle, hemisphere, getpoints, R: Optional parameters
    
    Output:
        traces: list of arrays or plots
    """
    #%normals = [x1,x2,...,xn;y1,y2,...,yn;z1,z2,...,zn];
    #%varargin{1} arclength in deg
    #normals = np.transpose(np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1]]))
    #
    if len(normals.shape)==1:
        normals = np.expand_dims(normals,axis=1)

    normals = normals.astype(float)
    normals /= np.sqrt((normals ** 2).sum(0))
    #print(normals)
    proj_normals = stereoprojection_directions(normals)

    idxs = np.where(abs(normals[0,:])+abs(normals[1,:])==0)[0]
    
    inplanedirs = np.vstack((-normals[1,:],normals[0,:],np.zeros(normals[0,:].shape)));
    inplanedirs[:,idxs] = np.vstack((np.zeros(normals[0,idxs].shape), -normals[2,idxs],normals[1,idxs]));
    
    inplanedirs /= np.sqrt((inplanedirs ** 2).sum(0))

    thirdaxis=np.cross(normals,inplanedirs,axisa=0,axisb=0,axisc=0)
#    thirdaxis = np.vstack((normals[1,:]*inplanedirs[2,:]-normals[2,:]*inplanedirs[1,:],
#                           -1*(normals[0,:]*inplanedirs[2,:]-normals[2,:]*inplanedirs[0,:]),
#                           normals[0,:]*inplanedirs[1,:]-normals[1,:]*inplanedirs[0,:]));
    t=np.linspace(iniangle,iniangle+arclength,180*2+1)*np.pi/180;
    basicarc = np.vstack((np.cos(t),np.sin(t),np.zeros(t.shape)));
    
    proj_planes=[];
    plane_points=[]
    for i in range(0,normals.shape[1]):
        Rot2Global = np.transpose(np.vstack((inplanedirs[:,i],thirdaxis[:,i],normals[:,i])));
        Ccp = R.dot(Rot2Global.dot(basicarc))
        if hemisphere == "both":            
            Ds = stereoprojection_directions(Ccp)
        elif hemisphere == "triangle":  
            #idxs = np.where(Ccp[2,:]>=0)[0]
            Ds=stereoprojection_intotriangle(Ccp)#[:,idxs])
        else:
            if hemisphere == "upper":
                idxs = np.where(Ccp[2,:]>=0)[0]
            elif hemisphere == "lower":
                idxs = np.where(Ccp[2,:]<=0)[0]
            Ds = stereoprojection_directions(Ccp[:,idxs])
            Ccp = Ccp[:,idxs]

        proj_planes.append(Ds)
        plane_points.append(Ccp)
        
    if len(proj_planes)==1:
        proj_planes = proj_planes[0]
        plane_points=plane_points[0]
    
    if getpoints:
        return  proj_planes, plane_points   
    else:
        return proj_planes

def iszero(a):
    """Check if value is approximately zero (|a| < 1e-9)."""
    return abs(a)<1e-9


def gcd(a,b):
    """Compute greatest common divisor of two numbers."""
    if iszero(b):
        return a
    return gcd(b,a%b)


def gcdarr(arr):
    """Compute GCD of all array elements."""
    return reduce(gcd,arr)


def vector2miller(arr, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3):
    """
    Convert vector to Miller indices using GCD reduction.
    
    Input:
        arr: array (3,) - Vector [x, y, z]
        MIN, Tol, tol, text, decimals: Optional parameters
    
    Output:
        miller: array (3,) - Miller indices [h, k, l]
    """
    #arguments are used just to keep compatibility... here are not used except decimals
    arrgcd = gcdarr(arr)
    vm = [(a/arrgcd) for a in arr]
    
    if (np.round(vm)==vm).all():
        return np.round(vm)
    else:
        return np.around(vm,decimals=decimals)
    


def vector2millerround(v, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3):
    """
    Convert vector to Miller indices with rounding.
    
    Input:
        v: array (3,) - Vector [x, y, z]
        MIN, Tol, tol, text, decimals: Optional parameters
    
    Output:
        miller: array (3,) - Miller indices
    """
    vm = np.round(v/Tol)*Tol
    #print(vm)
    #print((np.round(vm)==vm).all())
    #if (vm==np.array([ 2. ,-1. , 1.])).all():
        #print((np.round(vm)==vm).all())
        #print(vm)
    if (np.abs(vm)<=1).all():
        vm=np.round(vm/abs(min(vm[np.abs(vm)>1/tol]))*tol)/tol
        vm/=min(abs(vm[np.nonzero(vm)[0]]))   
    if not (np.round(vm)==vm).all():
        if MIN:
            #print((min(abs(vm[np.abs(vm)>1/tol]))))
            vm=np.round(vm/abs(min(np.abs(vm[np.abs(vm)>1/tol])))*tol)/tol
            #if not (np.round(vm)==vm).all():
            #    vm=np.round(vm/abs(max(np.abs(vm[np.abs(vm)>1/tol])))*tol)/tol
            #vm/=min(abs(vm[np.nonzero(vm)[0]]))
        else:
            #print(vm)
            vm=np.round(vm/abs(max(np.abs(vm[np.abs(vm)>1/tol])))*tol)/tol
            #if not (np.round(vm)==vm).all():
            #    vm=np.round(vm/abs(min(np.abs(vm[np.abs(vm)>1/tol])))*tol)/tol
            #vm/=max(abs(vm[np.nonzero(vm)[0]]))
    else:
        #print(vm)
        vm=vm.astype('int')
        #print(vm)
        gcd=math.gcd(math.gcd(vm[0],vm[1]),vm[2])
        vm=vm.astype('float')
        vm/=gcd
    #if (vm==np.array([-2.,-1., 1.])).all():
    #    print('====================================================')
    #    print(v)
    #    print('====================================================')
    vm=np.around(vm,decimals=decimals)
    if text:
        if (vm==vm.astype(int)).all():
            vm=f"$[{{{int(vm[0])}}}{{{int(vm[1])}}}{{{int(vm[2])}}}]$".replace('{-','\\overline{')
        else:
            f"$[{{{(vm[0])}}}{{{(vm[1])}}}{{{(vm[2])}}}]$".replace('{-','\\overline{')
    return(vm)


def vectors2miller(V, MIN=True, Tol=1e-9,tol=1e5,text=False):
    """
    Convert multiple vectors to Miller indices.
    
    Input:
        V: numpy array (3, N) - Vectors
        MIN, Tol, tol, text: Optional parameters
    
    Output:
        VM: numpy array (3, N) - Miller indices
    """
    VM=[]
    for v in V.T:
        VM.append(vector2miller(v,MIN=MIN,Tol=Tol,tol=tol,text=text))
    return np.array(VM).T
 


# Projlib.py Quick Reference Guide

## Most Commonly Used Functions

### 1. Generate Orientations

```python
# Generate orientation grid (5° resolution, upper hemisphere)
oris = genori(dangle=5.0, hemi='upper')
# Returns: (3, N) array of direction vectors

# Generate cubic triangle grid (1° resolution)
oris_tri = genoritri(resolution=1.0)
# Returns: (3, N) array for stereographic triangle
```

### 2. Project to 2D

```python
# Equal-area projection (Schmidt net) - RECOMMENDED for texture analysis
proj_ea = equalarea_directions(oris)
# Returns: (3, N) array, use proj_ea[0:2] for XY coordinates

# Stereographic projection (Wulff net) - preserves angles
proj_stereo = stereoprojection_directions(oris)
# Returns: (3, N) array, use proj_stereo[0:2] for XY coordinates
```

### 3. Create Density Grid

```python
# Generate regular grid with density data
grid_x, grid_y, grid_z, nummask, mask = genprojgrid(
    oris,                    # Orientation vectors
    gdata=density_values,    # Your density/intensity data
    nump=501,                # Grid resolution
    proj='equalarea',        # Projection type
    gdout=True              # Return gridded data
)
```

### 4. Generate Crystal Directions with Symmetry

```python
import numpy as np
from orilib import symmetry_elements

# Cubic crystal setup
a = 3.5  # Lattice parameter
L = a * np.eye(3)                # Direct basis
Lr = (2*np.pi/a) * np.eye(3)    # Reciprocal basis

# Important directions and planes
uvws = [[1,0,0], [1,1,0], [1,1,1]]
hkls = [[1,0,0], [1,1,0], [1,1,1]]

# Generate with cubic symmetry
symops = symmetry_elements('cubic')
dirs, normals = gen_dirs_norms(L, Lr, uvws, hkls, 
                                symops=symops, recsymops=symops)
```

### 5. Coordinate Transformations

```python
# Cartesian to spherical
theta, phi = xyz2spher(xyz_coords, deg=True)

# Spherical to Cartesian
xyz = spher2xyz(theta, phi, deg=True)
```

## Complete Workflow Examples

### Example 1: Basic Pole Figure

```python
import numpy as np
import matplotlib.pyplot as plt

# Step 1: Generate orientations
oris = genori(dangle=2.0, hemi='upper')

# Step 2: Create random density data
density = np.random.rand(oris.shape[1]) * 100

# Step 3: Generate grid
grid_x, grid_y, grid_z, _, _ = genprojgrid(
    oris, gdata=density, nump=301, proj='equalarea', gdout=True
)

# Step 4: Plot
fig, ax = plt.subplots(figsize=(8, 8))
contours = ax.contourf(grid_x, grid_y, grid_z, levels=20, cmap='jet')
plt.colorbar(contours, label='Density')
circle = plt.Circle((0, 0), np.sqrt(2), fill=False, color='k', linewidth=2)
ax.add_patch(circle)
ax.set_aspect('equal')
plt.title('Pole Figure')
plt.show()
```

### Example 2: Pole Figure with Crystal Directions

```python
import numpy as np
import matplotlib.pyplot as plt
from orilib import symmetry_elements

# Setup cubic crystal
a = 3.5
L = a * np.eye(3)
Lr = (2*np.pi/a) * np.eye(3)
uvws = [[1,0,0], [1,1,0], [1,1,1]]
hkls = [[1,0,0], [1,1,0], [1,1,1]]

# Generate symmetry-equivalent directions
symops = symmetry_elements('cubic')
dirs, normals = gen_dirs_norms(L, Lr, uvws, hkls, 
                                symops=symops, recsymops=symops)

# Plot
fig, ax = plt.subplots(figsize=(8, 8))

# Draw circle
circle = plt.Circle((0, 0), np.sqrt(2), fill=False, color='k', linewidth=2)
ax.add_patch(circle)

# Plot plane traces
for normal in normals:
    plane = normal['equalarea plane']
    ax.plot(plane[0], plane[1], 'k-', linewidth=0.5)

# Plot direction points
for dir in dirs:
    proj = dir['equalarea']
    ax.plot(proj[0], proj[1], 'ro', markersize=8)
    ax.text(proj[0]+0.05, proj[1]+0.05, dir['label'], fontsize=10)

ax.set_aspect('equal')
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
plt.title('Cubic Crystal Pole Figure')
plt.show()
```

### Example 3: Stereographic Triangle

```python
import matplotlib.pyplot as plt

# Generate triangle orientations
oris_tri = genoritri(resolution=1.0)

# Project
proj = equalarea_directions(oris_tri)

# Plot
fig, ax = stereotriangle(equalarea=True)
ax.scatter(proj[0], proj[1], s=1, alpha=0.5)
plt.title('Stereographic Triangle')
plt.show()
```

## Function Parameters Quick Reference

### genori()
- `dangle`: Angular resolution in degrees (default: 1.0)
- `hemi`: 'both', 'upper', 'lower' (default: 'both')
- `tol`: Hemisphere tolerance (default: 1e-2)
- `rot`: Rotation matrix (default: identity)
- `half`: Projection hemisphere filter (default: 'no')

### genprojgrid()
- `oris`: (3, N) orientation vectors
- `gdata`: (N,) data values (optional)
- `nump`: Grid points per dimension (default: 1001)
- `proj`: 'equalarea' or 'stereo' (default: 'equalarea')
- `method2`: 'linear', 'nearest', 'cubic' (default: 'linear')
- `gdout`: Return data grid (default: False)

### gen_dirs_norms()
- `L`: (3, 3) direct lattice basis
- `Lr`: (3, 3) reciprocal lattice basis
- `uvws`: List of [u,v,w] directions
- `hkls`: List of [h,k,l] planes
- `R2Proj`: Rotation to projection frame (default: identity)
- `symops`: Direct space symmetry operations (optional)
- `recsymops`: Reciprocal space symmetry operations (optional)

## Projection Formulas

### Equal-Area (Schmidt)
- **Formula**: r = √2 × sin(θ/2)
- **Max radius**: √2 (hemisphere)
- **Property**: Preserves area
- **Best for**: Texture analysis, pole figures

### Stereographic (Wulff)
- **Formula**: r = tan(θ/2)
- **Max radius**: ∞ (but typically 1 for hemisphere)
- **Property**: Preserves angles
- **Best for**: Angular relationships

## Common Pitfalls

1. **Array shapes**: Functions expect (3, N) not (N, 3)
   ```python
   # Wrong:
   oris = np.array([[1,0,0], [0,1,0]])  # Shape (2, 3)
   
   # Correct:
   oris = np.array([[1,0,0], [0,1,0]]).T  # Shape (3, 2)
   ```

2. **Projection coordinates**: Use [:2] or [0:2] to get X,Y
   ```python
   proj = equalarea_directions(oris)
   x_coords = proj[0]  # X coordinates
   y_coords = proj[1]  # Y coordinates
   xy_coords = proj[0:2]  # Both X and Y
   ```

3. **Grid data**: Set `gdout=True` when providing data
   ```python
   # With data:
   grid_x, grid_y, grid_z, _, _ = genprojgrid(oris, gdata=data, gdout=True)
   
   # Without data:
   grid_x, grid_y, _, _ = genprojgrid(oris)
   ```

## Performance Tips

1. **Lower resolution for preview**: Use dangle=10.0 or higher
2. **Higher resolution for publication**: Use dangle=1.0 or lower
3. **Grid resolution**: 301-501 for preview, 1001-2001 for publication
4. **Interpolation**: 'linear' is fastest, 'cubic' is smoothest
5. **Pre-compute projections**: Store proj = equalarea_directions(oris) if reusing

## Integration with Matplotlib

```python
import matplotlib.pyplot as plt
import numpy as np

# Setup figure
fig, ax = plt.subplots(figsize=(8, 8))

# Draw projection circle
if proj_type == 'equalarea':
    radius = np.sqrt(2)
else:  # stereo
    radius = 1.0
    
circle = plt.Circle((0, 0), radius, fill=False, color='k', linewidth=2)
ax.add_patch(circle)

# Plot data
ax.contourf(grid_x, grid_y, grid_z, levels=20, cmap='jet')

# Format
ax.set_aspect('equal')
ax.set_xlim(-1.5*radius, 1.5*radius)
ax.set_ylim(-1.5*radius, 1.5*radius)
plt.title('Pole Figure')
plt.show()
```

## For More Information

See the full documentation in `projlib_commented.py` where each function has:
- Comprehensive parameter descriptions
- Multiple usage examples
- Mathematical formulas
- Integration examples

Quick tips:
- Use equal-area for texture analysis (preserves density)
- Use stereographic for angular relationships
- Start with coarse grids, refine for final plots
- Always normalize your direction vectors (functions do this automatically)

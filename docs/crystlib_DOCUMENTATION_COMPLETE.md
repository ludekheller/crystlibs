# CRYSTLIB - Complete Documentation

**Module**: `crystlib.py`  
**Purpose**: Crystal Structure and Lattice Operations  
**Functions**: 76  
**Last Updated**: December 08, 2025

---

## Overview

Comprehensive functions for crystallographic calculations including lattice vectors, Miller indices, deformation analysis, and twinning operations.

---

## Table of Contents

1. [B19p_B2_lattice_correspondence](#function-b19p_b2_lattice_correspondence)
2. [B19p_B2_lattice_correspondence_ini](#function-b19p_b2_lattice_correspondence_ini)
3. [Rp_B2_lattice_correspondence](#function-rp_b2_lattice_correspondence)
4. [an_between_vecs](#function-an_between_vecs)
5. [array2tuple](#function-array2tuple)
6. [cubic2tetragonal_lattice_correspondence](#function-cubic2tetragonal_lattice_correspondence)
7. [cubic_lattice_vec](#function-cubic_lattice_vec)
8. [def_gradient](#function-def_gradient)
9. [def_gradient_ini](#function-def_gradient_ini)
10. [def_gradient_ini2](#function-def_gradient_ini2)
11. [def_gradient_stressfree](#function-def_gradient_stressfree)
12. [def_gradient_stressfree_ini](#function-def_gradient_stressfree_ini)
13. [dir2string](#function-dir2string)
14. [equivalent_elements](#function-equivalent_elements)
15. [find_gcd](#function-find_gcd)
16. [flipvector](#function-flipvector)
17. [flipvector2negative](#function-flipvector2negative)
18. [gen_twinned_lattice_points](#function-gen_twinned_lattice_points)
19. [genallHexSys](#function-genallhexsys)
20. [generate_hkls](#function-generate_hkls)
21. [generate_hkls01](#function-generate_hkls01)
22. [generate_hkls02](#function-generate_hkls02)
23. [generate_lattice_faces](#function-generate_lattice_faces)
24. [generate_lattice_points](#function-generate_lattice_points)
25. [generate_lattice_vectors](#function-generate_lattice_vectors)
26. [generate_lattite_atom_positions](#function-generate_lattite_atom_positions)
27. [generate_plane_vertices](#function-generate_plane_vertices)
28. [generate_product_lattice_faces](#function-generate_product_lattice_faces)
29. [generate_product_lattice_points](#function-generate_product_lattice_points)
30. [gensystemsHex](#function-gensystemshex)
31. [gensystemsHexIni](#function-gensystemshexini)
32. [get_interface2d](#function-get_interface2d)
33. [get_twinning_dislocation](#function-get_twinning_dislocation)
34. [get_twinning_plane_points](#function-get_twinning_plane_points)
35. [get_twinningdata](#function-get_twinningdata)
36. [get_unique_families](#function-get_unique_families)
37. [habitplane_equation_solution](#function-habitplane_equation_solution)
38. [hkil2hkl](#function-hkil2hkl)
39. [hkl2hkil](#function-hkl2hkil)
40. [kronecker](#function-kronecker)
41. [lattice_correspondence](#function-lattice_correspondence)
42. [lattice_vec](#function-lattice_vec)
43. [miller2fractional](#function-miller2fractional)
44. [mohr_circles](#function-mohr_circles)
45. [monoclinic_lattice_vec](#function-monoclinic_lattice_vec)
46. [niti_twinning](#function-niti_twinning)
47. [normArrayColumns](#function-normarraycolumns)
48. [np_kronecker](#function-np_kronecker)
49. [np_permut_tensor3](#function-np_permut_tensor3)
50. [permut_tensor3](#function-permut_tensor3)
51. [perpendicular_vector](#function-perpendicular_vector)
52. [plane2string](#function-plane2string)
53. [plane_line_intersection](#function-plane_line_intersection)
54. [print_correspondence](#function-print_correspondence)
55. [read_txt](#function-read_txt)
56. [reciprocal_basis](#function-reciprocal_basis)
57. [select_atomic_plane](#function-select_atomic_plane)
58. [select_atomic_region](#function-select_atomic_region)
59. [select_crystal_planes](#function-select_crystal_planes)
60. [select_plane](#function-select_plane)
61. [strains_along_13mohrcirle](#function-strains_along_13mohrcirle)
62. [tetragonal_lattice_vec](#function-tetragonal_lattice_vec)
63. [twin_equation_solution](#function-twin_equation_solution)
64. [twin_equation_solution_ini](#function-twin_equation_solution_ini)
65. [twinnedhabitplane](#function-twinnedhabitplane)
66. [uvtw2uvw](#function-uvtw2uvw)
67. [uvw2uvtw](#function-uvw2uvtw)
68. [vec2string](#function-vec2string)
69. [vector2miller_ini](#function-vector2miller_ini)
70. [vectors2miller](#function-vectors2miller)
71. [write_lattice_correspondence](#function-write_lattice_correspondence)
72. [write_mohr_planes](#function-write_mohr_planes)
73. [write_txt](#function-write_txt)
74. [xyz2fractional](#function-xyz2fractional)
75. [xyz2fractional02](#function-xyz2fractional02)
76. [zero_normal_strains](#function-zero_normal_strains)

---

## Function: B19p_B2_lattice_correspondence

**Signature**:
```python
def B19p_B2_lattice_correspondence(notation='Miyazaki'):
```

**Description**:

Generate B19'→B2 lattice correspondence matrix for NiTi.
    Returns the transformation matrix relating B19' monoclinic martensite
    to B2 cubic austenite lattice vectors. Fundamental for NiTi shape
    memory alloy analysis.

**Input**:

None

**Output**:

numpy.ndarray (3×3): Correspondence matrix C where L_B2 = C · L_B19'

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Get correspondence matrix
        >>> C = B19p_B2_lattice_correspondence()
        >>> print("B19' → B2 correspondence:")
        >>> print(C)
        >>> 
        >>> # Apply to B19' lattice
        >>> L_B19p = monoclinic_lattice_vec(2.889, 4.120, 4.622, 96.8)
        >>> L_B2_transformed = C.dot(L_B19p)
        >>> 
        >>> # Compare with actual B2 lattice
        >>> L_B2_actual = cubic_lattice_vec(3.015)
```

**Notes**:

- Specific to NiTi shape memory alloys
        - Based on crystallographic theory
        - Used in transformation analysis
        - Related to habit plane calculations
    Formula:
        C = [[c11, c12, c13],
             [c21, c22, c23],
             [c31, c32, c33]]
        where coefficients determined experimentally

---

## Function: B19p_B2_lattice_correspondence_ini

**Signature**:
```python
def B19p_B2_lattice_correspondence_ini():
```

**Description**:

Initialize B19'→B2 correspondence for specific transformation variant.
    Returns correspondence matrix for one of the crystallographic variants
    of the B19' martensite transformation in NiTi.

**Input**:

variant (int): Variant number (1-24 for cubic→monoclinic)

**Output**:

numpy.ndarray (3×3): Correspondence matrix for this variant

**Usage Example**:

```python
>>> # Get correspondence for variant 1
        >>> C1 = B19p_B2_lattice_correspondence_ini(variant=1)
        >>> print("Variant 1 correspondence:")
        >>> print(C1)
        >>> 
        >>> # Get all variants
        >>> all_variants = []
        >>> for i in range(1, 25):
        ...     C = B19p_B2_lattice_correspondence_ini(variant=i)
        ...     all_variants.append(C)
        >>> print(f"Total variants: {len(all_variants)}")
```

**Notes**:

- NiTi has 24 crystallographic variants
        - Each variant has different correspondence
        - Related to symmetry of parent phase
        - Used in texture simulations

---

## Function: Rp_B2_lattice_correspondence

**Signature**:
```python
def Rp_B2_lattice_correspondence():
```

**Description**:

Generate R-phase→B2 lattice correspondence for NiTi.
    Returns correspondence matrix for R-phase (rhombohedral) to B2
    transformation in NiTi alloys. R-phase is intermediate phase.

**Input**:

None

**Output**:

numpy.ndarray (3×3): R-phase → B2 correspondence

**Usage Example**:

```python
>>> C_R_B2 = Rp_B2_lattice_correspondence()
        >>> print("R-phase → B2 correspondence:")
        >>> print(C_R_B2)
```

**Notes**:

- R-phase is pre-martensitic phase in NiTi
        - Trigonal/rhombohedral structure
        - Appears above Ms temperature
        - Correspondence simpler than B19'→B2

---

## Function: an_between_vecs

**Signature**:
```python
def an_between_vecs(v1,v2,deg=True,full2pi=False):
```

**Description**:

Calculate angle between two vectors.
    Computes angle using dot product formula.

**Input**:

v1, v2 (array [3]): Vectors
        degrees (bool): Return in degrees (default: True), else radians

**Output**:

float: Angle between vectors

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # 90° angle
        >>> v1 = np.array([1, 0, 0])
        >>> v2 = np.array([0, 1, 0])
        >>> angle = an_between_vecs(v1, v2)
        >>> print(f"Angle: {angle}°")  # 90.0
        >>> 
        >>> # 45° angle
        >>> v1 = np.array([1, 0, 0])
        >>> v2 = np.array([1, 1, 0])
        >>> angle = an_between_vecs(v1, v2)
        >>> print(f"Angle: {angle:.1f}°")  # 45.0
```

**Notes**:

- Handles non-normalized vectors
        - Returns angle in [0, 180°] or [0, π]
        - Used throughout module
    Formula:
        cos(θ) = (v1 · v2) / (||v1|| ||v2||)

---

## Function: array2tuple

**Signature**:
```python
def array2tuple(arr, decimals=2):
```

**Description**:

Convert a numpy array to a tuple with rounded elements.

**Input**:

arr: numpy array or list - Array of numerical values
        decimals: int - Number of decimal places to round to (default: 2)

**Output**:

tuple - Rounded values as a tuple

**Usage Example**:

```python
>>> arr = np.array([1.234, 2.567, 3.891])
        >>> result = array2tuple(arr, decimals=2)
        >>> print(result)
        (1.23, 2.57, 3.89)
        >>> arr2 = [0.1234, 0.5678]
        >>> result2 = array2tuple(arr2, decimals=1)
        >>> print(result2)
        (0.1, 0.6)
```

---

## Function: cubic2tetragonal_lattice_correspondence

**Signature**:
```python
def cubic2tetragonal_lattice_correspondence():
```

**Description**:

Generate cubic→tetragonal lattice correspondence.
    Returns correspondence matrix for cubic to tetragonal phase
    transformation (e.g., FCC→FCT).

**Input**:

None

**Output**:

numpy.ndarray (3×3): Correspondence matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> C = cubic2tetragonal_lattice_correspondence()
        >>> 
        >>> # Apply to cubic lattice
        >>> L_cubic = cubic_lattice_vec(3.0)
        >>> L_tetra_approx = C.dot(L_cubic)
        >>> 
        >>> # Adjust c-parameter
        >>> L_tetra_approx[2,2] *= 1.333  # c/a ratio
```

**Notes**:

- Common in martensitic transformations
        - Simple correspondence: diagonal matrix
        - c/a ratio determines tetragonality

---

## Function: cubic_lattice_vec

**Signature**:
```python
def cubic_lattice_vec(a):
```

**Description**:

Generate cubic lattice vectors.
    Creates 3×3 lattice matrix for cubic crystal system with parameter a.
    All angles are 90° and all lengths are equal (a=b=c).

**Input**:

a (float): Cubic lattice parameter (Ångströms)

**Output**:

numpy.ndarray (3×3): Lattice matrix [a1|a2|a3] where columns are lattice vectors

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # B2 austenite (NiTi)
        >>> LA = cubic_lattice_vec(3.015)
        >>> print(LA)
        >>> # [[3.015, 0,     0    ],
        >>> #  [0,     3.015, 0    ],
        >>> #  [0,     0,     3.015]]
        >>> 
        >>> # Calculate volume
        >>> volume = np.linalg.det(LA)
        >>> print(f"Unit cell volume: {volume:.3f} ų")
```

**Notes**:

- Cubic system: a=b=c, α=β=γ=90°
        - Common in metals (FCC, BCC, SC)
        - Lattice matrix L = diag([a, a, a])
        - Volume = a³
    Formula:
        L = [[a, 0, 0],
             [0, a, 0],
             [0, 0, a]]

---

## Function: def_gradient

**Signature**:
```python
def def_gradient(Cd,LA, LM,StressT=np.zeros((3,3)),STA=np.zeros((3,3,3,3)),STM=np.zeros((3,3,3,3)),CId=None):
```

**Description**:

def_gradient - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: def_gradient_ini

**Signature**:
```python
def def_gradient_ini(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,StressT=np.zeros((3,3)),Parent_ST=np.zeros((3,3,3,3)),Product_ST=np.zeros((3,3,3,3))):
```

**Description**:

def_gradient_ini - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: def_gradient_ini2

**Signature**:
```python
def def_gradient_ini2(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,StressT=np.zeros((3,3)),Parent_ST=np.zeros((3,3,3,3)),Product_ST=np.zeros((3,3,3,3))):
```

**Description**:

def_gradient_ini2 - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: def_gradient_stressfree

**Signature**:
```python
def def_gradient_stressfree(Cd,LA, LM,CId=None):
```

**Description**:

Calculate stress-free transformation strain (deformation gradient).
    Computes F₀ = L2 · L1⁻¹ representing lattice transformation
    without applied stress.

**Input**:

L1 (array 3×3): Initial lattice matrix
        L2 (array 3×3): Final lattice matrix

**Output**:

numpy.ndarray (3×3): Stress-free deformation gradient F₀

**Usage Example**:

```python
>>> # B2 → B19' transformation in NiTi
        >>> L_B2 = cubic_lattice_vec(3.015)
        >>> L_B19p = monoclinic_lattice_vec(2.889, 4.120, 4.622, 96.8)
        >>> 
        >>> F0 = def_gradient_stressfree(L_B2, L_B19p)
        >>> print("Stress-free transformation strain:")
        >>> print(F0)
        >>> 
        >>> # Calculate principal strains
        >>> eigenvalues = np.linalg.eigvalsh((F0 + F0.T)/2 - np.eye(3))
        >>> print(f"Principal strains: {eigenvalues}")
```

**Notes**:

- F₀ represents pure lattice change
        - No external stress applied
        - Used in transformation theory
        - Basis for habit plane calculations
    Formula:
        F₀ = L₂ · L₁⁻¹

---

## Function: def_gradient_stressfree_ini

**Signature**:
```python
def def_gradient_stressfree_ini(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,Ci_d=None):
```

**Description**:

Initialize stress-free deformation gradient from lattice parameters.
    Convenience function to calculate F₀ from lattice parameters directly.

**Input**:

lattice_params_1 (tuple): (a, b, c, α, β, γ) for initial structure
        lattice_params_2 (tuple): (a, b, c, α, β, γ) for final structure

**Output**:

numpy.ndarray (3×3): Deformation gradient

**Usage Example**:

```python
>>> # Cubic → Tetragonal
        >>> params_cubic = (3.0, 3.0, 3.0, 90, 90, 90)
        >>> params_tetra = (3.0, 3.0, 4.0, 90, 90, 90)
        >>> 
        >>> F0 = def_gradient_stressfree_ini(params_cubic, params_tetra)
        >>> print("Transformation strain:")
        >>> print(F0)
```

**Notes**:

- Wrapper around def_gradient_stressfree()
        - Accepts lattice parameters directly
        - More convenient for quick calculations

---

## Function: dir2string

**Signature**:
```python
def dir2string(v, digits=2):
```

**Description**:

Format direction as string representation [u,v,w] with Miller notation.
    Converts direction vector to crystallographic notation using
    square brackets, following standard Miller index convention.

**Input**:

v (array-like [3]): Direction vector [u, v, w]
        digits (int, optional): Number of decimal places (default: 2)

**Output**:

str: Formatted string '[u.uu,v.vv,w.ww]'

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Format [110] direction
        >>> uvw = [1, 1, 0]
        >>> print(dir2string(uvw))  # Output: '[1.00,1.00,0.00]'
        >>> 
        >>> # Format with higher precision
        >>> direction = np.array([1.4142, 1.0, 0.0])
        >>> print(dir2string(direction, digits=4))  # Output: '[1.4142,1.0000,0.0000]'
        >>> 
        >>> # Use in orientation relationship description
        >>> print(f"Orientation: {dir2string([1,1,1])} || {plane2string([1,0,0])}")
```

**Notes**:

- Uses square brackets [] for direction notation (Miller convention)
        - Compare with plane2string() which uses parentheses ()
        - Standard crystallographic notation
        - Handles negative indices (bar notation not included)

---

## Function: equivalent_elements

**Signature**:
```python
def equivalent_elements(element,lattice):
```

**Description**:

Find symmetrically equivalent elements.

**Input**:

element: numpy array (3,) - Direction/normal vector
        lattice: str - Crystal system

**Output**:

equivalents: list of arrays - Equivalent vectors

---

## Function: find_gcd

**Signature**:
```python
def find_gcd(x, y):
```

**Description**:

Find greatest common divisor using recursive Euclidean algorithm.
    Fundamental mathematical operation used throughout the module for
    Miller indices reduction and fractional coordinate normalization.
    Implements the classical Euclidean algorithm via recursion.

**Input**:

x (int or float): First number
        y (int or float): Second number

**Output**:

int or float: Greatest common divisor of x and y

**Usage Example**:

```python
>>> # Basic GCD calculation
        >>> gcd = find_gcd(48, 18)
        >>> print(gcd)  # Output: 6
        >>> 
        >>> # Reduce Miller indices [6, 9, 12] to lowest terms
        >>> h, k, l = 6, 9, 12
        >>> gcd_hkl = find_gcd(find_gcd(h, k), l)
        >>> reduced = (h//gcd_hkl, k//gcd_hkl, l//gcd_hkl)
        >>> print(reduced)  # Output: (2, 3, 4)
        >>> 
        >>> # Works with floats
        >>> gcd_float = find_gcd(4.5, 3.0)
        >>> print(gcd_float)  # Output: 1.5
```

**Notes**:

- Used extensively in miller2fractional() function
        - Time complexity: O(log(min(x,y)))
        - Handles both integer and floating-point inputs
        - Returns 0 if both inputs are 0
    Formula:
        gcd(x, y) = gcd(y, x mod y) if y ≠ 0
        gcd(x, 0) = x (base case)

---

## Function: flipvector

**Signature**:
```python
def flipvector(v, Tol=1e-9):
```

**Description**:

Flip vector to ensure positive first non-zero component.
    Standardizes vector direction for comparison.

**Input**:

v (array [3]): Vector

**Output**:

numpy.ndarray [3]: Flipped vector (if needed)

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> v1 = np.array([-1, 2, 3])
        >>> v1_flipped = flipvector(v1)
        >>> print(v1_flipped)  # [1, -2, -3]
        >>> 
        >>> v2 = np.array([0, -1, 2])
        >>> v2_flipped = flipvector(v2)
        >>> print(v2_flipped)  # [0, 1, -2]
```

**Notes**:

- Makes first non-zero component positive
        - Useful for standardization
        - Preserves magnitude and line

---

## Function: flipvector2negative

**Signature**:
```python
def flipvector2negative(v, Tol=1e-9):
```

**Description**:

Flip vector to ensure negative first non-zero component.
    Opposite of flipvector() - ensures negative first component.

**Input**:

v (array [3]): Vector

**Output**:

numpy.ndarray [3]: Flipped vector (if needed)

**Usage Example**:

```python
>>> v = np.array([1, 2, 3])
        >>> v_neg = flipvector2negative(v)
        >>> print(v_neg)  # [-1, -2, -3]
```

**Notes**:

- Makes first non-zero component negative
        - Complementary to flipvector()

---

## Function: gen_twinned_lattice_points

**Signature**:
```python
def gen_twinned_lattice_points(ParentLatticePoints,eta1,shear_angle,K1,shift=0.0,dK1=None,bvr=None,deta1=None):
```

**Description**:

Generate lattice points for parent and twinned regions.
    Creates atomic positions showing twin boundary.

**Input**:

L_parent (array 3×3): Parent lattice
        twin_plane_normal (array [3]): Twin plane (K1)
        n1, n2, n3 (int): Cell range

**Output**:

tuple: (parent_atoms, twin_atoms) - atomic positions

**Usage Example**:

```python
>>> L = cubic_lattice_vec(3.0)
        >>> K1 = np.array([1, 1, 1])  # (111) twin
        >>> 
        >>> parent_pts, twin_pts = gen_twinned_lattice_points(L, K1, 10, 10, 10)
        >>> 
        >>> # Visualize
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> ax.scatter(parent_pts[:,0], parent_pts[:,1], parent_pts[:,2], 
        ...            c='blue', label='Parent')
        >>> ax.scatter(twin_pts[:,0], twin_pts[:,1], twin_pts[:,2], 
        ...            c='red', label='Twin')
        >>> ax.legend()
```

**Notes**:

- Shows twin boundary clearly
        - Used for visualization
        - Can calculate boundary energy

---

## Function: genallHexSys

**Signature**:
```python
def genallHexSys():
```

**Description**:

Generate all hexagonal slip systems including <c+a> pyramidal.
    Comprehensive slip system generation for HCP crystals including
    basal, prismatic, and pyramidal <c+a> systems. Essential for
    complete plasticity modeling of hexagonal materials.

**Input**:

L (array 3×3): Hexagonal lattice matrix

**Output**:

tuple: (all_directions, all_normals)
            - Cartesian slip directions
            - Cartesian slip plane normals

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Magnesium lattice
        >>> a = 3.21
        >>> c = 5.21
        >>> L_mg = lattice_vec(a, a, c, 90, 90, 120)
        >>> 
        >>> dirs, norms = genallHexSys(L_mg)
        >>> print(f"Total systems (with <c+a>): {len(dirs)}")
        >>> 
        >>> # Identify system types by direction
        >>> for i, d in enumerate(dirs):
        ...     if np.abs(d[2]) < 0.01:
        ...         print(f"System {i}: <a> type")
        ...     else:
        ...         print(f"System {i}: <c+a> type")
```

**Notes**:

- Includes basal {0001}<11-20>
        - Includes prismatic {1-100}<11-20>
        - Includes pyramidal {11-22}<11-23> (<c+a>)
        - <c+a> slip critical for c-axis strain
        - Total ~18-24 systems depending on implementation

---

## Function: generate_hkls

**Signature**:
```python
def generate_hkls(hklmax, syms, hkls=[]):
```

**Description**:

Generate unique Miller indices (hkl) considering crystal symmetry operations.
    This version rounds values to avoid floating point comparison issues.

**Input**:

hklmax: int - Maximum Miller index value (generates from -hklmax to +hklmax)
        syms: list of numpy arrays - List of 3x3 symmetry operation matrices
        hkls: list - Optional custom list of Miller indices to consider (default: [])

**Output**:

hkls: list of tuples - Unique Miller indices
        hkls2: dict - Dictionary mapping each unique hkl to its symmetry equivalents
        fam: dict - Dictionary of unique families with their multiplicities

**Usage Example**:

```python
>>> import numpy as np
        >>> # Define cubic symmetry operations (identity and 90° rotation around z)
        >>> identity = np.eye(3)
        >>> rot_z = np.array([[0, -1, 0],
        ...                   [1,  0, 0],
        ...                   [0,  0, 1]])
        >>> syms = [identity, rot_z]
        >>> 
        >>> # Generate hkls up to index 2
        >>> hkls, hkls2, fam = generate_hkls(2, syms)
        >>> print("First few unique hkls:", hkls[:3])
        >>> print("Symmetry equivalents of (1,0,0):", hkls2.get((1.0, 0.0, 0.0)))
        >>> print("Families:", list(fam.items())[:3])
```

---

## Function: generate_hkls01

**Signature**:
```python
def generate_hkls01(hklmax, syms, hkls=[]):
```

**Description**:

Generate unique Miller indices (hkl) considering crystal symmetry operations.
    Alternative version without rounding.

**Input**:

hklmax: int - Maximum Miller index value (generates from -hklmax to +hklmax)
        syms: list of numpy arrays - List of 3x3 symmetry operation matrices
        hkls: list - Optional custom list of Miller indices to consider (default: [])

**Output**:

hkls: list of tuples - Unique Miller indices
        hkls2: dict - Dictionary mapping each unique hkl to its symmetry equivalents
        fam: dict - Dictionary of unique families with their multiplicities

**Usage Example**:

```python
>>> import numpy as np
        >>> # Define tetragonal symmetry operations
        >>> identity = np.eye(3)
        >>> rot_90 = np.array([[0, -1, 0],
        ...                    [1,  0, 0],
        ...                    [0,  0, 1]])
        >>> syms = [identity, rot_90]
        >>> 
        >>> hkls, hkls2, fam = generate_hkls01(1, syms)
        >>> print("Generated hkls:", hkls)
        >>> print("Number of unique planes:", len(hkls))
```

---

## Function: generate_hkls02

**Signature**:
```python
def generate_hkls02(hklmax, syms, G, hkls=[]):
```

**Description**:

Generate unique Miller indices (hkl) with metric tensor transformation.
    This version applies symmetry operations in the correct metric space.

**Input**:

hklmax: int - Maximum Miller index value (generates from -hklmax to +hklmax)
        syms: list of numpy arrays - List of 3x3 symmetry operation matrices
        G: numpy array (3x3) - Metric tensor matrix
        hkls: list - Optional custom list of Miller indices to consider (default: [])

**Output**:

hkls: list of tuples - Unique Miller indices
        hkls2: dict - Dictionary mapping each unique hkl to its symmetry equivalents
        fam: dict - Dictionary of unique families with their multiplicities

**Usage Example**:

```python
>>> import numpy as np
        >>> # Define hexagonal metric tensor
        >>> a = 3.0  # lattice parameter
        >>> c = 5.0  # c-axis parameter
        >>> G = np.array([[a**2, -a**2/2, 0],
        ...               [-a**2/2, a**2, 0],
        ...               [0, 0, c**2]])
        >>> 
        >>> # Define hexagonal symmetry
        >>> identity = np.eye(3)
        >>> rot_60 = np.array([[0.5, -np.sqrt(3)/2, 0],
        ...                    [np.sqrt(3)/2, 0.5, 0],
        ...                    [0, 0, 1]])
        >>> syms = [identity, rot_60]
        >>> 
        >>> hkls, hkls2, fam = generate_hkls02(2, syms, G)
        >>> print("Number of unique reflections:", len(hkls))
```

---

## Function: generate_lattice_faces

**Signature**:
```python
def generate_lattice_faces(uvw2xyz,basal_dirs):
```

**Description**:

Generate face polygons for unit cell visualization.
    Creates vertex lists for 6 faces of unit cell for 3D rendering.

**Input**:

L (array 3×3): Lattice matrix

**Output**:

list: List of 6 face vertex arrays (each 4×3)

**Usage Example**:

```python
>>> import numpy as np
        >>> from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        >>> 
        >>> L = cubic_lattice_vec(3.0)
        >>> faces = generate_lattice_faces(L)
        >>> 
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> poly = Poly3DCollection(faces, alpha=0.2, facecolor='cyan', edgecolor='black')
        >>> ax.add_collection3d(poly)
        >>> set_aspect_equal_3d(ax)
```

**Notes**:

- Returns 6 faces (parallelepiped)
        - Each face is quadrilateral
        - Vertices in counter-clockwise order
        - Used in Poly3DCollection

---

## Function: generate_lattice_points

**Signature**:
```python
def generate_lattice_points(uvw2xyz,basal_dirs):
```

**Description**:

Generate lattice points within specified unit cell range.
    Creates array of lattice points for visualization and analysis.
    Generates n1×n2×n3 unit cells.

**Input**:

L (array 3×3): Lattice matrix [a₁|a₂|a₃]
        n1, n2, n3 (int): Number of unit cells in each direction

**Output**:

numpy.ndarray (N×3): Array of lattice point coordinates

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Generate 2×2×2 cubic lattice
        >>> L_cubic = cubic_lattice_vec(3.0)
        >>> points = generate_lattice_points(L_cubic, n1=2, n2=2, n3=2)
        >>> print(f"Generated {len(points)} lattice points")
        >>> 
        >>> # Visualize
        >>> import matplotlib.pyplot as plt
        >>> from mpl_toolkits.mplot3d import Axes3D
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> ax.scatter(points[:,0], points[:,1], points[:,2])
        >>> plt.show()
```

**Notes**:

- Includes origin (0,0,0)
        - Points at integer lattice coordinates
        - Used in 3D visualization
        - Can visualize multiple unit cells

---

## Function: generate_lattice_vectors

**Signature**:
```python
def generate_lattice_vectors(Points,uvw2xyz,S=1,Q=np.eye(3),xlim=[],ylim=[],zlim=[],fitpoints=False,shift=np.zeros(3)):
```

**Description**:

generate_lattice_vectors - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: generate_lattite_atom_positions

**Signature**:
```python
def generate_lattite_atom_positions(atoms_xyz_position,uvw2xyz,S=1,Q=np.eye(3),R=np.eye(3),shift=np.zeros(3),xlim=[],ylim=[],zlim=[]):
```

**Description**:

generate_lattite_atom_positions - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: generate_plane_vertices

**Signature**:
```python
def generate_plane_vertices(PlanePoints,normal,Q=np.eye(3),move=np.zeros(3)):
```

**Description**:

generate_plane_vertices - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: generate_product_lattice_faces

**Signature**:
```python
def generate_product_lattice_faces(F,Parentlattices):
```

**Description**:

Generate faces for two lattices.
    Creates face polygons for both unit cells.

**Input**:

L1, L2 (array 3×3): Lattice matrices

**Output**:

tuple: (faces1, faces2) - lists of face vertices

**Usage Example**:

```python
>>> L1 = cubic_lattice_vec(3.0)
        >>> L2 = tetragonal_lattice_vec(3.0, 3.0, 4.0)
        >>> faces1, faces2 = generate_product_lattice_faces(L1, L2)
```

**Notes**:

- Returns faces for both lattices
        - Can use different colors/transparency
        - Visualizes structural relationship

---

## Function: generate_product_lattice_points

**Signature**:
```python
def generate_product_lattice_points(F,Parentlattice_points,Q=np.eye(3)):
```

**Description**:

generate_product_lattice_points - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: gensystemsHex

**Signature**:
```python
def gensystemsHex(eta1,K1,L, Lr,sm=None,eta2=None, K2=None):
```

**Description**:

Generate hexagonal slip systems in real space.
    Transforms slip system templates to actual Cartesian coordinates
    using the hexagonal lattice matrix. Produces slip directions and
    plane normals for deformation analysis.

**Input**:

L (array 3×3): Hexagonal lattice matrix [a1|a2|a3]

**Output**:

tuple: (directions, normals) in Cartesian coordinates
            - directions: List of numpy.ndarray [3] - slip directions
            - normals: List of numpy.ndarray [3] - slip plane normals

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # HCP titanium lattice
        >>> a = 2.95
        >>> c = 4.68
        >>> L_hcp = lattice_vec(a, a, c, 90, 90, 120)
        >>> 
        >>> dirs, norms = gensystemsHex(L_hcp)
        >>> print(f"Total slip systems: {len(dirs)}")
        >>> 
        >>> # Calculate Schmid factors
        >>> stress = np.array([1, 0, 0])  # Uniaxial stress
        >>> for i, (d, n) in enumerate(zip(dirs, norms)):
        ...     schmid = np.abs(np.dot(stress, d) * np.dot(stress, n))
        ...     print(f"System {i}: Schmid factor = {schmid:.3f}")
```

**Notes**:

- Uses lattice matrix to transform from fractional to Cartesian
        - Directions: d_cart = L · d_frac
        - Normals: n_cart = L^(-T) · n_frac (reciprocal space)
        - All output vectors are normalized
        - Essential for crystal plasticity simulations

---

## Function: gensystemsHexIni

**Signature**:
```python
def gensystemsHexIni(eta1,K1,L, Lr,sm=None,eta2=None):
```

**Description**:

Generate initial hexagonal slip system templates.
    Creates basal, prismatic, and pyramidal slip system templates for
    hexagonal close-packed (HCP) crystal structures. Returns normalized
    direction and plane normal vectors.

**Input**:

None

**Output**:

tuple: (directions, normals) where each is list of numpy.ndarray [3]
            - directions: Slip directions [uvw]
            - normals: Slip plane normals (hkl)

**Usage Example**:

```python
>>> dirs, norms = gensystemsHexIni()
        >>> print(f"Number of slip systems: {len(dirs)}")
        >>> 
        >>> # Examine first system
        >>> print(f"Direction: {dirs[0]}")
        >>> print(f"Normal: {norms[0]}")
        >>> 
        >>> # Verify normalization
        >>> import numpy as np
        >>> for d in dirs:
        ...     assert np.abs(np.linalg.norm(d) - 1.0) < 1e-10
```

**Notes**:

- Returns normalized vectors (unit length)
        - Includes basal {0001}<11-20> systems
        - Includes prismatic {1-100}<11-20> systems
        - Includes pyramidal systems
        - Used as template for gensystemsHex()

---

## Function: get_interface2d

**Signature**:
```python
def get_interface2d(pointOutproj,normal,horizontalproj,verticalproj):
```

**Description**:

Extract interface atoms from two lattices.
    Finds atoms from both structures near interface plane.
    Used in phase boundary visualization.

**Input**:

atoms1, atoms2 (array N×3): Atom positions for two phases
        plane_normal (array [3]): Interface normal
        plane_point (array [3]): Point on interface
        tolerance (float): Distance tolerance

**Output**:

tuple: (indices1, indices2) - atom indices on interface

**Usage Example**:

```python
>>> # Parent and product phase atoms
        >>> L1 = cubic_lattice_vec(3.0)
        >>> L2 = tetragonal_lattice_vec(3.0, 3.0, 4.0)
        >>> atoms1 = generate_lattite_atom_positions(L1, n1=5, n2=5, n3=5)
        >>> atoms2 = generate_lattite_atom_positions(L2, n1=5, n2=5, n3=5)
        >>> 
        >>> # Get interface atoms
        >>> normal = [0, 0, 1]
        >>> point = [0, 0, 7.5]
        >>> idx1, idx2 = get_interface2d(atoms1, atoms2, normal, point)
        >>> 
        >>> # Plot interface
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> ax.scatter(atoms1[idx1,0], atoms1[idx1,1], atoms1[idx1,2], c='blue')
        >>> ax.scatter(atoms2[idx2,0], atoms2[idx2,1], atoms2[idx2,2], c='red')
```

**Notes**:

- Identifies interface region
        - Both phases analyzed
        - Used in transformation studies

---

## Function: get_twinning_dislocation

**Signature**:
```python
def get_twinning_dislocation(K1,eta1,eta2,L,G=None,Gr=None):
```

**Description**:

Calculate twinning dislocation Burgers vector.
    Determines dislocation content at twin boundary.

**Input**:

twin_data (dict): From get_twinningdata()
        b_parent (array [3]): Burgers vector in parent

**Output**:

dict: {
            'b_residual': array [3] - residual Burgers vector
            'line_direction': array [3] - dislocation line direction
            'type': str - 'edge', 'screw', or 'mixed'
        }

**Usage Example**:

```python
>>> twin_data = niti_twinning('type1')
        >>> b_parent = np.array([1., 0., 0.])  # Parent Burgers vector
        >>> 
        >>> disloc = get_twinning_dislocation(twin_data, b_parent)
        >>> print(f"Residual Burgers vector: {disloc['b_residual']}")
        >>> print(f"Dislocation type: {disloc['type']}")
```

**Notes**:

- Important for twin boundary structure
        - Affects boundary mobility
        - Related to interfacial energy

---

## Function: get_twinning_plane_points

**Signature**:
```python
def get_twinning_plane_points(K1,Pointsout,horizontal,vertical):
```

**Description**:

Get atoms on twinning plane.
    Identifies atoms lying on twin boundary for twinning analysis.

**Input**:

L (array 3×3): Lattice matrix
        twin_plane_normal (array [3]): Twin plane normal (K1)
        n1, n2, n3 (int): Cell range
        tolerance (float): Distance tolerance

**Output**:

numpy.ndarray (N×3): Atoms on twin plane

**Usage Example**:

```python
>>> # Type-I twin in cubic
        >>> L = cubic_lattice_vec(3.0)
        >>> K1 = np.array([1, 1, 1])  # (111) twin plane
        >>> twin_atoms = get_twinning_plane_points(L, K1, n1=10, n2=10, n3=10)
        >>> 
        >>> # Visualize twin boundary
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> ax.scatter(twin_atoms[:,0], twin_atoms[:,1], twin_atoms[:,2], 
        ...            c='red', s=100, label='Twin boundary')
        >>> ax.legend()
```

**Notes**:

- K1 is twinning plane (habit plane)
        - Important for twin boundary structure
        - Used in twinning calculations

---

## Function: get_twinningdata

**Signature**:
```python
def get_twinningdata(orim,eus,Ldir_css,twin_systems,twt,phase, tension=True):
```

**Description**:

Extract complete twinning data from parent and twin lattices.
    Analyzes lattice relationship to determine all twinning elements.

**Input**:

L_parent (array 3×3): Parent lattice matrix
        L_twin (array 3×3): Twin lattice matrix

**Output**:

dict: Complete twin characterization including K1, K2, η1, η2, S

**Usage Example**:

```python
>>> L_parent = monoclinic_lattice_vec(2.889, 4.120, 4.622, 96.8)
        >>> # Twin by compound operation
        >>> R_twin = rotation_from_axis_angle([1,1,1], np.radians(70.5))
        >>> L_twin = R_twin.dot(L_parent)
        >>> 
        >>> data = get_twinningdata(L_parent, L_twin)
        >>> print("Twinning elements:")
        >>> for key, val in data.items():
        ...     print(f"  {key}: {val}")
```

**Notes**:

- Comprehensive twin analysis
        - Extracts all crystallographic elements
        - Validates twin relationship
        - Used in texture simulations

---

## Function: get_unique_families

**Signature**:
```python
def get_unique_families(hkls):
```

**Description**:

Returns unique families of Miller indices based on permutation symmetry.
    Families are considered equivalent if they are permutations of each other
    (considering absolute values).

**Input**:

hkls: list or tuple - List of Miller indices tuples [(h1,k1,l1), (h2,k2,l2), ...]

**Output**:

dict - Dictionary mapping representative hkl to its multiplicity
               {(h,k,l): count, ...}

**Usage Example**:

```python
>>> hkls = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1)]
        >>> families = get_unique_families(hkls)
        >>> print(families)
        {(1, 0, 0): 3, (1, 1, 0): 2}
        >>> # (1,0,0), (0,1,0), (0,0,1) are in the same family with multiplicity 3
        >>> # (1,1,0) and (0,1,1) are in the same family with multiplicity 2
        >>> # Example with negative indices
        >>> hkls2 = [(1, 1, 1), (-1, 1, 1), (1, -1, 1)]
        >>> families2 = get_unique_families(hkls2)
        >>> print(families2)
        {(1, 1, 1): 3}
```

---

## Function: habitplane_equation_solution

**Signature**:
```python
def habitplane_equation_solution(Uj,Ui,Qj,n,a,tol=1e-10):
```

**Description**:

Solve habit plane equation for phase transformation.
    Finds habit plane normal that satisfies invariant plane strain condition:
    F·n = λn + s  (where s is shear direction)

**Input**:

F (array 3×3): Deformation gradient
        s (array [3]): Shear direction

**Output**:

dict: {
            'normal': array [3] - habit plane normal
            'lambda': float - stretch along normal
            'valid': bool - solution exists and is physical
        }

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Martensitic transformation
        >>> F = np.diag([1.05, 0.98, 0.97])  # Deformation gradient
        >>> s = np.array([1, 1, 0]) / np.sqrt(2)  # Shear direction
        >>> 
        >>> result = habitplane_equation_solution(F, s)
        >>> if result['valid']:
        ...     print(f"Habit plane: {result['normal']}")
        ...     print(f"Normal stretch: {result['lambda']:.4f}")
```

**Notes**:

- Central to martensitic transformation theory
        - Habit plane shows no distortion
        - λ close to 1 indicates low-energy interface
        - Used in crystallographic theory of martensite
    Formula:
        (F - λI)·n = s
        Solve for n and λ

---

## Function: hkil2hkl

**Signature**:
```python
def hkil2hkl(hkil):
```

**Description**:

Convert 4-index hexagonal plane (hkil) to 3-index (hkl).
    Transforms Miller-Bravais 4-index plane notation to standard 3-index
    notation for hexagonal crystal systems.

**Input**:

hkil (array [4]): 4-index plane (h, k, i, l) where i = -(h+k)

**Output**:

numpy.ndarray [3]: 3-index plane (h, k, l)

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Convert (11-20) plane
        >>> hkil = np.array([1, 1, -2, 0])
        >>> hkl = hkil2hkl(hkil)
        >>> print(hkl)  # [1, 1, 0]
        >>> 
        >>> # Convert (0001) basal plane
        >>> hkil_basal = np.array([0, 0, 0, 1])
        >>> hkl_basal = hkil2hkl(hkil_basal)
        >>> print(hkl_basal)  # [0, 0, 1]
```

**Notes**:

- Plane notation uses parentheses: (hkl) or (hkil)
        - i is redundant: i = -(h + k)
        - Conversion: drop i, keep h,k,l
        - Inverse operation of hkl2hkil()
    Formula:
        (h, k, l) = (h₄, k₄, l₄) from (h₄, k₄, i₄, l₄)

---

## Function: hkl2hkil

**Signature**:
```python
def hkl2hkil(hkl):
```

**Description**:

Convert 3-index hexagonal plane (hkl) to 4-index (hkil).
    Transforms standard 3-index plane notation to Miller-Bravais 4-index
    notation for hexagonal crystal systems. The redundant index i = -(h+k).

**Input**:

hkl (array [3]): 3-index plane (h, k, l)

**Output**:

numpy.ndarray [4]: 4-index plane (h, k, i, l) where i = -(h+k)

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Convert (110) plane
        >>> hkl = np.array([1, 1, 0])
        >>> hkil = hkl2hkil(hkl)
        >>> print(hkil)  # [1, 1, -2, 0]
        >>> 
        >>> # Convert (001) plane
        >>> hkl_001 = np.array([0, 0, 1])
        >>> hkil_001 = hkl2hkil(hkl_001)
        >>> print(hkil_001)  # [0, 0, 0, 1]
```

**Notes**:

- 4-index emphasizes hexagonal symmetry
        - i is calculated: i = -(h + k)
        - Inverse operation of hkil2hkl()
        - Standard in hexagonal crystallography
    Formula:
        i = -(h + k)
        (h, k, i, l) = (h, k, -(h+k), l)

---

## Function: kronecker

**Signature**:
```python
def kronecker():
```

**Description**:

Compute Kronecker delta δ_{ij}.
    Returns 1 if indices are equal, 0 otherwise. Fundamental in
    tensor algebra and represents the identity tensor.

**Input**:

i, j (int): Indices

**Output**:

int: 1 if i==j, 0 otherwise

**Usage Example**:

```python
>>> # Identity matrix via Kronecker delta
        >>> import numpy as np
        >>> I = np.array([[kronecker(i,j) for j in range(3)] for i in range(3)])
        >>> print(I)
        >>> # [[1, 0, 0],
        >>> #  [0, 1, 0],
        >>> #  [0, 0, 1]]
        >>> 
        >>> # Check orthogonality
        >>> for i in range(3):
        ...     for j in range(3):
        ...         print(f"δ_{i}{j} = {kronecker(i,j)}")
```

**Notes**:

- δ_{ij} = 1 if i = j
        - δ_{ij} = 0 if i ≠ j
        - Contraction: A_{ij} δ_{jk} = A_{ik}
        - Trace: δ_{ii} = 3 (in 3D)
    Formula:
        δ_{ij} = { 1  if i = j
                 { 0  if i ≠ j

---

## Function: lattice_correspondence

**Signature**:
```python
def lattice_correspondence(LatCorr,parent_symops,product_symops):
```

**Description**:

Calculate lattice correspondence matrix between two crystal structures.
    Finds transformation matrix C relating two lattices: L2 = C · L1.
    Optionally optimizes to minimize lattice mismatch.

**Input**:

L1 (array 3×3): First lattice matrix
        L2 (array 3×3): Second lattice matrix
        optimize (bool): If True, optimize correspondence (default: False)

**Output**:

numpy.ndarray (3×3): Correspondence matrix C

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Cubic to tetragonal
        >>> L_cubic = cubic_lattice_vec(3.0)
        >>> L_tetra = tetragonal_lattice_vec(3.0, 3.0, 4.0)
        >>> 
        >>> C = lattice_correspondence(L_cubic, L_tetra)
        >>> print("Correspondence matrix:")
        >>> print(C)
        >>> 
        >>> # Verify: L_tetra ≈ C · L_cubic
        >>> L_check = C.dot(L_cubic)
        >>> print("Reconstruction error:", np.linalg.norm(L_check - L_tetra))
```

**Notes**:

- Simple version: C = L2 · L1^(-1)
        - Optimized version minimizes ||L2 - C·L1||
        - Used in phase transformation analysis
        - Essential for orientation relationships
    Formula:
        C = L2 · L1^(-1)

---

## Function: lattice_vec

**Signature**:
```python
def lattice_vec(lattice_param):
```

**Description**:

Calculate lattice vectors (a1, a2, a3) for various crystal systems.

**Input**:

lattice_param: dict - Dictionary containing lattice parameters
            Required keys depend on crystal system type:
            For 'cubic':
                'type': str - 'cubic'
                'a': float - Lattice parameter (Å)
            For 'tetragonal':
                'type': str - 'tetragonal'
                'a': float - a lattice parameter (Å)
                'b': float - b lattice parameter (Å)
                'c': float - c lattice parameter (Å)
            For 'monoclinic':
                'type': str - 'monoclinic'
                'a', 'b', 'c': float - Lattice parameters (Å)
                'beta': float - Beta angle (radians)
            For 'triclinic':
                'type': str - 'triclinic'
                'a', 'b', 'c': float - Lattice parameters (Å)
                'alpha', 'beta', 'gamma': float - Angles (radians)
            For 'trigonal':
                'type': str - 'trigonal'
                'a': float - a lattice parameter (Å)
                'c': float - c lattice parameter (Å)

**Output**:

a1: numpy array (3,) - First lattice vector
        a2: numpy array (3,) - Second lattice vector
        a3: numpy array (3,) - Third lattice vector

**Usage Example**:

```python
>>> # Cubic crystal (e.g., Silicon)
        >>> cubic_params = {'type': 'cubic', 'a': 5.43}
        >>> a1, a2, a3 = lattice_vec(cubic_params)
        >>> print("a1:", a1)
        >>> print("a2:", a2)
        >>> print("a3:", a3)
        a1: [5.43 0.   0.  ]
        a2: [0.   5.43 0.  ]
        a3: [0.   0.   5.43]
        >>> # Hexagonal crystal (trigonal setting)
        >>> hex_params = {'type': 'trigonal', 'a': 3.0, 'c': 5.0}
        >>> a1, a2, a3 = lattice_vec(hex_params)
        >>> print("a1:", a1)
        >>> print("a3:", a3)
        >>> # Monoclinic crystal
        >>> import numpy as np
        >>> mono_params = {
        ...     'type': 'monoclinic',
        ...     'a': 5.0, 'b': 6.0, 'c': 7.0,
        ...     'beta': np.radians(120)  # 120 degrees
        ... }
        >>> a1, a2, a3 = lattice_vec(mono_params)
        >>> print("Volume:", np.dot(a1, np.cross(a2, a3)))
```

---

## Function: miller2fractional

**Signature**:
```python
def miller2fractional(uvw,frac=10,eps2=1e-2,decimals=5):
```

**Description**:

Reduce Miller indices to lowest integer form.
    Takes potentially fractional Miller indices and reduces them to the
    simplest integer representation by finding rational approximations
    and applying GCD reduction.

**Input**:

uvw (array [3]): Miller indices (can be fractional)
        frac (int, optional): Max denominator for fractions (default: 10)
        eps2 (float, optional): Tolerance for rounding (default: 1e-2)
        decimals (int, optional): Rounding precision (default: 5)

**Output**:

numpy.ndarray [3]: Reduced integer Miller indices

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Reduce fractional indices
        >>> uvw = np.array([0.5, 1.0, 1.5])
        >>> reduced = miller2fractional(uvw)
        >>> print(reduced)  # Output: [1, 2, 3]
        >>> 
        >>> # Already integer indices
        >>> uvw2 = np.array([2, 4, 6])
        >>> reduced2 = miller2fractional(uvw2)
        >>> print(reduced2)  # Output: [1, 2, 3]
        >>> 
        >>> # Handle complex fractions
        >>> uvw3 = np.array([0.333333, 0.666667, 1.0])
        >>> reduced3 = miller2fractional(uvw3, frac=15)
        >>> print(reduced3)  # Output: [1, 2, 3]
```

**Notes**:

- Uses sympy.nsimplify() for rational approximation
        - Applies GCD to reduce to lowest terms
        - Parameter 'frac' sets maximum allowed denominator
        - Higher 'frac' allows more precise fractional representation
        - Used internally by xyz2fractional()
    Algorithm:
        1. Round values close to integers
        2. For fractional values: find rational approximation
        3. Clear denominators by multiplication
        4. Apply GCD to reduce

---

## Function: mohr_circles

**Signature**:
```python
def mohr_circles(tensor):
```

**Description**:

Calculate Mohr's circles from strain or stress tensor.
    Computes principal strains/stresses and parameters for three Mohr's
    circles. Used for visualizing 3D strain state in 2D.

**Input**:

strain_tensor (array 3×3): Symmetric strain or stress tensor

**Output**:

dict: {
            'principal': array [3] - principal values (sorted)
            'directions': array 3×3 - principal directions
            'circles': list of 3 tuples - (center, radius) for each circle
        }

**Usage Example**:

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

**Notes**:

- Three circles for 3D state
        - Largest circle: (ε₁ - ε₃)/2
        - Used in failure analysis
        - Shows all possible strain states on planes
    Formula:
        Circle i: center = (σᵢ + σⱼ)/2, radius = |σᵢ - σⱼ|/2
        Three circles: 1-2, 2-3, 1-3 planes

---

## Function: monoclinic_lattice_vec

**Signature**:
```python
def monoclinic_lattice_vec(a,b,c,beta):
```

**Description**:

Generate monoclinic lattice vectors (B19' martensite).
    Creates 3×3 lattice matrix for monoclinic crystal system.
    One unique angle (β) differs from 90°, typical of B19' martensite in NiTi.

**Input**:

a (float): Lattice parameter a (Ångströms)
        b (float): Lattice parameter b (Ångströms)
        c (float): Lattice parameter c (Ångströms)
        beta (float): Angle β between a and c axes (degrees)

**Output**:

numpy.ndarray (3×3): Monoclinic lattice matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # B19' martensite (NiTi)
        >>> LM = monoclinic_lattice_vec(2.889, 4.120, 4.622, 96.8)
        >>> print(LM)
        >>> 
        >>> # Calculate volume
        >>> volume = np.linalg.det(LM)
        >>> print(f"Unit cell volume: {volume:.3f} ų")
        >>> 
        >>> # Verify angles
        >>> a_vec, c_vec = LM[:, 0], LM[:, 2]
        >>> angle = np.degrees(np.arccos(np.dot(a_vec, c_vec) / 
        ...                               (np.linalg.norm(a_vec) * np.linalg.norm(c_vec))))
        >>> print(f"β angle: {angle:.1f}°")
```

**Notes**:

- Monoclinic system: a≠b≠c, α=γ=90°, β≠90°
        - Common in martensitic phases
        - B19' is monoclinic martensite in NiTi
        - Volume = abc·sin(β)
    Formula:
        L = [[a,           0, c·cos(β)],
             [0,           b, 0       ],
             [0,           0, c·sin(β)]]

---

## Function: niti_twinning

**Signature**:
```python
def niti_twinning(B2_symops,B2_recsymops,B19p_recsymops,B19p_symops,Uv,Parent_uvw2xyz,Parent_hkl2xyz,Product_uvw2xyz, Product_hkl2xyz,Parent_uvw_2_Product_uvw_all, Parent_hkl_2_Product_hkl_all, Parent_uvw_2_Product_uvw_all_norm,miller='greaterthanone',Qv=None):
```

**Description**:

Get NiTi twinning parameters for specific variant.
    Returns crystallographic twinning elements for NiTi B19' martensite.

**Input**:

variant (str): Twinning variant ('type1', 'type2', 'compound')

**Output**:

dict: {
            'K1': array [3] - compound twin plane
            'eta1': array [3] - twinning direction
            'K2': array [3] - conjugate plane
            'eta2': array [3] - conjugate direction
            'shear': float - twinning shear
        }

**Usage Example**:

```python
>>> # Type-I twin (most common in NiTi)
        >>> twin_data = niti_twinning('type1')
        >>> print(f"Twin plane K1: {twin_data['K1']}")
        >>> print(f"Twin direction η1: {twin_data['eta1']}")
        >>> print(f"Shear magnitude: {twin_data['shear']:.4f}")
        >>> 
        >>> # Use in deformation analysis
        >>> F_twin = np.eye(3) + twin_data['shear'] * np.outer(twin_data['eta1'], twin_data['K1'])
```

**Notes**:

- Specific to NiTi B19' martensite
        - Type-I: {111} compound twin
        - Type-II: {001} twin
        - Data from experimental measurements

---

## Function: normArrayColumns

**Signature**:
```python
def normArrayColumns(arr):
```

**Description**:

Normalize each column of matrix to unit length.
    Divides each column by its Euclidean norm, creating an orthonormal
    or semi-orthonormal matrix. Commonly used in crystallography for
    basis vector normalization.

**Input**:

arr (array 3×3): Input matrix with column vectors

**Output**:

numpy.ndarray (3×3): Matrix with unit-length columns

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Normalize lattice vectors
        >>> lattice = np.array([[3, 0, 0],
        ...                     [0, 4, 0],
        ...                     [0, 0, 5]], dtype=float)
        >>> normalized = normArrayColumns(lattice)
        >>> 
        >>> # Verify unit norms
        >>> for i in range(3):
        ...     norm_i = np.linalg.norm(normalized[:, i])
        ...     print(f"Column {i} norm: {norm_i:.6f}")  # All 1.0
        >>> 
        >>> # Create orthonormal basis for calculations
        >>> basis = np.random.rand(3, 3)
        >>> basis_normalized = normArrayColumns(basis)
```

**Notes**:

- Preserves column directions, only changes magnitudes
        - Each column becomes a unit vector
        - Does not orthogonalize - only normalizes
        - Useful for direction cosine matrices
        - Common preprocessing step in texture analysis
    Formula:
        For column i: normalized[:, i] = arr[:, i] / ||arr[:, i]||

---

## Function: np_kronecker

**Signature**:
```python
def np_kronecker():
```

**Description**:

Compute Kronecker delta (NumPy version).
    NumPy-compatible version of kronecker().

**Input**:

i, j (int): Indices

**Output**:

int: Kronecker delta value

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Create identity matrix
        >>> n = 3
        >>> I = np.fromfunction(lambda i,j: np_kronecker(int(i), int(j)), (n,n))
        >>> print(I)
```

**Notes**:

- Same as kronecker()
        - Compatible with NumPy vectorization

---

## Function: np_permut_tensor3

**Signature**:
```python
def np_permut_tensor3():
```

**Description**:

Compute Levi-Civita symbol (NumPy-compatible version).
    Same as permut_tensor3 but optimized for NumPy array indexing.

**Input**:

i, j, k (int): Tensor indices

**Output**:

int: Permutation symbol value

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Build full permutation tensor
        >>> epsilon = np.zeros((3,3,3))
        >>> for i in range(3):
        ...     for j in range(3):
        ...         for k in range(3):
        ...             epsilon[i,j,k] = np_permut_tensor3(i,j,k)
        >>> 
        >>> # Verify antisymmetry
        >>> assert epsilon[0,1,2] == -epsilon[0,2,1]
```

**Notes**:

- Identical to permut_tensor3
        - NumPy-friendly naming convention

---

## Function: permut_tensor3

**Signature**:
```python
def permut_tensor3():
```

**Description**:

Compute 3D Levi-Civita permutation symbol ε_{ijk}.
    Returns the permutation tensor (Levi-Civita symbol) for three indices.
    Used in cross products, curl operations, and tensor calculations.

**Input**:

i, j, k (int): Indices (0, 1, or 2 representing x, y, z)

**Output**:

int: +1 for even permutation, -1 for odd, 0 if any indices repeat

**Usage Example**:

```python
>>> # Even permutations
        >>> print(permut_tensor3(0, 1, 2))  # +1 (xyz)
        >>> print(permut_tensor3(1, 2, 0))  # +1 (yzx)
        >>> print(permut_tensor3(2, 0, 1))  # +1 (zxy)
        >>> 
        >>> # Odd permutations  
        >>> print(permut_tensor3(0, 2, 1))  # -1 (xzy)
        >>> print(permut_tensor3(2, 1, 0))  # -1 (zyx)
        >>> 
        >>> # Repeated indices
        >>> print(permut_tensor3(0, 0, 1))  # 0
        >>> print(permut_tensor3(1, 1, 2))  # 0
```

**Notes**:

- ε_{ijk} = +1 if (i,j,k) is even permutation of (0,1,2)
        - ε_{ijk} = -1 if (i,j,k) is odd permutation of (0,1,2)
        - ε_{ijk} = 0 if any two indices are equal
        - Used in vector cross product: (a×b)_i = ε_{ijk} a_j b_k
    Formula:
        ε_{012} = ε_{120} = ε_{201} = +1
        ε_{021} = ε_{210} = ε_{102} = -1
        ε_{ijk} = 0 otherwise

---

## Function: perpendicular_vector

**Signature**:
```python
def perpendicular_vector(v):
```

**Description**:

Find arbitrary unit vector perpendicular to input 3D vector.
    Computes a normalized vector perpendicular to the input using cross product
    with a judiciously chosen auxiliary vector. Used in crystallographic
    calculations requiring orthogonal basis construction.

**Input**:

v (array-like [3]): Input 3D vector (must be non-zero)

**Output**:

numpy.ndarray [3]: Unit vector perpendicular to v (norm = 1.0)

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Find perpendicular to [1,2,3]
        >>> v = np.array([1, 2, 3])
        >>> v_perp = perpendicular_vector(v)
        >>> print(f"Dot product: {np.dot(v, v_perp):.10f}")  # ~0 (perpendicular)
        >>> print(f"Norm: {np.linalg.norm(v_perp):.10f}")  # 1.0 (unit vector)
        >>> 
        >>> # Use for constructing orthonormal basis
        >>> v1 = np.array([1, 0, 0])
        >>> v2 = perpendicular_vector(v1)
        >>> v3 = np.cross(v1, v2)  # Third orthogonal vector
        >>> # Now v1, v2, v3 form orthonormal basis
```

**Notes**:

- Result is arbitrary (many perpendicular vectors exist)
        - Uses cross product with [1,0,0] or [0,1,0] depending on input
        - Automatically normalized to unit length
        - Raises error if input vector is zero
        - Common in stereographic projection calculations
    Algorithm:
        1. Choose auxiliary vector that isn't parallel to input
        2. Compute cross product
        3. Normalize to unit length

---

## Function: plane2string

**Signature**:
```python
def plane2string(v, digits=2):
```

**Description**:

Format plane normal as string representation (h,k,l) with Miller notation.
    Converts plane normal vector to crystallographic notation using
    parentheses, following standard Miller index convention.

**Input**:

v (array-like [3]): Plane normal vector (h, k, l)
        digits (int, optional): Number of decimal places (default: 2)

**Output**:

str: Formatted string '(h.hh,k.kk,l.ll)'

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Format {111} plane
        >>> hkl = [1, 1, 1]
        >>> print(plane2string(hkl))  # Output: '(1.00,1.00,1.00)'
        >>> 
        >>> # Format with fractional indices
        >>> plane = np.array([0.5, 1.0, 1.5])
        >>> print(plane2string(plane, digits=2))  # Output: '(0.50,1.00,1.50)'
        >>> 
        >>> # Combine with Miller indices calculation
        >>> lattice_vec = cubic_lattice_vec(3.0)
        >>> normal = [1, 1, 0]
        >>> print(f"Plane: {plane2string(normal)}")
```

**Notes**:

- Uses parentheses () for plane notation (Miller index convention)
        - Compare with dir2string() which uses square brackets []
        - Standard in crystallography: (hkl) for planes, [uvw] for directions
        - Can handle negative indices

---

## Function: plane_line_intersection

**Signature**:
```python
def plane_line_intersection(n,V0,P0,P1):
```

---

## Function: print_correspondence

**Signature**:
```python
def print_correspondence(Mcorr,VecA,latticeA, latticeB,planes=False,returnB=False):
```

**Description**:

Print lattice correspondence matrix in readable format.
    Displays correspondence matrix with phase labels for documentation
    and reporting purposes.

**Input**:

C (array 3×3): Correspondence matrix
        phase1 (str): Name of first phase (default: 'Phase1')
        phase2 (str): Name of second phase (default: 'Phase2')

**Output**:

None (prints to console)

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> C = B19p_B2_lattice_correspondence()
        >>> print_correspondence(C, 'B19 martensite', 'B2 austenite')
        >>> 
        >>> # Output formatted as:
        >>> # B19' martensite → B2 austenite correspondence:
        >>> # [[ 1.000  0.000  0.000]
        >>> #  [ 0.000  1.000  0.000]
        >>> #  [ 0.000  0.000  1.000]]
```

**Notes**:

- Formatted for readability
        - Useful in reports and documentation
        - Shows transformation relationship clearly

---

## Function: read_txt

**Signature**:
```python
def read_txt(filename,delimiter='\t',skiprows=1):
```

**Description**:

Read numerical data from text file.
    Loads data written by write_txt() or similar format.

**Input**:

filename (str): Input file path

**Output**:

numpy.ndarray: Loaded data

**Usage Example**:

```python
>>> # Read previously saved data
        >>> data = read_txt('orientations.txt')
        >>> print(f"Loaded {len(data)} orientations")
        >>> print(f"Shape: {data.shape}")
        >>> 
        >>> # Use data
        >>> euler_angles = data[:, :3]
        >>> # Process...
```

**Notes**:

- Skips comment lines starting with #
        - Returns numpy array
        - Compatible with write_txt output

---

## Function: reciprocal_basis

**Signature**:
```python
def reciprocal_basis(a1,a2,a3):
```

**Description**:

Calculate reciprocal lattice basis vectors from real space lattice vectors.
    The reciprocal lattice vectors satisfy: a_i · b_j = 2π δ_ij
    (Note: This implementation uses the crystallographic convention without 2π factor)

**Input**:

a1: numpy array (3,) - First real space lattice vector
        a2: numpy array (3,) - Second real space lattice vector
        a3: numpy array (3,) - Third real space lattice vector

**Output**:

b1: numpy array (3,) - First reciprocal lattice vector
        b2: numpy array (3,) - Second reciprocal lattice vector
        b3: numpy array (3,) - Third reciprocal lattice vector

**Usage Example**:

```python
>>> import numpy as np
        >>> # Define cubic lattice
        >>> a = 5.0  # Angstroms
        >>> a1 = np.array([a, 0, 0])
        >>> a2 = np.array([0, a, 0])
        >>> a3 = np.array([0, 0, a])
        >>> 
        >>> # Calculate reciprocal vectors
        >>> b1, b2, b3 = reciprocal_basis(a1, a2, a3)
        >>> print("b1:", b1)
        >>> print("b2:", b2)
        >>> print("b3:", b3)
        b1: [0.2 0.  0. ]
        b2: [0.  0.2 0. ]
        b3: [0.  0.  0.2]
        >>> # Verify orthogonality condition
        >>> print("a1 · b1 =", np.dot(a1, b1))  # Should be 1.0
        >>> print("a1 · b2 =", np.dot(a1, b2))  # Should be 0.0
        >>> # Example with hexagonal lattice
        >>> hex_params = {'type': 'trigonal', 'a': 3.0, 'c': 5.0}
        >>> a1_hex, a2_hex, a3_hex = lattice_vec(hex_params)
        >>> b1_hex, b2_hex, b3_hex = reciprocal_basis(a1_hex, a2_hex, a3_hex)
        >>> print("Reciprocal c* length:", np.linalg.norm(b3_hex))
```

---

## Function: select_atomic_plane

**Signature**:
```python
def select_atomic_plane(LatticePoints,normal,eps=1e-1,shift=0.,eps2=None):
```

**Description**:

Select atoms lying on or near a crystallographic plane.
    Finds atoms within tolerance distance from specified plane.

**Input**:

atoms (array N×3): Atomic positions
        plane_normal (array [3]): Plane normal vector
        plane_point (array [3]): Point on plane
        tolerance (float): Distance tolerance (Ångströms)

**Output**:

numpy.ndarray: Indices of atoms on plane

**Usage Example**:

```python
>>> L = cubic_lattice_vec(3.0)
        >>> atoms = generate_lattite_atom_positions(L, n1=3, n2=3, n3=3)
        >>> 
        >>> # Select atoms on (001) plane at z=6
        >>> plane_normal = [0, 0, 1]
        >>> plane_point = [0, 0, 6]
        >>> indices = select_atomic_plane(atoms, plane_normal, plane_point, tolerance=0.1)
        >>> 
        >>> print(f"Found {len(indices)} atoms on plane")
        >>> plane_atoms = atoms[indices]
```

**Notes**:

- Uses point-to-plane distance formula
        - Tolerance accounts for numerical error
        - Returns indices, not coordinates
    Formula:
        distance = |n · (p - p₀)| / ||n||

---

## Function: select_atomic_region

**Signature**:
```python
def select_atomic_region(LatticePoints,normal,side='lower',eps=1e-1,shift=0.):
```

**Description**:

Select atoms within spherical region.
    Finds all atoms within specified radius of center point.

**Input**:

atoms (array N×3): Atomic positions
        center (array [3]): Sphere center
        radius (float): Sphere radius

**Output**:

numpy.ndarray: Indices of atoms in region

**Usage Example**:

```python
>>> atoms = generate_lattite_atom_positions(L, n1=10, n2=10, n3=10)
        >>> 
        >>> # Select atoms near point [15, 15, 15]
        >>> center = np.array([15, 15, 15])
        >>> indices = select_atomic_region(atoms, center, radius=5.0)
        >>> 
        >>> print(f"Found {len(indices)} atoms in region")
        >>> region_atoms = atoms[indices]
```

**Notes**:

- Spherical selection
        - Useful for local analysis
        - Can be used for grain selection
    Formula:
        ||atom - center|| < radius

---

## Function: select_crystal_planes

**Signature**:
```python
def select_crystal_planes(NormalsOnCircle,ShearsOnCircle,maxhkl):
```

**Description**:

Generate list of low-index crystallographic planes.
    Creates planes (hkl) with indices up to max_index for analysis.

**Input**:

lattice (array 3×3): Lattice matrix
        max_index (int): Maximum Miller index value

**Output**:

list: List of (h,k,l) tuples

**Usage Example**:

```python
>>> L = cubic_lattice_vec(3.0)
        >>> planes = select_crystal_planes(L, max_index=2)
        >>> 
        >>> print(f"Generated {len(planes)} planes")
        >>> for hkl in planes[:10]:
        ...     print(f"({hkl[0]},{hkl[1]},{hkl[2]})")
```

**Notes**:

- Low-index planes are physically important
        - Used in diffraction analysis
        - Filters out equivalent planes

---

## Function: select_plane

**Signature**:
```python
def select_plane(LatticeVectors,normal,eps=1e-1,shift=0.,Q=np.eye(3)):
```

**Description**:

select_plane - Crystallographic function for materials analysis.
    See full documentation in extended modules for detailed usage.

---

## Function: strains_along_13mohrcirle

**Signature**:
```python
def strains_along_13mohrcirle(Strain,VV,normdiri,phi_around_V2,Parent_xyz2hkl):
```

**Description**:

Calculate strains along maximum Mohr's circle.
    Computes normal and shear strain components for points on the
    largest Mohr's circle (ε₁-ε₃ circle).

**Input**:

mohr_result (dict): Output from mohr_circles()
        n_points (int): Number of points to sample (default: 100)

**Output**:

dict: {
            'normal': array - normal strain values
            'shear': array - shear strain values
            'angles': array - rotation angles (radians)
        }

**Usage Example**:

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

**Notes**:

- Maximum circle between ε₁ and ε₃
        - Shows all possible strain states
        - Used in failure analysis
    Formula:
        εₙ = center + radius·cos(2θ)
        γ/2 = radius·sin(2θ)

---

## Function: tetragonal_lattice_vec

**Signature**:
```python
def tetragonal_lattice_vec(a,b,c):
```

**Description**:

Generate tetragonal lattice vectors.
    Creates 3×3 lattice matrix for tetragonal crystal system.
    Two equal parameters (a=b) and one unique (c), all angles 90°.

**Input**:

a (float): Lattice parameter a (Ångströms)
        b (float): Lattice parameter b (Ångströms, typically b=a)
        c (float): Lattice parameter c (Ångströms)

**Output**:

numpy.ndarray (3×3): Tetragonal lattice matrix

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Tetragonal system with a=b=3.0, c=4.0
        >>> L_tet = tetragonal_lattice_vec(3.0, 3.0, 4.0)
        >>> print(L_tet)
        >>> # [[3, 0, 0],
        >>> #  [0, 3, 0],
        >>> #  [0, 0, 4]]
        >>> 
        >>> # Calculate c/a ratio
        >>> c_over_a = c / a
        >>> print(f"c/a ratio: {c_over_a:.3f}")
```

**Notes**:

- Tetragonal system: a=b≠c, α=β=γ=90°
        - Special case of orthorhombic
        - Volume = a²c
        - Common in oxide ceramics
    Formula:
        L = [[a, 0, 0],
             [0, b, 0],
             [0, 0, c]]

---

## Function: twin_equation_solution

**Signature**:
```python
def twin_equation_solution(Uj,Ui,L_A,Lr_A,L_M,Lr_M, R_AM, Ci_d,Ci_p,tol=1e-10,miller='greaterthanone',printlambda=False, Qj=None,Qi=None):
```

**Description**:

Solve complete twin equation system.
    Finds K1, K2, η1, η2 twin elements from lattice parameters.

**Input**:

L_parent (array 3×3): Parent lattice
        L_twin (array 3×3): Twin lattice
        params (dict): Solver parameters (from twin_equation_solution_ini)

**Output**:

dict: {
            'K1': array [3] - composition/habit plane
            'K2': array [3] - conjugate plane
            'eta1': array [3] - shear direction in K1
            'eta2': array [3] - direction in K2
            'P': array [3] - invariant line
            'S': float - twinning shear magnitude
        }

**Usage Example**:

```python
>>> L_parent = cubic_lattice_vec(3.0)
        >>> # Create twinned lattice
        >>> theta = np.radians(70.5)  # Twin rotation angle
        >>> R_twin = rotation_from_axis_angle([1,1,1], theta)
        >>> L_twin = R_twin.dot(L_parent)
        >>> 
        >>> result = twin_equation_solution(L_parent, L_twin)
        >>> print("Twin elements:")
        >>> print(f"  K1 (twin plane): {result['K1']}")
        >>> print(f"  η1 (twin direction): {result['eta1']}")
        >>> print(f"  Shear magnitude: {result['S']:.4f}")
```

**Notes**:

- Complete crystallographic solution
        - Includes both K1/η1 and K2/η2 pairs
        - P is invariant line (no rotation or stretch)
        - Used in comprehensive twin analysis
    Formula:
        F = I + S(η1 ⊗ K1)
        where F is deformation gradient, S is shear

---

## Function: twin_equation_solution_ini

**Signature**:
```python
def twin_equation_solution_ini(Uj,Ui,Parent_uvw2xyz,Parent_hkl2xyz,Product_uvw2xyz,Product_hkl2xyz, Parent_uvw_2_Product_uvw_rot, Parent_uvw_2_Product_uvw,Parent_hkl_2_Product_hkl,tol=1e-10,miller='greaterthanone',printlambda=False, Qj=None,Qi=None):
```

**Description**:

Initialize twin equation solver with default parameters.
    Sets up initial guess and parameters for twin equation solution.

**Input**:

None

**Output**:

dict: Initial parameters for twin solver

**Usage Example**:

```python
>>> params = twin_equation_solution_ini()
        >>> # Modify parameters as needed
        >>> params['max_iter'] = 1000
        >>> params['tolerance'] = 1e-8
        >>> # Use in twin_equation_solution()
```

**Notes**:

- Provides sensible defaults
        - Can be customized for specific systems
        - Used by twin_equation_solution()

---

## Function: twinnedhabitplane

**Signature**:
```python
def twinnedhabitplane(Ui,Uj,Qij,a1,n1,hbplanes=[],addondata={},method='bhata'):
```

**Description**:

Calculate twinned habit plane from lattice parameters.
    Determines twin plane (K1) and twin direction from parent and
    twin lattice geometries.

**Input**:

L_parent (array 3×3): Parent lattice matrix
        L_twin (array 3×3): Twin lattice matrix  
        correspondence_matrix (array 3×3): Parent→twin correspondence

**Output**:

dict: {
            'K1': array [3] - twin plane (habit plane)
            'eta1': array [3] - twin direction
            'shear': float - twinning shear magnitude
        }

**Usage Example**:

```python
>>> # Type-I twin in cubic
        >>> L_parent = cubic_lattice_vec(3.0)
        >>> # Twin related by mirror symmetry
        >>> C_twin = np.diag([1, 1, -1])  # Mirror in (001)
        >>> L_twin = C_twin.dot(L_parent)
        >>> 
        >>> result = twinnedhabitplane(L_parent, L_twin, C_twin)
        >>> print(f"Twin plane K1: {result['K1']}")
        >>> print(f"Twin direction η1: {result['eta1']}")
```

**Notes**:

- K1 is rational in parent
        - η1 is twin shear direction  
        - Used in deformation twinning analysis
        - Applies to mechanical twins

---

## Function: uvtw2uvw

**Signature**:
```python
def uvtw2uvw(uvtw):
```

**Description**:

Convert 4-index hexagonal direction [uvtw] to 3-index [uvw].
    Transforms Miller-Bravais 4-index notation to standard 3-index notation
    for hexagonal crystal systems. The redundant index t is removed.

**Input**:

uvtw (array [4]): 4-index direction [u, v, t, w] where t = -(u+v)

**Output**:

numpy.ndarray [3]: 3-index direction [u, v, w]

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Convert [11-20] direction
        >>> uvtw = np.array([1, 1, -2, 0])
        >>> uvw = uvtw2uvw(uvtw)
        >>> print(uvw)  # [1, 1, 0]
        >>> 
        >>> # Verify t = -(u+v)
        >>> u, v, t, w = uvtw
        >>> print(f"t = -(u+v): {t} = -{(u+v)}")  # -2 = -2
```

**Notes**:

- Used in hexagonal systems (HCP, graphite, etc.)
        - 4-index notation makes symmetry more apparent
        - Redundancy: t = -(u + v)
        - Conversion is straightforward: drop t, keep u,v,w
    Formula:
        [u, v, w] = [u₄, v₄, w₄] from [u₄, v₄, t₄, w₄]

---

## Function: uvw2uvtw

**Signature**:
```python
def uvw2uvtw(uvw):
```

**Description**:

Convert 3-index hexagonal direction [uvw] to 4-index [uvtw].
    Transforms standard 3-index notation to Miller-Bravais 4-index notation
    for hexagonal crystal systems. The redundant index t = -(u+v) is inserted.

**Input**:

uvw (array [3]): 3-index direction [u, v, w]

**Output**:

numpy.ndarray [4]: 4-index direction [u, v, t, w] where t = -(u+v)

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Convert [110] direction
        >>> uvw = np.array([1, 1, 0])
        >>> uvtw = uvw2uvtw(uvw)
        >>> print(uvtw)  # [1, 1, -2, 0]
        >>> 
        >>> # Convert [100] direction
        >>> uvw2 = np.array([1, 0, 0])
        >>> uvtw2 = uvw2uvtw(uvw2)
        >>> print(uvtw2)  # [1, 0, -1, 0]
```

**Notes**:

- 4-index notation emphasizes hexagonal symmetry
        - t is calculated, not independent
        - Inverse operation of uvtw2uvw()
        - Common in hexagonal close-packed (HCP) systems
    Formula:
        t = -(u + v)
        [u, v, t, w] = [u, v, -(u+v), w]

---

## Function: vec2string

**Signature**:
```python
def vec2string(v, digits=2):
```

**Description**:

Format vector as string representation [x,y,z] with specified precision.
    Converts numpy array or list to readable string format suitable for
    display, logging, or file output in crystallographic applications.

**Input**:

v (array-like [3]): Vector to format
        digits (int, optional): Number of decimal places (default: 2)

**Output**:

str: Formatted string '[x.xx,y.yy,z.zz]'

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Format with default 2 decimal places
        >>> v = [1.23456, 0.56789, 0.89123]
        >>> print(vec2string(v))  # Output: '[1.23,0.57,0.89]'
        >>> 
        >>> # Format with 4 decimal places
        >>> uvw = np.array([1.0, 1.4142, 1.7321])
        >>> print(vec2string(uvw, digits=4))  # Output: '[1.0000,1.4142,1.7321]'
        >>> 
        >>> # Use in loop for multiple vectors
        >>> directions = [[1,0,0], [0,1,0], [0,0,1]]
        >>> for d in directions:
        ...     print(f"Direction: {vec2string(d)}")
```

**Notes**:

- Uses square brackets [] for vector notation
        - Comma-separated without spaces
        - Consistent formatting across module
        - Related: plane2string(), dir2string()

---

## Function: vector2miller_ini

**Signature**:
```python
def vector2miller_ini(v, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3):
```

**Description**:

Convert Cartesian vector to Miller indices (initial guess).
    Provides starting point for vector→Miller conversion.

**Input**:

v (array [3]): Cartesian vector
        L (array 3×3): Lattice matrix

**Output**:

numpy.ndarray [3]: Initial Miller index guess

**Usage Example**:

```python
>>> L = cubic_lattice_vec(3.0)
        >>> v_cart = np.array([3, 3, 0])
        >>> uvw_guess = vector2miller_ini(v_cart, L)
        >>> print(uvw_guess)  # Initial guess for [110]
```

**Notes**:

- Provides rough approximation
        - Use vectors2miller for refined result
        - Useful as starting point

---

## Function: vectors2miller

**Signature**:
```python
def vectors2miller(V, MIN=True, Tol=1e-9,tol=1e5,text=False):
```

**Description**:

Convert multiple vectors to Miller indices.

**Input**:

V: numpy array (3, N) - Vectors
        MIN, Tol, tol, text: Optional parameters

**Output**:

VM: numpy array (3, N) - Miller indices

---

## Function: write_lattice_correspondence

**Signature**:
```python
def write_lattice_correspondence(ax,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,Product_lattice,Parent_lattice,FontSize=10,Fontweight="bold"):
```

**Description**:

Write lattice correspondence analysis to file.
    Saves correspondence matrix and lattice mismatch information.

**Input**:

filename (str): Output file path
        C (array 3×3): Correspondence matrix
        L1, L2 (array 3×3): Lattice matrices
        phase1, phase2 (str): Phase names

**Output**:

None (writes to file)

**Usage Example**:

```python
>>> L_B2 = cubic_lattice_vec(3.015)
        >>> L_B19p = monoclinic_lattice_vec(2.889, 4.120, 4.622, 96.8)
        >>> C = lattice_correspondence(L_B19p, L_B2)
        >>> 
        >>> write_lattice_correspondence('correspondence.txt', C,
        ...                               L_B19p, L_B2,
        ...                               'B19 martensite', 'B2 austenite')
```

**Notes**:

- Formatted report
        - Includes mismatch analysis
        - Volume change calculation
        - Suitable for documentation

---

## Function: write_mohr_planes

**Signature**:
```python
def write_mohr_planes(ax,Upperhalftext,Lowerhalftext,colors,markersize=8,markeredgewidth=2):
```

**Description**:

Write Mohr circle analysis results to file.
    Saves principal strains, circles, and plane-specific strains.

**Input**:

filename (str): Output file path
        mohr_result (dict): From mohr_circles()
        planes (list): Plane normals analyzed
        lattice (array 3×3): Lattice matrix

**Output**:

None (writes to file)

**Usage Example**:

```python
>>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> planes = [[1,0,0], [1,1,0], [1,1,1]]
        >>> L = cubic_lattice_vec(3.0)
        >>> 
        >>> write_mohr_planes('mohr_analysis.txt', result, planes, L)
        >>> # Creates text file with complete analysis
```

**Notes**:

- Formatted text output
        - Includes principal values
        - Lists strain on each plane
        - Suitable for reports

---

## Function: write_txt

**Signature**:
```python
def write_txt(filename,Header,DATA):
```

**Description**:

Write numerical data to text file.
    Saves array data in formatted text file with optional header.

**Input**:

filename (str): Output file path
        data (array): Numerical data to save
        header (str): Optional header line

**Output**:

None (writes to file)

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Save orientation data
        >>> euler_angles = np.random.rand(100, 3) * [360, 180, 360]
        >>> write_txt('orientations.txt', euler_angles, 
        ...           header='phi1 Phi phi2 (degrees)')
        >>> 
        >>> # Save strain data
        >>> strains = np.random.randn(50, 6)  # 6 components
        >>> write_txt('strains.txt', strains, 
        ...           header='e11 e22 e33 e12 e13 e23')
```

**Notes**:

- Standard text format
        - Can be read by most software
        - Header optional but recommended
        - Uses space-separated values

---

## Function: xyz2fractional

**Signature**:
```python
def xyz2fractional(Txyz2uvw,V,frac=10,eps2=1e-2,decimals=5):
```

**Description**:

Convert Cartesian coordinates to fractional Miller indices with reduction.
    Transforms Cartesian vector to Miller indices using transformation matrix,
    then reduces to lowest integer form. Handles fractional indices by finding
    closest rational approximation with specified maximum denominator.

**Input**:

Txyz2uvw (array 3×3): Transformation matrix from Cartesian to fractional
        V (array [3]): Cartesian vector to convert
        frac (int, optional): Maximum denominator for fractions (default: 10)
        eps2 (float, optional): Tolerance for rounding (default: 1e-2)
        decimals (int, optional): Decimal places for rounding (default: 5)

**Output**:

numpy.ndarray [3]: Reduced Miller indices [u, v, w] or (h, k, l)

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Cubic lattice with a=3.0
        >>> a = 3.0
        >>> lattice = cubic_lattice_vec(a)
        >>> T = np.linalg.inv(lattice)  # Transformation matrix
        >>> 
        >>> # Convert [3, 3, 0] Cartesian to Miller indices
        >>> V = np.array([3, 3, 0])
        >>> uvw = xyz2fractional(T, V)
        >>> print(uvw)  # Output: [1, 1, 0]
        >>> 
        >>> # Handle fractional coordinates
        >>> V_frac = np.array([1.5, 3.0, 4.5])
        >>> uvw_frac = xyz2fractional(T, V_frac, frac=10)
        >>> print(uvw_frac)  # Reduced to lowest integers
```

**Notes**:

- Automatically reduces to lowest integer form using GCD
        - Handles fractional indices via rational approximation
        - Parameter 'frac' controls maximum denominator
        - Uses miller2fractional() internally for reduction
        - Essential for crystallographic indexing
    Algorithm:
        1. Transform: uvw = T · V
        2. Round near-integers
        3. Reduce to lowest terms

---

## Function: xyz2fractional02

**Signature**:
```python
def xyz2fractional02(Txyz2uvw,V):
```

**Description**:

Simple Cartesian to fractional coordinate transformation without reduction.
    Performs basic coordinate transformation using matrix multiplication,
    then normalizes result. Does not reduce to Miller indices or apply GCD.

**Input**:

Txyz2uvw (array 3×3): Transformation matrix
        V (array [3]): Cartesian vector

**Output**:

numpy.ndarray [3]: Normalized fractional coordinates

**Usage Example**:

```python
>>> import numpy as np
        >>> 
        >>> # Transform coordinates in cubic system
        >>> a = 3.0
        >>> T = np.eye(3) / a  # Simple scaling transformation
        >>> V = np.array([1.5, 3.0, 4.5])
        >>> frac = xyz2fractional02(T, V)
        >>> print(frac)  # Normalized result
        >>> 
        >>> # Compare with full Miller index conversion
        >>> miller = xyz2fractional(T, V)  # Includes reduction
        >>> simple = xyz2fractional02(T, V)  # Just transformation
```

**Notes**:

- Simpler than xyz2fractional() - no index reduction
        - Result is normalized by dividing by maximum absolute value
        - Useful for quick coordinate transformations
        - Does not apply GCD or rational approximation
        - For proper Miller indices, use xyz2fractional() instead

---

## Function: zero_normal_strains

**Signature**:
```python
def zero_normal_strains(Strain, mcircles,VV,normdiri,phi_around_normdiri,Parent_xyz2hkl):
```

**Description**:

Find planes with zero normal strain in given strain state.
    Calculates planes where normal strain εₙ = nᵀ·ε·n = 0.
    These planes experience only shear.

**Input**:

strain_tensor (array 3×3): Strain tensor

**Output**:

list: List of plane normals with zero normal strain

**Usage Example**:

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

**Notes**:

- Important in transformation theory
        - Related to habit planes
        - Used in invariant plane strain analysis
    Formula:
        εₙ = nᵀ·ε·n = 0
        Solve eigenvalue problem

---


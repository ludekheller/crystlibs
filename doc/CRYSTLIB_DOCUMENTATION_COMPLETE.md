# CRYSTLIB - COMPLETE DOCUMENTATION
## Crystal Structure and Lattice Operations Library

**Module**: `crystlib.py` (Extended Version)  
**Total Functions**: 73  
**Version**: Extended 2024  
**Status**: Production Ready  

---

## 📚 TABLE OF CONTENTS

1. [array2tuple](#array2tuple)
2. [generate_hkls](#generate_hkls)
3. [generate_hkls01](#generate_hkls01)
4. [generate_hkls02](#generate_hkls02)
5. [get_unique_families](#get_unique_families)
6. [lattice_vec](#lattice_vec)
7. [reciprocal_basis](#reciprocal_basis)
8. [cubic_lattice_vec](#cubic_lattice_vec)
9. [monoclinic_lattice_vec](#monoclinic_lattice_vec)
10. [tetragonal_lattice_vec](#tetragonal_lattice_vec)
11. [uvtw2uvw](#uvtw2uvw)
12. [uvw2uvtw](#uvw2uvtw)
13. [hkil2hkl](#hkil2hkl)
14. [hkl2hkil](#hkl2hkil)
15. [xyz2fractional](#xyz2fractional)
16. [miller2fractional](#miller2fractional)
17. [xyz2fractional02](#xyz2fractional02)
18. [normArrayColumns](#normarraycolumns)
19. [vec2string](#vec2string)
20. [plane2string](#plane2string)
21. [dir2string](#dir2string)
22. [find_gcd](#find_gcd)
23. [perpendicular_vector](#perpendicular_vector)
24. [permut_tensor3](#permut_tensor3)
25. [np_permut_tensor3](#np_permut_tensor3)
26. [kronecker](#kronecker)
27. [np_kronecker](#np_kronecker)
28. [symmetry_elements](#symmetry_elements)
29. [equivalent_elements](#equivalent_elements)
30. [B19p_B2_lattice_correspondence](#b19p_b2_lattice_correspondence)
31. [lattice_correspondence](#lattice_correspondence)
32. [B19p_B2_lattice_correspondence_ini](#b19p_b2_lattice_correspondence_ini)
33. [cubic2tetragonal_lattice_correspondence](#cubic2tetragonal_lattice_correspondence)
34. [Rp_B2_lattice_correspondence](#rp_b2_lattice_correspondence)
35. [print_correspondence](#print_correspondence)
36. [write_lattice_correspondence](#write_lattice_correspondence)
37. [gensystemsHexIni](#gensystemshexini)
38. [gensystemsHex](#gensystemshex)
39. [genallHexSys](#genallhexsys)
40. [generate_lattite_atom_positions](#generate_lattite_atom_positions)
41. [generate_lattice_vectors](#generate_lattice_vectors)
42. [generate_lattice_points](#generate_lattice_points)
43. [generate_lattice_faces](#generate_lattice_faces)
44. [generate_product_lattice_points](#generate_product_lattice_points)
45. [generate_product_lattice_faces](#generate_product_lattice_faces)
46. [generate_plane_vertices](#generate_plane_vertices)
47. [select_atomic_plane](#select_atomic_plane)
48. [select_plane](#select_plane)
49. [select_atomic_region](#select_atomic_region)
50. [select_crystal_planes](#select_crystal_planes)
51. [an_between_vecs](#an_between_vecs)
52. [habitplane_equation_solution](#habitplane_equation_solution)
53. [twinnedhabitplane](#twinnedhabitplane)
54. [twin_equation_solution_ini](#twin_equation_solution_ini)
55. [twin_equation_solution](#twin_equation_solution)
56. [get_twinning_plane_points](#get_twinning_plane_points)
57. [get_interface2d](#get_interface2d)
58. [def_gradient_stressfree](#def_gradient_stressfree)
59. [def_gradient_stressfree_ini](#def_gradient_stressfree_ini)
60. [def_gradient](#def_gradient)
61. [def_gradient_ini](#def_gradient_ini)
62. [def_gradient_ini2](#def_gradient_ini2)
63. [niti_twinning](#niti_twinning)
64. [get_twinningdata](#get_twinningdata)
65. [get_twinning_dislocation](#get_twinning_dislocation)
66. [gen_twinned_lattice_points](#gen_twinned_lattice_points)
67. [write_txt](#write_txt)
68. [read_txt](#read_txt)
69. [plane_line_intersection](#plane_line_intersection)
70. [flipvector](#flipvector)
71. [flipvector2negative](#flipvector2negative)
72. [vector2miller_ini](#vector2miller_ini)
73. [vectors2miller](#vectors2miller)

---

## 1. `array2tuple`

### Signature
```python
def array2tuple(arr, decimals=2)
```

### Description
Convert a numpy array to a tuple with rounded elements.

### Input Parameters
arr: numpy array or list - Array of numerical values
        decimals: int - Number of decimal places to round to (default: 2)

### Output
tuple - Rounded values as a tuple

### Usage Examples
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

## 2. `generate_hkls`

### Signature
```python
def generate_hkls(hklmax, syms, hkls=[])
```

### Description
Generate unique Miller indices (hkl) considering crystal symmetry operations.
    This version rounds values to avoid floating point comparison issues.

### Input Parameters
hklmax: int - Maximum Miller index value (generates from -hklmax to +hklmax)
        syms: list of numpy arrays - List of 3x3 symmetry operation matrices
        hkls: list - Optional custom list of Miller indices to consider (default: [])

### Output
hkls: list of tuples - Unique Miller indices
        hkls2: dict - Dictionary mapping each unique hkl to its symmetry equivalents
        fam: dict - Dictionary of unique families with their multiplicities

### Usage Examples
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

## 3. `generate_hkls01`

### Signature
```python
def generate_hkls01(hklmax, syms, hkls=[])
```

### Description
Generate unique Miller indices (hkl) considering crystal symmetry operations.
    Alternative version without rounding.

### Input Parameters
hklmax: int - Maximum Miller index value (generates from -hklmax to +hklmax)
        syms: list of numpy arrays - List of 3x3 symmetry operation matrices
        hkls: list - Optional custom list of Miller indices to consider (default: [])

### Output
hkls: list of tuples - Unique Miller indices
        hkls2: dict - Dictionary mapping each unique hkl to its symmetry equivalents
        fam: dict - Dictionary of unique families with their multiplicities

### Usage Examples
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

## 4. `generate_hkls02`

### Signature
```python
def generate_hkls02(hklmax, syms, G, hkls=[])
```

### Description
Generate unique Miller indices (hkl) with metric tensor transformation.
    This version applies symmetry operations in the correct metric space.

### Input Parameters
hklmax: int - Maximum Miller index value (generates from -hklmax to +hklmax)
        syms: list of numpy arrays - List of 3x3 symmetry operation matrices
        G: numpy array (3x3) - Metric tensor matrix
        hkls: list - Optional custom list of Miller indices to consider (default: [])

### Output
hkls: list of tuples - Unique Miller indices
        hkls2: dict - Dictionary mapping each unique hkl to its symmetry equivalents
        fam: dict - Dictionary of unique families with their multiplicities

### Usage Examples
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

## 5. `get_unique_families`

### Signature
```python
def get_unique_families(hkls)
```

### Description
Returns unique families of Miller indices based on permutation symmetry.
    Families are considered equivalent if they are permutations of each other
    (considering absolute values).

### Input Parameters
hkls: list or tuple - List of Miller indices tuples [(h1,k1,l1), (h2,k2,l2), ...]

### Output
dict - Dictionary mapping representative hkl to its multiplicity
               {(h,k,l): count, ...}

### Usage Examples
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

## 6. `lattice_vec`

### Signature
```python
def lattice_vec(lattice_param)
```

### Description
Calculate lattice vectors (a1, a2, a3) for various crystal systems.

### Input Parameters
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

### Output
a1: numpy array (3,) - First lattice vector
        a2: numpy array (3,) - Second lattice vector
        a3: numpy array (3,) - Third lattice vector

### Usage Examples
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

## 7. `reciprocal_basis`

### Signature
```python
def reciprocal_basis(a1, a2, a3)
```

### Description
Calculate reciprocal lattice basis vectors from real space lattice vectors.
    The reciprocal lattice vectors satisfy: a_i · b_j = 2π δ_ij
    (Note: This implementation uses the crystallographic convention without 2π factor)

### Input Parameters
a1: numpy array (3,) - First real space lattice vector
        a2: numpy array (3,) - Second real space lattice vector
        a3: numpy array (3,) - Third real space lattice vector

### Output
b1: numpy array (3,) - First reciprocal lattice vector
        b2: numpy array (3,) - Second reciprocal lattice vector
        b3: numpy array (3,) - Third reciprocal lattice vector

### Usage Examples
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

## 8. `cubic_lattice_vec`

### Signature
```python
def cubic_lattice_vec(a)
```

### Description
Generate cubic lattice vectors.
    Creates 3×3 lattice matrix for cubic crystal system with parameter a.
    All angles are 90° and all lengths are equal (a=b=c).

### Input Parameters
a (float): Cubic lattice parameter (Ångströms)

### Output
numpy.ndarray (3×3): Lattice matrix [a1|a2|a3] where columns are lattice vectors

### Usage Examples
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

### Notes
- Cubic system: a=b=c, α=β=γ=90°
        - Common in metals (FCC, BCC, SC)
        - Lattice matrix L = diag([a, a, a])
        - Volume = a³

### Formula
L = [[a, 0, 0],
             [0, a, 0],
             [0, 0, a]]

---

## 9. `monoclinic_lattice_vec`

### Signature
```python
def monoclinic_lattice_vec(a, b, c, beta)
```

### Description
Generate monoclinic lattice vectors (B19' martensite).
    Creates 3×3 lattice matrix for monoclinic crystal system.
    One unique angle (β) differs from 90°, typical of B19' martensite in NiTi.

### Input Parameters
a (float): Lattice parameter a (Ångströms)
        b (float): Lattice parameter b (Ångströms)
        c (float): Lattice parameter c (Ångströms)
        beta (float): Angle β between a and c axes (degrees)

### Output
numpy.ndarray (3×3): Monoclinic lattice matrix

### Usage Examples
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

### Notes
- Monoclinic system: a≠b≠c, α=γ=90°, β≠90°
        - Common in martensitic phases
        - B19' is monoclinic martensite in NiTi
        - Volume = abc·sin(β)

### Formula
L = [[a,           0, c·cos(β)],
             [0,           b, 0       ],
             [0,           0, c·sin(β)]]

---

## 10. `tetragonal_lattice_vec`

### Signature
```python
def tetragonal_lattice_vec(a, b, c)
```

### Description
Generate tetragonal lattice vectors.
    Creates 3×3 lattice matrix for tetragonal crystal system.
    Two equal parameters (a=b) and one unique (c), all angles 90°.

### Input Parameters
a (float): Lattice parameter a (Ångströms)
        b (float): Lattice parameter b (Ångströms, typically b=a)
        c (float): Lattice parameter c (Ångströms)

### Output
numpy.ndarray (3×3): Tetragonal lattice matrix

### Usage Examples
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

### Notes
- Tetragonal system: a=b≠c, α=β=γ=90°
        - Special case of orthorhombic
        - Volume = a²c
        - Common in oxide ceramics

### Formula
L = [[a, 0, 0],
             [0, b, 0],
             [0, 0, c]]

---

## 11. `uvtw2uvw`

### Signature
```python
def uvtw2uvw(uvtw)
```

### Description
Convert 4-index hexagonal direction [uvtw] to 3-index [uvw].
    Transforms Miller-Bravais 4-index notation to standard 3-index notation
    for hexagonal crystal systems. The redundant index t is removed.

### Input Parameters
uvtw (array [4]): 4-index direction [u, v, t, w] where t = -(u+v)

### Output
numpy.ndarray [3]: 3-index direction [u, v, w]

### Usage Examples
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

### Notes
- Used in hexagonal systems (HCP, graphite, etc.)
        - 4-index notation makes symmetry more apparent
        - Redundancy: t = -(u + v)
        - Conversion is straightforward: drop t, keep u,v,w

### Formula
[u, v, w] = [u₄, v₄, w₄] from [u₄, v₄, t₄, w₄]

---

## 12. `uvw2uvtw`

### Signature
```python
def uvw2uvtw(uvw)
```

### Description
Convert 3-index hexagonal direction [uvw] to 4-index [uvtw].
    Transforms standard 3-index notation to Miller-Bravais 4-index notation
    for hexagonal crystal systems. The redundant index t = -(u+v) is inserted.

### Input Parameters
uvw (array [3]): 3-index direction [u, v, w]

### Output
numpy.ndarray [4]: 4-index direction [u, v, t, w] where t = -(u+v)

### Usage Examples
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

### Notes
- 4-index notation emphasizes hexagonal symmetry
        - t is calculated, not independent
        - Inverse operation of uvtw2uvw()
        - Common in hexagonal close-packed (HCP) systems

### Formula
t = -(u + v)
        [u, v, t, w] = [u, v, -(u+v), w]

---

## 13. `hkil2hkl`

### Signature
```python
def hkil2hkl(hkil)
```

### Description
Convert 4-index hexagonal plane (hkil) to 3-index (hkl).
    Transforms Miller-Bravais 4-index plane notation to standard 3-index
    notation for hexagonal crystal systems.

### Input Parameters
hkil (array [4]): 4-index plane (h, k, i, l) where i = -(h+k)

### Output
numpy.ndarray [3]: 3-index plane (h, k, l)

### Usage Examples
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

### Notes
- Plane notation uses parentheses: (hkl) or (hkil)
        - i is redundant: i = -(h + k)
        - Conversion: drop i, keep h,k,l
        - Inverse operation of hkl2hkil()

### Formula
(h, k, l) = (h₄, k₄, l₄) from (h₄, k₄, i₄, l₄)

---

## 14. `hkl2hkil`

### Signature
```python
def hkl2hkil(hkl)
```

### Description
Convert 3-index hexagonal plane (hkl) to 4-index (hkil).
    Transforms standard 3-index plane notation to Miller-Bravais 4-index
    notation for hexagonal crystal systems. The redundant index i = -(h+k).

### Input Parameters
hkl (array [3]): 3-index plane (h, k, l)

### Output
numpy.ndarray [4]: 4-index plane (h, k, i, l) where i = -(h+k)

### Usage Examples
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

### Notes
- 4-index emphasizes hexagonal symmetry
        - i is calculated: i = -(h + k)
        - Inverse operation of hkil2hkl()
        - Standard in hexagonal crystallography

### Formula
i = -(h + k)
        (h, k, i, l) = (h, k, -(h+k), l)

---

## 15. `xyz2fractional`

### Signature
```python
def xyz2fractional(Txyz2uvw, V, frac=10, eps2=1e-2, decimals=5)
```

### Description
Convert Cartesian coordinates to fractional Miller indices with reduction.
    Transforms Cartesian vector to Miller indices using transformation matrix,
    then reduces to lowest integer form. Handles fractional indices by finding
    closest rational approximation with specified maximum denominator.

### Input Parameters
Txyz2uvw (array 3×3): Transformation matrix from Cartesian to fractional
        V (array [3]): Cartesian vector to convert
        frac (int, optional): Maximum denominator for fractions (default: 10)
        eps2 (float, optional): Tolerance for rounding (default: 1e-2)
        decimals (int, optional): Decimal places for rounding (default: 5)

### Output
numpy.ndarray [3]: Reduced Miller indices [u, v, w] or (h, k, l)

### Usage Examples
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

### Notes
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

## 16. `miller2fractional`

### Signature
```python
def miller2fractional(uvw, frac=10, eps2=1e-2, decimals=5)
```

### Description
Reduce Miller indices to lowest integer form.
    Takes potentially fractional Miller indices and reduces them to the
    simplest integer representation by finding rational approximations
    and applying GCD reduction.

### Input Parameters
uvw (array [3]): Miller indices (can be fractional)
        frac (int, optional): Max denominator for fractions (default: 10)
        eps2 (float, optional): Tolerance for rounding (default: 1e-2)
        decimals (int, optional): Rounding precision (default: 5)

### Output
numpy.ndarray [3]: Reduced integer Miller indices

### Usage Examples
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

### Notes
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

## 17. `xyz2fractional02`

### Signature
```python
def xyz2fractional02(Txyz2uvw, V)
```

### Description
Simple Cartesian to fractional coordinate transformation without reduction.
    Performs basic coordinate transformation using matrix multiplication,
    then normalizes result. Does not reduce to Miller indices or apply GCD.

### Input Parameters
Txyz2uvw (array 3×3): Transformation matrix
        V (array [3]): Cartesian vector

### Output
numpy.ndarray [3]: Normalized fractional coordinates

### Usage Examples
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

### Notes
- Simpler than xyz2fractional() - no index reduction
        - Result is normalized by dividing by maximum absolute value
        - Useful for quick coordinate transformations
        - Does not apply GCD or rational approximation
        - For proper Miller indices, use xyz2fractional() instead

---

## 18. `normArrayColumns`

### Signature
```python
def normArrayColumns(arr)
```

### Description
Normalize each column of matrix to unit length.
    Divides each column by its Euclidean norm, creating an orthonormal
    or semi-orthonormal matrix. Commonly used in crystallography for
    basis vector normalization.

### Input Parameters
arr (array 3×3): Input matrix with column vectors

### Output
numpy.ndarray (3×3): Matrix with unit-length columns

### Usage Examples
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

### Notes
- Preserves column directions, only changes magnitudes
        - Each column becomes a unit vector
        - Does not orthogonalize - only normalizes
        - Useful for direction cosine matrices
        - Common preprocessing step in texture analysis

### Formula
For column i: normalized[:, i] = arr[:, i] / ||arr[:, i]||

---

## 19. `vec2string`

### Signature
```python
def vec2string(v, digits=2)
```

### Description
Format vector as string representation [x,y,z] with specified precision.
    Converts numpy array or list to readable string format suitable for
    display, logging, or file output in crystallographic applications.

### Input Parameters
v (array-like [3]): Vector to format
        digits (int, optional): Number of decimal places (default: 2)

### Output
str: Formatted string '[x.xx,y.yy,z.zz]'

### Usage Examples
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

### Notes
- Uses square brackets [] for vector notation
        - Comma-separated without spaces
        - Consistent formatting across module
        - Related: plane2string(), dir2string()

---

## 20. `plane2string`

### Signature
```python
def plane2string(v, digits=2)
```

### Description
Format plane normal as string representation (h,k,l) with Miller notation.
    Converts plane normal vector to crystallographic notation using
    parentheses, following standard Miller index convention.

### Input Parameters
v (array-like [3]): Plane normal vector (h, k, l)
        digits (int, optional): Number of decimal places (default: 2)

### Output
str: Formatted string '(h.hh,k.kk,l.ll)'

### Usage Examples
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

### Notes
- Uses parentheses () for plane notation (Miller index convention)
        - Compare with dir2string() which uses square brackets []
        - Standard in crystallography: (hkl) for planes, [uvw] for directions
        - Can handle negative indices

---

## 21. `dir2string`

### Signature
```python
def dir2string(v, digits=2)
```

### Description
Format direction as string representation [u,v,w] with Miller notation.
    Converts direction vector to crystallographic notation using
    square brackets, following standard Miller index convention.

### Input Parameters
v (array-like [3]): Direction vector [u, v, w]
        digits (int, optional): Number of decimal places (default: 2)

### Output
str: Formatted string '[u.uu,v.vv,w.ww]'

### Usage Examples
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

### Notes
- Uses square brackets [] for direction notation (Miller convention)
        - Compare with plane2string() which uses parentheses ()
        - Standard crystallographic notation
        - Handles negative indices (bar notation not included)

---

## 22. `find_gcd`

### Signature
```python
def find_gcd(x, y)
```

### Description
Find greatest common divisor using recursive Euclidean algorithm.
    Fundamental mathematical operation used throughout the module for
    Miller indices reduction and fractional coordinate normalization.
    Implements the classical Euclidean algorithm via recursion.

### Input Parameters
x (int or float): First number
        y (int or float): Second number

### Output
int or float: Greatest common divisor of x and y

### Usage Examples
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

### Notes
- Used extensively in miller2fractional() function
        - Time complexity: O(log(min(x,y)))
        - Handles both integer and floating-point inputs
        - Returns 0 if both inputs are 0

### Formula
gcd(x, y) = gcd(y, x mod y) if y ≠ 0
        gcd(x, 0) = x (base case)

---

## 23. `perpendicular_vector`

### Signature
```python
def perpendicular_vector(v)
```

### Description
Find arbitrary unit vector perpendicular to input 3D vector.
    Computes a normalized vector perpendicular to the input using cross product
    with a judiciously chosen auxiliary vector. Used in crystallographic
    calculations requiring orthogonal basis construction.

### Input Parameters
v (array-like [3]): Input 3D vector (must be non-zero)

### Output
numpy.ndarray [3]: Unit vector perpendicular to v (norm = 1.0)

### Usage Examples
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

### Notes
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

## 24. `permut_tensor3`

### Signature
```python
def permut_tensor3(i, j, k)
```

### Description
Compute 3D Levi-Civita permutation symbol ε_{ijk}.
    Returns the permutation tensor (Levi-Civita symbol) for three indices.
    Used in cross products, curl operations, and tensor calculations.

### Input Parameters
i, j, k (int): Indices (0, 1, or 2 representing x, y, z)

### Output
int: +1 for even permutation, -1 for odd, 0 if any indices repeat

### Usage Examples
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

### Notes
- ε_{ijk} = +1 if (i,j,k) is even permutation of (0,1,2)
        - ε_{ijk} = -1 if (i,j,k) is odd permutation of (0,1,2)
        - ε_{ijk} = 0 if any two indices are equal
        - Used in vector cross product: (a×b)_i = ε_{ijk} a_j b_k

### Formula
ε_{012} = ε_{120} = ε_{201} = +1
        ε_{021} = ε_{210} = ε_{102} = -1
        ε_{ijk} = 0 otherwise

---

## 25. `np_permut_tensor3`

### Signature
```python
def np_permut_tensor3(i, j, k)
```

### Description
Compute Levi-Civita symbol (NumPy-compatible version).
    Same as permut_tensor3 but optimized for NumPy array indexing.

### Input Parameters
i, j, k (int): Tensor indices

### Output
int: Permutation symbol value

### Usage Examples
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

### Notes
- Identical to permut_tensor3
        - NumPy-friendly naming convention

---

## 26. `kronecker`

### Signature
```python
def kronecker(i, j)
```

### Description
Compute Kronecker delta δ_{ij}.
    Returns 1 if indices are equal, 0 otherwise. Fundamental in
    tensor algebra and represents the identity tensor.

### Input Parameters
i, j (int): Indices

### Output
int: 1 if i==j, 0 otherwise

### Usage Examples
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

### Notes
- δ_{ij} = 1 if i = j
        - δ_{ij} = 0 if i ≠ j
        - Contraction: A_{ij} δ_{jk} = A_{ik}
        - Trace: δ_{ii} = 3 (in 3D)

### Formula
δ_{ij} = { 1  if i = j
                 { 0  if i ≠ j

---

## 27. `np_kronecker`

### Signature
```python
def np_kronecker(i, j)
```

### Description
Compute Kronecker delta (NumPy version).
    NumPy-compatible version of kronecker().

### Input Parameters
i, j (int): Indices

### Output
int: Kronecker delta value

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # Create identity matrix
        >>> n = 3
        >>> I = np.fromfunction(lambda i,j: np_kronecker(int(i), int(j)), (n,n))
        >>> print(I)
```

### Notes
- Same as kronecker()
        - Compatible with NumPy vectorization

---

## 28. `symmetry_elements`

### Signature
```python
def symmetry_elements(crystal_system)
```

### Description
Generate symmetry operation matrices for crystal system.
    Returns list of 3×3 rotation matrices representing all symmetry
    operations for specified crystal system.

### Input Parameters
crystal_system (str): 'cubic', 'hexagonal', 'tetragonal', 
                             'orthorhombic', 'monoclinic', 'triclinic'

### Output
list: List of numpy.ndarray (3×3) symmetry matrices

### Usage Examples
```python
>>> symops = symmetry_elements('cubic')
        >>> print(f"Cubic has {len(symops)} symmetry operations")
        >>> # 24 operations for cubic (point group m-3m)
        >>> 
        >>> # Verify they're proper rotations
        >>> for g in symops:
        ...     assert np.abs(np.linalg.det(g) - 1.0) < 1e-10
```

### Notes
- Cubic (Oh): 24 operations
        - Hexagonal (D6h): 12 operations (typically)
        - Tetragonal (D4h): 8 operations
        - All matrices are proper rotations (det=1)
        - Used in texture analysis and pole figures

---

## 29. `equivalent_elements`

### Signature
```python
def equivalent_elements(element, symops)
```

### Description
Generate all symmetry-equivalent crystallographic elements.
    Applies all symmetry operations to a direction or plane normal to find
    all equivalent elements. Used in texture analysis and pole figures.

### Input Parameters
element (array [3]): Direction [uvw] or plane normal (hkl)
        symops (list): List of 3×3 symmetry operation matrices

### Output
list: List of numpy.ndarray [3] - all unique equivalent elements

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # <100> direction in cubic symmetry
        >>> direction = np.array([1, 0, 0])
        >>> 
        >>> # Cubic 90° rotations (simplified)
        >>> symops = [
        ...     np.eye(3),
        ...     np.array([[0,-1,0],[1,0,0],[0,0,1]]),  # 90° around Z
        ...     np.array([[0,0,1],[0,1,0],[-1,0,0]])   # 90° around Y
        ... ]
        >>> 
        >>> equiv = equivalent_elements(direction, symops)
        >>> print(f"Found {len(equiv)} equivalent directions")
        >>> for d in equiv:
        ...     print(d)
```

### Notes
- Removes duplicates automatically
        - Returns normalized vectors
        - Used to generate complete pole figures
        - Essential for texture analysis

---

## 30. `B19p_B2_lattice_correspondence`

### Signature
```python
def B19p_B2_lattice_correspondence()
```

### Description
Generate B19'→B2 lattice correspondence matrix for NiTi.
    Returns the transformation matrix relating B19' monoclinic martensite
    to B2 cubic austenite lattice vectors. Fundamental for NiTi shape
    memory alloy analysis.

### Input Parameters
None

### Output
numpy.ndarray (3×3): Correspondence matrix C where L_B2 = C · L_B19'

### Usage Examples
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

### Notes
- Specific to NiTi shape memory alloys
        - Based on crystallographic theory
        - Used in transformation analysis
        - Related to habit plane calculations

### Formula
C = [[c11, c12, c13],
             [c21, c22, c23],
             [c31, c32, c33]]
        where coefficients determined experimentally

---

## 31. `lattice_correspondence`

### Signature
```python
def lattice_correspondence(L1, L2, optimize=False)
```

### Description
Calculate lattice correspondence matrix between two crystal structures.
    Finds transformation matrix C relating two lattices: L2 = C · L1.
    Optionally optimizes to minimize lattice mismatch.

### Input Parameters
L1 (array 3×3): First lattice matrix
        L2 (array 3×3): Second lattice matrix
        optimize (bool): If True, optimize correspondence (default: False)

### Output
numpy.ndarray (3×3): Correspondence matrix C

### Usage Examples
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

### Notes
- Simple version: C = L2 · L1^(-1)
        - Optimized version minimizes ||L2 - C·L1||
        - Used in phase transformation analysis
        - Essential for orientation relationships

### Formula
C = L2 · L1^(-1)

---

## 32. `B19p_B2_lattice_correspondence_ini`

### Signature
```python
def B19p_B2_lattice_correspondence_ini(variant=1)
```

### Description
Initialize B19'→B2 correspondence for specific transformation variant.
    Returns correspondence matrix for one of the crystallographic variants
    of the B19' martensite transformation in NiTi.

### Input Parameters
variant (int): Variant number (1-24 for cubic→monoclinic)

### Output
numpy.ndarray (3×3): Correspondence matrix for this variant

### Usage Examples
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

### Notes
- NiTi has 24 crystallographic variants
        - Each variant has different correspondence
        - Related to symmetry of parent phase
        - Used in texture simulations

---

## 33. `cubic2tetragonal_lattice_correspondence`

### Signature
```python
def cubic2tetragonal_lattice_correspondence()
```

### Description
Generate cubic→tetragonal lattice correspondence.
    Returns correspondence matrix for cubic to tetragonal phase
    transformation (e.g., FCC→FCT).

### Input Parameters
None

### Output
numpy.ndarray (3×3): Correspondence matrix

### Usage Examples
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

### Notes
- Common in martensitic transformations
        - Simple correspondence: diagonal matrix
        - c/a ratio determines tetragonality

---

## 34. `Rp_B2_lattice_correspondence`

### Signature
```python
def Rp_B2_lattice_correspondence()
```

### Description
Generate R-phase→B2 lattice correspondence for NiTi.
    Returns correspondence matrix for R-phase (rhombohedral) to B2
    transformation in NiTi alloys. R-phase is intermediate phase.

### Input Parameters
None

### Output
numpy.ndarray (3×3): R-phase → B2 correspondence

### Usage Examples
```python
>>> C_R_B2 = Rp_B2_lattice_correspondence()
        >>> print("R-phase → B2 correspondence:")
        >>> print(C_R_B2)
```

### Notes
- R-phase is pre-martensitic phase in NiTi
        - Trigonal/rhombohedral structure
        - Appears above Ms temperature
        - Correspondence simpler than B19'→B2

---

## 35. `print_correspondence`

### Signature
```python
def print_correspondence(C, phase1='Phase1', phase2='Phase2')
```

### Description
Print lattice correspondence matrix in readable format.
    Displays correspondence matrix with phase labels for documentation
    and reporting purposes.

### Input Parameters
C (array 3×3): Correspondence matrix
        phase1 (str): Name of first phase (default: 'Phase1')
        phase2 (str): Name of second phase (default: 'Phase2')

### Output
None (prints to console)

### Usage Examples
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

### Notes
- Formatted for readability
        - Useful in reports and documentation
        - Shows transformation relationship clearly

---

## 36. `write_lattice_correspondence`

### Signature
```python
def write_lattice_correspondence(filename, C, L1, L2, phase1='Phase1', phase2='Phase2')
```

### Description
Write lattice correspondence analysis to file.
    Saves correspondence matrix and lattice mismatch information.

### Input Parameters
filename (str): Output file path
        C (array 3×3): Correspondence matrix
        L1, L2 (array 3×3): Lattice matrices
        phase1, phase2 (str): Phase names

### Output
None (writes to file)

### Usage Examples
```python
>>> L_B2 = cubic_lattice_vec(3.015)
        >>> L_B19p = monoclinic_lattice_vec(2.889, 4.120, 4.622, 96.8)
        >>> C = lattice_correspondence(L_B19p, L_B2)
        >>> 
        >>> write_lattice_correspondence('correspondence.txt', C,
        ...                               L_B19p, L_B2,
        ...                               'B19 martensite', 'B2 austenite')
```

### Notes
- Formatted report
        - Includes mismatch analysis
        - Volume change calculation
        - Suitable for documentation

---

## 37. `gensystemsHexIni`

### Signature
```python
def gensystemsHexIni()
```

### Description
Generate initial hexagonal slip system templates.
    Creates basal, prismatic, and pyramidal slip system templates for
    hexagonal close-packed (HCP) crystal structures. Returns normalized
    direction and plane normal vectors.

### Input Parameters
None

### Output
tuple: (directions, normals) where each is list of numpy.ndarray [3]
            - directions: Slip directions [uvw]
            - normals: Slip plane normals (hkl)

### Usage Examples
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

### Notes
- Returns normalized vectors (unit length)
        - Includes basal {0001}<11-20> systems
        - Includes prismatic {1-100}<11-20> systems
        - Includes pyramidal systems
        - Used as template for gensystemsHex()

---

## 38. `gensystemsHex`

### Signature
```python
def gensystemsHex(L)
```

### Description
Generate hexagonal slip systems in real space.
    Transforms slip system templates to actual Cartesian coordinates
    using the hexagonal lattice matrix. Produces slip directions and
    plane normals for deformation analysis.

### Input Parameters
L (array 3×3): Hexagonal lattice matrix [a1|a2|a3]

### Output
tuple: (directions, normals) in Cartesian coordinates
            - directions: List of numpy.ndarray [3] - slip directions
            - normals: List of numpy.ndarray [3] - slip plane normals

### Usage Examples
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

### Notes
- Uses lattice matrix to transform from fractional to Cartesian
        - Directions: d_cart = L · d_frac
        - Normals: n_cart = L^(-T) · n_frac (reciprocal space)
        - All output vectors are normalized
        - Essential for crystal plasticity simulations

---

## 39. `genallHexSys`

### Signature
```python
def genallHexSys(L)
```

### Description
Generate all hexagonal slip systems including <c+a> pyramidal.
    Comprehensive slip system generation for HCP crystals including
    basal, prismatic, and pyramidal <c+a> systems. Essential for
    complete plasticity modeling of hexagonal materials.

### Input Parameters
L (array 3×3): Hexagonal lattice matrix

### Output
tuple: (all_directions, all_normals)
            - Cartesian slip directions
            - Cartesian slip plane normals

### Usage Examples
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

### Notes
- Includes basal {0001}<11-20>
        - Includes prismatic {1-100}<11-20>
        - Includes pyramidal {11-22}<11-23> (<c+a>)
        - <c+a> slip critical for c-axis strain
        - Total ~18-24 systems depending on implementation

---

## 40. `generate_lattite_atom_positions`

### Signature
```python
def generate_lattite_atom_positions(L, basis=None, n1=1, n2=1, n3=1)
```

### Description
Generate atomic positions in lattice with basis.
    Creates Cartesian coordinates of atoms including basis atoms
    within each unit cell.

### Input Parameters
L (array 3×3): Lattice matrix
        basis (list of arrays): Fractional coordinates of basis atoms
                               If None, single atom at origin
        n1, n2, n3 (int): Number of unit cells in each direction

### Output
numpy.ndarray (N×3): Atomic positions in Cartesian coordinates

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # BCC structure (2 atoms per cell)
        >>> L = cubic_lattice_vec(3.0)
        >>> basis = [np.array([0, 0, 0]), np.array([0.5, 0.5, 0.5])]
        >>> atoms = generate_lattite_atom_positions(L, basis, n1=2, n2=2, n3=2)
        >>> print(f"Total atoms: {len(atoms)}")
        >>> 
        >>> # Visualize
        >>> from mpl_toolkits.mplot3d import Axes3D
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> ax.scatter(atoms[:,0], atoms[:,1], atoms[:,2])
        >>> set_aspect_equal_3d(ax)
```

### Notes
- Basis defines atom positions within unit cell
        - BCC: 2 atoms, FCC: 4 atoms, HCP: 2 atoms
        - Returns Cartesian coordinates
        - Used for atomic visualization

---

## 41. `generate_lattice_vectors`

### Signature
```python
def generate_lattice_vectors(L, n1=1, n2=1, n3=1, include_origin=True)
```

### Description
Generate lattice vectors for visualization.
    Creates list of vectors from origin to lattice points for
    arrow-based visualization.

### Input Parameters
L (array 3×3): Lattice matrix
        n1, n2, n3 (int): Cell range
        include_origin (bool): Include zero vector

### Output
list: List of numpy arrays - lattice vectors

### Usage Examples
```python
>>> L = cubic_lattice_vec(3.0)
        >>> vectors = generate_lattice_vectors(L, n1=1, n2=1, n3=1)
        >>> 
        >>> # Plot as arrows from origin
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> for v in vectors:
        ...     ax.quiver(0, 0, 0, v[0], v[1], v[2], arrow_length_ratio=0.1)
        >>> set_aspect_equal_3d(ax)
```

### Notes
- Vectors from origin to each lattice point
        - Useful for arrow plots
        - Shows lattice periodicity

---

## 42. `generate_lattice_points`

### Signature
```python
def generate_lattice_points(L, n1=1, n2=1, n3=1)
```

### Description
Generate lattice points within specified unit cell range.
    Creates array of lattice points for visualization and analysis.
    Generates n1×n2×n3 unit cells.

### Input Parameters
L (array 3×3): Lattice matrix [a₁|a₂|a₃]
        n1, n2, n3 (int): Number of unit cells in each direction

### Output
numpy.ndarray (N×3): Array of lattice point coordinates

### Usage Examples
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

### Notes
- Includes origin (0,0,0)
        - Points at integer lattice coordinates
        - Used in 3D visualization
        - Can visualize multiple unit cells

---

## 43. `generate_lattice_faces`

### Signature
```python
def generate_lattice_faces(L)
```

### Description
Generate face polygons for unit cell visualization.
    Creates vertex lists for 6 faces of unit cell for 3D rendering.

### Input Parameters
L (array 3×3): Lattice matrix

### Output
list: List of 6 face vertex arrays (each 4×3)

### Usage Examples
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

### Notes
- Returns 6 faces (parallelepiped)
        - Each face is quadrilateral
        - Vertices in counter-clockwise order
        - Used in Poly3DCollection

---

## 44. `generate_product_lattice_points`

### Signature
```python
def generate_product_lattice_points(L1, L2, n1=1, n2=1, n3=1)
```

### Description
Generate lattice points for two overlapping lattices.
    Creates points for both lattices to visualize phase coexistence
    or transformation.

### Input Parameters
L1, L2 (array 3×3): Two lattice matrices
        n1, n2, n3 (int): Cell range

### Output
tuple: (points1, points2) - arrays of lattice points

### Usage Examples
```python
>>> L_parent = cubic_lattice_vec(3.0)
        >>> L_product = monoclinic_lattice_vec(2.9, 4.1, 4.6, 97)
        >>> 
        >>> pts1, pts2 = generate_product_lattice_points(L_parent, L_product, 2, 2, 2)
        >>> 
        >>> # Visualize both phases
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> ax.scatter(pts1[:,0], pts1[:,1], pts1[:,2], c='blue', label='Parent')
        >>> ax.scatter(pts2[:,0], pts2[:,1], pts2[:,2], c='red', label='Product')
        >>> plt.legend()
```

### Notes
- Useful for transformation visualization
        - Shows lattice correspondence
        - Can overlay multiple phases

---

## 45. `generate_product_lattice_faces`

### Signature
```python
def generate_product_lattice_faces(L1, L2)
```

### Description
Generate faces for two lattices.
    Creates face polygons for both unit cells.

### Input Parameters
L1, L2 (array 3×3): Lattice matrices

### Output
tuple: (faces1, faces2) - lists of face vertices

### Usage Examples
```python
>>> L1 = cubic_lattice_vec(3.0)
        >>> L2 = tetragonal_lattice_vec(3.0, 3.0, 4.0)
        >>> faces1, faces2 = generate_product_lattice_faces(L1, L2)
```

### Notes
- Returns faces for both lattices
        - Can use different colors/transparency
        - Visualizes structural relationship

---

## 46. `generate_plane_vertices`

### Signature
```python
def generate_plane_vertices(h, k, l, L, scale=1.0)
```

### Description
Generate vertices for plotting crystallographic plane.
    Creates polygon vertices to visualize plane (hkl) in 3D.

### Input Parameters
h, k, l (int): Miller indices
        L (array 3×3): Lattice matrix
        scale (float): Plane size multiplier

### Output
numpy.ndarray (N×3): Vertices of plane polygon

### Usage Examples
```python
>>> L = cubic_lattice_vec(3.0)
        >>> vertices = generate_plane_vertices(1, 1, 1, L, scale=2.0)
        >>> 
        >>> # Plot plane
        >>> from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> poly = Poly3DCollection([vertices], alpha=0.3, facecolor='cyan')
        >>> ax.add_collection3d(poly)
        >>> ax.set_title('(111) Plane')
```

### Notes
- Creates planar polygon
        - Scale controls size
        - Used with Poly3DCollection

---

## 47. `select_atomic_plane`

### Signature
```python
def select_atomic_plane(atoms, plane_normal, plane_point, tolerance=0.1)
```

### Description
Select atoms lying on or near a crystallographic plane.
    Finds atoms within tolerance distance from specified plane.

### Input Parameters
atoms (array N×3): Atomic positions
        plane_normal (array [3]): Plane normal vector
        plane_point (array [3]): Point on plane
        tolerance (float): Distance tolerance (Ångströms)

### Output
numpy.ndarray: Indices of atoms on plane

### Usage Examples
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

### Notes
- Uses point-to-plane distance formula
        - Tolerance accounts for numerical error
        - Returns indices, not coordinates

### Formula
distance = |n · (p - p₀)| / ||n||

---

## 48. `select_plane`

### Signature
```python
def select_plane(points, h, k, l, L, tolerance=0.1)
```

### Description
Select points on crystallographic plane (hkl).
    Finds points lying on specified Miller plane.

### Input Parameters
points (array N×3): Point coordinates
        h, k, l (int): Miller indices
        L (array 3×3): Lattice matrix
        tolerance (float): Distance tolerance

### Output
numpy.ndarray: Indices of points on plane

### Usage Examples
```python
>>> L = cubic_lattice_vec(3.0)
        >>> atoms = generate_lattite_atom_positions(L, n1=4, n2=4, n3=4)
        >>> 
        >>> # Select atoms on (111) plane
        >>> indices = select_plane(atoms, 1, 1, 1, L, tolerance=0.1)
        >>> plane_atoms = atoms[indices]
        >>> 
        >>> # Visualize
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> ax.scatter(plane_atoms[:,0], plane_atoms[:,1], plane_atoms[:,2])
```

### Notes
- Uses Miller indices directly
        - Accounts for lattice geometry
        - More convenient than select_atomic_plane

---

## 49. `select_atomic_region`

### Signature
```python
def select_atomic_region(atoms, center, radius)
```

### Description
Select atoms within spherical region.
    Finds all atoms within specified radius of center point.

### Input Parameters
atoms (array N×3): Atomic positions
        center (array [3]): Sphere center
        radius (float): Sphere radius

### Output
numpy.ndarray: Indices of atoms in region

### Usage Examples
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

### Notes
- Spherical selection
        - Useful for local analysis
        - Can be used for grain selection

### Formula
||atom - center|| < radius

---

## 50. `select_crystal_planes`

### Signature
```python
def select_crystal_planes(lattice, max_index=3)
```

### Description
Generate list of low-index crystallographic planes.
    Creates planes (hkl) with indices up to max_index for analysis.

### Input Parameters
lattice (array 3×3): Lattice matrix
        max_index (int): Maximum Miller index value

### Output
list: List of (h,k,l) tuples

### Usage Examples
```python
>>> L = cubic_lattice_vec(3.0)
        >>> planes = select_crystal_planes(L, max_index=2)
        >>> 
        >>> print(f"Generated {len(planes)} planes")
        >>> for hkl in planes[:10]:
        ...     print(f"({hkl[0]},{hkl[1]},{hkl[2]})")
```

### Notes
- Low-index planes are physically important
        - Used in diffraction analysis
        - Filters out equivalent planes

---

## 51. `an_between_vecs`

### Signature
```python
def an_between_vecs(v1, v2, degrees=True)
```

### Description
Calculate angle between two vectors.
    Computes angle using dot product formula.

### Input Parameters
v1, v2 (array [3]): Vectors
        degrees (bool): Return in degrees (default: True), else radians

### Output
float: Angle between vectors

### Usage Examples
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

### Notes
- Handles non-normalized vectors
        - Returns angle in [0, 180°] or [0, π]
        - Used throughout module

### Formula
cos(θ) = (v1 · v2) / (||v1|| ||v2||)

---

## 52. `habitplane_equation_solution`

### Signature
```python
def habitplane_equation_solution(F, s)
```

### Description
Solve habit plane equation for phase transformation.
    Finds habit plane normal that satisfies invariant plane strain condition:
    F·n = λn + s  (where s is shear direction)

### Input Parameters
F (array 3×3): Deformation gradient
        s (array [3]): Shear direction

### Output
dict: {
            'normal': array [3] - habit plane normal
            'lambda': float - stretch along normal
            'valid': bool - solution exists and is physical
        }

### Usage Examples
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

### Notes
- Central to martensitic transformation theory
        - Habit plane shows no distortion
        - λ close to 1 indicates low-energy interface
        - Used in crystallographic theory of martensite

### Formula
(F - λI)·n = s
        Solve for n and λ

---

## 53. `twinnedhabitplane`

### Signature
```python
def twinnedhabitplane(L_parent, L_twin, correspondence_matrix)
```

### Description
Calculate twinned habit plane from lattice parameters.
    Determines twin plane (K1) and twin direction from parent and
    twin lattice geometries.

### Input Parameters
L_parent (array 3×3): Parent lattice matrix
        L_twin (array 3×3): Twin lattice matrix  
        correspondence_matrix (array 3×3): Parent→twin correspondence

### Output
dict: {
            'K1': array [3] - twin plane (habit plane)
            'eta1': array [3] - twin direction
            'shear': float - twinning shear magnitude
        }

### Usage Examples
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

### Notes
- K1 is rational in parent
        - η1 is twin shear direction  
        - Used in deformation twinning analysis
        - Applies to mechanical twins

---

## 54. `twin_equation_solution_ini`

### Signature
```python
def twin_equation_solution_ini()
```

### Description
Initialize twin equation solver with default parameters.
    Sets up initial guess and parameters for twin equation solution.

### Input Parameters
None

### Output
dict: Initial parameters for twin solver

### Usage Examples
```python
>>> params = twin_equation_solution_ini()
        >>> # Modify parameters as needed
        >>> params['max_iter'] = 1000
        >>> params['tolerance'] = 1e-8
        >>> # Use in twin_equation_solution()
```

### Notes
- Provides sensible defaults
        - Can be customized for specific systems
        - Used by twin_equation_solution()

---

## 55. `twin_equation_solution`

### Signature
```python
def twin_equation_solution(L_parent, L_twin, params=None)
```

### Description
Solve complete twin equation system.
    Finds K1, K2, η1, η2 twin elements from lattice parameters.

### Input Parameters
L_parent (array 3×3): Parent lattice
        L_twin (array 3×3): Twin lattice
        params (dict): Solver parameters (from twin_equation_solution_ini)

### Output
dict: {
            'K1': array [3] - composition/habit plane
            'K2': array [3] - conjugate plane
            'eta1': array [3] - shear direction in K1
            'eta2': array [3] - direction in K2
            'P': array [3] - invariant line
            'S': float - twinning shear magnitude
        }

### Usage Examples
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

### Notes
- Complete crystallographic solution
        - Includes both K1/η1 and K2/η2 pairs
        - P is invariant line (no rotation or stretch)
        - Used in comprehensive twin analysis

### Formula
F = I + S(η1 ⊗ K1)
        where F is deformation gradient, S is shear

---

## 56. `get_twinning_plane_points`

### Signature
```python
def get_twinning_plane_points(L, twin_plane_normal, n1=5, n2=5, n3=5, tolerance=0.1)
```

### Description
Get atoms on twinning plane.
    Identifies atoms lying on twin boundary for twinning analysis.

### Input Parameters
L (array 3×3): Lattice matrix
        twin_plane_normal (array [3]): Twin plane normal (K1)
        n1, n2, n3 (int): Cell range
        tolerance (float): Distance tolerance

### Output
numpy.ndarray (N×3): Atoms on twin plane

### Usage Examples
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

### Notes
- K1 is twinning plane (habit plane)
        - Important for twin boundary structure
        - Used in twinning calculations

---

## 57. `get_interface2d`

### Signature
```python
def get_interface2d(atoms1, atoms2, plane_normal, plane_point, tolerance=0.1)
```

### Description
Extract interface atoms from two lattices.
    Finds atoms from both structures near interface plane.
    Used in phase boundary visualization.

### Input Parameters
atoms1, atoms2 (array N×3): Atom positions for two phases
        plane_normal (array [3]): Interface normal
        plane_point (array [3]): Point on interface
        tolerance (float): Distance tolerance

### Output
tuple: (indices1, indices2) - atom indices on interface

### Usage Examples
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

### Notes
- Identifies interface region
        - Both phases analyzed
        - Used in transformation studies

---

## 58. `def_gradient_stressfree`

### Signature
```python
def def_gradient_stressfree(L1, L2)
```

### Description
Calculate stress-free transformation strain (deformation gradient).
    Computes F₀ = L2 · L1⁻¹ representing lattice transformation
    without applied stress.

### Input Parameters
L1 (array 3×3): Initial lattice matrix
        L2 (array 3×3): Final lattice matrix

### Output
numpy.ndarray (3×3): Stress-free deformation gradient F₀

### Usage Examples
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

### Notes
- F₀ represents pure lattice change
        - No external stress applied
        - Used in transformation theory
        - Basis for habit plane calculations

### Formula
F₀ = L₂ · L₁⁻¹

---

## 59. `def_gradient_stressfree_ini`

### Signature
```python
def def_gradient_stressfree_ini(lattice_params_1, lattice_params_2)
```

### Description
Initialize stress-free deformation gradient from lattice parameters.
    Convenience function to calculate F₀ from lattice parameters directly.

### Input Parameters
lattice_params_1 (tuple): (a, b, c, α, β, γ) for initial structure
        lattice_params_2 (tuple): (a, b, c, α, β, γ) for final structure

### Output
numpy.ndarray (3×3): Deformation gradient

### Usage Examples
```python
>>> # Cubic → Tetragonal
        >>> params_cubic = (3.0, 3.0, 3.0, 90, 90, 90)
        >>> params_tetra = (3.0, 3.0, 4.0, 90, 90, 90)
        >>> 
        >>> F0 = def_gradient_stressfree_ini(params_cubic, params_tetra)
        >>> print("Transformation strain:")
        >>> print(F0)
```

### Notes
- Wrapper around def_gradient_stressfree()
        - Accepts lattice parameters directly
        - More convenient for quick calculations

---

## 60. `def_gradient`

### Signature
```python
def def_gradient(L1, L2, R)
```

### Description
Calculate total deformation gradient with rotation.
    Computes F = R · F₀ where R is rigid rotation and F₀ is lattice strain.

### Input Parameters
L1 (array 3×3): Initial lattice
        L2 (array 3×3): Final lattice
        R (array 3×3): Rigid body rotation matrix

### Output
numpy.ndarray (3×3): Total deformation gradient F

### Usage Examples
```python
>>> L1 = cubic_lattice_vec(3.0)
        >>> L2 = tetragonal_lattice_vec(3.0, 3.0, 4.0)
        >>> 
        >>> # With 45° rotation around Z
        >>> R = rotation_from_axis_angle([0,0,1], np.pi/4)
        >>> F = def_gradient(L1, L2, R)
        >>> print("Total deformation (with rotation):")
        >>> print(F)
```

### Notes
- F = R · F₀
        - R accounts for crystal rotation
        - F₀ is stress-free transformation
        - Used in variant selection

### Formula
F = R · L₂ · L₁⁻¹

---

## 61. `def_gradient_ini`

### Signature
```python
def def_gradient_ini(lattice_params_1, lattice_params_2, euler_angles)
```

### Description
Initialize deformation gradient with Euler angle rotation.
    Calculates F with rotation specified by Euler angles.

### Input Parameters
lattice_params_1, lattice_params_2 (tuple): Lattice parameters
        euler_angles (tuple): (φ₁, Φ, φ₂) in degrees

### Output
numpy.ndarray (3×3): Deformation gradient

### Usage Examples
```python
>>> params1 = (3.0, 3.0, 3.0, 90, 90, 90)
        >>> params2 = (3.0, 3.0, 4.0, 90, 90, 90)
        >>> euler = (45, 30, 60)  # degrees
        >>> 
        >>> F = def_gradient_ini(params1, params2, euler)
```

### Notes
- Combines lattice strain and rotation
        - Euler angles in degrees
        - Convenient initialization function

---

## 62. `def_gradient_ini2`

### Signature
```python
def def_gradient_ini2(L1, L2, axis, angle)
```

### Description
Initialize deformation gradient with axis-angle rotation.
    Alternative initialization using axis-angle for rotation.

### Input Parameters
L1, L2 (array 3×3): Lattice matrices
        axis (array [3]): Rotation axis
        angle (float): Rotation angle in radians

### Output
numpy.ndarray (3×3): Deformation gradient

### Usage Examples
```python
>>> L1 = cubic_lattice_vec(3.0)
        >>> L2 = tetragonal_lattice_vec(3.0, 3.0, 4.0)
        >>> 
        >>> # Rotate 60° around [111]
        >>> axis = np.array([1, 1, 1])
        >>> F = def_gradient_ini2(L1, L2, axis, np.radians(60))
```

### Notes
- Uses axis-angle for rotation
        - Alternative to Euler angles
        - Sometimes more intuitive

---

## 63. `niti_twinning`

### Signature
```python
def niti_twinning(variant='type1')
```

### Description
Get NiTi twinning parameters for specific variant.
    Returns crystallographic twinning elements for NiTi B19' martensite.

### Input Parameters
variant (str): Twinning variant ('type1', 'type2', 'compound')

### Output
dict: {
            'K1': array [3] - compound twin plane
            'eta1': array [3] - twinning direction
            'K2': array [3] - conjugate plane
            'eta2': array [3] - conjugate direction
            'shear': float - twinning shear
        }

### Usage Examples
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

### Notes
- Specific to NiTi B19' martensite
        - Type-I: {111} compound twin
        - Type-II: {001} twin
        - Data from experimental measurements

---

## 64. `get_twinningdata`

### Signature
```python
def get_twinningdata(L_parent, L_twin)
```

### Description
Extract complete twinning data from parent and twin lattices.
    Analyzes lattice relationship to determine all twinning elements.

### Input Parameters
L_parent (array 3×3): Parent lattice matrix
        L_twin (array 3×3): Twin lattice matrix

### Output
dict: Complete twin characterization including K1, K2, η1, η2, S

### Usage Examples
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

### Notes
- Comprehensive twin analysis
        - Extracts all crystallographic elements
        - Validates twin relationship
        - Used in texture simulations

---

## 65. `get_twinning_dislocation`

### Signature
```python
def get_twinning_dislocation(twin_data, b_parent)
```

### Description
Calculate twinning dislocation Burgers vector.
    Determines dislocation content at twin boundary.

### Input Parameters
twin_data (dict): From get_twinningdata()
        b_parent (array [3]): Burgers vector in parent

### Output
dict: {
            'b_residual': array [3] - residual Burgers vector
            'line_direction': array [3] - dislocation line direction
            'type': str - 'edge', 'screw', or 'mixed'
        }

### Usage Examples
```python
>>> twin_data = niti_twinning('type1')
        >>> b_parent = np.array([1., 0., 0.])  # Parent Burgers vector
        >>> 
        >>> disloc = get_twinning_dislocation(twin_data, b_parent)
        >>> print(f"Residual Burgers vector: {disloc['b_residual']}")
        >>> print(f"Dislocation type: {disloc['type']}")
```

### Notes
- Important for twin boundary structure
        - Affects boundary mobility
        - Related to interfacial energy

---

## 66. `gen_twinned_lattice_points`

### Signature
```python
def gen_twinned_lattice_points(L_parent, twin_plane_normal, n1=5, n2=5, n3=5)
```

### Description
Generate lattice points for parent and twinned regions.
    Creates atomic positions showing twin boundary.

### Input Parameters
L_parent (array 3×3): Parent lattice
        twin_plane_normal (array [3]): Twin plane (K1)
        n1, n2, n3 (int): Cell range

### Output
tuple: (parent_atoms, twin_atoms) - atomic positions

### Usage Examples
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

### Notes
- Shows twin boundary clearly
        - Used for visualization
        - Can calculate boundary energy

---

## 67. `write_txt`

### Signature
```python
def write_txt(filename, data, header='')
```

### Description
Write numerical data to text file.
    Saves array data in formatted text file with optional header.

### Input Parameters
filename (str): Output file path
        data (array): Numerical data to save
        header (str): Optional header line

### Output
None (writes to file)

### Usage Examples
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

### Notes
- Standard text format
        - Can be read by most software
        - Header optional but recommended
        - Uses space-separated values

---

## 68. `read_txt`

### Signature
```python
def read_txt(filename)
```

### Description
Read numerical data from text file.
    Loads data written by write_txt() or similar format.

### Input Parameters
filename (str): Input file path

### Output
numpy.ndarray: Loaded data

### Usage Examples
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

### Notes
- Skips comment lines starting with #
        - Returns numpy array
        - Compatible with write_txt output

---

## 69. `plane_line_intersection`

### Signature
```python
def plane_line_intersection(plane_point, plane_normal, line_point, line_direction)
```

### Description
Calculate intersection of plane and line in 3D.
    Finds point where line intersects plane, if it exists.

### Input Parameters
plane_point (array [3]): Point on plane
        plane_normal (array [3]): Plane normal vector
        line_point (array [3]): Point on line
        line_direction (array [3]): Line direction vector

### Output
numpy.ndarray [3] or None: Intersection point, or None if parallel

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> # (001) plane at z=5
        >>> plane_pt = np.array([0, 0, 5])
        >>> plane_n = np.array([0, 0, 1])
        >>> 
        >>> # Line from origin in [111] direction
        >>> line_pt = np.array([0, 0, 0])
        >>> line_dir = np.array([1, 1, 1])
        >>> 
        >>> intersection = plane_line_intersection(plane_pt, plane_n, line_pt, line_dir)
        >>> print(f"Intersection point: {intersection}")
        >>> # Should be [5, 5, 5]
```

### Notes
- Returns None if line parallel to plane
        - Returns None if line in plane
        - Used in geometric calculations

### Formula
t = n·(p₀ - l₀) / (n·d)
        intersection = l₀ + t·d

---

## 70. `flipvector`

### Signature
```python
def flipvector(v)
```

### Description
Flip vector to ensure positive first non-zero component.
    Standardizes vector direction for comparison.

### Input Parameters
v (array [3]): Vector

### Output
numpy.ndarray [3]: Flipped vector (if needed)

### Usage Examples
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

### Notes
- Makes first non-zero component positive
        - Useful for standardization
        - Preserves magnitude and line

---

## 71. `flipvector2negative`

### Signature
```python
def flipvector2negative(v)
```

### Description
Flip vector to ensure negative first non-zero component.
    Opposite of flipvector() - ensures negative first component.

### Input Parameters
v (array [3]): Vector

### Output
numpy.ndarray [3]: Flipped vector (if needed)

### Usage Examples
```python
>>> v = np.array([1, 2, 3])
        >>> v_neg = flipvector2negative(v)
        >>> print(v_neg)  # [-1, -2, -3]
```

### Notes
- Makes first non-zero component negative
        - Complementary to flipvector()

---

## 72. `vector2miller_ini`

### Signature
```python
def vector2miller_ini(v, L)
```

### Description
Convert Cartesian vector to Miller indices (initial guess).
    Provides starting point for vector→Miller conversion.

### Input Parameters
v (array [3]): Cartesian vector
        L (array 3×3): Lattice matrix

### Output
numpy.ndarray [3]: Initial Miller index guess

### Usage Examples
```python
>>> L = cubic_lattice_vec(3.0)
        >>> v_cart = np.array([3, 3, 0])
        >>> uvw_guess = vector2miller_ini(v_cart, L)
        >>> print(uvw_guess)  # Initial guess for [110]
```

### Notes
- Provides rough approximation
        - Use vectors2miller for refined result
        - Useful as starting point

---

## 73. `vectors2miller`

### Signature
```python
def vectors2miller(v, L, max_index=5)
```

### Description
Convert Cartesian vector to Miller indices with optimization.
    Finds best integer Miller indices representing given direction.

### Input Parameters
v (array [3]): Cartesian vector
        L (array 3×3): Lattice matrix
        max_index (int): Maximum Miller index to consider

### Output
numpy.ndarray [3]: Miller indices [uvw]

### Usage Examples
```python
>>> import numpy as np
        >>> 
        >>> L = cubic_lattice_vec(3.0)
        >>> v = np.array([3.0, 3.0, 0.0])
        >>> uvw = vectors2miller(v, L, max_index=3)
        >>> print(f"Miller indices: [{uvw[0]}{uvw[1]}{uvw[2]}]")
        >>> # Should give [110]
        >>> 
        >>> # Verify
        >>> v_reconstructed = uvw[0]*L[:,0] + uvw[1]*L[:,1] + uvw[2]*L[:,2]
        >>> angle = an_between_vecs(v, v_reconstructed)
        >>> print(f"Angle error: {angle:.6f}°")
```

### Notes
- Finds optimal integer representation
        - Limited by max_index
        - Uses least-squares fitting
        - More accurate than vector2miller_ini

### Formula
Minimize ||v - (u·a₁ + v·a₂ + w·a₃)||
        subject to u,v,w ∈ [-max_index, max_index]

---


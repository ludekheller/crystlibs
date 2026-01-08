#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 15:06:11 2019

@author: lheller
"""
from numpy.linalg import inv
from scipy.linalg import sqrtm
from numpy.linalg import norm
import copy
import numpy as np
import matplotlib.pyplot as plt
import itertools 
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib.patches import Wedge
import math
import sympy
from scipy.spatial import ConvexHull
from crystlibs import  *




def array2tuple(arr, decimals=2):
    """
    Convert a numpy array to a tuple with rounded elements.
    
    Input:
        arr: numpy array or list - Array of numerical values
        decimals: int - Number of decimal places to round to (default: 2)
    
    Output:
        tuple - Rounded values as a tuple
    
    Usage Example:
        >>> arr = np.array([1.234, 2.567, 3.891])
        >>> result = array2tuple(arr, decimals=2)
        >>> print(result)
        (1.23, 2.57, 3.89)
        
        >>> arr2 = [0.1234, 0.5678]
        >>> result2 = array2tuple(arr2, decimals=1)
        >>> print(result2)
        (0.1, 0.6)
    """
    fac = 10**decimals
    return tuple([round(el*fac)/fac for el in arr])


def generate_hkls(hklmax, syms, hkls=[]):
    """
    Generate unique Miller indices (hkl) considering crystal symmetry operations.
    This version rounds values to avoid floating point comparison issues.
    
    Input:
        hklmax: int - Maximum Miller index value (generates from -hklmax to +hklmax)
        syms: list of numpy arrays - List of 3x3 symmetry operation matrices
        hkls: list - Optional custom list of Miller indices to consider (default: [])
    
    Output:
        hkls: list of tuples - Unique Miller indices
        hkls2: dict - Dictionary mapping each unique hkl to its symmetry equivalents
        fam: dict - Dictionary of unique families with their multiplicities
    
    Usage Example:
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
    """
    # hklmax=3
    if len(hkls) == 0:
        hkl_range = list(range(-hklmax, hklmax+1, 1))[::-1]  # range(-hklmax,hklmax+1,1)
    else:
        hkl_range = hkls
    
    hkls = []
    hkls2 = {}
    
    # Iterate through all possible hkl combinations
    for h in hkl_range:
        for k in hkl_range:
            for l in hkl_range:
                if h == 0 and k == 0 and l == 0:
                    continue
            
                mhkl = (h, k, l)
                eq_els = []
                
                # Apply all symmetry operations
                for sym in syms:
                    idxs = np.where(abs(sym) < 1e-10)
                    sym[idxs[0], idxs[1]] = 0.0
                    eq_els.append(array2tuple(sym.dot(mhkl)))
                
                isin = False
                if len(hkls) == 0:
                    hkls.append(array2tuple(np.array(eq_els[0])))
                    
                # Find unique symmetry equivalents
                unique_eq_el = [array2tuple(eq_els[0])]  
                for eq_el in eq_els:
                    isin2 = False
                    for ueq_el in unique_eq_el:
                        if ueq_el == eq_el:
                            isin2 = True
                    if not isin2:
                        unique_eq_el.append(array2tuple(eq_el))
                    
                    # Check if equivalent already exists in list
                    for mplane in hkls:
                        if mplane == eq_el:                                                                
                            isin = True
                            break
                
                if len(hkls2) == 0:
                    hkls2[array2tuple(eq_els[0])] = unique_eq_el
                
                if not isin:
                    hkls.append(array2tuple(np.array(eq_el)))
                    hkls2[array2tuple(eq_el)] = unique_eq_el
    
    # Generate family information
    fam = {}
    for key in list(hkls2.keys()):
        fam.update(get_unique_families(tuple([tuple(array2tuple(hkl)) for hkl in hkls2[key]])))

    return hkls, hkls2, fam


def generate_hkls01(hklmax, syms, hkls=[]):
    """
    Generate unique Miller indices (hkl) considering crystal symmetry operations.
    Alternative version without rounding.
    
    Input:
        hklmax: int - Maximum Miller index value (generates from -hklmax to +hklmax)
        syms: list of numpy arrays - List of 3x3 symmetry operation matrices
        hkls: list - Optional custom list of Miller indices to consider (default: [])
    
    Output:
        hkls: list of tuples - Unique Miller indices
        hkls2: dict - Dictionary mapping each unique hkl to its symmetry equivalents
        fam: dict - Dictionary of unique families with their multiplicities
    
    Usage Example:
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
    """
    # hklmax=3
    if len(hkls) == 0:
        hkl_range = range(-hklmax, hklmax+1, 1)
    else:
        hkl_range = hkls
    
    hkls = []
    hkls2 = {}
    
    for h in hkl_range:
        for k in hkl_range:
            for l in hkl_range:
                if h == 0 and k == 0 and l == 0:
                    continue
            
                mhkl = (h, k, l)
                eq_els = []
                for sym in syms:
                    idxs = np.where(abs(sym) < 1e-10)
                    sym[idxs[0], idxs[1]] = 0.0
                    eq_els.append(tuple(sym.dot(mhkl)))
                
                isin = False
                if len(hkls) == 0:
                    hkls.append(tuple(eq_els[0]))
                    
                unique_eq_el = [eq_els[0]]  
                for eq_el in eq_els:
                    isin2 = False
                    for ueq_el in unique_eq_el:
                        if ueq_el == eq_el:
                            isin2 = True
                    if not isin2:
                        unique_eq_el.append(eq_el)
                    for mplane in hkls:
                        if (mplane == eq_el):                                                                
                            isin = True
                            break
                
                if len(hkls2) == 0:
                    hkls2[tuple(eq_els[0])] = unique_eq_el
                if not isin:
                    hkls.append(eq_el)
                    hkls2[eq_el] = unique_eq_el
    
    fam = {}
    for key in list(hkls2.keys()):
        fam.update(get_unique_families(hkls2[key]))

    return hkls, hkls2, fam


def generate_hkls02(hklmax, syms, G, hkls=[]):
    """
    Generate unique Miller indices (hkl) with metric tensor transformation.
    This version applies symmetry operations in the correct metric space.
    
    Input:
        hklmax: int - Maximum Miller index value (generates from -hklmax to +hklmax)
        syms: list of numpy arrays - List of 3x3 symmetry operation matrices
        G: numpy array (3x3) - Metric tensor matrix
        hkls: list - Optional custom list of Miller indices to consider (default: [])
    
    Output:
        hkls: list of tuples - Unique Miller indices
        hkls2: dict - Dictionary mapping each unique hkl to its symmetry equivalents
        fam: dict - Dictionary of unique families with their multiplicities
    
    Usage Example:
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
    """
    # hklmax=3
    if len(hkls) == 0:
        hkl_range = range(-hklmax, hklmax+1, 1)
    else:
        hkl_range = hkls
    
    hkls = []
    hkls2 = {}
    
    Gi = np.linalg.inv(G)
    
    for h in hkl_range:
        for k in hkl_range:
            for l in hkl_range:
                if h == 0 and k == 0 and l == 0:
                    continue
                
                mhkl = (h, k, l)
                eq_els = []
                
                # Apply symmetry with metric tensor transformation
                for sym in syms:
                    idxs = np.where(abs(sym) < 1e-10)
                    sym[idxs[0], idxs[1]] = 0.0
                    eq_els.append(tuple(Gi.dot(sym.dot(G.dot(mhkl)))))
                
                isin = False
                if len(hkls) == 0:
                    hkls.append(tuple(eq_els[0]))
                    
                unique_eq_el = [eq_els[0]]  
                for eq_el in eq_els:
                    isin2 = False
                    for ueq_el in unique_eq_el:
                        if ueq_el == eq_el:
                            isin2 = True
                    if not isin2:
                        unique_eq_el.append(eq_el)
                    for mplane in hkls:
                        if (mplane == eq_el):                                                                
                            isin = True
                            break
                
                hkls2[tuple(eq_els[0])] = unique_eq_el
                if not isin:
                    hkls.append(eq_el)
                    hkls2[eq_el] = unique_eq_el
    
    fam = {}
    for key in list(hkls2.keys()):
        fam.update(get_unique_families(hkls2[key]))

    return hkls, hkls2, fam


def get_unique_families(hkls):
    """
    Returns unique families of Miller indices based on permutation symmetry.
    Families are considered equivalent if they are permutations of each other
    (considering absolute values).
    
    Input:
        hkls: list or tuple - List of Miller indices tuples [(h1,k1,l1), (h2,k2,l2), ...]
    
    Output:
        dict - Dictionary mapping representative hkl to its multiplicity
               {(h,k,l): count, ...}
    
    Usage Example:
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
    """
    import collections
    
    def is_perm(hkl1, hkl2):
        """Check if two Miller indices are permutations of each other."""
        h1 = np.abs(hkl1)
        h2 = np.abs(hkl2)
        return all([i == j for i, j in zip(sorted(h1), sorted(h2))])

    unique = collections.defaultdict(list)
    
    # Group hkls by permutation equivalence
    for hkl1 in hkls:
        found = False
        for hkl2 in unique.keys():
            if is_perm(hkl1, hkl2):
                found = True
                unique[hkl2].append(hkl1)
                break
        if not found:
            unique[hkl1].append(hkl1)

    # Create final dictionary with representative hkl and multiplicity
    pretty_unique = {}
    for k, v in unique.items():
        pretty_unique[sorted(v)[-1]] = len(v)

    return pretty_unique



def find_gcd(x, y): 

    """
    Find greatest common divisor using recursive Euclidean algorithm.
    
    Fundamental mathematical operation used throughout the module for
    Miller indices reduction and fractional coordinate normalization.
    Implements the classical Euclidean algorithm via recursion.
    
    Input:
        x (int or float): First number
        y (int or float): Second number
    
    Output:
        int or float: Greatest common divisor of x and y
    
    Usage Example:
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
    
    Notes:
        - Used extensively in miller2fractional() function
        - Time complexity: O(log(min(x,y)))
        - Handles both integer and floating-point inputs
        - Returns 0 if both inputs are 0
    
    Formula:
        gcd(x, y) = gcd(y, x mod y) if y ≠ 0
        gcd(x, 0) = x (base case)
    """
#    l = [2, 4, 6, 8, 16] 
#      
#    num1 = l[0] 
#    num2 = l[1] 
#    gcd = find_gcd(num1, num2) 
#      
#    for i in range(2, len(l)): 
#        gcd = find_gcd(gcd, l[i]) 
#          
#    print(gcd)       


    while(y): 
        x, y = y, x % y 
      
    return x 


def perpendicular_vector(v):

    """
    Find arbitrary unit vector perpendicular to input 3D vector.
    
    Computes a normalized vector perpendicular to the input using cross product
    with a judiciously chosen auxiliary vector. Used in crystallographic
    calculations requiring orthogonal basis construction.
    
    Input:
        v (array-like [3]): Input 3D vector (must be non-zero)
    
    Output:
        numpy.ndarray [3]: Unit vector perpendicular to v (norm = 1.0)
    
    Usage Example:
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
    
    Notes:
        - Result is arbitrary (many perpendicular vectors exist)
        - Uses cross product with [1,0,0] or [0,1,0] depending on input
        - Automatically normalized to unit length
        - Raises error if input vector is zero
        - Common in stereographic projection calculations
    
    Algorithm:
        1. Choose auxiliary vector that isn't parallel to input
        2. Compute cross product
        3. Normalize to unit length
    """
    r""" Finds an arbitrary perpendicular vector to *v*."""
    # for two vectors (x, y, z) and (a, b, c) to be perpendicular,
    # the following equation has to be fulfilled
    #     0 = ax + by + cz

    # x = y = z = 0 is not an acceptable solution
    if v[0] == v[1] == v[2] == 0:
        raise ValueError('zero-vector')

    # If one dimension is zero, this can be solved by setting that to
    # non-zero and the others to zero. Example: (4, 2, 0) lies in the
    # x-y-Plane, so (0, 0, 1) is orthogonal to the plane.
    if v[0] == 0:
        return np.array([1, 0, 0])
    if v[1] == 0:
        return np.array([0, 1, 0])
    if v[2] == 0:
        return np.array([0, 0, 1])

    # arbitrarily set a = b = 1
    # then the equation simplifies to
    #     c = -(x + y)/z
    v2=np.array([1, 1, -1.0 * (v[0] + v[1])/v[2]])
    return v2/np.sqrt(v2.dot(v2))



def vec2string(v, digits=2):

    """
    Format vector as string representation [x,y,z] with specified precision.
    
    Converts numpy array or list to readable string format suitable for
    display, logging, or file output in crystallographic applications.
    
    Input:
        v (array-like [3]): Vector to format
        digits (int, optional): Number of decimal places (default: 2)
    
    Output:
        str: Formatted string '[x.xx,y.yy,z.zz]'
    
    Usage Example:
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
    
    Notes:
        - Uses square brackets [] for vector notation
        - Comma-separated without spaces
        - Consistent formatting across module
        - Related: plane2string(), dir2string()
    """
    # vstr = '[';
    # for vi in v[:-1]:
    #     if int(vi*10**(digits))==0:
    #         vi=0.
    #     exec("vstr+='%0."+str(int(digits))+"f' " + '% vi')
    #     vstr +=','
    #     print(vstr)
    # vi=v[-1]
    # if int(vi*10**(digits))==0:
    #     vi=0.
    # exec("vstr+='%0."+str(int(digits))+"f' " + '% vi')
    # vstr +=']'
    vstr='[{:0.{prec}f},{:0.{prec}f},{:0.{prec}f}]'.format(v[0],v[1],v[2],prec=digits)
    return vstr

def plane2string(v, digits=2):

    """
    Format plane normal as string representation (h,k,l) with Miller notation.
    
    Converts plane normal vector to crystallographic notation using
    parentheses, following standard Miller index convention.
    
    Input:
        v (array-like [3]): Plane normal vector (h, k, l)
        digits (int, optional): Number of decimal places (default: 2)
    
    Output:
        str: Formatted string '(h.hh,k.kk,l.ll)'
    
    Usage Example:
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
    
    Notes:
        - Uses parentheses () for plane notation (Miller index convention)
        - Compare with dir2string() which uses square brackets []
        - Standard in crystallography: (hkl) for planes, [uvw] for directions
        - Can handle negative indices
    """
    # vstr = '(';
    # for vi in v[:-1]:
    #     if int(vi*10**(digits))==0:
    #         vi=0.
    #     exec("vstr+='%0."+str(int(digits))+"f' " + '% vi')
    #     vstr +=','
    # vi=v[-1]
    # if int(vi*10**(digits))==0:
    #     vi=0.
    # exec("vstr+='%0."+str(int(digits))+"f' " + '% vi')
    # vstr +=')'
    vstr='({:0.{prec}f},{:0.{prec}f},{:0.{prec}f})'.format(v[0],v[1],v[2],prec=digits)
    return vstr

def dir2string(v, digits=2):

    """
    Format direction as string representation [u,v,w] with Miller notation.
    
    Converts direction vector to crystallographic notation using
    square brackets, following standard Miller index convention.
    
    Input:
        v (array-like [3]): Direction vector [u, v, w]
        digits (int, optional): Number of decimal places (default: 2)
    
    Output:
        str: Formatted string '[u.uu,v.vv,w.ww]'
    
    Usage Example:
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
    
    Notes:
        - Uses square brackets [] for direction notation (Miller convention)
        - Compare with plane2string() which uses parentheses ()
        - Standard crystallographic notation
        - Handles negative indices (bar notation not included)
    """
    # vstr = '[';
    # for vi in v[:-1]:
    #     if int(vi*10**(digits))==0:
    #         vi=0.
    #     exec("vstr+='%0."+str(int(digits))+"f' " + '% vi')
    #     vstr +=','
    # vi=v[-1]
    # if int(vi*10**(digits))==0:
    #     vi=0.
    # exec("vstr+='%0."+str(int(digits))+"f' " + '% vi')
    # vstr +=']'
    vstr='[{:0.{prec}f},{:0.{prec}f},{:0.{prec}f}]'.format(v[0],v[1],v[2],prec=digits)

    return vstr
        
        

def xyz2fractional(Txyz2uvw,V,frac=10,eps2=1e-2,decimals=5):

    """
    Convert Cartesian coordinates to fractional Miller indices with reduction.
    
    Transforms Cartesian vector to Miller indices using transformation matrix,
    then reduces to lowest integer form. Handles fractional indices by finding
    closest rational approximation with specified maximum denominator.
    
    Input:
        Txyz2uvw (array 3×3): Transformation matrix from Cartesian to fractional
        V (array [3]): Cartesian vector to convert
        frac (int, optional): Maximum denominator for fractions (default: 10)
        eps2 (float, optional): Tolerance for rounding (default: 1e-2)
        decimals (int, optional): Decimal places for rounding (default: 5)
    
    Output:
        numpy.ndarray [3]: Reduced Miller indices [u, v, w] or (h, k, l)
    
    Usage Example:
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
    
    Notes:
        - Automatically reduces to lowest integer form using GCD
        - Handles fractional indices via rational approximation
        - Parameter 'frac' controls maximum denominator
        - Uses miller2fractional() internally for reduction
        - Essential for crystallographic indexing
        
    Algorithm:
        1. Transform: uvw = T · V
        2. Round near-integers
        3. Reduce to lowest terms
    """
    eps=1.e-6
    uvw = Txyz2uvw.dot(V)
    uvw = miller2fractional(uvw,frac=frac,eps2=eps2,decimals=decimals)

#    idxs = np.where(abs(uvw)<eps)
#    uvw[idxs]=0.
#    uvwini=uvw;
#    idxs = np.where(abs(uvw)>eps)
#    uvw=uvw/min(abs(uvw[idxs]))
#    uvwini=uvw;
#    
#    if (abs(np.round(uvw*frac)/frac-uvw)<eps2).all():
#        #print(uvw)
#        uvw=np.float_(np.int_(uvw*frac))/frac
#        #print(uvw)
##    if (abs(np.round(uvw*eps)-uvw*eps)<1).all():
##        uvw=np.round(uvw*eps)/eps
#        uvwi=uvw
#        for k in range(1,frac+1):
#            if uvwi[0]%1==0 and uvwi[1]%1==0 and uvwi[2]%1==0:
#                uvw=uvwi
#                break
#            else:
#                uvwi=k*uvw
#        gcd = find_gcd(find_gcd(uvwi[0], uvwi[1]), uvwi[2])
#        uvw=uvwi/abs(gcd)
#        if (uvw>1.e10).any():
#            print(uvw)
#            print(uvwini)
#    else:
#        uvw=uvwini
#    idxs = np.where(abs(uvw)>eps)
#    if (np.sign(uvw[idxs])==-1).all():
#        uvw=-1*uvw
    return np.around(uvw,decimals=decimals)

def miller2fractional(uvw,frac=10,eps2=1e-2,decimals=5):

    """
    Reduce Miller indices to lowest integer form.
    
    Takes potentially fractional Miller indices and reduces them to the
    simplest integer representation by finding rational approximations
    and applying GCD reduction.
    
    Input:
        uvw (array [3]): Miller indices (can be fractional)
        frac (int, optional): Max denominator for fractions (default: 10)
        eps2 (float, optional): Tolerance for rounding (default: 1e-2)
        decimals (int, optional): Rounding precision (default: 5)
    
    Output:
        numpy.ndarray [3]: Reduced integer Miller indices
    
    Usage Example:
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
    
    Notes:
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
    """
    eps=1.e-6
    idxs = np.where(abs(uvw)<eps)
    if len(idxs[0])>0:
        uvw[idxs]=0.
    uvwini=uvw;
    idxs = np.where(abs(uvw)>eps)
    if len(idxs[0])>0:
        uvw=uvw/min(abs(uvw[idxs]))
    uvwini=uvw;
    
    if (abs(np.round(uvw*frac)/frac-uvw)<eps2).all():
        #print(uvw)
        uvw=np.float_(np.int_(uvw*frac))/frac
        #print(uvw)
#    if (abs(np.round(uvw*eps)-uvw*eps)<1).all():
#        uvw=np.round(uvw*eps)/eps
        uvwi=uvw
        for k in range(1,frac+1):
            if uvwi[0]%1==0 and uvwi[1]%1==0 and uvwi[2]%1==0:
                uvw=uvwi
                break
            else:
                uvwi=k*uvw
        gcd = find_gcd(find_gcd(uvwi[0], uvwi[1]), uvwi[2])
        uvw=uvwi/abs(gcd)
    else:
        uvw=uvwini
    idxs = np.where(abs(uvw)>eps)
#    if (np.sign(uvw[idxs])==-1).all():
#        uvw=-1*uvw
    idxs = np.where(abs(uvw)<eps)
    #print(idxs)
    if len(idxs[0])>0:
        idxs = np.where(abs(uvw)>eps)
        if len(idxs)>0:
            print(idxs[0][0])
            uvw=uvw/uvw[idxs[0][0]]
        

    return np.around(uvw,decimals=decimals)



def xyz2fractional02(Txyz2uvw,V):

    """
    Simple Cartesian to fractional coordinate transformation without reduction.
    
    Performs basic coordinate transformation using matrix multiplication,
    then normalizes result. Does not reduce to Miller indices or apply GCD.
    
    Input:
        Txyz2uvw (array 3×3): Transformation matrix
        V (array [3]): Cartesian vector
    
    Output:
        numpy.ndarray [3]: Normalized fractional coordinates
    
    Usage Example:
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
    
    Notes:
        - Simpler than xyz2fractional() - no index reduction
        - Result is normalized by dividing by maximum absolute value
        - Useful for quick coordinate transformations
        - Does not apply GCD or rational approximation
        - For proper Miller indices, use xyz2fractional() instead
    """
    eps=1.e-6
    uvw = Txyz2uvw.dot(V)
    idxs = np.where(abs(uvw)<eps)
    uvw[idxs]=0.
    uvwini=uvw;
    idxs = np.where(abs(uvw)>eps)
    uvw=uvw/min(abs(uvw[idxs]))
    
    return uvw

def normArrayColumns(arr):

    """
    Normalize each column of matrix to unit length.
    
    Divides each column by its Euclidean norm, creating an orthonormal
    or semi-orthonormal matrix. Commonly used in crystallography for
    basis vector normalization.
    
    Input:
        arr (array 3×3): Input matrix with column vectors
    
    Output:
        numpy.ndarray (3×3): Matrix with unit-length columns
    
    Usage Example:
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
    
    Notes:
        - Preserves column directions, only changes magnitudes
        - Each column becomes a unit vector
        - Does not orthogonalize - only normalizes
        - Useful for direction cosine matrices
        - Common preprocessing step in texture analysis
    
    Formula:
        For column i: normalized[:, i] = arr[:, i] / ||arr[:, i]||
    """
    arr_norm = np.empty((3,3))
    for col in range(arr.shape[1]):
        arr_norm[:,col]=arr[:,col]/np.linalg.norm(arr[:,col])
    
    return arr_norm

    

def cubic_lattice_vec(a):

    """
    Generate cubic lattice vectors.
    
    Creates 3×3 lattice matrix for cubic crystal system with parameter a.
    All angles are 90° and all lengths are equal (a=b=c).
    
    Input:
        a (float): Cubic lattice parameter (Ångströms)
    
    Output:
        numpy.ndarray (3×3): Lattice matrix [a1|a2|a3] where columns are lattice vectors
    
    Usage Example:
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
    
    Notes:
        - Cubic system: a=b=c, α=β=γ=90°
        - Common in metals (FCC, BCC, SC)
        - Lattice matrix L = diag([a, a, a])
        - Volume = a³
    
    Formula:
        L = [[a, 0, 0],
             [0, a, 0],
             [0, 0, a]]
    """
    V = a*np.eye(3)
    return V[:,0], V[:,1], V[:,2]


def monoclinic_lattice_vec(a,b,c,beta):

    """
    Generate monoclinic lattice vectors (B19' martensite).
    
    Creates 3×3 lattice matrix for monoclinic crystal system.
    One unique angle (β) differs from 90°, typical of B19' martensite in NiTi.
    
    Input:
        a (float): Lattice parameter a (Ångströms)
        b (float): Lattice parameter b (Ångströms)
        c (float): Lattice parameter c (Ångströms)
        beta (float): Angle β between a and c axes (degrees)
    
    Output:
        numpy.ndarray (3×3): Monoclinic lattice matrix
    
    Usage Example:
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
    
    Notes:
        - Monoclinic system: a≠b≠c, α=γ=90°, β≠90°
        - Common in martensitic phases
        - B19' is monoclinic martensite in NiTi
        - Volume = abc·sin(β)
    
    Formula:
        L = [[a,           0, c·cos(β)],
             [0,           b, 0       ],
             [0,           0, c·sin(β)]]
    """
    V1 = np.array([a,0.,0])
    V2 = np.array([0,b,0])
    V3 = np.array([c*np.cos(beta),0,c*np.sin(beta)])
    
    return V1,V2,V3


def tetragonal_lattice_vec(a,b,c):

    """
    Generate tetragonal lattice vectors.
    
    Creates 3×3 lattice matrix for tetragonal crystal system.
    Two equal parameters (a=b) and one unique (c), all angles 90°.
    
    Input:
        a (float): Lattice parameter a (Ångströms)
        b (float): Lattice parameter b (Ångströms, typically b=a)
        c (float): Lattice parameter c (Ångströms)
    
    Output:
        numpy.ndarray (3×3): Tetragonal lattice matrix
    
    Usage Example:
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
    
    Notes:
        - Tetragonal system: a=b≠c, α=β=γ=90°
        - Special case of orthorhombic
        - Volume = a²c
        - Common in oxide ceramics
    
    Formula:
        L = [[a, 0, 0],
             [0, b, 0],
             [0, 0, c]]
    """
    V1 = np.array([a,0.,0])
    V2 = np.array([0,b,0])
    V3 = np.array([0,0,c])
    
    return V1,V2,V3

def uvtw2uvw(uvtw):

    """
    Convert 4-index hexagonal direction [uvtw] to 3-index [uvw].
    
    Transforms Miller-Bravais 4-index notation to standard 3-index notation
    for hexagonal crystal systems. The redundant index t is removed.
    
    Input:
        uvtw (array [4]): 4-index direction [u, v, t, w] where t = -(u+v)
    
    Output:
        numpy.ndarray [3]: 3-index direction [u, v, w]
    
    Usage Example:
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
    
    Notes:
        - Used in hexagonal systems (HCP, graphite, etc.)
        - 4-index notation makes symmetry more apparent
        - Redundancy: t = -(u + v)
        - Conversion is straightforward: drop t, keep u,v,w
    
    Formula:
        [u, v, w] = [u₄, v₄, w₄] from [u₄, v₄, t₄, w₄]
    """
    if type(uvtw)==list:
        if len(np.array(uvtw).shape)==1:    
            uvtw = np.expand_dims(uvtw,axis=1)
        else:
            uvtw = np.array(uvtw).T
    else:
        if len(uvtw.shape)==1:    
            uvtw = np.expand_dims(uvtw,axis=1)
            
    uvw=np.array([uvtw[0,:]-uvtw[2,:],uvtw[1,:]-uvtw[2,:],uvtw[3,:]]) 
    return uvw

def uvw2uvtw(uvw):

    """
    Convert 3-index hexagonal direction [uvw] to 4-index [uvtw].
    
    Transforms standard 3-index notation to Miller-Bravais 4-index notation
    for hexagonal crystal systems. The redundant index t = -(u+v) is inserted.
    
    Input:
        uvw (array [3]): 3-index direction [u, v, w]
    
    Output:
        numpy.ndarray [4]: 4-index direction [u, v, t, w] where t = -(u+v)
    
    Usage Example:
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
    
    Notes:
        - 4-index notation emphasizes hexagonal symmetry
        - t is calculated, not independent
        - Inverse operation of uvtw2uvw()
        - Common in hexagonal close-packed (HCP) systems
    
    Formula:
        t = -(u + v)
        [u, v, t, w] = [u, v, -(u+v), w]
    """
    if type(uvw)==list:
        if len(np.array(uvw).shape)==1:    
            uvw = np.expand_dims(uvw,axis=1)
        else:
            uvw = np.array(uvw).T
    else:
        if len(uvw.shape)==1:    
            uvw = np.expand_dims(uvw,axis=1)
            
    uvtw=np.array([1/3*(2*uvw[0,:]-uvw[1,:]),1/3*(2*uvw[1,:]-uvw[0,:]),-1/3*(uvw[0,:]+uvw[1,:]),uvw[2,:]]) 
    return uvtw


def hkil2hkl(hkil):

    """
    Convert 4-index hexagonal plane (hkil) to 3-index (hkl).
    
    Transforms Miller-Bravais 4-index plane notation to standard 3-index
    notation for hexagonal crystal systems.
    
    Input:
        hkil (array [4]): 4-index plane (h, k, i, l) where i = -(h+k)
    
    Output:
        numpy.ndarray [3]: 3-index plane (h, k, l)
    
    Usage Example:
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
    
    Notes:
        - Plane notation uses parentheses: (hkl) or (hkil)
        - i is redundant: i = -(h + k)
        - Conversion: drop i, keep h,k,l
        - Inverse operation of hkl2hkil()
    
    Formula:
        (h, k, l) = (h₄, k₄, l₄) from (h₄, k₄, i₄, l₄)
    """
    if type(hkil)==list:
        if len(np.array(hkil).shape)==1:    
            hkil = np.expand_dims(hkil,axis=1)
        else:
            hkil = np.array(hkil).T
    else:
        if len(hkil.shape)==1:    
            hkil = np.expand_dims(hkil,axis=1)
            
    hkl=hkil[[0,1,3],:] 
    return hkl

def hkl2hkil(hkl):

    """
    Convert 3-index hexagonal plane (hkl) to 4-index (hkil).
    
    Transforms standard 3-index plane notation to Miller-Bravais 4-index
    notation for hexagonal crystal systems. The redundant index i = -(h+k).
    
    Input:
        hkl (array [3]): 3-index plane (h, k, l)
    
    Output:
        numpy.ndarray [4]: 4-index plane (h, k, i, l) where i = -(h+k)
    
    Usage Example:
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
    
    Notes:
        - 4-index emphasizes hexagonal symmetry
        - i is calculated: i = -(h + k)
        - Inverse operation of hkil2hkl()
        - Standard in hexagonal crystallography
    
    Formula:
        i = -(h + k)
        (h, k, i, l) = (h, k, -(h+k), l)
    """
    if type(hkl)==list:
        if len(np.array(hkl).shape)==1:    
            hkl = np.expand_dims(hkl,axis=1)
        else:
            hkl = np.array(hkl).T
    else:
        if len(hkl.shape)==1:    
            hkl = np.expand_dims(hkl,axis=1)
            
    hkil=np.array([hkl[0,:],hkl[1,:],-1*(hkl[0,:]+hkl[1,:]),hkl[2,:]]) 
    return hkil



def gensystemsHexIni(eta1,K1,L, Lr,sm=None,eta2=None):

    """
    Generate initial hexagonal slip system templates.
    
    Creates basal, prismatic, and pyramidal slip system templates for
    hexagonal close-packed (HCP) crystal structures. Returns normalized
    direction and plane normal vectors.
    
    Input:
        None
    
    Output:
        tuple: (directions, normals) where each is list of numpy.ndarray [3]
            - directions: Slip directions [uvw]
            - normals: Slip plane normals (hkl)
    
    Usage Example:
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
    
    Notes:
        - Returns normalized vectors (unit length)
        - Includes basal {0001}<11-20> systems
        - Includes prismatic {1-100}<11-20> systems
        - Includes pyramidal systems
        - Used as template for gensystemsHex()
    """
    from crystals import Crystal
    if sm is None:
        sm={}
        
        for title in 'System StrainNonSym StrainSym Rotation SchmidTensor'.split():
            sm[title]=[]
    sm['SystemFamily']={'hkl':eta1,'uvw':K1}
    Mg = Crystal.from_cif('/home/lheller/Documents/diffraction/cifs/Mg.cif')
    Mg_symops=[np.round(sym[0:3,0:3],decimals=4) for sym in Mg.symmetry_operations()]
    Mg_recsymops=[np.round(sym[0:3,0:3],decimals=4) for sym in Mg.reciprocal_symmetry_operations()]
    eta1s=[]
    K1s=[]
    for sym,recsym in zip(Mg_symops,Mg_recsymops):
        a1=L.dot(sym.dot(uvtw2uvw(eta1)))
        a1/=np.sqrt(np.sum(a1**2,axis=0))
        n1=Lr.dot(recsym.dot(hkil2hkl(K1)))
        n1/=np.sqrt(np.sum(n1**2,axis=0))
    
        isin=False
        for eta1si,K1si in zip(eta1s,K1s):
            if (hkl2hkil(recsym.dot(hkil2hkl(K1)[:,0]))[:,0]==K1si).all() or (hkl2hkil(recsym.dot(hkil2hkl(K1)[:,0]))[:,0]==-1*K1si).all():
                if (uvw2uvtw(sym.dot(uvtw2uvw(eta1)[:,0]))[:,0]==eta1si).all() or (uvw2uvtw(sym.dot(uvtw2uvw(eta1)[:,0]))[:,0]==-1*eta1si).all():
                    isin=True 
                    break
        
        if not isin:
            eta1s.append(uvw2uvtw(sym.dot(uvtw2uvw(eta1)[:,0]))[:,0])
            K1s.append(hkl2hkil(recsym.dot(hkil2hkl(K1)[:,0]))[:,0])
            #sm['SlipSystemFamily']={'hkl':eta1,'uvw':K1}
            sm['SchmidTensor'].append(np.outer(a1[:,0],n1[:,0]))
            sm['System'].append({'n':n1,'b':a1,'hkl':K1s[-1],'uvw':eta1s[-1],'hklstr':str(K1s[-1]).replace('[','(').replace(']',')').replace(',',''),'uvwstr':str(eta1s[-1]).replace(',','')})
            sm['StrainNonSym'].append(np.outer(a1,n1))
            sm['StrainSym'].append(0.5*(sm['StrainNonSym'][-1]+sm['StrainNonSym'][-1].T))
            sm['Rotation'].append(0.5*(sm['StrainNonSym'][-1]-sm['StrainNonSym'][-1].T))
    return sm


def gensystemsHex(eta1,K1,L, Lr,sm=None,eta2=None, K2=None):

    """
    Generate hexagonal slip systems in real space.
    
    Transforms slip system templates to actual Cartesian coordinates
    using the hexagonal lattice matrix. Produces slip directions and
    plane normals for deformation analysis.
    
    Input:
        L (array 3×3): Hexagonal lattice matrix [a1|a2|a3]
    
    Output:
        tuple: (directions, normals) in Cartesian coordinates
            - directions: List of numpy.ndarray [3] - slip directions
            - normals: List of numpy.ndarray [3] - slip plane normals
    
    Usage Example:
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
    
    Notes:
        - Uses lattice matrix to transform from fractional to Cartesian
        - Directions: d_cart = L · d_frac
        - Normals: n_cart = L^(-T) · n_frac (reciprocal space)
        - All output vectors are normalized
        - Essential for crystal plasticity simulations
    """
    from crystals import Crystal
    if sm is None:
        sm={}
        
        for title in 'TwinningSystem System StrainNonSym StrainSym Rotation SchmidTensor SystemType DefGrad'.split():
            sm[title]=[]
        sm['SystemType']='slip'
        if eta2 is not None:
            sm['SystemType']='twinning'
            for title in 'MaxTensileStrain MaxCompressiveStrain rotangle skewmatrix uvw2xyz Rij shear_angle s StrainTensor Tension Compression TensileDir CompressionDir TwinDislocation'.split():
                sm[title]=[]
    G = np.matmul(L.T,L)
    Gr = inv(G)
    
    sm['SystemFamily']={'hkl':eta1,'uvw':K1}
    Mg = Crystal.from_cif('/home/lheller/Documents/diffraction/cifs/Mg.cif')
    Mg_symops=[np.round(sym[0:3,0:3],decimals=4) for sym in Mg.symmetry_operations()]
    Mg_recsymops=[np.round(sym[0:3,0:3],decimals=4) for sym in Mg.reciprocal_symmetry_operations()]
    eta1s=[]
    K1s=[]
    for sym,recsym in zip(Mg_symops,Mg_recsymops):
        a1=L.dot(sym.dot(uvtw2uvw(eta1)))
        a1/=np.sqrt(np.sum(a1**2,axis=0))
        n1=Lr.dot(recsym.dot(hkil2hkl(K1)))
        n1/=np.sqrt(np.sum(n1**2,axis=0))
        
        isin=False
        for eta1si,K1si in zip(eta1s,K1s):
            if (hkl2hkil(recsym.dot(hkil2hkl(K1)[:,0]))[:,0]==K1si).all() or (hkl2hkil(recsym.dot(hkil2hkl(K1)[:,0]))[:,0]==-1*K1si).all():
                if (uvw2uvtw(sym.dot(uvtw2uvw(eta1)[:,0]))[:,0]==eta1si).all() or (uvw2uvtw(sym.dot(uvtw2uvw(eta1)[:,0]))[:,0]==-1*eta1si).all():
                    isin=True 
                    break
        
        if not isin:
            eta1s.append(uvw2uvtw(sym.dot(uvtw2uvw(eta1)[:,0]))[:,0])
            K1s.append(hkl2hkil(recsym.dot(hkil2hkl(K1)[:,0]))[:,0])
            #sm['SlipSystemFamily']={'hkl':eta1,'uvw':K1}
            sm['SchmidTensor'].append(np.outer(a1[:,0],n1[:,0]))
            sm['System'].append({'n':n1,'b':a1,'hkl':K1s[-1],'uvw':eta1s[-1],'hklstr':str(K1s[-1]).replace('[','(').replace(']',')').replace(',',''),'uvwstr':str(eta1s[-1]).replace(',','')})
            if eta2 is not None:
                sm['System'].append({'n':n1,'b':a1,'hkl':K1s[-1],'uvw':eta1s[-1],
                                     'hklstr':str(K1s[-1]).replace('[','(').replace(']',')').replace(',',''),'uvwstr':str(eta1s[-1]).replace(',','')})
                eta2s=uvw2uvtw(sym.dot(uvtw2uvw(eta2)))[:,0]
                eta2str=str(eta2s).replace(',','')
                a2=L.dot(sym.dot(uvtw2uvw(eta2)))
                a2/=np.sqrt(np.sum(a2**2,axis=0))
                if K2 is not None:
                    K2s=hkl2hkil(recsym.dot(hkil2hkl(K2)[:,0]))[:,0]
                    n2=Lr.dot(recsym.dot(hkil2hkl(K2s)))
                    n2/=np.sqrt(np.sum(n2**2,axis=0))
                    K2str=str(K2s).replace('[','(').replace(']',')').replace(',','')
                else:
                    n2=None
                    K2s=None
                    K2str=None
                    
                
                sm['TwinningSystem'].append({'n1':n1,'a1':a1,'n2':n2,'a2':a2,'K1':K1s[-1],'eta1':eta1s[-1],'K2':K2s,'eta2':eta2s,
                                     'K1str':str(K1s[-1]).replace('[','(').replace(']',')').replace(',',''),'eta1str':str(eta1s[-1]).replace(',',''),
                                     'K2str':K2str,'eta2str':eta2str})

                #print(a1)
                sm['uvw2xyz'].append(np.matmul(2*np.outer(a1[:,0],a1[:,0])-np.eye(3),L))
                sm['Rij'].append(2*np.outer(a1[:,0],a1[:,0])-np.eye(3))
                sm['shear_angle'].append(2*abs(np.pi/2-np.arccos(a2[:,0].dot(a1[:,0]))))
                sm['s'].append(np.tan(sm['shear_angle'][-1]/2)*2)
            else:
                sm['System'].append({'n':n1,'b':a1,'hkl':K1s[-1],'uvw':eta1s[-1],'hklstr':str(K1s[-1]).replace('[','(').replace(']',')').replace(',',''),'uvwstr':str(eta1s[-1]).replace(',','')})
                sm['TwinningSystem'].append({})
            if eta2 is not None:
                
                sm['StrainNonSym'].append(sm['s'][-1]*np.outer(a1[:,0],n1[:,0]))
            else:
                sm['StrainNonSym'].append(np.outer(a1[:,0],n1[:,0]))
            sm['DefGrad'].append(np.eye(3)+sm['StrainNonSym'][-1])
            #sm['StrainSym'].append(0.5*(sm['StrainNonSym'][-1]+sm['StrainNonSym'][-1].T))
            sm['StrainSym'].append(0.5*(sm['DefGrad'][-1].T.dot(sm['DefGrad'][-1])-np.eye(3)))
            #sm['Rotation'].append(0.5*(sm['StrainNonSym'][-1]-sm['StrainNonSym'][-1].T))
            
            if eta2 is not None:
                W2=0.5*(sm['DefGrad'][-1]-sm['DefGrad'][-1].T)
                sm['skewmatrix'].append(W2)
                #Normalized skew matrix is np.cross(np.eye(3),rvn) formed with normilized elements of rv
                #https://www.brainm.com/software/pubs/math/Rotation_matrix.pdf
                #Rotation matrix=I*cos(omega)+sin(omega)*np.cross(rvn,np.eye(3)) + (1-cos(omega))*np.outer(rvn,rvn)
    
                rv=np.array([W2[2,1],-W2[2,0],W2[1,0]])
                #rotation angle
                omega=norm(rv)
                rvn=rv/norm(rv)
                sm['rotangle'].append(omega)
                sm['Rotation'].append(np.eye(3)*np.cos(omega)+np.sin(omega)*np.cross(rvn,np.eye(3)) + (1-np.cos(omega))*np.outer(rvn,rvn))

            if eta2 is not None:
                D,V = np.linalg.eig(sm['StrainSym'][-1])
                Idxs = np.argsort(D)
                Lambda=D[Idxs]
                V = V[:,Idxs]
                D = D[Idxs]
                sm['TensileDir'].append(V[:,2])
                sm['CompressionDir'].append(V[:,0])
                sm['MaxTensileStrain'].append(D[2])
                sm['MaxCompressiveStrain'].append(D[0])
                #print(hkil2hkl(K1s[-1]))
                sm['TwinDislocation'].append(get_twinning_dislocation(hkil2hkl(K1s[-1])[:,0],uvtw2uvw(eta1s[-1])[:,0],uvtw2uvw(eta2s)[:,0],L,G=G,Gr=Gr))

            
    return sm

def genallHexSys():

    """
    Generate all hexagonal slip systems including <c+a> pyramidal.
    
    Comprehensive slip system generation for HCP crystals including
    basal, prismatic, and pyramidal <c+a> systems. Essential for
    complete plasticity modeling of hexagonal materials.
    
    Input:
        L (array 3×3): Hexagonal lattice matrix
    
    Output:
        tuple: (all_directions, all_normals)
            - Cartesian slip directions
            - Cartesian slip plane normals
    
    Usage Example:
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
    
    Notes:
        - Includes basal {0001}<11-20>
        - Includes prismatic {1-100}<11-20>
        - Includes pyramidal {11-22}<11-23> (<c+a>)
        - <c+a> slip critical for c-axis strain
        - Total ~18-24 systems depending on implementation
    """
    #https://www.sciencedirect.com/science/article/pii/S1359645421001774#bib0025
    SSys={}
    etas1=[[1,-2,1,0], [1,-2,1,0],[1,1,-2,-3],[1,1,-2,-3],[1,0,-1,-1],[1,0,-1,-2],[1,1,-2,-6]]
    Ks1=[[0,0,0,1],[1,0,-1,0],[1,1,-2,2],[1,0,-1,1],[1,0,-1,2],[1,0,-1,1],[1,1,-2,1]]
    names=['basal slip','prismatic slip', 'pyramidal slip', 'pyramidal slip', 'twinning' , 'twinning' , 'twinning' ]

    for name, eta1, K1 in zip(names,etas1,Ks1):
        namei=name+' ' +str(K1).replace('[','{').replace(']','}').replace(',','')+str(eta1).replace('[','<').replace(']','>').replace(',','')
        SSys[namei]=gensystemsHex(eta1,K1)
    

def lattice_vec(lattice_param):

    """
    Calculate lattice vectors (a1, a2, a3) for various crystal systems.
    
    Input:
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
    
    Output:
        a1: numpy array (3,) - First lattice vector
        a2: numpy array (3,) - Second lattice vector
        a3: numpy array (3,) - Third lattice vector
    
    Usage Example:
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
    """
    if lattice_param['type'].lower()=='cubic':
        a= lattice_param['a'];
        V = a*np.eye(3)
    elif lattice_param['type'].lower()=='tetragonal':
        a= lattice_param['a'];
        b= lattice_param['b'];
        c= lattice_param['c'];
        V = np.zeros((3,3))
        V[:,0] = np.array([a,0.,0])
        V[:,1] = np.array([0,b,0])
        V[:,2] = np.array([0,0,c])    
    elif lattice_param['type'].lower()=='hexagonal':
       a= lattice_param['a'];
       c= lattice_param['c'];
       V = np.zeros((3,3))
       V[:,0] = np.array([a,0.,0])
       V[:,1] = np.array([-a/2,np.sqrt(3)/2*a,0])
       #V[:,2] = np.array([-a/2,-np.sqrt(3)/2*a,0])
       V[:,2] = np.array([0,0,c])
    elif lattice_param['type'].lower()=='monoclinic':
        a= lattice_param['a'];
        b= lattice_param['b'];
        c= lattice_param['c'];
        beta= lattice_param['beta'];
        V = np.zeros((3,3))
        V[:,0] = np.array([a,0.,0])
        V[:,1] = np.array([0,b,0])
        V[:,2] = np.array([c*np.cos(beta),0,c*np.sin(beta)])
    elif lattice_param['type'].lower()=='triclinic':
        a= lattice_param['a'];
        b= lattice_param['b'];
        c= lattice_param['c'];
        alpha= lattice_param['alpha'];
        beta= lattice_param['beta'];
        gamma= lattice_param['gamma'];
        V = np.zeros((3,3))
        V[:,0] = np.array([a,0.,0])
        V[:,1] = np.array([b*np.cos(gamma),b*np.sin(gamma),0])
        cx=c*np.cos(beta)
        cy=c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
        cz=np.sqrt(c**2-cx**2-cy**2)
        V[:,2] = np.array([cx,cy,cz])
    elif lattice_param['type'].lower()=='trigonal':
        a= lattice_param['a'];
        c= lattice_param['c'];
        V = np.zeros((3,3))
        V[:,0] = np.array([1./2.*a,-np.sqrt(3)/2.*a,0])
        V[:,1] = np.array([1./2.*a,np.sqrt(3)/2.*a,0])
        V[:,2] = np.array([0,0,c])
        
        

        
    return V[:,0], V[:,1], V[:,2]


def reciprocal_basis(a1,a2,a3):
    

    """
    Calculate reciprocal lattice basis vectors from real space lattice vectors.
    
    The reciprocal lattice vectors satisfy: a_i · b_j = 2π δ_ij
    (Note: This implementation uses the crystallographic convention without 2π factor)
    
    Input:
        a1: numpy array (3,) - First real space lattice vector
        a2: numpy array (3,) - Second real space lattice vector
        a3: numpy array (3,) - Third real space lattice vector
    
    Output:
        b1: numpy array (3,) - First reciprocal lattice vector
        b2: numpy array (3,) - Second reciprocal lattice vector
        b3: numpy array (3,) - Third reciprocal lattice vector
    
    Usage Example:
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
    """
    b1 = np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
    b2 = np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
    b3 = np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
    return b1,b2,b3


def permut_tensor3():

    """
    Compute 3D Levi-Civita permutation symbol ε_{ijk}.
    
    Returns the permutation tensor (Levi-Civita symbol) for three indices.
    Used in cross products, curl operations, and tensor calculations.
    
    Input:
        i, j, k (int): Indices (0, 1, or 2 representing x, y, z)
    
    Output:
        int: +1 for even permutation, -1 for odd, 0 if any indices repeat
    
    Usage Example:
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
    
    Notes:
        - ε_{ijk} = +1 if (i,j,k) is even permutation of (0,1,2)
        - ε_{ijk} = -1 if (i,j,k) is odd permutation of (0,1,2)
        - ε_{ijk} = 0 if any two indices are equal
        - Used in vector cross product: (a×b)_i = ε_{ijk} a_j b_k
    
    Formula:
        ε_{012} = ε_{120} = ε_{201} = +1
        ε_{021} = ε_{210} = ε_{102} = -1
        ε_{ijk} = 0 otherwise
    """
    epsilon=[];
    for i in range(0,3):
        pj=[]
        for j in range(0,3):
            pk=[]
            for k in range(0,3):
                if i==j or j==k or k==i:
                    pk.append(0)
                elif (i,j,k)==(0,1,2) or (i,j,k)==(1,2,0) or (i,j,k)==(2,0,1):
                    pk.append(1)
                else:
                    pk.append(-1)
            pj.append(pk)
        epsilon.append(pj)
    return epsilon


def np_permut_tensor3():

    """
    Compute Levi-Civita symbol (NumPy-compatible version).
    
    Same as permut_tensor3 but optimized for NumPy array indexing.
    
    Input:
        i, j, k (int): Tensor indices
    
    Output:
        int: Permutation symbol value
    
    Usage Example:
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
    
    Notes:
        - Identical to permut_tensor3
        - NumPy-friendly naming convention
    """
    epsilon=np.empty((3,3,3));
    for i in range(0,3):
        for j in range(0,3):
            for k in range(0,3):
                if i==j or j==k or k==i:
                    val=0
                elif (i,j,k)==(0,1,2) or (i,j,k)==(1,2,0) or (i,j,k)==(2,0,1):
                    val=1
                else:
                    val=-1
                epsilon[i,j,k]=val
    return epsilon


def kronecker():

    """
    Compute Kronecker delta δ_{ij}.
    
    Returns 1 if indices are equal, 0 otherwise. Fundamental in
    tensor algebra and represents the identity tensor.
    
    Input:
        i, j (int): Indices
    
    Output:
        int: 1 if i==j, 0 otherwise
    
    Usage Example:
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
    
    Notes:
        - δ_{ij} = 1 if i = j
        - δ_{ij} = 0 if i ≠ j
        - Contraction: A_{ij} δ_{jk} = A_{ik}
        - Trace: δ_{ii} = 3 (in 3D)
    
    Formula:
        δ_{ij} = { 1  if i = j
                 { 0  if i ≠ j
    """
    delta=[];
    for i in range(0,3):
        pj=[]
        for j in range(0,3):
            if i==j:
                pj.append(1.)
            else:
                pj.append(0)
        delta.append(pj)
    return delta


def np_kronecker():

    """
    Compute Kronecker delta (NumPy version).
    
    NumPy-compatible version of kronecker().
    
    Input:
        i, j (int): Indices
    
    Output:
        int: Kronecker delta value
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Create identity matrix
        >>> n = 3
        >>> I = np.fromfunction(lambda i,j: np_kronecker(int(i), int(j)), (n,n))
        >>> print(I)
    
    Notes:
        - Same as kronecker()
        - Compatible with NumPy vectorization
    """
    delta=np.empty((3,3));
    for i in range(0,3):
        for j in range(0,3):
            if i==j:
                val=1
            else:
                val=0
            delta[i,j]=val
    return delta




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
    
  
    

# def wulffnet(ax=None,basedirs=False,baseplanes=False,basedirscol='k',baseplanescol='tab:gray', description=True,R2Proj=np.eye(3),eps=1e-5,eps2=-1e-5,rotategrid=True):
#     if ax==None:
#         fig, ax = plt.subplots()

#         fig.patch.set_alpha(0)
#     else:
#         fig=[]
#     if basedirs:
#         basicdirections = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1]]);
#         #basicdirections = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];
#         basicdirectionstext = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1]]);
#         #basicdirectionstext = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];
#         basicdirections = np.array([[1,0,0],[0,1,0],[0,0,1],\
#                                         [-1,0,0],[0,-1,0],\
#                                         [1,1,0],[1,-1,0],[-1,1,0],[-1,-1,0],\
#                                         [0,1,1],[0,-1,1],[0,1,-1],[0,-1,-1],\
#                                         [1,0,1],[-1,0,1],[1,0,-1],[-1,0,-1],\
#                                         [-1,1,1],[1,-1,1],[1,1,-1],\
#                                         [1,1,1],[-1,-1,1],[1,-1,-1],[-1,1,-1],[-1,-1,-1]]);
#         basicdirectionstext = np.array([[1,0,0],[0,1,0],[0,0,1],\
#                                         [-1,0,0],[0,-1,0],\
#                                         [1,1,0],[1,-1,0],[-1,1,0],[-1,-1,0],\
#                                         [0,1,1],[0,-1,1],[0,1,-1],[0,-1,-1],\
#                                         [1,0,1],[-1,0,1],[1,0,-1],[-1,0,-1],\
#                                         [-1,1,1],[1,-1,1],[1,1,-1],\
#                                         [1,1,1],[-1,-1,1],[1,-1,-1],[-1,1,-1],[-1,-1,-1]]);

#         basicplanes = np.array([[1,0,0],[0,1,0],[0,0,1],\
#                                         [-1,0,0],[0,-1,0],\
#                                         [1,1,0],[1,-1,0],[-1,1,0],[-1,-1,0],\
#                                         [0,1,1],[0,-1,1],\
#                                         [1,0,1],[-1,0,1]]);

#         # basicdirections = np.array([[1,0,0],[0,1,0],[0,0,1],\
#         #                                 [-1,0,0],[0,-1,0],\
#         #                                 [1,1,0],[1,-1,0],[-1,1,0],[-1,-1,0],\
#         #                                 [0,1,1],[0,-1,1],\
#         #                                 [1,0,1],[-1,0,1],\
#         #                                 [-1,1,1],[1,-1,1],\
#         #                                 [1,1,1],[-1,-1,1]]);
#         # basicdirectionstext = np.array([[1,0,0],[0,1,0],[0,0,1],\
#         #                                 [-1,0,0],[0,-1,0],\
#         #                                 [1,1,0],[1,-1,0],[-1,1,0],[-1,-1,0],\
#         #                                 [0,1,1],[0,-1,1],\
#         #                                 [1,0,1],[-1,0,1],\
#         #                                 [-1,1,1],[1,-1,1],\
#         #                                 [1,1,1],[-1,-1,1]]);

#     #longitude lines
#     #fig, ax = plt.subplots()
#     ax.tick_params(
#         axis='both',
#         which='both',
#         bottom=False,
#         top=False,
#         left=False,
#         labelbottom=False,
#         labelleft=False)
#     ax.plot(0, 0, 'k+')
#     circ = plt.Circle((0, 0), 1.0, facecolor=(210./255.,235./255.,255./255.), edgecolor='black')
#     ax.add_patch(circ)

#     ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
#   # equal aspect ratio
#     ax.axis('off')  # remove the box
#     #plt.show()
    
#     t=np.linspace(0,180,180*2+1)*np.pi/180;
#     xc = np.sin(t);
#     yc = np.cos(t);
#     AltitudeAngle = np.linspace(0,180,37)*np.pi/180;
#     AltitudeAngle = np.linspace(0,360,73)*np.pi/180;
#     for an in AltitudeAngle:
#         RotY = passive_rotation(an, 'y')
#         if rotategrid:
#             Ccp = R2Proj.dot(np.matmul(RotY,np.vstack((xc,yc,np.zeros(yc.shape)))))
#         else:
#             Ccp = np.matmul(RotY,np.vstack((xc,yc,np.zeros(yc.shape))))
#         Upper=np.where(Ccp[2,:]>eps)[0]
#         proj_dirs = stereoprojection_directions(Ccp)
#         idxs=np.where((np.diff(Upper)>1)==True)[0]
#         if idxs.shape[0]==0:
#             ax.plot(proj_dirs[0,Upper],proj_dirs[1,Upper],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
#         else:
#             idxini=0
#             for idx in idxs:
#                 ax.plot(proj_dirs[0,Upper[idxini:idx+1]],proj_dirs[1,Upper[idxini:idx+1]],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
#                 idxini=idx+1
#             ax.plot(proj_dirs[0,Upper[idxini:]],proj_dirs[1,Upper[idxini:]],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')

#         #ax.plot(proj_dirs[0,Upper],proj_dirs[1,Upper],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')

#     #Latitude Lines
#     LatitudeAngle = np.linspace(-90,90,37)*np.pi/180; 
#     LatitudeAngle = np.linspace(-180,180,73)*np.pi/180; 
#     #LatitudeAngle = np.linspace(125,135,2)*np.pi/180; 
    
#     #t=np.linspace(0,180,360*2+1)*np.pi/180;
#     zc = np.sin(t);
#     xc = np.cos(t);
#     for an in LatitudeAngle:#[0.]:#LatitudeAngle:
#         Rmeridian = np.cos(an);
#         px = Rmeridian*xc;
#         py = np.sin(an)*np.ones(t.shape);
#         pz = Rmeridian*zc;
#         if rotategrid:
#             Ccp = R2Proj.dot(np.vstack((px,py,pz)))
#         else:
#             Ccp = np.vstack((px,py,pz))
#         Upper=np.where(Ccp[2,:]>eps)[0]
#         proj_dirs = stereoprojection_directions(Ccp)
#         #print(Upper)
#         idxs=np.where((np.diff(Upper)>1)==True)[0]
#         if idxs.shape[0]==0:
#             ax.plot(proj_dirs[0,Upper],proj_dirs[1,Upper],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
#         else:
#             idxini=0
#             for idx in idxs:
#                 ax.plot(proj_dirs[0,Upper[idxini:idx+1]],proj_dirs[1,Upper[idxini:idx+1]],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
#                 idxini=idx+1
#             ax.plot(proj_dirs[0,Upper[idxini:]],proj_dirs[1,Upper[idxini:]],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
        
        
    
#     if basedirs:
#         an=45.*np.pi/180.;
#         an=0.*np.pi/180.;
#         Rotz = active_rotation(an, 'z')
#         an=np.arccos(1/np.sqrt(3));
#         an=0.*np.pi/180.;
#         Rotx = active_rotation(an, 'x')
#         dirs = R2Proj.dot(np.matmul(np.matmul(Rotx,Rotz),np.transpose(basicdirections)))
#         Upper=np.where(dirs[2,:]>eps2)[0]
#         proj_dirs = stereoprojection_directions(dirs)
#         ax.plot(proj_dirs[0,Upper],proj_dirs[1,Upper],color=basedirscol,marker='o',linestyle='',zorder=10010)
#         if description:
#             for diri,proj_diri in zip(basicdirectionstext.T[:,Upper].T,np.transpose(proj_dirs[:,Upper])):
#                 ax.text(0.03+proj_diri[0],0.03+proj_diri[1],str(diri),zorder=10010)

#     if baseplanes:
#         an=45.*np.pi/180.;
#         an=0.*np.pi/180.;
#         Rotz = active_rotation(an, 'z')
#         an=np.arccos(1/np.sqrt(3));
#         an=0.*np.pi/180.;
#         Rotx = active_rotation(an, 'x')
#         planes = R2Proj.dot(np.matmul(np.matmul(Rotx,Rotz),np.transpose(basicplanes)))
#         proj_planes, points=stereoprojection_planes(planes)
#         for proj_plane,point in zip(proj_planes,points):
#             Upper=np.where(point[2,:]>eps)[0]
#             ax.plot(proj_plane[0,Upper],proj_plane[1,Upper],color=baseplanescol,linestyle='-',zorder=10000)

#     return fig,ax

# def schmidtnet(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.)):
#     if ax==None:
#         fig, ax = plt.subplots()

#         fig.patch.set_alpha(0)
#     else:
#         fig=[]
#     if basedirs:
#         basicdirections = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1]]);
#         #basicdirections = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];
#         basicdirectionstext = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1]]);
#         #basicdirectionstext = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];



#     #longitude lines
#     #fig, ax = plt.subplots()
#     ax.tick_params(
#         axis='both',
#         which='both',
#         bottom=False,
#         top=False,
#         left=False,
#         labelbottom=False,
#         labelleft=False)
#     ax.plot(0, 0, 'k+')
#     equaarea_factor = 2./np.sqrt(2)
#     circ = plt.Circle((0, 0), equaarea_factor*1.0, facecolor=facecolor, edgecolor='black')
#     ax.add_patch(circ)

#     ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
#   # equal aspect ratio
#     ax.axis('off')  # remove the box
#     #plt.show()
    
#     t=np.linspace(0,180,180*2+1)*np.pi/180;
#     xc = np.sin(t);
#     yc = np.cos(t);
#     AltitudeAngle = np.linspace(0,180,37)*np.pi/180;
#     for an in AltitudeAngle:
#         RotY = passive_rotation(an, 'y')
#         Ccp = np.matmul(RotY,np.vstack((xc,yc,np.zeros(yc.shape))))

#         proj_dirs = equalarea_directions(Ccp)
#         ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')

#     #Latitude Lines
#     LatitudeAngle = np.linspace(-90,90,37)*np.pi/180; 
#     #t=np.linspace(0,180,360*2+1)*np.pi/180;
#     zc = np.sin(t);
#     xc = np.cos(t);
#     for an in LatitudeAngle:#[0.]:#LatitudeAngle:
#         Rmeridian = np.cos(an);
#         px = Rmeridian*xc;
#         py = np.sin(an)*np.ones(t.shape);
#         pz = Rmeridian*zc;

#         proj_dirs = equalarea_directions(np.vstack((px,py,pz)))

#         ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
    
#     if basedirs:
#         an=45.*np.pi/180.;
#         an=0.*np.pi/180.;
#         Rotz = active_rotation(an, 'z')
#         an=np.arccos(1/np.sqrt(3));
#         an=0.*np.pi/180.;
#         Rotx = active_rotation(an, 'x')
#         dirs = np.matmul(np.matmul(Rotx,Rotz),np.transpose(basicdirections))
#         proj_dirs = equalarea_directions(dirs)
#         ax.plot(proj_dirs[0,:],proj_dirs[1,:],color='b',marker='o',linestyle='')
#         for diri,proj_diri in zip(basicdirectionstext,np.transpose(proj_dirs)):
#             ax.text(0.03+proj_diri[0],0.03+proj_diri[1],str(diri))

#     return fig,ax

# def stereotriangle02(ax=None,basedirs=False):
#     if ax==None:
#         fig, ax = plt.subplots()
#     else:
#         fig=[]
#     ax.tick_params(
#         axis='both',
#         which='both',
#         bottom=False,
#         top=False,
#         left=False,
#         labelbottom=False,
#         labelleft=False)
#     ax.plot(0, 0, 'k+')
#     #circ = plt.Circle((0, 0), 1.0, facecolor=(210./255.,235./255.,255./255.), edgecolor='black')
#     #ax.add_patch(circ)

#     ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
#   # equal aspect ratio
#     ax.axis('off')  # remove the box
#     #plt.show()
    
#     normals = np.array([1,-1,0]);
#     arclength = 90#-np.arccos(np.sqrt(2)/np.sqrt(3))*180/np.pi;
#     proj_normals, points = stereoprojection_planes(normals,arclength=arclength)
#     ax.plot(proj_normals[0,:], proj_normals[1,:], 'k')
    
#     t=np.linspace(0,45,180*2)*np.pi/180;
#     xc = np.cos(t);
#     yc = np.sin(t);
#     ax.plot(xc,yc, 'k')
#     ax.plot([0,1],[0,0], 'k')
    
#     dirs = np.column_stack([[0,0,1],[1,1,0],[1,0,0]]);
# #    dirs = np.column_stack([[0,0,1],[1,1,1],[0,1,1]]);
#     proj_dirs = stereoprojection_directions(dirs)
#     ax.plot(proj_dirs[0,:], proj_dirs[1,:], 'ro')
#     if basedirs:
#         for diri,proj_diri in zip(np.transpose(dirs),np.transpose(proj_dirs)):
#             if diri==[0,0,1]:
#                 ax.text(-0.03+proj_diri[0],0.03+proj_diri[1],str(diri))
#             else:
#                 ax.text(0.03+proj_diri[0],0.03+proj_diri[1],str(diri))
    
#     return fig,ax

# def stereotriangle_ini(ax=None,basedirs=False):
#     #print('ok')
#     if ax==None:
#         fig, ax = plt.subplots()
#     else:
#         fig=[]
#     ax.tick_params(
#         axis='both',
#         which='both',
#         bottom=False,
#         top=False,
#         left=False,
#         labelbottom=False,
#         labelleft=False)
#     ax.plot(0, 0, 'k+')
#     #circ = plt.Circle((0, 0), 1.0, facecolor=(210./255.,235./255.,255./255.), edgecolor='black')
#     #ax.add_patch(circ)

#     ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
#   # equal aspect ratio
#     ax.axis('off')  # remove the box
#     #plt.show()





    
#     normals = np.array([1,-1,0]);
# #    normals = np.array([1,1,1]);
#     proj_111=stereoprojection_directions(np.array([1,1,1]))
#     arclength = 55#-np.arccos(np.sqrt(2)/np.sqrt(3))*180/np.pi;
#     proj_normals, points = stereoprojection_planes(normals,arclength=arclength,iniangle=35)
#     ax.plot(proj_normals[0,:], proj_normals[1,:], 'k')

#     normals = np.array([-1,0,1]);
#     arclength = 35#-np.arccos(np.sqrt(2)/np.sqrt(3))*180/np.pi;
#     proj_normals, points = stereoprojection_planes(normals,arclength=arclength,iniangle=90)
#     ax.plot(proj_normals[0,:], proj_normals[1,:], 'k')
    
#     R = stereoprojection_directions(np.array([1,0,1]))[0]
#     ax.plot([0,R],[0,0], 'k')
# #    dirs = np.column_stack([[0,0,1],[1,1,0],[1,0,0]]);
#     dirs = np.column_stack([[0,0,1],[1,1,1],[1,0,1]]);
#     proj_dirs = stereoprojection_directions(dirs)
#     ax.plot(proj_dirs[0,:], proj_dirs[1,:], 'ro')
#     if basedirs:
#         for diri,proj_diri in zip(np.transpose(dirs),np.transpose(proj_dirs)):
#             ax.text(0.01+proj_diri[0],0.01+proj_diri[1],str(diri))
    
#     return fig,ax

# def stereotriangle(ax=None,basedirs=False,equalarea=False):

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
#     if ax==None:
#         fig, ax = plt.subplots()
#     else:
#         fig=[]
#     ax.tick_params(
#         axis='both',
#         which='both',
#         bottom=False,
#         top=False,
#         left=False,
#         labelbottom=False,
#         labelleft=False)
#     ax.plot(0, 0, 'k+')
#     #circ = plt.Circle((0, 0), 1.0, facecolor=(210./255.,235./255.,255./255.), edgecolor='black')
#     #ax.add_patch(circ)

#     ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
  # equal aspect ratio
#     ax.axis('off')  # remove the box
#     #plt.show()

#     # t=np.linspace(0,180,180*2+1)*np.pi/180;
#     # xc = np.sin(t);
#     # yc = np.cos(t);
#     # AltitudeAngle = np.linspace(0,180,37)*np.pi/180;
#     # for an in AltitudeAngle:
#     #     RotY = passive_rotation(an, 'y')
#     #     Ccp = np.matmul(RotY,np.vstack((xc,yc,np.zeros(yc.shape))))

#     #     proj_dirs = stereoprojection_directions(Ccp)
#     #     ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')




    
#     normals = np.array([1,-1,0]);
# #    normals = np.array([1,1,1]);
#     if equalarea:
#         proj_111=equalarea_directions(np.array([1,1,1]))
#     else:
#         proj_111=stereoprojection_directions(np.array([1,1,1]))
#     arclength = 55#-np.arccos(np.sqrt(2)/np.sqrt(3))*180/np.pi;
    
#     if equalarea:
#         proj_normals = equalarea_planes(normals,arclength=arclength,iniangle=35)
#     else:
#         proj_normals,points = stereoprojection_planes(normals,arclength=arclength,iniangle=35)
        
#     ax.plot(proj_normals[0,:], proj_normals[1,:], 'k')

#     normals = np.array([-1,0,1]);
#     arclength = 35#-np.arccos(np.sqrt(2)/np.sqrt(3))*180/np.pi;
#     if equalarea:
#         proj_normals = equalarea_planes(normals,arclength=arclength,iniangle=90)
#     else:
#         proj_normals,points = stereoprojection_planes(normals,arclength=arclength,iniangle=90)
#     ax.plot(proj_normals[0,:], proj_normals[1,:], 'k')
    
#     if equalarea:
#         R = equalarea_directions(np.array([1,0,1]))[0]
#     else:
#         R = stereoprojection_directions(np.array([1,0,1]))[0]
#     ax.plot([0,R],[0,0], 'k')
# #    dirs = np.column_stack([[0,0,1],[1,1,0],[1,0,0]]);

#     dirs = np.column_stack([[0,0,1],[1,1,1],[1,0,1]]);
    
    
#     if equalarea:
#         proj_dirs = equalarea_directions(dirs)
#     else:
#         proj_dirs = stereoprojection_directions(dirs)
#     ax.plot(proj_dirs[0,:], proj_dirs[1,:], 'ro')
#     if basedirs:
#         for diri,proj_diri in zip(np.transpose(dirs),np.transpose(proj_dirs)):
#             if diri==[0,0,1]:
#                 ax.text(-0.05+proj_diri[0],0.01+proj_diri[1],str(diri))
#             else:
#                 ax.text(0.01+proj_diri[0],0.01+proj_diri[1],str(diri))
#     # resolution=1
#     # grid_cub = get_beam_directions_grid("cubic", resolution, mesh="spherified_cube_edge")
#     # grid_stereo = Rotation.from_euler(np.deg2rad(grid_cub))*Vector3d.zvector()

#     # proj_Ds = stereoprojection_directions(grid_stereo.data.T)
#     # #Colors=stereotriangle_colors(proj_Ds)
#     # #fig,ax=stereotriangle(ax=None,basedirs=basedirs)
#     # ax.scatter(proj_Ds[0,:],proj_Ds[1,:],c='k',s=5)#,'.',color='r',markersize=1)
#     #plt.show()
    
#     return fig,ax

# def wulffnet_half(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.)):
#     if ax==None:
#         fig, ax = plt.subplots()
#     else:
#         fig=[]
#     if basedirs:
#         basicdirections = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1],[1,0,2]]);
#         #basicdirections = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];
#         basicdirectionstext = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1],[1,0,2]]);
#         #basicdirectionstext = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];



#     #longitude lines
#     #fig, ax = plt.subplots()
#     ax.tick_params(
#         axis='both',
#         which='both',
#         bottom=False,
#         top=False,
#         left=False,
#         labelbottom=False,
#         labelleft=False)
#     ax.plot(0, 0, 'k+')
    
    
    
    
#     w1 = Wedge((0,0), 1.0, 0, 180, fc=facecolor, edgecolor='black')
#     ax.add_artist(w1)
    
            
    
    
# #    circ = plt.Circle((0, 0), 1.0, facecolor=(210./255.,235./255.,255./255.), edgecolor='black')
# #    ax.add_patch(circ)

#     ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
#   # equal aspect ratio
#     ax.axis('off')  # remove the box
#     #plt.show()
    
#     t=np.linspace(0,90,180*2+1)*np.pi/180;
#     xc = np.sin(t);
#     yc = np.cos(t);
#     AltitudeAngle = np.linspace(0,180,37)*np.pi/180;
#     for an in AltitudeAngle:
#         RotY = passive_rotation(an, 'y')
#         Ccp = np.matmul(RotY,np.vstack((xc,yc,np.zeros(yc.shape))))

#         proj_dirs = stereoprojection_directions(Ccp)
#         ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')

#     #Latitude Lines
#     LatitudeAngle = np.linspace(0,90,37)*np.pi/90; 
#     #t=np.linspace(0,180,360*2+1)*np.pi/180;
#     zc = np.sin(t);
#     xc = np.cos(t);
#     for an in LatitudeAngle:#[0.]:#LatitudeAngle:
#         Rmeridian = np.cos(an);
#         px = Rmeridian*xc;
#         py = np.sin(an)*np.ones(t.shape);
#         pz = Rmeridian*zc;

#         proj_dirs = stereoprojection_directions(np.vstack((px,py,pz)))

#         ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
    
#     if basedirs:
#         an=45.*np.pi/180.;
#         an=0.*np.pi/180.;
#         Rotz = active_rotation(an, 'z')
#         an=np.arccos(1/np.sqrt(3));
#         an=0.*np.pi/180.;
#         Rotx = active_rotation(an, 'x')
#         dirs = np.matmul(np.matmul(Rotx,Rotz),np.transpose(basicdirections))
#         proj_dirs = stereoprojection_directions(dirs)
#         ax.plot(proj_dirs[0,:],proj_dirs[1,:],color='b',marker='o',linestyle='')
#         for diri,proj_diri in zip(basicdirectionstext,np.transpose(proj_dirs)):
#             ax.text(0.03+proj_diri[0],0.03+proj_diri[1],str(diri))

#     return fig,ax
# def schmidtnet_half(ax=None,basedirs=False,facecolor=(210./255.,235./255.,255./255.)):
#     if ax==None:
#         fig, ax = plt.subplots()
#     else:
#         fig=[]
#     if basedirs:
#         basicdirections = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1],[1,0,2]]);
#         #basicdirections = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];
#         basicdirectionstext = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[0,1,1],[1,0,1],[1,0,2]]);
#         #basicdirectionstext = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1;0,1,1;1,0,1;1,1,-2;-1,-1,2;1,-1,0;-1,1,0];



#     #longitude lines
#     #fig, ax = plt.subplots()
#     ax.tick_params(
#         axis='both',
#         which='both',
#         bottom=False,
#         top=False,
#         left=False,
#         labelbottom=False,
#         labelleft=False)
#     ax.plot(0, 0, 'k+')
#     equaarea_factor = 2./np.sqrt(2)

#     w1 = Wedge((0,0), equaarea_factor*1.0, 0, 180, fc=facecolor, edgecolor='black')
#     ax.add_artist(w1)

#     #circ = plt.Circle((0, 0), equaarea_factor*1.0, facecolor=facecolor, edgecolor='black')
#     #ax.add_patch(circ)

#     ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
#   # equal aspect ratio
#     ax.axis('off')  # remove the box
#     #plt.show()
    
#     t=np.linspace(0,90,180*2+1)*np.pi/180;
#     xc = np.sin(t);
#     yc = np.cos(t);
#     AltitudeAngle = np.linspace(0,180,37)*np.pi/180;
#     for an in AltitudeAngle:
#         RotY = passive_rotation(an, 'y')
#         Ccp = np.matmul(RotY,np.vstack((xc,yc,np.zeros(yc.shape))))

#         proj_dirs = equalarea_directions(Ccp)
#         ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')

#     #Latitude Lines
#     LatitudeAngle = np.linspace(0,180,37)*np.pi/180; 
#     #t=np.linspace(0,180,360*2+1)*np.pi/180;
#     zc = np.sin(t);
#     xc = np.cos(t);
#     for an in LatitudeAngle:#[0.]:#LatitudeAngle:
#         Rmeridian = np.cos(an);
#         px = Rmeridian*xc;
#         py = np.sin(an)*np.ones(t.shape);
#         pz = Rmeridian*zc;

#         proj_dirs = equalarea_directions(np.vstack((px,py,pz)))

#         ax.plot(proj_dirs[0,:],proj_dirs[1,:],color=(0.5,0.5,0.5),linewidth=0.5,linestyle='--')
    
#     if basedirs:
#         an=45.*np.pi/180.;
#         an=0.*np.pi/180.;
#         Rotz = active_rotation(an, 'z')
#         an=np.arccos(1/np.sqrt(3));
#         an=0.*np.pi/180.;
#         Rotx = active_rotation(an, 'x')
#         dirs = np.matmul(np.matmul(Rotx,Rotz),np.transpose(basicdirections))
#         proj_dirs = equalarea_directions(dirs)
#         ax.plot(proj_dirs[0,:],proj_dirs[1,:],color='b',marker='o',linestyle='')
#         for diri,proj_diri in zip(basicdirectionstext,np.transpose(proj_dirs)):
#             ax.text(0.03+proj_diri[0],0.03+proj_diri[1],str(diri))

#     return fig,ax

# def wulffnet_regular_grid_ini(ax,dangle,markersize=1):
#     #dphi=10.deg
#     #dtheta=10.deg
#     #dangle = 10.
    
#     Phi1=np.linspace(0.,360.-dangle,int(360./dangle))
#     Phi2=np.linspace(0.,180.-dangle,int(180./dangle))
#     GridX=[];
#     GridY=[];
#     Dc=[1.,0.,0.];
#     for phi1 in Phi1:
#         for phi2 in Phi2:        
#             RotZ = active_rotation(phi1, 'z', deg=True) 
#             RotY = active_rotation(phi2, 'y', deg=True)
#             Ds = np.matmul(RotY,RotZ).dot(Dc)
#             proj_Ds = stereoprojection_directions(Ds)
#             GridX.append(proj_Ds[0,0])
#             GridY.append(proj_Ds[1,0])
                
    

#     #fig,ax = wulffnet()
#     ax.plot(GridX,GridY,'.',color='r',markersize=markersize)
    
#     return GridX,GridY
# def wulffnet_regular_grid(ax,dangle,dirout=False, plot=True):

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
#     #dphi=10.deg
#     #dtheta=10.deg
#     #dangle = 10.
    
#     Phi1=np.linspace(0.,360.-dangle,int(360./dangle))
#     Phi2=np.linspace(0.,180.-dangle,int(180./dangle))
#     GridX=[];
#     GridY=[];
#     Dc=[1.,0.,0.];
#     dirs=[]
#     for phi1 in Phi1:
#         for phi2 in Phi2:        
#             RotZ = active_rotation(phi1, 'z', deg=True) 
#             RotY = active_rotation(phi2, 'y', deg=True)
#             Ds = np.matmul(RotY,RotZ).dot(Dc)
#             dirs.append(Ds)
#             proj_Ds = stereoprojection_directions(Ds)
#             GridX.append(proj_Ds[0,0])
#             GridY.append(proj_Ds[1,0])
                
    

#     #fig,ax = wulffnet()
#     if plot:
#         ax.plot(GridX,GridY,'.',color='r',markersize=1)
#     if dirout:
#         return GridX,GridY,dirs
#     else:
#         return GridX,GridY

# def schmidt_regular_grid_ini(ax,dangle):
#     #dphi=10.deg
#     #dtheta=10.deg
#     #dangle = 10.
    
#     Phi1=np.linspace(0.,360.-dangle,int(360./dangle))
#     Phi2=np.linspace(0.,180.-dangle,int(180./dangle))
#     GridX=[0.];
#     GridY=[0.];
#     Dc=[1.,0.,0.];
#     for phi1 in Phi1:
#         for phi2 in Phi2:        
#             RotZ = active_rotation(phi1, 'z', deg=True) 
#             RotY = active_rotation(phi2, 'y', deg=True)
#             Ds = np.matmul(RotY,RotZ).dot(Dc)
#             proj_Ds = equalarea_directions(Ds)
#             GridX.append(proj_Ds[0,0])
#             GridY.append(proj_Ds[1,0])
                
    

#     #fig,ax = wulffnet()
#     ax.plot(GridX,GridY,'.',color='r',markersize=1)
    
#     return GridX,GridY

# def schmidt_regular_grid(ax,Na=72,Nr=20,plot=True):

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
#     dphi1=360/Na
#     phi1=np.linspace(0,360-dphi1,Na)
#     R=equalarea_directions(np.array([1,0,0]))[0,0]
#     TotalArea=np.pi*R**2
#     dr=R/(Nr+0.5)
#     r=np.linspace(0,R-dr/2,Nr+1)
# #    GridX=[0.];
# #    GridY=[0.];
# #    Weight= np.pi*(r[1]/2)**2/TotalArea
#     AreaRatio=[]
#     for ri in r:
#         Nari=int(Na*(ri/r[-1]))
        
#         if Nari<8:
#             Nari=8
#         else:
#             Nari=Nari-Nari%8    
#         #print(Nari)
#         #Nari=8
#         dphi1=360./(Nari)
#         phi1=np.linspace(0,360-dphi1,Nari)
#         phi1=np.linspace(dphi1/2,360-dphi1/2,Nari)
#         #phi1=phi1[0:-1]
#         #print(2*np.pi*((ri+dr/2)**2/2-(ri-dr/2)**2/2))
#         #tot=0
#         for phi1i in phi1:        
#             if ri==0:
#                 GridX=[0.];
#                 GridPhi=[0.]
#                 GridY=[0.];
#                 GridR=[r[1]-dr/2]
#                 GridR=[ri]
#                 AreaRatio= [np.pi*(r[1]-dr/2.)**2/TotalArea]
#             else:
#                 GridX.append(ri*np.cos(phi1i*np.pi/180.))
#                 GridY.append(ri*np.sin(phi1i*np.pi/180.))
#                 AreaRatio.append(dphi1*np.pi/180.*((ri+dr/2)**2/2.-(ri-dr/2)**2/2.)/TotalArea)
#                 GridR.append(ri)
#                 phi=np.arctan2(GridY[-1],GridX[-1])*180./np.pi
#                 GridPhi.append(phi)
#                 #tot+=dphi1*np.pi/180.*((ri+dr/2)**2/2-(ri-dr/2)**2/2)
#         #print(tot)
#     #sum(AreaRatio)
#     if plot:
#         ax.plot(GridX,GridY,'.',color='r',markersize=1)
    
#     return GridX,GridY,GridR,GridPhi,AreaRatio

# def schmidt_triangle_regular_grid(ax,Na=72,Nr=20,plot=True):
#     amax =45
#     dphi1=amax/Na
#     phi1=np.linspace(0,amax-dphi1,Na)
#     R=np.linalg.norm(equalarea_directions(np.array([1,1,1]))[:,0])
#     Rw=equalarea_directions(np.array([1,0,0]))[0,0]
#     TotalArea=np.pi*Rw**2
#     dr=R/(Nr+0.5)
#     r=np.linspace(0,R-dr/2,Nr+1)
#     normals = np.array([-1,0,1]);
#     arclength = 35#-np.arccos(np.sqrt(2)/np.sqrt(3))*180/np.pi;
#     proj_normals = equalarea_planes(normals,arclength=arclength,iniangle=90)
#     #ax.plot(proj_normals[0,:], proj_normals[1,:], 'k')
#     RN=np.sqrt(proj_normals[0,:]**2+proj_normals[1,:]**2)
#     AN=np.arctan(proj_normals[1,:]/proj_normals[0,:])*180/np.pi
# #    GridX=[0.];
# #    GridY=[0.];
# #    Weight= np.pi*(r[1]/2)**2/TotalArea
#     AreaRatio=[]
#     for ri in r:
#         Nari=int(Na*(ri/r[-1]))
        
#         if Nari<8:
#             Nari=8
#         else:
#             Nari=Nari-Nari%8    
#         #print(Nari)
#         #Nari=8
#         dphi1=amax/(Nari)
#         phi1=np.linspace(0,amax-dphi1,Nari)
#         phi1=np.linspace(dphi1/2,amax-dphi1/2,Nari)
#         #phi1=phi1[0:-1]
#         #print(2*np.pi*((ri+dr/2)**2/2-(ri-dr/2)**2/2))
#         #tot=0
#         for phi1i in phi1: 
            
#             rmax=RN[np.where(abs(AN-phi1i)==min(abs(AN-phi1i)))[0]]
            
#             #print(rmax)
#             if ri==0:
#                 GridX=[0.];
#                 GridPhi=[0.]
#                 GridY=[0.];
#                 GridR=[r[1]-dr/2]
#                 GridR=[ri]
#                 AreaRatio= [np.pi*(r[1]-dr/2.)**2/TotalArea]
#             elif ri<=rmax:
#                 GridX.append(ri*np.cos(phi1i*np.pi/180.))
#                 GridY.append(ri*np.sin(phi1i*np.pi/180.))
#                 AreaRatio.append(dphi1*np.pi/180.*((ri+dr/2)**2/2.-(ri-dr/2)**2/2.)/TotalArea)
#                 GridR.append(ri)
#                 phi=np.arctan2(GridY[-1],GridX[-1])*180./np.pi
#                 GridPhi.append(phi)
#                 #tot+=dphi1*np.pi/180.*((ri+dr/2)**2/2-(ri-dr/2)**2/2)
#         #print(tot)
#     #sum(AreaRatio)
#     if plot:
#         ax.plot(GridX,GridY,'.',color='r',markersize=1)
    
#     return GridX,GridY,GridR,GridPhi,AreaRatio

# def pf(gPhi1,gPHI,gPhi2,Dc,lattice,Na=72,Nr=20):
#     #fig,ax = schmidtnet()
#     GridX,GridY,GridR,GridPhi,AreaRatio=schmidt_regular_grid([],Na=Na,Nr=Nr,plot=False)
    
#     Intensity = np.array(GridX)*0
#     inc=0
#     #Dc=[1,0,0]
#     Syms = symmetry_elements(lattice)
    
#     for Phi1,PHI,Phi2 in zip(gPhi1,gPHI,gPhi2):
#         inc+=1
#         print(str(inc)+'/'+str(len(gPhi1)))
#         R = np.array(np_euler_matrix(Phi1, PHI,Phi2))
#         for Sym in Syms:
#             Ri=np.matmul(Sym,R)
#             Phi1,PHI,Phi2=euler_angles_from_matrix(Ri)
#             Phi1,PHI,Phi2 = euler_angles_reduction(Phi1,PHI,Phi2)
#             U = np_inverse_euler_matrix(Phi1, PHI,Phi2)
#             Ds = np.array(U).dot(Dc)
#             proj_Ds = equalarea_directions(Ds)
#             phi = np.arctan2(proj_Ds[1],proj_Ds[0])[0]*180./np.pi
#             r=np.sqrt(proj_Ds[:,0].dot(proj_Ds[:,0]))
    
#             dr=np.array(GridR)-r
#             idxmin=np.where(abs(dr)==min(abs(dr)))
#             ri=GridR[idxmin[0][0]]
    
#             idxr=np.where(np.array(GridR)==ri)[0]
#             dtan = np.array(GridPhi)[idxr]-phi
#             idxtan = np.where(abs(dtan)==min(abs(dtan)))
#             idxmin3=idxr[idxtan[0]][0]
            
#     #        dx = np.array(GridX)-proj_Ds[0];
#     #        dy = np.array(GridY)-proj_Ds[1];
#     #        dr = dx**2+dy**2
#     #        idxmin = np.where(dr==min(dr))[0][0]
#             Intensity[idxmin3]+=1/AreaRatio[idxmin3]
#     #        GridX[idxmin]
#     #        GridY[idxmin]
    
    
    
        
#     #fig, ax = plt.subplots()
#     fig,ax = schmidtnet()
#     GridX,GridY,GridR,GridPhi,AreaRatio=schmidt_regular_grid(ax,Na=Na,Nr=Nr,plot=True)

#     plt.scatter(GridX,GridY, c=Intensity, s=50, edgecolor='',zorder=10,cmap='jet')
#     cb=plt.colorbar()
#     #plt.show()

def B19p_B2_lattice_correspondence(notation='Miyazaki'):

    """
    Generate B19'→B2 lattice correspondence matrix for NiTi.
    
    Returns the transformation matrix relating B19' monoclinic martensite
    to B2 cubic austenite lattice vectors. Fundamental for NiTi shape
    memory alloy analysis.
    
    Input:
        None
    
    Output:
        numpy.ndarray (3×3): Correspondence matrix C where L_B2 = C · L_B19'
    
    Usage Example:
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
    
    Notes:
        - Specific to NiTi shape memory alloys
        - Based on crystallographic theory
        - Used in transformation analysis
        - Related to habit plane calculations
    
    Formula:
        C = [[c11, c12, c13],
             [c21, c22, c23],
             [c31, c32, c33]]
        where coefficients determined experimentally
    """
    #check
    #Variant = 8
    #B19puvw_2_B2uvw_all[:,:,Variant]
    #testB19v=[1,2,3]
    #testB2v = B19puvw_2_B2uvw_all[:,:,Variant].dot(testB2v)            
    #testB19v-np.linalg.inv(B19puvw_2_B2uvw_all[:,:,Variant]).dot(testB2v)
    #testB19p=[1,2,1]
    #testB2p = B19phkl_2_B2hkl_all[:,:,Variant].dot(testB19p)            
    #testB19p-np.linalg.inv(B19phkl_2_B2hkl_all[:,:,Variant]).dot(testB2p)
    #testB2p=[1,-3,1]

    #correspondance matrix betwen uvw of B19p and that of B2 /B19[uvw]->B2[uvw]/
    B19puvw_2_B2uvw_all = np.empty((3,3,12))
    B19phkl_2_B2hkl_all = np.empty((3,3,12))
    B2uvw_2_B19puvw_all = np.empty((3,3,12))
    B2hkl_2_B19phkl_all = np.empty((3,3,12))
    
    for Variant in range(B19puvw_2_B2uvw_all.shape[2]):
        
        if (Variant+1)==1:
            
            if notation=='Miyazaki':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 1,0,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 0,1,1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,1];
            elif notation=='Waitz':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 1,0,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 0,1,-1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,1,1];

            
            B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;          
                
    #            B19puvw_2_B2uvw_all[:,0,Variant] = [-1, 0, 0];
    #            B19puvw_2_B2uvw_all[:,1,Variant] = [ 0, 1,-1];
    #            B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,-1];
                
    #            B19puvw_2_B2uvw_all[:,0,Variant] = [1, 0, 0];
    #            B19puvw_2_B2uvw_all[:,1,Variant] = [ 0, 1,1];
    #            B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,1];
    #            
    
        if (Variant+1)==2:
            if notation=='Miyazaki':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ -1,0,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 0,-1,-1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,1];
            elif notation=='Waitz':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 1,0,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 0,-1,1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,-1];
                
            B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
                
    #            B19puvw_2_B2uvw_all[:,0,Variant] = [-1, 0, 0];
    #            B19puvw_2_B2uvw_all[:,1,Variant] = [ 0, -1,-1];
    #            B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,1];
    #
        if (Variant+1)==3:
            if notation=='Miyazaki':
                B19puvw_2_B2uvw_all[:,0,Variant] = [1, 0, 0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 0, -1,1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,-1];
            elif notation=='Waitz':
                B19puvw_2_B2uvw_all[:,0,Variant] = [1, 0, 0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 0, 1,1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,1];
    
            B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
                
    #            B19puvw_2_B2uvw_all[:,0,Variant] = [1, 0, 0];
    #            B19puvw_2_B2uvw_all[:,1,Variant] = [ 0, -1,1];
    #            B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,-1];
    
        if (Variant+1)==4:
            if notation=='Miyazaki':
                B19puvw_2_B2uvw_all[:,0,Variant] = [-1,0,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 0,1,-1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,-1];
            elif notation=='Waitz':
                B19puvw_2_B2uvw_all[:,0,Variant] = [1,0,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 0,-1,-1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,1,-1];
            
            B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==5:
            if notation=='Miyazaki':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,1,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [1,0,1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [1,0,-1];
            elif notation=='Waitz':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,-1,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [-1,0,1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [-1,0,-1];
                
            B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==6:
            if notation=='Miyazaki':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,-1,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ -1,0,-1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [1,0,-1];
            elif notation=='Waitz':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,-1,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [1,0,-1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [1,0,1];
        
            B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==7:
            if notation=='Miyazaki':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,1,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 1,0,-1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ -1,0,-1];
            elif notation=='Waitz':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,1,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 1,0,1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 1,0,-1];
        
            B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==8:
            if notation=='Miyazaki':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,-1,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [-1,0,1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ -1,0,-1];
            elif notation=='Waitz':        
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,1,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [-1,0,-1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ -1,0,1];
            B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
                
        if (Variant+1)==9:
            if notation=='Miyazaki':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,0,1];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 1,1,0];
                B19puvw_2_B2uvw_all[:,2,Variant] = [-1,1,0];
            elif notation=='Waitz':        
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,0,-1];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 1,-1,0];
                B19puvw_2_B2uvw_all[:,2,Variant] = [-1,-1,0];
        
            B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==10:
            if notation=='Miyazaki':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,0,-1];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ -1,-1,0];
                B19puvw_2_B2uvw_all[:,2,Variant] = [-1,1,0];
            elif notation=='Waitz':        
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,0,-1];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ -1,1,0];
                B19puvw_2_B2uvw_all[:,2,Variant] = [1,1,0];
        
            B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==11:
            if notation=='Miyazaki':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,0,1];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ -1,1,0];
                B19puvw_2_B2uvw_all[:,2,Variant] = [-1,-1,0];
            elif notation=='Waitz':        
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,0,1];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ -1,-1,0];
                B19puvw_2_B2uvw_all[:,2,Variant] = [1,-1,0];
    
        
            B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==12:
            if notation=='Miyazaki':
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,0,-1];
                B19puvw_2_B2uvw_all[:,1,Variant] = [1,-1,0];
                B19puvw_2_B2uvw_all[:,2,Variant] = [-1,-1,0];
            elif notation=='Waitz':        
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,0,1];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 1,1,0];
                B19puvw_2_B2uvw_all[:,2,Variant] = [-1,1,0];
          
            B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        
        B2uvw_2_B19puvw_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]);
        B2hkl_2_B19phkl_all[:,:,Variant] = inv(B19phkl_2_B2hkl_all[:,:,Variant]);         

    return B19puvw_2_B2uvw_all,B2uvw_2_B19puvw_all,B19phkl_2_B2hkl_all,B2hkl_2_B19phkl_all


def lattice_correspondence(LatCorr,parent_symops,product_symops):

    """
    Calculate lattice correspondence matrix between two crystal structures.
    
    Finds transformation matrix C relating two lattices: L2 = C · L1.
    Optionally optimizes to minimize lattice mismatch.
    
    Input:
        L1 (array 3×3): First lattice matrix
        L2 (array 3×3): Second lattice matrix
        optimize (bool): If True, optimize correspondence (default: False)
    
    Output:
        numpy.ndarray (3×3): Correspondence matrix C
    
    Usage Example:
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
    
    Notes:
        - Simple version: C = L2 · L1^(-1)
        - Optimized version minimizes ||L2 - C·L1||
        - Used in phase transformation analysis
        - Essential for orientation relationships
    
    Formula:
        C = L2 · L1^(-1)
    """
    #LatCorr=np.empty((3,3,1))
    #Variant=0
    #LatCorr[:,0,Variant] = [ 1,0,0];
    #LatCorr[:,1,Variant] = [ 0,1,-1];
    #LatCorr[:,2,Variant] = [ 0,1,1];
    count=-1
    #print(LatCorr[:,:,0])
    for syms in parent_symops:
        count+=1
        isin=False
        for var in range(LatCorr.shape[2]):#check if not yet included
            for symsB19p in product_symops:#take into account symmetry of product lattice
                if (abs(LatCorr[:,:,var]-syms.dot(LatCorr[:,:,0].dot(symsB19p)))<1e-10).all():# or (abs(-1*LatCorr[:,:,var]-syms.dot(LatCorr[:,:,0]))<1e-10).all():
                    isin=True
    
        if not isin and (syms.dot(LatCorr[:,:,0])[:,0]>=0).all() and np.linalg.det(syms.dot(LatCorr[:,:,0]))>0:#select correspondence that preserves right-handed system between [100]_a, [010]_a [001]_a
            #Variant+=1
            #LatCorr[:,:,Variant]=syms.dot(LatCorr[:,:,0])
            LatCorr=np.concatenate((LatCorr,np.expand_dims(syms.dot(LatCorr[:,:,0]),axis=(2))),axis=2)
            #print(f'count={count}')
            #print(LatCorr[:,:,-1])
            #print


def B19p_B2_lattice_correspondence_ini():

    """
    Initialize B19'→B2 correspondence for specific transformation variant.
    
    Returns correspondence matrix for one of the crystallographic variants
    of the B19' martensite transformation in NiTi.
    
    Input:
        variant (int): Variant number (1-24 for cubic→monoclinic)
    
    Output:
        numpy.ndarray (3×3): Correspondence matrix for this variant
    
    Usage Example:
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
    
    Notes:
        - NiTi has 24 crystallographic variants
        - Each variant has different correspondence
        - Related to symmetry of parent phase
        - Used in texture simulations
    """
    #check
    #Variant = 8
    #B19puvw_2_B2uvw_all[:,:,Variant]
    #testB19v=[1,2,3]
    #testB2v = B19puvw_2_B2uvw_all[:,:,Variant].dot(testB2v)            
    #testB19v-np.linalg.inv(B19puvw_2_B2uvw_all[:,:,Variant]).dot(testB2v)
    #testB19p=[1,2,1]
    #testB2p = B19phkl_2_B2hkl_all[:,:,Variant].dot(testB19p)            
    #testB19p-np.linalg.inv(B19phkl_2_B2hkl_all[:,:,Variant]).dot(testB2p)
    #testB2p=[1,-3,1]

    #correspondance matrix betwen uvw of B19p and that of B2 /B19[uvw]->B2[uvw]/
    B19puvw_2_B2uvw_all = np.empty((3,3,12))
    B19phkl_2_B2hkl_all = np.empty((3,3,12))
    B2uvw_2_B19puvw_all = np.empty((3,3,12))
    B2hkl_2_B19phkl_all = np.empty((3,3,12))
    
    for Variant in range(B19puvw_2_B2uvw_all.shape[2]):
        
        if (Variant+1)==1:
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 1,0,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 0,-1,1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,1,1];
                B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;          
                
    #            B19puvw_2_B2uvw_all[:,0,Variant] = [-1, 0, 0];
    #            B19puvw_2_B2uvw_all[:,1,Variant] = [ 0, 1,-1];
    #            B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,-1];
                
    #            B19puvw_2_B2uvw_all[:,0,Variant] = [1, 0, 0];
    #            B19puvw_2_B2uvw_all[:,1,Variant] = [ 0, 1,1];
    #            B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,1];
    #            
    
        if (Variant+1)==2:
                B19puvw_2_B2uvw_all[:,0,Variant] = [ -1,0,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 0,1,-1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,1,1];
                B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
                
    #            B19puvw_2_B2uvw_all[:,0,Variant] = [-1, 0, 0];
    #            B19puvw_2_B2uvw_all[:,1,Variant] = [ 0, -1,-1];
    #            B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,1];
    #
        if (Variant+1)==3:
                B19puvw_2_B2uvw_all[:,0,Variant] = [1, 0, 0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 0, -1,-1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,1];
    
                B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
                
    #            B19puvw_2_B2uvw_all[:,0,Variant] = [1, 0, 0];
    #            B19puvw_2_B2uvw_all[:,1,Variant] = [ 0, -1,1];
    #            B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,-1];
    
        if (Variant+1)==4:
                B19puvw_2_B2uvw_all[:,0,Variant] = [-1,0,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 0,1,1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 0,-1,1];
                B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==5:
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,1,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [1,0,-1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [1,0,1];
                B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==6:
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,-1,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ -1,0,1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [1,0,1];
                B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==7:
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,-1,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 1,0,1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 1,0,-1];
                B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==8:
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,1,0];
                B19puvw_2_B2uvw_all[:,1,Variant] = [-1,0,-1];
                B19puvw_2_B2uvw_all[:,2,Variant] = [ 1,0,-1];
                B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
                
        if (Variant+1)==9:
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,0,1];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ -1,1,0];
                B19puvw_2_B2uvw_all[:,2,Variant] = [1,1,0];
                B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==10:
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,0,-1];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 1,-1,0];
                B19puvw_2_B2uvw_all[:,2,Variant] = [1,1,0];
                B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==11:
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,0,-1];
                B19puvw_2_B2uvw_all[:,1,Variant] = [ 1,1,0];
                B19puvw_2_B2uvw_all[:,2,Variant] = [-1,1,0];
    
                B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        if (Variant+1)==12:
                B19puvw_2_B2uvw_all[:,0,Variant] = [ 0,0,1];
                B19puvw_2_B2uvw_all[:,1,Variant] = [-1,-1,0];
                B19puvw_2_B2uvw_all[:,2,Variant] = [-1,1,0];
                B19phkl_2_B2hkl_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]).T;
        
        B2uvw_2_B19puvw_all[:,:,Variant] = inv(B19puvw_2_B2uvw_all[:,:,Variant]);
        B2hkl_2_B19phkl_all[:,:,Variant] = inv(B19phkl_2_B2hkl_all[:,:,Variant]);         

    return B19puvw_2_B2uvw_all,B2uvw_2_B19puvw_all,B19phkl_2_B2hkl_all,B2hkl_2_B19phkl_all
#nn=0
#for ii in range(0,12):
#    Cd=B19puvw_2_B2uvw_all[:,:,ii];
#    for i in range(0,2):
#        for j in range(0,2):
#            if Cd[i,j]==round(Cd[i,j]):
#                Cd[i,j]=round(Cd[i,j])
#    nn+=1
#    #print("\\begin{{pmatrix}}{} & {} &{}\\\\{} &{} &{}\\\\{} &{} &{}\\end{{pmatrix}}".format(Cd[0,0],Cd[0,1],Cd[0,2],Cd[1,0],Cd[1,1],Cd[1,2],Cd[2,0],Cd[2,1],Cd[2,2]))
#    print("\\Bigl(\\begin{{smallmatrix}}{} & {} &{}\\\\{} &{} &{}\\\\{} &{} &{}\\end{{smallmatrix}}\\Bigr)".format(Cd[0,0],Cd[0,1],Cd[0,2],Cd[1,0],Cd[1,1],Cd[1,2],Cd[2,0],Cd[2,1],Cd[2,2]))
#    print('&')
#    Cd=B19phkl_2_B2hkl_all[:,:,ii]
#    #print("\\begin{{pmatrix}}{} & {} &{}\\\\{} &{} &{}\\\\{} &{} &{}\\end{{pmatrix}}".format(Cd[0,0],Cd[0,1],Cd[0,2],Cd[1,0],Cd[1,1],Cd[1,2],Cd[2,0],Cd[2,1],Cd[2,2]))
#    print("\\Bigl(\\begin{{smallmatrix}}{} & {} &{}\\\\{} &{} &{}\\\\{} &{} &{}\\end{{smallmatrix}}\\Bigr)".format(Cd[0,0],Cd[0,1],Cd[0,2],Cd[1,0],Cd[1,1],Cd[1,2],Cd[2,0],Cd[2,1],Cd[2,2]))
#    if nn<4:
#        print('&')
#    else:
#        print('\\\\')
#        print('\hline')
#        nn=0
##    print(B19puvw_2_B2uvw_all[:,:,ii])
#    print(B19phkl_2_B2hkl_all[:,:,ii])
#    print('============================================')
    
    

def cubic2tetragonal_lattice_correspondence():

    """
    Generate cubic→tetragonal lattice correspondence.
    
    Returns correspondence matrix for cubic to tetragonal phase
    transformation (e.g., FCC→FCT).
    
    Input:
        None
    
    Output:
        numpy.ndarray (3×3): Correspondence matrix
    
    Usage Example:
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
    
    Notes:
        - Common in martensitic transformations
        - Simple correspondence: diagonal matrix
        - c/a ratio determines tetragonality
    """
    #check
    #Variant = 8
    #B19puvw_2_B2uvw_all[:,:,Variant]
    #testB19v=[1,2,3]
    #testB2v = B19puvw_2_B2uvw_all[:,:,Variant].dot(testB2v)            
    #testB19v-np.linalg.inv(B19puvw_2_B2uvw_all[:,:,Variant]).dot(testB2v)
    #testB19p=[1,2,1]
    #testB2p = B19phkl_2_B2hkl_all[:,:,Variant].dot(testB19p)            
    #testB19p-np.linalg.inv(B19phkl_2_B2hkl_all[:,:,Variant]).dot(testB2p)
    #testB2p=[1,-3,1]

    #correspondance matrix betwen uvw of B19p and that of B2 /B19[uvw]->B2[uvw]/
    Productuvw_2_Parentuvw_all = np.empty((3,3,3))
    Producthkl_2_Parenthkl_all = np.empty((3,3,3))
    Parentuvw_2_Productuvw_all = np.empty((3,3,3))
    Parenthkl_2_Producthkl_all = np.empty((3,3,3))
    
    for Variant in range(Productuvw_2_Parentuvw_all.shape[2]):
        
        if (Variant+1)==1:
                Productuvw_2_Parentuvw_all[:,0,Variant] = [ 0,1,1];
                Productuvw_2_Parentuvw_all[:,1,Variant] = [ 0,-1,1];
                Productuvw_2_Parentuvw_all[:,2,Variant] = [ 1,0,0];
                Producthkl_2_Parenthkl_all[:,:,Variant] = inv(Productuvw_2_Parentuvw_all[:,:,Variant]).T;          
                
        if (Variant+1)==2:
                Productuvw_2_Parentuvw_all[:,0,Variant] = [ 1,0,-1];
                Productuvw_2_Parentuvw_all[:,1,Variant] = [ 1,0,1];
                Productuvw_2_Parentuvw_all[:,2,Variant] = [ 0,1,0];
                Producthkl_2_Parenthkl_all[:,:,Variant] = inv(Productuvw_2_Parentuvw_all[:,:,Variant]).T;          

        if (Variant+1)==3:
                Productuvw_2_Parentuvw_all[:,0,Variant] = [ 1,1,0];
                Productuvw_2_Parentuvw_all[:,1,Variant] = [ -1,1,0];
                Productuvw_2_Parentuvw_all[:,2,Variant] = [ 0,0,1];
                Producthkl_2_Parenthkl_all[:,:,Variant] = inv(Productuvw_2_Parentuvw_all[:,:,Variant]).T;          


        
        Parentuvw_2_Productuvw_all[:,:,Variant] = inv(Productuvw_2_Parentuvw_all[:,:,Variant]);
        Parenthkl_2_Producthkl_all[:,:,Variant] = inv(Producthkl_2_Parenthkl_all[:,:,Variant]);         

    return Parentuvw_2_Productuvw_all,Productuvw_2_Parentuvw_all,Parenthkl_2_Producthkl_all,Producthkl_2_Parenthkl_all

def Rp_B2_lattice_correspondence():

    """
    Generate R-phase→B2 lattice correspondence for NiTi.
    
    Returns correspondence matrix for R-phase (rhombohedral) to B2
    transformation in NiTi alloys. R-phase is intermediate phase.
    
    Input:
        None
    
    Output:
        numpy.ndarray (3×3): R-phase → B2 correspondence
    
    Usage Example:
        >>> C_R_B2 = Rp_B2_lattice_correspondence()
        >>> print("R-phase → B2 correspondence:")
        >>> print(C_R_B2)
    
    Notes:
        - R-phase is pre-martensitic phase in NiTi
        - Trigonal/rhombohedral structure
        - Appears above Ms temperature
        - Correspondence simpler than B19'→B2
    """
    #correspondance matrix betwen uvw of Rp and that of B2 /B19[uvw]->B2[uvw]/
    Rpuvw_2_B2uvw_all = np.empty((3,3,4))
    Rphkl_2_B2hkl_all = np.empty((3,3,4))
    B2uvw_2_Rpuvw_all = np.empty((3,3,4))
    B2hkl_2_Rphkl_all = np.empty((3,3,4))
    
    for Variant in range(Rpuvw_2_B2uvw_all.shape[2]):
        
        if (Variant+1)==1:
                Rpuvw_2_B2uvw_all[:,0,Variant] = [ 1,-2,1];
                Rpuvw_2_B2uvw_all[:,1,Variant] = [ 1,1,-2];
                Rpuvw_2_B2uvw_all[:,2,Variant] = [ 1,1,1];
                Rphkl_2_B2hkl_all[:,:,Variant] = inv(Rpuvw_2_B2uvw_all[:,:,Variant]).T;          
                
    #            
    
        if (Variant+1)==2:
                Rpuvw_2_B2uvw_all[:,0,Variant] = [ 2,1,1];
                Rpuvw_2_B2uvw_all[:,1,Variant] = [ -1,1,-2];
                Rpuvw_2_B2uvw_all[:,2,Variant] = [ -1,1,1];
                Rphkl_2_B2hkl_all[:,:,Variant] = inv(Rpuvw_2_B2uvw_all[:,:,Variant]).T;
                
    #
        if (Variant+1)==3:
                Rpuvw_2_B2uvw_all[:,0,Variant] = [1,2,1];
                Rpuvw_2_B2uvw_all[:,0,Variant] = [1,-2,-1];
                Rpuvw_2_B2uvw_all[:,1,Variant] = [1,1,2];
                #Rpuvw_2_B2uvw_all[:,1,Variant] = [1,-1,-2];
                Rpuvw_2_B2uvw_all[:,2,Variant] = [ -1,-1,1];
    
                Rphkl_2_B2hkl_all[:,:,Variant] = inv(Rpuvw_2_B2uvw_all[:,:,Variant]).T;
                
    
        if (Variant+1)==4:
                Rpuvw_2_B2uvw_all[:,0,Variant] = [2,1,1];
                Rpuvw_2_B2uvw_all[:,0,Variant] = [2,1,-1];
                Rpuvw_2_B2uvw_all[:,1,Variant] = [ 1,1,-2];
                Rpuvw_2_B2uvw_all[:,1,Variant] = [ -1,1,2];
                Rpuvw_2_B2uvw_all[:,2,Variant] = [ 1,-1,1];
                Rphkl_2_B2hkl_all[:,:,Variant] = inv(Rpuvw_2_B2uvw_all[:,:,Variant]).T;

        
        B2uvw_2_Rpuvw_all[:,:,Variant] = inv(Rpuvw_2_B2uvw_all[:,:,Variant]);
        B2hkl_2_Rphkl_all[:,:,Variant] = inv(Rphkl_2_B2hkl_all[:,:,Variant]);         

    return Rpuvw_2_B2uvw_all,B2uvw_2_Rpuvw_all,Rphkl_2_B2hkl_all,B2hkl_2_Rphkl_all


def print_correspondence(Mcorr,VecA,latticeA, latticeB,planes=False,returnB=False):

    """
    Print lattice correspondence matrix in readable format.
    
    Displays correspondence matrix with phase labels for documentation
    and reporting purposes.
    
    Input:
        C (array 3×3): Correspondence matrix
        phase1 (str): Name of first phase (default: 'Phase1')
        phase2 (str): Name of second phase (default: 'Phase2')
    
    Output:
        None (prints to console)
    
    Usage Example:
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
    
    Notes:
        - Formatted for readability
        - Useful in reports and documentation
        - Shows transformation relationship clearly
    """
    if type(VecA[0])==list:
        VecA=np.array(VecA)
    if not VecA.ndim==2:
        VecA=np.expand_dims(VecA, axis=0)    
    if not VecA.shape[0]==3:
        VecA=VecA.T
    #print basal directions
    if planes:
        strfunction='plane2string'
    else:
        strfunction='dir2string'
        
    BasalDirs = np.empty((3,3,Mcorr.shape[2]))
    VecB = np.empty((3,VecA.shape[1],Mcorr.shape[2]))
    print('================================================================================')
    print('Basal directions')
    print('================================================================================')
    for var in range(0,Mcorr.shape[2]):
        print("\t"+'--------------------------------------------------------------------------------')
        print("\t"+'Variant '+str(var+1))
        print("\t"+'--------------------------------------------------------------------------------')
        for ei,i in zip(np.eye(3),range(0,3)):
            strvari=eval(strfunction+'(ei)')+"_"+latticeA+"=="+eval(strfunction+'(Mcorr[:,i,var])')+"_"+latticeB;
            BasalDirs[:,i,var]=Mcorr[:,i,var]
            print("\t"+strvari)

    print('================================================================================')
    print('Input directions')
    print('================================================================================')
    #print(VecA.shape)
    for var in range(0,Mcorr.shape[2]):
        print("\t"+'--------------------------------------------------------------------------------')
        print("\t"+'Variant '+str(var+1))
        print("\t"+'--------------------------------------------------------------------------------')
        inc=-1
        for va in VecA.T:
            inc+=1
#            print(va)
#            print('================================================================================')
#            print('Direction '+dir2string(va)+"_"+latticeA)
#            print('================================================================================')
#            print(va)
            vb = Mcorr[:,:,var].dot(va)
            if not planes:
                vb=miller2fractional(vb)
            str2print=eval(strfunction+'(va)')+"_"+latticeA+"=="+eval(strfunction+'(vb)')+"_"+latticeB;
            print("\t"+str2print)
            #print(vb)
            VecB[:,inc,var]=vb
 
    if returnB:
        return VecB
            #    var2=dir2string([0,1,0])+"_B19'=="+dir2string(B19puvw_2_B2uvw_all[:,1,var])+'_B2'
#    var3=dir2string([0,0,1])+"_B19'=="+dir2string(B19puvw_2_B2uvw_all[:,2,var])+'_B2'
#    #print('Varianta '+str(var+1)+':'+var1+','+var2+','+var3)
#    print('Varianta '+str(var+1))
#    print(plane2string(b19pvar, digits=2))


def mohr_circles(tensor):

    """
    Calculate Mohr's circles from strain or stress tensor.
    
    Computes principal strains/stresses and parameters for three Mohr's
    circles. Used for visualizing 3D strain state in 2D.
    
    Input:
        strain_tensor (array 3×3): Symmetric strain or stress tensor
    
    Output:
        dict: {
            'principal': array [3] - principal values (sorted)
            'directions': array 3×3 - principal directions
            'circles': list of 3 tuples - (center, radius) for each circle
        }
    
    Usage Example:
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
    
    Notes:
        - Three circles for 3D state
        - Largest circle: (ε₁ - ε₃)/2
        - Used in failure analysis
        - Shows all possible strain states on planes
    
    Formula:
        Circle i: center = (σᵢ + σⱼ)/2, radius = |σᵢ - σⱼ|/2
        Three circles: 1-2, 2-3, 1-3 planes
    """
    DD,VV = np.linalg.eig(tensor)
    Idxs = np.argsort(DD)[::-1]      
    DD=DD[Idxs]
    #print(DD)
    VV=VV[:,Idxs]    
    #EpsT =  np.matmul(matmul(VV.T,tensor),VV)
    
    
    #plot mohr cicrles
    R13=(DD[0]-DD[2])/2
    C13=DD[0]-R13
    R12=(DD[0]-DD[1])/2
    C12=DD[0]-R12#(DD[0]+DD[1])/2
    R23=(DD[1]-DD[2])/2
    C23=DD[1]-R23

    mohr_circles={'C12':C12,'C13':C13,'C23':C23,'R12':R12,'R13':R13,'R23':R23,\
                  'D1':DD[0],'D2':DD[1],'D3':DD[2],\
                  'V1':VV[:,0],'V2':VV[:,1],'V3':VV[:,2]}
    return mohr_circles,VV,DD


def generate_lattice_points(uvw2xyz,basal_dirs):

    """
    Generate lattice points within specified unit cell range.
    
    Creates array of lattice points for visualization and analysis.
    Generates n1×n2×n3 unit cells.
    
    Input:
        L (array 3×3): Lattice matrix [a₁|a₂|a₃]
        n1, n2, n3 (int): Number of unit cells in each direction
    
    Output:
        numpy.ndarray (N×3): Array of lattice point coordinates
    
    Usage Example:
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
    
    Notes:
        - Includes origin (0,0,0)
        - Points at integer lattice coordinates
        - Used in 3D visualization
        - Can visualize multiple unit cells
    """
    Lattice=[]
    #basic lattice
    for idx in range(0,basal_dirs.shape[1]):
        v=uvw2xyz.dot(basal_dirs[:,idx])
        Lattice.append([[0,v[0]],[0,v[1]],[0,v[2]]])
        v3=np.array([0,0,0]);
        for idx2 in range(0,basal_dirs.shape[1]):
            if not idx2==idx:
                v2=uvw2xyz.dot(basal_dirs[:,idx2])
                Lattice.append([[0,v[0]]+v2[0],[0,v[1]]+v2[1],[0,v[2]]+v2[2]])
                v3=v2+v3
        Lattice.append([[0,v[0]]+v3[0],[0,v[1]]+v3[1],[0,v[2]]+v3[2]])
        Lattice2=copy.deepcopy(Lattice)
    #Translation of basic lattice
    #Lattice=[]
    for basal_dir in basal_dirs.T:
        TransDir = -1.
        Trans_vec = TransDir*uvw2xyz.dot(basal_dir)
        for point in Lattice2:
            trpoint=[];
            for p,trv in zip(point,Trans_vec):
                trpoint.append([p[0]+trv,p[1]+trv])
            Lattice.append(trpoint)
    for idxs in [[0,1],[0,2],[1,2],[0,1,2]]:
        TransDir = -1.
        Trans_vec=np.array([0.,0.,0.]);
        for idx in idxs:
            Trans_vec+=TransDir*uvw2xyz.dot(basal_dirs[:,idx])
        for point in Lattice2:
            trpoint=[];
            for p,trv in zip(point,Trans_vec):
                trpoint.append([p[0]+trv,p[1]+trv])
            Lattice.append(trpoint)
            
    #remove duplicates            
    Lattice2=[];
    eps=1.e-6
    for points in Lattice:
        if len(Lattice2)==0:
            Lattice2.append(points)
        else:
            isin=False
            for points2 in Lattice2:
                for point2,point in zip(points2,points):
                    if abs(point2[0]-point[0])<eps and abs(point2[1]-point[1])<eps:
                        isin=True
                    else:False
            if not isin:
                Lattice2.append(points)
#    Lattice=[]
#    for phix in [0,180]:    
#        Rx=passive_rotation(phix, 'x', deg=True)     
#
#        for phiz in [0,90,180,270]:    
#            Rz=passive_rotation(phiz, 'z', deg=True)     
#            for point in Lattice2:
#                p0=[]
#                p1=[]
#                for p in point:
#                    p0.append(p[0])    
#                    p1.append(p[1])    
#                p0r=Rz.dot(Rx.dot(p0))
#                p1r=Rz.dot(Rx.dot(p1))
#                B2point=[]
#                for p0ri,p1ri in zip(p0r,p1r):
#                    B2point.append([p0ri,p1ri])
#                Lattice.append(B2point)
    return Lattice

def generate_lattice_faces(uvw2xyz,basal_dirs):


    """
    Generate face polygons for unit cell visualization.
    
    Creates vertex lists for 6 faces of unit cell for 3D rendering.
    
    Input:
        L (array 3×3): Lattice matrix
    
    Output:
        list: List of 6 face vertex arrays (each 4×3)
    
    Usage Example:
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
    
    Notes:
        - Returns 6 faces (parallelepiped)
        - Each face is quadrilateral
        - Vertices in counter-clockwise order
        - Used in Poly3DCollection
    """
#uvw2xyz=Parent_uvw2xyz,basal_dirs=Parent
    LatticeFaces=[]
    #basic lattice
    for idx in range(0,basal_dirs.shape[1]):
        idxs=list(range(basal_dirs.shape[1]))#[0,1,2];
        idxs.remove(idx)
        face=[]
        face.append(np.array([0.,0.,0.]))
        for idx2 in idxs:
            face.append(face[-1]+uvw2xyz.dot(basal_dirs[:,idx2]))
        face.append(face[-1]-uvw2xyz.dot(basal_dirs[:,idxs[0]]))    
        face.append(face[0])    
        LatticeFaces.append(face)
        face2=copy.deepcopy(face)
        for i in range(0,len(face2)):
            face2[i]=face2[i]+uvw2xyz.dot(basal_dirs[:,idx])
        LatticeFaces.append(face2)
    
    LatticeFaces2=copy.deepcopy(LatticeFaces)
    #LatticeFaces=[]
    Lattices=[]
    Lattices.append(LatticeFaces)        
    #Translation of basic lattice
    #Lattice=[]
    for basal_dir in basal_dirs.T:
        LatticeFaces=[]
        TransDir = -1.
        Trans_vec = TransDir*uvw2xyz.dot(basal_dir)
        for face in LatticeFaces2:            
            trface=[];
            for edge in face:
                trface.append(edge+Trans_vec)
            LatticeFaces.append(trface)
        Lattices.append(LatticeFaces)
    for idxs in [[0,1],[0,2],[1,2],[0,1,2]]:
        LatticeFaces=[]
        TransDir = -1.
        Trans_vec=np.array([0.,0.,0.]);
        for idx in idxs:
            Trans_vec+=TransDir*uvw2xyz.dot(basal_dirs[:,idx])
        for face in LatticeFaces2:
            trface=[];
            for edge in face:
                trface.append(edge+Trans_vec)
            LatticeFaces.append(trface)
        Lattices.append(LatticeFaces)

            
    #remove duplicates            
    Lattices2=[];
    eps=1.e-6
    for lattice in Lattices:
        if len(Lattices2)==0:
            Lattices2.append(lattice)
        else:
            Lat=[]
            for face in lattice:
                ISIN=False
                for lattice2 in Lattices2:
                    for face2 in lattice2:
                        isin=False
                        for edge2,edge in zip(face2,face):
                            if abs(edge.dot(edge2)-edge.dot(edge))<eps:
                                isin=True
                            else:
                                isin=False
                        if isin:
                            ISIN=True
                        
                if not ISIN:
                    Lat.append(face)
                    
            Lattices2.append(Lat)
    return Lattices2


def generate_product_lattice_points(F,Parentlattice_points,Q=np.eye(3)):
    """
    generate_product_lattice_points - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """
    #F-deformation gradient
    Productlattice_points=[]
    
    for v in Parentlattice_points:
        p0=np.array([vi[0] for vi in v])
        p1=np.array([vi[1] for vi in v])
#        if not (p0==0).all():
#            pn0=p0/np.sqrt(p0.dot(p0))
#            p0B19p = np.ndarray.tolist(pn0.dot(F.dot(p0))*pn0)
#        else:
#            p0B19p = np.ndarray.tolist(p0);
#            
#        if not (p1==0).all():
#            pn1=p1/np.sqrt(p1.dot(p1))
#            p1B19p = np.ndarray.tolist(pn1.dot(F.dot(p1))*pn1)
#        else:
#            p1B19p = np.ndarray.tolist(p1);
            
        p0B19p=Q.dot(F.dot(p0))
        p1B19p=Q.dot(F.dot(p1))
        Productlattice_points.append([[p0B19pi,p1B19pi] for p0B19pi,p1B19pi in zip(p0B19p,p1B19p) ])
    
    return Productlattice_points

def generate_product_lattice_faces(F,Parentlattices):

    """
    Generate faces for two lattices.
    
    Creates face polygons for both unit cells.
    
    Input:
        L1, L2 (array 3×3): Lattice matrices
    
    Output:
        tuple: (faces1, faces2) - lists of face vertices
    
    Usage Example:
        >>> L1 = cubic_lattice_vec(3.0)
        >>> L2 = tetragonal_lattice_vec(3.0, 3.0, 4.0)
        >>> faces1, faces2 = generate_product_lattice_faces(L1, L2)
    
    Notes:
        - Returns faces for both lattices
        - Can use different colors/transparency
        - Visualizes structural relationship
    """
    #F-deformation gradient
    Productlattices=[]
    for lattice in Parentlattices:
        plattice=[]
        for face in lattice:
            pface=[]
            for edge in face:
                pface.append(F.dot(edge))
            plattice.append(pface)
        Productlattices.append(plattice)
    
    return Productlattices



def zero_normal_strains(Strain, mcircles,VV,normdiri,phi_around_normdiri,Parent_xyz2hkl):

    """
    Find planes with zero normal strain in given strain state.
    
    Calculates planes where normal strain εₙ = nᵀ·ε·n = 0.
    These planes experience only shear.
    
    Input:
        strain_tensor (array 3×3): Strain tensor
    
    Output:
        list: List of plane normals with zero normal strain
    
    Usage Example:
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
    
    Notes:
        - Important in transformation theory
        - Related to habit planes
        - Used in invariant plane strain analysis
    
    Formula:
        εₙ = nᵀ·ε·n = 0
        Solve eigenvalue problem
    """
    Re =  mcircles['R13']
    dRe = mcircles['C13']
    #phi_around_normdiri=np.linspace(0.,360,361);
    sheardiri=np.cross(normdiri,[0.,1.,0.])
    an = np.arccos(dRe/Re)/2.*180./np.pi        

    #strains, planes for orientation where normal strain==0 and \
    #for 90deg oriented plane 
    Shear={'rotabout':'V2','angle':[],'strainmag':[],'inprincipaldir':[],'inlatticedir':[],'hkldir':[]}
    Normal={'rotabout':'V2','angle':[],'strainmag':[],'inprincipaldir':[],'inlatticedir':[],'hkldir':[]}
    #strains in these planes rotated around normdiri by phi_around_normdiri
    Shears={'rotabout':'V1','angle':[],'angleaboutV2':[],'strainmag':[],'inprincipaldir':[],'inlatticedir':[],'hkldir':[]}
    Normals={'rotabout':'V1','angle':[],'angleaboutV2':[],'strainmag':[],'inprincipaldir':[],'inlatticedir':[],'hkldir':[]}
    InPlaneShears={'rotabout':[],'angle':[],'angleaboutV2':[],'strainmag':[],'inprincipaldir':[],'inlatticedir':[],'hkldir':[]}
    InPlaneNormals={'rotabout':[],'angle':[],'angleaboutV2':[],'strainmag':[],'inprincipaldir':[],'inlatticedir':[],'hkldir':[]}


    Shear['strainmag'].append([])
    Normal['strainmag'].append([])
    Shear['inprincipaldir'].append([])
    Normal['inprincipaldir'].append([])
    Shear['inlatticedir'].append([])
    Normal['inlatticedir'].append([])
    Shear['hkldir'].append([])
    Normal['hkldir'].append([])
    Shear['angle'].append([])
    Normal['angle'].append([])

    
    #sheardiri=[0,0,1]

    for angle in [an,an+90.]:
#        normdiri=[0,0,1]
#        sheardiri=[1,0,0]
        
        R=passive_rotation(angle,'y',deg=True)      

        Shear['inprincipaldir'][-1].append(R.dot(sheardiri))      
        Normal['inprincipaldir'][-1].append(R.dot(normdiri)) 
        Shear['inlatticedir'][-1].append(np.matmul(VV,R).dot(sheardiri))
        Normal['inlatticedir'][-1].append(np.matmul(VV,R).dot(normdiri)) 
        Shear['hkldir'][-1].append(xyz2fractional(Parent_xyz2hkl,np.matmul(VV,R).dot(sheardiri)))
        Normal['hkldir'][-1].append(xyz2fractional(Parent_xyz2hkl,np.matmul(VV,R).dot(normdiri)))
        
        Shear['strainmag'][-1].append(Shear['inlatticedir'][-1][-1].dot(Strain.dot(Normal['inlatticedir'][-1][-1])))
        Normal['strainmag'][-1].append(Normal['inlatticedir'][-1][-1].dot(Strain.dot(Normal['inlatticedir'][-1][-1])))
        
        Shear['angle'][-1].append(angle)
        Normal['angle'][-1].append(angle)
        InPlaneShears['angleaboutV2'].append(angle)
        InPlaneNormals['angleaboutV2'].append(angle)

        InPlaneShears['strainmag'].append([])
        InPlaneNormals['strainmag'].append([])
        InPlaneShears['inprincipaldir'].append([])
        InPlaneNormals['inprincipaldir'].append([])
        InPlaneShears['inlatticedir'].append([])
        InPlaneNormals['inlatticedir'].append([])
        InPlaneShears['hkldir'].append([])
        InPlaneNormals['hkldir'].append([])
        InPlaneShears['angle'].append([])
        InPlaneNormals['angle'].append([])

        if normdiri==[1,0,0]:
            InPlaneShears['rotabout']='V1'
            InPlaneNormals['rotabout']='V1'
        else:
            InPlaneShears['rotabout']='V3'
            InPlaneNormals['rotabout']='V3'

        StrainPrincipal=np.matmul(np.matmul(VV.T,Strain),VV)
        Strain_in_zero = np.matmul(np.matmul(R,StrainPrincipal),R.T)
        eps=1e-10;
        if Strain_in_zero[0,0]<eps:
            rotaxisidx = 2
#            InPlaneShears['rotabout'].append('V3')
#            InPlaneNormals['rotabout'].append('V3')
            inplanesheardiri=[1.,0.,0.]  
        else:
            rotaxisidx = 0
#            InPlaneShears['rotabout'].append('V1')
#            InPlaneNormals['rotabout'].append('V1')
            inplanesheardiri=[0.,0.,1.]  

        inplanenormdiri=[0.,1.,0.]   
        for phii in phi_around_normdiri:

            InPlaneShears['angle'][-1].append([phii])
            InPlaneNormals['angle'][-1].append([phii])

#            if normdiri==[1,0,0]:
#                R2=passive_rotation(phii,'x',deg=True) 
#            else:
#                R2=passive_rotation(phii,'z',deg=True)
            if rotaxisidx == 0:
                R2=passive_rotation(phii,'x',deg=True) 
                Norm=R.T.dot([1.,0.,0.])
            else:
                R2=passive_rotation(phii,'z',deg=True)
                Norm=R.T.dot([0.,0.,1.])

            ShDir = R2.dot(R.dot(inplanesheardiri))
            NoDir = R2.dot(R.dot(inplanenormdiri))
                
            ShDirg = VV.dot(ShDir)
            NoDirg = VV.dot(NoDir)

#            ShDirg = np.matmul(VV,R).dot(R2.dot(sheardiri))   
#            NoDirg = np.matmul(VV,R).dot(NoDir)
#            Strain_in_zero.dot()
#            Strain_in_zero_rot = np.matmul(np.matmul(R2,Strain_in_zero),R2.T)
#            vv=np.array([0,1,0])
#            vv.dot(Strain_in_zero_rot.dot([0,0,1]))
#            NoDir.dot(StrainPrincipal.dot(Norm))
#                                ShDir.dot(StrainPrincipal.dot(Norm))
#            
#            np.matmul(np.matmul(R2,Strain_in_zero),R2.T)
            
            InPlaneShears['strainmag'][-1].append(ShDirg.dot(Strain.dot(NoDirg)))
            InPlaneNormals['strainmag'][-1].append(NoDirg.dot(Strain.dot(NoDirg)))
            InPlaneShears['inprincipaldir'][-1].append(ShDir)
            InPlaneNormals['inprincipaldir'][-1].append(NoDir)
            InPlaneShears['inlatticedir'][-1].append(ShDirg)
            InPlaneNormals['inlatticedir'][-1].append(NoDirg)
            InPlaneShears['hkldir'][-1].append(xyz2fractional(Parent_xyz2hkl,ShDirg))
            print(NoDirg)
            InPlaneNormals['hkldir'][-1].append(xyz2fractional(Parent_xyz2hkl,NoDirg))

    return     Shear,Normal,InPlaneShears,InPlaneNormals,[an,an+90.]


def strains_along_13mohrcirle(Strain,VV,normdiri,phi_around_V2,Parent_xyz2hkl):

    """
    Calculate strains along maximum Mohr's circle.
    
    Computes normal and shear strain components for points on the
    largest Mohr's circle (ε₁-ε₃ circle).
    
    Input:
        mohr_result (dict): Output from mohr_circles()
        n_points (int): Number of points to sample (default: 100)
    
    Output:
        dict: {
            'normal': array - normal strain values
            'shear': array - shear strain values
            'angles': array - rotation angles (radians)
        }
    
    Usage Example:
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
    
    Notes:
        - Maximum circle between ε₁ and ε₃
        - Shows all possible strain states
        - Used in failure analysis
    
    Formula:
        εₙ = center + radius·cos(2θ)
        γ/2 = radius·sin(2θ)
    """
    ShearsOnCircle={'rotabout':'V2','angle':[],'strainmag':[],'inprincipaldir':[],'inlatticedir':[],'hkldir':[]}
    NormalsOnCircle={'rotabout':'V2','angle':[],'strainmag':[],'inprincipaldir':[],'inlatticedir':[],'hkldir':[]}

    
    sheardiri=np.cross(normdiri,[0.,1.,0.])
    
    for angle in phi_around_V2:
#        ShearsOnCircle['inprincipaldir'].append([])
#        NormalsOnCircle['inprincipaldir'].append([])
#        ShearsOnCircle['inlatticedir'].append([])
#        NormalsOnCircle['inlatticedir'].append([])
#        ShearsOnCircle['hkldir'].append([])
#        NormalsOnCircle['hkldir'].append([])
        
        ShDir = []
        NoDir = []
        ShDirg = []
        NoDirg = []
        ShDirhkl = []
        NoDirhkl = []
        a=[]
        for halfcircle in [0,90.]:
        
            anglei=angle+halfcircle
    
            R=passive_rotation(anglei,'y',deg=True)      
    
            ShDir.append(R.dot(sheardiri))
            NoDir.append(R.dot(normdiri))
            
                
            ShDirg.append(np.matmul(VV,R).dot(sheardiri))   
            NoDirg.append(np.matmul(VV,R).dot(normdiri))
            ShDirhkl.append(xyz2fractional(Parent_xyz2hkl,ShDirg[-1]))
            NoDirhkl.append(xyz2fractional(Parent_xyz2hkl,NoDirg[-1]))
    
        isin=False
        for nn in NormalsOnCircle['hkldir']:
            if (nn[1]==NoDirhkl[1]).all() and (nn[0]==NoDirhkl[0]).all():
                isin=True
                #print(nn[0])
        if not isin and NoDirhkl[1].dot(NoDirhkl[0])==0.:
        #if not isin:
            
            ShearsOnCircle['strainmag'].append([])
            NormalsOnCircle['strainmag'].append([])
            ShearsOnCircle['angle'].append([])
            NormalsOnCircle['angle'].append([])

            ShearsOnCircle['angle'][-1].append([angle,angle+halfcircle])
            NormalsOnCircle['angle'][-1].append([angle,angle+halfcircle])
            for sh,no in zip(ShDirg,NoDirg): 
                ShearsOnCircle['strainmag'][-1].append(sh.dot(Strain.dot(no)))
                NormalsOnCircle['strainmag'][-1].append(no.dot(Strain.dot(no)))
            ShearsOnCircle['inprincipaldir'].append(ShDir)
            NormalsOnCircle['inprincipaldir'].append(NoDir)
            ShearsOnCircle['inlatticedir'].append(ShDirg)
            NormalsOnCircle['inlatticedir'].append(NoDirg)
            ShearsOnCircle['hkldir'].append(ShDirhkl)
            NormalsOnCircle['hkldir'].append(NoDirhkl)
            
            
    return ShearsOnCircle, NormalsOnCircle


def select_crystal_planes(NormalsOnCircle,ShearsOnCircle,maxhkl):
    

    """
    Generate list of low-index crystallographic planes.
    
    Creates planes (hkl) with indices up to max_index for analysis.
    
    Input:
        lattice (array 3×3): Lattice matrix
        max_index (int): Maximum Miller index value
    
    Output:
        list: List of (h,k,l) tuples
    
    Usage Example:
        >>> L = cubic_lattice_vec(3.0)
        >>> planes = select_crystal_planes(L, max_index=2)
        >>> 
        >>> print(f"Generated {len(planes)} planes")
        >>> for hkl in planes[:10]:
        ...     print(f"({hkl[0]},{hkl[1]},{hkl[2]})")
    
    Notes:
        - Low-index planes are physically important
        - Used in diffraction analysis
        - Filters out equivalent planes
    """
    #NormalsAlli=[]
    idxmore=[]
    inc=-1
    #maxhkl=5
    for norm in NormalsOnCircle['hkldir']:
        inc+=1
        uvwi=norm[0]
        if uvwi[0]%1==0 and uvwi[1]%1==0 and uvwi[2]%1==0 and max(abs(uvwi))<=maxhkl:
            idxmore.append(inc)
            #NormalsAlli.append([np.round(norm[0]),np.round(norm[1])])
        else:
            uvwi=norm[1]
            if uvwi[0]%1==0 and uvwi[1]%1==0 and uvwi[2]%1==0 and max(abs(uvwi))<=maxhkl:
                idxmore.append(inc)
                #NormalsAlli.append([np.round(norm[0]),np.round(norm[1])])
            
    SelShearsOnCircle={}
    SelNormalsOnCircle={}
    for key in ShearsOnCircle.keys():
        if type(ShearsOnCircle[key])==list:
            SelShearsOnCircle[key]=[]
            SelNormalsOnCircle[key]=[]
        else:
            SelShearsOnCircle[key]=ShearsOnCircle[key]
            SelNormalsOnCircle[key]=NormalsOnCircle[key]
            
    for key in SelShearsOnCircle.keys():
        if type(SelShearsOnCircle[key])==list:
            for i in idxmore:
                SelShearsOnCircle[key].append(ShearsOnCircle[key][i])
                SelNormalsOnCircle[key].append(NormalsOnCircle[key][i])
    return  SelShearsOnCircle, SelNormalsOnCircle           
                


def write_mohr_planes(ax,Upperhalftext,Lowerhalftext,colors,markersize=8,markeredgewidth=2):

    """
    Write Mohr circle analysis results to file.
    
    Saves principal strains, circles, and plane-specific strains.
    
    Input:
        filename (str): Output file path
        mohr_result (dict): From mohr_circles()
        planes (list): Plane normals analyzed
        lattice (array 3×3): Lattice matrix
    
    Output:
        None (writes to file)
    
    Usage Example:
        >>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> planes = [[1,0,0], [1,1,0], [1,1,1]]
        >>> L = cubic_lattice_vec(3.0)
        >>> 
        >>> write_mohr_planes('mohr_analysis.txt', result, planes, L)
        >>> # Creates text file with complete analysis
    
    Notes:
        - Formatted text output
        - Includes principal values
        - Lists strain on each plane
        - Suitable for reports
    """
    #write planes
    uppery=1
    lowery=0.45
    xall=0    
    inc=-1;
    for ut,lt in zip(Upperhalftext,Lowerhalftext):          
        if type(colors) is list:
            if len(colors)>=len(Upperhalftext):
                inc+=1
            else:
                colors=[colors[0]]
                inc=0
        else:
            colors=[colors]
            inc=0
        ax.text(xall,uppery,ut,color=colors[inc],fontsize=10,fontweight="bold")
        ax.text(xall,lowery,lt,color=colors[inc],fontsize=10,fontweight="bold")
        uppery-=0.035
        lowery-=0.035

#    uppery-=0.035
#    lowery-=0.035
    xall+= 1./len(Upperhalftext)/4   
    inc=-1;
    for ut,lt in zip(Upperhalftext,Lowerhalftext):     
        if type(colors) is list:
            if len(colors)>=len(Upperhalftext):
                inc+=1
            else:
                colors=[colors[0]]
                inc=0
        else:
            colors=[colors]
            inc=0
        markerfacecolor='None'
        ax.plot(xall,uppery,'o',markerfacecolor=markerfacecolor,markeredgecolor=colors[inc],\
                     markeredgewidth=markeredgewidth,markersize=markersize+4)
        markerfacecolor=colors[inc]
        ax.plot(xall,lowery,'o',markerfacecolor=markerfacecolor,markeredgecolor=colors[inc],\
                     markeredgewidth=markeredgewidth,markersize=markersize)
        xall+=1./len(Upperhalftext)
    ax.plot([0,1],[0.5,0.5],'k')
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])



def write_lattice_correspondence(ax,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,Product_lattice,Parent_lattice,FontSize=10,Fontweight="bold"):

    """
    Write lattice correspondence analysis to file.
    
    Saves correspondence matrix and lattice mismatch information.
    
    Input:
        filename (str): Output file path
        C (array 3×3): Correspondence matrix
        L1, L2 (array 3×3): Lattice matrices
        phase1, phase2 (str): Phase names
    
    Output:
        None (writes to file)
    
    Usage Example:
        >>> L_B2 = cubic_lattice_vec(3.015)
        >>> L_B19p = monoclinic_lattice_vec(2.889, 4.120, 4.622, 96.8)
        >>> C = lattice_correspondence(L_B19p, L_B2)
        >>> 
        >>> write_lattice_correspondence('correspondence.txt', C,
        ...                               L_B19p, L_B2,
        ...                               'B19 martensite', 'B2 austenite')
    
    Notes:
        - Formatted report
        - Includes mismatch analysis
        - Volume change calculation
        - Suitable for documentation
    """
    #write lattice correspondence
    tt=''
    for vm in [[1,0,0],[0,1,0],[0,0,1]]:
        va = Product_uvw_2_Parent_uvw_all_norm.dot(vm)
        vmreal = Product_uvw2xyz.dot(vm)
        
        tt += r''+vec2string(va, digits=0)+r'$^{\mathbf{uvw}}_{\mathbf{'+Parent_lattice+'}}$'+r' -> '+ \
        vec2string(vm, digits=0)+r"$^{\mathbf{uvw}}_{\mathbf{"+Product_lattice+"}}$="+\
        vec2string(vmreal,digits=3)+r"$^{\mathbf{xyz}}_{\mathbf{"+Product_lattice+"}}$\n"
            
    tt=tt[:-1]
    ax.text(0,0,tt,color='k',fontsize=FontSize,fontweight=Fontweight)
    
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])


def generate_lattite_atom_positions(atoms_xyz_position,uvw2xyz,S=1,Q=np.eye(3),R=np.eye(3),shift=np.zeros(3),xlim=[],ylim=[],zlim=[]):
    """
    generate_lattite_atom_positions - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """
    S_range = list(range(-S,S+1)) 
    Points=[]
    eps=1e-5
    for atoms in atoms_xyz_position:
        Points1=[]
        for atom in atoms:
            triplets = list(itertools.product(S_range, repeat=3)) 
            triplets = np.array(triplets) 
            triplets = triplets.T
            for ii,lim in zip(range(0,3),[xlim,ylim,zlim]):
                if len(lim)>0:
                    idxs = np.where(triplets[ii,:]<lim[0])[0]
                    triplets = np.delete(triplets,idxs,1)
                    if (atom==[0.,0.,0.]).all():
                        idxs = np.where(triplets[ii,:]>lim[1])[0]
                        triplets = np.delete(triplets,idxs,1) 
                    elif atom[ii]==0:
                        idxs = np.where(triplets[ii,:]>lim[1])[0]
                        triplets = np.delete(triplets,idxs,1) 
                    else:
                        idxs = np.where(triplets[ii,:]>lim[1]-1)[0]
                        triplets = np.delete(triplets,idxs,1)

            points = R.dot(uvw2xyz.dot(triplets)+np.repeat([atom],triplets.shape[1],axis=0).T)
            #print(points.shape)
            points+=np.tile(np.array([shift]).T,(1,points.shape[1]))
            #print(points.shape)
#            for ii,lim,ei in zip(range(0,3),[xlim,ylim,zlim],np.eye(3)):
#                if len(lim)>0:
#                    aei = sum(uvw2xyz.dot(ei))
#                    idxs = np.where(points[ii,:]<aei*lim[0])[0]
#                    points = np.delete(points,idxs,1)
#                    idxs = np.where(points[ii,:]>aei*lim[1])[0]
#                    points = np.delete(points,idxs,1)
            points=np.matmul(Q,points)
                
        #    if len(Max)>0:
        #        sum1=np.sum((points>np.repeat([Max],triplets.shape[1],axis=0).T).astype(int),axis=0)
        #        sum1=sum1+np.sum((points<np.repeat([Min],triplets.shape[1],axis=0).T).astype(int),axis=0)
        #        idx = np.where(sum1==0)
        #        points=points[:,idx]
#            Max=[max(p) for p in points]
#            Min=[min(p) for p in points]
#            points = np.matmul(R,points)
            Points1.append(points)
        Points.append(Points1)
    return Points

def generate_lattice_vectors(Points,uvw2xyz,S=1,Q=np.eye(3),xlim=[],ylim=[],zlim=[],fitpoints=False,shift=np.zeros(3)):
    """
    generate_lattice_vectors - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """
    S_range = list(range(-S,S+1)) 
    atom=np.array([0., 0., 0.]);
    if fitpoints:
        allPoints=np.hstack([p for points in Points for p in points])+np.array([shift]).T
    for points in Points[0][0]*0.:
        triplets = list(itertools.product(S_range, repeat=3)) 
        triplets = np.array(triplets) 
        triplets = triplets.T
        for ii,lim in zip(range(0,3),[xlim,ylim,zlim]):
            
            if len(lim)>0:
                #print('ok')
                idxs = np.where(triplets[ii,:]<lim[0])[0]
                triplets = np.delete(triplets,idxs,1)
                if (atom==[0.,0.,0.]).all():
                    idxs = np.where(triplets[ii,:]>lim[1])[0]
                    #print(ii)
                    #print(triplets[ii,idxs])
                    triplets = np.delete(triplets,idxs,1)                        
        
        #print('==============')
        points = uvw2xyz.dot(triplets)+np.repeat([atom],triplets.shape[1],axis=0).T
        #points+=np.tile(np.array([shift]).T,(1,points.shape[1]))
        points2 = np.eye(3).dot(triplets)+np.repeat([atom],triplets.shape[1],axis=0).T
    Max=[max(p) for p in points]
    Min=[min(p) for p in points]
    
    LatticeVectors=[]        
    for point,point2 in zip(points.T,points2.T):   
        LatticeVectors1=[]
        for ii,lattice_vec,lim in zip([0,1,2],uvw2xyz.T,[xlim,ylim,zlim]):
            vec=np.array([lattice_vec*0,lattice_vec]).T+np.repeat([point],2,axis=0).T
            #print('------------------')
            #print(np.array([shift]).T.shape)
            #sum1=np.sum((vec>np.repeat([Max],2,axis=0).T).astype(int),axis=0)
            #sum1=sum1+np.sum((vec<np.repeat([Min],2,axis=0).T).astype(int),axis=0)
            #idx = np.where(sum1>0)[0]
            #print(lim)
            if len(lim)==0:
                lim=[0,point2[ii]+10]
            if point2[ii]<lim[1]:#True:#not len(idx)>0:
                vec=np.matmul(Q,vec)
                if fitpoints:
                    fit=False
                    for idx2 in [0,1]:
                        #print(np.min(np.linalg.norm(allPoints-np.array([vec[:,idx2]]).T,axis=0)))
                        #print(vec)
                        if np.min(np.linalg.norm(allPoints-np.array([vec[:,idx2]]).T,axis=0))<1e-1:
                            fit=True
                            
                    #fit=True
                    if fit:
                        LatticeVectors1.append(vec)
                    #else:
                        #print(vec)
                        #print('======================================')
                        #print(allPoints)
                        #print('======================================')
                else:
                    LatticeVectors1.append(vec)
                #print(vec)
                #print('===============')
        LatticeVectors.append(LatticeVectors1)
        


    return LatticeVectors
    

def select_atomic_plane(LatticePoints,normal,eps=1e-1,shift=0.,eps2=None):

    """
    Select atoms lying on or near a crystallographic plane.
    
    Finds atoms within tolerance distance from specified plane.
    
    Input:
        atoms (array N×3): Atomic positions
        plane_normal (array [3]): Plane normal vector
        plane_point (array [3]): Point on plane
        tolerance (float): Distance tolerance (Ångströms)
    
    Output:
        numpy.ndarray: Indices of atoms on plane
    
    Usage Example:
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
    
    Notes:
        - Uses point-to-plane distance formula
        - Tolerance accounts for numerical error
        - Returns indices, not coordinates
    
    Formula:
        distance = |n · (p - p₀)| / ||n||
    """
    if not isinstance(normal, np.ndarray):
        normal=np.array(normal);
    normal=normal/np.sqrt(normal.dot(normal))
    PointsOut=[]    
    for Points1 in LatticePoints:
        PointsOut1=[]
        #print(Points1)
        for points in Points1:   
            #point_proj_n = normal.dot(points/np.sqrt(np.sum(points**2,axis=0)))
            point_proj_n = normal.dot(points)
            
            idxs = np.where(abs(point_proj_n-shift)<(eps))[0]
            #print(point_proj_n)
            if len(idxs)>0:
                PointsOut1.append(points[:,idxs])
        PointsOut.append(PointsOut1)
    
    return PointsOut

def get_interface2d(pointOutproj,normal,horizontalproj,verticalproj):

    """
    Extract interface atoms from two lattices.
    
    Finds atoms from both structures near interface plane.
    Used in phase boundary visualization.
    
    Input:
        atoms1, atoms2 (array N×3): Atom positions for two phases
        plane_normal (array [3]): Interface normal
        plane_point (array [3]): Point on interface
        tolerance (float): Distance tolerance
    
    Output:
        tuple: (indices1, indices2) - atom indices on interface
    
    Usage Example:
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
    
    Notes:
        - Identifies interface region
        - Both phases analyzed
        - Used in transformation studies
    """
    vertices=[]
    for points in pointOutproj:
        for p in points:
            for px,py in zip(p[0],p[1]):
                vertices.append([px,py])
    vertices=np.array(vertices)
    hull = ConvexHull(vertices)
    verts=vertices[hull.vertices,:].T
    projnormal=np.array([horizontalproj.dot(normal),verticalproj.dot(normal)])
    twpoints=verts[:,np.abs(projnormal.dot(verts))<1e-10]
    return twpoints

def select_plane(LatticeVectors,normal,eps=1e-1,shift=0.,Q=np.eye(3)):
    """
    select_plane - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """
    if not isinstance(normal, np.ndarray):
        normal=np.array(normal);
    normal=normal/np.sqrt(normal.dot(normal))
    normal=Q.dot(normal)
    PointsOut=[]    
    for LatticeVectors1 in LatticeVectors:        
        for vec1 in LatticeVectors1:
            vec=Q.dot(vec1)
            for points in vec.T:   
                #point_proj_n = normal.dot(points/np.sqrt(np.sum(points**2,axis=0)))
                point_proj_n = normal.dot(points)
                
                idxs = np.where(abs(point_proj_n-shift)<(eps))[0]
                #print(point_proj_n)
                if len(idxs)>0:
                    isin=False
                    for PointOut in PointsOut:
                        if (PointOut==points).all():
                            isin=True
                            break
                    if not isin:
                        if not np.isnan(points).any():
                            PointsOut.append(points)
                        #else:
                        #    print(points)
            
    for LatticeVectors1 in LatticeVectors:        
        for vec1 in LatticeVectors1:
            vec=Q.dot(vec1)
#            point1_proj_n = normal.dot(point1)
#            point2_proj_n = normal.dot(point2)
            point_proj_n = normal.dot(vec)-shift
            #print(point_proj_n)
            #print(point_proj_n)
            if np.sign(point_proj_n[0])==-1*np.sign(point_proj_n[1]):
                points=plane_line_intersection(normal,shift,vec[:,0],vec[:,1])
                if not np.isnan(points).any():
                    PointsOut.append(points)
                else:
                    PointsOut.append(list(vec[:,0]))
                    PointsOut.append(list(vec[:,1]))
                    #print(normal)
                    #print(vec)
                    #print(points)
                    
            elif abs(point_proj_n[0])<1e-10 and abs(point_proj_n[1])<1e-10:
                #print("ooo")
                PointsOut.append(list(vec[:,0]))
                PointsOut.append(list(vec[:,1]))
                

    
    return np.array(PointsOut).T


def generate_plane_vertices(PlanePoints,normal,Q=np.eye(3),move=np.zeros(3)):
    """
    generate_plane_vertices - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """
    verts = list(zip(PlanePoints[0,:], PlanePoints[1,:],PlanePoints[2,:]))
    normal=normal/np.sqrt(normal.dot(normal))
    v2=perpendicular_vector(normal)
    v3=np.cross(normal,v2)
    Qr=np.vstack((normal,v2,v3))
    Qr.dot(v3)
    x=[]
    y=[]
    z=[]
    for points in verts:   
        v4=Qr.dot(points)
#        v4=v4/np.linalg.norm(v4)
        x.append(v4[0])
        y.append(v4[1])
        z.append(v4[2])
    
    yc=np.mean(y)
    zc=np.mean(z)
    phase=[]
    for yi,zi in zip(y,z):
        phase.append(np.arctan2(zi-zc,yi-yc))
        
        
#    normal2=normal-shift*normal
#    normal2=normal2/np.sqrt(normal2.dot(normal2))
#    v2=np.array(verts[0])-shift*normal
#    v2=v2/np.linalg.norm(v2)
#    v3=np.cross(normal2,v2)
#    Qr=np.vstack((normal2,v2,v3))
#    Qr.dot(v3)
#    phase=[]
#    x=[]
#    y=[]
#    z=[]
#    for points in verts:   
#        v4=np.array(points)-shift*normal
##        v4=v4/np.linalg.norm(v4)
#        x.append(v4[0])
#        y.append(v4[1])
#        z.append(v4[2])
#    rc=np.array([np.mean(x),np.mean(y),np.mean(z)])    
#    for points in verts:   
#        #v4r=Qr.dot(v4)+       
#        v4=np.array(points)-shift*normal-rc
#        v4=v4/np.linalg.norm(v4)
#        phase.append(np.arctan2(Qr.dot(v4)[2],Qr.dot(v4)[1]))
    
    phase=np.unwrap(phase)
    np.argsort(phase)
    verts2=np.array(verts)[np.argsort(phase)[::1],:]
    #verts3=np.vstack((verts2,verts2[0,:]))
    x = list(np.matmul(Q,verts2.T)[0,:]+move[0])
    y = list(np.matmul(Q,verts2.T)[1,:]+move[1])
    z = list(np.matmul(Q,verts2.T)[2,:]+move[2])
    verts = [list(zip(x,y,z))]
    return verts,phase


def select_atomic_region(LatticePoints,normal,side='lower',eps=1e-1,shift=0.):

    """
    Select atoms within spherical region.
    
    Finds all atoms within specified radius of center point.
    
    Input:
        atoms (array N×3): Atomic positions
        center (array [3]): Sphere center
        radius (float): Sphere radius
    
    Output:
        numpy.ndarray: Indices of atoms in region
    
    Usage Example:
        >>> atoms = generate_lattite_atom_positions(L, n1=10, n2=10, n3=10)
        >>> 
        >>> # Select atoms near point [15, 15, 15]
        >>> center = np.array([15, 15, 15])
        >>> indices = select_atomic_region(atoms, center, radius=5.0)
        >>> 
        >>> print(f"Found {len(indices)} atoms in region")
        >>> region_atoms = atoms[indices]
    
    Notes:
        - Spherical selection
        - Useful for local analysis
        - Can be used for grain selection
    
    Formula:
        ||atom - center|| < radius
    """
    if not isinstance(normal, np.ndarray):
        normal=np.array(normal);
    normal=normal/np.sqrt(normal.dot(normal))
    PointsOut=[]    
    for Points1 in LatticePoints:
        PointsOut1=[]
        for points in Points1:   
            #point_proj_n = normal.dot(points/np.sqrt(np.sum(points**2,axis=0)))
            #print(points[:,0]+shift)
            
            point_proj_n = normal.dot(points)
            
            if side=='lower':
                idxs = np.where((point_proj_n-shift)<(eps))[0]
            else:
                idxs = np.where((point_proj_n-shift)>(-eps))[0]
                #print(np.where((point_proj_n-shift)<(-eps))[0])
            #print(point_proj_n)
            PointsOut1.append(points[:,idxs])
        PointsOut.append(PointsOut1)
    
    return PointsOut
    

def get_twinning_plane_points(K1,Pointsout,horizontal,vertical):

    """
    Get atoms on twinning plane.
    
    Identifies atoms lying on twin boundary for twinning analysis.
    
    Input:
        L (array 3×3): Lattice matrix
        twin_plane_normal (array [3]): Twin plane normal (K1)
        n1, n2, n3 (int): Cell range
        tolerance (float): Distance tolerance
    
    Output:
        numpy.ndarray (N×3): Atoms on twin plane
    
    Usage Example:
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
    
    Notes:
        - K1 is twinning plane (habit plane)
        - Important for twin boundary structure
        - Used in twinning calculations
    """
    K1/=norm(K1)
    horizontal/=norm(horizontal)
    vertical/=norm(vertical)
    projnormal=np.array([horizontal.dot(K1),vertical.dot(K1)])
    allPoints=np.hstack([p for points in Pointsout for p in points])
    hull = ConvexHull(allPoints.T)
    verts=allPoints[:,hull.vertices]
    #print(verts)
    twpoints=verts[:,np.abs(projnormal.dot(verts))<1e-10]
    return twpoints
    

def an_between_vecs(v1,v2,deg=True,full2pi=False):

    """
    Calculate angle between two vectors.
    
    Computes angle using dot product formula.
    
    Input:
        v1, v2 (array [3]): Vectors
        degrees (bool): Return in degrees (default: True), else radians
    
    Output:
        float: Angle between vectors
    
    Usage Example:
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
    
    Notes:
        - Handles non-normalized vectors
        - Returns angle in [0, 180°] or [0, π]
        - Used throughout module
    
    Formula:
        cos(θ) = (v1 · v2) / (||v1|| ||v2||)
    """
    #output 0-360 deg or -180-180 angle from v1 to v2
    #if an>0 counter clockwise otherwire clockwise
    v1=v1/np.sqrt(v1.dot(v1))
    v2=v2/np.sqrt(v2.dot(v2))
    cosa=v1.dot(v2)
    crossp=np.cross(v1,v2)
    an=np.arccos(cosa)
    if np.cross(crossp,v1).dot(v2)<0:
        an*=-1
    if deg:
        an*=180./np.pi
    if full2pi and an<0:
        an+=360.
    return an

def habitplane_equation_solution(Uj,Ui,Qj,n,a,tol=1e-10):

    """
    Solve habit plane equation for phase transformation.
    
    Finds habit plane normal that satisfies invariant plane strain condition:
    F·n = λn + s  (where s is shear direction)
    
    Input:
        F (array 3×3): Deformation gradient
        s (array [3]): Shear direction
    
    Output:
        dict: {
            'normal': array [3] - habit plane normal
            'lambda': float - stretch along normal
            'valid': bool - solution exists and is physical
        }
    
    Usage Example:
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
    
    Notes:
        - Central to martensitic transformation theory
        - Habit plane shows no distortion
        - λ close to 1 indicates low-energy interface
        - Used in crystallographic theory of martensite
    
    Formula:
        (F - λI)·n = s
        Solve for n and λ
    """
    #Qj*Uj-Ui=nxa

    #habit plane according to Bhattacharya
    delta=a.T.dot(Ui).dot(inv(Ui.dot(Ui)-np.eye(3)).dot(n))
    #print(delta)
    eta=np.trace(Ui.dot(Ui))-np.linalg.det(Ui.dot(Ui))-2+norm(a)**2/2/delta
    #print(eta)
    lam=0.5*(1-np.sqrt(1+2/delta))
    #print('Volum fraction of Uj: {}'.format(lam))
    C=(Ui+lam*np.outer(n,a)).dot(Ui+lam*np.outer(a,n))
    #print(C)
    D,V = np.linalg.eig(C)
    Idxs = np.argsort(D)
        
    Lambda=D[Idxs]
    V = V[:,Idxs]

    if Lambda[0]<1 and abs(Lambda[1]-1)<tol and  Lambda[0]>1:
        HB=True
    else:
        HB=False
        
    if HB:
        habitpdata={}
        keys=['shear_angle','s','m_a','b_a','Q_a','m_m','b_m','Q_m']
        for key in keys:
            habitpdata[key]=[]
    
        for k in [-1,1]:
            hpd={}
            
            m1_a=(np.sqrt(Lambda[2])-np.sqrt(Lambda[0]))/np.sqrt(Lambda[2]-Lambda[0])*(-1*np.sqrt(1-Lambda[0])*V[:,0]+k*np.sqrt(Lambda[2]-1)*V[:,2]);
            rho =1* norm(m1_a)
            m1_a=1/rho*m1_a
            b1_a=rho*(np.sqrt(Lambda[2]*(1-Lambda[0])/(Lambda[2]-Lambda[0]))*V[:,0]+k*np.sqrt(Lambda[0]*(Lambda[2]-1)/(Lambda[2]-Lambda[0]))*V[:,2]);
            b1_a/=norm(b1_a)
            Q1_a=(np.eye(3)+np.outer(b1_a,m1_a)).dot(inv(lam*Qj.dot(Uj)+(1-lam)*Ui))

            # hpd['m_a']=m1_a#np.linalg.inv(Ui).dot(n1_a)
            # hpd['b_a']=b1_a
            # hpd['Q_a']=Q_a
            # hpd['s']=s1
            # hpd['shear_angle']=np.arctan(s1/2)*2


        # k=-1
        # m2_a=(np.sqrt(Lambda[2])-np.sqrt(Lambda[0]))/np.sqrt(Lambda[2]-Lambda[0])*(-1*np.sqrt(1-Lambda[0])*V[:,0]+k*np.sqrt(Lambda[2]-1)*V[:,2]);
        # rho =1* norm(m2_a)
        # m2_a=1/rho*m2_a
        # b2_a=rho*(np.sqrt(Lambda[2]*(1-Lambda[0])/(Lambda[2]-Lambda[0]))*V[:,0]+k*np.sqrt(Lambda[0]*(Lambda[2]-1)/(Lambda[2]-Lambda[0]))*V[:,2]);
        # b2_a/=norm(b2_a)
        # Q2_a=(np.eye(3)+np.outer(b1_a,m1_a)).dot(inv(lam*Qij.dot(Uj)+(1-lam)*Ui))
        #         twind['s']=s1
        #         twind['shear_angle']=np.arctan(s1/2)*2


def twinnedhabitplane(Ui,Uj,Qij,a1,n1,hbplanes=[],addondata={},method='bhata'):        

    """
    Calculate twinned habit plane from lattice parameters.
    
    Determines twin plane (K1) and twin direction from parent and
    twin lattice geometries.
    
    Input:
        L_parent (array 3×3): Parent lattice matrix
        L_twin (array 3×3): Twin lattice matrix  
        correspondence_matrix (array 3×3): Parent→twin correspondence
    
    Output:
        dict: {
            'K1': array [3] - twin plane (habit plane)
            'eta1': array [3] - twin direction
            'shear': float - twinning shear magnitude
        }
    
    Usage Example:
        >>> # Type-I twin in cubic
        >>> L_parent = cubic_lattice_vec(3.0)
        >>> # Twin related by mirror symmetry
        >>> C_twin = np.diag([1, 1, -1])  # Mirror in (001)
        >>> L_twin = C_twin.dot(L_parent)
        >>> 
        >>> result = twinnedhabitplane(L_parent, L_twin, C_twin)
        >>> print(f"Twin plane K1: {result['K1']}")
        >>> print(f"Twin direction η1: {result['eta1']}")
    
    Notes:
        - K1 is rational in parent
        - η1 is twin shear direction  
        - Used in deformation twinning analysis
        - Applies to mechanical twins
    """
#    if len(hbplanes.keys())==0:
#        hbplanes={}
#        titles='m_a b_a F_a Q_a Lambda epsilon gamma s alpha Li Lj var_i var_j var_i2 var_j2'
#        for key in titles.split():
#            hbplanes[key]=[]
    titles='m_a b_a b_an F_a U_a Q_a Lambda epsilon gamma s alpha Li Lj sgn'#' RotAngle RotMat'    
    if method=='bhata':
        #solution to :
        #1. QijUj - Ui =a1 x n1
        #2. Q_a(Lambda*QijUj+(1-Lambda)Ui)= I + b_a x m_a
        #Rotation of the habit plane with respect to variant i Q_a*Qi^-1 (Ui=Qi^-1*Fi)
        SOL=False
        #habit plane according to Bhattacharya: Microstructure of martensite page 113
        delta=a1.T.dot(Ui).dot(inv(Ui.dot(Ui)-np.eye(3)).dot(n1))
        eta=np.trace(Ui.dot(Ui))-np.linalg.det(Ui.dot(Ui))-2+norm(a1)**2/2/delta
        if delta<=-2 and eta >=0:
            SOL=True
            f=0.5*(1-np.sqrt(1+2/delta))
            for lam in [f,1-f]:
            #print('Volum fraction of Ui: {}, Uj: {}'.format(1-lam,lam))
                C=(Ui+lam*np.outer(n1,a1)).dot(Ui+lam*np.outer(a1,n1))
                D,V = np.linalg.eig(C)
                Idxs = np.argsort(D)  
                #print(D)
                #print(V)
                Lambda=D[Idxs]
                V = V[:,Idxs]
                #print(Lambda)
                #print(V)
                
                for kk in [-1,1]:
                    m_a=(np.sqrt(Lambda[2])-np.sqrt(Lambda[0]))/np.sqrt(Lambda[2]-Lambda[0])*(-1*np.sqrt(1-Lambda[0])*V[:,0]+kk*np.sqrt(Lambda[2]-1)*V[:,2]);
                    rho =1* norm(m_a)
                    m_a=1/rho*m_a
                    b_a=rho*(np.sqrt(Lambda[2]*(1-Lambda[0])/(Lambda[2]-Lambda[0]))*V[:,0]+kk*np.sqrt(Lambda[0]*(Lambda[2]-1)/(Lambda[2]-Lambda[0]))*V[:,2]);
                    Q_a=(np.eye(3)+np.outer(b_a,m_a)).dot(inv(lam*Qij.dot(Uj)+(1-lam)*Ui))
                    F_a=Q_a.dot(lam*Qij.dot(Uj)+(1-lam)*Ui)
                    U_a=Q_a.dot(lam*Qij.dot(Uj)+(1-lam)*Ui)
                    if False:
                        W2=0.5*(F_a-F_a)
                        rv=np.array([W2[2,1],-W2[2,0],W2[1,0]])#W=np.cross(np.eye(3),rv)
                        #rotation angle
                        omega=norm(rv)
                        #rotation axis normalized
                        rvn=rv/norm(rv)
                        #Normalized skew matrix is np.cross(np.eye(3),rvn) formed with normilized elements of rv
                        #https://www.brainm.com/software/pubs/math/Rotation_matrix.pdf
                        #Rotation matrix=I*cos(omega)+sin(omega)*np.cross(rvn,np.eye(3)) + (1-cos(omega))*np.outer(rvn,rvn)
                        R=np.eye(3)*np.cos(omega)+np.sin(omega)*np.cross(np.eye(3),rvn)+(1-np.cos(omega))*np.outer(rvn,rvn)
                        RotAngle=omega
                        RotMat=R                    
                        
                    
                    #F=np.eye(3)+np.outer(b_a,m_a)
                    #F=Ainv=np.eye(3)+gamma*np.outer(s,m_a)+alpha*np.outer(m_a,m_a)

                    b_an=b_a/norm(b_a)
                    epsilon=(Lambda[2]-Lambda[0])/2
                    gamma=epsilon*norm(b_an-(b_an.dot(m_a))*m_a)
                    s=(b_an-(b_an.dot(m_a))*m_a)/norm(b_an-(b_an.dot(m_a))*m_a)
                    alpha=epsilon*b_an.dot(m_a)
                    #b_a/=norm(b_a)
                    #Li=1-f
                    #Lj=f
                    Li=1-lam
                    Lj=lam
                    hbplane={}
                    sgn=kk
                    for key in titles.split():#hbplanes.keys():
                        exec('hbplane["{}"]={}'.format(key,key))
                        #exec('hbplane["{}"].append({})'.format(key,key))
                    for key in addondata.keys():
                        hbplane[key]=addondata[key]
                    hbplanes.append(hbplane)
    elif  method=='Yongmei':
        #habit plane according to Yongmei https://www.sciencedirect.com/science/article/pii/S1359645402001234
        SOL=False
        x = sympy.Symbol('x')
        
        #Asw=(1-x)*Ui+x*Qij.dot(Uj)
        Asw=(x)*Ui+(1-x)*Qij.dot(Uj)
        fx=sympy.Matrix(Asw.T.dot(Asw)-np.eye(3)).det()
        all_terms=sympy.Poly(fx, x).all_terms()
        fx=sum(x**(len(all_terms)-n-1) * term[1] for n,term in enumerate(all_terms) if abs(term[1])>1e-10)
        sol=np.roots(fx.as_poly().coeffs())
        if np.isreal(sol[0]) and np.isreal(sol[1]):
            SOL=True
            for soli in sol:
                Aswi=np.array(sympy.Matrix(Asw).subs(x,soli)).astype(np.float64)
                F2=Aswi.T.dot(Aswi)
                D,V = np.linalg.eig(F2)
                Idxs = np.argsort(D)
                Lambda=D[Idxs]
                V = V[:,Idxs]
                for kk in [-1,1]:
                    m_a=np.sqrt((Lambda[2]-1)/(Lambda[2]-Lambda[0]))*V[:,2]-kk*np.sqrt((1-Lambda[0])/(Lambda[2]-Lambda[0]))*V[:,0]
                    b_a=np.sqrt(Lambda[0]*(Lambda[2]-1)/(Lambda[2]-Lambda[0]))*V[:,2]+kk*np.sqrt(Lambda[2]*(1-Lambda[0])/(Lambda[2]-Lambda[0]))*V[:,0]
                    b_an=np.sqrt(Lambda[0]*(Lambda[2]-1)/(Lambda[2]-Lambda[0]))*V[:,2]+kk*np.sqrt(Lambda[2]*(1-Lambda[0])/(Lambda[2]-Lambda[0]))*V[:,0]
                    epsilon=(Lambda[2]-Lambda[0])/2
                    gamma=epsilon*norm(b_a-(b_a.dot(m_a))*m_a)
                    s=(b_a-(b_a.dot(m_a))*m_a)/norm(b_a-(b_a.dot(m_a))*m_a)
                    alpha=epsilon*b_a.dot(m_a)
                    F=V*0
                    for i in range(0,3):
                        F+=Lambda[i]**0.5*np.outer(V[:,i],V[:,i])
                    
                    Ainv=np.eye(3)+gamma*np.outer(s,m_a)+alpha*np.outer(m_a,m_a)
                    Q_a=Ainv.dot(np.linalg.inv(Aswi))#Rinv
                    F_a=np.eye(3)+np.outer(epsilon*b_a, m_a)
                    U_a=np.eye(3)+np.outer(epsilon*b_a, m_a)
                    Lj=1-soli
                    Li=soli
                    hbplane={}
                    sgn=kk
                    for key in titles.split():#hbplanes.keys():
                        exec('hbplane["{}"]={}'.format(key,key))
                        #exec('hbplanes["{}"].append({})'.format(key,key))
                    for key in addondata.keys():
                        hbplane[key]=addondata[key]
                    hbplanes.append(hbplane)
    
        

    return hbplanes

def twin_equation_solution_ini(Uj,Ui,Parent_uvw2xyz,Parent_hkl2xyz,Product_uvw2xyz,Product_hkl2xyz, 
                           Parent_uvw_2_Product_uvw_rot, Parent_uvw_2_Product_uvw,Parent_hkl_2_Product_hkl,tol=1e-10,miller='greaterthanone',printlambda=False,
                           Qj=None,Qi=None):

    """
    Initialize twin equation solver with default parameters.
    
    Sets up initial guess and parameters for twin equation solution.
    
    Input:
        None
    
    Output:
        dict: Initial parameters for twin solver
    
    Usage Example:
        >>> params = twin_equation_solution_ini()
        >>> # Modify parameters as needed
        >>> params['max_iter'] = 1000
        >>> params['tolerance'] = 1e-8
        >>> # Use in twin_equation_solution()
    
    Notes:
        - Provides sensible defaults
        - Can be customized for specific systems
        - Used by twin_equation_solution()
    """
    from numpy import matmul
    from numpy import sqrt
    from numpy.linalg import inv
    from scipy.linalg import sqrtm
    from numpy.linalg import norm
    if miller=='greaterthanone':
        MIN=True
    else:
        MIN=False
    #Solves equation Q*Uj-Ui=axn
    #C=matmul(matmul(inv(Ui).T,Uj.T),matmul(Uj,inv(Ui)))
    
    #C=matmul(matmul(Uj,inv(Ui)).T,matmul(Uj,inv(Ui)))
    #print('ok')
    C=np.linalg.multi_dot([inv(Ui.T),Uj.T,Uj,inv(Ui)])
    
    
    D,V = np.linalg.eig(C)
    #V[:,1]=-1*V[:,1]
    
    matmul(C,V) - matmul(V,D*np.eye(3))
    
    #sorting according to eigenvalues Lambda[0]<Lambda[1]<Lambda[2]
    Idxs = np.argsort(D)
    
    Lambda=D[Idxs]
    #print(Lambda)
    if printlambda:
        print(Lambda)
    V = V[:,Idxs]
    
    twindata={}
    keys=['shear_angle','s','b_a', 'Rij_a', 'k','Type','n_a','a_a','Q_a','R_a','eta_a','K_a','eta_a_type','K_a_type','n_m','a_m','Q_m','R_m','eta_m','K_m','eta_m_type','K_m_type','C_m','C_a']
    for key in keys:
        twindata[key]=[]
        
    # shear_angle=[]
    # s=[]
    # n_a=[]
    # a_a=[]

    # Q_a=[]
    # R_a=[]
    # eta_a=[]
    # K_a=[]

    # n_m=[]
    # a_m=[]
    # Q_m=[]
    # R_m=[]
    # eta_m=[]
    # K_m=[]
    TWINDATA=[]
    #print(Lambda[2]-Lambda[1])
    #solution exists if Lambda[0]<1,Lambda[1]=1,Lambda[2]>1
    if Lambda[0]<1 and np.abs(Lambda[1]-1)<tol and Lambda[2]>1:
        
        #solution in the frame where Ui Uj are defioned
        #solutions are in terms of
        #twinning shear direction a
        #twinning plane normal n
        #rotation Q so that Q*Uj-Ui=axn
        #180 twinning rotation R so that  R1_a*Ui*R1_a.T=Uj, R1_a*Ui*R1_a.T=Uj

        for k in [-1,1]:
            twind={}
            
            n1_a=(sqrt(Lambda[2])-sqrt(Lambda[0]))/sqrt(Lambda[2]-Lambda[0])*(-1*sqrt(1-Lambda[0])*Ui.T.dot(V[:,0])+k*sqrt(Lambda[2]-1)*Ui.T.dot(V[:,2]));
            rho =1* norm(n1_a)
            n1_a=1/rho*n1_a
            a1_a=rho*(sqrt(Lambda[2]*(1-Lambda[0])/(Lambda[2]-Lambda[0]))*V[:,0]+k*sqrt(Lambda[0]*(Lambda[2]-1)/(Lambda[2]-Lambda[0]))*V[:,2]);
            
            
            n1_a,c = flipvector(n1_a);
            
            #c=1
            #a1_a=c*np.linalg.inv(Ui).dot(a1_a);
            a1_a=c*a1_a;
            #twind['n_a']=n1_a
            #twind['a_a']=a1_a
            
            twind['n_a']=n1_a#np.linalg.inv(Ui).dot(n1_a)
            twind['a_a']=a1_a
            twind['k']=k
            
            
            
            #np.linalg.inv(Ui).dot(a1_a)
            
            
            
            eta1_a=vector2miller(np.linalg.inv(Parent_uvw2xyz).dot(np.linalg.inv(Ui).dot(a1_a)),MIN=MIN)
            #eta1_a=vector2miller(np.linalg.inv(Parent_uvw2xyz).dot(a1_a))
            twind['eta_a']=eta1_a
            twind['eta_a_type']=eta1_a
            K1_a=vector2miller(inv(Parent_hkl2xyz).dot(n1_a),MIN=MIN)
            #K1_a=vector2miller(inv(Parent_hkl2xyz).dot(np.linalg.inv(Ui).dot(n1_a)))
            #print(K1_a)
            twind['K_a']=K1_a
            twind['K_a_type']=K1_a

            
            eta1_m = vector2miller(Parent_uvw_2_Product_uvw.dot(eta1_a),MIN=MIN)
            
            K1_m=vector2miller(Parent_hkl_2_Product_hkl.dot(K1_a),MIN=MIN)
            
            
            n1_m = Product_hkl2xyz.dot(K1_m)
            a1_m = Product_uvw2xyz.dot(eta1_m)
            #n1_m = Parent_uvw_2_Product_uvw_rot.dot(n1_a)
            #a1_m = Parent_uvw_2_Product_uvw_rot.dot(a1_a)
            n1_m=n1_m/norm(n1_m)
            a1_m=a1_m/norm(a1_m)
            #n1_m,c=flipvector2negative(n1_m)
            #a1_m=c*a1_m
            c=1.0            
            #K1_m=vector2miller(inv(Product_hkl2xyz).dot(n1_m))
            #eta1_m=vector2miller(inv(Product_uvw2xyz).dot(a1_m))
            twind['eta_m']=c*eta1_m
            twind['K_m']=c*K1_m

            twind['eta_m_type']=c*eta1_m
            twind['K_m_type']=c*K1_m

            
            # n1_m = Parent_uvw_2_Product_uvw_rot.dot(twind['n_a'])
            # a1_m = Parent_uvw_2_Product_uvw_rot.dot(twind['a_a'])
            # n1_m=n1_m/norm(n1_m)
            # a1_m=a1_m/norm(a1_m)
            # n1_m,c=flipvector2negative(n1_m)
            # a1_m=c*a1_m
            
            # eta1_m = vector2miller(np.linalg.inv(Product_uvw2xyz).dot(a1_m))
            # twind['eta_m']=eta1_m
            # print(a1_m) 
            # print(eta1_m)
            # K1_m=vector2miller(np.linalg.inv(Product_hkl2xyz).dot(n1_m))
            # twind['K_m']=K1_m
            # print(n1_m) 
            # print(K1_m)
            # print('============================================================================')
            #print(Parent_hkl_2_Product_hkl.dot(K1_a))
            #print(twind['K_m'])
            
            #eta1_m = vector2miller(Parent_uvw_2_Product_uvw.dot(eta1_a))
            #twind['eta_m']=eta1_m
            #K1_m=vector2miller(Parent_hkl_2_Product_hkl.dot(K1_a))
            #twind['K_m']=K1_m
            
            #n1_m = Product_hkl2xyz.dot(K1_m)
            #a1_m = Parent_uvw_2_Product_uvw_rot.dot(a1_a)
            #n1_m=n1_m/norm(n1_m)
            #a1_m=a1_m/norm(a1_m)
            
            # eta1_m = vector2miller(Parent_uvw_2_Product_uvw.dot(eta1_a))
            # twind['eta_m']=eta1_m
            # K1_m=vector2miller(Parent_hkl_2_Product_hkl.dot(K1_a))
            # twind['K_m']=K1_m
            
            # n1_m = Product_hkl2xyz.dot(K1_m)
            # a1_m = Product_uvw2xyz.dot(eta1_m)
            # n1_m=n1_m/norm(n1_m)
            # a1_m=a1_m/norm(a1_m)


            twind['n_m']=n1_m
            twind['a_m']=a1_m
            Q1_a=matmul(np.outer(a1_a,n1_a)+Ui,inv(Uj))
            
            #(matmul(Q1_a,Uj)-Ui)-np.outer(a1_a,n1_a)
            R1_a=-np.eye(3)+2*np.outer(a1_a,a1_a)
            twind['Q_a']=Q1_a
            twind['R_a']=R1_a
            if Qj is not None and Qi is not None:
                twind['b_a']=Qi.dot(twind['a_a'])
                twind['Rij_a']=Qi.dot(twind['Q_a']).dot(Qj.T)
                #print(twind['Rij_a'])
            else:
                twind['b_a']=None
                twind['Rij_a']=None
                
            
            R1_m=Parent_uvw_2_Product_uvw_rot.dot(R1_a.dot(inv(Parent_uvw_2_Product_uvw_rot)))
            #np.outer(a1_m,a1_m).dot(v) - it is 0 if v is perpendicular to a1_m, otherwise it is a projection of the vector into a1_m
            #therefore if v is a vector along a1_m R1_m*v+I*v=v+v=2*np.outer(a1_m,a1_m)*v=2v
            R1_m=-np.eye(3)+2*np.outer(a1_m,a1_m)
            #or R1_m=np.eye(3)-2*np.outer(n1_m,n1_m)
            #if v is perpendicular to twinning plane thne R1_m*v=-v=>-v-I*v=-2*np.outer(n1_m,n1_m)*v=-2v
            #
            #change of lattice directions in the twin
            #if vector v lies in the twinning plane np.outer(n1_m,n1_m) is zeros and directions change signs C1_m*v=-I*v+0, the signes along n1 dop not change sign
            if (K1_m != np.round(K1_m)).any() and (eta1_m == np.round(eta1_m)).all():
                #for Type II twins - eta1 is preserved - rotation around eta1
                C1_m=-np.eye(3)+2*np.outer(a1_m,a1_m)
                Type='Type II'
                #print(K1_a)
            else:
                #for other K1 and eta1 is preserved - rotation 180 around K1
                #C1_m=np.eye(3)-2*np.outer(n1_m,n1_m)
                if (twind['eta_m'] == np.round(twind['eta_m'])).all():
                    Type='compound'
                else:
                    Type='Type I'
                C1_m=-np.eye(3)+2*np.outer(a1_m,a1_m)
            
            Q1_m=Parent_uvw_2_Product_uvw_rot.dot(Q1_a.dot(inv(Parent_uvw_2_Product_uvw_rot)))
            
            twind['Q_m']=Q1_m
            twind['R_m']=C1_m
            twind['C_m']=C1_m
            twind['C_a']=C1_m
            twind['Type']=Type
            s1=norm(a1_a)*norm(inv(Ui).dot(n1_a))
            twind['s']=s1
            twind['shear_angle']=np.arctan(s1/2)*2
            
            
            
            for key in keys:
                twindata[key].append(twind[key])
            # s.append(s1)
            
            # n_a.append(n1_a)
            # a_a.append(a1_a)
            # eta_a.append(eta1_a)
            # K_a.append(K1_a)
            # Q_a.append(Q1_a)
            # R_a.append(R1_a)
            
            # n_m.append(n1_m)
            # a_m.append(a1_m)
            # eta_m.append(eta1_m)
            # K_m.append(K1_m)
            # Q_m.append(Q1_m)
            # R_m.append(R1_m)
    
        KEYS=['shear_angle','s','b_a','Rij_a','k','Type','n1_a','a1_a','n2_a','a2_a','Q_a','R_a','eta1_a','K1_a','eta1_a_type',
              'K1_a_type','eta2_a','K2_a','eta2_a_type','K2_a_type','n1_m','a1_m','n2_m','a2_m','Q_m','R_m',
              'eta1_m','K1_m','eta2_m','K2_m','C_m','C_a']
        
        for ij in [[1,0],[0,1]]:
            TWIN={}
            for key in KEYS:
                if '1' in key:
                    TWIN[key]=twindata[key.replace('1','')][ij[0]]
                elif '2' in key:
                    TWIN[key]=twindata[key.replace('2','')][ij[1]]
                else:
                    TWIN[key]=twindata[key][ij[0]]
            TWINDATA.append(TWIN)   
        
    return TWINDATA



def twin_equation_solution(Uj,Ui,L_A,Lr_A,L_M,Lr_M, 
                           R_AM, Ci_d,Ci_p,tol=1e-10,miller='greaterthanone',printlambda=False,
                           Qj=None,Qi=None):

    """
    Solve complete twin equation system.
    
    Finds K1, K2, η1, η2 twin elements from lattice parameters.
    
    Input:
        L_parent (array 3×3): Parent lattice
        L_twin (array 3×3): Twin lattice
        params (dict): Solver parameters (from twin_equation_solution_ini)
    
    Output:
        dict: {
            'K1': array [3] - composition/habit plane
            'K2': array [3] - conjugate plane
            'eta1': array [3] - shear direction in K1
            'eta2': array [3] - direction in K2
            'P': array [3] - invariant line
            'S': float - twinning shear magnitude
        }
    
    Usage Example:
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
    
    Notes:
        - Complete crystallographic solution
        - Includes both K1/η1 and K2/η2 pairs
        - P is invariant line (no rotation or stretch)
        - Used in comprehensive twin analysis
    
    Formula:
        F = I + S(η1 ⊗ K1)
        where F is deformation gradient, S is shear
    """
    #Uj transformation stretches of a possible twinned variant related to corresponding transformation gradient Fj=Qj*Uj
    #Ui transformation matrix of the matrix variant related to corresponding transformation gradient Fi=Qi*Ui
    #Qj rotational part of the transformation gradien Fj=Qj*Uj
    #Qi rotational part of the transformation gradien Fi=Qi*Ui
    #L_A/L_M matrix converting UVW->XYZ for parent phase, austenite/product phase, martensite
    #Lr_A/Lr_M matrix converting HKL->XYZ for parent phase, austenite/product phase, martensite
    #R_AM rotation matrix rotating the reference space of the austenite phase aligned with basal directions x=[1,0,0],y=[0,1,0]
    #into the space where x,y,z are aligned with austenite directions corresponding to martensite directions [1,0,0]m,[0,1,0]m,[0,0,1]m of variant i.
    #Basically R_AM.dot(inv(Ci_d))= diagonal matrix
    #Ci_d lattice correspondence for directions converting UVW of austenite with lattice corresponding UVW of martensite
    #Ci_p lattice correspondence for planes converting HKL of austenite with lattice corresponding HKL of martensite
    
    from numpy import matmul
    from numpy import sqrt
    from numpy.linalg import inv
    from scipy.linalg import sqrtm
    from numpy.linalg import norm
    if miller=='greaterthanone':
        MIN=True
    else:
        MIN=False
    #Solves equation Q*Uj-Ui=axn for stretches so the reference space Ui rotated Qi^(-1) with respect to real parent space
    #C=matmul(matmul(inv(Ui).T,Uj.T),matmul(Uj,inv(Ui)))
    
    #C=matmul(matmul(Uj,inv(Ui)).T,matmul(Uj,inv(Ui)))
    #print('ok')
    C=np.linalg.multi_dot([inv(Ui.T),Uj.T,Uj,inv(Ui)])
    
    
    D,V = np.linalg.eig(C)
    #V[:,1]=-1*V[:,1]
    
    matmul(C,V) - matmul(V,D*np.eye(3))
    
    #sorting according to eigenvalues Lambda[0]<Lambda[1]<Lambda[2]
    Idxs = np.argsort(D)
    
    Lambda=D[Idxs]
    #print(Lambda)
    if printlambda:
        print(Lambda)
    V = V[:,Idxs]
    
    twindata={}
    keys=['shear_angle','s','b_a', 'Rij_a', 'k','Type','n_a','a_a','Q_a','eta_a','K_a','eta_a_type','K_a_type','n_m','a_m','Q_m','R_m','eta_m','K_m','eta_m_type','K_m_type','C_m','C_a']
    for key in keys:
        twindata[key]=[]
        
    TWINDATA=[]
    #print(Lambda[2]-Lambda[1])
    #solution exists if Lambda[0]<1,Lambda[1]=1,Lambda[2]>1
    if Lambda[0]<1 and np.abs(Lambda[1]-1)<tol and Lambda[2]>1:
        
        #solution in the frame where Ui Uj are defioned
        #solutions are in terms of
        #twinning shear direction a
        #twinning plane normal n
        #rotation Q so that Q*Uj-Ui=axn
        #180 twinning rotation R so that  R1_a*Ui*R1_a.T=Uj, R1_a*Ui*R1_a.T=Uj

        for k in [-1,1]:
            twind={}
            
            n1_a=(sqrt(Lambda[2])-sqrt(Lambda[0]))/sqrt(Lambda[2]-Lambda[0])*(-1*sqrt(1-Lambda[0])*Ui.T.dot(V[:,0])+k*sqrt(Lambda[2]-1)*Ui.T.dot(V[:,2]));
            rho =1* norm(n1_a)
            n1_a=1/rho*n1_a
            a1_a=rho*(sqrt(Lambda[2]*(1-Lambda[0])/(Lambda[2]-Lambda[0]))*V[:,0]+k*sqrt(Lambda[0]*(Lambda[2]-1)/(Lambda[2]-Lambda[0]))*V[:,2]);
            
            #n1_a is twinning normal undistorted by transformation!! It is real normal but only lattice corresponding to martensite twinning normal            
            n1_a,c = flipvector(n1_a);
            
            #a1_a is real twinning direction in parent phase but with subtracted 
            a1_a=c*a1_a;
            
            twind['n_a']=n1_a#this is twinning normal undistorted by transformation!! It is real normal lattice corresponding to martensite twinning normal 
            twind['a_a']=a1_a#This is real twinning direction rotated Qi^(-1) with respect to real parent space
            twind['k']=k#auxiliary
            
            
            
            #np.linalg.inv(Ui).dot(a1_a)
            eta1_a=vector2miller(np.linalg.inv(L_A).dot(np.linalg.inv(Ui).dot(a1_a)),MIN=MIN)##this is twinning direction undistorted by 
            #transformation!! See np.linalg.inv(Ui).dot(a1_a). It is real twinning shear direction lattice corresponding to martensite twinning direction 
            
            #eta1_a=vector2miller(np.linalg.inv(L_A).dot(a1_a))
            twind['eta_a']=eta1_a
            twind['eta_a_type']=eta1_a
            K1_a=vector2miller(inv(Lr_A).dot(n1_a),MIN=MIN)#Miller indexes of twinning normal undistorted by transformation!! It is HKL lattice corresponding to martensite twinning normal 
            twind['K_a']=K1_a
            twind['K_a_type']=K1_a

            #now for martensite
            #variant A based on lattice correspondence
            c=1.0            
            eta1_m = vector2miller(Ci_d.dot(eta1_a),MIN=MIN)#UVW of twinning direction in the martensite system          
            K1_m=vector2miller(Ci_p.dot(K1_a),MIN=MIN)#HKL of twinning direction in the martensite system           
            twind['eta_m']=c*eta1_m
            twind['K_m']=c*K1_m

            twind['eta_m_type']=c*eta1_m
            twind['K_m_type']=c*K1_m

            
            n1_m = Lr_M.dot(K1_m)#Real space of twinning normal in the martensite system
            a1_m = L_M.dot(eta1_m)#Real space twinning direction in the martensite system
            n1_m=n1_m/norm(n1_m)
            a1_m=a1_m/norm(a1_m)
            twind['n_m']=n1_m
            twind['a_m']=a1_m
            
            #shear of the twinning
            s1=norm(a1_a)*norm(inv(Ui).dot(n1_a))
            twind['s']=s1
            twind['shear_angle']=np.arctan(s1/2)*2


            #rigid body rotation necessary to apply to streteches of j variant to make twin with variant i
            #however there are other rotations Qj Qi already subtracted from deformation gradients Fj, Fi.
            #therefore complete rotation Rij_a is calculated below if Qj, Qi provided
            Q1_a=matmul(np.outer(a1_a,n1_a)+Ui,inv(Uj))
            
            #(matmul(Q1_a,Uj)-Ui)-np.outer(a1_a,n1_a)
            #R1_a=-np.eye(3)+2*np.outer(a1_a,a1_a)
            
            twind['Q_a']=Q1_a
            #twind['R_a']=R1_a
            if Qj is not None and Qi is not None: 
                #this is real twinning shear direction accounting for the subtracted rotation
                twind['b_a']=Qi.dot(twind['a_a'])
                twind['Rij_a']=Qi.dot(twind['Q_a']).dot(Qj.T)
                #print(twind['Rij_a'])
                #the other way to reach twinning elements in martensite
                #get real twinning normal (distorted by transformation)
                n_m=twind['n_a'].dot(np.linalg.inv(Qi.dot(Ui)))
                #rotate it into martensite frame, normalize, convert to HKL and make miller
                n_m=R_AM.dot(n_m)
                n_m=n_m/norm(n_m)
                twind['n_m']=n_m
                K_m=vector2miller(inv(Lr_M).dot(n_m),MIN=MIN)
                twind['K_m']=K_m
                twind['K_m_type']=K_m
                #eta_m corresponds to real twinning shear direction  
                a_m=Qi.dot(twind['a_a'])
                a_m=R_AM.dot(a_m)
                a_m/=norm(a_m)
                twind['a_m']=a_m
                eta_m=vector2miller(inv(L_M).dot(a_m),MIN=MIN)
                twind['eta_m']=eta_m
                twind['eta_m_type']=eta_m
            else:
                twind['b_a']=None
                twind['Rij_a']=None
                
            
            #change of lattice directions in the twin
            #this is passive and active, it is symmetric since it is 180 deg rotation
            #if vector v lies in the twinning plane np.outer(n1_m,n1_m) is zeros and directions change signs C1_m*v=-I*v+0, the signes along n1 dop not change sign
            if (K1_m != np.round(K1_m)).any() and (eta1_m == np.round(eta1_m)).all():
                #for Type II twins - eta1 is preserved - rotation around eta1
                C1_m=-np.eye(3)+2*np.outer(a1_m,a1_m)
                Type='Type II'
                #print(K1_a)
            else:
                #for other K1 and eta1 is preserved - rotation 180 around K1
                C1_m=np.eye(3)-2*np.outer(n1_m,n1_m)
                if (twind['eta_m'] == np.round(twind['eta_m'])).all():
                    Type='compound'
                else:
                    Type='Type I'
                #C1_m=-np.eye(3)+2*np.outer(a1_m,a1_m)
            
            Q1_m=R_AM.dot(Q1_a.dot(inv(R_AM)))
            
            twind['Q_m']=Q1_m
            twind['R_m']=C1_m
            twind['C_m']=C1_m
            twind['C_a']=C1_m
            twind['Type']=Type
            
            
            
            for key in keys:
                twindata[key].append(twind[key])
    
        KEYS=['shear_angle','s','b_a','Rij_a','k','Type','n1_a','a1_a','n2_a','a2_a','Q_a','eta1_a','K1_a','eta1_a_type',
              'K1_a_type','eta2_a','K2_a','eta2_a_type','K2_a_type','n1_m','a1_m','n2_m','a2_m','Q_m','R_m',
              'eta1_m','K1_m','eta2_m','K2_m','C_m','C_a']
        
        for ij in [[1,0],[0,1]]:
            TWIN={}
            for key in KEYS:
                if '1' in key:
                    TWIN[key]=twindata[key.replace('1','')][ij[0]]
                elif '2' in key:
                    TWIN[key]=twindata[key.replace('2','')][ij[1]]
                else:
                    TWIN[key]=twindata[key][ij[0]]
            TWINDATA.append(TWIN)   
        
    return TWINDATA




def def_gradient_stressfree(Cd,LA, LM,CId=None):


    """
    Calculate stress-free transformation strain (deformation gradient).
    
    Computes F₀ = L2 · L1⁻¹ representing lattice transformation
    without applied stress.
    
    Input:
        L1 (array 3×3): Initial lattice matrix
        L2 (array 3×3): Final lattice matrix
    
    Output:
        numpy.ndarray (3×3): Stress-free deformation gradient F₀
    
    Usage Example:
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
    
    Notes:
        - F₀ represents pure lattice change
        - No external stress applied
        - Used in transformation theory
        - Basis for habit plane calculations
    
    Formula:
        F₀ = L₂ · L₁⁻¹
    """
    T_MA = np.empty((3,3,Cd.shape[2]))
    T_AM = np.empty((3,3,Cd.shape[2]))
    for Variant in range(Cd.shape[2]):        
        T_MA[:,:,Variant] = Cd[:,:,Variant].dot(inv(sqrtm((Cd[:,:,Variant].T.dot(Cd[:,:,Variant])))))
        T_AM[:,:,Variant] = inv(T_MA[:,:,Variant])
        
        
    #Bain strain calculation
    #Deformation gradient
    F_AM = np.empty((3,3,Cd.shape[2]))
    
    #stretch - right Cauchy-Green deformation tensor
    U_AM = np.empty((3,3,Cd.shape[2]))
    #rotation
    Q_M = np.empty((3,3,Cd.shape[2]))
    for Variant in range(Cd.shape[2]):
        if CId is not None:
            F_AM[:,:,Variant]=T_MA[:,:,Variant].dot(LM.dot(CId[:,:,Variant].dot(inv(LA))))
            #print("Using Ci_d - correspondense for austenite directions")
        else:
            F_AM[:,:,Variant]=T_MA[:,:,Variant].dot(LM.dot(inv(Cd[:,:,Variant]).dot(inv(LA))))
            #print("Using C_d - correspondense for martensite directions")

    
        #stretch matrix
        U_AM[:,:,Variant] = sqrtm(F_AM[:,:,Variant].T.dot(F_AM[:,:,Variant]));
        
        #rotation matrix
        Q_M[:,:,Variant] = F_AM[:,:,Variant].dot(inv(U_AM[:,:,Variant]));
        

    return F_AM, U_AM, Q_M, T_MA, T_AM





def def_gradient_stressfree_ini(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,Ci_d=None):


    """
    Initialize stress-free deformation gradient from lattice parameters.
    
    Convenience function to calculate F₀ from lattice parameters directly.
    
    Input:
        lattice_params_1 (tuple): (a, b, c, α, β, γ) for initial structure
        lattice_params_2 (tuple): (a, b, c, α, β, γ) for final structure
    
    Output:
        numpy.ndarray (3×3): Deformation gradient
    
    Usage Example:
        >>> # Cubic → Tetragonal
        >>> params_cubic = (3.0, 3.0, 3.0, 90, 90, 90)
        >>> params_tetra = (3.0, 3.0, 4.0, 90, 90, 90)
        >>> 
        >>> F0 = def_gradient_stressfree_ini(params_cubic, params_tetra)
        >>> print("Transformation strain:")
        >>> print(F0)
    
    Notes:
        - Wrapper around def_gradient_stressfree()
        - Accepts lattice parameters directly
        - More convenient for quick calculations
    """
    Parent_uvw2xyz_stressfree = np.array(lattice_vec(parent_lattice_param)).T

    
    Parent_uvw2xyz=Parent_uvw2xyz_stressfree
    


    parent_lattice_param_ortho=parent_lattice_param
    parent_lattice_param_ortho['alpha']=np.pi/2.
    parent_lattice_param_ortho['beta']=np.pi/2.
    parent_lattice_param_ortho['gamma']=np.pi/2.
    Parent_uvw2xyz_ortho = np.array(lattice_vec(parent_lattice_param_ortho)).T


    Product_uvw2xyz_stressfree = np.array(lattice_vec(product_lattice_param)).T;
    
    product_lattice_param_ortho=product_lattice_param.copy()
    product_lattice_param_ortho['alpha']=np.pi/2.
    product_lattice_param_ortho['beta']=np.pi/2.
    product_lattice_param_ortho['gamma']=np.pi/2.
    Product_uvw2xyz_ortho = np.array(lattice_vec(product_lattice_param_ortho)).T

    Product_uvw_2_Parent_uvw_all_norm = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    Parent_uvw_2_Product_uvw_all_norm = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    for Variant in range(Product_uvw_2_Parent_uvw_all.shape[2]):
        F_ortho = np.matmul(Product_uvw2xyz_ortho,\
               np.matmul(np.linalg.inv(Product_uvw_2_Parent_uvw_all[:,:,Variant]),np.linalg.inv(Parent_uvw2xyz_ortho)))
#        F_ortho = np.matmul(Product_uvw2xyz_ortho,\
#                            np.linalg.inv(np.matmul(Parent_uvw2xyz_ortho,Parent_uvw_2_Product_uvw_all[:,:,Variant])))                
#        F_ortho = np.matmul(np.matmul(Product_uvw2xyz_ortho,Product_uvw_2_Parent_uvw_all[:,:,Variant]),\
#                            np.linalg.inv(Parent_uvw2xyz_ortho))                
        #Stretch = sqrtm(np.matmul(F_ortho.T,F_ortho));
        #Parent_uvw_2_Product_uvw_all_norm[:,:,Variant] = np.matmul(F_ortho,inv(Stretch))
        #Product_uvw_2_Parent_uvw_all_norm[:,:,Variant] = np.linalg.inv(Parent_uvw_2_Product_uvw_all_norm[:,:,Variant])

        PP=Product_uvw_2_Parent_uvw_all[:,:,Variant]
        Product_uvw_2_Parent_uvw_all_norm[:,:,Variant] = np.matmul(PP,inv(sqrtm(np.matmul(PP.T,PP))))
        Parent_uvw_2_Product_uvw_all_norm[:,:,Variant] = np.linalg.inv(Product_uvw_2_Parent_uvw_all_norm[:,:,Variant])
        
        
    #Bain strain calculation
    #Deformation gradient
    Fv = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    
    #stretch
    Uv = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    #rotation
    Qv = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    for Variant in range(Product_uvw_2_Parent_uvw_all.shape[2]):
        #real coords of Parent in the system aligned with Parent
        #Variant=5
        #In case of deformation Product_hkl2xyz changes
        #rotation of martensite compliance to to the system aligned with Parent
        A=Parent_uvw_2_Product_uvw_all_norm[:,:,Variant]
        #Parent_uvw_2_Product_uvw_all[:,:,Variant].dot(loaddir)
        #Product_ST_rot = np.einsum('ia,jb,kc,ld,abcd->ijkl',A, A, A, A,Product_ST)
        #deformation gradient
        Product_uvw2xyz = Product_uvw2xyz_stressfree#np.matmul(Product_strain +np.eye(3),Product_e_in_parent[:,:,Variant])
    
        
        #Deformation gradient    
#        Fv[:,:,Variant] = np.matmul(Product_uvw2xyz,\
#               np.matmul(np.linalg.inv(Product_uvw_2_Parent_uvw_all[:,:,Variant]),np.linalg.inv(Parent_uvw2xyz)))
#        Fv[:,:,Variant] = np.matmul(Product_uvw2xyz,\
#               np.linalg.inv(np.matmul(Parent_uvw2xyz,Product_uvw_2_Parent_uvw_all[:,:,Variant])))

        # Fv[:,:,Variant] = np.matmul(Product_uvw2xyz,\
        #        np.linalg.inv(np.matmul(Parent_uvw2xyz,Product_uvw_2_Parent_uvw_all[:,:,Variant])))
        Fv[:,:,Variant] = np.matmul(np.matmul(Product_uvw_2_Parent_uvw_all_norm[:,:,Variant],Product_uvw2xyz),\
               np.linalg.inv(np.matmul(Parent_uvw2xyz,Product_uvw_2_Parent_uvw_all[:,:,Variant])))
        #print("Using C_d - correspondense for martensite directions")

        if Ci_d is not None:
            Fv[:,:,Variant] = Product_uvw_2_Parent_uvw_all_norm[:,:,Variant].dot(Product_uvw2xyz.dot(Ci_d[:,:,Variant])).dot(\
                   np.linalg.inv(Parent_uvw2xyz))
            #print("Using Ci_d - correspondense for austenite directions")



#        Fv[:,:,Variant] = np.matmul(np.matmul(Product_uvw2xyz,np.linalg.inv(Product_uvw_2_Parent_uvw_all[:,:,Variant])),\
#               np.linalg.inv(np.matmul(Parent_uvw2xyz,Product_uvw_2_Parent_uvw_all[:,:,Variant])))
#        Fv[:,:,Variant] = np.matmul(np.matmul(Product_uvw2xyz,np.linalg.inv(Product_uvw_2_Parent_uvw_all[:,:,Variant])),\
#               np.linalg.inv(np.matmul(Parent_uvw2xyz,np.eye(3))))
#        Fv[:,:,Variant] = np.matmul(np.matmul(Product_uvw_2_Parent_uvw_all_norm[:,:,Variant],Product_uvw2xyz),\
#               np.linalg.inv(Parent_uvw2xyz))
    
    
        #stretch matrix
        Uv[:,:,Variant] = sqrtm(np.matmul(Fv[:,:,Variant].T,Fv[:,:,Variant]));
        
        #rotation matrix
        Qv[:,:,Variant] = np.matmul(Fv[:,:,Variant],inv(Uv[:,:,Variant]));
        

    return Fv, Uv, Qv, Product_uvw_2_Parent_uvw_all_norm, Parent_uvw_2_Product_uvw_all_norm, Parent_uvw2xyz, Product_uvw2xyz


def def_gradient(Cd,LA, LM,StressT=np.zeros((3,3)),STA=np.zeros((3,3,3,3)),STM=np.zeros((3,3,3,3)),CId=None):
    """
    def_gradient - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """

 

    Parent_strain = np.einsum('ijkl,kl',STA,StressT)
    F_parent=Parent_strain+np.eye(3)
    LAStress=F_parent.dot(LA)
    
    T_MA = np.empty((3,3,Cd.shape[2]))
    T_AM = np.empty((3,3,Cd.shape[2]))
    for Variant in range(Cd.shape[2]):        
        T_MA[:,:,Variant] = Cd[:,:,Variant].dot(inv(sqrtm((Cd[:,:,Variant].T.dot(Cd[:,:,Variant])))))
        T_AM[:,:,Variant] = inv(T_MA[:,:,Variant])
        
        
    #Bain strain calculation
    #Deformation gradient
    F_AM = np.empty((3,3,Cd.shape[2]))
    
    #stretch - right Cauchy-Green deformation tensor
    U_AM = np.empty((3,3,Cd.shape[2]))
    #rotation
    Q_M = np.empty((3,3,Cd.shape[2]))
    
    Product_strain= np.empty((3,3,Cd.shape[2]))
    F_product= np.empty((3,3,Cd.shape[2]))
    Product_e_in_parent = np.empty((3,3,Cd.shape[2]))
    for Variant in range(Cd.shape[2]):        
        #transformation of the stress into martensite coordinate system
        StressT_in_Product=np.einsum('ia,jb,ab->ij',T_AM[:,:,Variant], T_AM[:,:,Variant],StressT)
        
        
        Product_strain[:,:,Variant] = np.einsum('ijkl,kl',STM,StressT_in_Product)
        #deformation gradient
        F_product[:,:,Variant] =Product_strain[:,:,Variant] +np.eye(3)
        LMStress = F_product[:,:,Variant].dot(LM)
        
        #print(LM)
        #Deformation gradient   
        Product_e_in_parent[:,:,Variant]  = T_MA[:,:,Variant].dot(LMStress)
        #print(inv(Cd[:,:,Variant]))
        if CId is not None:
            F_AM[:,:,Variant] = Product_e_in_parent[:,:,Variant] .dot(CId[:,:,Variant].dot(inv(LAStress)))
            #print("Using Ci_d - correspondense for austenite directions")
        else:
            F_AM[:,:,Variant] = Product_e_in_parent[:,:,Variant] .dot(inv(Cd[:,:,Variant]).dot(inv(LAStress)))
            #print("Using C_d - correspondense for martensite directions")
        #print(F_AM[:,:,Variant])
        #stretch matrix
        U_AM[:,:,Variant] = sqrtm(F_AM[:,:,Variant].T.dot(F_AM[:,:,Variant]));
        
        #rotation matrix
        Q_M[:,:,Variant] = F_AM[:,:,Variant].dot(inv(U_AM[:,:,Variant]));
        

    return F_AM, U_AM, Q_M, T_MA, T_AM, LAStress, LMStress, Parent_strain, Product_strain, F_parent, F_product




def def_gradient_ini(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,StressT=np.zeros((3,3)),Parent_ST=np.zeros((3,3,3,3)),Product_ST=np.zeros((3,3,3,3))):
    """
    def_gradient_ini - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """

    Parent_uvw2xyz_stressfree = np.array(lattice_vec(parent_lattice_param)).T

    Parent_strain = np.einsum('ijkl,kl',Parent_ST,StressT)
    
    F_parent=Parent_strain+np.eye(3)
    Parent_uvw2xyz=np.matmul(F_parent,Parent_uvw2xyz_stressfree)
    



    Product_uvw2xyz_stressfree = np.array(lattice_vec(product_lattice_param)).T;
    

    Product_uvw_2_Parent_uvw_all_norm = np.empty((3,3,12))
    Parent_uvw_2_Product_uvw_all_norm = np.empty((3,3,12))
    for Variant in range(Product_uvw_2_Parent_uvw_all.shape[2]):
        PP=Product_uvw_2_Parent_uvw_all[:,:,Variant]
        Product_uvw_2_Parent_uvw_all_norm[:,:,Variant] = np.matmul(PP,inv(sqrtm(np.matmul(PP.T,PP))))
        Parent_uvw_2_Product_uvw_all_norm[:,:,Variant] = np.linalg.inv(Product_uvw_2_Parent_uvw_all_norm[:,:,Variant])
        #print(F_ortho)
        #print(Product_uvw_2_Parent_uvw_all_norm[:,:,Variant] )        
        #print(np.matmul(PP,inv(sqrtm(np.matmul(PP.T,PP)))))
        #print("======================================================")
        #print(np.matmul((Stretch),inv(F_ortho)) )
        
        
    #Bain strain calculation
    #Deformation gradient
    Fv = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    
    #stretch
    Uv = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    #rotation
    Qv = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    Product_e_in_parent = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    Parent_e = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    for Variant in range(Product_uvw_2_Parent_uvw_all.shape[2]):
        #real coords of Parent in the system aligned with Parent
        #Variant=5
        #In case of deformation Product_hkl2xyz changes
        #rotation of martensite compliance to to the system aligned with Parent
        A=Parent_uvw_2_Product_uvw_all_norm[:,:,Variant]
        #Parent_uvw_2_Product_uvw_all[:,:,Variant].dot(loaddir)
        #Product_ST_rot = np.einsum('ia,jb,kc,ld,abcd->ijkl',A, A, A, A,Product_ST)
        StressT_in_Product=np.einsum('ia,jb,ab->ij',A, A,StressT)
        
        
        Product_strain = np.einsum('ijkl,kl',Product_ST,StressT_in_Product)
        #deformation gradient
        F_product =Product_strain +np.eye(3)
        Product_uvw2xyz = np.matmul(F_product,Product_uvw2xyz_stressfree)#np.matmul(Product_strain +np.eye(3),Product_e_in_parent[:,:,Variant])
        #print(Product_uvw2xyz_stressfree)
        Product_e = np.matmul(Product_uvw2xyz,np.eye(3));
        Product_e_in_parent[:,:,Variant] =  np.matmul(Product_uvw_2_Parent_uvw_all_norm[:,:,Variant],Product_e);
        
        Parent_e[:,:,Variant] =np.matmul(Parent_uvw2xyz,np.matmul(Product_uvw_2_Parent_uvw_all[:,:,Variant],np.eye(3)));
        #Deformation gradient    
        Fv[:,:,Variant] = np.matmul(Product_uvw2xyz,\
               np.matmul(np.linalg.inv(Product_uvw_2_Parent_uvw_all[:,:,Variant]),np.linalg.inv(Parent_uvw2xyz)))
        Fv[:,:,Variant] = np.matmul(Product_uvw2xyz,\
               np.linalg.inv(np.matmul(Parent_uvw2xyz,Product_uvw_2_Parent_uvw_all[:,:,Variant])))
    
        Fv[:,:,Variant] = np.matmul(np.matmul(Product_uvw_2_Parent_uvw_all_norm[:,:,Variant],Product_uvw2xyz),\
               np.linalg.inv(np.matmul(Parent_uvw2xyz,Product_uvw_2_Parent_uvw_all[:,:,Variant])))
        #Calculation of deformation gradient
        Fv[:,:,Variant] = np.matmul(Product_e_in_parent[:,:,Variant],inv(Parent_e[:,:,Variant]));    
        #stretch matrix
        Uv[:,:,Variant] = sqrtm(np.matmul(Fv[:,:,Variant].T,Fv[:,:,Variant]));
        
        #rotation matrix
        Qv[:,:,Variant] = np.matmul(Fv[:,:,Variant],inv(Uv[:,:,Variant]));
        

    return Fv, Uv, Qv, Product_uvw_2_Parent_uvw_all_norm, Parent_uvw_2_Product_uvw_all_norm, Parent_uvw2xyz, Product_uvw2xyz, Parent_strain, Product_strain, F_parent, F_product





def def_gradient_ini2(Product_uvw_2_Parent_uvw_all,parent_lattice_param, product_lattice_param,StressT=np.zeros((3,3)),Parent_ST=np.zeros((3,3,3,3)),Product_ST=np.zeros((3,3,3,3))):
    """
    def_gradient_ini2 - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """

    Parent_uvw2xyz_stressfree = np.array(lattice_vec(parent_lattice_param)).T

    Parent_strain = np.einsum('ijkl,kl',Parent_ST,StressT)
    
    F_parent=Parent_strain+np.eye(3)
    Parent_uvw2xyz=np.matmul(F_parent,Parent_uvw2xyz_stressfree)
    


    parent_lattice_param_ortho=parent_lattice_param
    parent_lattice_param_ortho['alpha']=np.pi/2.
    parent_lattice_param_ortho['beta']=np.pi/2.
    parent_lattice_param_ortho['gamma']=np.pi/2.
    Parent_uvw2xyz_ortho = np.array(lattice_vec(parent_lattice_param_ortho)).T


    Product_uvw2xyz_stressfree = np.array(lattice_vec(product_lattice_param)).T;
    
    product_lattice_param_ortho=product_lattice_param
    product_lattice_param_ortho['alpha']=np.pi/2.
    product_lattice_param_ortho['beta']=np.pi/2.
    product_lattice_param_ortho['gamma']=np.pi/2.
    Product_uvw2xyz_ortho = np.array(lattice_vec(product_lattice_param_ortho)).T

    Product_uvw_2_Parent_uvw_all_norm = np.empty((3,3,12))
    Parent_uvw_2_Product_uvw_all_norm = np.empty((3,3,12))
    for Variant in range(Product_uvw_2_Parent_uvw_all.shape[2]):
#        F_ortho = np.matmul(Product_uvw2xyz_ortho,\
#               np.matmul(np.linalg.inv(Product_uvw_2_Parent_uvw_all[:,:,Variant]),np.linalg.inv(Parent_uvw2xyz_ortho)))
#        F_ortho = np.matmul(np.eye(3,3),\
#               np.matmul(np.linalg.inv(Product_uvw_2_Parent_uvw_all[:,:,Variant]),np.linalg.inv(np.eye(3,3))))
#        Stretch = sqrtm(np.matmul(F_ortho.T,F_ortho));
#        Parent_uvw_2_Product_uvw_all_norm[:,:,Variant] = np.matmul(F_ortho,inv(Stretch))
#        #print(np.matmul(F_ortho,inv(Stretch)))
#        #print(F_ortho.T)
#        Product_uvw_2_Parent_uvw_all_norm[:,:,Variant] = np.linalg.inv(Parent_uvw_2_Product_uvw_all_norm[:,:,Variant])
        PP=Product_uvw_2_Parent_uvw_all[:,:,Variant]
        Product_uvw_2_Parent_uvw_all_norm[:,:,Variant] = np.matmul(PP,inv(sqrtm(np.matmul(PP.T,PP))))
        Parent_uvw_2_Product_uvw_all_norm[:,:,Variant] = np.linalg.inv(Product_uvw_2_Parent_uvw_all_norm[:,:,Variant])
        #print(F_ortho)
        print(Product_uvw_2_Parent_uvw_all_norm[:,:,Variant] )        
        print(np.matmul(PP,inv(sqrtm(np.matmul(PP.T,PP)))))
        print("======================================================")
        #print(np.matmul((Stretch),inv(F_ortho)) )
        
        
    #Bain strain calculation
    #Deformation gradient
    Fv = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    
    #stretch
    Uv = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    #rotation
    Qv = np.empty((3,3,Product_uvw_2_Parent_uvw_all.shape[2]))
    for Variant in range(Product_uvw_2_Parent_uvw_all.shape[2]):
        #real coords of Parent in the system aligned with Parent
        #Variant=5
        #In case of deformation Product_hkl2xyz changes
        #rotation of martensite compliance to to the system aligned with Parent
        A=Parent_uvw_2_Product_uvw_all_norm[:,:,Variant]
        #Parent_uvw_2_Product_uvw_all[:,:,Variant].dot(loaddir)
        #Product_ST_rot = np.einsum('ia,jb,kc,ld,abcd->ijkl',A, A, A, A,Product_ST)
        StressT_in_Product=np.einsum('ia,jb,ab->ij',A, A,StressT)
        
        
        Product_strain = np.einsum('ijkl,kl',Product_ST,StressT_in_Product)
        #deformation gradient
        F_product =Product_strain +np.eye(3)
        Product_uvw2xyz = np.matmul(F_product,Product_uvw2xyz_stressfree)#np.matmul(Product_strain +np.eye(3),Product_e_in_parent[:,:,Variant])
    
        
        #Deformation gradient    
        Fv[:,:,Variant] = np.matmul(Product_uvw2xyz,\
               np.matmul(np.linalg.inv(Product_uvw_2_Parent_uvw_all[:,:,Variant]),np.linalg.inv(Parent_uvw2xyz)))
        Fv[:,:,Variant] = np.matmul(Product_uvw2xyz,\
               np.linalg.inv(np.matmul(Parent_uvw2xyz,Product_uvw_2_Parent_uvw_all[:,:,Variant])))
    
    
        #stretch matrix
        Uv[:,:,Variant] = sqrtm(np.matmul(Fv[:,:,Variant].T,Fv[:,:,Variant]));
        
        #rotation matrix
        Qv[:,:,Variant] = np.matmul(Fv[:,:,Variant],inv(Uv[:,:,Variant]));
        

    return Fv, Uv, Qv, Product_uvw_2_Parent_uvw_all_norm, Parent_uvw_2_Product_uvw_all_norm, Parent_uvw2xyz, Product_uvw2xyz, Parent_strain, Product_strain, F_parent, F_product









def niti_twinning(B2_symops,B2_recsymops,B19p_recsymops,B19p_symops,Uv,Parent_uvw2xyz,Parent_hkl2xyz,Product_uvw2xyz,
                  Product_hkl2xyz,Parent_uvw_2_Product_uvw_all, Parent_hkl_2_Product_hkl_all,
                  Parent_uvw_2_Product_uvw_all_norm,miller='greaterthanone',Qv=None):

    """
    Get NiTi twinning parameters for specific variant.
    
    Returns crystallographic twinning elements for NiTi B19' martensite.
    
    Input:
        variant (str): Twinning variant ('type1', 'type2', 'compound')
    
    Output:
        dict: {
            'K1': array [3] - compound twin plane
            'eta1': array [3] - twinning direction
            'K2': array [3] - conjugate plane
            'eta2': array [3] - conjugate direction
            'shear': float - twinning shear
        }
    
    Usage Example:
        >>> # Type-I twin (most common in NiTi)
        >>> twin_data = niti_twinning('type1')
        >>> print(f"Twin plane K1: {twin_data['K1']}")
        >>> print(f"Twin direction η1: {twin_data['eta1']}")
        >>> print(f"Shear magnitude: {twin_data['shear']:.4f}")
        >>> 
        >>> # Use in deformation analysis
        >>> F_twin = np.eye(3) + twin_data['shear'] * np.outer(twin_data['eta1'], twin_data['K1'])
    
    Notes:
        - Specific to NiTi B19' martensite
        - Type-I: {111} compound twin
        - Type-II: {001} twin
        - Data from experimental measurements
    """
    L_A=Parent_uvw2xyz
    L_M=Product_uvw2xyz
    #Meteric tensor ||[x,y,x]||^2=[uvw]*G_A*[uvw]^T
    G_A=np.matmul(L_A.T,L_A)
    G_M=np.matmul(L_M.T,L_M)
    #Reciprocal Meteric tensor ||[x,y,x]||^2=d_hkl^2=1/[hkl]*Gr_A*[hkl]^T
    Gr_A = np.linalg.inv(G_A)
    Gr_M = np.linalg.inv(G_M)
    
    i=4
    j=5
    if Qv is None:
        twindata=twin_equation_solution(Uv[:,:,i],Uv[:,:,j],
                               Parent_uvw2xyz,Parent_hkl2xyz,
                               Product_uvw2xyz,Product_hkl2xyz, 
                               Parent_uvw_2_Product_uvw_all_norm[:,:,j], 
                               Parent_uvw_2_Product_uvw_all[:,:,j],
                               Parent_hkl_2_Product_hkl_all[:,:,j],tol=1e-10)
    else:
        twindata=twin_equation_solution(Uv[:,:,i],Uv[:,:,j],
                               Parent_uvw2xyz,Parent_hkl2xyz,
                               Product_uvw2xyz,Product_hkl2xyz, 
                               Parent_uvw_2_Product_uvw_all_norm[:,:,j], 
                               Parent_uvw_2_Product_uvw_all[:,:,j],
                               Parent_hkl_2_Product_hkl_all[:,:,j],tol=1e-10,miller=miller,Qj=Qv[:,:,i],Qi=Qv[:,:,j])
        #print(list(twindata[0].keys()))
        
    varnotation={0:'1',1:'1\'',2:'2',3:'2\'',4:'3',5:'3\'',6:'4',7:'4\'',8:'5',9:'5\'',10:'6',11:'6\''}
    twin_systems={}
    for twintype in ['100','001','Type I','Type II']:
        twin_systems[twintype]={}
        for key in list(twindata[0].keys()):
            twin_systems[twintype][key]=[]
        twin_systems[twintype]['vars']=[]
        twin_systems[twintype]['uvw2xyz_m']=[]
        twin_systems[twintype]['uvw2xyz_a']=[]
        twin_systems[twintype]['StrainTensor_m']=[]
        twin_systems[twintype]['F_m']=[]
        twin_systems[twintype]['F_a']=[]
        twin_systems[twintype]['RotMat']=[]
        twin_systems[twintype]['RotAngle']=[]

        twin_systems[twintype]['Tension_m']=[]
        twin_systems[twintype]['Compression_m']=[]
        twin_systems[twintype]['StrainTensor_a']=[]
        twin_systems[twintype]['Tension_a']=[]
        twin_systems[twintype]['Compression_a']=[]
        twin_systems[twintype]['Variant pairs']=[]
        twin_systems[twintype]['T_AM']=[]
        twin_systems[twintype]['Variant pairs 2']=[]
        twin_systems[twintype]['dm_a']=[]
        twin_systems[twintype]['dm_m']=[]
    #print(twin_systems[twintype])
    i=5
    j=4
    i=4
    j=5
    
    #Trnaformation twins
    Twins = np.zeros((Parent_hkl_2_Product_hkl_all.shape[2],Parent_hkl_2_Product_hkl_all.shape[2]), dtype=object)
    numtype2=0
    numtype1=0
    num100=0
    num001=0
    numall=0
    for j in range(0,Parent_hkl_2_Product_hkl_all.shape[2]): 
        for i in range(0,Parent_hkl_2_Product_hkl_all.shape[2]): 
            
            if i!=j:
                #print((i,j))
                #twin_equation_solution(Uj,Ui,...) - Solves equation Q*Uj-Ui=axn
                if Qv is None:
                    twindata=twin_equation_solution(Uv[:,:,i],Uv[:,:,j],
                                           Parent_uvw2xyz,Parent_hkl2xyz,
                                           Product_uvw2xyz,Product_hkl2xyz, 
                                           Parent_uvw_2_Product_uvw_all_norm[:,:,j], 
                                           Parent_uvw_2_Product_uvw_all[:,:,j],
                                           Parent_hkl_2_Product_hkl_all[:,:,j],tol=1e-10,miller=miller)
                    RotMat=(np.eye(3),np.eye(3))
                    RotAngle=(0.,0.)
                else:
                    twindata=twin_equation_solution(Uv[:,:,i],Uv[:,:,j],
                                           Parent_uvw2xyz,Parent_hkl2xyz,
                                           Product_uvw2xyz,Product_hkl2xyz, 
                                           Parent_uvw_2_Product_uvw_all_norm[:,:,j], 
                                           Parent_uvw_2_Product_uvw_all[:,:,j],
                                           Parent_hkl_2_Product_hkl_all[:,:,j],tol=1e-10,miller=miller,Qj=Qv[:,:,i],Qi=Qv[:,:,j])
                    RotAngle=[]
                    RotMat=[]
                    for ji in [j,i]:
                        F=Qv[:,:,ji].dot(Uv[:,:,ji])
                        W2=0.5*(F-F.T)
                        rv=np.array([W2[2,1],-W2[2,0],W2[1,0]])#W=np.cross(np.eye(3),rv)
                        #rotation angle
                        omega=norm(rv)
                        #rotation axis normalized
                        rvn=rv/norm(rv)
                        #Normalized skew matrix is np.cross(np.eye(3),rvn) formed with normilized elements of rv
                        #https://www.brainm.com/software/pubs/math/Rotation_matrix.pdf
                        #Rotation matrix=I*cos(omega)+sin(omega)*np.cross(rvn,np.eye(3)) + (1-cos(omega))*np.outer(rvn,rvn)
                        R=np.eye(3)*np.cos(omega)+np.sin(omega)*np.cross(np.eye(3),rvn)+(1-np.cos(omega))*np.outer(rvn,rvn)
                        RotAngle.append(omega)
                        RotMat.append(R)
                    RotAngle=tuple(RotAngle)
                    RotMat = tuple(RotMat)
                    
                    #print(twin_systems[twintype])
                #print(len(twindata))
                #print(i)
                if len(twindata)>0:

                    K1_arounddigit=30
                    eta1_arounddigit=30
                    arounddigit=0
                    for twin in  twindata:
                        #if twintype == '001':
                        #    print(list(twin.keys()))

                        twintype='no'
                        if (twin['eta1_m'] == np.round(twin['eta1_m'])).all():
                            uvw2xyz_m = np.matmul((2*np.outer(twin['a1_m'],twin['a1_m'])-np.eye(3)),Product_uvw2xyz)
                            uvw2xyz_a = np.matmul((2*np.outer(twin['a1_a'],twin['a1_a'])-np.eye(3)),np.matmul(np.linalg.inv(Parent_uvw_2_Product_uvw_all_norm[:,:,i]),Product_uvw2xyz))
                            if (twin['K1_m'] == np.round(twin['K1_m'])).all():
                                Twins[i,j]='Compound'
                                if (np.abs(twin['K1_m'])-np.array([0,0,1])==0).all(): #or (K_m[0]-np.array([0,0,1])==0).all():
                                    twintype='001'
                                    num001+=1
                                    #print('001: {}'.format(num001))
                                elif (np.abs(twin['K1_m'])-np.array([1,0,0])==0).all():
                                    twintype='100'
                                    num100+=1
                                    #print('100: {}'.format(num100))
                                else:
                                    numall+=1
                                    print('noncategorized twin')
                                    print('K1_m: {}, eta1_m:{}'.format(twin['K1_m'],twin['eta1_m']))        
                                    print((twin['eta1_m'] == np.around(twin['eta1_m'],decimals=arounddigit)))
                                    print(np.around(twin['eta1_m'],decimals=arounddigit))
                            else:
                                #print(twin['eta1_m'])
                                Twins[i,j]='Type II'
                                twintype='Type II'
                                K1_arounddigit=1
                                numtype2+=1
                                #print('Type II: {}'.format(numtype2))
                                #print(twin['K1_a'])
                        elif (twin['K1_m'] == np.round(twin['K1_m'])).all():
                            uvw2xyz_m = np.matmul((-2*np.outer(twin['n1_m'],twin['n1_m'])+np.eye(3)),Product_uvw2xyz)
                            uvw2xyz_a = np.matmul((-2*np.outer(twin['n1_a'],twin['n1_a'])+np.eye(3)),np.matmul(np.linalg.inv(Parent_uvw_2_Product_uvw_all_norm[:,:,i]),Product_uvw2xyz))
                            twintype='Type I'
                            Twins[i,j]='Type I'
                            eta1_arounddigit=1
                            numtype1+=1
                            #print('Type I: {}'.format(numtype1))
                            #if (twin['K1_a']==[0.,0.,1.0]).all():
                                #print('{}-{}'.format(i,j))
                        #else:
                            #numall+=1
                            #print('numall: {}'.format(numall))                        
                            
                        #if  twintype =='no':
                        #   numall+=1
                        #    print('numall: {}'.format(numall))                        
                        if twintype !='no':
                            #isin=False
                            #print(twin_systems[twintype])
                            #if twintype=='Type II':
                            #    numtype2+=1
                            #    print(numtype2)
                            twels=['K1_a_type','eta1_a_type','K2_a_type','eta2_a_type']
                            twels=['K1_m_type','eta1_m_type']
                            twels=['K1_a_type','eta1_a_type']
                            SymO=[B2_recsymops,B2_symops]
                            
                            ISIN=False
                            isin=[]
                            if len(twin_systems[twintype][twels[0]])>0:
                                isin=[]
                                for twi in range(0,len(twin_systems[twintype][twels[0]])):
                                    isin.append(False)
                                    alldiff=[]
                                    for twel,Sym in zip(twels,SymO):
                                        alldiff.append(0)
                                        for symops in Sym:
                                            twitwel = twin_systems[twintype][twel][twi]
                                            twintwel=twin[twel]                                            
                                            if 'K1' in twel:
                                                twitwel=np.around(twitwel,decimals=K1_arounddigit)
                                                twintwel=np.around(twintwel,decimals=K1_arounddigit)
                                            elif 'eta1' in twel:
                                                twitwel=np.around(twitwel,decimals=eta1_arounddigit)
                                                twintwel=np.around(twintwel,decimals=eta1_arounddigit)   

                                            if (abs(twitwel-symops.dot(twintwel))>1e-10).any():
                                                # if twintype=='Type I' and len(twin_systems[twintype][twels[0]])==1:
                                                #     print(twin_systems[twintype][twel][twi])
                                                #     print(symops.dot(twin[twel]))
                                                #     print('================================================================================')
                                                alldiff[-1]+=1
                                    if  alldiff[0]==len(SymO[0]) or alldiff[1]==len(SymO[0]):
                                        isin[-1]=False
                                    else:
                                        isin[-1]=True
                                    #K1in=[1 for K1Mi in [twin_systems[twintype]['K1_m'][twi]] if (twin['K1_m'] == K1Mi).all() or (-1*twin['K1_m'] == K1Mi).all()]
                                    #E1in=[1 for E1Mi in [twin_systems[twintype]['eta1_m'][twi]] if (twin['eta1_m'] == E1Mi).all() or (-1*twin['eta1_m'] == E1Mi).all()]
                                    #K1E1in = [1 for K1Mi,E1Mi in zip([twin_systems[twintype]['K1_m'][twi]],[twin_systems[twintype]['eta1_m'][twi]]) 
                                    #                                 if ((twin['K1_m'] == K1Mi).all() and (twin['eta1_m'] == E1Mi).all()) or 
                                    #                                 ((-1*twin['K1_m'] == K1Mi).all() and (-1*twin['eta1_m'] == E1Mi).all())] 
                                    #print(twin_systems[twintype]['K1_m'][twi])
                                    #print('=================')
                                    #if len(K1in)!=0 and len(E1in)!=0:
                                    #if len(K1E1in)!=0:
                                    #    isin[-1]=True
                                    #if twintype=='Type II':
                                    #    print(twin_systems[twintype]['eta1_m'][twi])
                                    #    print('================================================')
                            
                            keys2list=['Q_a','b_a','Rij_a','a1_a','n1_a','a2_a','n2_a','K1_a','eta1_a','eta2_a','K2_a','eta1_m','K1_m','eta2_m','K2_m']#,'eta1_a','K1_a']            
                            #keys2list=['Q_a','a1_a','n1_a','a2_a','n2_a']  
                            try:
                                idx=isin.index(True)
                                #if twintype=='Type II':
                                #    print(twin[twels[0]])
                                #    print(twin_systems[twintype][twels[0]][idx])
                                #    print(twin[twels[1]])
                                #    print(twin_systems[twintype][twels[1]][idx])
                                #    print('=============================================')
                                ISIN=True
                                #twin_systems[twintype]['Variant pairs'][idx].append((i,j))
                                #twin_systems[twintype]['Variant pairs 2'][idx].append((varnotation[i],varnotation[j]))
                                twin_systems[twintype]['Variant pairs'][idx].append((j,i))
                                twin_systems[twintype]['Variant pairs 2'][idx].append((varnotation[j],varnotation[i]))
                                twin_systems[twintype]['T_AM'][idx].append((Parent_uvw_2_Product_uvw_all_norm[:,:,j],Parent_uvw_2_Product_uvw_all_norm[:,:,i]))
                                twin_systems[twintype]['RotMat'][idx].append(RotMat)
                                twin_systems[twintype]['RotAngle'][idx].append(RotAngle)
                                
                                
                                for key in keys2list:
                                    twin_systems[twintype][key][idx].append(twin[key])
                                #twin_systems[twintype]['Q_a'][idx].append(twin['Q_a'])
                                #twin_systems[twintype]['a1_a'][idx].append(twin['a1_a'])
                                #twin_systems[twintype]['n1_a'][idx].append(twin['n1_a'])
                                #twin_systems[twintype]['a2_a'][idx].append(twin['a2_a'])
                                #twin_systems[twintype]['n2_a'][idx].append(twin['n2_a'])
                                #twin_systems[twintype]['vars'].append((i,j))
                                twin_systems[twintype]['vars'].append((j,i))
                            except:
                                pass
                            #ISIN=False        
                            
                            if not ISIN:
                                for key in list(twin.keys()):                                    
                                    if key in keys2list:
                                        twin_systems[twintype][key].append([twin[key]])
                                    else:
                                        twin_systems[twintype][key].append(twin[key])
                                #twin_systems[twintype]['vars'].append((i,j))
                                twin_systems[twintype]['vars'].append((j,i))
                                twin_systems[twintype]['uvw2xyz_m'].append(uvw2xyz_m)
                                #twin_systems[twintype]['uvw2xyz_a'].append(uvw2xyz_a)
                                
                                n2eta1K1=np.cross(twin_systems[twintype]['a1_m'][-1],twin_systems[twintype]['n1_m'][-1]);
                                Qrot=np.vstack((twin_systems[twintype]['a1_m'][-1],n2eta1K1,twin_systems[twintype]['n1_m'][-1]))
                                Qrot.dot(twin_systems[twintype]['a1_m'][-1])
                                Qri=np.linalg.inv(Qrot)
                                StrainTensor=np.array([[0,0,twin_systems[twintype]['s'][-1]/2],[0,0,0],[twin_systems[twintype]['s'][-1]/2,0,0]])
                                DefGrad=np.array([[1,0,twin_systems[twintype]['s'][-1]],[0,1,0],[0,0,1]])
                                
                                twin_systems[twintype]['StrainTensor_m'].append(np.matmul(np.matmul(Qri,StrainTensor),Qrot))
                                twin_systems[twintype]['F_m'].append(np.matmul(np.matmul(Qri,DefGrad),Qrot))
                                #if twintype=='001':
                                #    print(np.matmul(np.matmul(Qri,DefGrad),Qrot))
                                #    print(Qrot)
                                D,V = np.linalg.eig(twin_systems[twintype]['StrainTensor_m'][-1])
                                Idxs = np.argsort(D)[::-1]
                                Lambda=D[Idxs]
                                V = V[:,Idxs]
                                #if twintype == '001':
                                #    print('ok')
                                twin_systems[twintype]['Tension_m'].append(V[:,0])
                                twin_systems[twintype]['Compression_m'].append(V[:,2])
                                #Qa=np.linalg.inv(Parent_uvw_2_Product_uvw_all_norm[:,:,i])
                                #twin_systems[twintype]['StrainTensor_a'].append(Qa.dot(twin_systems[twintype]['StrainTensor_m'][-1]).dot(Qa.T))
                                #twin_systems[twintype]['Tension_a'].append(Qa.dot(twin_systems[twintype]['Tension_m'][-1]))
                                #twin_systems[twintype]['Compression_a'].append(Qa.dot(twin_systems[twintype]['Compression_m'][-1]))
                                
                                #twin_systems[twintype]['Variant pairs'].append([(i,j)])
                                #twin_systems[twintype]['Variant pairs 2'].append([(varnotation[i],varnotation[j])])
                                twin_systems[twintype]['Variant pairs'].append([(j,i)])
                                twin_systems[twintype]['T_AM'].append([(Parent_uvw_2_Product_uvw_all_norm[:,:,j],Parent_uvw_2_Product_uvw_all_norm[:,:,i])])
                                
                                twin_systems[twintype]['RotMat'].append([RotMat])
                                twin_systems[twintype]['RotAngle'].append([RotAngle])
                                
                                twin_systems[twintype]['Variant pairs 2'].append([(varnotation[j],varnotation[i])])
                                #print(twin_systems[twintype]['K1_m'][-1])
                                
                                
                                twin_systems[twintype]['dm_m'].append(get_twinning_dislocation(twin_systems[twintype]['K1_m'][-1][0],
                                                                                        twin_systems[twintype]['eta1_m'][-1][0],
                                                                                        twin_systems[twintype]['eta2_m'][-1][0],
                                                                                        L_M,G=G_M,Gr=Gr_M))
                           
                            twin_systems[twintype]['uvw2xyz_a'].append(uvw2xyz_a)
                            Qa=np.linalg.inv(Parent_uvw_2_Product_uvw_all_norm[:,:,i])
                            twin_systems[twintype]['StrainTensor_a'].append(Qa.dot(twin_systems[twintype]['StrainTensor_m'][-1]).dot(Qa.T))
                            twin_systems[twintype]['F_a'].append(Qa.dot(twin_systems[twintype]['F_m'][-1]).dot(Qa.T))
                            twin_systems[twintype]['Tension_a'].append(Qa.dot(twin_systems[twintype]['Tension_m'][-1]))
                            twin_systems[twintype]['Compression_a'].append(Qa.dot(twin_systems[twintype]['Compression_m'][-1]))                           
                            twin_systems[twintype]['dm_a'].append(get_twinning_dislocation(twin_systems[twintype]['K1_a'][-1][0],
                                                                                     twin_systems[twintype]['eta1_a'][-1][0],
                                                                                     twin_systems[twintype]['eta2_a'][-1][0],
                                                                                     L_A,G=G_A,Gr=Gr_A))   

    
    #deformation twins in martensite
    #20-1
    
    twintype='20-1'
    twin_systems[twintype]={}
    for key in list(twindata[0].keys()):
        twin_systems[twintype][key]=[]
    twin_systems[twintype]['vars']=[]
    twin_systems[twintype]['uvw2xyz_m']=[]
    twin_systems[twintype]['uvw2xyz_a']=[]
    twin_systems[twintype]['StrainTensor_m']=[]
    twin_systems[twintype]['Tension_m']=[]
    twin_systems[twintype]['Compression_m']=[]
    twin_systems[twintype]['Variant pairs']=[]
    twin_systems[twintype]['dm_a']=[]
    twin_systems[twintype]['dm_m']=[]
 
    twin_systems[twintype]['eta1_m'].append(np.array([-1,0,-2]))
    twin_systems[twintype]['K1_m'].append(np.array([2,0,-1]))
    twin_systems[twintype]['eta2_m'].append(np.array([1,0,0]))
    twin_systems[twintype]['K2_m'].append(np.array([0,0,-1]))
    uvw2xyz=Product_uvw2xyz
    hkl2xyz=Product_hkl2xyz
    a1=uvw2xyz.dot(twin_systems[twintype]['eta1_m'][0])
    twin_systems[twintype]['a1_m'].append(a1/np.linalg.norm(a1))
    a2=uvw2xyz.dot(twin_systems[twintype]['eta2_m'][0])
    twin_systems[twintype]['a2_m'].append(a2/np.linalg.norm(a2))
    n1=hkl2xyz.dot(twin_systems[twintype]['K1_m'][0])
    twin_systems[twintype]['n1_m'].append(n1/np.linalg.norm(n1))
    n2=hkl2xyz.dot(twin_systems[twintype]['K2_m'][0])
    twin_systems[twintype]['n2_m'].append(n2/np.linalg.norm(n2))
    
    twin_systems[twintype]['uvw2xyz_m'].append(np.matmul((2*np.outer(twin_systems[twintype]['a1_m'][0],twin_systems[twintype]['a1_m'][0])-np.eye(3)),uvw2xyz))
    twin_systems[twintype]['R_m'].append(twin_systems[twintype]['uvw2xyz_m'][-1].dot(inv(uvw2xyz)))
    #C1_m=-np.eye(3)+2*np.outer(n1_m,n1_m)
    #twin_systems[twintype]['C_m'].append(np.eye(3)-2*np.outer(twin_systems[twintype]['n1_m'][0],twin_systems[twintype]['n1_m'][0]))
    twin_systems[twintype]['C_m'].append(2*np.outer(twin_systems[twintype]['a1_m'][0],twin_systems[twintype]['a1_m'][0])-np.eye(3))
    twin_systems[twintype]['shear_angle'].append(2*abs(np.pi/2-np.arccos((twin_systems[twintype]['uvw2xyz_m'][0].dot(twin_systems[twintype]['eta2_m'][0])
                                         /norm(twin_systems[twintype]['uvw2xyz_m'][0].dot(twin_systems[twintype]['eta2_m'][0]))).dot(
                                             twin_systems[twintype]['a1_m'][0]))))
    twin_systems[twintype]['s'].append(np.tan(twin_systems[twintype]['shear_angle'][0]/2)*2)
    
    
    
    n2eta1K1=np.cross(twin_systems[twintype]['a1_m'][-1],twin_systems[twintype]['n1_m'][-1]);
    Qrot=np.vstack((twin_systems[twintype]['a1_m'][-1],n2eta1K1,twin_systems[twintype]['n1_m'][-1]))
    Qrot.dot(twin_systems[twintype]['a1_m'][-1])
    Qri=np.linalg.inv(Qrot)
    StrainTensor=np.array([[0,0,twin_systems[twintype]['s'][-1]/2],[0,0,0],[twin_systems[twintype]['s'][-1]/2,0,0]])
    twin_systems[twintype]['StrainTensor_m'].append(np.matmul(np.matmul(Qri,StrainTensor),Qrot))
    
    D,V = np.linalg.eig(twin_systems[twintype]['StrainTensor_m'][-1])
    Idxs = np.argsort(D)
    Lambda=D[Idxs]
    V = V[:,Idxs]
    twin_systems[twintype]['Tension_m'].append(V[:,0])
    twin_systems[twintype]['Compression_m'].append(V[:,2])
    twin_systems[twintype]['dm_m'].append(get_twinning_dislocation(twin_systems[twintype]['K1_m'][-1],
                                                            twin_systems[twintype]['eta1_m'][-1],
                                                            twin_systems[twintype]['eta2_m'][-1],
                                                            L_M,G=G_M,Gr=Gr_M))
    
    
    for symhkl,symuvw in zip( B19p_recsymops,B19p_symops):
        K1=np.round(symhkl[0:3,0:3].dot(twin_systems[twintype]['K1_m'][0]),5)
        K2=np.round(symhkl[0:3,0:3].dot(twin_systems[twintype]['K2_m'][0]),5)
        eta1=np.round(symuvw[0:3,0:3].dot(twin_systems[twintype]['eta1_m'][0]),5)
        eta2=np.round(symuvw[0:3,0:3].dot(twin_systems[twintype]['eta2_m'][0]),5)
        isin=False
        for i in range(0,len(twin_systems[twintype]['eta1_m'])):
            if (twin_systems[twintype]['eta1_m'][i]==eta1).all():
                if (twin_systems[twintype]['eta2_m'][i]==eta2).all():
                    if (twin_systems[twintype]['K1_m'][i]==K1).all():
                        if (twin_systems[twintype]['K2_m'][i]==K2).all():
                            isin=True  
        K1in=[1 for K1Mi in twin_systems[twintype]['K1_m'] if (K1 == K1Mi).all() or (-1*K1 == K1Mi).all()]
        E1in=[1 for E1Mi in twin_systems[twintype]['eta1_m'] if (eta1 == E1Mi).all() or (-1*eta1 == E1Mi).all()]
        if len(K1in)==0 or len(E1in)==0:

#        if not isin:
            twin_systems[twintype]['eta1_m'].append(eta1)
            twin_systems[twintype]['K1_m'].append(K1)
            twin_systems[twintype]['eta2_m'].append(eta2)
            twin_systems[twintype]['K2_m'].append(K2)
            a1=uvw2xyz.dot(twin_systems[twintype]['eta1_m'][-1])
            twin_systems[twintype]['a1_m'].append(a1/np.linalg.norm(a1))
            a2=uvw2xyz.dot(twin_systems[twintype]['eta2_m'][-1])
            twin_systems[twintype]['a2_m'].append(a2/np.linalg.norm(a2))
            n1=hkl2xyz.dot(twin_systems[twintype]['K1_m'][-1])
            twin_systems[twintype]['n1_m'].append(n1/np.linalg.norm(n1))
            n2=hkl2xyz.dot(twin_systems[twintype]['K2_m'][-1])
            twin_systems[twintype]['n2_m'].append(n2/np.linalg.norm(n2))
            
            twin_systems[twintype]['uvw2xyz_m'].append(np.matmul((2*np.outer(twin_systems[twintype]['a1_m'][-1],twin_systems[twintype]['a1_m'][-1])-np.eye(3)),uvw2xyz))
            twin_systems[twintype]['R_m'].append(twin_systems[twintype]['uvw2xyz_m'][-1].dot(inv(uvw2xyz)))
            
            #C1_m=-np.eye(3)+2*np.outer(n1_m,n1_m)
            #twin_systems[twintype]['C_m'].append(-np.eye(3)+2*np.outer(twin_systems[twintype]['n1_m'][0],twin_systems[twintype]['n1_m'][0]))
            twin_systems[twintype]['C_m'].append(2*np.outer(twin_systems[twintype]['a1_m'][-1],twin_systems[twintype]['a1_m'][-1])-np.eye(3))
            twin_systems[twintype]['shear_angle'].append(2*abs(np.pi/2-np.arccos((twin_systems[twintype]['uvw2xyz_m'][-1].dot(twin_systems[twintype]['eta2_m'][-1])
                                                 /norm(twin_systems[twintype]['uvw2xyz_m'][-1].dot(twin_systems[twintype]['eta2_m'][-1]))).dot(
                                                     twin_systems[twintype]['a1_m'][-1]))))
            twin_systems[twintype]['s'].append(np.tan(twin_systems[twintype]['shear_angle'][-1]/2)*2)
            n2eta1K1=np.cross(twin_systems[twintype]['a1_m'][-1],twin_systems[twintype]['n1_m'][-1]);
            Qrot=np.vstack((twin_systems[twintype]['a1_m'][-1],n2eta1K1,twin_systems[twintype]['n1_m'][-1]))
            Qrot.dot(twin_systems[twintype]['a1_m'][-1])
            Qri=np.linalg.inv(Qrot)
            StrainTensor=np.array([[0,0,twin_systems[twintype]['s'][-1]/2],[0,0,0],[twin_systems[twintype]['s'][-1]/2,0,0]])
            twin_systems[twintype]['StrainTensor_m'].append(np.matmul(np.matmul(Qri,StrainTensor),Qrot))
            
            D,V = np.linalg.eig(twin_systems[twintype]['StrainTensor_m'][-1])
            Idxs = np.argsort(D)
            Lambda=D[Idxs]
            V = V[:,Idxs]
            twin_systems[twintype]['Tension_m'].append(V[:,0])
            twin_systems[twintype]['Compression_m'].append(V[:,2])
            
            twin_systems[twintype]['dm_m'].append(get_twinning_dislocation(twin_systems[twintype]['K1_m'][-1],
                                                                    twin_systems[twintype]['eta1_m'][-1],
                                                                    twin_systems[twintype]['eta2_m'][-1],
                                                                    L_M,G=G_M,Gr=Gr_M))
            
        
    
    
    
    #deformation twins in austenite
    #114
    twintype='114'
    twin_systems[twintype]={}
    for key in list(twindata[0].keys()):
        twin_systems[twintype][key]=[]
    twin_systems[twintype]['vars']=[]
    twin_systems[twintype]['uvw2xyz_m']=[]
    twin_systems[twintype]['uvw2xyz_a']=[]
    twin_systems[twintype]['StrainTensor_a']=[]
    twin_systems[twintype]['Tension_a']=[]
    twin_systems[twintype]['Compression_a']=[]
    twin_systems[twintype]['Variant pairs']=[]
    twin_systems[twintype]['dm_a']=[]
    twin_systems[twintype]['dm_m']=[]
 
    # twin_systems[twintype]['eta1_a'].append([-1,-2,-2])
    # twin_systems[twintype]['K1_a'].append([-1,-1,4])
    # twin_systems[twintype]['eta2_a'].append([0,0,1])
    # twin_systems[twintype]['K2_a'].append([1,1,0])
    
    twin_systems[twintype]['eta1_a'].append([1,2,2])
    twin_systems[twintype]['K1_a'].append([-4,1,1])
    twin_systems[twintype]['eta2_a'].append([-1,0,0])
    twin_systems[twintype]['K2_a'].append([0,-1,1])
    uvw2xyz=Parent_uvw2xyz
    hkl2xyz=Parent_hkl2xyz
    a1=uvw2xyz.dot(twin_systems[twintype]['eta1_a'][0])
    twin_systems[twintype]['a1_a'].append(a1/np.linalg.norm(a1))
    a2=uvw2xyz.dot(twin_systems[twintype]['eta2_a'][0])
    twin_systems[twintype]['a2_a'].append(a2/np.linalg.norm(a2))
    n1=hkl2xyz.dot(twin_systems[twintype]['K1_a'][0])
    twin_systems[twintype]['n1_a'].append(n1/np.linalg.norm(n1))
    n2=hkl2xyz.dot(twin_systems[twintype]['K2_a'][0])
    twin_systems[twintype]['n2_a'].append(n2/np.linalg.norm(n2))
    
    twin_systems[twintype]['uvw2xyz_a'].append(np.matmul((2*np.outer(twin_systems[twintype]['a1_a'][0],twin_systems[twintype]['a1_a'][0])-np.eye(3)),uvw2xyz))
    
    #twin_systems[twintype]['R_a'].append(twin_systems[twintype]['uvw2xyz_a'][-1].dot(inv(uvw2xyz)))
    #C1_m=-np.eye(3)+2*np.outer(n1_m,n1_m)
    #np.eye(3)-2*np.outer(n1_m,n1_m)
    #twin_systems[twintype]['C_a'].append(np.eye(3)-2*np.outer(twin_systems[twintype]['n1_a'][0],twin_systems[twintype]['n1_a'][0]))
    twin_systems[twintype]['C_a'].append(2*np.outer(twin_systems[twintype]['a1_a'][0],twin_systems[twintype]['a1_a'][0])-np.eye(3))
    
    twin_systems[twintype]['shear_angle'].append(2*abs(np.pi/2-np.arccos((twin_systems[twintype]['uvw2xyz_a'][0].dot(twin_systems[twintype]['eta2_a'][0])
                                         /norm(twin_systems[twintype]['uvw2xyz_a'][0].dot(twin_systems[twintype]['eta2_a'][0]))).dot(
                                             twin_systems[twintype]['a1_a'][0]))))
    twin_systems[twintype]['s'].append(np.tan(twin_systems[twintype]['shear_angle'][0]/2)*2)
    
    n2eta1K1=np.cross(twin_systems[twintype]['a1_a'][-1],twin_systems[twintype]['n1_a'][-1]);
    Qrot=np.vstack((twin_systems[twintype]['a1_a'][-1],n2eta1K1,twin_systems[twintype]['n1_a'][-1]))
    Qrot.dot(twin_systems[twintype]['a1_a'][-1])
    Qri=np.linalg.inv(Qrot)
    StrainTensor=np.array([[0,0,twin_systems[twintype]['s'][-1]/2],[0,0,0],[twin_systems[twintype]['s'][-1]/2,0,0]])
    twin_systems[twintype]['StrainTensor_a'].append(np.matmul(np.matmul(Qri,StrainTensor),Qrot))
    
    D,V = np.linalg.eig(twin_systems[twintype]['StrainTensor_a'][-1])
    Idxs = np.argsort(D)[::-1]
    Lambda=D[Idxs]
    V = V[:,Idxs]
    twin_systems[twintype]['Tension_a'].append(V[:,0])
    twin_systems[twintype]['Compression_a'].append(V[:,2])
    
    twin_systems[twintype]['dm_a'].append(get_twinning_dislocation(twin_systems[twintype]['K1_a'][-1],
                                                            twin_systems[twintype]['eta1_a'][-1],
                                                            twin_systems[twintype]['eta2_a'][-1],
                                                            L_A,G=G_A,Gr=Gr_A))
    
    
    
    
    
    for symhkl,symuvw in zip(B2_recsymops,B2_symops):
        K1=np.round(symhkl[0:3,0:3].dot(twin_systems[twintype]['K1_a'][0]),5)
        K2=np.round(symhkl[0:3,0:3].dot(twin_systems[twintype]['K2_a'][0]),5)
        eta1=np.round(symuvw[0:3,0:3].dot(twin_systems[twintype]['eta1_a'][0]),5)
        eta2=np.round(symuvw[0:3,0:3].dot(twin_systems[twintype]['eta2_a'][0]),5)
        isin=False
        for i in range(0,len(twin_systems[twintype]['eta1_a'])):
            if (twin_systems[twintype]['eta1_a'][i]==eta1).all():
                if (twin_systems[twintype]['eta2_a'][i]==eta2).all():
                    if (twin_systems[twintype]['K1_a'][i]==K1).all():
                        if (twin_systems[twintype]['K2_a'][i]==K2).all():
                            isin=True  
        if not isin:
            K1in=[1 for K1Mi in twin_systems[twintype]['K1_a'] if (K1 == K1Mi).all() or (-1*K1 == K1Mi).all()]
            E1in=[1 for E1Mi in twin_systems[twintype]['eta1_a'] if (eta1 == E1Mi).all() or (-1*eta1 == E1Mi).all()]
            if len(K1in)==0 or len(E1in)==0:
                twin_systems[twintype]['eta1_a'].append(eta1)
                twin_systems[twintype]['K1_a'].append(K1)
                twin_systems[twintype]['eta2_a'].append(eta2)
                twin_systems[twintype]['K2_a'].append(K2)
                a1=uvw2xyz.dot(twin_systems[twintype]['eta1_a'][-1])
                twin_systems[twintype]['a1_a'].append(a1/np.linalg.norm(a1))
                a2=uvw2xyz.dot(twin_systems[twintype]['eta2_a'][-1])
                twin_systems[twintype]['a2_a'].append(a2/np.linalg.norm(a2))
                n1=uvw2xyz.dot(twin_systems[twintype]['K1_a'][-1])
                twin_systems[twintype]['n1_a'].append(n1/np.linalg.norm(n1))
                n2=uvw2xyz.dot(twin_systems[twintype]['K2_a'][-1])
                twin_systems[twintype]['n2_a'].append(n2/np.linalg.norm(n2))
                
                twin_systems[twintype]['uvw2xyz_a'].append(np.matmul((2*np.outer(twin_systems[twintype]['a1_a'][-1],twin_systems[twintype]['a1_a'][-1])-np.eye(3)),uvw2xyz))
                #twin_systems[twintype]['R_a'].append(twin_systems[twintype]['uvw2xyz_a'][-1].dot(inv(uvw2xyz)))
                #twin_systems[twintype]['C_a'].append(np.eye(3)-2*np.outer(twin_systems[twintype]['n1_a'][-1],twin_systems[twintype]['n1_a'][-1]))
                twin_systems[twintype]['C_a'].append(2*np.outer(twin_systems[twintype]['a1_a'][-1],twin_systems[twintype]['a1_a'][-1])-np.eye(3))
                twin_systems[twintype]['shear_angle'].append(2*abs(np.pi/2-np.arccos((twin_systems[twintype]['uvw2xyz_a'][-1].dot(twin_systems[twintype]['eta2_a'][-1])
                                                     /norm(twin_systems[twintype]['uvw2xyz_a'][-1].dot(twin_systems[twintype]['eta2_a'][-1]))).dot(
                                                         twin_systems[twintype]['a1_a'][-1]))))
                twin_systems[twintype]['s'].append(np.tan(twin_systems[twintype]['shear_angle'][-1]/2)*2)
                
                n2eta1K1=np.cross(twin_systems[twintype]['a1_a'][-1],twin_systems[twintype]['n1_a'][-1]);
                Qrot=np.vstack((twin_systems[twintype]['a1_a'][-1],n2eta1K1,twin_systems[twintype]['n1_a'][-1]))
                Qrot.dot(twin_systems[twintype]['a1_a'][-1])
                Qri=np.linalg.inv(Qrot)
                StrainTensor=np.array([[0,0,twin_systems[twintype]['s'][-1]/2],[0,0,0],[twin_systems[twintype]['s'][-1]/2,0,0]])
                twin_systems[twintype]['StrainTensor_a'].append(np.matmul(np.matmul(Qri,StrainTensor),Qrot))
                
                D,V = np.linalg.eig(twin_systems[twintype]['StrainTensor_a'][-1])
                Idxs = np.argsort(D)[::-1]
                Lambda=D[Idxs]
                V = V[:,Idxs]
                twin_systems[twintype]['Tension_a'].append(V[:,0])
                twin_systems[twintype]['Compression_a'].append(V[:,2])

                twin_systems[twintype]['dm_a'].append(get_twinning_dislocation(twin_systems[twintype]['K1_a'][-1],
                                                                        twin_systems[twintype]['eta1_a'][-1],
                                                                        twin_systems[twintype]['eta2_a'][-1],
                                                                        L_A,G=G_A,Gr=Gr_A))


    #114
    twintypes=['111']
    eta1s_a=[[2,1,1]]
    K1s_a=[[-1,1,1]]
    eta2s_a=[[-1,0,0]]
    K2s_a=[[0,-1,1]]


    twintypes.append('112')
    eta1s_a.append([1,1,1])
    K1s_a.append([-2,1,1])
    eta2s_a.append([-1,0,0])
    K2s_a.append([0,-1,1])
    
    twintypes.append('113')
    eta1s_a.append([2,3,3])
    K1s_a.append([-3,1,1])
    eta2s_a.append([-1,0,0])
    K2s_a.append([0,-1,1])

    twintypes.append('115')
    eta1s_a.append([2,5,5])
    K1s_a.append([-5,1,1])
    eta2s_a.append([-1,0,0])
    K2s_a.append([0,-1,1])

    twintypes.append('116')
    eta1s_a.append([1,3,3])
    K1s_a.append([-6,1,1])
    eta2s_a.append([-1,0,0])
    K2s_a.append([0,-1,1])

    twintypes.append('117')
    eta1s_a.append([2,7,7 ])
    K1s_a.append([-7,1,1])
    eta2s_a.append([-1,0,0])
    K2s_a.append([0,-1,1])

    #twintypes.append('111')
    #eta1s_a.append([-1,-1,-2])
    #K1s_a.append([-1,-1,1])
    #eta2s_a.append([-1,-1,2])
    #K2s_a.append([-1,-1,-1])

    

    for twintype,eta1_a,K1_a,eta2_a,K2_a in zip(twintypes,eta1s_a,K1s_a,eta2s_a,K2s_a):   
        twin_systems[twintype]={}
        for key in list(twindata[0].keys()):
            twin_systems[twintype][key]=[]
        twin_systems[twintype]['vars']=[]
        twin_systems[twintype]['uvw2xyz_m']=[]
        twin_systems[twintype]['uvw2xyz_a']=[]
        twin_systems[twintype]['StrainTensor_a']=[]
        twin_systems[twintype]['Tension_a']=[]
        twin_systems[twintype]['Compression_a']=[]
        twin_systems[twintype]['Variant pairs']=[]
        twin_systems[twintype]['dm_a']=[]
        twin_systems[twintype]['dm_m']=[]
         
        # twin_systems[twintype]['eta1_a'].append([-1,-2,-2])
        # twin_systems[twintype]['K1_a'].append([-1,-1,4])
        # twin_systems[twintype]['eta2_a'].append([0,0,1])
        # twin_systems[twintype]['K2_a'].append([1,1,0])
        
        twin_systems[twintype]['eta1_a'].append(eta1_a)
        twin_systems[twintype]['K1_a'].append(K1_a)
        twin_systems[twintype]['eta2_a'].append(eta2_a)
        twin_systems[twintype]['K2_a'].append(K2_a)
        uvw2xyz=Parent_uvw2xyz
        hkl2xyz=Parent_hkl2xyz
        a1=uvw2xyz.dot(twin_systems[twintype]['eta1_a'][0])
        twin_systems[twintype]['a1_a'].append(a1/np.linalg.norm(a1))
        a2=uvw2xyz.dot(twin_systems[twintype]['eta2_a'][0])
        twin_systems[twintype]['a2_a'].append(a2/np.linalg.norm(a2))
        n1=hkl2xyz.dot(twin_systems[twintype]['K1_a'][0])
        twin_systems[twintype]['n1_a'].append(n1/np.linalg.norm(n1))
        n2=hkl2xyz.dot(twin_systems[twintype]['K2_a'][0])
        twin_systems[twintype]['n2_a'].append(n2/np.linalg.norm(n2))
        
        twin_systems[twintype]['uvw2xyz_a'].append(np.matmul((2*np.outer(twin_systems[twintype]['a1_a'][0],twin_systems[twintype]['a1_a'][0])-np.eye(3)),uvw2xyz))
        #twin_systems[twintype]['R_a'].append(twin_systems[twintype]['uvw2xyz_a'][-1].dot(inv(uvw2xyz)))
        #C1_m=-np.eye(3)+2*np.outer(n1_m,n1_m)
        #np.eye(3)-2*np.outer(n1_m,n1_m)
        #twin_systems[twintype]['C_a'].append(np.eye(3)-2*np.outer(twin_systems[twintype]['n1_a'][0],twin_systems[twintype]['n1_a'][0]))
        twin_systems[twintype]['C_a'].append(2*np.outer(twin_systems[twintype]['a1_a'][0],twin_systems[twintype]['a1_a'][0])-np.eye(3))
        twin_systems[twintype]['shear_angle'].append(2*abs(np.pi/2-np.arccos((twin_systems[twintype]['uvw2xyz_a'][0].dot(twin_systems[twintype]['eta2_a'][0])
                                             /norm(twin_systems[twintype]['uvw2xyz_a'][0].dot(twin_systems[twintype]['eta2_a'][0]))).dot(
                                                 twin_systems[twintype]['a1_a'][0]))))
        twin_systems[twintype]['s'].append(np.tan(twin_systems[twintype]['shear_angle'][0]/2)*2)
        
        n2eta1K1=np.cross(twin_systems[twintype]['a1_a'][-1],twin_systems[twintype]['n1_a'][-1]);
        Qrot=np.vstack((twin_systems[twintype]['a1_a'][-1],n2eta1K1,twin_systems[twintype]['n1_a'][-1]))
        Qrot.dot(twin_systems[twintype]['a1_a'][-1])
        Qri=np.linalg.inv(Qrot)
        StrainTensor=np.array([[0,0,twin_systems[twintype]['s'][-1]/2],[0,0,0],[twin_systems[twintype]['s'][-1]/2,0,0]])
        twin_systems[twintype]['StrainTensor_a'].append(np.matmul(np.matmul(Qri,StrainTensor),Qrot))
        
        D,V = np.linalg.eig(twin_systems[twintype]['StrainTensor_a'][-1])
        Idxs = np.argsort(D)[::-1]
        Lambda=D[Idxs]
        V = V[:,Idxs]
        twin_systems[twintype]['Tension_a'].append(V[:,0])
        twin_systems[twintype]['Compression_a'].append(V[:,2])
        
        
        twin_systems[twintype]['dm_a'].append(get_twinning_dislocation(twin_systems[twintype]['K1_a'][-1],
                                                                twin_systems[twintype]['eta1_a'][-1],
                                                                twin_systems[twintype]['eta2_a'][-1],
                                                                L_A,G=G_A,Gr=Gr_A))
            
        #print(twin_systems[twintype]['dm_a'])
        #print(twin_systems[twintype]['eta2_a'][-1])
        #print(G_A)
        
        for symhkl,symuvw in zip( B2_recsymops,B2_symops):
            K1=np.round(symhkl[0:3,0:3].dot(twin_systems[twintype]['K1_a'][0]),5)
            K2=np.round(symhkl[0:3,0:3].dot(twin_systems[twintype]['K2_a'][0]),5)
            eta1=np.round(symuvw[0:3,0:3].dot(twin_systems[twintype]['eta1_a'][0]),5)
            eta2=np.round(symuvw[0:3,0:3].dot(twin_systems[twintype]['eta2_a'][0]),5)
            isin=False
            for i in range(0,len(twin_systems[twintype]['eta1_a'])):
                if (twin_systems[twintype]['eta1_a'][i]==eta1).all():
                    if (twin_systems[twintype]['eta2_a'][i]==eta2).all():
                        if (twin_systems[twintype]['K1_a'][i]==K1).all():
                            if (twin_systems[twintype]['K2_a'][i]==K2).all():
                                isin=True  
            if not isin:
                K1in=[1 for K1Mi in twin_systems[twintype]['K1_a'] if (K1 == K1Mi).all() or (-1*K1 == K1Mi).all()]
                E1in=[1 for E1Mi in twin_systems[twintype]['eta1_a'] if (eta1 == E1Mi).all() or (-1*eta1 == E1Mi).all()]
                if len(K1in)==0 or len(E1in)==0:
                    twin_systems[twintype]['eta1_a'].append(eta1)
                    twin_systems[twintype]['K1_a'].append(K1)
                    twin_systems[twintype]['eta2_a'].append(eta2)
                    twin_systems[twintype]['K2_a'].append(K2)
                    a1=uvw2xyz.dot(twin_systems[twintype]['eta1_a'][-1])
                    twin_systems[twintype]['a1_a'].append(a1/np.linalg.norm(a1))
                    a2=uvw2xyz.dot(twin_systems[twintype]['eta2_a'][-1])
                    twin_systems[twintype]['a2_a'].append(a2/np.linalg.norm(a2))
                    n1=uvw2xyz.dot(twin_systems[twintype]['K1_a'][-1])
                    twin_systems[twintype]['n1_a'].append(n1/np.linalg.norm(n1))
                    n2=uvw2xyz.dot(twin_systems[twintype]['K2_a'][-1])
                    twin_systems[twintype]['n2_a'].append(n2/np.linalg.norm(n2))
                    
                    twin_systems[twintype]['uvw2xyz_a'].append(np.matmul((2*np.outer(twin_systems[twintype]['a1_a'][-1],twin_systems[twintype]['a1_a'][-1])-np.eye(3)),uvw2xyz))
                    #twin_systems[twintype]['R_a'].append(twin_systems[twintype]['uvw2xyz_a'][-1].dot(inv(uvw2xyz)))
                    #twin_systems[twintype]['C_a'].append(np.eye(3)-2*np.outer(twin_systems[twintype]['n1_a'][-1],twin_systems[twintype]['n1_a'][-1]))
                    twin_systems[twintype]['C_a'].append(2*np.outer(twin_systems[twintype]['a1_a'][-1],twin_systems[twintype]['a1_a'][-1])-np.eye(3))
                    twin_systems[twintype]['shear_angle'].append(2*abs(np.pi/2-np.arccos((twin_systems[twintype]['uvw2xyz_a'][-1].dot(twin_systems[twintype]['eta2_a'][-1])
                                                         /norm(twin_systems[twintype]['uvw2xyz_a'][-1].dot(twin_systems[twintype]['eta2_a'][-1]))).dot(
                                                             twin_systems[twintype]['a1_a'][-1]))))
                    twin_systems[twintype]['s'].append(np.tan(twin_systems[twintype]['shear_angle'][-1]/2)*2)
                    
                    n2eta1K1=np.cross(twin_systems[twintype]['a1_a'][-1],twin_systems[twintype]['n1_a'][-1]);
                    Qrot=np.vstack((twin_systems[twintype]['a1_a'][-1],n2eta1K1,twin_systems[twintype]['n1_a'][-1]))
                    Qrot.dot(twin_systems[twintype]['a1_a'][-1])
                    Qri=np.linalg.inv(Qrot)
                    StrainTensor=np.array([[0,0,twin_systems[twintype]['s'][-1]/2],[0,0,0],[twin_systems[twintype]['s'][-1]/2,0,0]])
                    twin_systems[twintype]['StrainTensor_a'].append(np.matmul(np.matmul(Qri,StrainTensor),Qrot))
                    
                    D,V = np.linalg.eig(twin_systems[twintype]['StrainTensor_a'][-1])
                    Idxs = np.argsort(D)[::-1]
                    Lambda=D[Idxs]
                    V = V[:,Idxs]
                    twin_systems[twintype]['Tension_a'].append(V[:,0])
                    twin_systems[twintype]['Compression_a'].append(V[:,2])
    
                    twin_systems[twintype]['dm_a'].append(get_twinning_dislocation(twin_systems[twintype]['K1_a'][-1],
                                                                            twin_systems[twintype]['eta1_a'][-1],
                                                                            twin_systems[twintype]['eta2_a'][-1],
                                                                            L_A,G=G_A,Gr=Gr_A))
    
    

    for twintype in twin_systems.keys():
        for key in twin_systems[twintype]:
            if type(twin_systems[twintype][key])==list:
                if len(twin_systems[twintype][key])>0:
                    if len(twin_systems[twintype][key])==1:
                        if type(twin_systems[twintype][key][0])==list:
                            #print('reducing dimension for ',twintype,key)
                            twin_systems[twintype][key] = twin_systems[twintype][key][0]
                            #print(twin_systems[twintype][key])
    #print(twin_systems['001']['eta1_m'])        
    return twin_systems           


def get_twinningdata(orim,eus,Ldir_css,twin_systems,twt,phase, tension=True):

    """
    Extract complete twinning data from parent and twin lattices.
    
    Analyzes lattice relationship to determine all twinning elements.
    
    Input:
        L_parent (array 3×3): Parent lattice matrix
        L_twin (array 3×3): Twin lattice matrix
    
    Output:
        dict: Complete twin characterization including K1, K2, η1, η2, S
    
    Usage Example:
        >>> L_parent = monoclinic_lattice_vec(2.889, 4.120, 4.622, 96.8)
        >>> # Twin by compound operation
        >>> R_twin = rotation_from_axis_angle([1,1,1], np.radians(70.5))
        >>> L_twin = R_twin.dot(L_parent)
        >>> 
        >>> data = get_twinningdata(L_parent, L_twin)
        >>> print("Twinning elements:")
        >>> for key, val in data.items():
        ...     print(f"  {key}: {val}")
    
    Notes:
        - Comprehensive twin analysis
        - Extracts all crystallographic elements
        - Validates twin relationship
        - Used in texture simulations
    """
    #inputs:
    #orim: list of orientations as list of orientation matrices csc=orim[gi]*css (css - coord.sys of sample, csc - coord sys of lattice)
    #eus - list of Euler angles 
    #Ldir_css - unit vector of loading
    #twin_systems: dictionary of twinning systems
    #twt: string with type of twinning to be considered - will be used as twin_systems[twt]
    #phase: 'a' for austenite twinning system, 'm' for martensite twinning system
    #outputs:
    #Dictionary containing under key ['n_csl'] list of the normal of the twinning system with higesthest probability to by activated
    #in individual grains (in the coordinate system of the lattice,
    #under key ['n_css'] list of normals the coordinate system of the sample,
    #under key ['angle_n_css'] list of angles between normals and loading direction
    #and under key ['SF'] list of Schmid factors ranging between -0.5 to 0.5, i.e. propensity to twinning
    #and under key ['twsimax'] index of twinning system with highest SF
    #and under key ['neworim'] list of orientations updated by rotation due to twinning
    #and under key ['neweus'] euler angles of updated orientations
    TwinnigData={}
    TwinnigData['n_csl']=[]
    TwinnigData['n_css']=[]
    TwinnigData['angle_n_css']=[]
    TwinnigData['SF']=[]
    TwinnigData['neworim']=[]
    TwinnigData['neweus']=[]
    TwinnigData['eus']=[]
    TwinnigData['twsimax']=[]
    TwinnigData['n1']=[]
    TwinnigData['a1']=[]
    titles='StrainNonSym DefGrad StrainSymEng StrainSymGL DefGradTwin StrainSymEngTwin StrainSymGLTwin StrainLdirSymGl StrainLdirSymEng'.split()
    for title in titles:
        TwinnigData[title]=[]
    #grain index
    gi=0
    #Unit vector of the loading direction in coordinate system of the sample
    #Ldir_css=np.array([0,0,1])
    for gi in range(0,len(orim)):
        # Loading direction in coordinate system of the lattice
        Ldir_csl=orim[gi].dot(Ldir_css)
        SFgi=[]
        #index over all the twinning systems within twt family and get all their Schmid factors
        for twsi in range(0,len(twin_systems[twt]['n1_'+phase])):
            #Unit vector of the normal to twinning plane (correspond to K1)
            n1=twin_systems[twt]['n1_'+phase][twsi]
            #Unit vector of the shear direction (correspond to eta1)
            a1=twin_systems[twt]['a1_'+phase][twsi]
            #Get propensity to twinning
            #SFgi.append(np.sign(a1.dot(Ldir_csl))*n1.dot(Ldir_csl)*a1.dot(Ldir_csl))
            SFgi.append(n1.dot(Ldir_csl)*a1.dot(Ldir_csl))

        SFgi=np.array(SFgi)
        if tension:
            twsimax=np.argmax(SFgi)
        else:
            twsimax=np.argmax(SFgi)
        TwinnigData['n_csl'].append(twin_systems[twt]['n1_'+phase][twsimax])
        TwinnigData['n_css'].append(orim[gi].T.dot(twin_systems[twt]['n1_'+phase][twsimax]))
        TwinnigData['angle_n_css'].append(np.arccos(TwinnigData['n_csl'][-1].dot(Ldir_csl))*180/np.pi)
        TwinnigData['SF'].append(SFgi[twsimax])
        TwinnigData['neworim'].append(twin_systems[twt]['C_'+phase][twsimax].dot(orim[gi]))
        TwinnigData['neweus'].append(euler_angles_from_matrix(TwinnigData['neworim'][-1]))
        TwinnigData['eus'].append(eus[gi])
        TwinnigData['twsimax'].append(twsimax)
        TwinnigData['n1'].append(twin_systems[twt]['n1_'+phase][twsimax])
        TwinnigData['a1'].append(twin_systems[twt]['a1_'+phase][twsimax])
        a1=twin_systems[twt]['a1_'+phase][twsimax]
        n1=twin_systems[twt]['n1_'+phase][twsimax]
        s=twin_systems[twt]['s'][twsimax]
        z=np.cross(a1,n1)
        #transformation matrix into system aligned with a1, n1
        T=np.array([n1,a1,z])
        #transformation matrix into system aligned with basal directions of the twinned lattice
        Tij=2*np.outer(a1,a1)-np.eye(3)
        #Strains and def. gradient in the system aligned with basal directions of the matrix
        StrainNonSym=s*np.outer(a1,n1)
        DefGrad=np.eye(3)+StrainNonSym
        StrainSymEng=0.5*(StrainNonSym+StrainNonSym.T)#0.5*(DefGrad+DefGrad.T) - np.eye(3)
        #Green Lagrange takes into account the elongation of the unit side [010] of the reference square that is sheared along [100] by s
        #so that [010] becomes [s10] (relatively longer (sqrt(s**2+1)) compared to [010])
        StrainSymGL=0.5*(DefGrad.T.dot(DefGrad)-np.eye(3))
        StrainLdirSymGl=np.array(Ldir_csl).dot(StrainSymGL.dot(Ldir_csl))
        StrainLdirSymEng=np.array(Ldir_csl).dot(StrainSymEng.dot(Ldir_csl))
        #Strains and def. gradient in the system aligned with basal directions of the twin
        DefGradTwin=Tij.dot(DefGrad.dot(Tij.T))
        StrainSymEngTwin=0.5*(DefGradTwin+DefGradTwin.T) - np.eye(3)
        StrainSymGLTwin=0.5*(DefGradTwin.T.dot(DefGradTwin)-np.eye(3))
        for title in titles:
            exec(f'TwinnigData["{title}"].append({title})')


    return TwinnigData



def get_twinning_dislocation(K1,eta1,eta2,L,G=None,Gr=None):

    """
    Calculate twinning dislocation Burgers vector.
    
    Determines dislocation content at twin boundary.
    
    Input:
        twin_data (dict): From get_twinningdata()
        b_parent (array [3]): Burgers vector in parent
    
    Output:
        dict: {
            'b_residual': array [3] - residual Burgers vector
            'line_direction': array [3] - dislocation line direction
            'type': str - 'edge', 'screw', or 'mixed'
        }
    
    Usage Example:
        >>> twin_data = niti_twinning('type1')
        >>> b_parent = np.array([1., 0., 0.])  # Parent Burgers vector
        >>> 
        >>> disloc = get_twinning_dislocation(twin_data, b_parent)
        >>> print(f"Residual Burgers vector: {disloc['b_residual']}")
        >>> print(f"Dislocation type: {disloc['type']}")
    
    Notes:
        - Important for twin boundary structure
        - Affects boundary mobility
        - Related to interfacial energy
    """
    #Checked for cubic lattice and K1=-1-11,eta1=-1-1-2,eta2=-1-12 =1/6
    #https://www.sciencedirect.com/science/article/pii/S0966979514000892?via%3Dihub
    #L -latice tensor converting uvw->xyz
    #if type(K1) is list:
        #K1=K1[0]   
    #    print(K1)        
    #if type(eta1) is list:
    #    eta1=eta1[0]
    #if type(eta2) is list:
    #    eta2=eta2[0]
    #print(K1)
    #print(eta1)
    #print(eta2)
    if True:
        if G is None:
            #Meteric tensor ||[x,y,x]||^2=[uvw]*G_A*[uvw]^T
            G = np.matmul(L.T,L)
        if Gr is None:
            #Reciprocal Meteric tensor ||[x,y,x]||^2=d_hkl^2=1/[hkl]*Gr_A*[hkl]^T
            Gr = np.linalg.inv(G)   
        if type(K1) is not np.ndarray:
            K1=np.array(K1)
        if type(eta1) is not np.ndarray:
            eta1=np.array(eta1)
        if type(eta2) is not np.ndarray:
            eta2=np.array(eta2)  
        #rotation angle of the twinning = shear is defined as 2* tan(shear angle/2)
        shear_angle=2*abs(np.pi/2-np.arccos((L.dot(eta2)/np.linalg.norm(L.dot(eta2))).dot(L.dot(eta1)/np.linalg.norm(L.dot(eta1)))))
        #print(shear_angle
        #d-spacing of twinning plane K1
        dK1=1/np.sqrt(K1.dot(Gr).dot(K1))
        #length of eta1
        deta1=np.sqrt(eta1.dot(G).dot(eta1))
        #slip of the firs K1 plane above the twinning plane
        s=2*np.tan(shear_angle/2)*dK1
        #twinning dislocation magnitude as fraction of eta1 (dm*eta1)
        dm=s/deta1
    return dm


def gen_twinned_lattice_points(ParentLatticePoints,eta1,shear_angle,K1,shift=0.0,dK1=None,bvr=None,deta1=None):

    """
    Generate lattice points for parent and twinned regions.
    
    Creates atomic positions showing twin boundary.
    
    Input:
        L_parent (array 3×3): Parent lattice
        twin_plane_normal (array [3]): Twin plane (K1)
        n1, n2, n3 (int): Cell range
    
    Output:
        tuple: (parent_atoms, twin_atoms) - atomic positions
    
    Usage Example:
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
    
    Notes:
        - Shows twin boundary clearly
        - Used for visualization
        - Can calculate boundary energy
    """
    #bvr
    #bvr - burgers vector length in terms of fraction of eta1 direction - deta1
    #dK1 - d-spacing of K1
    #deta1 - length of eta1 direction
    TwinnedPoints=[]
    for Points in ParentLatticePoints:
        twpoints=[]
        for points in Points:
            #Spoints=copy.deepcopy(points)*0
            #for idx in range(3):
                #Spoints[idx,:]=shift[idx]
            if bvr is None:
                twpoints.append(points+2*np.array([eta1]).T*np.tan(shear_angle/2)*(K1.dot(points)+shift))
            else:
                twpoints.append(points+bvr*deta1*np.array([eta1]).T*(np.modf(((K1.dot(points)+shift*np.sign(K1.dot(points)))/dK1))[1]))
        
        TwinnedPoints.append(twpoints)
    return TwinnedPoints


def write_txt(filename,Header,DATA):

    """
    Write numerical data to text file.
    
    Saves array data in formatted text file with optional header.
    
    Input:
        filename (str): Output file path
        data (array): Numerical data to save
        header (str): Optional header line
    
    Output:
        None (writes to file)
    
    Usage Example:
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
    
    Notes:
        - Standard text format
        - Can be read by most software
        - Header optional but recommended
        - Uses space-separated values
    """
    f=open(filename, "w")
    f.write('\t'.join(Header)+'\n')
    
    inc=0;
    for data in DATA:
        inc+=1
        print(str(inc)+'/'+str(len(DATA)))
        formatted = ['%.8f'% item  for item in data]
        f.write('\t'.join(formatted)+'\n')
    
    f.close()


def read_txt(filename,delimiter='\t',skiprows=1):

    """
    Read numerical data from text file.
    
    Loads data written by write_txt() or similar format.
    
    Input:
        filename (str): Input file path
    
    Output:
        numpy.ndarray: Loaded data
    
    Usage Example:
        >>> # Read previously saved data
        >>> data = read_txt('orientations.txt')
        >>> print(f"Loaded {len(data)} orientations")
        >>> print(f"Shape: {data.shape}")
        >>> 
        >>> # Use data
        >>> euler_angles = data[:, :3]
        >>> # Process...
    
    Notes:
        - Skips comment lines starting with #
        - Returns numpy array
        - Compatible with write_txt output
    """
    f=open(filename, "r")
    Header=f.readline().strip().split(delimiter)
    f.close()
    Data = np.loadtxt(filename,  skiprows=skiprows, delimiter=delimiter);
    return Header,Data


def plane_line_intersection(n,V0,P0,P1):
    # n: normal vector of the Plane 
    # V0: any point that belongs to the Plane 
    # P0: end point 1 of the segment P0P1
    # P1:  end point 2 of the segment P0P1


    w = P0 - V0;
    u = P1-P0;
    N = -np.dot(n,w);
    D = np.dot(n,u)
    sI = N / D
    I = P0+ sI*u
    return I



def flipvector(v, Tol=1e-9):

    """
    Flip vector to ensure positive first non-zero component.
    
    Standardizes vector direction for comparison.
    
    Input:
        v (array [3]): Vector
    
    Output:
        numpy.ndarray [3]: Flipped vector (if needed)
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> v1 = np.array([-1, 2, 3])
        >>> v1_flipped = flipvector(v1)
        >>> print(v1_flipped)  # [1, -2, -3]
        >>> 
        >>> v2 = np.array([0, -1, 2])
        >>> v2_flipped = flipvector(v2)
        >>> print(v2_flipped)  # [0, 1, -2]
    
    Notes:
        - Makes first non-zero component positive
        - Useful for standardization
        - Preserves magnitude and line
    """
    vm = np.round(v/Tol)*Tol
    c=1.0
#    if len(vm[vm<0])==len(vm[vm<>0]):
    if len(vm[vm<0])==len(vm[vm!=0]):
        vm=-1*vm
        c=-1.0
    return vm,c


def flipvector2negative(v, Tol=1e-9):

    """
    Flip vector to ensure negative first non-zero component.
    
    Opposite of flipvector() - ensures negative first component.
    
    Input:
        v (array [3]): Vector
    
    Output:
        numpy.ndarray [3]: Flipped vector (if needed)
    
    Usage Example:
        >>> v = np.array([1, 2, 3])
        >>> v_neg = flipvector2negative(v)
        >>> print(v_neg)  # [-1, -2, -3]
    
    Notes:
        - Makes first non-zero component negative
        - Complementary to flipvector()
    """
    vm = np.round(v/Tol)*Tol
    c=1.0
#    if len(vm[vm<0])==len(vm[vm<>0]):
    if len(vm[vm>0])==len(vm[vm!=0]):
        vm=-1*vm
        c=-1.0
    return vm,c
  

def vector2miller_ini(v, MIN=True, Tol=1e-9,tol=1e5,text=False,decimals=3):

    """
    Convert Cartesian vector to Miller indices (initial guess).
    
    Provides starting point for vector→Miller conversion.
    
    Input:
        v (array [3]): Cartesian vector
        L (array 3×3): Lattice matrix
    
    Output:
        numpy.ndarray [3]: Initial Miller index guess
    
    Usage Example:
        >>> L = cubic_lattice_vec(3.0)
        >>> v_cart = np.array([3, 3, 0])
        >>> uvw_guess = vector2miller_ini(v_cart, L)
        >>> print(uvw_guess)  # Initial guess for [110]
    
    Notes:
        - Provides rough approximation
        - Use vectors2miller for refined result
        - Useful as starting point
    """
    vm = np.round(v/Tol)*Tol
    #print('ok')
    #print((np.round(vm)==vm).all())
    #if (vm==np.array([ 2. ,-1. , 1.])).all():
        #print((np.round(vm)==vm).all())
        #print(vm)
    if (np.abs(vm)<=1).all():
        vm=np.round(vm/abs(min(vm[np.abs(vm)>1/tol]))*tol)/tol
        vm/=min(abs(vm[np.nonzero(vm)[0]]))   
    if not (np.round(vm)==vm).all():
        if MIN:
            vm=np.round(vm/abs(min(vm[np.abs(vm)>1/tol]))*tol)/tol
            #print(vm)
            vm/=min(abs(vm[np.nonzero(vm)[0]]))
            
            if not (np.round(vm)==vm).all():
                vm=np.round(vm/abs(max(vm[np.abs(vm)>1/tol]))*tol)/tol                
                vm/=max(abs(vm[np.nonzero(vm)[0]]))
        else:
            #print(vm)
            vm=np.round(vm/abs(max(vm[np.abs(vm)>1/tol]))*tol)/tol
            
            vm/=max(abs(vm[np.nonzero(vm)[0]]))
            if not (np.round(vm)==vm).all():
                vm=np.round(vm/abs(min(vm[np.abs(vm)>1/tol]))*tol)/tol
                #print(vm)
                vm/=min(abs(vm[np.nonzero(vm)[0]]))
                
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
 

    


"""
CRYSTLIB.PY - Crystal Structure and Lattice Operations Module

This module provides comprehensive functions for crystallographic calculations
and lattice operations in materials science research, with focus on crystal
structures, lattice transformations, and phase relationships.

Main Categories:
================

1. Lattice Vector Generation
   - cubic_lattice_vec, monoclinic_lattice_vec, tetragonal_lattice_vec
   - lattice_vec (general), reciprocal_basis

2. Miller Index Operations
   - uvtw2uvw, uvw2uvtw (4-index ↔ 3-index conversions)
   - hkil2hkl, hkl2hkil (hexagonal notations)
   - miller2fractional, xyz2fractional (coordinate transformations)
   - vector2miller, vectors2miller (index determination)

3. Lattice Correspondence
   - B19p_B2_lattice_correspondence (NiTi transformations)
   - cubic2tetragonal_lattice_correspondence
   - Rp_B2_lattice_correspondence (R-phase)
   - lattice_correspondence (general)

4. Deformation & Strain Analysis
   - def_gradient, def_gradient_stressfree
   - mohr_circles, zero_normal_strains
   - strains_along_13mohrcirle

5. Twinning Analysis
   - niti_twinning, get_twinningdata
   - habitplane_equation_solution
   - twin_equation_solution
   - get_twinning_plane_points

6. Lattice Point Generation
   - generate_lattice_points, generate_lattice_faces
   - generate_product_lattice_points
   - generate_lattite_atom_positions
   - gen_twinned_lattice_points

7. Crystal Plane Operations
   - select_plane, select_atomic_plane
   - select_crystal_planes
   - plane_line_intersection
   - generate_plane_vertices

8. Hexagonal Systems
   - gensystemsHex, gensystemsHexIni
   - genallHexSys

9. Tensor Operations
   - permut_tensor3, np_permut_tensor3 (Levi-Civita)
   - kronecker, np_kronecker (Kronecker delta)

10. Utility Functions
    - find_gcd, perpendicular_vector
    - vec2string, plane2string, dir2string
    - flipvector, normArrayColumns
    - write_txt, read_txt

Key Applications:
- Shape memory alloy analysis (NiTi, Cu-based)
- Martensitic phase transformations
- Crystallographic texture analysis
- Deformation mechanics and twinning
- EBSD data processing

Author: lheller
Created: 2019-2023
"""

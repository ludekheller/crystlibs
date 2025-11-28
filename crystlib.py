#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Crystallography Library Module

This module provides utility functions for crystallographic calculations including:
- Miller indices generation and manipulation
- Lattice vector calculations for various crystal systems
- Reciprocal basis vector calculations

Created on Mon Sep  7 09:06:48 2020
@author: lheller
"""
import numpy as np


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
    if lattice_param['type'].lower() == 'cubic':
        a = lattice_param['a']
        V = a * np.eye(3)
        
    elif lattice_param['type'].lower() == 'tetragonal':
        a = lattice_param['a']
        b = lattice_param['b']
        c = lattice_param['c']
        V = np.zeros((3, 3))
        V[:, 0] = np.array([a, 0., 0])
        V[:, 1] = np.array([0, b, 0])
        V[:, 2] = np.array([0, 0, c])
        
    elif lattice_param['type'].lower() == 'monoclinic':
        a = lattice_param['a']
        b = lattice_param['b']
        c = lattice_param['c']
        beta = lattice_param['beta']
        V = np.zeros((3, 3))
        V[:, 0] = np.array([a, 0., 0])
        V[:, 1] = np.array([0, b, 0])
        V[:, 2] = np.array([c*np.cos(beta), 0, c*np.sin(beta)])
        
    elif lattice_param['type'].lower() == 'triclinic':
        a = lattice_param['a']
        b = lattice_param['b']
        c = lattice_param['c']
        alpha = lattice_param['alpha']
        beta = lattice_param['beta']
        gamma = lattice_param['gamma']
        V = np.zeros((3, 3))
        V[:, 0] = np.array([a, 0., 0])
        V[:, 1] = np.array([b*np.cos(gamma), b*np.sin(gamma), 0])
        cx = c * np.cos(beta)
        cy = c * (np.cos(alpha) - np.cos(beta)*np.cos(gamma)) / np.sin(gamma)
        cz = np.sqrt(c**2 - cx**2 - cy**2)
        V[:, 2] = np.array([cx, cy, cz])
        
    elif lattice_param['type'].lower() == 'trigonal':
        a = lattice_param['a']
        c = lattice_param['c']
        V = np.zeros((3, 3))
        V[:, 0] = np.array([1./2.*a, -np.sqrt(3)/2.*a, 0])
        V[:, 1] = np.array([1./2.*a, np.sqrt(3)/2.*a, 0])
        V[:, 2] = np.array([0, 0, c])
        
    return V[:, 0], V[:, 1], V[:, 2]


def reciprocal_basis(a1, a2, a3):
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
    b1 = np.cross(a2, a3) / np.dot(a1, np.cross(a2, a3))
    b2 = np.cross(a3, a1) / np.dot(a2, np.cross(a3, a1))
    b3 = np.cross(a1, a2) / np.dot(a3, np.cross(a1, a2))
    
    return b1, b2, b3

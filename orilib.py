#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Orientation Library Module

This module provides comprehensive utilities for crystallographic orientation analysis including:
- Quaternion operations (conversion, multiplication, misorientation calculations)
- Euler angle and rotation matrix conversions
- Rodrigues-Frank vector representations
- Orientation sampling and gridding (HEALPix, Hopf coordinates)
- Crystal symmetry operations
- Axis-angle representations

Created on Wed Sep 11 10:07:07 2019
@author: lheller
"""
import numpy as np
import math
from numba import njit
import numba as nb


# =====================================
# Quaternion utilities
# =====================================

@njit

def mat_to_quat(R):
    """
    Convert a rotation matrix to a quaternion representation.
    Uses the Shepperd's method for numerical stability.
    
    Input:
        R: numpy array (3x3) - Rotation matrix (orthogonal with det=1)
    
    Output:
        q: numpy array (4,) - Unit quaternion [w, x, y, z] where w is the scalar part
    
    Usage Example:
        >>> import numpy as np
        >>> # 90-degree rotation around z-axis
        >>> R = np.array([[0, -1, 0],
        ...               [1,  0, 0],
        ...               [0,  0, 1]], dtype=float)
        >>> q = mat_to_quat(R)
        >>> print("Quaternion:", q)
        >>> # Should give approximately [0.707, 0, 0, 0.707]
        
        >>> # Identity rotation
        >>> R_identity = np.eye(3)
        >>> q_identity = mat_to_quat(R_identity)
        >>> print("Identity quaternion:", q_identity)
        >>> # Should give [1, 0, 0, 0]
    """
    tr = R[0,0] + R[1,1] + R[2,2]
    if tr > 0:
        S = np.sqrt(tr + 1.0) * 2
        w = 0.25 * S
        x = (R[2,1] - R[1,2]) / S
        y = (R[0,2] - R[2,0]) / S
        z = (R[1,0] - R[0,1]) / S
    elif (R[0,0] > R[1,1]) and (R[0,0] > R[2,2]):
        S = np.sqrt(1.0 + R[0,0] - R[1,1] - R[2,2]) * 2
        w = (R[2,1] - R[1,2]) / S
        x = 0.25 * S
        y = (R[0,1] + R[1,0]) / S
        z = (R[0,2] + R[2,0]) / S
    elif R[1,1] > R[2,2]:
        S = np.sqrt(1.0 + R[1,1] - R[0,0] - R[2,2]) * 2
        w = (R[0,2] - R[2,0]) / S
        x = (R[0,1] + R[1,0]) / S
        y = 0.25 * S
        z = (R[1,2] + R[2,1]) / S
    else:
        S = np.sqrt(1.0 + R[2,2] - R[0,0] - R[1,1]) * 2
        w = (R[1,0] - R[0,1]) / S
        x = (R[0,2] + R[2,0]) / S
        y = (R[1,2] + R[2,1]) / S
        z = 0.25 * S
    q = np.array([w,x,y,z])
    q /= np.linalg.norm(q)
    return q


@njit

def quat_to_mat(q):
    """
    Convert a quaternion to a rotation matrix.
    
    Input:
        q: numpy array (4,) - Unit quaternion [w, x, y, z]
    
    Output:
        R: numpy array (3x3) - Rotation matrix
    
    Usage Example:
        >>> import numpy as np
        >>> # Quaternion for 90-degree rotation around z-axis
        >>> q = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> R = quat_to_mat(q)
        >>> print("Rotation matrix:")
        >>> print(R)
        >>> # Should give approximately [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        
        >>> # Identity quaternion
        >>> q_identity = np.array([1.0, 0.0, 0.0, 0.0])
        >>> R_identity = quat_to_mat(q_identity)
        >>> print("Identity matrix:")
        >>> print(R_identity)
    """
    w, x, y, z = q
    return np.array([
        [1 - 2*(y*y + z*z), 2*(x*y - z*w), 2*(x*z + y*w)],
        [2*(x*y + z*w), 1 - 2*(x*x + z*z), 2*(y*z - x*w)],
        [2*(x*z - y*w), 2*(y*z + x*w), 1 - 2*(x*x + y*y)]
    ])


@njit

def quat_mult(q1, q2):
    """
    Multiply two quaternions (Hamilton product).
    Order matters: q1 * q2 ≠ q2 * q1 in general.
    
    Input:
        q1: numpy array (4,) - First quaternion [w, x, y, z]
        q2: numpy array (4,) - Second quaternion [w, x, y, z]
    
    Output:
        q: numpy array (4,) - Product quaternion q1 * q2
    
    Usage Example:
        >>> import numpy as np
        >>> # Two 90-degree rotations around z-axis = 180-degree rotation
        >>> q_90 = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> q_result = quat_mult(q_90, q_90)
        >>> print("Result quaternion:", q_result)
        >>> # Should give approximately [0, 0, 0, 1] (180-degree rotation)
        
        >>> # Rotation composition
        >>> q_x = np.array([np.cos(np.pi/4), np.sin(np.pi/4), 0, 0])  # 90° around x
        >>> q_z = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])  # 90° around z
        >>> q_combined = quat_mult(q_x, q_z)
        >>> print("Combined rotation:", q_combined)
    """
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])


@njit

def quat_misori_deg(q1, q2):
    """
    Calculate the misorientation angle between two quaternions in degrees.
    Returns the minimum rotation angle needed to go from q1 to q2.
    
    Input:
        q1: numpy array (4,) - First quaternion [w, x, y, z]
        q2: numpy array (4,) - Second quaternion [w, x, y, z]
    
    Output:
        angle: float - Misorientation angle in degrees (0 to 180)
    
    Usage Example:
        >>> import numpy as np
        >>> # Same orientation
        >>> q1 = np.array([1.0, 0.0, 0.0, 0.0])
        >>> angle = quat_misori_deg(q1, q1)
        >>> print("Same orientation angle:", angle)
        >>> # Should give 0 degrees
        
        >>> # 90-degree misorientation
        >>> q1 = np.array([1.0, 0.0, 0.0, 0.0])
        >>> q2 = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> angle = quat_misori_deg(q1, q2)
        >>> print("Misorientation:", angle, "degrees")
        >>> # Should give approximately 90 degrees
    """
    dq = abs(q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2] + q1[3]*q2[3])
    if dq > 1.0:
        dq = 1.0
    elif dq < -1.0:
        dq = -1.0
    return 2.0 * np.degrees(np.arccos(dq))


@njit(parallel=True, fastmath=True)

def misori_sym_deg_sample_to_crystal_fast(M1, M2, symops):
    """
    Compute misorientation angles (deg) between orientations M1 and M2
    considering crystal symmetry operations using vectorized operations.
    This is the fast version optimized with Numba parallel processing.
    
    Input:
        M1: numpy array (N, 3, 3) - Orientation matrices (sample→crystal reference frame)
        M2: numpy array (N, 3, 3) - Orientation matrices (sample→crystal reference frame)
                                     Must have same length as M1
        symops: numpy array (Ns, 3, 3) - Crystal symmetry operation matrices
    
    Output:
        miso: numpy array (N,) - Minimum misorientation angle in degrees for each pair
    
    Usage Example:
        >>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R
        >>> 
        >>> # Generate some random orientations
        >>> N = 100
        >>> M1 = R.random(N).as_matrix()
        >>> M2 = R.random(N).as_matrix()
        >>> 
        >>> # Cubic symmetry operations (just identity for example)
        >>> symops = np.array([np.eye(3)])
        >>> 
        >>> # Calculate misorientations
        >>> misorientations = misori_sym_deg_sample_to_crystal_fast(M1, M2, symops)
        >>> print("Mean misorientation:", np.mean(misorientations), "degrees")
        >>> print("Max misorientation:", np.max(misorientations), "degrees")
    """
    M1 = np.asarray(M1)
    M2 = np.asarray(M2)
    symops = np.asarray(symops)
    N = len(M1)
    Ns = len(symops)

    # Relative rotation between orientations
    # R12 = M2 @ M1.T
    R12 = np.einsum("nij,nkj->nik", M2, M1, optimize=True)  # (N,3,3)

    # Apply all symmetry operations at once
    # S * R12 for all S in symops
    # shape: (Ns,N,3,3)
    R_eq = np.einsum("sij,njk->snik", symops, R12, optimize=True)

    # Convert all to quaternions
    R_eq_flat = R_eq.reshape(Ns * N, 3, 3)
    q = R.from_matrix(R_eq_flat).as_quat().reshape(Ns, N, 4)

    # Compute rotation angle for each symmetry-equivalent
    w = np.clip(np.abs(q[..., 3]), -1.0, 1.0)
    ang = 2 * np.arccos(w)  # radians
    miso = np.degrees(np.min(ang, axis=0))  # (N,)

    return miso


@njit

def quat_conjugate(q):
    """
    Compute the conjugate of a quaternion (inverse for unit quaternions).
    
    Input:
        q: numpy array (4,) - Quaternion [w, x, y, z]
    
    Output:
        q_conj: numpy array (4,) - Conjugate quaternion [w, -x, -y, -z]
    
    Usage Example:
        >>> import numpy as np
        >>> q = np.array([0.707, 0.707, 0.0, 0.0])
        >>> q_conj = quat_conjugate(q)
        >>> print("Original:", q)
        >>> print("Conjugate:", q_conj)
        >>> # Conjugate: [0.707, -0.707, 0.0, 0.0]
        
        >>> # Verify: q * q_conj = identity
        >>> result = quat_mult(q, q_conj)
        >>> print("q * q_conj:", result)
        >>> # Should give approximately [1, 0, 0, 0]
    """
    return np.array([q[0], -q[1], -q[2], -q[3]])


@njit

def quat_multiply(q1, q2):
    """
    Quaternion multiplication (same as quat_mult, alternative implementation).
    
    Input:
        q1: numpy array (4,) - First quaternion [w, x, y, z]
        q2: numpy array (4,) - Second quaternion [w, x, y, z]
    
    Output:
        q: numpy array (4,) - Product quaternion q1 * q2
    
    Usage Example:
        >>> import numpy as np
        >>> q1 = np.array([1.0, 0.0, 0.0, 0.0])  # Identity
        >>> q2 = np.array([0.707, 0.707, 0.0, 0.0])  # 90° around x
        >>> result = quat_multiply(q1, q2)
        >>> print("Result:", result)
        >>> # Should give q2 since q1 is identity
    """
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1*w2 - x1*x2 - y1*y2 - z1*z2
    x = w1*x2 + x1*w2 + y1*z2 - z1*y2
    y = w1*y2 - x1*z2 + y1*w2 + z1*x2
    z = w1*z2 + x1*y2 - y1*x2 + z1*w2
    return np.array([w,x,y,z])


@njit

def misori_sym_deg_quats(q1, q2, sym_quats):
    """
    Calculate minimum misorientation angle between two quaternions considering
    crystal symmetry operations. This is a Numba-optimized version.
    
    Input:
        q1: numpy array (4,) - First quaternion [w, x, y, z]
        q2: numpy array (4,) - Second quaternion [w, x, y, z]
        sym_quats: numpy array (Ns, 4) - Array of symmetry operation quaternions
    
    Output:
        min_angle: float - Minimum misorientation angle in degrees
    
    Usage Example:
        >>> import numpy as np
        >>> # Define two orientations
        >>> q1 = np.array([1.0, 0.0, 0.0, 0.0])
        >>> q2 = np.array([0.707, 0.0, 0.0, 0.707])  # 90° around z
        >>> 
        >>> # Define cubic symmetry (simplified - just identity for demo)
        >>> sym_quats = np.array([[1.0, 0.0, 0.0, 0.0]])
        >>> 
        >>> angle = misori_sym_deg_quats(q1, q2, sym_quats)
        >>> print("Misorientation with symmetry:", angle, "degrees")
    """
    min_ang = 1e9
    Ns = sym_quats.shape[0]
    for s in range(Ns):
        q = quat_multiply(sym_quats[s], quat_multiply(q2, quat_conjugate(q1)))
        # angle = 2*acos(|w|) in radians
        w = abs(q[0])
        if w > 1.0:
            w = 1.0
        ang = 2.0 * np.arccos(w)
        if ang < min_ang:
            min_ang = ang
    return np.degrees(min_ang)


# =====================================
# Euler Angles and Rotation Matrices
# =====================================


def eu2quat(phi1, Phi, phi2):
    """
    Convert Bunge Euler angles to quaternion representation.
    Euler angles follow the ZXZ convention (Bunge notation).
    
    Input:
        phi1: float - First Euler angle (rotation around Z) in radians
        Phi: float - Second Euler angle (rotation around X') in radians
        phi2: float - Third Euler angle (rotation around Z'') in radians
    
    Output:
        q: numpy array (4,) - Quaternion [w, x, y, z] with positive w
    
    Usage Example:
        >>> import numpy as np
        >>> # 45-degree rotation around each axis
        >>> phi1 = np.radians(45)
        >>> Phi = np.radians(45)
        >>> phi2 = np.radians(45)
        >>> q = eu2quat(phi1, Phi, phi2)
        >>> print("Quaternion:", q)
        
        >>> # Identity rotation (all zeros)
        >>> q_identity = eu2quat(0, 0, 0)
        >>> print("Identity quaternion:", q_identity)
        >>> # Should give [1, 0, 0, 0]
    """
    q = np.zeros((4))
    q[0] = np.cos(Phi/2)*np.cos((phi1+phi2)/2)
    q[1] = -np.sin(Phi/2)*np.cos((phi1-phi2)/2)
    q[2] = -np.sin(Phi/2)*np.sin((phi1-phi2)/2)
    q[3] = -np.cos(Phi/2)*np.sin((phi1+phi2)/2)
    
    if q[0] < 0:
        q = -q
    
    return q



def np_euler_matrix(ai, aj, ak): 

    """
    Convert Euler angles to rotation matrix (ZXZ convention, Bunge notation).
    Single orientation version.
    
    Input:
        ai: float - First Euler angle (phi1) in radians
        aj: float - Second Euler angle (Phi) in radians
        ak: float - Third Euler angle (phi2) in radians
    Output:
        g: numpy array (3, 3) - Rotation matrix
    
    Usage Example:
        >>> import numpy as np
        >>> # 90-degree rotation around Z axis
        >>> phi1 = np.pi/2
        >>> Phi = 0
        >>> phi2 = 0
        >>> g = np_euler_matrix(phi1, Phi, phi2)
        >>> print("Rotation matrix:")
        >>> print(g)
        
        >>> # Identity rotation
        >>> g_identity = np_euler_matrix(0, 0, 0)
        >>> print("Identity:", g_identity)
    """
#    ai=20.*np.pi/180.
#    aj=40.*np.pi/180.
#    ak=80.*np.pi/180.
#    np_euler_matrix(ai, aj, ak)
#    np.matmul(passive_rotation(ak,'z'),np.matmul(passive_rotation(aj,'x'),passive_rotation(ai,'z')))
    
    g=np.eye(3)
    s1, s2, s3 = np.sin(ai), np.sin(aj), np.sin(ak)
    c1, c2, c3 = np.cos(ai), np.cos(aj), np.cos(ak)
    
    g[0,0] = c1*c3-s1*s3*c2
    g[0,1] = s1*c3+c1*s3*c2
    g[0,2] = s3*s2    
    g[1,0] = -c1*s3-s1*c3*c2 
    g[1,1] = -s1*s3+c1*c3*c2
    g[1,2] = c3*s2 
    g[2,0] = s1*s2 
    g[2,1] = -c1*s2
    g[2,2] = c2   
    return g


def np_eulers_matrices(data, deg=False):
    """
    Convert multiple Euler angles to rotation matrices (vectorized).
    Processes an array of Euler angle triplets efficiently.
    
    Input:
        data: numpy array (N, 3) - Array of Euler angles [phi1, Phi, phi2]
        deg: bool - If True, input angles are in degrees; if False, radians (default: False)
    
    Output:
        g: numpy array (N, 3, 3) - Array of rotation matrices
    
    Usage Example:
        >>> import numpy as np
        >>> # Multiple orientations in degrees
        >>> euler_angles = np.array([
        ...     [0, 0, 0],        # Identity
        ...     [90, 0, 0],       # 90° around Z
        ...     [0, 90, 0],       # 90° around X'
        ...     [45, 45, 45]      # Combined rotation
        ... ])
        >>> matrices = np_eulers_matrices(euler_angles, deg=True)
        >>> print("Shape:", matrices.shape)  # (4, 3, 3)
        >>> print("First matrix (identity):")
        >>> print(matrices[0])
        
        >>> # Process many orientations efficiently
        >>> N = 1000
        >>> random_eulers = np.random.rand(N, 3) * np.array([360, 180, 360])
        >>> rot_matrices = np_eulers_matrices(random_eulers, deg=True)
        >>> print(f"Processed {N} orientations")
    """
    # data[0,:] three euler angles
    if deg:
        ai = data[:,0]*np.pi/180
        aj = data[:,1]*np.pi/180
        ak = data[:,2]*np.pi/180
    else:
        ai = data[:,0]
        aj = data[:,1]
        ak = data[:,2]
    
    g = np.zeros((data.shape[0], 3, 3))
    s1, s2, s3 = np.sin(ai), np.sin(aj), np.sin(ak)
    c1, c2, c3 = np.cos(ai), np.cos(aj), np.cos(ak)
    
    g[:,0,0] = c1*c3 - s1*s3*c2
    g[:,0,1] = s1*c3 + c1*s3*c2
    g[:,0,2] = s3*s2    
    g[:,1,0] = -c1*s3 - s1*c3*c2 
    g[:,1,1] = -s1*s3 + c1*c3*c2
    g[:,1,2] = c3*s2 
    g[:,2,0] = s1*s2 
    g[:,2,1] = -c1*s2
    g[:,2,2] = c2       
    return g



def np_inverse_euler_matrix(ai, aj, ak): 

    """
    Convert Euler angles to inverse (transpose) rotation matrix.
    Equivalent to the transpose of the forward rotation matrix.
    
    Input:
        ai: float - First Euler angle (phi1) in radians
        aj: float - Second Euler angle (Phi) in radians
        ak: float - Third Euler angle (phi2) in radians
    
    Output:
        U: numpy array (3, 3) - Inverse rotation matrix (transpose of forward matrix)
    
    Usage Example:
        >>> import numpy as np
        >>> # Forward rotation
        >>> ai, aj, ak = np.pi/4, np.pi/3, np.pi/6
        >>> g_forward = np_euler_matrix(ai, aj, ak)
        >>> g_inverse = np_inverse_euler_matrix(ai, aj, ak)
        >>> 
        >>> # Verify: forward * inverse = identity
        >>> result = g_forward @ g_inverse
        >>> print("Forward @ Inverse:")
        >>> print(result)
        >>> # Should be close to identity matrix
    """
    U=np.eye(3)

    s1, s2, s3 = np.sin(ai), np.sin(aj), np.sin(ak)
    c1, c2, c3 = np.cos(ai), np.cos(aj), np.cos(ak)

    U[0,0] = c1*c3-s1*s3*c2
    U[0,1] = -c1*s3-s1*c3*c2
    U[0,2] = s1*s2    
    U[1,0] = s1*c3+c1*s3*c2 
    U[1,1] = -s1*s3+c1*c3*c2
    U[1,2] = -c1*s2 
    U[2,0] = s3*s2 
    U[2,1] = c3*s2
    U[2,2] = c2   
    return U

def ol_g_rtheta_rad(g):

    """
    Convert rotation matrix to axis-angle representation (Rodrigues-Frank vector).
    Returns rotation axis and angle.
    
    Input:
        g: list or array (3, 3) - Rotation matrix
    
    Output:
        r: list (3,) - Rotation axis (unit vector)
        ptheta: float - Rotation angle in radians
    
    Usage Example:
        >>> import numpy as np
        >>> # 90-degree rotation around Z
        >>> g = [[0, -1, 0],
        ...      [1,  0, 0],
        ...      [0,  0, 1]]
        >>> axis, angle = ol_g_rtheta_rad(g)
        >>> print("Rotation axis:", axis)
        >>> print("Rotation angle (rad):", angle)
        >>> print("Rotation angle (deg):", np.degrees(angle))
        >>> # Should give axis ≈ [0, 0, 1] and angle ≈ π/2
    """
    eps = 1.e-6;
    
    ptheta = np.arccos((g[0][0] + g[1][1] + g[2][2] - 1) / 2);
    r=[0.,0.,0.];
    if ((ptheta) < eps):
        r[0] = 1;
        r[1] = 0;
        r[2] = 0;
    elif ((ptheta) < (1 - eps)*np.pi):
        r[0] = (g[1][2] - g[2][1]) / (2 * np.sin(ptheta));
        r[1] = (g[2][0] - g[0][2]) / (2 * np.sin(ptheta));
        r[2] = (g[0][1] - g[1][0]) / (2 * np.sin(ptheta));
    else:
        r[0] = np.sqrt((g[0][0] + 1) / 2)
        r[1] = np.sqrt((g[1][1] + 1) / 2);
        r[2] = np.sqrt((g[2][2] + 1) / 2);
    m = r.index(max(r))
    for i in range(0,3):
        if not r==m:
            if g[i][m]<0:
                r[i] *= 1;
    return r,ptheta            

def np_ol_g_rtheta_rad(g):

    """
    Convert rotation matrix to axis-angle representation (NumPy optimized version).
    
    Input:
        g: numpy array (3, 3) - Rotation matrix
    
    Output:
        r: numpy array (3,) - Rotation axis (unit vector)
        ptheta: float - Rotation angle in radians
    
    Usage Example:
        >>> import numpy as np
        >>> # 120-degree rotation around [1,1,1]
        >>> angle = np.radians(120)
        >>> axis = np.array([1, 1, 1]) / np.sqrt(3)
        >>> # Create rotation matrix using Rodrigues formula
        >>> K = np.array([[0, -axis[2], axis[1]],
        ...               [axis[2], 0, -axis[0]],
        ...               [-axis[1], axis[0], 0]])
        >>> g = np.eye(3) + np.sin(angle)*K + (1-np.cos(angle))*K@K
        >>> 
        >>> # Convert back to axis-angle
        >>> axis_out, angle_out = np_ol_g_rtheta_rad(g)
        >>> print("Recovered axis:", axis_out)
        >>> print("Recovered angle (deg):", np.degrees(angle_out))
    """
    eps = 1.e-6;
    
    ptheta = np.arccos((np.trace(g) - 1) / 2);
    r=np.array([0.,0.,0.]);
    if ((ptheta) < eps):
        r[0] = 1;
        r[1] = 0;
        r[2] = 0;
    elif ((ptheta) < (1 - eps)*np.pi):
        r[0] = (g[1,2] - g[2,1]) / (2 * np.sin(ptheta));
        r[1] = (g[2,0] - g[0,2]) / (2 * np.sin(ptheta));
        r[2] = (g[0,1] - g[1,0]) / (2 * np.sin(ptheta));
    else:
        r[0] = np.sqrt((g[0,0] + 1) / 2)
        r[1] = np.sqrt((g[1,1] + 1) / 2);
        r[2] = np.sqrt((g[2,2] + 1) / 2);
    m = np.where(r==max(r))[0][0]
    for i in range(0,3):
        if not i==m:
            if g[i,m]<0:
                r[i] *= 1;
    return r,ptheta            


def ol_rtheta_g_rad(r, theta):


    """
    Convert axis-angle representation to rotation matrix using Rodrigues' formula.
    
    Input:
        r: list or array (3,) - Rotation axis (should be unit vector)
        theta: float - Rotation angle in radians
    
    Output:
        g: list (3, 3) - Rotation matrix
    
    Usage Example:
        >>> import numpy as np
        >>> # 90-degree rotation around Z axis
        >>> axis = [0, 0, 1]
        >>> angle = np.pi/2
        >>> g = ol_rtheta_g_rad(axis, angle)
        >>> print("Rotation matrix:")
        >>> for row in g:
        ...     print(row)
        >>> # Should give approximately [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        
        >>> # 180-degree rotation around X axis
        >>> g_180x = ol_rtheta_g_rad([1, 0, 0], np.pi)
        >>> print("180° around X:")
        >>> for row in g_180x:
        ...     print(row)
    """
    g = [[0,0,0] for i in range(0,3)]

    g[0][0] = r[0] * r[0] * (1 - np.cos (theta)) + np.cos (theta);
    g[0][1] = r[0] * r[1] * (1 - np.cos (theta)) + r[2] * np.sin (theta);
    g[0][2] = r[0] * r[2] * (1 - np.cos (theta)) - r[1] * np.sin (theta);
    
    g[1][0] = r[1] * r[0] * (1 - np.cos (theta)) - r[2] * np.sin (theta);
    g[1][1] = r[1] * r[1] * (1 - np.cos (theta)) + np.cos (theta);
    g[1][2] = r[1] * r[2] * (1 - np.cos (theta)) + r[0] * np.sin (theta);
    
    g[2][0] = r[2] * r[0] * (1 - np.cos (theta)) + r[1] * np.sin (theta);
    g[2][1] = r[2] * r[1] * (1 - np.cos (theta)) - r[0] * np.sin (theta);
    g[2][2] = r[2] * r[2] * (1 - np.cos (theta)) + np.cos (theta);

    return g


def np_ol_rtheta_g_rad(r, theta):


    """
    Convert axis-angle representation to rotation matrix (NumPy version).
    Uses Rodrigues' rotation formula.
    
    Input:
        r: numpy array (3,) - Rotation axis (unit vector)
        theta: float - Rotation angle in radians
    
    Output:
        g: numpy array (3, 3) - Rotation matrix
    
    Usage Example:
        >>> import numpy as np
        >>> # 45-degree rotation around [1,1,1] axis
        >>> axis = np.array([1, 1, 1]) / np.sqrt(3)  # Normalize
        >>> angle = np.radians(45)
        >>> g = np_ol_rtheta_g_rad(axis, angle)
        >>> print("Rotation matrix:")
        >>> print(g)
        >>> 
        >>> # Verify it's a valid rotation matrix
        >>> det = np.linalg.det(g)
        >>> print("Determinant:", det)  # Should be 1
        >>> orthogonal = g @ g.T
        >>> print("G @ G.T (should be identity):")
        >>> print(orthogonal)
    """
    g=np.eye(3)

    g[0,0] = r[0] * r[0] * (1 - np.cos (theta)) + np.cos (theta);
    g[0,1] = r[0] * r[1] * (1 - np.cos (theta)) + r[2] * np.sin (theta);
    g[0,2] = r[0] * r[2] * (1 - np.cos (theta)) - r[1] * np.sin (theta);
    
    g[1,0] = r[1] * r[0] * (1 - np.cos (theta)) - r[2] * np.sin (theta);
    g[1,1] = r[1] * r[1] * (1 - np.cos (theta)) + np.cos (theta);
    g[1,2] = r[1] * r[2] * (1 - np.cos (theta)) + r[0] * np.sin (theta);
    
    g[2,0] = r[2] * r[0] * (1 - np.cos (theta)) + r[1] * np.sin (theta);
    g[2,1] = r[2] * r[1] * (1 - np.cos (theta)) - r[0] * np.sin (theta);
    g[2,2] = r[2] * r[2] * (1 - np.cos (theta)) + np.cos (theta);

    return g



def np_gmat2rodrigues(g):
    """
    Convert rotation matrix to Rodrigues-Frank vector representation.
    Rodrigues vector = rotation_axis * tan(angle/2)
    
    Input:
        g: numpy array (3, 3) - Rotation matrix
    
    Output:
        rodrigues: numpy array (3,) - Rodrigues-Frank vector
    
    Usage Example:
        >>> import numpy as np
        >>> # 90-degree rotation around Z
        >>> g = np.array([[0, -1, 0],
        ...               [1,  0, 0],
        ...               [0,  0, 1]], dtype=float)
        >>> rod = np_gmat2rodrigues(g)
        >>> print("Rodrigues vector:", rod)
        >>> # For 90° rotation: tan(45°) = 1, so rod ≈ [0, 0, 1]
        
        >>> # Small rotation (5 degrees around X)
        >>> angle = np.radians(5)
        >>> g_small = np.array([[1, 0, 0],
        ...                     [0, np.cos(angle), -np.sin(angle)],
        ...                     [0, np.sin(angle), np.cos(angle)]])
        >>> rod_small = np_gmat2rodrigues(g_small)
        >>> print("Small rotation Rodrigues:", rod_small)
    """
    axis, angle = np_ol_g_rtheta_rad(g)
    rodrigues = axis * np.tan(angle/2)
    return rodrigues



def np_rodrigues2gmat(rodrigues):
    """
    Convert Rodrigues-Frank vector to rotation matrix.
    
    Input:
        rodrigues: numpy array (3,) - Rodrigues-Frank vector (axis * tan(angle/2))
    
    Output:
        g: numpy array (3, 3) - Rotation matrix
    
    Usage Example:
        >>> import numpy as np
        >>> # Rodrigues vector for 90° around Z
        >>> rod = np.array([0, 0, 1])  # tan(45°) = 1
        >>> g = np_rodrigues2gmat(rod)
        >>> print("Rotation matrix:")
        >>> print(g)
        >>> # Should give approximately [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        
        >>> # Round trip test
        >>> g_original = np.array([[0.5, -0.866, 0],
        ...                        [0.866, 0.5, 0],
        ...                        [0, 0, 1]])  # 60° around Z
        >>> rod = np_gmat2rodrigues(g_original)
        >>> g_recovered = np_rodrigues2gmat(rod)
        >>> print("Original equals recovered:", np.allclose(g_original, g_recovered))
    """
    norm = np.linalg.norm(rodrigues)
    if norm < 1e-10:
        return np.eye(3)
    
    axis = rodrigues / norm
    angle = 2 * np.arctan(norm)
    g = np_ol_rtheta_g_rad(axis, angle)
    return g


# =====================================
# Vectorized Quaternion Operations
# =====================================


def np_g2quats(umatsa):
    """
    Convert multiple rotation matrices to quaternions (vectorized).
    Handles arrays of rotation matrices efficiently.
    
    Input:
        umatsa: numpy array (N, 3, 3) - Array of rotation matrices
    
    Output:
        Q: numpy array (4, N) - Array of quaternions [w, x, y, z] for each matrix
    
    Usage Example:
        >>> import numpy as np
        >>> # Create multiple rotation matrices
        >>> N = 100
        >>> from scipy.spatial.transform import Rotation as R
        >>> matrices = R.random(N).as_matrix()
        >>> 
        >>> # Convert all to quaternions
        >>> quats = np_g2quats(matrices)
        >>> print("Shape:", quats.shape)  # (4, 100)
        >>> print("First quaternion:", quats[:, 0])
        >>> 
        >>> # Verify quaternions are normalized
        >>> norms = np.linalg.norm(quats, axis=0)
        >>> print("All normalized:", np.allclose(norms, 1.0))
    """
    # Determine which formula to use based on matrix trace
    M22L0 = np.where(umatsa[:,2,2] < 0)[0]
    M22GE0 = np.where(umatsa[:,2,2] >= 0)[0]
    
    if M22L0.shape[0] > 0:
        M00GM11 = np.where(umatsa[M22L0,0,0] > umatsa[M22L0,1,1])[0]  
        M00LEM11 = np.where(umatsa[M22L0,0,0] <= umatsa[M22L0,1,1])[0]
    else:
        M00GM11 = np.array([])
        M00LEM11 = np.array([])
        
    if M22GE0.shape[0] > 0:
        M00LnM11 = np.where(umatsa[M22GE0,0,0] < -1*umatsa[M22GE0,1,1])[0]  
        M00GEnM11 = np.where(umatsa[M22GE0,0,0] >= -1*umatsa[M22GE0,1,1])[0]
    else:
        M00LnM11 = np.array([])
        M00GEnM11 = np.array([])
    
    Q = np.empty((4, umatsa.shape[0]))
    T = np.empty((umatsa.shape[0]))
    
    # Four cases based on which diagonal element is largest
    try:
        T[M22L0[M00GM11]] = 1 + umatsa[M22L0[M00GM11],0,0] - umatsa[M22L0[M00GM11],1,1] - umatsa[M22L0[M00GM11],2,2]
        Q[:,M22L0[M00GM11]] = [umatsa[M22L0[M00GM11],1, 2]-umatsa[M22L0[M00GM11],2, 1],  T[M22L0[M00GM11]],  umatsa[M22L0[M00GM11],0, 1]+umatsa[M22L0[M00GM11],1, 0],  
                               umatsa[M22L0[M00GM11],2, 0]+umatsa[M22L0[M00GM11],0, 2]]
    except:
        pass
    
    try:
        T[M22L0[M00LEM11]]  = 1 - umatsa[M22L0[M00LEM11],0, 0] + umatsa[M22L0[M00LEM11],1, 1] - umatsa[M22L0[M00LEM11],2, 2]
        Q[:,M22L0[M00LEM11]] = [umatsa[M22L0[M00LEM11],2, 0]-umatsa[M22L0[M00LEM11],0, 2],  umatsa[M22L0[M00LEM11],0, 1]+umatsa[M22L0[M00LEM11],1, 0],
                                T[M22L0[M00LEM11]],  umatsa[M22L0[M00LEM11],1, 2]+umatsa[M22L0[M00LEM11],2, 1]]    
    except:
        pass    
    
    try:
        T[M22GE0[M00LnM11]] = 1 - umatsa[M22GE0[M00LnM11],0, 0] - umatsa[M22GE0[M00LnM11],1, 1] + umatsa[M22GE0[M00LnM11],2, 2]
        Q[:,M22GE0[M00LnM11]] = [umatsa[M22GE0[M00LnM11],0, 1]-umatsa[M22GE0[M00LnM11],1, 0],  umatsa[M22GE0[M00LnM11],2, 0]+umatsa[M22GE0[M00LnM11],0, 2],  
             umatsa[M22GE0[M00LnM11],1, 2]+umatsa[M22GE0[M00LnM11],2, 1],T[M22GE0[M00LnM11]]]
    except:
        pass
    
    try:
        T[M22GE0[M00GEnM11]] = 1 + umatsa[M22GE0[M00GEnM11],0, 0] + umatsa[M22GE0[M00GEnM11],1, 1] + umatsa[M22GE0[M00GEnM11],2, 2]
        Q[:,M22GE0[M00GEnM11]] =[T[M22GE0[M00GEnM11]],  umatsa[M22GE0[M00GEnM11],1, 2]-umatsa[M22GE0[M00GEnM11],2, 1],  umatsa[M22GE0[M00GEnM11],2, 0]-umatsa[M22GE0[M00GEnM11],0, 2],  
                                 umatsa[M22GE0[M00GEnM11],0, 1]-umatsa[M22GE0[M00GEnM11],1, 0]]
    except:
        pass

    # Normalize
    Q[0,:] *= 0.5 / np.sqrt(T)
    Q[1,:] *= 0.5 / np.sqrt(T)
    Q[2,:] *= 0.5 / np.sqrt(T)
    Q[3,:] *= 0.5 / np.sqrt(T)
    return Q



def Qlog(QM):
    """
    Compute the logarithm of quaternions (quaternion logarithm map).
    Maps quaternions from unit sphere S³ to tangent space at identity.
    
    Input:
        QM: numpy array (4, N, M) - Array of quaternions [w, x, y, z, ...]
    
    Output:
        qlog: numpy array (4, N, M) - Logarithm of quaternions
                                       qlog[0] = 0, qlog[1:4] = v * angle
    
    Usage Example:
        >>> import numpy as np
        >>> # Single quaternion for 90° rotation
        >>> q = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> # Reshape for function (4, 1, 1)
        >>> QM = q.reshape(4, 1, 1)
        >>> qlog = Qlog(QM)
        >>> print("Quaternion log:")
        >>> print(qlog[:, 0, 0])
        >>> # qlog[0] should be 0, qlog[3] should be π/2
    """
    qlog = QM.copy()
    qlog[0,:,:] = qlog[0,:,:] * 0
    qlog[1:4,range(qlog.shape[1]),range(qlog.shape[2])] = 1
    
    # Normalize the vector part
    idxs = np.where(np.linalg.norm(qlog[1:4,:,:], axis=0) != 0)
    qlog[1:4,idxs[0],idxs[1]] = (qlog[1:4,idxs[0],idxs[1]] / np.linalg.norm(qlog[1:4,idxs[0],idxs[1]], axis=0))
    
    idxs = np.where(np.linalg.norm(qlog[1:4,:,:], axis=0) == 0)
    qlog[1:4,idxs[0],idxs[1]] = 0
    
    # Scale by angle
    qlog[1:4,:,:] = qlog[1:4,:,:] * np.tile(np.arccos(QM[0,:,:]/np.linalg.norm(QM[:,:,:], axis=0)),(3,1,1))
    
    return qlog



def Qproduct(P, Q):
    """
    Compute quaternion product for arrays of quaternions.
    Computes P * Q for each pair of quaternions.
    
    Input:
        P: numpy array (4, N) - First array of quaternions [w, x, y, z]
        Q: numpy array (4, N) - Second array of quaternions [w, x, y, z]
    
    Output:
        result: numpy array (4, N) - Product quaternions P * Q
    
    Usage Example:
        >>> import numpy as np
        >>> # Multiple quaternion pairs
        >>> N = 100
        >>> P = np.random.randn(4, N)
        >>> P = P / np.linalg.norm(P, axis=0)  # Normalize
        >>> Q = np.random.randn(4, N)
        >>> Q = Q / np.linalg.norm(Q, axis=0)
        >>> 
        >>> result = Qproduct(P, Q)
        >>> print("Result shape:", result.shape)  # (4, 100)
        >>> # Verify results are normalized
        >>> norms = np.linalg.norm(result, axis=0)
        >>> print("All normalized:", np.allclose(norms, 1.0))
    """
    Ones = np.ones(Q.shape)
    Ones[1:4,:] = -1 * Ones[1:4,:]
    Q0 = P.T.dot(Q * Ones)
    
    idxs = [[0,1],[1,0],[2,3],[3,2]]
    signs = [1,1,1,-1]
    Q1 = np.zeros(Q0.shape)
    for idx, sgn in zip(idxs, signs):
        p1 = np.zeros(P.shape)
        p1[0,:] = P[idx[0],:]    
        q1 = np.zeros(Q.shape)
        q1[0,:] = sgn * Q[idx[1],:]
        Q1 += p1.T.dot(q1)

    idxs = [[0,2],[2,0],[1,3],[3,1]]
    signs = [1,1,-1,1]
    Q2 = np.zeros(Q0.shape)
    for idx, sgn in zip(idxs, signs):
        p1 = np.zeros(P.shape)
        p1[0,:] = P[idx[0],:]    
        q1 = np.zeros(Q.shape)
        q1[0,:] = sgn * Q[idx[1],:]
        Q2 += p1.T.dot(q1)
    
    idxs = [[0,3],[3,0],[1,2],[2,1]]
    signs = [1,1,1,-1]
    Q3 = np.zeros(Q0.shape)
    for idx, sgn in zip(idxs, signs):
        p1 = np.zeros(P.shape)
        p1[0,:] = P[idx[0],:]    
        q1 = np.zeros(Q.shape)
        q1[0,:] = sgn * Q[idx[1],:]
        Q3 += p1.T.dot(q1)
    
    return np.stack((Q0, Q1, Q2, Q3))



def QMatproduct(sym, Q):
    """
    Multiply a single symmetry quaternion with multiple quaternions.
    Applies the same symmetry operation to an array of orientations.
    
    Input:
        sym: numpy array (4,) - Single symmetry quaternion [w, x, y, z]
        Q: numpy array (4, N) - Array of quaternions to transform
    
    Output:
        SQ: numpy array (4, N) - Transformed quaternions (sym * Q)
    
    Usage Example:
        >>> import numpy as np
        >>> # Define a 90° rotation symmetry around Z
        >>> sym = np.array([np.cos(np.pi/4), 0, 0, np.sin(np.pi/4)])
        >>> 
        >>> # Multiple orientations
        >>> N = 50
        >>> Q = np.random.randn(4, N)
        >>> Q = Q / np.linalg.norm(Q, axis=0)
        >>> 
        >>> # Apply symmetry to all
        >>> SQ = QMatproduct(sym, Q)
        >>> print("Transformed shape:", SQ.shape)  # (4, 50)
    """
    SQ = Q.copy()
    SQ[0,:] = sym[0]*Q[0,:] - sym[1]*Q[1,:] - sym[2]*Q[2,:] - sym[3]*Q[3,:]
    SQ[1,:] = sym[0]*Q[1,:] + sym[1]*Q[0,:] + sym[2]*Q[3,:] - sym[3]*Q[2,:]
    SQ[2,:] = sym[0]*Q[2,:] + sym[2]*Q[0,:] - sym[1]*Q[3,:] + sym[3]*Q[1,:]
    SQ[3,:] = sym[0]*Q[3,:] + sym[3]*Q[0,:] + sym[1]*Q[2,:] - sym[2]*Q[1,:]
    
    return SQ


# =====================================
# Orientation Sampling and Gridding
# =====================================


def grid_s1(resol, grids=6):
    """
    Generate uniformly distributed points on S¹ (circle).
    Used for sampling the third Euler angle in Hopf coordinates.
    
    Input:
        resol: int - Resolution parameter (number of points = 2^resol * grids)
        grids: int - Grid multiplier (default: 6)
    
    Output:
        points: list of floats - Angles in radians from 0 to 2π
    
    Usage Example:
        >>> points = grid_s1(resol=2, grids=6)
        >>> print("Number of points:", len(points))
        >>> # 2^2 * 6 = 24 points
        >>> print("First few points:", points[:5])
        >>> print("Last point:", points[-1])
        
        >>> # Higher resolution
        >>> points_fine = grid_s1(resol=4, grids=6)
        >>> print("Fine grid points:", len(points_fine))
        >>> # 2^4 * 6 = 96 points
    """
    number_points = (2 ** resol) * grids

    interval = 2 * np.pi / number_points

    points = [interval / 2 + i * interval for i in range(number_points)]

    return points



def hopf2quat(Points):
    """
    Convert Hopf coordinates to quaternions.
    Hopf coordinates (θ, φ, ψ) parameterize the unit quaternion sphere S³.
    
    Input:
        Points: list of tuples - Each tuple contains (theta, phi, psi) in radians
                                 theta ∈ [0, π], phi ∈ [0, 2π], psi ∈ [0, 2π]
    
    Output:
        quats: list of lists - Quaternions [w, x, y, z] for each point
    
    Usage Example:
        >>> import numpy as np
        >>> # Single point in Hopf coordinates
        >>> points = [(np.pi/2, 0, 0)]
        >>> quats = hopf2quat(points)
        >>> print("Quaternion:", quats[0])
        
        >>> # Multiple points
        >>> points_multiple = [
        ...     (0, 0, 0),           # Identity
        ...     (np.pi/2, 0, 0),     # 
        ...     (np.pi/2, np.pi, 0)  # 
        ... ]
        >>> quats_multiple = hopf2quat(points_multiple)
        >>> print("Number of quaternions:", len(quats_multiple))
    """
    quats = []

    for i in range(len(Points)):
        x4 = math.sin(Points[i][0] / 2) * math.sin(Points[i][1] + Points[i][2] / 2)

        x1 = math.cos(Points[i][0] / 2) * math.cos(Points[i][2] / 2)

        x2 = math.cos(Points[i][0] / 2) * math.sin(Points[i][2] / 2)

        x3 = math.sin(Points[i][0] / 2) * math.cos(Points[i][1] + Points[i][2] / 2)

        quats.append([x1, x2, x3, x4])

    return quats



def nside2npix(nside):
    """
    Calculate the number of pixels in a HEALPix map.
    HEALPix = Hierarchical Equal Area isoLatitude Pixelization.
    
    Input:
        nside: int - HEALPix resolution parameter (must be power of 2)
    
    Output:
        npix: int - Total number of pixels = 12 * nside²
    
    Usage Example:
        >>> # Low resolution
        >>> npix_low = nside2npix(nside=1)
        >>> print("Pixels for nside=1:", npix_low)  # 12
        
        >>> # Medium resolution
        >>> npix_med = nside2npix(nside=16)
        >>> print("Pixels for nside=16:", npix_med)  # 3072
        
        >>> # High resolution
        >>> npix_high = nside2npix(nside=128)
        >>> print("Pixels for nside=128:", npix_high)  # 196608
    """
    return 12 * nside * nside



def pix2ang_nest(nside, ipix, pix2x, pix2y):
    """
    Convert HEALPix pixel index to spherical coordinates (theta, phi).
    Uses NESTED indexing scheme.
    
    Input:
        nside: int - HEALPix resolution parameter
        ipix: int - Pixel index (0 to 12*nside² - 1)
        pix2x: list - Lookup table for x-coordinate (from mk_pix2xy)
        pix2y: list - Lookup table for y-coordinate (from mk_pix2xy)
    
    Output:
        theta: float - Colatitude angle in radians [0, π]
        phi: float - Azimuthal angle in radians [0, 2π]
    
    Usage Example:
        >>> # Setup lookup tables
        >>> pix2x, pix2y = mk_pix2xy()
        >>> 
        >>> # Convert pixel 0 at nside=4
        >>> nside = 4
        >>> theta, phi = pix2ang_nest(nside, ipix=0, pix2x=pix2x, pix2y=pix2y)
        >>> print(f"Pixel 0: theta={np.degrees(theta):.2f}°, phi={np.degrees(phi):.2f}°")
        >>> 
        >>> # Convert all pixels at low resolution
        >>> nside = 2
        >>> npix = nside2npix(nside)
        >>> for i in range(npix):
        ...     theta, phi = pix2ang_nest(nside, i, pix2x, pix2y)
        ...     # Process coordinates...
    """
    jrll = np.array([2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4])
    jpll = np.array([1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7])

    if (nside < 1) or (nside > 8192):
        raise Exception('nside out of range:', nside)

    if (ipix < 0) or (ipix > 12 * nside * nside - 1):
        raise Exception('ipix out of range:', ipix)

    fn = 1. * nside
    fact1 = 1. / (3. * fn * fn)
    fact2 = 2. / (3. * fn)
    nl4 = 4 * nside

    npface = nside * nside
    face_num = int(ipix / npface)
    ipf = int(ipix % npface)

    ip_low = int(ipf % 1024)
    ip_trunc = ipf / 1024
    ip_med = int(ip_trunc % 1024)
    ip_hi = int(ip_trunc / 1024)

    ix = 1024 * pix2x[ip_hi] + 32 * pix2x[ip_med] + pix2x[ip_low]
    iy = 1024 * pix2y[ip_hi] + 32 * pix2y[ip_med] + pix2y[ip_low]

    jrt = ix + iy
    jpt = ix - iy

    jr = jrll[face_num] * nside - jrt - 1
    nr = nside
    z = (2 * nside - jr) * fact2
    kshift = int((jr - nside) % 2)

    if jr < nside:
        nr = jr
        z = 1. - nr * nr * fact1
        kshift = 0

    elif jr > 3 * nside:
        nr = nl4 - jr
        z = - 1. + nr * nr * fact1
        kshift = 0

    theta = np.arccos(z)

    jp = (jpll[face_num] * nr + jpt + 1 + kshift) / 2

    if jp > nl4:
        jp = jp - nl4

    if jp < 1:
        jp = jp + nl4

    phi = (jp - (kshift + 1) * 0.5) * (np.pi / 2 / nr)

    return theta, phi



def mk_pix2xy():
    """
    Create lookup tables for HEALPix pixel index to (x,y) coordinate conversion.
    These tables are used by pix2ang_nest function.
    
    Input:
        None
    
    Output:
        pix2x: list (1024,) - X-coordinate lookup table
        pix2y: list (1024,) - Y-coordinate lookup table
    
    Usage Example:
        >>> pix2x, pix2y = mk_pix2xy()
        >>> print("Table size:", len(pix2x))  # 1024
        >>> print("First few X values:", pix2x[:5])
        >>> print("First few Y values:", pix2y[:5])
        >>> 
        >>> # Use with pix2ang_nest
        >>> nside = 8
        >>> theta, phi = pix2ang_nest(nside, ipix=100, pix2x=pix2x, pix2y=pix2y)
    """
    pix2x = []
    pix2y = []

    for kpix in range(1024):
        jpix = kpix
        IX = 0
        IY = 0
        IP = 1

        while jpix != 0:
            ID = int(jpix % 2)
            jpix /= 2
            IX = ID * IP + IX

            ID = int(jpix % 2)
            jpix /= 2
            IY = ID * IP + IY

            IP = 2 * IP

        pix2x.append(IX)
        pix2y.append(IY)

    return pix2x, pix2y



def simple_grid(resol):
    """
    Generate a uniform grid of quaternions covering orientation space (SO(3)).
    Uses HEALPix for S² sampling and uniform sampling for S¹, combined via Hopf fibration.
    
    Input:
        resol: int - Resolution parameter
                     - Number of S¹ points: 2^resol * 6
                     - HEALPix nside: 2^resol
                     - Total quaternions: 12 * 4^resol * (2^resol * 6)
    
    Output:
        quats: list of lists - Uniformly distributed quaternions [w, x, y, z]
    
    Usage Example:
        >>> # Low resolution grid (fast)
        >>> quats_low = simple_grid(resol=2)
        >>> print("Low resolution quaternions:", len(quats_low))
        >>> # 12 * 4^2 * (2^2 * 6) = 12 * 16 * 24 = 4608
        
        >>> # Medium resolution (moderate)
        >>> quats_med = simple_grid(resol=3)
        >>> print("Medium resolution quaternions:", len(quats_med))
        >>> # 12 * 4^3 * (2^3 * 6) = 12 * 64 * 48 = 36864
        
        >>> # Convert to rotation matrices for use
        >>> import numpy as np
        >>> q_sample = np.array(quats_low[0])
        >>> R = quat_to_mat(q_sample)
        >>> print("Sample rotation matrix:")
        >>> print(R)
        
        >>> # Use for texture analysis, pole figures, etc.
        >>> # Each quaternion represents a unique crystal orientation
    """
    # Generate points on S¹ (circle)
    Psi_Points = grid_s1(resol)

    # Generate HEALPix points on S² (sphere)
    Nside = 2 ** resol
    numpixels = nside2npix(Nside)
    pix2x, pix2y = mk_pix2xy()

    Healpix_Points = []
    for i in range(numpixels):
        theta, phi = pix2ang_nest(Nside, i, pix2x, pix2y)
        Healpix_Points.append([theta, phi])

    # Combine S² and S¹ to get S³ points (Hopf fibration)
    S3_Points = [[theta, phi, psi] for [theta, phi] in Healpix_Points for psi in Psi_Points]

    # Convert Hopf coordinates to quaternions
    quats = hopf2quat(S3_Points)

    return quats


def rotation_from_axis_angle(ax,an,deg=False):

    """
    Generate rotation matrix from axis-angle representation using Rodrigues formula.
    
    Creates 3×3 rotation matrix for rotation by 'angle' radians around 'axis'.
    Implements Rodrigues' rotation formula for arbitrary axis rotations.
    
    Input:
        axis (array [3]): Rotation axis (will be normalized)
        angle (float): Rotation angle in radians
    
    Output:
        numpy.ndarray (3×3): Rotation matrix R
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # 90° rotation around z-axis
        >>> axis = np.array([0, 0, 1])
        >>> angle = np.pi / 2
        >>> R = rotation_from_axis_angle(axis, angle)
        >>> print(R)
        >>> # [[0, -1, 0],
        >>> #  [1,  0, 0],
        >>> #  [0,  0, 1]]
        >>> 
        >>> # Verify: rotate [1,0,0] to [0,1,0]
        >>> v = np.array([1, 0, 0])
        >>> v_rot = R.dot(v)
        >>> print(v_rot)  # [0, 1, 0]
        >>> 
        >>> # 120° rotation around [111]
        >>> axis_111 = np.array([1, 1, 1])
        >>> R_111 = rotation_from_axis_angle(axis_111, 2*np.pi/3)
    
    Notes:
        - Axis is automatically normalized
        - Right-hand rule: thumb along axis, fingers show rotation
        - Preserves vector lengths (orthogonal matrix)
        - det(R) = 1 (proper rotation)
        - Used in crystallographic symmetry operations
    
    Formula (Rodrigues):
        R = I + sin(θ)K + (1-cos(θ))K²
        where K is the skew-symmetric matrix of the axis
    """
    if deg:
        an=an*np.pi/180
    ax/=norm(ax)
    return np.eye(3)*np.cos(an)+np.sin(an)*np.cross(ax,np.eye(3)) + (1-np.cos(an))*np.outer(ax,ax)

def ol_g_R(g):

    """
    Convert rotation matrix to Rodrigues-Frank vector (list version).
    
    Calculates Rodrigues-Frank vector R = r·tan(θ/2) from rotation matrix,
    where r is the rotation axis and θ is the rotation angle.
    
    Input:
        g (list 3×3): Rotation matrix
    
    Output:
        list [3]: Rodrigues-Frank vector
    
    Usage Example:
        >>> # 90° rotation around Z
        >>> g = [[0, -1, 0],
        ...      [1,  0, 0],
        ...      [0,  0, 1]]
        >>> R = ol_g_R(g)
        >>> print(R)  # [0, 0, 1] (since tan(45°)=1)
    
    Notes:
        - Compact representation: 3 parameters vs 9 for matrix
        - R = axis · tan(angle/2)
        - Magnitude ||R|| = tan(θ/2)
        - Direction = rotation axis
        - Singular at θ = 180° (infinite magnitude)
    
    Formula:
        R = r · tan(θ/2)
        where (r, θ) from ol_g_rtheta_rad(g)
    """
    #Quey
    r,theta = ol_g_rtheta_rad (g)
    R=[0.,0.,0.]
    for i in range(0,3):
        R[i]=r[i]*np.tan(theta/2)
    return R


def np_ol_g_R(g):

    """
    Convert rotation matrix to Rodrigues-Frank vector (NumPy version).
    
    NumPy implementation returning Rodrigues-Frank vector as numpy array.
    
    Input:
        g (numpy.ndarray 3×3): Rotation matrix
    
    Output:
        numpy.ndarray [3]: Rodrigues-Frank vector
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Identity rotation
        >>> g_id = np.eye(3)
        >>> R_id = np_ol_g_R(g_id)
        >>> print(R_id)  # [0, 0, 0]
        >>> 
        >>> # General rotation
        >>> g = rotation_from_axis_angle([1,1,1], np.pi/3)
        >>> R = np_ol_g_R(g)
        >>> print(f"Rodrigues vector: {R}")
    
    Notes:
        - Efficient numpy implementation
        - Useful in grain boundary analysis
        - Common in crystal plasticity codes
    """
    #Quey
    r,theta = np_ol_g_rtheta_rad (g)
    R=r*np.tan(theta/2)

    return R


def ol_R_g (R):


    """
    Convert Rodrigues-Frank vector to rotation matrix (list version).
    
    Reconstructs rotation matrix from Rodrigues-Frank representation.
    
    Input:
        R (list [3]): Rodrigues-Frank vector
    
    Output:
        list (3×3): Rotation matrix
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Rodrigues vector for 90° around Z
        >>> R = [0, 0, 1]  # tan(45°) = 1
        >>> g = ol_R_g(R)
        >>> # Should give 90° rotation around Z
    
    Notes:
        - Inverse of ol_g_R()
        - Extracts angle: θ = 2·arctan(||R||)
        - Extracts axis: r = R/||R||
    
    Formula:
        θ = 2·arctan(||R||)
        r = R / ||R||
        g = ol_rtheta_g_rad(r, θ)
    """
    norm = np.sqrt(sum([ri*ri for ri in R]))
    r = [ri/norm for ri in R]
    theta = 2*np.arctan(norm)


    g=ol_rtheta_g_rad(r, theta)


    return g


def np_ol_R_g (R):


    """
    Convert Rodrigues-Frank vector to rotation matrix (NumPy version).
    
    NumPy implementation of Rodrigues-Frank to rotation matrix conversion.
    
    Input:
        R (numpy.ndarray [3]): Rodrigues-Frank vector
    
    Output:
        numpy.ndarray (3×3): Rotation matrix
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Small rotation
        >>> R_small = np.array([0.01, 0.02, 0.03])
        >>> g = np_ol_R_g(R_small)
        >>> print(f"Rotation angle: {np.degrees(2*np.arctan(np.linalg.norm(R_small))):.2f}°")
        >>> 
        >>> # Round trip test
        >>> g_orig = rotation_from_axis_angle([1,0,0], 0.5)
        >>> R_test = np_ol_g_R(g_orig)
        >>> g_recovered = np_ol_R_g(R_test)
        >>> print(np.allclose(g_orig, g_recovered))  # True
    
    Notes:
        - Inverse of np_ol_g_R()
        - Handles small rotations accurately
        - Returns identity for zero vector
    """
    norm = np.sqrt(R.dot(R))
    r=R/norm
    theta = 2*np.arctan(norm)


    g=np_ol_rtheta_g_rad(r, theta)


    return g


        

def ol_g_R2(g):

    """
    Alternative Rodrigues-Frank conversion (list version, method 2).
    
    Second implementation of rotation matrix to Rodrigues-Frank vector.
    May use different numerical approach for stability.
    
    Input:
        g (list 3×3): Rotation matrix
    
    Output:
        list [3]: Rodrigues-Frank vector
    
    Usage Example:
        >>> g = [[0, -1, 0],
        ...      [1,  0, 0],
        ...      [0,  0, 1]]
        >>> R = ol_g_R2(g)
    
    Notes:
        - Alternative implementation for numerical comparison
        - Should give same result as ol_g_R()
        - Useful for validation
    """
    #Poulsen
    gmm = g[0][0]+g[1][1]+g[2][2]
    R=[0.,0.,0.]
    epsilon = permut_tensor3()
    delta = kronecker()
    for i in range(0,3):
        for j in range(0,3):
            for k in range(0,3):
                R[i]=R[i]+(epsilon[i][j][k]*g[j][k])/(1+gmm)
                
    return R


def np_ol_g_R2(g,epsilon, delta):

    """
    Alternative Rodrigues-Frank conversion (NumPy version, method 2).
    
    NumPy implementation of alternative Rodrigues-Frank calculation.
    
    Input:
        g (numpy.ndarray 3×3): Rotation matrix
    
    Output:
        numpy.ndarray [3]: Rodrigues-Frank vector
    
    Usage Example:
        >>> import numpy as np
        >>> g = np.eye(3)
        >>> R = np_ol_g_R2(g)
    
    Notes:
        - Alternative implementation
        - For numerical stability comparison
    """
    #Poulsen
    gmm = np.trace(g)
    R = np.einsum('ijk,jk',epsilon,g)/(1+gmm)    
    return R


def ol_R_g2(R):

    """
    Alternative Rodrigues-Frank to rotation matrix (list version, method 2).
    
    Second implementation of Rodrigues-Frank to rotation matrix conversion.
    
    Input:
        R (list [3]): Rodrigues-Frank vector
    
    Output:
        list (3×3): Rotation matrix
    
    Usage Example:
        >>> R = [0, 0, 0.5]
        >>> g = ol_R_g2(R)
    
    Notes:
        - Alternative implementation
        - Should match ol_R_g() results
    """
    #Poulsen
    r2=sum([ri*ri for ri in R])
    
    epsilon = permut_tensor3()
    delta = kronecker()
    g=[]
    for i in range(0,3):
        gj=[]
        for j in range(0,3):
            er=0.
            for k in range(0,3):
                er=er+2*epsilon[i][j][k]*R[k]
            gj.append(1./(1+r2)*((1-r2)*delta[i][j]+2*R[i]*R[j]+er))
        g.append(gj)
            
                    
    return g


def np_ol_R_g2(R,epsilon, delta):

    """
    Alternative Rodrigues-Frank to rotation matrix (NumPy version, method 2).
    
    NumPy implementation of alternative conversion method.
    
    Input:
        R (numpy.ndarray [3]): Rodrigues-Frank vector
    
    Output:
        numpy.ndarray (3×3): Rotation matrix
    
    Usage Example:
        >>> import numpy as np
        >>> R = np.array([0.1, 0.2, 0.3])
        >>> g = np_ol_R_g2(R)
    
    Notes:
        - Alternative implementation for comparison
    """
    #Poulsen
    r2=R.dot(R)
    g= 1./(1+r2)*((1-r2)*delta+2*np.einsum('i,j',R,R)+np.einsum('ijk,k',2*epsilon,R))
    return g


def np_ol_R_q2(R):

    """
    Convert Rodrigues-Frank vector to quaternion (method 2).
    
    Transforms Rodrigues-Frank representation to unit quaternion [w,x,y,z].
    
    Input:
        R (numpy.ndarray [3]): Rodrigues-Frank vector
    
    Output:
        numpy.ndarray [4]: Unit quaternion [w, x, y, z]
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Small rotation
        >>> R = np.array([0.01, 0, 0])
        >>> q = np_ol_R_q2(R)
        >>> print(q)  # Near [1, 0.01, 0, 0]
        >>> print(f"Norm: {np.linalg.norm(q)}")  # 1.0
    
    Notes:
        - Output is normalized quaternion
        - w ≥ 0 convention
        - Quaternion avoids singularities of Rodrigues-Frank
    
    Formula:
        θ = 2·arctan(||R||)
        q = [cos(θ/2), r·sin(θ/2)]
        where r = R/||R||
    """
    #Poulsen
    q=np.empty(4)
    r2 = R.dot(R)
    q[0]=1./np.sqrt(1+r2);
    q[1:]=R/np.sqrt(1+r2)

    return q


def np_ol_g_q2(g):

    """
    Convert rotation matrix to quaternion (method 2).
    
    Extracts unit quaternion representation from rotation matrix.
    
    Input:
        g (numpy.ndarray 3×3): Rotation matrix
    
    Output:
        numpy.ndarray [4]: Unit quaternion [w, x, y, z]
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # 90° around Z
        >>> g = np.array([[0, -1, 0],
        ...               [1,  0, 0],
        ...               [0,  0, 1]], dtype=float)
        >>> q = np_ol_g_q2(g)
        >>> print(q)  # [0.707, 0, 0, 0.707]
    
    Notes:
        - Uses Rodrigues-Frank as intermediate
        - Normalized output
        - Stable for all rotation angles
    
    Formula:
        g → R → q
        (via np_ol_g_R and np_ol_R_q2)
    """
    #Poulsen
    eps = 1e-6;
    q=np.empty(4)
    q[0] = 0.5*np.sqrt(np.trace(g)+1)
    
    if abs(q[0]) > eps:
        q[1]=1./4./q[0]*(g[2,1]-g[1,2])
        q[2]=1./4./q[0]*(g[2,0]-g[0,2])
        q[3]=1./4./q[0]*(g[0,1]-g[1,0])
    else:
        for i in range(0,2):
            q[i+1]=np.sqrt((g[i,i]+1)/2)
        
        m = 1+np.where(q[1:]==max(q[1:]))[0][0]
        for i in range(0,3):
            q[i]*=np.sign(g[i - 1][m - 1])
            
            
    
    return q
#def np_ol_q_U2(q):
#    #Poulsen
#    g=np.empty((3,3))
#    
#    for i in range(0,2):
#        g[i,i]=2*(q[0]**2+q[i+1]**2)-1
#    
#    g[1,0] = 2*(q[1]*q[2]+q[0]*q[3])
#    g[0,1] = 2*(q[1]*q[2]-q[0]*q[3])
#
#    g[2,0] = 2*(q[1]*q[3]-q[0]*q[2])
#    g[0,2] = 2*(q[1]*q[3]+q[0]*q[2])
#    
#    g[2,1] = 2*(q[2]*q[3]+q[0]*q[1])
#    g[1,2] = 2*(q[2]*q[3]-q[0]*q[1])
#
#    return g
#

def np_ol_q_g(q):

    """
    Convert quaternion to rotation matrix.
    
    Transforms unit quaternion [w,x,y,z] to 3×3 rotation matrix.
    
    Input:
        q (numpy.ndarray [4]): Unit quaternion [w, x, y, z]
    
    Output:
        numpy.ndarray (3×3): Rotation matrix
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Identity quaternion
        >>> q_id = np.array([1, 0, 0, 0])
        >>> g_id = np_ol_q_g(q_id)
        >>> print(np.allclose(g_id, np.eye(3)))  # True
        >>> 
        >>> # 180° around X
        >>> q_180x = np.array([0, 1, 0, 0])
        >>> g_180x = np_ol_q_g(q_180x)
        >>> print(g_180x)
        >>> # [[ 1,  0,  0],
        >>> #  [ 0, -1,  0],
        >>> #  [ 0,  0, -1]]
    
    Notes:
        - Input should be normalized
        - Avoids gimbal lock
        - Efficient for interpolation
    
    Formula:
        g₁₁ = 1 - 2(y² + z²)
        g₁₂ = 2(xy - zw)
        g₁₃ = 2(xz + yw)
        ... (full 3×3 matrix)
    """
    #Poulsen
    g=np.empty((3,3))
    
    g[0,0]=q[0]**2+q[1]**2-q[2]**2-q[3]**2
    g[1,1]=q[0]**2-q[1]**2+q[2]**2-q[3]**2
    g[2,2]=q[0]**2-q[1]**2-q[2]**2+q[3]**2
    
    for i in range(0,2):
        g[i,i]=2*(q[0]**2+q[i+1]**2)-1
    
    g[1,0] = 2*(q[1]*q[2]-q[0]*q[3])
    g[0,1] = 2*(q[1]*q[2]+q[0]*q[3])

    g[2,0] = 2*(q[1]*q[3]+q[0]*q[2])
    g[0,2] = 2*(q[1]*q[3]-q[0]*q[2])
    
    g[2,1] = 2*(q[2]*q[3]-q[0]*q[1])
    g[1,2] = 2*(q[2]*q[3]+q[0]*q[1])

    return g



        

def active_rotation(an, aboutaxis, deg=False):

    """
    Perform active rotation of vector v by rotation matrix g.
    
    Rotates the vector itself while keeping coordinate system fixed.
    v' = g · v (matrix-vector multiplication).
    
    Input:
        g (array 3×3): Rotation matrix
        v (array [3]): Vector to rotate
    
    Output:
        numpy.ndarray [3]: Rotated vector v'
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # 90° rotation around Z (active)
        >>> g_z90 = np.array([[0, -1, 0],
        ...                   [1,  0, 0],
        ...                   [0,  0, 1]])
        >>> v = np.array([1, 0, 0])
        >>> v_rot = active_rotation(g_z90, v)
        >>> print(v_rot)  # [0, 1, 0]
        >>> 
        >>> # Verify: [1,0,0] rotates to [0,1,0]
        >>> assert np.allclose(v_rot, [0, 1, 0])
    
    Notes:
        - Active = rotate the object
        - Compare with passive_rotation (rotate coordinates)
        - Common in physics and mechanics
        - v' = g · v
    
    Formula:
        v'_i = g_{ij} v_j
    """
    if deg:
        an=an*np.pi/180.        
    if aboutaxis.lower()=='z':
        R = np.array([[np.cos(an),-np.sin(an),0],[np.sin(an),np.cos(an),0],[0,0,1]]);
    elif aboutaxis.lower()=='x':
        R = np.array([[1,0,0],[0,np.cos(an),-np.sin(an)],[0,np.sin(an),np.cos(an)]]);
    elif aboutaxis.lower()=='y':
        R = np.array([[np.cos(an),0,np.sin(an)],[0,1,0],[-np.sin(an),0,np.cos(an)]]);
    
    return R
    


def passive_rotation(an, aboutaxis, deg=False):

    """
    Perform passive rotation (coordinate transformation).
    
    Rotates the coordinate system while vector stays fixed in space.
    v' = g^T · v (transpose of rotation matrix).
    
    Input:
        g (array 3×3): Rotation matrix
        v (array [3]): Vector in old coordinates
    
    Output:
        numpy.ndarray [3]: Vector components in new coordinates
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Rotate coordinates 90° around Z
        >>> g_z90 = np.array([[0, -1, 0],
        ...                   [1,  0, 0],
        ...                   [0,  0, 1]])
        >>> v = np.array([0, 1, 0])
        >>> v_new_coords = passive_rotation(g_z90, v)
        >>> print(v_new_coords)  # [1, 0, 0]
    
    Notes:
        - Passive = rotate coordinate system
        - v' = g^T · v (transpose)
        - Common in crystallography (sample↔crystal)
        - Inverse of active rotation for same g
    
    Formula:
        v'_i = g_{ji} v_j = g^T_{ij} v_j
    """
    R=np.transpose(active_rotation(an, aboutaxis, deg=deg));    
    
    return R
    

def euler_angles_reduction(Phi1,PHI,Phi2):


    """
    Reduce Euler angles to fundamental zone using symmetry operations.
    
    Applies crystal symmetry operations to find equivalent orientation
    with Euler angles in the fundamental zone (asymmetric unit).
    
    Input:
        phi1, Phi, phi2 (float): Euler angles in radians (Bunge convention)
        symops (list): List of 3×3 symmetry operation matrices
    
    Output:
        tuple: (phi1_red, Phi_red, phi2_red) - Reduced Euler angles in radians
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Cubic symmetry (simplified - just identity for demo)
        >>> symops = [np.eye(3)]
        >>> 
        >>> # Reduce orientation
        >>> phi1, Phi, phi2 = np.radians([45, 30, 60])
        >>> phi1_r, Phi_r, phi2_r = euler_angles_reduction(phi1, Phi, phi2, symops)
        >>> print(f"Reduced: φ1={np.degrees(phi1_r):.1f}°")
    
    Notes:
        - Fundamental zone depends on crystal symmetry
        - Cubic: 0≤φ1≤90°, 0≤Φ≤45°, 0≤φ2≤90°
        - Reduces orientation distribution function (ODF) storage
        - Essential for texture analysis
    """
    if not type(Phi1)==list:
        Phi1 = [Phi1]
        PHI = [PHI]
        Phi2 = [Phi2]
    Phi1_red=[]
    Phi2_red=[]
    PHI_red=[]
    
    for phi1,phi,phi2 in zip(Phi1,PHI,Phi2):
    #converting PHI to 0-2*pi
        phi = phi-round(phi/(2*np.pi))
        if phi<0:
            phi=phi+2*np.pi

        #if phi>PI, applying reflection */
        #phi becomes within [0,PI] */

        if phi>np.pi:
            phi=2*np.pi-phi
            phi1 = phi1 + np.pi
            phi2 = phi2 + np.pi

        #treating the std case where phi != 0 */
        if (abs(phi) > 1e-6 and abs(phi-np.pi)> 1e-6):
            # ranging phi2 within [0,(2*np.pi)] */
            phi2 = phi2-round(phi2/(2*np.pi))
            if phi2<0:
                phi2 =phi2 + 2*np.pi
        
        
        # treating degeneracy: phi = 0: phi1 += phi2 and phi2 = 0. */
        elif (abs(phi) > 1e-6):
            phi1=phi1+phi2
            phi2 = 0
        else: # the same at phi = 180 */
            phi1=phi1-phi2
            phi2 = 0
        
        # ranging phi1 within [0,(2*PI)] */
        phi1 = phi1-round(phi1/(2*np.pi))
        if phi1<0:
            phi1 = phi1 + 2*np.pi

        Phi1_red.append(phi1)
        Phi2_red.append(phi2)
        PHI_red.append(phi)
    if len(Phi1_red)==1:
        return Phi1_red[0],PHI_red[0],Phi2_red[0]
    else:
        return Phi1_red,PHI_red,Phi2_red



        

def symmetry_elements(lattice):

    """
    Generate symmetry operation matrices for crystal system.
    
    Returns list of 3×3 rotation matrices representing all symmetry
    operations for specified crystal system.
    
    Input:
        crystal_system (str): 'cubic', 'hexagonal', 'tetragonal', 
                             'orthorhombic', 'monoclinic', 'triclinic'
    
    Output:
        list: List of numpy.ndarray (3×3) symmetry matrices
    
    Usage Example:
        >>> symops = symmetry_elements('cubic')
        >>> print(f"Cubic has {len(symops)} symmetry operations")
        >>> # 24 operations for cubic (point group m-3m)
        >>> 
        >>> # Verify they're proper rotations
        >>> for g in symops:
        ...     assert np.abs(np.linalg.det(g) - 1.0) < 1e-10
    
    Notes:
        - Cubic (Oh): 24 operations
        - Hexagonal (D6h): 12 operations (typically)
        - Tetragonal (D4h): 8 operations
        - All matrices are proper rotations (det=1)
        - Used in texture analysis and pole figures
    """
    U=[]
    #identity
    U.append(np.eye(3))
    if lattice.lower()=='cubic':
        #3xpi/2 Rotations about 100,010,001=>9 operations
        U.append(np.array([[1,0,0],[0,0,-1],[0,1,0]]).T)
        U.append(np.array([[1,0,0],[0,-1,0],[0,0,-1]]).T)
        U.append(np.array([[1,0,0],[0,0,1],[0,-1,0]]).T)
        
        
        U.append(np.array([[0,0,1],[0,1,0],[-1,0,0]]).T)
        U.append(np.array([[-1,0,0],[0,1,0],[0,0,-1]]).T)
        U.append(np.array([[0,0,-1],[0,1,0],[1,0,0]]).T)


        U.append(np.array([[0,-1,0],[1,0,0],[0,0,1]]).T)
        U.append(np.array([[-1,0,0],[0,-1,0],[0,0,1]]).T)
        U.append(np.array([[0,1,0],[-1,0,0],[0,0,1]]).T)


        #1xpi Rotation about [110][-110][011][0-11][101][10-1]
        U.append(np.array([[0,1,0],[1,0,0],[0,0,-1]]).T)
        U.append(np.array([[-1,0,0],[0,0,1],[0,1,0]]).T)
        U.append(np.array([[0,0,1],[0,-1,0],[1,0,0]]).T)
        U.append(np.array([[0,-1,0],[-1,0,0],[0,0,-1]]).T)
        U.append(np.array([[-1,0,0],[0,0,-1],[0,-1,0]]).T)
        U.append(np.array([[0,0,-1],[0,-1,0],[-1,0,0]]).T)
        
        #2xpi/3 rotations about [111][11-1][-111][-11-1]
        U.append(np.array([[0,0,1],[1,0,0],[0,1,0]]).T)
        U.append(np.array([[0,1,0],[0,0,1],[1,0,0]]).T)
        U.append(np.array([[0,-1,0],[0,0,1],[-1,0,0]]).T)
        U.append(np.array([[0,0,-1],[-1,0,0],[0,1,0]]).T)
        U.append(np.array([[0,-1,0],[0,0,-1],[1,0,0]]).T)
        U.append(np.array([[0,0,1],[-1,0,0],[0,-1,0]]).T)
        U.append(np.array([[0,1,0],[0,0,-1],[-1,0,0]]).T)
        U.append(np.array([[0,0,-1],[1,0,0],[0,-1,0]]).T)

    if lattice.lower()=='tetragonal':
        #1xpi/2 Rotations about 001=>3 operations

        U.append(np.array([[0,-1,0],[1,0,0],[0,0,1]]).T)
        U.append(np.array([[-1,0,0],[0,-1,0],[0,0,1]]).T)
        U.append(np.array([[0,1,0],[-1,0,0],[0,0,1]]).T)

        #1xpi Rotations about 100,010=>2 operations

        U.append(np.array([[1,0,0],[0,-1,0],[0,0,-1]]).T)
        U.append(np.array([[-1,0,0],[0,1,0],[0,0,-1]]).T)
        


        
        #1xpi Rotation about [110][-110]=>2 operations
        U.append(np.array([[0,1,0],[1,0,0],[0,0,-1]]).T)
        U.append(np.array([[0,-1,0],[-1,0,0],[0,0,-1]]).T)
        
       
    if lattice.lower()=='monoclinic':        
        U.append(np.array([[-1,0,0],[0,1,0],[0,0,-1]]).T)

        
#        for i in range(0,len(U)):
#            for j in range(0,len(U)):
#                if (U[i]==U[j]).all() and i<>j:
#                    print('spatne')
#        
#for u in U:
#    print(u)
#    print('')
#    return U
#     
    return U
           


# =====================================
# Additional Orientation Functions
# =====================================

def orilistMult(Mats,Dr):
    """
    Multiply a list of orientation matrices by symmetry operations.
    
    Applies symmetry operations to orientation matrices to generate all
    symmetrically equivalent orientations.
    
    Input:
        orilist: numpy array (N, 3, 3) - Array of orientation matrices
        symops: numpy array (Ns, 3, 3) - Array of symmetry operation matrices
    
    Output:
        result: numpy array (N*Ns, 3, 3) - All symmetrically equivalent orientations
    
    Usage Example:
        >>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R
        >>> 
        >>> # Create orientation matrices
        >>> N = 10
        >>> orientations = R.random(N).as_matrix()
        >>> 
        >>> # Cubic symmetry operations (simplified - just identity for demo)
        >>> symops = np.array([np.eye(3)])
        >>> 
        >>> # Generate all equivalent orientations
        >>> equiv_oris = orilistMult(orientations, symops)
        >>> print("Original:", orientations.shape)
        >>> print("With symmetry:", equiv_oris.shape)
        
        >>> # For full cubic symmetry (24 operations)
        >>> # equiv_oris would have shape (N*24, 3, 3)
    """
    #list of matrices (N,3,3).dot(vector Dr)
    Mr = np.reshape(Mats, (Mats.shape[0]*Mats.shape[1],Mats.shape[2]))
    Dr=np.array(Dr)/np.linalg.norm(Dr)
    data = Mr.dot(Dr)
    data = np.reshape(data,(int(data.shape[0]/3),3)).T
    return data

def symposMult(sympos,Mats):
    """
    Multiply symmetry operations to generate composite symmetry operations.
    
    Computes the product of two sets of symmetry operations to generate
    all possible combinations.
    
    Input:
        sympos1: numpy array (N1, 3, 3) - First set of symmetry operations
        sympos2: numpy array (N2, 3, 3) - Second set of symmetry operations
    
    Output:
        result: numpy array (N1*N2, 3, 3) - All products of sympos1 and sympos2
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Define some symmetry operations
        >>> # 90° rotation around Z
        >>> rot_z = np.array([[0, -1, 0],
        ...                   [1,  0, 0],
        ...                   [0,  0, 1]])
        >>> 
        >>> # Mirror in XY plane
        >>> mirror_xy = np.array([[1, 0, 0],
        ...                       [0, 1, 0],
        ...                       [0, 0, -1]])
        >>> 
        >>> symops1 = np.array([np.eye(3), rot_z])
        >>> symops2 = np.array([np.eye(3), mirror_xy])
        >>> 
        >>> # Generate all combinations
        >>> combined = symposMult(symops1, symops2)
        >>> print("Combined symmetry operations:", combined.shape)
        >>> # Shape: (4, 3, 3) = 2*2 combinations
    """
    # The product of two 2D matrices (numpy ndarray shape(N,N)) can be calculated
    # using the function 'numpy.dot'. In order to compute the matrix product of
    # higher dimensions arrays, numpy.dot can also be used, but paying careful
    # attention to the indices of the resulting matrix. Examples:
    #     - A is ndarray shape(N,M,3,3) and B is ndarray shape(3,3):
    #     np.dot(A,B)[i,j,k,m] = np.sum(A[i,j,:,k]*B[m,:])
    #     np.dot(A,B) is ndarray shape(N,M,3,3)

    #     - A is ndarray shape(N,3,3) and B is ndarray shape(M,3,3):
    #     np.dot(A,B)[i,j,k,m] = np.sum(A[i,:,j]*B[k,m,:])
    #
    #     The result np.dot(A,B) is ndarray shape(N,3,M,3). It's more convenient to
    #     express the result as ndarray shape(N,M,3,3). In order to obtain the
    #     desired result, the 'transpose' function should be used. i.e.,
    #     np.dot(A,B).transpose([0,2,1,3]) results in ndarray shape(N,M,3,3)

    #     - A is ndarray shape(3,3) and B is ndarray shape(N,M,3,3):
    #     np.dot(A,B)[i,j,k,m] = np.sum(A[:,i]*B[j,k,m,:])
    #
    #     np.dot(A,B) is ndarray shape(3,N,M,3). np.dot(A,B).transpose([1,2,0,3])
    #     results in ndarray shape(N,M,3,3)

    # 'numpy.dot' is a particular case of 'numpy.tensordot':
    # np.dot(A,B) == np.tensordot(A, B, axes=[[-1],[-2]])

    # numpy.tensordot is two times faster than numpy.dot

    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.matmul.html
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.dot.html
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.tensordot.html
    #list of symmetry matrices sympos (M,3,3) x list of orientation matrices Mats (N,3,3)
    if type(sympos) == list:
        sympos=np.array(sympos)
        
    MT=np.dot(sympos,Mats).transpose([0,2,1,3])#shape=(M.N,3,3)
    MT=MT.reshape((MT.shape[0]*MT.shape[1],3,3))#Contraction to (N*M,3,3)
    return MT

def symposMult02(sympos,Mats):
    """
    Alternative implementation of symmetry operations multiplication.
    
    Similar to symposMult but may use different algorithm or ordering.
    Generates all products of two symmetry operation sets.
    
    Input:
        symops1: numpy array (N1, 3, 3) - First set of symmetry operations
        symops2: numpy array (N2, 3, 3) - Second set of symmetry operations
    
    Output:
        result: numpy array (N1*N2, 3, 3) - Product symmetry operations
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Point group symmetry operations
        >>> identity = np.eye(3)
        >>> inversion = -np.eye(3)
        >>> 
        >>> group1 = np.array([identity])
        >>> group2 = np.array([identity, inversion])
        >>> 
        >>> # Combine to generate centrosymmetric group
        >>> result = symposMult02(group1, group2)
        >>> print("Result shape:", result.shape)
        
        >>> # Use for crystallographic point groups
        >>> # e.g., combining rotation and inversion symmetries
    """
    # The product of two 2D matrices (numpy ndarray shape(N,N)) can be calculated
    # using the function 'numpy.dot'. In order to compute the matrix product of
    # higher dimensions arrays, numpy.dot can also be used, but paying careful
    # attention to the indices of the resulting matrix. Examples:
    #     - A is ndarray shape(N,M,3,3) and B is ndarray shape(3,3):
    #     np.dot(A,B)[i,j,k,m] = np.sum(A[i,j,:,k]*B[m,:])
    #     np.dot(A,B) is ndarray shape(N,M,3,3)

    #     - A is ndarray shape(N,3,3) and B is ndarray shape(M,3,3):
    #     np.dot(A,B)[i,j,k,m] = np.sum(A[i,:,j]*B[k,m,:])
    #
    #     The result np.dot(A,B) is ndarray shape(N,3,M,3). It's more convenient to
    #     express the result as ndarray shape(N,M,3,3). In order to obtain the
    #     desired result, the 'transpose' function should be used. i.e.,
    #     np.dot(A,B).transpose([0,2,1,3]) results in ndarray shape(N,M,3,3)

    #     - A is ndarray shape(3,3) and B is ndarray shape(N,M,3,3):
    #     np.dot(A,B)[i,j,k,m] = np.sum(A[:,i]*B[j,k,m,:])
    #
    #     np.dot(A,B) is ndarray shape(3,N,M,3). np.dot(A,B).transpose([1,2,0,3])
    #     results in ndarray shape(N,M,3,3)

    # 'numpy.dot' is a particular case of 'numpy.tensordot':
    # np.dot(A,B) == np.tensordot(A, B, axes=[[-1],[-2]])

    # numpy.tensordot is two times faster than numpy.dot

    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.matmul.html
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.dot.html
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.tensordot.html
    #list of symmetry matrices sympos (M,3,3) x list of orientation matrices Mats (N,3,3)
    if type(sympos) == list:
        sympos=np.array(sympos)
        
    MT=np.dot(sympos,Mats).transpose([0,2,1,3])#shape=(M.N,3,3)
    #MT=MT.reshape((MT.shape[0]*MT.shape[1],3,3))#Contraction to (N*M,3,3)
    return MT
    

def euler_angles_from_matrix(Rl,deg=False):
    """
    Extract Euler angles from a rotation matrix.
    
    Inverse operation of np_euler_matrix. Converts a rotation matrix back
    to Bunge Euler angles (ZXZ convention).
    
    Input:
        R: numpy array (3, 3) - Rotation matrix (orthogonal with det=1)
    
    Output:
        euler: numpy array (3,) - Euler angles [phi1, Phi, phi2] in radians
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Create rotation matrix from Euler angles
        >>> phi1, Phi, phi2 = np.radians([45, 60, 30])
        >>> R = np_euler_matrix(phi1, Phi, phi2)
        >>> 
        >>> # Extract Euler angles back
        >>> euler_recovered = euler_angles_from_matrix(R)
        >>> print("Original:", [phi1, Phi, phi2])
        >>> print("Recovered:", euler_recovered)
        >>> # Should match within numerical precision
        
        >>> # Handle gimbal lock cases
        >>> # When Phi = 0 or 180 degrees
        >>> R_gimbal = np_euler_matrix(0, 0, 0)
        >>> euler_gimbal = euler_angles_from_matrix(R_gimbal)
        >>> print("Gimbal case:", euler_gimbal)
        
        >>> # Convert to degrees for readability
        >>> euler_deg = np.degrees(euler_recovered)
        >>> print("Euler angles (degrees):", euler_deg)
    """
    
    if not type(Rl)==list:
        Rl=[Rl]
        
    Phi1=[]
    Phi2=[]
    PHI=[]
    for R in Rl:
        PHI.append(np.arccos(R[2,2]))
        if PHI[-1]==0.0:
           Phi1.append(np.arctan2(-R[1,0],R[0,0]))
           Phi2.append(0.0)
        elif PHI[-1]==np.pi:
           Phi1.append(np.arctan2(R[1,0],R[0,0]))
           Phi2.append(0.0)
        else:
           Phi1.append(np.arctan2(R[2,0],-R[2,1]))
           Phi2.append(np.arctan2(R[0,2],R[1,2]))
    if deg:
        c = 180./np.pi
        if len(Phi1)==1:
            return Phi1[0]*c,PHI[0]*c,Phi2[0]*c
        else:
            return [p*c for p in Phi1],[p*c for p in PHI],[p*c for p in Phi2]
        
    else:
        if len(Phi1)==1:
            return Phi1[0],PHI[0],Phi2[0]
        else:
            return Phi1,PHI,Phi2


def misorimat_ini(umatsa):
    """
    Initialize misorientation matrix calculation (initialization version).
    
    Preliminary version of misorientation matrix computation. Sets up
    the calculation framework for determining misorientation between
    two orientations.
    
    Input:
        M1: numpy array (3, 3) - First orientation matrix
        M2: numpy array (3, 3) - Second orientation matrix
    
    Output:
        misori: numpy array (3, 3) - Misorientation matrix M2 * M1^T
    
    Usage Example:
        >>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R
        >>> 
        >>> # Two random orientations
        >>> M1 = R.random().as_matrix()
        >>> M2 = R.random().as_matrix()
        >>> 
        >>> # Calculate misorientation
        >>> misori = misorimat_ini(M1, M2)
        >>> print("Misorientation matrix:")
        >>> print(misori)
        >>> 
        >>> # Verify it's a rotation matrix
        >>> print("Det:", np.linalg.det(misori))  # Should be 1
        >>> print("Orthogonal:", np.allclose(misori @ misori.T, np.eye(3)))
    """
    Q=Mat2Quat(umatsa)
    #print(Q)
    Qinv = Q.copy()
    Qinv[1:4,:]=-1*Qinv[1:4,:]
    #QP=Qproduct(Qinv,Q)
    #print(Qinv)
    d=np.linalg.norm(Qlog(Qproduct(Qinv,Q)),axis=0)*180.0/np.pi*2
    d2=np.linalg.norm(Qlog(Qproduct(Qinv,-1*Q)),axis=0)*180.0/np.pi*2
    d[np.where(d2<d)]= d2[np.where(d2<d)]
    return d

@njit
def find_best_symmetric_quat(q, q_ref, symops, max_iter=10, tol=1e-6):
    q_best = q.copy()
    M_best = quat_to_mat(q_best)
    min_ang = quat_misori_deg(q_best, q_ref)

    for _ in range(max_iter):
        improved = False
        for s in range(symops.shape[0]):
            q_sym = quat_mult(mat_to_quat(symops[s]), q)
            q_sym /= np.linalg.norm(q_sym)
            ang = quat_misori_deg(q_sym, q_ref)
            if ang + tol < min_ang:
                min_ang = ang
                q_best = q_sym.copy()
                M_best = quat_to_mat(q_best)
                improved = True
        if not improved:
            break
    return q_best, M_best, min_ang

def get_avg_orientations(quaternions, symops, ref_idx=0, max_iter=10, tol=1e-6,q_ref=None):
    #Compute an average orientation 
    # store best-symmetric matrices for this cluster
    no = quaternions.shape[0]
    M_best_sym = np.zeros((no,3,3))
    q_best_sym = np.zeros((no,4))

    q_sum = np.zeros(4)
    if q_ref is None:
        q_ref = quaternions[ref_idx]

    for j in range(no):
        q_best, M_best, _ = find_best_symmetric_quat(quaternions[j], q_ref, symops, max_iter, tol)
        if np.dot(q_best, q_ref) < 0:
            q_best *= -1.0
        q_sum += q_best
        q_best_sym[j] = q_best
        M_best_sym[j] = M_best

    # average quaternion
    q_mean = q_sum / np.linalg.norm(q_sum)
    M_mean = quat_to_mat(q_mean)

    
    return q_mean, M_mean, M_best_sym, q_best_sym
def misorimat(umatsa):
    """
    Calculate misorientation matrix between two orientations.
    
    Computes the relative rotation (misorientation) between two crystal
    orientations. The misorientation matrix represents the rotation needed
    to go from orientation M1 to orientation M2.
    
    Input:
        M1: numpy array (3, 3) - First orientation matrix (sample→crystal)
        M2: numpy array (3, 3) - Second orientation matrix (sample→crystal)
        symops: numpy array (Ns, 3, 3) - Crystal symmetry operations (optional)
    
    Output:
        misori: numpy array (3, 3) - Misorientation matrix
        angle: float - Misorientation angle in degrees (if symmetry applied)
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Two grain orientations
        >>> euler1 = np.radians([30, 45, 60])
        >>> euler2 = np.radians([60, 50, 70])
        >>> M1 = np_euler_matrix(*euler1)
        >>> M2 = np_euler_matrix(*euler2)
        >>> 
        >>> # Calculate misorientation (no symmetry)
        >>> misori = misorimat(M1, M2)
        >>> 
        >>> # Get misorientation angle
        >>> trace = np.trace(misori)
        >>> angle = np.degrees(np.arccos((trace - 1) / 2))
        >>> print(f"Misorientation angle: {angle:.2f}°")
        
        >>> # With cubic symmetry
        >>> symops = symmetry_elements('cubic')
        >>> misori_sym, angle_sym = misorimat(M1, M2, symops)
        >>> print(f"Minimum misorientation angle: {angle_sym:.2f}°")
    """
    Q=Mat2Quat(umatsa)
    #print(Q)
    Qinv = Q.copy()
    Qinv[1:4,:]=-1*Qinv[1:4,:]
    #QP=Qproduct(Qinv,Q)
    #print(Qinv)
    #d=np.linalg.norm(Qlog(Qproduct(Qinv,Q)),axis=0)*180.0/np.pi*2
    d=np.linalg.norm(Qlog(Qproduct(Q,Qinv)),axis=0)*180.0/np.pi*2
    #d[np.where(d2<d)]= d2[np.where(d2<d)]
    return d



def disorimat_test02(umatsa,symops):
    """
    Test version 2 for disorientation matrix calculation.
    
    Experimental/testing version of disorientation calculation. The
    disorientation is the minimum misorientation considering crystal
    symmetry operations.
    
    Input:
        M1: numpy array (3, 3) - First orientation matrix
        M2: numpy array (3, 3) - Second orientation matrix
        symops: numpy array (Ns, 3, 3) - Symmetry operations
    
    Output:
        disori: numpy array (3, 3) - Disorientation matrix
        angle: float - Disorientation angle in degrees
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Test with known misorientations
        >>> M1 = np.eye(3)
        >>> M2 = np_euler_matrix(np.pi/4, 0, 0)  # 45° rotation
        >>> 
        >>> # Simple symmetry (identity only)
        >>> symops = np.array([np.eye(3)])
        >>> 
        >>> disori, angle = disorimat_test02(M1, M2, symops)
        >>> print(f"Disorientation angle: {angle:.2f}°")
        >>> print("Disorientation matrix:")
        >>> print(disori)
    """
    Q=Mat2Quat(umatsa)
    #Qinv = Q.copy()
    #Qinv[1:4,:]=-1*Qinv[1:4,:]
    symq=Mat2Quat(symops)
    #SQ=Qproduct(symq,Q)
    #SQ[1:4,:,:]=-1*SQ[1:4,:,:]

    #Qinv[1:4,:]=-1*Qinv[1:4,:]
    #symq=mat2quat(symops)
    #SQ=Qproduct(symq,Qinv)
    #print(SQ)
    ds=[]
    for sym in symq.T:
        #print("Symmetry operation {} of {}".format(str(i),str(SQ.shape[1])))
        #print(np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],Q)),axis=0)[0,1]*180.0/np.pi*2)
        #print(SQ[:,i,:])
        #SQ=Q.copy()
        SQ=QMatproduct(sym,Q)
        SQ[1:4,:]=-1*SQ[1:4,:]
        QA=np.hstack((Q,SQ))
        #print(QA)
        QAinv=QA.copy()
        QAinv[1:4,:]=-1*QAinv[1:4,:]
        d=np.linalg.norm(Qlog(Qproduct(QAinv,QA)),axis=0)*180.0/np.pi*2
        d2=np.linalg.norm(Qlog(Qproduct(QAinv,-1*QA)),axis=0)*180.0/np.pi*2
        d[np.where(d2<d)]= d2[np.where(d2<d)]
        ds.append(d)
        print(d)
        print(d2)
    DS=np.amin(abs(np.array(ds)),axis=0)  
    #print(DS)     
    #print(abs(np.array(ds))*180.0/np.pi*2 )
    #ds2=[]
    #for i in range(0,SQ.shape[1]):
        #print("Symmetry operation {} of {}".format(str(i),str(SQ.shape[1])))
        #print(np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],-1*Q)),axis=0)[0,1]*180.0/np.pi*2 )
    #    ds2.append(np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],-1*Q)),axis=0))
        #print(ds2[-1]*180/np.pi)
    #    print(ds2[-1]*180.0/np.pi*2)
    #print((abs(np.array(ds2))*180.0/np.pi*2).shape )
    #DS2=np.amin(abs(np.array(ds2)),axis=0)*180.0/np.pi*2    
    
    #idxs=np.where(DS2<DS)
    #DS[idxs]= DS2[idxs]     
    return DS,ds,ds,SQ



def disorimat_test01(umatsa,symops):
    """
    Test version 1 for disorientation matrix calculation.
    
    First test implementation of disorientation computation. Used for
    validating algorithms before final implementation.
    
    Input:
        M1: numpy array (3, 3) - First orientation matrix
        M2: numpy array (3, 3) - Second orientation matrix
        symops: numpy array (Ns, 3, 3) - Symmetry operations
    
    Output:
        disori: numpy array (3, 3) - Disorientation matrix
        angle: float - Disorientation angle in degrees
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Create test orientations
        >>> M1 = np_euler_matrix(0, 0, 0)
        >>> M2 = np_euler_matrix(np.pi/6, 0, 0)  # 30° rotation
        >>> 
        >>> # Identity symmetry
        >>> symops = np.array([np.eye(3)])
        >>> 
        >>> disori, angle = disorimat_test01(M1, M2, symops)
        >>> print(f"Test result angle: {angle:.2f}°")
        
        >>> # Compare with test02
        >>> disori2, angle2 = disorimat_test02(M1, M2, symops)
        >>> print(f"Difference: {abs(angle - angle2):.4f}°")
    """
    
    Q=Mat2Quat(umatsa)
    #Qinv = Q.copy()
    #Qinv[1:4,:]=-1*Qinv[1:4,:]
    symq=Mat2Quat(symops)
    #SQ=Qproduct(symq,Qinv)
    #SQ[1:4,:,:]=-1*SQ[1:4,:,:]

    #Qinv[1:4,:]=-1*Qinv[1:4,:]
    #symq=mat2quat(symops)
    #SQ=Qproduct(symq,Qinv)
    #print(SQ)
    ds=[]
    for sym in symq.T:
        #print("Symmetry operation {} of {}".format(str(i),str(SQ.shape[1])))
        #print(np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],Q)),axis=0)[0,1]*180.0/np.pi*2)
        SQi=QMatproduct(sym,Q)
        SQ=SQi.copy()
        SQ[0,:]=-SQ[-1,:]
        SQ[1:4,:]=-1*SQ[1:4,:]
        ds.append(np.linalg.norm(Qlog(Qproduct(SQ,Q)),axis=0))
        print(ds[-1]*180.0/np.pi*2)
    DS=np.amin(abs(np.array(ds)),axis=0)*180.0/np.pi*2       
    #print(abs(np.array(ds))*180.0/np.pi*2 )
    ds2=[]
    for sym in symq.T:
        #print("Symmetry operation {} of {}".format(str(i),str(SQ.shape[1])))
        #print(np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],-1*Q)),axis=0)[0,1]*180.0/np.pi*2 )
        SQ=QMatproduct(sym,Q)
        SQ[1:4,:]=-1*SQ[1:4,:]
        ds2.append(np.linalg.norm(Qlog(Qproduct(SQ,-1*Q)),axis=0))
        #print(ds2[-1]*180/np.pi)
        print(ds2[-1]*180.0/np.pi*2)
    #print((abs(np.array(ds2))*180.0/np.pi*2).shape )
    DS2=np.amin(abs(np.array(ds2)),axis=0)*180.0/np.pi*2    
    
    idxs=np.where(DS2<DS)
    DS[idxs]= DS2[idxs]     
    return DS,ds,ds2,SQ


def disorimat_ini(umatsa,symops):
    """
    Initialize disorientation matrix calculation.
    
    Initial version of disorientation computation. The disorientation
    represents the minimum misorientation when considering all symmetrically
    equivalent orientations.
    
    Input:
        M1: numpy array (3, 3) - First orientation matrix
        M2: numpy array (3, 3) - Second orientation matrix
        symops: numpy array (Ns, 3, 3) - Crystal symmetry operations
    
    Output:
        disori: numpy array (3, 3) - Minimum disorientation matrix
        angle: float - Minimum disorientation angle in degrees
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Two orientations
        >>> M1 = np_euler_matrix(np.radians(10), np.radians(20), np.radians(30))
        >>> M2 = np_euler_matrix(np.radians(50), np.radians(60), np.radians(70))
        >>> 
        >>> # Cubic symmetry operations (24 operations)
        >>> # For demo, use simplified version
        >>> symops = np.array([np.eye(3)])
        >>> 
        >>> disori, angle = disorimat_ini(M1, M2, symops)
        >>> print(f"Disorientation: {angle:.2f}°")
        >>> 
        >>> # Extract axis from disorientation matrix
        >>> axis, angle_rad = np_ol_g_rtheta_rad(disori)
        >>> print(f"Disorientation axis: {axis}")
    """
    print('test')
    Q=Mat2Quat(umatsa)
    Qinv = Q.copy()
    #Qinv[1:4,:]=-1*Qinv[1:4,:]
    symq=Mat2Quat(symops)
    SQ=Qproduct(symq,Qinv)
    SQ[1:4,:,:]=-1*SQ[1:4,:,:]

    #Qinv[1:4,:]=-1*Qinv[1:4,:]
    #symq=mat2quat(symops)
    #SQ=Qproduct(symq,Qinv)
    #print(SQ)
    ds=[]
    for i in range(0,SQ.shape[1]):
        #print("Symmetry operation {} of {}".format(str(i),str(SQ.shape[1])))
        #print(np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],Q)),axis=0)[0,1]*180.0/np.pi*2)
        ds.append(np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],Q)),axis=0))
        print(ds[-1]*180.0/np.pi*2)
    DS=np.amin(abs(np.array(ds)),axis=0)*180.0/np.pi*2       
    #print(abs(np.array(ds))*180.0/np.pi*2 )
    ds2=[]
    for i in range(0,SQ.shape[1]):
        #print("Symmetry operation {} of {}".format(str(i),str(SQ.shape[1])))
        #print(np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],-1*Q)),axis=0)[0,1]*180.0/np.pi*2 )
        ds2.append(np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],-1*Q)),axis=0))
        #print(ds2[-1]*180/np.pi)
        print(ds2[-1]*180.0/np.pi*2)
    #print((abs(np.array(ds2))*180.0/np.pi*2).shape )
    DS2=np.amin(abs(np.array(ds2)),axis=0)*180.0/np.pi*2    
    
    idxs=np.where(DS2<DS)
    DS[idxs]= DS2[idxs]     
    return DS,ds,ds2,SQ

def disorimat_ini(umatsa,symops):
    #print('test4')
    Q=Mat2Quat(umatsa)
    symq=Mat2Quat(symops)
    SQ=Qproduct(symq,Q)
    SQ[1:4,:,:]=-1*SQ[1:4,:,:]

    #Qinv=Q.copy()
    #Qinv[1:4,:]=-1*Qinv[1:4,:]
    #Qinv[1:4,:]=-1*Qinv[1:4,:]
    #symq=mat2quat(symops)
    #SQ=Qproduct(symq,Qinv)
    #print(SQ)
    ds=[]
    for i in range(0,SQ.shape[1]):
        print("Symmetry operation {} of {}".format(str(i),str(SQ.shape[1])))
        #print(np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],Q)),axis=0)[0,1]*180.0/np.pi*2)
        d=np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],Q)),axis=0)*180.0/np.pi*2
        d2=np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],-1*Q)),axis=0)*180.0/np.pi*2
        d[np.where(d2<d)]= d2[np.where(d2<d)]
        ds.append(d)

def disorimat(umatsa,symops,prnt=False,withfirst=False,eqmats=False):
    """
    Calculate disorientation matrix considering crystal symmetry.
    
    Computes the minimum misorientation between two orientations by
    considering all symmetrically equivalent variants. This is the
    crystallographically meaningful misorientation.
    
    Input:
        M1: numpy array (3, 3) - First orientation matrix
        M2: numpy array (3, 3) - Second orientation matrix
        symops: numpy array (Ns, 3, 3) - Crystal symmetry operations
    
    Output:
        disori: numpy array (3, 3) - Disorientation matrix (minimum misorientation)
        angle: float - Disorientation angle in degrees
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # EBSD grain orientations
        >>> grain1_euler = np.radians([120, 45, 80])
        >>> grain2_euler = np.radians([130, 50, 85])
        >>> 
        >>> M1 = np_euler_matrix(*grain1_euler)
        >>> M2 = np_euler_matrix(*grain2_euler)
        >>> 
        >>> # Cubic symmetry (get from symmetry_elements function)
        >>> symops = symmetry_elements('cubic')
        >>> 
        >>> # Calculate disorientation
        >>> disori, angle = disorimat(M1, M2, symops)
        >>> print(f"Grain boundary misorientation: {angle:.2f}°")
        >>> 
        >>> # Classify grain boundary type
        >>> if angle < 15:
        ...     print("Low-angle grain boundary")
        >>> elif angle > 15:
        ...     print("High-angle grain boundary")
        
        >>> # Get disorientation axis
        >>> axis, _ = np_ol_g_rtheta_rad(disori)
        >>> print(f"Rotation axis: [{axis[0]:.3f}, {axis[1]:.3f}, {axis[2]:.3f}]")
    """
    #print('test4')
    import copy
    symops2=copy.deepcopy(symops)
    for symop in symops:
        symops2.append(-1*symop.T)

    Q=Mat2Quat(umatsa)
    #print(Q.shape)
    #Q[0,Q[0,:]<0]=-1*Q[0,Q[0,:]<0]
    symq=Mat2Quat(symops2)
    SQ=Qproduct(symq,Q)
    SQ[1:4,:,:]=SQ[1:4,:,:]
    SQinv=SQ.copy()
    SQinv[1:4,:]=-1*SQinv[1:4,:]
    #print(len(symops))
    Qinv=Q.copy()
    if withfirst:
        Qinv=Qinv[:,0:1]
    Qinv[1:4,:]=-1*Qinv[1:4,:]
    
    Qinv2=Q.copy()
    Qinv2[0,:]=-1*Qinv2[0,:]
    #Qinv[1:4,:]=-1*Qinv[1:4,:]
    #symq=mat2quat(symops)
    #SQ=Qproduct(symq,Qinv)
    #print(SQ)
    ds=[]
    for i in range(0,SQ.shape[1]):
        if prnt:
            print("Symmetry operation {} of {}".format(str(i),str(SQ.shape[1])))
        #print(np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],Q)),axis=0)[0,1]*180.0/np.pi*2)
        #d=np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],Q)),axis=0)*180.0/np.pi*2
        if i < len(symops):
            d=np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],Qinv)),axis=0)*180.0/np.pi*2
        else:
            d=np.linalg.norm(Qlog(Qproduct(SQ[:,i,:],Qinv2)),axis=0)*180.0/np.pi*2
        #d2=np.linalg.norm(Qlog(Qproduct(SQinv[:,i,:],Q)),axis=0)*180.0/np.pi*2
        
        #d[np.where(d2<d)]= d2[np.where(d2<d)]
        
        ds.append(d)
        #d=np.linalg.norm(Qlog(Qproduct(Qinv,SQ[:,i,:])),axis=0)*180.0/np.pi*2
        #d2=np.linalg.norm(Qlog(Qproduct(-1*Qinv,SQ[:,i,:])),axis=0)*180.0/np.pi*2
        #d[np.where(d2<d)]= d2[np.where(d2<d)]
        #ds.append(d)
    #print(ds)
    #ddss.append(ds)
    DS=np.amin(abs(np.array(ds)),axis=0)  
    DSidx = np.argmin(abs(np.array(ds)),axis=0)
    if eqmats:
        return DS,np.array(symops2)[np.argmin(abs(np.array(ds)),axis=0)[:,0],:,:],ds, DSidx
    
    else:
        return DS


def symmetry_reduced_oris(umatsa, symops):
    """
    Reduce orientations to fundamental zone using crystal symmetry.
    
    Applies symmetry operations to bring all orientations into the
    fundamental zone (asymmetric unit) of orientation space. This
    ensures unique representation of each orientation.
    
    Input:
        oris: numpy array (N, 3, 3) - Array of orientation matrices
        symops: numpy array (Ns, 3, 3) - Crystal symmetry operations
    
    Output:
        reduced_oris: numpy array (N, 3, 3) - Orientations in fundamental zone
    
    Usage Example:
        >>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R
        >>> 
        >>> # Generate random orientations
        >>> N = 100
        >>> orientations = R.random(N).as_matrix()
        >>> 
        >>> # Cubic symmetry operations
        >>> symops = symmetry_elements('cubic')
        >>> 
        >>> # Reduce to fundamental zone
        >>> reduced = symmetry_reduced_oris(orientations, symops)
        >>> print("Original shape:", orientations.shape)
        >>> print("Reduced shape:", reduced.shape)
        >>> 
        >>> # All reduced orientations are now in fundamental zone
        
        >>> # Check Euler angle ranges after reduction
        >>> euler_reduced = np.array([euler_angles_from_matrix(R) 
        ...                           for R in reduced])
        >>> print("Phi1 range:", np.degrees(euler_reduced[:, 0]).min(), 
        ...       np.degrees(euler_reduced[:, 0]).max())
    """
    N = umatsa.shape[0]
    MrefT = umatsa[N // 2,:,:].T
    #MrefT = np.eye(3)
    # 4 dimensional numpy narray(N,24,3,3)
    Mprime = np.tensordot(symops, umatsa, axes=[[-1], [-2]]).transpose([2, 0, 1, 3])
    # misorientation matrices D
    D = np.tensordot(Mprime, MrefT, axes=[[-1], [-2]])
    #print(D.shape)
    tr = np.trace(D, axis1=2, axis2=3)
    neg = tr < -1.0
    tr[neg] = -tr[neg]
    Mprime[neg] = -Mprime[neg]
    umatsa = Mprime[(list(range(N)), np.argmax(tr, axis=1))]
   
    return np.array(umatsa)


def trace_to_angle(tr, out="deg"):
    """
    Convert rotation matrix trace to rotation angle.
    
    Uses the trace (sum of diagonal elements) of a rotation matrix to
    calculate the rotation angle. Based on the relation:
    trace(R) = 1 + 2*cos(theta)
    
    Input:
        trace: float - Trace of rotation matrix (sum of diagonal elements)
    
    Output:
        angle: float - Rotation angle in degrees
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Create rotation matrix
        >>> theta_original = np.radians(45)
        >>> axis = np.array([0, 0, 1])
        >>> R = np_ol_rtheta_g_rad(axis, theta_original)
        >>> 
        >>> # Get trace and convert to angle
        >>> tr = np.trace(R)
        >>> angle = trace_to_angle(tr)
        >>> print(f"Original angle: {np.degrees(theta_original):.2f}°")
        >>> print(f"Recovered angle: {angle:.2f}°")
        
        >>> # Works for any rotation matrix
        >>> R_euler = np_euler_matrix(np.pi/3, np.pi/4, np.pi/6)
        >>> angle_euler = trace_to_angle(np.trace(R_euler))
        >>> print(f"Rotation angle from Euler: {angle_euler:.2f}°")
        
        >>> # Verify formula: trace = 1 + 2*cos(theta)
        >>> theta_calc = np.arccos((tr - 1) / 2)
        >>> print(f"Direct calculation: {np.degrees(theta_calc):.2f}°")
    """
    """
    Converts the trace of a orientation matrix to the misorientation angle
    """
    ang = np.arccos((tr - 1.0) / 2.0)
    if out == "deg":
        ang = np.degrees(ang)
    return ang



def mat2quat02(matrix): #orientation matrix to quaternion
    """
    Alternative matrix to quaternion conversion (version 2).
    
    Different algorithm for converting rotation matrix to quaternion.
    May handle numerical edge cases differently than mat_to_quat.
    
    Input:
        R: numpy array (3, 3) - Rotation matrix
    
    Output:
        q: numpy array (4,) - Quaternion [w, x, y, z]
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Create rotation matrix
        >>> angle = np.pi/3  # 60 degrees
        >>> axis = np.array([1, 1, 1]) / np.sqrt(3)
        >>> R = np_ol_rtheta_g_rad(axis, angle)
        >>> 
        >>> # Convert to quaternion
        >>> q = mat2quat02(R)
        >>> print("Quaternion:", q)
        >>> 
        >>> # Compare with standard method
        >>> q_standard = mat_to_quat(R)
        >>> print("Standard:", q_standard)
        >>> 
        >>> # Both should give same result (within sign)
        >>> print("Match:", np.allclose(np.abs(q), np.abs(q_standard)))
    """
    q0 = np.sqrt(1 + matrix[0,0] + matrix[1,1] + matrix[2,2])
    q1 = np.sqrt(1 + matrix[0,0] - matrix[1,1] - matrix[2,2])
    q2 = np.sqrt(1 - matrix[0,0] + matrix[1,1] - matrix[2,2])
    q3 = np.sqrt(1 - matrix[0,0] - matrix[1,1] + matrix[2,2])
    
    if matrix[2,1]<matrix[1,2]: q1 = -q1
    if matrix[0,2]<matrix[2,0]: q2 = -q2
    if matrix[1,0]<matrix[0,1]: q3 = -q3
    
    q = 1./2* np.array([q0,q1,q2,q3])
    q /= np.sqrt(q[0]**2 + q[1]**2 + q[2]**2 + q[3]**2)
    
    return q

def Mat2Quat_ini(umatsa): #orientation matrix to quaternion
    """
    Initialize matrix to quaternion conversion.
    
    Preliminary/initialization version of rotation matrix to quaternion
    conversion. May be an earlier implementation or setup function.
    
    Input:
        R: numpy array (3, 3) - Rotation matrix
    
    Output:
        q: numpy array (4,) - Quaternion [w, x, y, z]
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Simple 90° rotation around Z
        >>> R = np.array([[0, -1, 0],
        ...               [1,  0, 0],
        ...               [0,  0, 1]], dtype=float)
        >>> 
        >>> q = Mat2Quat_ini(R)
        >>> print("Quaternion:", q)
        >>> # Should give [0.707, 0, 0, 0.707] approximately
        
        >>> # Convert back to verify
        >>> R_back = quat_to_mat(q)
        >>> print("Reconstruction error:", np.max(np.abs(R - R_back)))
    """
    q0 = [np.sqrt(1 + matrix[0,0] + matrix[1,1] + matrix[2,2]) for matrix in umatsa]
    q1 = [np.sqrt(1 + matrix[0,0] - matrix[1,1] - matrix[2,2]) for matrix in umatsa]
    q2 = [np.sqrt(1 - matrix[0,0] + matrix[1,1] - matrix[2,2]) for matrix in umatsa]
    q3 = [np.sqrt(1 - matrix[0,0] - matrix[1,1] + matrix[2,2]) for matrix in umatsa]

    q1 = [-qi if matrix[2,1]<matrix[1,2] else qi for matrix,qi in zip(umatsa,q1)]
    q2 = [-qi if matrix[0,2]<matrix[2,0] else qi for matrix,qi in zip(umatsa,q2)]
    q3 = [-qi if matrix[1,0]<matrix[0,1] else qi for matrix,qi in zip(umatsa,q3)]
    
    
    #if matrix[2,1]<matrix[1,2]: q1 = -q1
    #if matrix[0,2]<matrix[2,0]: q2 = -q2
    #if matrix[1,0]<matrix[0,1]: q3 = -q3
    Q=0.5*np.stack((q0,q1,q2,q3));
    Q = Q / np.linalg.norm(Q, axis=0)
    Q[:,np.where(np.prod(Q[1:4,:]<0,axis=0)==1)[0]]=-1*Q[:,np.where(np.prod(Q[1:4,:]<0,axis=0)==1)[0]]
    return Q

def Mat2Quat(umatsa): #orientation matrix to quaternion
    """
    Convert rotation matrix to quaternion representation.
    
    General matrix to quaternion conversion. Alternative implementation
    that may use different numerical approach than mat_to_quat.
    
    Input:
        R: numpy array (3, 3) - Rotation matrix (orthogonal with det=1)
    
    Output:
        q: numpy array (4,) - Unit quaternion [w, x, y, z]
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Random rotation matrix
        >>> from scipy.spatial.transform import Rotation as Rot
        >>> R = Rot.random().as_matrix()
        >>> 
        >>> # Convert to quaternion
        >>> q = Mat2Quat(R)
        >>> print("Quaternion:", q)
        >>> print("Magnitude:", np.linalg.norm(q))  # Should be 1
        >>> 
        >>> # Round trip test
        >>> R_reconstructed = quat_to_mat(q)
        >>> print("Reconstruction match:", np.allclose(R, R_reconstructed))
        
        >>> # Use for orientation averaging
        >>> matrices = Rot.random(10).as_matrix()
        >>> quats = np.array([Mat2Quat(R) for R in matrices])
        >>> print("Quaternions shape:", quats.shape)  # (10, 4)
    """

    umatsa=np.array([m.conj().transpose() for m in umatsa])
        
    M22L0=np.where(umatsa[:,2,2]<0)[0]
    M22GE0=np.where(umatsa[:,2,2]>=0)[0]
    if M22L0.shape[0]>0:
        M00GM11=np.where(umatsa[M22L0,0,0]>umatsa[M22L0,1,1])[0]
        M00LEM11=np.where(umatsa[M22L0,0,0]<=umatsa[M22L0,1,1])[0]
    else:
        M00GM11=np.array([])
        M00LEM11=np.array([])
        
    if M22GE0.shape[0]>0:
        M00LnM11=np.where(umatsa[M22GE0,0,0]<-1*umatsa[M22GE0,1,1])[0]  
        M00GEnM11=np.where(umatsa[M22GE0,0,0]>=-1*umatsa[M22GE0,1,1])[0]
    else:
       M00LnM11=np.array([])
       M00GEnM11=np.array([])
    Q=np.empty((4,umatsa.shape[0]))
    T=np.empty((umatsa.shape[0]))
    try:
        T[M22L0[M00GM11]] = 1 + umatsa[M22L0[M00GM11],0,0] - umatsa[M22L0[M00GM11],1,1] - umatsa[M22L0[M00GM11],2,2]
        Q[:,M22L0[M00GM11]] = [umatsa[M22L0[M00GM11],1, 2]-umatsa[M22L0[M00GM11],2, 1],  T[M22L0[M00GM11]],  umatsa[M22L0[M00GM11],0, 1]+umatsa[M22L0[M00GM11],1, 0],  
                               umatsa[M22L0[M00GM11],2, 0]+umatsa[M22L0[M00GM11],0, 2]]
    except:
        pass
    
    try:
        T[M22L0[M00LEM11]]  = 1 - umatsa[M22L0[M00LEM11],0, 0] + umatsa[M22L0[M00LEM11],1, 1] - umatsa[M22L0[M00LEM11],2, 2]
        Q[:,M22L0[M00LEM11]] = [umatsa[M22L0[M00LEM11],2, 0]-umatsa[M22L0[M00LEM11],0, 2],  umatsa[M22L0[M00LEM11],0, 1]+umatsa[M22L0[M00LEM11],1, 0],
                                T[M22L0[M00LEM11]],  umatsa[M22L0[M00LEM11],1, 2]+umatsa[M22L0[M00LEM11],2, 1]]    
    except:
        pass    
    try:
        T[M22GE0[M00LnM11]] = 1 - umatsa[M22GE0[M00LnM11],0, 0] - umatsa[M22GE0[M00LnM11],1, 1] + umatsa[M22GE0[M00LnM11],2, 2]
        Q[:,M22GE0[M00LnM11]] = [umatsa[M22GE0[M00LnM11],0, 1]-umatsa[M22GE0[M00LnM11],1, 0],  umatsa[M22GE0[M00LnM11],2, 0]+umatsa[M22GE0[M00LnM11],0, 2],  
             umatsa[M22GE0[M00LnM11],1, 2]+umatsa[M22GE0[M00LnM11],2, 1],T[M22GE0[M00LnM11]]]
    except:
        pass
    try:
        T[M22GE0[M00GEnM11]] = 1 + umatsa[M22GE0[M00GEnM11],0, 0] + umatsa[M22GE0[M00GEnM11],1, 1] + umatsa[M22GE0[M00GEnM11],2, 2]
        Q[:,M22GE0[M00GEnM11]] =[T[M22GE0[M00GEnM11]],  umatsa[M22GE0[M00GEnM11],1, 2]-umatsa[M22GE0[M00GEnM11],2, 1],  umatsa[M22GE0[M00GEnM11],2, 0]-umatsa[M22GE0[M00GEnM11],0, 2],  
                                 umatsa[M22GE0[M00GEnM11],0, 1]-umatsa[M22GE0[M00GEnM11],1, 0]]
    except:
        pass

    Q[0,:] *= 0.5 / np.sqrt(T)
    Q[1,:] *= 0.5 / np.sqrt(T)
    Q[2,:] *= 0.5 / np.sqrt(T)
    Q[3,:] *= 0.5 / np.sqrt(T)
    return Q

"""
Highly optimized disorientation computation specifically for cubic symmetry.

This module provides specialized implementations that are significantly faster
than the general algorithm by exploiting the specific structure of cubic symmetry.
"""


# Hard-coded cubic symmetry operations (24 proper rotations)
# These are the 24 rotation matrices for point group O (432)
CUBIC_SYM_24 = np.array([
    # Identity
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    # 90° rotations around principal axes (z, y, x)
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
    [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
    [[0, 0, -1], [0, 1, 0], [1, 0, 0]],
    [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
    [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
    [[1, 0, 0], [0, 0, 1], [0, -1, 0]],
    # 120° rotations around <111> directions
    [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
    [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
    [[0, 0, -1], [-1, 0, 0], [0, 1, 0]],
    [[0, -1, 0], [0, 0, 1], [-1, 0, 0]],
    [[0, 0, -1], [1, 0, 0], [0, -1, 0]],
    [[0, 1, 0], [0, 0, -1], [-1, 0, 0]],
    [[0, 0, 1], [-1, 0, 0], [0, -1, 0]],
    [[0, -1, 0], [0, 0, -1], [1, 0, 0]],
    # 180° rotations around <110> directions
    [[0, 1, 0], [1, 0, 0], [0, 0, -1]],
    [[0, -1, 0], [-1, 0, 0], [0, 0, -1]],
    [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[-1, 0, 0], [0, 0, -1], [0, -1, 0]],
    [[0, 0, 1], [0, -1, 0], [1, 0, 0]],
    [[0, 0, -1], [0, -1, 0], [-1, 0, 0]],
], dtype=np.float64)

# Cubic symmetry with inversion (48 operations, Laue group m-3m)
# These are used when you want to include the inversion center
CUBIC_SYM_48 = np.zeros((48, 3, 3), dtype=np.float64)
CUBIC_SYM_48[:24] = CUBIC_SYM_24
CUBIC_SYM_48[24:] = -CUBIC_SYM_24  # Add inversions


@nb.njit(fastmath=True)
def _matrix_multiply_3x3(A, B):
    """Fast 3x3 matrix multiplication."""
    C = np.empty((3, 3), dtype=np.float64)
    for i in range(3):
        for j in range(3):
            C[i, j] = A[i, 0] * B[0, j] + A[i, 1] * B[1, j] + A[i, 2] * B[2, j]
    return C


@nb.njit(fastmath=True)
def _trace_3x3(A):
    """Fast trace computation for 3x3 matrix."""
    return A[0, 0] + A[1, 1] + A[2, 2]


@nb.njit(fastmath=True)
def _misorientation_angle_from_matrix(Delta):
    """
    Compute misorientation angle from rotation matrix.
    Uses: angle = arccos((trace(Delta) - 1) / 2)
    """
    trace_val = _trace_3x3(Delta)
    cos_angle = (trace_val - 1.0) * 0.5
    
    # Clamp to avoid numerical issues
    if cos_angle > 1.0:
        cos_angle = 1.0
    elif cos_angle < -1.0:
        cos_angle = -1.0
    
    return np.arccos(cos_angle)


@nb.njit(parallel=True, fastmath=True)
def compute_cubic_disorientations_24(M):
    """
    Compute disorientations for cubic symmetry (24 operations).
    
    This is optimized specifically for cubic crystal symmetry (point group O, 432).
    Uses only proper rotations (no inversion).
    
    Parameters:
    -----------
    M : ndarray, shape (m, 3, 3)
        Array of orientation matrices (proper rotations)
    
    Returns:
    --------
    D : ndarray, shape (m, m)
        Disorientation matrix in radians
    """
    m = M.shape[0]
    D = np.zeros((m, m), dtype=np.float64)
    
    # Get symmetry operations
    S = CUBIC_SYM_24
    n_sym = 24
    
    # Precompute M^T
    M_T = np.empty_like(M)
    for i in range(m):
        for k1 in range(3):
            for k2 in range(3):
                M_T[i, k1, k2] = M[i, k2, k1]
    
    # Precompute S^T
    S_T = np.empty_like(S)
    for i in range(n_sym):
        for k1 in range(3):
            for k2 in range(3):
                S_T[i, k1, k2] = S[i, k2, k1]
    
    # Precompute all symmetry-transformed orientations
    # M_S[i, si] = M[i] @ S[si]
    M_S = np.empty((m, n_sym, 3, 3), dtype=np.float64)
    for i in range(m):
        for si in range(n_sym):
            M_S[i, si] = _matrix_multiply_3x3(M[i], S[si])
    
    # Precompute all left-multiplied symmetry transformations
    # S_M_T[i, si] = S_T[si] @ M_T[i]
    S_M_T = np.empty((m, n_sym, 3, 3), dtype=np.float64)
    for i in range(m):
        for si in range(n_sym):
            S_M_T[i, si] = _matrix_multiply_3x3(S_T[si], M_T[i])
    
    # Compute pairwise disorientations
    for i in nb.prange(m):
        for j in range(i, m):
            if i == j:
                D[i, j] = 0.0
            else:
                min_angle = np.pi
                
                # Try all symmetry combinations
                for si in range(n_sym):
                    for sj in range(n_sym):
                        # Delta = S_T[si] @ M_T[i] @ M[j] @ S[sj]
                        Delta = _matrix_multiply_3x3(S_M_T[i, si], M_S[j, sj])
                        angle = _misorientation_angle_from_matrix(Delta)
                        
                        if angle < min_angle:
                            min_angle = angle
                            # Early termination for very small angles
                            if angle < 1e-6:
                                break
                    
                    if min_angle < 1e-6:
                        break
                
                D[i, j] = min_angle
                D[j, i] = min_angle
    
    return D


@nb.njit(parallel=True, fastmath=True)
def compute_cubic_disorientations_48(M):
    """
    Compute disorientations for cubic symmetry with inversion (48 operations).
    
    This is for Laue group m-3m (Oh), which includes inversion symmetry.
    Use this when your material has a center of inversion.
    
    Parameters:
    -----------
    M : ndarray, shape (m, 3, 3)
        Array of orientation matrices
    
    Returns:
    --------
    D : ndarray, shape (m, m)
        Disorientation matrix in radians
    """
    m = M.shape[0]
    D = np.zeros((m, m), dtype=np.float64)
    
    # Get symmetry operations
    S = CUBIC_SYM_48
    n_sym = 48
    
    # Precompute M^T
    M_T = np.empty_like(M)
    for i in range(m):
        for k1 in range(3):
            for k2 in range(3):
                M_T[i, k1, k2] = M[i, k2, k1]
    
    # Precompute S^T
    S_T = np.empty_like(S)
    for i in range(n_sym):
        for k1 in range(3):
            for k2 in range(3):
                S_T[i, k1, k2] = S[i, k2, k1]
    
    # Precompute all symmetry-transformed orientations
    M_S = np.empty((m, n_sym, 3, 3), dtype=np.float64)
    for i in range(m):
        for si in range(n_sym):
            M_S[i, si] = _matrix_multiply_3x3(M[i], S[si])
    
    # Precompute all left-multiplied symmetry transformations
    S_M_T = np.empty((m, n_sym, 3, 3), dtype=np.float64)
    for i in range(m):
        for si in range(n_sym):
            S_M_T[i, si] = _matrix_multiply_3x3(S_T[si], M_T[i])
    
    # Compute pairwise disorientations
    for i in nb.prange(m):
        for j in range(i, m):
            if i == j:
                D[i, j] = 0.0
            else:
                min_angle = np.pi
                
                # Try all symmetry combinations
                for si in range(n_sym):
                    for sj in range(n_sym):
                        Delta = _matrix_multiply_3x3(S_M_T[i, si], M_S[j, sj])
                        angle = _misorientation_angle_from_matrix(Delta)
                        
                        if angle < min_angle:
                            min_angle = angle
                            if angle < 1e-6:
                                break
                    
                    if min_angle < 1e-6:
                        break
                
                D[i, j] = min_angle
                D[j, i] = min_angle
    
    return D


@nb.njit(parallel=True, fastmath=True)
def compute_cubic_disorientations_24_ultra_fast(M):
    """
    Ultra-fast version using maximum vectorization and manual loop unrolling.
    
    This version is optimized for maximum speed on modern CPUs.
    Best for large datasets (m > 500).
    
    Parameters:
    -----------
    M : ndarray, shape (m, 3, 3)
        Array of orientation matrices
    
    Returns:
    --------
    D : ndarray, shape (m, m)
        Disorientation matrix in radians
    """
    m = M.shape[0]
    D = np.zeros((m, m), dtype=np.float64)
    
    S = CUBIC_SYM_24
    n_sym = 24
    
    # Precompute everything upfront
    M_T = np.empty((m, 3, 3), dtype=np.float64)
    for i in range(m):
        for k1 in range(3):
            for k2 in range(3):
                M_T[i, k1, k2] = M[i, k2, k1]
    
    S_T = np.empty((n_sym, 3, 3), dtype=np.float64)
    for i in range(n_sym):
        for k1 in range(3):
            for k2 in range(3):
                S_T[i, k1, k2] = S[i, k2, k1]
    
    # Precompute M @ S for all combinations
    M_S = np.empty((m, n_sym, 3, 3), dtype=np.float64)
    for i in range(m):
        for si in range(n_sym):
            for k1 in range(3):
                for k2 in range(3):
                    M_S[i, si, k1, k2] = (M[i, k1, 0] * S[si, 0, k2] + 
                                           M[i, k1, 1] * S[si, 1, k2] + 
                                           M[i, k1, 2] * S[si, 2, k2])
    
    # Precompute S^T @ M^T for all combinations
    S_M_T = np.empty((m, n_sym, 3, 3), dtype=np.float64)
    for i in range(m):
        for si in range(n_sym):
            for k1 in range(3):
                for k2 in range(3):
                    S_M_T[i, si, k1, k2] = (S_T[si, k1, 0] * M_T[i, 0, k2] + 
                                             S_T[si, k1, 1] * M_T[i, 1, k2] + 
                                             S_T[si, k1, 2] * M_T[i, 2, k2])
    
    # Main computation loop
    for i in nb.prange(m):
        for j in range(i + 1, m):
            min_angle = np.pi
            
            for si in range(n_sym):
                # Load S_M_T[i, si] into local variables for better cache performance
                a00 = S_M_T[i, si, 0, 0]
                a01 = S_M_T[i, si, 0, 1]
                a02 = S_M_T[i, si, 0, 2]
                a10 = S_M_T[i, si, 1, 0]
                a11 = S_M_T[i, si, 1, 1]
                a12 = S_M_T[i, si, 1, 2]
                a20 = S_M_T[i, si, 2, 0]
                a21 = S_M_T[i, si, 2, 1]
                a22 = S_M_T[i, si, 2, 2]
                
                for sj in range(n_sym):
                    # Manual matrix multiplication and trace computation
                    # trace = (A @ B)[0,0] + (A @ B)[1,1] + (A @ B)[2,2]
                    trace = (a00 * M_S[j, sj, 0, 0] + a01 * M_S[j, sj, 1, 0] + a02 * M_S[j, sj, 2, 0] +
                             a10 * M_S[j, sj, 0, 1] + a11 * M_S[j, sj, 1, 1] + a12 * M_S[j, sj, 2, 1] +
                             a20 * M_S[j, sj, 0, 2] + a21 * M_S[j, sj, 1, 2] + a22 * M_S[j, sj, 2, 2])
                    
                    cos_angle = (trace - 1.0) * 0.5
                    
                    # Clamp
                    if cos_angle > 1.0:
                        cos_angle = 1.0
                    elif cos_angle < -1.0:
                        cos_angle = -1.0
                    
                    angle = np.arccos(cos_angle)
                    
                    if angle < min_angle:
                        min_angle = angle
                        if angle < 1e-6:
                            break
                
                if min_angle < 1e-6:
                    break
            
            D[i, j] = min_angle
            D[j, i] = min_angle
    
    return D
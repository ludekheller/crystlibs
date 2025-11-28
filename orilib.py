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
    g = np.eye(3)
    s1, s2, s3 = np.sin(ai), np.sin(aj), np.sin(ak)
    c1, c2, c3 = np.cos(ai), np.cos(aj), np.cos(ak)
    
    g[0,0] = c1*c3 - s1*s3*c2
    g[0,1] = s1*c3 + c1*s3*c2
    g[0,2] = s3*s2    
    g[1,0] = -c1*s3 - s1*c3*c2 
    g[1,1] = -s1*s3 + c1*c3*c2
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
    U = np.eye(3)

    s1, s2, s3 = np.sin(ai), np.sin(aj), np.sin(ak)
    c1, c2, c3 = np.cos(ai), np.cos(aj), np.cos(ak)

    U[0,0] = c1*c3 - s1*s3*c2
    U[0,1] = -c1*s3 - s1*c3*c2
    U[0,2] = s1*s2    
    U[1,0] = s1*c3 + c1*s3*c2 
    U[1,1] = -s1*s3 + c1*c3*c2
    U[1,2] = -c1*s2 
    U[2,0] = s3*s2 
    U[2,1] = c3*s2
    U[2,2] = c2   
    return U


# =====================================
# Rodrigues-Frank and Axis-Angle Representations
# =====================================

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
    eps = 1.e-6
    
    ptheta = np.arccos((g[0][0] + g[1][1] + g[2][2] - 1) / 2)
    r = [0., 0., 0.]
    
    if ptheta < eps:
        # Small angle - axis is arbitrary
        r[0] = 1
        r[1] = 0
        r[2] = 0
    elif ptheta < (1 - eps)*np.pi:
        # Normal case
        r[0] = (g[1][2] - g[2][1]) / (2 * np.sin(ptheta))
        r[1] = (g[2][0] - g[0][2]) / (2 * np.sin(ptheta))
        r[2] = (g[0][1] - g[1][0]) / (2 * np.sin(ptheta))
    else:
        # Near 180 degrees
        r[0] = np.sqrt((g[0][0] + 1) / 2)
        r[1] = np.sqrt((g[1][1] + 1) / 2)
        r[2] = np.sqrt((g[2][2] + 1) / 2)
    
    m = r.index(max(r))
    for i in range(0, 3):
        if not r == m:
            if g[i][m] < 0:
                r[i] *= 1
    
    return r, ptheta


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
    eps = 1.e-6
    
    ptheta = np.arccos((np.trace(g) - 1) / 2)
    r = np.array([0., 0., 0.])
    
    if ptheta < eps:
        r[0] = 1
        r[1] = 0
        r[2] = 0
    elif ptheta < (1 - eps)*np.pi:
        r[0] = (g[1,2] - g[2,1]) / (2 * np.sin(ptheta))
        r[1] = (g[2,0] - g[0,2]) / (2 * np.sin(ptheta))
        r[2] = (g[0,1] - g[1,0]) / (2 * np.sin(ptheta))
    else:
        r[0] = np.sqrt((g[0,0] + 1) / 2)
        r[1] = np.sqrt((g[1,1] + 1) / 2)
        r[2] = np.sqrt((g[2,2] + 1) / 2)
    
    m = np.where(r == max(r))[0][0]
    for i in range(0, 3):
        if not i == m:
            if g[i,m] < 0:
                r[i] *= 1
    
    return r, ptheta


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

    g[0][0] = r[0] * r[0] * (1 - np.cos(theta)) + np.cos(theta)
    g[0][1] = r[0] * r[1] * (1 - np.cos(theta)) + r[2] * np.sin(theta)
    g[0][2] = r[0] * r[2] * (1 - np.cos(theta)) - r[1] * np.sin(theta)
    
    g[1][0] = r[1] * r[0] * (1 - np.cos(theta)) - r[2] * np.sin(theta)
    g[1][1] = r[1] * r[1] * (1 - np.cos(theta)) + np.cos(theta)
    g[1][2] = r[1] * r[2] * (1 - np.cos(theta)) + r[0] * np.sin(theta)
    
    g[2][0] = r[2] * r[0] * (1 - np.cos(theta)) + r[1] * np.sin(theta)
    g[2][1] = r[2] * r[1] * (1 - np.cos(theta)) - r[0] * np.sin(theta)
    g[2][2] = r[2] * r[2] * (1 - np.cos(theta)) + np.cos(theta)
    
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
    g = np.zeros((3,3))

    g[0,0] = r[0] * r[0] * (1 - np.cos(theta)) + np.cos(theta)
    g[0,1] = r[0] * r[1] * (1 - np.cos(theta)) + r[2] * np.sin(theta)
    g[0,2] = r[0] * r[2] * (1 - np.cos(theta)) - r[1] * np.sin(theta)
    
    g[1,0] = r[1] * r[0] * (1 - np.cos(theta)) - r[2] * np.sin(theta)
    g[1,1] = r[1] * r[1] * (1 - np.cos(theta)) + np.cos(theta)
    g[1,2] = r[1] * r[2] * (1 - np.cos(theta)) + r[0] * np.sin(theta)
    
    g[2,0] = r[2] * r[0] * (1 - np.cos(theta)) + r[1] * np.sin(theta)
    g[2,1] = r[2] * r[1] * (1 - np.cos(theta)) - r[0] * np.sin(theta)
    g[2,2] = r[2] * r[2] * (1 - np.cos(theta)) + np.cos(theta)
    
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

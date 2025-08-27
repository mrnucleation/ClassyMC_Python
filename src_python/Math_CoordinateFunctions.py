import numpy as np
from math import pi, sqrt, sin, cos, acos
import random

# Assuming dp is like float64, but in Python we use float
# grnd is random.uniform(0,1)
# norm2 is np.linalg.norm**2 but actually norm2(v) is sqrt(sum v**2) no, in Fortran norm2 is sqrt(dot), so magnitude
# error stop -> raise ValueError

def GetPerpendicularVector(v1, phi=None):
    # v1, v2 are lists or np arrays of 3 elements
    # assume np.array for ease
    v1 = np.array(v1)
    v2 = np.zeros(3)
    two_pi = 2 * pi
    if phi is not None:
        tors_angle = phi
    else:
        tors_angle = two_pi * random.random()
    
    s_term = sin(tors_angle)
    c_term = cos(tors_angle)      
    r_proj = sqrt(v1[0]**2 + v1[1]**2)
    
    # coeff2 = 1.0 / r_proj
    coeff2 = 1.0 / r_proj
    coeff3 = coeff2 / np.linalg.norm(v1)
    
    v2[0] = -coeff2 * c_term * v1[1] - coeff3 * s_term * v1[0] * v1[2]
    v2[1] = +coeff2 * c_term * v1[0] - coeff3 * s_term * v1[1] * v1[2]
    v2[2] = + coeff3 * s_term * (r_proj * r_proj)
    
    return v2
import numpy as np
from math import pi, sin, cos, acos
import random

# Assuming 'dp' corresponds to double precision, which in Python is float
# 'grnd()' is replaced with random.random() for uniform [0,1)
# 'norm2(v)' is replaced with np.linalg.norm(v)
# 'error stop' is replaced with raise ValueError
# Vectors are assumed to be numpy arrays of shape (3,) for simplicity
# External modules like VarPrecision, ClassyConstants, RandomGen, CoordinateTypes are not present;
#   we inline constants like pi, two_pi and use built-ins/random

def GetPerpendicularVector(v1, phi=None):
    """
    Produces a vector, v2, that is perpendicular to the input v1 vector.
    
    Parameters:
    v1 : numpy.ndarray of shape (3,)
        The Input Reference Vector which the user needs a perpendicular vector towards
    phi : float, optional
        Rotational angle. If not provided, a random angle is used.
    
    Returns:
    numpy.ndarray of shape (3,)
        The Output Vector that will be perpendicular to v1
    """
    v1 = np.array(v1)
    two_pi = 2 * pi
    if phi is not None:
        tors_angle = phi
    else:
        tors_angle = two_pi * random.random()
    
    s_term = sin(tors_angle)
    c_term = cos(tors_angle)      
    r_proj = np.sqrt(v1[0]**2 + v1[1]**2)
    
    # Since the bond_ang is fixed at pi/2 the first coefficient drops to 0 and the
    # second's sine term is fixed at 1.
    coeff2 = 1.0 / r_proj
    r1 = np.linalg.norm(v1)
    coeff3 = coeff2 / r1
    
    v2 = np.zeros(3)
    v2[0] = -coeff2 * c_term * v1[1] - coeff3 * s_term * v1[0] * v1[2]
    v2[1] = +coeff2 * c_term * v1[0] - coeff3 * s_term * v1[1] * v1[2]
    v2[2] = + coeff3 * s_term * (r_proj * r_proj)
    
    return v2

def Generate_RelativeVector(v1, r2, bond_ang, phi):
    """
    The purpose of this function is to generate a vector for an atom
    given a fixed bond dist, angle, and torsion with a given vector (v1).
    The coordinate is created using a relative orthonormal framework given by these vectors
    w1=(x1,y1,z1)   w2=(-y1,x1,0)  w3=(-x1*z1, -y1*z1, x1^2 + y1^2)
    Where w1 is a normalized v1, w2 is perpendicular to v1, and w3 is generated
    by the cross product w3 = (w1 x w2).
    Using these vectors, the new vector(v2) is calculated using a rotational matrix
    
    Parameters:
    v1 : numpy.ndarray of shape (3,)
        Input vector
    r2 : float
        Bond distance
    bond_ang : float
        Bond angle
    phi : float
        Torsion angle
    
    Returns:
    numpy.ndarray of shape (3,)
        Output vector v2
    """
    v1 = np.array(v1)
    r1 = np.linalg.norm(v1)
    s_term = sin(phi)
    c_term = cos(phi)      
    r_proj = np.sqrt(v1[0]**2 + v1[1]**2)
    
    coeff1 = (r2 / r1) * cos(bond_ang)
    coeff2 = (r2 / r_proj) * sin(bond_ang)
    coeff3 = coeff2 / r1

    v2 = np.zeros(3)
    v2[0] = coeff1 * v1[0] - coeff2 * c_term * v1[1] - coeff3 * s_term * v1[0] * v1[2]
    v2[1] = coeff1 * v1[1] + coeff2 * c_term * v1[0] - coeff3 * s_term * v1[1] * v1[2]
    v2[2] = coeff1 * v1[2] + coeff3 * s_term * (r_proj * r_proj)
    
    return v2

def AngleFromVectors(v1, v2):
    """
    Standard dot product angle calculation for computing the angle between vectors.
    
    Parameters:
    v1, v2 : numpy.ndarray
        Input vectors (any shape, but compatible for dot product)
    
    Returns:
    float
        The angle theta in radians
    """
    r1 = np.linalg.norm(v1)
    if r1 <= 0.0:
        raise ValueError("Vector of zero length has been passed to AngleFromVectors function!")
    r2 = np.linalg.norm(v2)
    if r2 <= 0.0:
        raise ValueError("Vector of zero length has been passed to AngleFromVectors function!")
    dotprod = np.dot(v1, v2)
    theta = dotprod / (r1 * r2)
    theta = acos(theta)
    return theta

def DihedralAngle(ang1, ang2, ang3):
    """
    Computes the dihedral angle made by three different bond angles.
    
    Parameters:
    ang1, ang2, ang3 : float
        Input angles
    
    Returns:
    float
        The dihedral angle
    """
    dihed = acos((cos(ang3) - cos(ang1) * cos(ang2)) / (sin(ang1) * sin(ang2)))
    return dihed

def DihedralAngle_FromVecs(v1, v2, v3):
    """
    Computes the dihedral angle made by three different vectors
    
    Parameters:
    v1, v2, v3 : numpy.ndarray of shape (3,)
        Input vectors
    
    Returns:
    float
        The dihedral angle
    """
    theta12 = AngleFromVectors(v1, v2)
    theta13 = AngleFromVectors(v1, v3)
    theta23 = AngleFromVectors(v2, v3)
    dihed = DihedralAngle(theta12, theta13, theta23)
    return dihed
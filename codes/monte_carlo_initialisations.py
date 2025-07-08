import numpy as np


def transform_direction(direction, axis):
    axis = axis / np.linalg.norm(axis)
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(z, axis)
    s = np.linalg.norm(v)
    c = np.dot(z, axis)

    if s < 1e-8:
        if c > 0:
            return direction
        else:
            return -direction

    vx = np.array([[0, -v[2], v[1]],
                   [v[2], 0, -v[0]],
                   [-v[1], v[0], 0]])

    R = np.eye(3) + vx + vx @ vx * ((1 - c) / s**2)
    return R @ direction


def isotropic_direction_in_angle(angle):
    
    if angle < 0 or angle > np.pi:
        raise ValueError("Angle must be between 0 and pi radians.")
    
    n = np.zeros(3)
    n[2] = np.random.rand()*(1 - np.cos(angle)) + np.cos(angle)

    theta = np.arccos(n[2])
    beta = 2 * np.pi * np.random.rand()

    n[0] = np.sin(theta) * np.cos(beta)
    n[1] = np.sin(theta) * np.sin(beta)


    return n


def isotropic_direction_in_cone(angle, source_position, detector_height, detector_radius):
    axis = -source_position / np.linalg.norm(source_position)

    if np.sqrt(detector_radius**2 + detector_height**2 / 4) >= np.linalg.norm(source_position):
        n = isotropic_direction_in_angle(np.pi)

    else:
        n = isotropic_direction_in_angle(angle)
        n = transform_direction(n, axis)
        
    return n


def photon_angle(energy_in):
    alpha = energy_in / 0.511  
    
    while True:

        epsilon = 1 + 2*alpha  # Maximum energy ratio (for backscatter)
        r1, r2 = np.random.rand(2)
        
        if r1 < (alpha + 1)/(9*alpha + 1):
            # Sample from high-acceptance region
            epsilon = 1 + 2*alpha*r2
            g = (1 + 2*alpha/epsilon + (1/epsilon)**2)/4
        else:
            # Sample from low-rejection region
            epsilon = (1 + 2*alpha)/(1 + 2*alpha*r2)
            g = (epsilon + 1/epsilon - (1 - ((epsilon - 1)/(alpha*epsilon))**2)/2)
        
        # Acceptance test
        if np.random.rand() <= g:
            break
    
    costheta = 1 - (epsilon - 1)/alpha
    theta = np.arccos(np.clip(costheta, -1, 1))
    energy_out = energy_in/epsilon
    
    return theta, energy_out


def photon_direction(angle):

    n = np.zeros(3)
    n[2] = np.cos(angle)

    rho = np.sin(angle)
    phi = np.random.rand() * 2 * np.pi

    n[0] = rho * np.cos(phi)
    n[1] = rho * np.sin(phi)

    return n


def compton_scatter_photon(energy_in, direction):
    
    angle, energy_out = photon_angle(energy_in)
    direction_around_z = photon_direction(angle)

    direction = transform_direction(direction_around_z, direction)

    detector_energy = energy_in - energy_out

    return direction, energy_out, detector_energy
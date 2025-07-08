import numpy as np


def calc_angle_of_cone(source_position, detector_height, detector_radius):
    R_g = np.sqrt(detector_radius**2 + detector_height**2 / 4)

    if np.linalg.norm(source_position) <= R_g:
        alpha = np.pi
    else:
        alpha = np.arcsin(R_g / np.linalg.norm(source_position))

    return alpha


def intersect_plane(position, direction, planeZ):
    
    if np.sqrt(direction[0]**2 + direction[1]**2 + direction[2]**2) < 1e-3:
        dist = np.inf
        return dist
    
    else:
        dist = (planeZ - position[2]) / direction[2]

        return dist
    

def intersect_cylinder(position, direction, radius):

    dist = np.zeros(2)

    a = direction[0]**2 + direction[1]**2
    b = 2 * (position[0] * direction[0] + position[1] * direction[1])
    c = position[0]**2 + position[1]**2 - radius**2
    d = b**2 - 4 * a * c

    if d < 0:

        dist[0] = np.inf
        dist[1] = np.inf
        return dist
    
    else:
        dist[0] = (-b + np.sqrt(d)) / (2 * a)
        dist[1] = (-b - np.sqrt(d)) / (2 * a)

        return dist



def intersect_cylinder_in(position, direction, pztop, pzbottom, radius):
    d1, d2 = intersect_cylinder(position, direction, radius)
    d_plus = intersect_plane(position, direction, pztop)
    d_minus = intersect_plane(position, direction, pzbottom)

    dist = np.min([np.max([d1, d2]), np.max([d_plus, d_minus])])

    return dist           


def intersect_cylinder_out(position, direction, pztop, pzbottom, radius):
    d1, d2 = intersect_cylinder(position, direction, radius)
    d_min = np.min([d1, d2])
    d_max = np.max([d1, d2])

    if d_min == np.inf or d_max < 0 or ((position[2] > pztop and direction[2] > 0) or (position[2] < pzbottom and direction[2] < 0)):
        dist = np.inf
    else:
        if d1*d2 > 0:
            d_cyl = d_min

        else:
            d_cyl = d_max

        z_cyl = position[2] + d_cyl * direction[2]

        if z_cyl > pztop or z_cyl < pzbottom:
            d_cyl = np.inf
        
        dtop = intersect_plane(position, direction, pztop)
        xtop = position[0] + dtop * direction[0]
        ytop = position[1] + dtop * direction[1]

        if xtop**2 + ytop**2 > radius**2:
            dtop = np.inf

        dbottom = intersect_plane(position, direction, pzbottom)
        xbottom = position[0] + dbottom * direction[0]
        ybottom = position[1] + dbottom * direction[1]


        if xbottom**2 + ybottom**2 > radius**2:
            dbottom = np.inf


        dist = np.min([d_cyl, dtop, dbottom])

    return dist


def intersect_cylinder_starting_points(position, direction, pztop, pzbottom, radius):
    dist = intersect_cylinder_out(position, direction, pztop, pzbottom, radius)
    if dist != np.inf:
        point = position + dist * direction
        return point
    else:
        return None
    

def goes_outside(position, direction, pztop, pzbottom, radius, lambd):
    dist = intersect_cylinder_in(position, direction, pztop, pzbottom, radius)
    out = dist < lambd
    return out
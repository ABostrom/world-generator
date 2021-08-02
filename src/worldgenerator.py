



from random import random, uniform
from map import Map, voronoi_polygons_2D, in_box
from scipy.spatial import Voronoi, Delaunay, voronoi_plot_2d, delaunay_plot_2d

import numpy as np

def create_points(width, height, spacing):
    radius = spacing / 2
    return np.array([j_point(x, y, radius) for y in range(-1, height+1, spacing) for x in range(-1, width+1, spacing)])

def j_point(x, y, radius):
    return [j_val(x, radius), j_val(y, radius)]

def j_val(x, radius):
    return x + (uniform(-radius, radius) / 2)



def centroid_region(vertices):
    # Polygon's signed area
    A = 0
    # Centroid's x
    C_x = 0
    # Centroid's y
    C_y = 0
    for i in range(0, len(vertices) - 1):
        s = (vertices[i, 0] * vertices[i + 1, 1] - vertices[i + 1, 0] * vertices[i, 1])
        A = A + s
        C_x = C_x + (vertices[i, 0] + vertices[i + 1, 0]) * s
        C_y = C_y + (vertices[i, 1] + vertices[i + 1, 1]) * s
    A = 0.5 * A
    C_x = (1.0 / (6.0 * A)) * C_x
    C_y = (1.0 / (6.0 * A)) * C_y
    return np.array([[C_x, C_y]])

# Generates a bounded vornoi diagram with finite regions
def bounded_voronoi(points, bounding_box):
    # # Select towers inside the bounding box
    i = in_box(points, bounding_box)

    # Mirror points left, right, above, and under to provide finite regions for the edge regions of the bounding box
    points_center = points[i, :]

    points_left = np.copy(points_center)
    points_left[:, 0] = bounding_box[0] - (points_left[:, 0] - bounding_box[0])

    points_right = np.copy(points_center)
    points_right[:, 0] = bounding_box[1] + (bounding_box[1] - points_right[:, 0])

    points_down = np.copy(points_center)
    points_down[:, 1] = bounding_box[2] - (points_down[:, 1] - bounding_box[2])

    points_up = np.copy(points_center)
    points_up[:, 1] = bounding_box[3] + (bounding_box[3] - points_up[:, 1])

    points = np.append(points_center,
                       np.append(np.append(points_left,
                                           points_right,
                                           axis=0),
                                 np.append(points_down,
                                           points_up,
                                           axis=0),
                                 axis=0),
                       axis=0)
    # Compute Voronoi
    vor = Voronoi(points)

    vor.filtered_points = points_center # creates a new attibute for points that form the diagram within the region
    vor.filtered_regions = np.array(vor.regions)[vor.point_region[:vor.npoints//5]] # grabs the first fifth of the regions, which are the original regions
    vor.point_region = vor.point_region[:vor.npoints//5]

    return vor

# Performs x iterations of loyd's algorithm to calculate a centroidal vornoi diagram
def loyds_algorithm(points, iterations, bounding_box):
    for i in range(iterations):
        vor = bounded_voronoi(points, bounding_box)
        centroids = []

        for region in vor.filtered_regions:
        #for region in vor.regions:
            vertices = vor.vertices[region + [region[0]], :] # grabs vertices for the region and adds a duplicate of the first one to the end
            centroid = centroid_region(vertices)
            centroids.append(list(centroid[0, :]))

        points = np.array(centroids)

    vor = bounded_voronoi(points, bounding_box)
    vor.point_region[:vor.npoints//5]
    return vor



width = 10
height = 10
bounding_box = np.array([0., width, 0., height])
points = create_points(width, height, 1)

voronoi = loyds_algorithm(points, 10, bounding_box)

map = Map(bounding_box, voronoi)

map.draw_land()
input()

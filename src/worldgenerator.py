



from random import random, uniform
from map import Map, voronoi_polygons_2D
from scipy.spatial import Voronoi, Delaunay, voronoi_plot_2d, delaunay_plot_2d

import numpy as np

def create_points(width, height, spacing):
    radius = int(spacing / 2)
    return np.array([j_point(x, y, radius) for y in range(0, height, spacing) for x in range(0, width, spacing)])

def j_point(x, y, radius):
    return [j_val(x, radius), j_val(y, radius)]

def j_val(x, radius):
    return x + (uniform(-radius, radius) / 2)


def loyds_algorithm(points, iters):

    if iters == 1:
        return Voronoi(points, qhull_options=f'Qbb Qc')

    for i in range(iters):
        print(len(points))
        voronoi = Voronoi(points, qhull_options=f'Qbb Qc')
        vertices = np.asarray(voronoi.vertices)
        def strip(reg):
            reg = np.asarray(reg)
            return reg[np.where(reg >= 0)]

        points = np.asarray([vertices[strip(region)].mean(axis=0) for region in voronoi.regions if len(region) > 0])

    return voronoi


def inside_bounds(pt):
    print(pt)


width = 50
height = 50
points = create_points(width, height, 4)


print(points)
#voronoi = Voronoi(points, qhull_options='Qbb Qc Qz')
voronoi = loyds_algorithm(points, 1)
delaunay = Delaunay(points)

map = Map(width, height, points, voronoi, delaunay)
#map.draw_vor_centres()

map.draw_land()
input()





from random import random, uniform
from map import Map
from scipy.spatial import Voronoi, Delaunay, voronoi_plot_2d, delaunay_plot_2d

import numpy as np

def create_points(width, height, spacing):
    radius = int(spacing / 2)
    return np.array([j_point(x, y, radius) for y in range(0, height, spacing) for x in range(0, width, spacing)])

def j_point(x, y, radius):
    return [j_val(x, radius), j_val(y, radius)]

def j_val(x, radius):
    return x + (uniform(-radius, radius) / 2)



width = 100
height = 100
points = create_points(width, height, 4)
delaunay = Delaunay(points)

print(points)
voronoi = Voronoi(points, qhull_options='Qbb Qc Qz')
map = Map(width, height, points, voronoi, delaunay)
#map.draw_vor_centres()

map.draw_land()
input()

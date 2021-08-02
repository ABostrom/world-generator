import enum
from matplotlib.collections import LineCollection, PolyCollection
from scipy.spatial import Voronoi, Delaunay, voronoi_plot_2d, delaunay_plot_2d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
from collections import defaultdict

def voronoi_polygons_2D(vor, return_ridges=False):
    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    radius = vor.points.ptp().max()
    #radius = max(self.width, self.height)
    # Construct a map containing all ridges for a given point
    all_ridges = defaultdict(list)
    #construct map of all regions and all ridge points.
    #foreach pair of points that form a ridge. extract the ridge points and vertices.
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges[p1].append((p2, v1, v2))
        all_ridges[p2].append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

        new_regions, np.asarray(new_vertices)

    if return_ridges:
        return new_regions, np.asarray(new_vertices), all_ridges
    else:
        return new_regions, np.asarray(new_vertices)
    

class Cell:

    def __init__(self, pt,vert, region, region_verts, neighbours) -> None:
        self.index = pt
        self.centre = np.asarray(vert)
        self.region = region        #this is an 
        self.corners = region_verts #this is an index into the Corners array.
        self.neighbours = neighbours

        #setup the borders and the corners

    def __repr__(self):
        return f'index:{self.index}'

#TODO: This class contains the Edge data
# Should have the Edge index
# The Two indexes for the verts
# The vert data
# The two polygons that it borders
class Edge:
    
    def __init__(self, v0, v1, d0, d1) -> None:
        
        self.v0 = np.asarray(v0)
        self.v1 = np.asarray(v1)
        self.d0 = d0
        self.d1 = d1
        

#TODO: This class contains the corner data.
# Should have the corners index.
# Should have the edges it connects too
# Should have the polygons that it meets
class Corner:
    
    def __init__(self, vert, touches, protrudes, adjecent) -> None:
        self.vert = vert
        self.touches = touches # polygons this corner touches
        self.protrudes = protrudes #edges touching the corner
        self.adjacent = adjecent #corners adjacent to this one.


class Map:

    def __init__(self, width, height, points, voronoi, delaunay) -> None:
        
        self.points = voronoi.points
        self.vor = voronoi
        self.tri = delaunay
        self.width = width
        self.height = height
        self.regions = voronoi.regions
        self.vertices = voronoi.vertices

        connections = defaultdict(list)
        #calculate polygon neighbours
        for pt1, pt2 in voronoi.ridge_points:
            connections[pt1].append(pt2)
            connections[pt2].append(pt1)


        # create the cell graph.

        #TODO: self.vertices should be corners.
        self.cells = {pt: Cell(pt, self.points[pt], voronoi.regions[index], self.vertices[voronoi.regions[index]], connections[pt]) for pt, index in enumerate(voronoi.point_region) if -1 not in voronoi.regions[index]}
        self.num_cells = len(self.cells)
        self.elevation = np.random.normal(size=self.num_cells)

        print(self.cells)
        print(voronoi.ridge_vertices)


        #TODO: given a vert (index), get the cells
        #TODO: given a vert (index), get the edges
        vert_polys = defaultdict(set)
        vert_edges = defaultdict(list)
        vert_connections = defaultdict(set)

        # build set of polys for each vertex
        for (p1, p2), (v0, v1) in zip(voronoi.ridge_points, voronoi.ridge_vertices):
            #skip -1 indices  
            if -1  in (v0, v1):
                continue

            #given a vert index get the set of polys we're connected too.
            vert_polys[v0].update((p1, p2))
            vert_polys[v1].update((p1, p2))
            #poly pairs that we can use to index into the edges dictionary.
            vert_edges[v0].append((p1, p2))
            vert_edges[v1].append((p1, p2))
            #verts connected to this one.
            vert_connections[v0].add(v1)
            vert_connections[v1].add(v0)

        self.corners = [Corner(vert, vert_polys[i], vert_edges[i], vert_connections[i]) for i, vert in enumerate(voronoi.vertices)]

        #Data structure: given two poly indexes, get the edge
        #this will contain each edge, and the cells it bisects, the cells that touch one end of the corner, and the other.
        #given a pair of polygons, get the corresponding edge.
        #IE two adjacent poly's should contain a shared pair of vertices
        self.edges = {(p1, p2): Edge(self.vertices[v0], self.vertices[v1], self.cells.get(p1, None), self.cells.get(p2, None)) 
                        for (p1, p2), (v0, v1) in zip(voronoi.ridge_points, voronoi.ridge_vertices)
                        if -1  not in (v0, v1)}

        

    def compute_line_segments(self):
        vor = self.vor
        pts = self.pts

        #calculate the centre point of all points.
        center = vor.points.mean(axis=0)

        #calculate peak to peak for the points the get a boundary
        ptp_bound = vor.points.ptp(axis=0)

        #find the line segments from the ridge vertexes
        finite_segments = []
        for pointidx, simplex in zip(vor.ridge_points,vor.ridge_vertices):
            simplex = np.asarray(simplex)
            if np.all(simplex >= 0):
                finite_segments.append(vor.vertices[simplex])
            else:
                far_point, i= self.calculate_point(vor.points, vor.vertices, pointidx, simplex, center, ptp_bound, vor.furthest_site)
                finite_segments.append([vor.vertices[i], far_point])

        return finite_segments

    def generate_colours(self, data, color_pallette):
        # find min/max values for normalization
        minima = min(data)
        maxima = max(data)

        # normalize chosen colormap
        norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
        mapper = cm.ScalarMappable(norm=norm, cmap=color_pallette)

        return mapper



    def draw(self, data, color, neighbours=False, show_text=False, corner_connection=20):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        # for i, pt in enumerate(self.points):
        #     if self.cells.__contains__(i):
        #         ax.text(*(pt + (np.ones(2) / 4)), str(i), color = 'red')

        #plot the cells
        for i, (pt, cell) in enumerate(self.cells.items()):
            if show_text:
                ax.text(*cell.centre, str(cell.index), color='white')

            ax.fill(*zip(*cell.corners), color=color.to_rgba(data[i]))

            if neighbours:
                for neighbour in cell.neighbours:
                    if self.cells.__contains__(neighbour):
                        #TODO:
                        #find the line segment for the boundary. 
                        #find the midpoint.
                        #get the normal. 
                        #come out from the normal by 0.25
                        dir = (self.cells[neighbour].centre - cell.centre) * 0.25
                        if show_text:
                            ax.text(*(cell.centre + dir), str(self.cells[neighbour].index))

        for (p1, p2), edge in self.edges.items():
            mid = (edge.v0 + edge.v1) / 2
            ax.plot([edge.v0[0], edge.v1[0]], [edge.v0[1], edge.v1[1]], marker = 'o', color='black')
            #TODO: could normalise direction vector to make more consistent
            if show_text:
                ax.text(*(mid + ((self.cells[p1].centre - mid)*0.25)), str(self.cells[p1].index)) 
                ax.text(*(mid + ((self.cells[p2].centre - mid)*0.25)), str(self.cells[p2].index)) 


        #draw a specific corner, and it's associated connectivity.

        
        if corner_connection >= 0:
            #get all the corners from the cells.
            for cor in self.cells[corner_connection].region:
                corner = self.corners[cor]
                
                print(corner.vert)
                ax.plot(*corner.vert, marker = 'x', color='green')

                #draw connected edges.
                for e in corner.protrudes:
                    if self.edges.__contains__(e):
                        edge = self.edges[e]
                        ax.plot([edge.v0[0], edge.v1[0]], [edge.v0[1], edge.v1[1]], marker = 'o', color='red')

                #draw adjacent corners.
                for c in corner.adjacent:
                    ax.plot(*self.corners[c].vert, marker = 'x', color='green')

                
                for p in corner.touches:
                    if self.cells.__contains__(p):
                        ax.plot(*self.cells[p].centre, *corner.vert, marker='v', color='cyan')

        

        plt.xlim(self.vor.min_bound[0] - 0.1, self.vor.max_bound[0] - 0.1)
        plt.ylim(self.vor.min_bound[1] - 0.1, self.vor.max_bound[1] - 0.1)

        fig.show()


    def draw_land(self, **kwargs):
        self.draw(self.elevation, self.generate_colours(self.elevation, cm.terrain), kwargs)

    def calculate_point(self, points, verts, pointidx, simplex, center, ptp_bound, furthest_site):
        i = simplex[simplex >= 0][0]  # finite end Voronoi vertex
        t = points[pointidx[1]] - points[pointidx[0]]  # tangent
        t /= np.linalg.norm(t)
        # -dy / dx
        n = np.array([-t[1], t[0]])  # normal

        midpoint = points[pointidx].mean(axis=0)
        direction = np.sign(np.dot(midpoint - center, n)) * n
        if (furthest_site):
            direction = -direction

        far_point = verts[i] + direction * ptp_bound.max()

        return far_point, i
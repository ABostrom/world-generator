from matplotlib.collections import LineCollection, PolyCollection
from scipy.spatial import Voronoi, Delaunay, voronoi_plot_2d, delaunay_plot_2d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
from collections import defaultdict

class Cell:

    def __init__(self, pt,vert, region, region_verts) -> None:
        self.pt = pt
        self.vert = vert
        self.region = region
        self.region_verts = region_verts


class Map:

    def __init__(self, width, height, points, voronoi, delaunay) -> None:
        self.pts = points
        self.num_points = len(points)
        self.vor = voronoi
        self.tri = delaunay
        self.width = width
        self.height = height


        self.line_segments = self.compute_line_segments()
        regions, vertices = self.finite_polygons_2D(voronoi)

        self.cells = [Cell(pt, points[pt], region, vertices[region]) for pt, region in enumerate(regions)]


    def finite_polygons_2D(self, vor):
        new_regions = []
        new_vertices = vor.vertices.tolist()

        center = vor.points.mean(axis=0)
        radius = vor.points.ptp().max()
        #radius = max(self.width, self.height)
        # Construct a map containing all ridges for a given point
        all_ridges = defaultdict(list)
        #construct map of all regions and all ridge points.
        #foreach pair of points that form a ridge. extract the ridge points and vertices.
        for (p1, p2), (v1, v2) in zip(self.vor.ridge_points, self.vor.ridge_vertices):
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

        return new_regions, np.asarray(new_vertices)

    def draw(self):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        voronoi_plot_2d(self.vor, ax=ax)   
        #delaunay_plot_2d(self.tri, ax=ax)

        fig.show()

    def draw_vor_centres(self):
        #plt.triplot(self.pts[:,0], self.pts[:,1], self.tri.simplices) 
        plt.plot(self.pts[:,0], self.pts[:,1], 'o')
        plt.show()


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

    def draw_voronoi(self):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        #plot the points
        ax.plot(self.pts[:,0], self.pts[:,1], 'o')

        #create a line collection of the line segments.
        ax.add_collection(LineCollection(self.line_segments,
                                    colors='k',
                                    lw=1.0,
                                    alpha=1.0,
                                    linestyle='solid'))

        fig.show()


    def draw_land(self):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        
        vor = self.vor

        height = np.random.uniform(low=0.0, high=5.0, size=self.num_points)

        # find min/max values for normalization
        minima = min(height)
        maxima = max(height)

        # normalize chosen colormap
        norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
        mapper = cm.ScalarMappable(norm=norm, cmap=cm.Blues_r)

        
        #create a line collection of the line segments.
        ax.add_collection(LineCollection(self.line_segments,
                                    colors='k',
                                    lw=1.0,
                                    alpha=1.0,
                                    linestyle='solid'))

        for i, cell in enumerate(self.cells):
            ax.fill(*zip(*cell.region_verts), color=mapper.to_rgba(height[i]))

        plt.xlim(vor.min_bound[0] - 0.1, vor.max_bound[0] + 0.1)
        plt.ylim(vor.min_bound[1] - 0.1, vor.max_bound[1] + 0.1)
        fig.show()


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
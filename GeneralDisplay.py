__author__ = 'Leo Cao'

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay

# significance level for permutation test on median velocity of culvert and surrounding area
PERM_ALPHA = 0.01
TRI_THRESHOLD = 2.0


# ############# plots based on rating #############################################################
def culvert_info_display(results):
    inner_results = results["inner"]
    surrounding_results = results["surrounding"]
    culverts = results["culverts"]
    search_radius = results["radius"]
    counter = 0
    total_num_points_in_search_radius = 0
    for culvert in culverts:
        circle = plt.Circle((culvert.x, culvert.y), search_radius, color='r')
        figure1 = plt.figure(1)
        axe = figure1.gca()

        if culvert.rate > 5:
            circle.set_color('g')

        elif culvert.rate > 4:
            circle.set_color('y')

        else:
            circle.set_color('r')

        if 0 < culvert.perm_test_p_value < PERM_ALPHA:
            axe.annotate('Incoherent with Surrounding',
                         xy=(culvert.x, culvert.y),
                         xytext=(culvert.x + 1.5 * search_radius, culvert.y + 2 * search_radius)
                         )
            print culvert.federal_structure_id, "Incoherent Movement compared to Surrounding "

        if culvert.tri_result > TRI_THRESHOLD:
            axe.annotate('Surface Change',
                         xy=(culvert.x, culvert.y),
                         xytext=(culvert.x + 1.5 * search_radius, culvert.y + 3 * search_radius)
                         )
            print culvert.federal_structure_id, " Surface Fluctuation Warning"

        axe.add_artist(circle)
        inner_result = inner_results[culvert.federal_structure_id]
        surrounding_result = surrounding_results[culvert.federal_structure_id]
        if inner_result.size != 0:
            counter += 1

            total_num_points_in_search_radius += inner_result.size

            inner_xs = [p[0] for p in inner_result]
            inner_ys = [p[1] for p in inner_result]
            axe.plot(inner_xs, inner_ys, "r.")

            surrounding_xs = [p[0] for p in surrounding_result]
            surrounding_ys = [p[1] for p in surrounding_result]
            axe.plot(surrounding_xs, surrounding_ys, "k.")

    plt.show()

    print "There are ", len(culverts), " Culverts"
    print counter, "culverts has data points within ", search_radius, " feet"
    print "On average, each culvert has ", total_num_points_in_search_radius/counter, " nearby data points"


def triangulation_result_display(culvert, inner_result, all_points, features, radius, plot_num=10):
    background_color = 0
    ft2mm = 304.8

    inner_points = all_points[inner_result]

    culvert_mean_location = inner_points.mean(axis=0)
    points = inner_points - culvert_mean_location
    # print points
    tri = Delaunay(points)
    x = [p[0] for p in points]
    y = [p[1] for p in points]

    tri_num = tri.simplices.shape[0]
    zfaces = np.array([i for i in range(tri_num)])
    # for each triangle calculate area change
    # calculation unit: mm
    area_sum = 0
    area_extreme_sum = 0
    area_sum_ts = [0] * tri_num
    for i in range(0, tri_num):
        # zfaces.append(0)
        # get 3 vertices of each tri simplice
        tri_i = tri.simplices[i, :]
        p0 = points[tri_i[0]]
        p1 = points[tri_i[1]]
        p2 = points[tri_i[2]]

        # restore to original location
        p0_ori = p0 + culvert_mean_location
        p1_ori = p1 + culvert_mean_location
        p2_ori = p2 + culvert_mean_location

        feature0 = features[(p0_ori[0], p0_ori[1])]
        feature1 = features[(p1_ori[0], p1_ori[1])]
        feature2 = features[(p2_ori[0], p2_ori[1])]

        p0 = np.append(p0, feature0["HEIGHT"]*1000)
        p1 = np.append(p1, feature1["HEIGHT"]*1000)
        p2 = np.append(p2, feature2["HEIGHT"]*1000)

        p0[0] *= ft2mm
        p0[1] *= ft2mm
        p1[0] *= ft2mm
        p1[1] *= ft2mm
        p2[0] *= ft2mm
        p2[1] *= ft2mm

        height_ts0 = []
        height_ts1 = []
        height_ts2 = []
        for col in feature0.keys():
            if col[0] == 'D' and col[1] == "2" and col[2] == "0":
                height_ts0.append(feature0[col])
                height_ts1.append(feature1[col])
                height_ts2.append(feature2[col])

        # edge length
        e0 = np.linalg.norm(p0 - p1)
        e1 = np.linalg.norm(p0 - p2)
        e2 = np.linalg.norm(p1 - p2)

        # calculate init area
        # calculated with (x, y, HEIGHT)  in feet^2
        area_init = np.sqrt((e0 + e1 + e2) * (e0 + e1 - e2) * (e0 + e2 - e1) * (e1 + e2 - e0)) / 4

        # get max deviation position of each vertex p0, p1, p2
        d0 = max(height_ts0) if abs(max(height_ts0)) > abs(min(height_ts0)) else min(height_ts0)
        d1 = max(height_ts1) if abs(max(height_ts1)) > abs(min(height_ts1)) else min(height_ts1)
        d2 = max(height_ts2) if abs(max(height_ts2)) > abs(min(height_ts2)) else min(height_ts2)

        p0_extreme = p0
        p1_extreme = p1
        p2_extreme = p2

        p0_extreme[2] += d0
        p1_extreme[2] += d1
        p2_extreme[2] += d2

        e0_extreme = np.linalg.norm(p0 - p1)
        e1_extreme = np.linalg.norm(p0 - p2)
        e2_extreme = np.linalg.norm(p1 - p2)

        area_extreme = np.sqrt((e0_extreme + e1_extreme + e2_extreme) * (e0_extreme + e1_extreme - e2_extreme)
                               * (e0_extreme + e2_extreme - e1_extreme) *
                               (e1_extreme + e2_extreme - e0_extreme)) / 4

        area_sum += area_init
        area_extreme_sum += area_extreme

        # calculate a ts of area
        area_diff_ts = [0] * len(height_ts0)
        # print "hh", height_ts0, "h"
        for j in range(0, len(height_ts0)):
            p0_j = p0
            p1_j = p1
            p2_j = p2

            d0_j = height_ts0[j]
            d1_j = height_ts1[j]
            d2_j = height_ts2[j]

            p0_j[2] += d0_j
            p1_j[2] += d1_j
            p2_j[2] += d2_j

            e0_j = np.linalg.norm(p0 - p1)
            e1_j = np.linalg.norm(p0 - p2)
            e2_j = np.linalg.norm(p1 - p2)

            area_j = np.sqrt((e0_j + e1_j + e2_j) * (e0_j + e1_j - e2_j) * (e0_j + e2_j - e1_j) *
                             (e1_j + e2_j - e0_j)) / 4
            area_diff_ts[j] = area_j-area_init

        temp = max(area_diff_ts) if abs(max(area_diff_ts)) > abs(min(area_diff_ts)) else min(area_diff_ts)
        area_sum_ts[i] = temp
        zfaces[i] = 100 * temp/area_init

    change_percentage = 100 * sum(area_sum_ts) / area_sum
    if 100 * sum(area_sum_ts) / area_sum > TRI_THRESHOLD:
        fig, ax = plt.subplots()
        ax.triplot(x, y, tri.simplices.copy(), "go-")

        color = np.clip(zfaces, 0, 50)
        plt.tripcolor(x, y, tri.simplices.copy(), facecolors=color, edgecolors="k")

        levels = np.arange(0, 50, 5)
        plt.colorbar(ticks=levels)

        ax.set_xlabel('Feet')
        ax.set_ylabel('Feet')

        plt.title("warning on Structure:" + str(int(culvert.federal_structure_id))
                  + " Current rating:" + str(culvert.rate) + " Surface Change " + str(change_percentage))

    plt.draw()

    return change_percentage, tri


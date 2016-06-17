__author__ = 'Leo Cao'
import ogr
import osr
from ProcessCulvertExcel import get_culverts_from_excel
from scipy import spatial
import os
import math
from GeneralDisplay import *

# ######################### Control Variables ################################
OUTER_RADIUS_RATIO = 6
SHOW_MESSAGE = False

# the number of sampling for permutation test
PERM_NUM = 10000

# ################# Helper function for permutation test using monte-carlo method #################
def exact_mc_perm_test(xs, ys, nmc):
    n, k = len(xs), 0.0
    diff = np.abs(np.median(xs) - np.median(ys))
    zs = np.concatenate([xs, ys])
    for j in range(nmc):
        np.random.shuffle(zs)
        k += diff < np.abs(np.median(zs[:n]) - np.median(zs[n:]))
    return k / nmc


def data_summary(shapefile_path, displacement_date):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    # input data source
    input_data = driver.Open(shapefile_path)
    in_layer = input_data.GetLayer()

    vel_list = []
    v_var_list = []
    displacement_list = []

    for feature in in_layer:
        if feature["EFF_AREA"] != 0:
            continue
        else:
            vel_list.append(feature["VEL"])
            v_var_list.append(feature["V_STDEV"] * feature["V_STDEV"])
            displacement_list.append(feature[displacement_date])

    avg_vel = sum(vel_list)/len(vel_list)
    avg_v_std = math.sqrt(sum(v_var_list)/len(v_var_list))
    avg_displacement = sum(displacement_list)/len(displacement_list)
    return avg_vel, avg_v_std, avg_displacement


def extract_points_in_culvert(culvert_excel_path, shapefile_path, out_shape_file_name, tri_shape_file="tri.shp",
                              input_search_radius=50, input_time_range=-1, input_displacement_threshhold=0.5):

    search_radius = input_search_radius

    # #################### process input excel ################################
    culverts = get_culverts_from_excel(culvert_excel_path)

    # ######################## input shapefile #########################
    # message contains shapeFile information
    message = ""

    # create input driver
    driver = ogr.GetDriverByName('ESRI Shapefile')
    # input data source
    input_data = driver.Open(shapefile_path)
    in_layer = input_data.GetLayer()
    message += "input shape file contains " + str(in_layer.GetFeatureCount()) + " features\n"

    # get current coordinate system
    source_spacial_reference = in_layer.GetSpatialRef()
    message += "------ get_points_from_shape_file ------\n"
    message += "the Spatial Reference of the shapes is: " + str(source_spacial_reference) + "\n"

    # new spatial reference in lat/long ---- geographical coordinate system
    lat_long_sr = osr.SpatialReference()
    # lat_long_sr.SetWellKnownGeogCS("EPSG:4326")
    lat_long_sr.ImportFromEPSG(4326)

    # Create a Coordinate Transformation
    # (From the projected coordinate system:"NAD_1983_StatePlane_Virginia_North_FIPS_4501_Feet"
    # to geographic coordinate system: WGS84/EPSG:4326 )
    # projection_coordinate_transformation = osr.CreateCoordinateTransformation(source_spacial_reference, lat_long_sr)

    # from lat/long to feet in "NAD_1983_StatePlane_Virginia_North_FIPS_4501_Feet"
    lat_long_to_feet_transformation = osr.CreateCoordinateTransformation(lat_long_sr, source_spacial_reference)

    # ##################### create output shapefile #################################
    out_driver = ogr.GetDriverByName("ESRI Shapefile")

    out_shape_file = out_shape_file_name
    if __name__ == "__main__":
        print "out_shape_file", out_shape_file

    # Remove output shape file if it already exists
    if os.path.exists(out_shape_file):
        out_driver.DeleteDataSource(out_shape_file)

    # Create the output shape file
    out_data_source = out_driver.CreateDataSource(out_shape_file)

    out_lyr_name = os.path.splitext(os.path.split(out_shape_file)[1])[0]
    out_layer = out_data_source.CreateLayer(out_lyr_name, srs=source_spacial_reference, geom_type=ogr.wkbPoint)

    # Add input Layer Fields to the output Layer if it is the one we want
    in_layer_defn = in_layer.GetLayerDefn()
    for i in range(0, in_layer_defn.GetFieldCount()):
        out_field_defn = in_layer_defn.GetFieldDefn(i)
        out_layer.CreateField(out_field_defn)

    id_field = ogr.FieldDefn("FED_ID", ogr.OFTInteger)
    out_layer.CreateField(id_field)

    # ######################## create shapefile for culvert surface Triangles #########################
    tri_driver = ogr.GetDriverByName("ESRI Shapefile")

    if __name__ == "__main__":
        print "tri_shape_file", tri_shape_file

    # Remove output shape file if it already exists
    if os.path.exists(tri_shape_file):
        tri_driver.DeleteDataSource(tri_shape_file)

    # Create the output shape file
    tri_data_source = tri_driver.CreateDataSource(tri_shape_file)

    tri_lyr_name = os.path.splitext(os.path.split(tri_shape_file)[1])[0]
    tri_layer = tri_data_source.CreateLayer(tri_lyr_name, srs=source_spacial_reference, geom_type=ogr.wkbPolygon)
    tri_layer.CreateField(id_field)

    area_change_field = ogr.FieldDefn("AreaChange", ogr.OFTReal)
    tri_layer.CreateField(area_change_field)

    # ######################################## Creating FeaturePoint Objects ###########################################
    features = {}
    all_features = {}
    counter = 0

    feature_grouped_by_culvert_inner = {}
    feature_grouped_by_culvert_surrounding = {}

    for feature in in_layer:

        geo = feature.GetGeometryRef()

        if geo is None:
            continue
        # ignore DS points in output shapefile
        if feature["EFF_AREA"] != 0:
            continue
        new_point = (geo.GetX(), geo.GetY())
        features[new_point] = dict(HEIGHT=feature["HEIGHT"], VEL=feature["VEL"],
                                   ACC=feature["ACC"], V_STD=feature["V_STDEV"], Feature=feature)
        all_features[new_point] = feature
        counter += 1

    # ########################################### Construct KDT ########################################################
    all_points = np.array([list(i) for i in features.keys()])
    tree = spatial.KDTree(all_points)

    message += "------ End of get_points_from_shape_file------\n"
    if SHOW_MESSAGE:
        print message

    # ############################# Getting General Statistics for points in Culverts ##################################

    # avg_vel, avg_v_std, avg_displacement = data_summary(shapefile_path, "D20141124")

    # ##################################### Getting Nearest Points for Culverts ########################################

    counter = 0
    total_num_points_in_search_radius = 0
    num_culvert_has_point = 0

    # for each culvert in input file
    total_culvert_sum_vel = 0
    total_culvert_sum_v_var = 0
    total_outer_sum_vel = 0
    total_outer_sum_v_var = 0
    outer_num_data = 0

    change_percent_summary = []
    perm_p_value_summary = []

    for culvert in culverts:

        culvert_sum_vel = 0
        culvert_sum_v_var = 0
        # permutation test result on median velocity of culvert and surrounding area
        # -1: test cannot be conducted
        # positive number: the p-value of permutation test result
        perm_result = -1

        # Transform the lat/long of culvert into VA's feet coordinate system
        pt = ogr.Geometry(ogr.wkbPoint)
        pt.AddPoint(culvert.lg, culvert.lat)
        pt.Transform(lat_long_to_feet_transformation)
        culvert.x = pt.GetX()
        culvert.y = pt.GetY()
        # perform query with input search radius in KDT
        # and perform another search with twice radius, for comparison
        inner_result = tree.query_ball_point([pt.GetX(), pt.GetY()], search_radius)
        outer_result = tree.query_ball_point([pt.GetX(), pt.GetY()], search_radius * OUTER_RADIUS_RATIO)

        surrounding_result = list(set(outer_result) - set(inner_result))

        feature_grouped_by_culvert_inner[culvert.federal_structure_id] = all_points[inner_result]
        total_num_points_in_search_radius += len(inner_result)
        feature_grouped_by_culvert_surrounding[culvert.federal_structure_id] = all_points[surrounding_result]

        # get output layer's field definition
        out_layer_defn = out_layer.GetLayerDefn()

        # outer circle stat
        outer_sum_vel = 0
        outer_sum_var = 0
        if outer_result:
            for p in all_points[outer_result]:
                in_feature = features[tuple(p)]["Feature"]
                outer_sum_vel += in_feature["VEL"]
                outer_sum_var += in_feature["V_STDEV"] ** 2
            outer_num_data += len(outer_result)
            outer_avg_vel = outer_sum_vel/len(outer_result)
            outer_v_std = math.sqrt(outer_sum_var/len(outer_result))
            total_outer_sum_vel += outer_sum_vel
            total_outer_sum_v_var += outer_sum_var

            # conduct permutation test on median velocity of culvert and surrounding area only if there are enough data
            if len(surrounding_result) > 3 and len(inner_result) > 3:
                culvert_vel_list = []
                surrounding_vel_list = []
                for p in all_points[inner_result]:
                    feature_i = features[tuple(p)]["Feature"]
                    culvert_vel_list.append(feature_i["VEL"])
                for p in all_points[surrounding_result]:
                    feature_j = features[tuple(p)]["Feature"]
                    surrounding_vel_list.append(feature_j["VEL"])
                perm_result = exact_mc_perm_test(culvert_vel_list, surrounding_vel_list, PERM_NUM)
                culvert.perm_test_p_value = perm_result
                perm_p_value_summary.append(perm_result)

        # retrieve and store all feature points with in search radius into output layer
        if inner_result:
            num_culvert_has_point += 1
            features_in_culvert = {}
            for p in all_points[inner_result]:
                # add all result InSAR feature points withing culvert radius to output shape file
                out_feature = ogr.Feature(out_layer_defn)
                in_feature = features[tuple(p)]["Feature"]
                features_in_culvert[tuple(p)] = in_feature
                out_feature.SetGeometry(in_feature.GetGeometryRef())

                for i in range(0, in_layer_defn.GetFieldCount()):
                    out_feature.SetField(out_layer_defn.GetFieldDefn(i).GetNameRef(), in_feature.GetField(i))

                out_feature["FED_ID"] = culvert.federal_structure_id
                out_layer.CreateFeature(out_feature)
                culvert_sum_vel += in_feature["VEL"]
                total_culvert_sum_vel += in_feature["VEL"]
                culvert_sum_v_var += in_feature["V_STDEV"] ** 2
                total_culvert_sum_v_var += in_feature["V_STDEV"] ** 2

            # ###################### Triangulation Analysis ###########################################################
            if len(inner_result) > 4:
                counter += 1
                change_percent, tri = \
                    triangulation_result_display(culvert=culvert, inner_result=inner_result, all_points=all_points,
                                                 features=features_in_culvert, radius=input_search_radius,
                                                 plot_num=counter)
                change_percent_summary.append(change_percent)

                inner_points = all_points[inner_result]
                culvert_mean_location = inner_points.mean(axis=0)
                points = inner_points - culvert_mean_location
                tri_num = tri.simplices.shape[0]

                for i in range(0, tri_num):
                    # get 3 vertices of each tri simplice
                    tri_i = tri.simplices[i, :]
                    p0 = points[tri_i[0]]
                    p1 = points[tri_i[1]]
                    p2 = points[tri_i[2]]

                    # restore to original location
                    p0_ori = p0 + culvert_mean_location
                    p1_ori = p1 + culvert_mean_location
                    p2_ori = p2 + culvert_mean_location

                    ring = ogr.Geometry(ogr.wkbLinearRing)
                    ring.AddPoint(p0_ori[0], p0_ori[1])
                    ring.AddPoint(p1_ori[0], p1_ori[1])
                    ring.AddPoint(p2_ori[0], p2_ori[1])
                    poly = ogr.Geometry(ogr.wkbPolygon)
                    poly.AddGeometry(ring)
                    tri_feature = ogr.Feature(tri_layer.GetLayerDefn())
                    tri_feature.SetGeometry(poly)
                    tri_feature["FED_ID"] = culvert.federal_structure_id
                    tri_feature["AreaChange"] = change_percent
                    tri_layer.CreateFeature(tri_feature)

                # update culvert attribute
                culvert.tri_result = change_percent

    out_data_source.Destroy()
    print "total number of points within radius ", total_num_points_in_search_radius
    print "num of culverts: ", num_culvert_has_point
    print "on average each culvert has ", total_num_points_in_search_radius/num_culvert_has_point, " points"

    plt.figure()
    plt.hist(change_percent_summary)
    plt.figure()
    plt.hist(perm_p_value_summary)

    plt.show()

    return {"all_points": all_points, "inner": feature_grouped_by_culvert_inner,
            "surrounding": feature_grouped_by_culvert_surrounding, "culverts": culverts,
            "radius": input_search_radius, "features": features}

if __name__ == "__main__":
    in_shapefile_path = \
        r"D:/VIVA/Data_from_VIVA_Server/Data/TRE Data/RITARS-14-H-UVA/Site 1/" \
        r"RITA2_Site1_Staunton_SqueeSAR_Nov2014/SqueeSAR_Data/SITE1_STAUNTON_RITA2_CSK_H40B_A_T213_" \
        r"C1_Nov14-TSR.shp"

    radius = 50
    results = extract_points_in_culvert(r"Cleaned_STAUNT&NOVA CULV.xlsx",
                                        in_shapefile_path, r"Output.shp", input_search_radius=radius)

    culvert_info_display(results)



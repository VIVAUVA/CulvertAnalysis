__author__ = 'Leo Cao'
import ogr


class Culvert:

    def __init__(self, lat, lg, rate, age=-1, federal_structure_id=-1):
        self.federal_structure_id = federal_structure_id
        self.lat = lat
        self.lg = lg
        self.age = age
        self.rate = rate
        self.x = 0
        self.y = 0
        self.tri_result = 0
        self.perm_test_p_value = -1
        # create a Geometry object of type Point based on lat/long
        self.point = ogr.Geometry(ogr.wkbPoint)
        self.point.AddPoint(self.lat, self.lg)
        self.neighbor = []

    def __str__(self):
        return "culvertID: " + str(self.federal_structure_id) + ", lat: " + str(self.lat) + \
               ", long: " + str(self.lg) + ", age: " + str(self.age) + ", rate: " + str(self.rate) + \
               ", # points near by: " + str(self.neighbor)

    @property
    def lat(self):
        return self.lat

    @property
    def lg(self):
        return self.lg

    @property
    def point(self):
        return self.point


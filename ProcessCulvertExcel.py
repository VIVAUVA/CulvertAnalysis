__author__ = 'Leo Cao'
from Culvert import Culvert
import xlrd


def get_culverts_from_excel(excel_name):
    # message contains excel information
    message = ""

    # the function reads in data from excel spreadsheet with name passed in as parameter excel_name
    # the function return a list of Culvert objects defined in Culvert.py
    work_book = xlrd.open_workbook(excel_name)
    sheet_names = work_book.sheet_names()
    # the culverts' information is expected to be on the first spread sheet
    culverts_sheet = work_book.sheet_by_name(sheet_names[0])
    # the culvert sheet are expected to have header on the first row
    header = culverts_sheet.row_values(0)

    message += "------ get_culverts_from_excel -----\n"
    message += "shread sheet header: \n"
    message += str(header) + "\n"

    # the headers of culvert id, latitude, longitude, age and rating are expected to be same as shown below
    id_col = header.index('Federal Structure ID')
    lat_col = header.index('Latitude')
    long_col = header.index('Longitude')
    rate_col = header.index('Culvert Rating')

    culverts = []

    for i in range(1, culverts_sheet.nrows-1):
        entry = culverts_sheet.row_values(i)
        federal_structure_id = entry[id_col]
        lat = entry[lat_col]
        lg = entry[long_col]
        # age = entry[age_col]
        rate = entry[rate_col]
        # check if lat or long is null
        if lat and lg:
            new_culvert = Culvert(lat, lg, rate, federal_structure_id=federal_structure_id)
            culverts.append(new_culvert)

    message += "------ End of get_culverts_from_excel ------\n"
    # print message
    return culverts


def get_culverts_from_excel_old(excel_name):
    # message contains excel information
    message = ""

    # the function reads in data from excel spreadsheet with name passed in as parameter excel_name
    # the function return a list of Culvert objects defined in Culvert.py
    work_book = xlrd.open_workbook(excel_name)
    sheet_names = work_book.sheet_names()
    # the culverts' information is expected to be on the first spread sheet
    culverts_sheet = work_book.sheet_by_name(sheet_names[0])
    # the culvert sheet are expected to have header on the first row
    header = culverts_sheet.row_values(0)

    message += "------ get_culverts_from_excel -----\n"
    message += "shread sheet header: \n"
    message += str(header) + "\n"

    # the headers of culvert id, latitude, longitude, age and rating are expected to be same as shown below
    id_col = header.index('Federal Structure ID')
    lat_col = header.index('lat')
    long_col = header.index('long')
    # age_col = header.index('Age (yrs to 2015)')
    rate_col = header.index('Culvert Rating')

    culverts = []

    for i in range(1, culverts_sheet.nrows-1):
        entry = culverts_sheet.row_values(i)
        culvert_id = entry[id_col]
        lat = entry[lat_col]
        lg = entry[long_col]
        # age = entry[age_col]
        rate = entry[rate_col]
        # check if lat or long is null
        if lat and lg:
            new_culvert = Culvert(culvert_id, lat, lg, rate)
            culverts.append(new_culvert)

    message += "------ End of get_culverts_from_excel ------\n"
    # print message
    return culverts
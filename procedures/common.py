import csv
from os.path import exists


# open an existing csv file for appending or create
# it with a given header
# returns a handler to the file and csv_data
def create_or_append_file(file_name: str, header: list, force_create=False):
    csv_data = None
    file = None
    if exists(file_name) and not force_create:
        file = open(file_name, "a")
        csv_data = csv.writer(file, delimiter=";", quotechar="'")
    else:
        file = open(file_name, "w")
        csv_data = csv.writer(file, delimiter=";", quotechar="'")
        csv_data.writerow(header)
    return file, csv_data

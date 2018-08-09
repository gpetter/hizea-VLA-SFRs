import os


def return_galaxy_list():
    data_path = '/users/gpetter/DATA/data_v2'
    flag_list = ['J090133.42', 'J094417.84']

    # Go to directory, get list of galaxies
    names = os.listdir(data_path)
    os.chdir(data_path)

    # Sort by RA
    stripped_names = []
    sorted_names = []

    for name in names:
        new_name = name.split('J')[1]
        stripped_names.append(new_name)

    sorted_nums = sorted(stripped_names)

    for thing in sorted_nums:
        sorted_names.append('J' + thing)

    for x in range(len(flag_list)):
        sorted_names.remove(flag_list[x])

    return sorted_names

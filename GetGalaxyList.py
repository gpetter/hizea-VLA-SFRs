import os

def return_galaxy_list():
    #data_path = '/lustre/aoc/students/gpetter/imaging'
    data_path = '/Users/graysonpetter/Desktop/mac_copy'
    flag_list = ['J090133.42']

    # Go to directory, get list of galaxies
    names = os.listdir(data_path)
    os.chdir(data_path)
    sorted_names = sorted([f for f in os.listdir(data_path) if not f.startswith('.')], key=lambda f: f.lower())

    for x in range(len(flag_list)):
	if flag_list[x] in sorted_names:
        	sorted_names.remove(flag_list[x])

    return sorted_names

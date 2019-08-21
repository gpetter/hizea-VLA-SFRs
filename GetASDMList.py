import os

names2 = []

def return_galaxy_list():
    data_path = '/lustre/aoc/students/gpetter/flagging/'

    # Go to directory, get list of galaxies
    names = os.listdir(data_path)
    for name in names:
	if name[-4:]=='pipe':
		names2.append(data_path+name)
	

    return names2

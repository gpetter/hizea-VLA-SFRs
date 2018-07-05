import os


def get_imfit(gal_name):

    os.chdir(gal_name)

    with open('summary.log', 'r') as f:
        lines = f.readlines()
        val_line = lines[2].split(' ')
        flux = val_line[8]
        flux_err = val_line[9]
        max_val = val_line[10]
        max_err = val_line[11]

    output = []
    output.append(flux)
    output.append(flux_err)
    output.append(max_val)
    output.append(max_err)

    os.chdir('..')

    return output
import os
import stat
import GetASDMList
reload(GetASDMList)

vises = []

current_dir = os.getcwd()
names = GetASDMList.return_galaxy_list()






# generates the pipeline script
os.chdir(current_dir)
with open('flaggingrun', 'w') as f:
    for x in range(len(names)):
        f.write("""cd %s; rm -rf oussid*; rm -rf *.last; rm -rf *.ms*; rm -rf *.g*; rm -rf *.b*; rm -rf *.k*; rm -rf *.log; rm -rf *.png; rm -rf *.csv; rm -rf *.flag*; rm -rf pipeline*; xvfb-run -d casa-pipe --nogui -c %s\n""" % (names[x], (names[x]+'/casa_pipescript.py')))

st = os.stat('flaggingrun')
os.chmod('flaggingrun', st.st_mode | 0111)

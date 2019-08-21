import os
import glob

cwd = os.getcwd()


with open('plotms_command.txt', 'r') as f:
	command = f.read()

os.chdir('/lustre/aoc/students/gpetter/flagging/')

split1 = command.split('(')[1]
split2 = split1.split(')')[0]
split3 = split2.split(' ')

msname = split3[0].split('=')[1].strip(',').strip("""'""")[:30]
goob = glob.glob('%s*' % (msname))[0]
path = os.path.abspath(goob)

vis = split3[0]
index = vis.find("""'""")
line = vis[:index] + path + '/' + vis[index:][1:]


index2 = line.find('/')
line2 = line[:index2] + """'""" + line[index2:]
print(line2)

split3[0] = line2



#mspath = '/lustre/aoc/students/gpetter/flagging/%s' % ()

for x in range(len(split3)):
	if str(split3[x]).startswith('spw'):
		split3[x] = """spw='0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15',"""
for x in range(len(split3)):
	if str(split3[x]).startswith('showgui'):
		split3[x] = 'showgui=True,'

shorter = [x for x in split3 if not (x.startswith('plotfile'))]


linkup = " ".join(shorter)
final = 'plotms(' + linkup + ')'

print(final)

os.chdir(cwd)

exec(final)

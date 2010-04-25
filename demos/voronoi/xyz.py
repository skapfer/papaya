
import fileio


def read_xyz_with_comment (infilename):
	fp = open (infilename)
	#numpoints = fp.readline ()
	size = float(fp.readline () )
	comment = fp.readline ()
	label = 0
	ret = []
	for l in fp.readlines ():
		bl = fileio.break_line_into_string (l)
		x = float (bl[1])
		y = float (bl[2])
		R = 0.
		ret.append ([label, x, y, R])
		label += 1
	fp.close ()
	return (ret, comment)

def read_xyz (infilename):
    return read_xyz_with_comment (infilename)[0]

def split_comment (c):
    ca = c.split (',')
    ret = {}
    for x in ca:
        x = x.split ('=')
        if len (x) == 2:
            ret[x[0].strip ()] = x[1].strip ()
    return ret

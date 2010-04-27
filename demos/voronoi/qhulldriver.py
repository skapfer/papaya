import os
import tempfile
import poly

def get_temp_file_name (suffix = ''):
    (number, name) = tempfile.mkstemp (suffix = suffix)
    os.close (number)
    return name

def qvoronoi (points):

    point_fn = get_temp_file_name ()
    fp = open (point_fn, "w")
    fp.write ("2 rbox %i\n%i D2\n" % (len(points), len(points)))
    for x in points:
        fp.write ("%.16f %.16f\n" % (x[0],x[1]))
    fp.close ()

    off_fn = get_temp_file_name ();
    print off_fn
    cmd = 'cat %s | qvoronoi s o >%s' % (point_fn, off_fn)
    print cmd
    if os.system (cmd) != 0:
        die ("qvoronoi failed")

    tessel = poly.poly ()
    tessel.read_from_off (off_fn)

    os.remove (off_fn)
    #os.remove (point_fn)

    return tessel
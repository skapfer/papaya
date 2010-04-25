def add (x, y):
    assert len (x) == len (y)
    if len (x) == 2:
        return (x[0]+y[0], x[1]+y[1])
    else:
        die ("NYI")

def add_periodic_copies (points, box, lattice_vectors, overlap_factor):
    max_dist = overlap_factor * max (box)
    ret = []
    for L in lattice_vectors:
        for P in points:
            np = add (L, P)
            if not (np[0] > -max_dist and np[1] > -max_dist and\
                np[0]-box[0] < max_dist and np[1]-box[1] < max_dist):

                continue

            ret.append (np)

    return ret

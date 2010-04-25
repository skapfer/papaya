# ***************************************************************************
# break line into strings
#
# SPECIAL VERSION 
#
# this function breaks a string 'line' as read by
# line=file.readline()
# into separate strings where it recognises blanc spaces, '(', ')' and  ':' as
# separators between these strings.
# "asdfsd 234 dfsdf   234" -> {"asdfsd","234","dfsdf","234"}
#
# A line read from a file probably ends with a carriage return,
# this will not be returned
#
def break_line_into_string(line):
    list=[]
    length=len(line)
    i=0
    if line[len(line)-1] == '\012':
        line = line[:-1]
    while i < len(line):
        while i<len(line):
            break_now = 1
            if line[i]==' ' or line[i] == ':' or line[i] == '\t' or \
	       line[i]=='(' or line[i] == ')' or line[i] == ',' or\
	       line[i]=='[' or line[i] == ']' or line[i] == '=':
                break_now = 0
                i=i+1
            if break_now == 1:
                break
        word=""
        while i<len(line):
            if line[i] == ' ' or line[i] == ':'  or line[i] == '\t' or\
	       line[i] == '(' or line[i] == ')' or line[i] == ',' or\
	       line[i]=='[' or line[i] == ']' or line[i] == '=':
                break
            word=word+line[i]
            i=i+1
        if word != "":
            list.append(word)
    return list
# ***************************************************************************

from string import split
from string import atoi
from string import atof
from sys import argv
from os import system

# ***************************************************************************
# READ POLY FILE
def ReadPoly(filename,factor=1):
    po=open(filename,'r')
    # read first line, should be points
    line = po.readline()
    if line!='POINTS\012':
	error('I dont know what this file is')
    lines = po.readlines()
    # analyse first line of file and determine what information is contained
    # in each line
    line=lines[0]
    s = break_line_into_string(line)
    j = 0
    uv_i = -1
    n_i  = -1
    co_i  = -1
    while j < len(s):
        if s[j] == "uv":
            uv_i = j+1
        if s[j] == "n":
            n_i  = j+1
        if s[j] == "c":
            co_i = j+1
        j = j+1
    points_coord=[]
    points_index=[]
    findex2index={}
    uv=[]
    normals = []
    co = []
    # read points until 'POLYS' is found
    for i in range(0,len(lines)):
	line=lines[i]
	s = break_line_into_string(line)
	if s[0]=='POLYS':
		pos_poly = i
		break
        # assign point coordinates
	points_coord.append([atof(s[1]),atof(s[2]),atof(s[3])])
	points_index.append(atoi(s[0]))
        findex2index[atoi(s[0])]=len(points_coord)-1
        if n_i != -1:
            normals.append([atof(s[n_i]),atof(s[n_i+1]),atof(s[n_i+2])])
        if uv_i != -1:
            uv.append([atof(s[uv_i]),atof(s[uv_i+1])])
        if co_i != -1:
            co.append([atof(s[co_i]),atof(s[co_i+1]),\
                       atof(s[co_i+2]),atof(s[co_i+3])])
    # divide by factor
    if factor != 1:
        for i in range(len(points_coord)):
            for j in range(3):
                points_coord[i][j] = points_coord[i][j]/factor
                if n_i != -1:
                    normals[i][j] = normals[i][j]/factor
            for j in range(4):
                if co_i != -1:
                    co[i][j] = co[i][j]/factor
            for j in range(2):
                if uv_i != -1:
                    uv[i][j] = uv[i][j]/factor
                
    # read polygons
    polys_verts = []
    polys_index = []
    for i in range(pos_poly+1,len(lines)):
	line=lines[i]
	s=break_line_into_string(line)
	if s[0]=="END":
		pos_comments = i+1
		break
	# determine wether line ends with '<' (=closed)
	last_vertex=len(s)
        while s[len(s)-1] == " ":
            del s[len(s)-1]
	if s[len(s)-1]=="<":
		last_vertex = len(s)-1
        polys_verts.append([])
        polyindex=len(polys_verts)-1
        polys_index.append(atoi(s[0]))
        for j in range(1,last_vertex):
            if j > len(s):
                print "Error in ReadPoly: j > len(s). "
                print "j, len(s), i : ", j, ", ", len(s), ", ", i
            if findex2index.has_key(atoi(s[j])) == 0:
                print "Error in ReadPoly: findex2index.has_key(atoi(s[j])) == 0"
                print "atoi(s[j]), len(findex2index), s[j], i : ", atoi(s[j]), ", ", len(findex2index), \
                      ", ", s[j], ", ", i
            polys_verts[polyindex].append(findex2index[atoi(s[j])])
    # read comments
    if pos_comments < len(lines):
        comments=lines[pos_comments:len(lines)]
    else:
        comments=[]
    # return list
    po.close()
    return([points_index,points_coord,polys_index,polys_verts,normals,uv,co,comments])
# ***************************************************************************



# ***************************************************************************
def PrintOff(filename,coords,polys,corners=[]):
    outfile=open(filename,'w')
    outfile.write("OFF \n")
    outfile.write("%s %s 0 \n" % (len(coords)+2,len(polys)))
    # print vertices
    for coord in coords:
        x=coord[0]
        y=coord[1]
        z=coord[2]
        outfile.write("%f %f %f\n" % (x,y,z))
    # print two extra vertices in corners
    if corners==[]:
        [xmm,ymm,zmm]=coords[0]
        [xmp,ymp,zmp]=coords[0]
    outfile.write("%f %f %f\n" % (xmm,ymm,zmm))
    outfile.write("%f %f %f\n" % (xmp,ymp,zmp))
    # print polygons
    for polygon in polys:
        outfile.write("%s " % len(polygon))
        for vertex in polygon:
            outfile.write("%s " % vertex)
        outfile.write("\n")
    outfile.close()
# ***************************************************************************


# ***************************************************************************
def PrintPoly(filename,coords,polys,nnc=[],uvc=[],coc=[],factor=1,comments=[]):
    outfile=open(filename,'w')
    outfile.write("POINTS\n")
    # print vertices
    for i in range(len(coords)):
        x=coords[i][0]
        y=coords[i][1]
        z=coords[i][2]
        if factor != 1:
            x=int(x*factor)
            y=int(y*factor)
            z=int(z*factor)
        outfile.write("%s: %f %f %f " % (i+1,x,y,z))
        if len(coc) == len(coords):
            [cr,cg,cb,al]=coc[i]
            if factor!=1:
                [cr,cg,cb,al]=[int(cr*factor),int(cg*factor),int(cb*factor),int(al*factor)]
            outfile.write("c(%f,%f,%f,%f) " % (cr,cg,cb,al))
        if len(uvc) == len(coords):
            [u,v]=uvc[i]
            if factor != 1:
                [u,v]=[int(u*factor),int(v*factor)]
            outfile.write("uv(%f,%f) " % (u,v))
        if len(nnc) == len(coords):
            [nx,ny,nz]=nnc[i]
            if factor != 1:
                [nx,ny,nz]=[int(factor*nx),int(factor*ny),int(factor*nz)]
            outfile.write("n(%f,%f,%f) " % (nx,ny,nz))
        outfile.write("\n")
    # print polygons
    outfile.write("POLYS\n")
    for i in range(len(polys)):
	polygon=polys[i]
        outfile.write("%s: " % (i+1))
        for vertex in polygon:
            outfile.write("%s " % (vertex+1))
        outfile.write("<\n")
    outfile.write("END\n")
    # print comments
    for comment in comments:
        outfile.write(comment+"\n")
    outfile.close()
# ***************************************************************************

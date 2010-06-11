from shared import break_line_into_string

class poly:

    # constructor 
       def __init__(self):
              self.coords =[]                            # (label, x,y,z)
              self.polys = []                            # (label, first_index, ... , last_index)
       
       def read_from_off(self,infilename):
              from string import atoi, atof              
              offfile=open(infilename,'r')
              lines = offfile.readlines()
              s = break_line_into_string(lines[1])
              no_of_points=atoi(s[0])
              no_of_polys=atoi(s[1])
              for i in range(2,no_of_points+2):
                    s = break_line_into_string(lines[i])
                    x = [atof (f) for f in s]
                    self.coords.append (x)
              counter=0
              for i in range(no_of_points+2,no_of_points+no_of_polys+2):
                     label=None
                     s = break_line_into_string(lines[i])
                     poly=[]
                     poly.append(label)
                     for j in range(1,len(s)):
                            poly.append(atoi(s[j]))
                     self.polys.append(poly)
                     counter=counter+1
              offfile.close()

       def reindex (self):
            new_p = []
            for i in range (len (self.polys)):
                for j in range (1, len (self.polys[i])):
                    new_p.append (self.coords[self.polys[i][j]][:])
                    self.polys[i][j] = len (new_p) - 1
            print len (self.coords)
            print len (new_p)
            self.coords = new_p



       def save_poly(self,outfilename):
              polyfile = open (outfilename,'w')
              polyfile.write ("POINTS\n")
              for i in range (0,len(self.coords)):
                     if len (self.coords[i]) == 3:
                         polyfile.write ("%s: %f %f %f\n" % (i+1,self.coords[i][0],self.coords[i][1],self.coords[i][2]))
                     else:
                         polyfile.write ("%s: %f %f %f\n" % (i+1,self.coords[i][0],self.coords[i][1],0.0))
              polyfile.write ("POLYS\n")       
              for i in range(0,len(self.polys)):
                     polyfile.write("%s: " % int(i+1))
                     for j in range(1,len(self.polys[i])):
                            polyfile.write("%s " % (self.polys[i][j]+1) )
                     polyfile.write("< c(0.0, 0.0, 0.0, %s)\n" % (self.polys[i][0]) )
              print "File >",outfilename,"< saved with",len(self.coords),"points and" , len(self.polys), "polys"              
              polyfile.write("END\n")              
              polyfile.close()
                     
       def point_in_poly(self,unlabeled_point,poly_idx):
              # test if a point is inside a closed (compact) polygon
              # Ref: http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/
              #       Solution 1 
              #       Algorithm by Randolph Franklin
              N=len(self.polys[poly_idx])
              poly=[]
              for i in range(1,N):
                     poly.append(self.polys[poly_idx][i])
              n=N-1
              x=unlabeled_point[0]
              y=unlabeled_point[1]
              i=0
              j=n-1
              inside=False #False = outside
                     
              while i<n:                         
                     xpi=self.coords[poly[i]][0]
                     xpj=self.coords[poly[j]][0]
                     ypi=self.coords[poly[i]][1]
                     ypj=self.coords[poly[j]][1]

                     if ( ( (( ypi <= y) and (y < ypj)) or ((ypj <= y) and (y < ypi)) ) \
                            and (x<(xpj-xpi)*(y-ypi)/(ypj-ypi)+xpi ) ):
                            if inside==True:
                                   inside=False
                            else:
                                   inside=True
                     
                     j=i
                     i=i+1
              
              return inside

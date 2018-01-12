#!/usr/bin/env python
#pbm plotting module
import numpy
def myplot(A,file):
    n,m=A.shape
    f=open(file,'w')
    f.write("P1 \n")
    f.write("#myplot module \n")
    st=str(n) + " " +str(m) + "\n"
    f.write(st)
   
    B=A.astype(bool)
    numpy.savetxt(f,B,fmt="%i")



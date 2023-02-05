import matplotlib
import matplotlib.pyplot as pylab
import numpy
from scipy.optimize import curve_fit

# data load
t,c1,dc1,c2,dc2,c3,dc3,c4,dc4,c5,dc5,m,dm,e,de=numpy.loadtxt('2wolf_3.txt',unpack=True)
T,C1,dC1,C2,dC2,C3,dC3,C4,dC4,C5,dC5,M,dM,E,dE=numpy.loadtxt('2wolfcri_3.txt',unpack=True)

L=80

t1=numpy.append(t,T)
u=numpy.append(c1,C1)
du=numpy.append(dc1,dC1)
d=numpy.append(c2,C2)
dd=numpy.append(dc2,dC2)
v=numpy.append(c3,C3)
dv=numpy.append(dc3,dC3)
q=numpy.append(c4,C4)
dq=numpy.append(dc4,dC4)
y=numpy.append(c5,C5)
dy=numpy.append(dc5,dC5)
m1=numpy.append(m,M)
dm1=numpy.append(dm,dM)
e1=numpy.append(e,E)
de1=numpy.append(de,dE)


g1=u-m1*m1+0.2
dg1=du+2*dm1*m1
g2=d-m1*m1+0.2
dg2=dd+2*dm1*m1
g3=v-m1*m1+0.2
dg3=dv+2*dm1*m1
g4=q-m1*m1+0.2
dg4=dq+2*dm1*m1
g5=y-m1*m1+0.2
dg5=dy+2*dm1*m1




pylab.subplot(2,1,1)
pylab.rc('font',size=18)
pylab.xlabel('T')
pylab.ylabel('M')
pylab.minorticks_on()
pylab.errorbar(t1, m1, dm1, linestyle='', marker='.', color= 'blue')
pylab.ylim(0,1.1)

pylab.subplot(2,1,2)
pylab.rc('font',size=18)
pylab.xlabel('T')
pylab.ylabel('g1')
pylab.minorticks_on()
pylab.errorbar(t1, g1, dg1, linestyle='', marker='.', color= 'blue')
pylab.errorbar(t1, g2, dg2, linestyle='', marker='.', color= 'green')
pylab.errorbar(t1, g3, dg3, linestyle='', marker='.', color= 'black')
pylab.errorbar(t1, g4, dg4, linestyle='', marker='.', color= 'yellow')
pylab.errorbar(t1, g5, dg5, linestyle='', marker='.', color= 'red')

pylab.ylim(0,1)

# show the plot
pylab.show()

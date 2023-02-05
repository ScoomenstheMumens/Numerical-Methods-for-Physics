import matplotlib
import matplotlib.pyplot as pylab
import numpy
from scipy.optimize import curve_fit

# data load
t1,m1,dm1,e1,de1,x1,dx1,h1,dh1=numpy.loadtxt('wolf_10.txt',unpack=True)
t2,m2,dm2,e2,de2,x2,dx2,h2,dh2=numpy.loadtxt('wolf_20.txt',unpack=True)
t3,m3,dm3,e3,de3,x3,dx3,h3,dh3=numpy.loadtxt('wolf_30.txt',unpack=True)
t4,m4,dm4,e4,de4,x4,dx4,h4,dh4=numpy.loadtxt('wolfcri_40.txt',unpack=True)
L=80

dm1=dm1/100
dm2=dm2/100
dm3=dm3/100
dm4=dm4/100
de1=de1/100
de2=de2/100
de3=de3/100
de4=de4/100

pylab.subplot(221)
pylab.rc('font',size=18)
pylab.xlabel('T')
pylab.ylabel('C')
pylab.minorticks_on()
pylab.errorbar(t1, h1, dh1, linestyle='-', marker='.', color= 'blue')
pylab.errorbar(t2, h2, dh2, linestyle='-', marker='.', color= 'green')
pylab.errorbar(t3, h3, dh3, linestyle='-', marker='.', color= 'red')
pylab.errorbar(t4, h4, dh4, linestyle='-', marker='.', color= 'violet')
pylab.ylim(0,1.5)

pylab.subplot(222)
pylab.rc('font',size=18)
pylab.xlabel('T')
pylab.ylabel('M')
pylab.minorticks_on()
pylab.errorbar(t1, m1, dm1, linestyle='-', marker='.', color= 'blue')
pylab.errorbar(t2, m2, dm2, linestyle='-', marker='.', color= 'green')
pylab.errorbar(t3, m3, dm3, linestyle='-', marker='.', color= 'red')
pylab.errorbar(t4, m4, dm4, linestyle='-', marker='.', color= 'violet')
pylab.ylim(0,1)


pylab.subplot(223)
pylab.rc('font',size=18)
pylab.xlabel('T')
pylab.ylabel('E')
pylab.minorticks_on()
pylab.errorbar(t1, e1, de1, linestyle='-', marker='.', color= 'blue')
pylab.errorbar(t2, e2, de2, linestyle='-', marker='.', color= 'green')
pylab.errorbar(t3, e3, de3, linestyle='-', marker='.', color= 'red')
pylab.errorbar(t4, e4, de4, linestyle='-', marker='.', color= 'violet')
pylab.ylim(1.1,2.1)


pylab.subplot(224)
pylab.rc('font',size=18)
pylab.xlabel('T')
pylab.ylabel('X')
pylab.minorticks_on()
pylab.errorbar(t1, x1, dx1, linestyle='-', marker='.', color= 'blue')
pylab.errorbar(t2, x2, dx2, linestyle='-', marker='.', color= 'green')
pylab.errorbar(t3, x3, dx3, linestyle='-', marker='.', color= 'red')
pylab.errorbar(t4, x4, dx4, linestyle='-', marker='.', color= 'violet')
pylab.ylim(0,40)


# show the plot
pylab.show()

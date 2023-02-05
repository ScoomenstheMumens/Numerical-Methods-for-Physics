import matplotlib
import matplotlib.pyplot as pylab
import numpy
from scipy.optimize import curve_fit

# data load
t,c1,dc1,c2,dc2,c3,dc3,c4,dc4,c5,dc5,m,dm,e,de,x,dx,h,dh=numpy.loadtxt('2wolfcri_2.txt',unpack=True)
L=100


def f1(x,b,f,c):
    return b*pow((f-x),c)
L=50
x1=t
y1=m
xx1=numpy.linspace(0.1,1.15,20000)
init1=(1,1.20,0.12)
# set the error
sigma1=dm+0.0001
w1=1/sigma1**2

def f2(x,b,tc,c,d):
    return b*pow((tc-x),c)+d
x2=t
y2=h
xx2=numpy.linspace(0.1,1.15,20000)
init2=(0.1,1.20,-1.75,1.0)
# set the error
sigma2=dh+0.0001
w2=1/sigma2**2
# call the minimization routine
pars1,covm1=curve_fit(f1,x1,y1,init1,sigma1, absolute_sigma=False)
pars2,covm2=curve_fit(f2,x2,y2,init2,sigma2, absolute_sigma=False)

# calculate the chisquare for the best-fit function
chi21 = ((w1*(y1-f1(x1,*pars1))**2)).sum()
chi22 = ((w2*(y2-f2(x2,*pars2))**2)).sum()

# determine the ndof
ndof1=len(x1)-len(init1)
ndof2=len(x2)-len(init2)
# print results on the console
print('pars1:',pars1)
print('covm1:',covm1)
print ('chi21, ndof1:',chi21, ndof1)

print('pars2:',pars2)
print('covm2:',covm2)
print ('chi22, ndof2:',chi22, ndof2)

pylab.subplot(211)
pylab.rc('font',size=18)
pylab.xlabel('x')
pylab.ylabel('y')
pylab.minorticks_on()
pylab.errorbar(x1, y1, sigma1, linestyle='', marker='+', color= 'black')
pylab.plot(xx1,f1(xx1,*pars1), color='red')
pylab.grid(2, ls='--', alpha=0.5, lw=1)

pylab.subplot(212)
pylab.rc('font',size=18)
pylab.xlabel('x1')
pylab.ylabel('y1')
pylab.minorticks_on()
pylab.errorbar(x2, y2, sigma2, linestyle='', marker='+', color= 'black')
pylab.plot(xx2,f2(xx2,*pars2), color='blue')
pylab.grid(2, ls='--', alpha=0.5, lw=1)
pylab.ylim(0,100)
# show the plot
pylab.show()


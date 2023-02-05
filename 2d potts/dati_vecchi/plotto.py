import matplotlib
import matplotlib.pyplot as plt
import numpy as np

m,e=np.loadtxt('2wolfprova_2.txt',unpack=True)

def autocorr(x,lags):
    '''manualy compute, non partial'''

    mean=np.mean(x)
    var=np.var(x)
    xp=x-mean
    corr=[1. if l==0 else np.sum(xp[l:]*xp[:-l])/len(x)/var for l in lags]

    return np.array(corr)

t=range(len(m))
corr=autocorr(m,t)
dc=0.01

plt.subplot(2,1,1)
plt.rc('font',size=18)
plt.xlabel('T')
plt.ylabel('correlazione M')
plt.minorticks_on()
plt.errorbar(t, corr, dc, linestyle='-', marker='.', color= 'blue')
plt.ylim(-0.2,1.1)


plt.subplot(2,1,2)
plt.rc('font',size=18)
plt.xlabel('T')
plt.ylabel('correlazione M')
plt.minorticks_on()
plt.errorbar(t, m, dc, linestyle='-', marker='.', color= 'blue')
plt.ylim(0,1.1)
    
plt.show()

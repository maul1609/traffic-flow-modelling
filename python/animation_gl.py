"""
    this script animates the results from the traffic model
"""
import numpy as np
import matplotlib.pyplot as plt

(r,c)=np.shape(s.y)
num1=int(c/2)
ypos=s.y[:,0:num1]
ymax=np.max(ypos[:])

plt.figure()
ax=plt.subplot(121)

for i in range(len(s.t)):
    
    ax.clear()

    ind,=np.where((s.y[i,0:num1]>0) & (s.vehicle[i,0:num1]==1))
    if(len(ind)):
        p1=plt.plot(np.zeros(len(ind)), s.y[i,ind],'.',markersize=6,color='red')

    ind,=np.where((s.y[i,0:num1]>0) & (s.vehicle[i,0:num1]==2))
    if(len(ind)):
        p1=plt.plot(np.zeros(len(ind)), s.y[i,ind],'.',markersize=20,color='green')

    ind,=np.where((s.y[i,0:num1]>0) & (s.vehicle[i,0:num1]==3))
    if(len(ind)):
        p1=plt.plot(np.zeros(len(ind)), s.y[i,ind],'.',markersize=6,color='green')
    
    
    plt.text(-8,0,'Start')
    plt.text(-8,450,'M-S')
    plt.text(-8,780,'UPS')
    plt.text(-8,1000,'Traffic light')
    plt.text(-8,1250,'Traffic light')
    
    
    if(s.lights[i,0]==0):
        plt.plot(-1,1000.,'gp',markersize=8)
    else:    
        plt.plot(-1,1000.,'rp',markersize=8)
    
    if(s.lights[i,1]==0):
        plt.plot(-1,1250,'gp',markersize=8)
    else:    
        plt.plot(-1,1250,'rp',markersize=8)
    
    plt.ylim((0,ymax))
    plt.xlim((-10,5))
    plt.title('Time: ' +str(s.t[i]) + ' s')
    plt.pause(0.25)
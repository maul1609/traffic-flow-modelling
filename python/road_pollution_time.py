import numpy as np
import matplotlib.pyplot as plt
import gaussian_plume_model as gp
import scipy.interpolate as scint
from tqdm import tqdm

def road_pollution_time01(t,y,num_pass, vehicle, timer, \
    xrec=[200,800,1000],yrec=[100,100,100],zlevel=0.,emission1=[10.,100,15.]):
    """
        this scrip just calculates the distribution of pollution using a 
        gaussian plume model
    """
    
    nrec=len(yrec)
    (r,c)=np.shape(y)
    num1=int(c/2)
    
    ypos=y[:,0:num1]
    ymax=np.max(ypos[:])
    
    # now loop over all vehicles and their positions to calculate the map of pollution
    (r,c)=np.shape(vehicle)
    pollution_rec=np.zeros((r,nrec))

    # call just to define xs    
    (pollution,xs,ys)=gp.gaussian_plume_model([0,500],10.,zlevel,ymax)

    
    pollution_store=np.zeros(np.shape(pollution))
    
    indxarr=np.mgrid[0:len(xs[0,:])]
    indyarr=np.mgrid[0:len(ys[:,0])]
    # loop over time-steps
    for i in tqdm(range(r)):
#         print('Time: ',str(t[i]))
        
        ncalcs=0
        
        pollution_store[:]=0.
        
        # loop over each vehicle
        for j in range(c):
            if vehicle[i,j]==0:
                continue
                
            # use gaussian plume model to calculate emission based on position
            posx=y[i,j]
            speedx=y[i,j+c] # m/s
            
            # g/s if it is emission1 g/mile
            emission=emission1[int(vehicle[i,j])-1]*speedx/1609.344 
            
            (pollution,xs,ys)= \
                gp.gaussian_plume_model([posx,500],emission,zlevel,ymax)
            
            pollution_store[:,:] += pollution # add up for every car
            ncalcs += 1
            
        # find the position of the receptors        
        f1=scint.interp1d(xs[0,:],indxarr,kind='nearest')
        f2=scint.interp1d(ys[:,0],indyarr,kind='nearest')
        for k in range(nrec):
            ipos = int(f1(xrec[k]))
            jpos = int(f2(yrec[k]))
            pollution_rec[i,k]=pollution_store[jpos,ipos]

    return (t,pollution_rec)
        
(t,pollution_rec)=road_pollution_time01(s.t,s.y,s.num_pass, s.vehicle, s.timer)    

plt.figure()
plt.ion()

plt.plot(t,pollution_rec)
plt.xlabel('time')
plt.ylabel('Pollution level (g m$^{-3}$)')
plt.yscale('log') 
plt.ylim((1e-8,1e-4))   
plt.legend(['R1','R2','R3'])   
plt.show()

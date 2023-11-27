import numpy as np
import matplotlib.pyplot as plt
import gaussian_plume_model as gp
import scipy.interpolate as scint
import matplotlib.colors as colors
from tqdm import tqdm

def road_pollution01(t,y,num_pass, vehicle, timer, zlevel=2.,emission1=[10.,100,15.]):
    """
        this scrip just calculates the distribution of pollution using a 
        gaussian plume model
    """
    
    (r,c)=np.shape(y)
    num1=int(c/2)
    
    ypos=y[:,0:num1]
    ymax=np.max(ypos[:])
    
    # now loop over all vehicles and their positions to calculate the map of pollution
    (r,c)=np.shape(vehicle)

    (pollution,xs,ys)= \
                gp.gaussian_plume_model([0.,500],10.,zlevel,ymax)
    pollution_store=np.zeros(np.shape(pollution))
    
    # loop over time-steps
    ncalcs=0
    for i in tqdm(range(r)):
        #print('Time: ',str(t[i]))
                
        # loop over each vehicle
        for j in range(c):
            if vehicle[i,j]==0:
                continue
                
            # use gaussian plume model to calculate emission based on position
            posx=y[i,j]
            speedx=y[i,j+c] # m/s
            
            # g/s if it is emission1 g/mile
            emission=emission1[int(vehicle[i,j])-1] #*speedx/1609.344 
            
            (pollution,xs,ys)= \
                gp.gaussian_plume_model([posx,500],emission,zlevel,ymax)
                            
            pollution_store[:,:] += pollution*(t[1]-t[0]) # add up for every car
            ncalcs += 1
            
    pollution_store /= (t[-1]-t[0])
    
    return (xs,ys,pollution_store)
        
(xs,ys,pollution_store)=road_pollution01(s.t,s.y,s.num_pass, s.vehicle, s.timer)    

plt.ion()
plt.figure()
plt.pcolor(xs,ys,pollution_store,norm=colors.LogNorm(vmin=pollution_store.min(), vmax=pollution_store.max()))
plt.xlabel('x-position along road (m)')
plt.ylabel('y-position perpendicular to road (m)')
plt.colorbar()
plt.show()
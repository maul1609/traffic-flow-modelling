import numpy as np
import matplotlib.pyplot as plt

def cars_decibel(x):
    """
        figure 6 of steven et al. 2005
        https://www.umweltbundesamt.de/sites/default/files/medien/publikation/long/3092.pdf
    """
#     return 60.*np.ones(np.shape(x))
    return 13.257*np.log(x)+20.624
    
def hgvs_decibel(x):
    """
        figure 13 of steven et al. 2005
        https://www.umweltbundesamt.de/sites/default/files/medien/publikation/long/3092.pdf
    """
#     return 75.*np.ones(np.shape(x))
#     return 0.1233*x+75.769
    return 0.2194*x+75.556
    

def noise_pollution01(t,y,vehicle,house_location1=450,house_location2=20):
    """
        This script simulates the noise levels from simulated vehicles
        travelling down a road
        
        in order to estimate the increase in noise that hgvs have on green lane
        assume that the sound always adds constructively
        
        Call as follows (note you must have run a simulation first)
        
        >> sound_level=noise_pollution01(t,y,vehicle)
        the sound_level variables contains noise level in dBs
    """
    sound_level=np.zeros((len(t)))
    cars_reference=7.5
    hgvs_reference=7.5
    pref=20e-6
    
    (xx,yy)=np.shape(y)
    num1=int(yy/2)
    # loop through all vehicles and calculate the decibel level at the house
    for i in range(len(t)):
        ind1,=np.where(vehicle[i,:]==1) # cars
        ind2,=np.where(vehicle[i,:]==2) # hgvs
        ind3,=np.where(vehicle[i,:]==3) # ups
        
        ind1a=ind1
        ind2a=ind2
        ind3a=ind3

        if(len(ind1a)):
            ind1a=ind1a+num1
        if(len(ind2a)):
            ind2a=ind2a+num1
        if(len(ind3a)):
            ind3a=ind3a+num1
        
        # distance of sound sources from house - Pythagoras' theorem - these
        # are root mean square values
        distance1 = np.sqrt((y[i,ind1]-house_location1)**2+house_location2**2)
        distance2 = np.sqrt((y[i,ind2]-house_location1)**2+house_location2**2)
        distance3 = np.sqrt((y[i,ind3]-house_location1)**2+house_location2**2)
        
        cars_db=np.maximum(cars_decibel(y[i,ind1a]*3600./1000.),60.)
        hgvs_db=np.maximum(hgvs_decibel(y[i,ind2a]*3600./1000.),75.)
        ups_db=np.maximum(cars_decibel(y[i,ind3a]*3600./1000.),60.)
        try:
            cars_db=np.nanmax(cars_db)
        except:
            cars_db=60.
        try:
            hgvs_db=np.nanmax(hgvs_db)
        except:
            hgvs_db=75.
        try:
            ups_db=np.nanmax(ups_db)
        except:
            ups_db=60.
            
        cars_pref=pref*10**(cars_db/20.)
        hgvs_pref=pref*10**(hgvs_db/20.)
        ups_pref=pref*10**(ups_db/20.)
        """
            now add up all cars and hgvs, using superposition principle
        """
        pressure_amplitudes=np.append(cars_pref*cars_reference/distance1, \
            hgvs_pref*hgvs_reference/distance2)
        pressure_amplitudes=np.append(pressure_amplitudes, \
            ups_pref*cars_reference/distance3)
        
        pressure_amplitudes = pressure_amplitudes*np.sqrt(2.)

        # randomise the phases of each wave arriving at house
        rand1=2.*np.pi*np.random.rand(1000,len(pressure_amplitudes))

        # find the average of the total pressure wave
        total_pressure=np.sqrt(np.mean( \
            np.sum(np.repeat(np.transpose(np.expand_dims(pressure_amplitudes,1)),[1000],axis=0) \
                *np.sin(rand1),axis=1 )**2)  )
            
        pressure_level=total_pressure
    
        sound_level[i]=20.*np.log10(pressure_level/pref)
        
    return sound_level
        
sound_level=noise_pollution01(s.t,s.y,s.vehicle) 
plt.figure()
plt.ion()

plt.plot(s.t,sound_level)
plt.xlabel('time')
plt.ylabel('Noise level (dB)')
plt.show()

 
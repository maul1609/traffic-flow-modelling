import numpy as np
import scipy.stats as scs
from tqdm import tqdm
from scipy.integrate import RK45, RK23, BDF, LSODA, odeint, ode, solve_ivp

class State1:
    def __init__(self,runtime=3600,dt=10,rands01=np.random.rand(10000,5)):
#         self.acc_fac=5./7.
        self.acc_fac_car=5./7.
        self.acc_fac_hgv=5./30.
        self.acc_fac=self.acc_fac_car
        self.hgv_turn_flag = 1
        self.ups_turn_flag = 1
        self.n_car=0
        self.rands01=rands01
        self.count_car=0
        self.count_hgv=0     
        self.count_ups=0     
        self.length=1250         # length of green lane in metres
        self.m_s=490             # distance from roundabout to mitchell-shackleton
        self.ups=780             # distance from roundabout to UPS
        self.cromwell=940        # distance from roundabout to cromwell
        self.trafficlight1=1000  # distance from roundabout to traffic light 1
        self.trafficlight2=1250  # distance from roundabout to traffic light 2
        self.traffic_light_switch=1 # to switch them on and off
    
        self.n_cars_rate=500*24  # number of cars per day
        self.n_ups_rate=50*24    # number of UPS vans per day
        self.n_hgv_rate=80*24*8  # number of HGVs going in and out of mitchell-shackleton each day
        self.n_tot=int(np.ceil((self.n_cars_rate+self.n_ups_rate+self.n_hgv_rate) / 24))
    
        self.acc_car=0.4 # acceleration of car
        self.dec_car=0.8 # deceleration of car
        self.acc_hgv=0.2 # acceleration of large truck
        self.dec_hgv=0.6 # deceleration of large truck
        self.acc_ups=0.4 # acceleration of ups
        self.dec_ups=0.8 # deceleration of ups
    
        self.top_speed=30*1609.344/3600.    # speed limit along road - 30 mph
    
        # rules
        # lights change every 1 minute, standard deviation of 5 seconds
        self.traffic_change_mean=75
        self.traffic_change_std=5
        self.traffic_inter1=scs.norm.ppf(self.rands01[:,0],self.traffic_change_mean, \
                        self.traffic_change_std)
        self.traffic_inter2=scs.norm.ppf(self.rands01[:,1],self.traffic_change_mean, \
                        self.traffic_change_std)
        self.traffic_ind1=0
        self.traffic_ind2=0
    
        # cars arrive so that mean interarrival time is n /day and 
        # standard deviation of 10%
        self.car_arrive_mean=86400/self.n_cars_rate
        self.car_arrive_std=self.car_arrive_mean*0.1
    
        # UPS vans arrive so that mean interarrival time is n /day and 
        # standard deviation of 40%
        self.ups_arrive_mean=86400/self.n_ups_rate
        self.ups_arrive_std=self.ups_arrive_mean*0.4
    
        # HGVs arrive so that mean interarrival time is n /day and 
        # standard deviation of 40%
        self.hgv_arrive_mean=86400/self.n_hgv_rate
        self.hgv_arrive_std=self.hgv_arrive_mean*0.4
        
        self.roundabout_waiting_time=10. # roundabout waiting time (s)
        
        # model stuff
        self.dt=dt
        self.t=np.mgrid[0:runtime+dt:dt]  
        self.num_pass=np.zeros(np.shape(self.t))
        self.y=np.zeros((len(self.t),self.n_tot*2))
        self.timer=np.zeros((self.n_tot,2))
        self.j_start=0
        self.j_stop=0
        self.tt=0.
        self.tt1=0.
        self.thgv=0.
        self.tups=0.
        self.tr1=0.;self.red_light1=0.;
        self.tr2=0.;self.red_light2=0.;
        self.lights=np.zeros((len(self.t),2))
        self.vehicle=np.zeros((len(self.t),self.n_tot))
        self.i=0
        self.flag=np.zeros((self.n_tot,1))
        
    
    def alter_HGVs(self, fracHGV,fracUPS):
        # numebr of cars in one hour
        num_car1=3600./self.car_arrive_mean
        # if supplied, argument fracHGV is the fraction of traffic that is HGVs
        # and num_car1 is the total number of vehicles
        if((fracHGV+fracUPS)==1.):
            self.car_arrive_mean=1000000.
        else:
            self.car_arrive_mean=3600./((1.-fracHGV-fracUPS) * num_car1)  

        if(fracHGV>0.): 
            self.hgv_arrive_mean=3600./((fracHGV) * num_car1)     
        else:
            self.hgv_arrive_mean=1000000.

        if(fracUPS>0.): 
            self.ups_arrive_mean=3600./((fracUPS) * num_car1)     
        else:
            self.ups_arrive_mean=1000000.

        self.car_arrive_std=self.car_arrive_mean*0.1 
        self.ups_arrive_std=self.ups_arrive_mean*0.1 
        self.hgv_arrive_std=self.hgv_arrive_mean*0.1
        
        
    def update(self):
        self.waiting_time=scs.norm.ppf(self.rands01[0,2],self.car_arrive_mean,\
                self.car_arrive_std)
        self.count_car=self.count_car+1
        
        self.waiting_time_hgv=scs.norm.ppf(self.rands01[0,3],self.hgv_arrive_mean, \
                self.hgv_arrive_std)
        self.count_hgv=self.count_hgv+1
        
        self.waiting_time_ups=scs.norm.ppf(self.rands01[0,4],self.ups_arrive_mean, \
                self.ups_arrive_std)
        self.count_ups=self.count_ups+1
        
        self.traffic_time1=self.traffic_inter1[self.traffic_ind1]
        self.traffic_ind1 += 1
        self.traffic_time2=self.traffic_inter2[self.traffic_ind2]
        self.traffic_ind2 += 1
        
def go1(t, y):
    """
        calculates the velocities and accelerations
    """
    ydot = np.zeros(((s.n_tot)*2))
    t1=t-s.tt1

    if(s.traffic_light_switch==1):
        # set traffic light 1
        if((t1+s.tt1-s.tr1)>=s.traffic_time1):
            s.traffic_time1=s.traffic_inter1[s.traffic_ind1];s.traffic_ind1 += 1
            s.tr1=s.tt1+t1
            if(s.red_light1==1):
                s.red_light1=0
            else:
                s.red_light1=1
            
    
        # set traffic light 2
        if((t1+s.tt1-s.tr2)>=s.traffic_time2):
            s.traffic_time2=s.traffic_inter2[s.traffic_ind2];s.traffic_ind2 += 1
            s.tr2=s.tt1+t1
            if(s.red_light2==1):
                s.red_light2=0
            else:
                s.red_light2=1

    # decide what vehicle type is coming from the roundabout
    # time it starts journey
    if((t1+s.tt1-s.tt) >= s.waiting_time):
        s.vehicle[s.i,s.n_car]=1
        s.n_car += 1
        if(s.n_car>s.n_cars_rate):
            print('n_cars > n_cars_rate')
        s.waiting_time=scs.norm.ppf(s.rands01[s.count_car,2], \
            s.car_arrive_mean, s.car_arrive_std)
        s.count_car += 1
        s.tt=s.tt1+t1
        s.timer[s.j_start,0]=s.tt
        s.j_start += 1
        
    # time it reaches end of road
    for i in range(s.n_car):
        if((y[i]>=s.length) & (s.flag[i,0]==0) & (s.vehicle[s.i,i]==1)):
            # during the outer loop y[i] will be set back to 0, and be available
            s.timer[s.j_stop,1]=t1+s.tt1
            s.j_stop += 1
            s.flag[i,0]=1
            
    if((t1+s.tt1-s.thgv)>=s.waiting_time_hgv):
        s.vehicle[s.i,s.n_car]=2
        s.n_car += 1
        if(s.n_car>s.n_cars_rate):
            print('n_car > n_cars_rate')
        s.waiting_time_hgv=scs.norm.ppf(s.rands01[s.count_hgv,3], \
            s.hgv_arrive_mean, s.hgv_arrive_std)
        s.count_hgv += 1
        s.thgv=s.tt1+t1
        
    if((t1+s.tt1-s.tups)>=s.waiting_time_ups):
        s.vehicle[s.i,s.n_car]=3
        s.n_car += 1
        if(s.n_car>s.n_cars_rate):
            print('n_car > n_cars_rate')
        s.waiting_time_ups=scs.norm.ppf(s.rands01[s.count_ups,4], \
            s.ups_arrive_mean, s.ups_arrive_std)
        s.count_ups += 1
        s.tups=s.tt1+t1
        
    # initial definition of accelerations
    for i in range(s.n_car):
        ydot[i+s.n_tot]=acceleration(y[i+s.n_tot],s.top_speed, s.vehicle[s.i,i])
        # velocity
        ydot[i] = y[i+s.n_tot]
    
    # now adjust acceleration so that cars do not collide
    for i in range(0,s.n_car-1):
        if((s.flag[i,0] == 1) |(s.vehicle[s.i,i]==0) ):  # if at end of road
            continue
        if(s.vehicle[s.i,i]==1):
            # if the stopping distance is greater than 5 m
            if((y[i+1]-5.-y[i]) < (y[i+s.n_tot]/(s.acc_fac*s.dec_car))):
                ydot[i+s.n_tot]=acceleration(y[i+s.n_tot],y[i+1+s.n_tot],s.vehicle[s.i,i])
                ydot[i]=y[i+s.n_tot]
        elif(s.vehicle[s.i,i]==2):
            # if the stopping distance is greater than 5 m
            if((y[i+1]-5.-y[i]) < (y[i+s.n_tot]/(s.acc_fac*s.dec_hgv))):
                ydot[i+s.n_tot]=acceleration(y[i+s.n_tot],y[i+1+s.n_tot],s.vehicle[s.i,i])
                ydot[i]=y[i+s.n_tot]
        elif(s.vehicle[s.i,i]==3):
            # if the stopping distance is greater than 5 m
            if((y[i+1]-5.-y[i]) < (y[i+s.n_tot]/(s.acc_fac*s.dec_ups))):
                ydot[i+s.n_tot]=acceleration(y[i+s.n_tot],y[i+1+s.n_tot],s.vehicle[s.i,i])
                ydot[i]=y[i+s.n_tot]
          
    # slow HGVs so they can turn  
    if(s.hgv_turn_flag==1):
        for i in range(s.n_car):
            if(s.vehicle[s.i,i]==2):
                if((s.m_s-y[i]) < (y[i+s.n_tot]/(s.acc_fac_hgv*s.dec_hgv))):
                    s.acc_fac=s.acc_fac_hgv
                    ydot[i+s.n_tot]=acceleration(y[i+s.n_tot],0.,s.vehicle[s.i,i])
                    ydot[i]=y[i+s.n_tot]
                    s.acc_fac=s.acc_fac_car
                    
    # slow UPS so they can turn  
    if(s.ups_turn_flag==1):
        for i in range(s.n_car):
            if(s.vehicle[s.i,i]==3):
                if((s.ups-y[i]) < (y[i+s.n_tot]/(s.acc_fac*s.dec_ups))):
                    ydot[i+s.n_tot]=acceleration(y[i+s.n_tot],0.,s.vehicle[s.i,i])
                    ydot[i]=y[i+s.n_tot]
                    
    # slow vehicles for traffic light 1
    if(s.red_light1==1):
        for i in range(s.n_car):
            if(s.vehicle[s.i,i]==0):
                continue
            elif(s.vehicle[s.i,i]==1):
                dec1=s.dec_car
            elif(s.vehicle[s.i,i]==2):
                dec1=s.dec_hgv
            elif(s.vehicle[s.i,i]==3):
                dec1=s.dec_ups
                
            if(((s.trafficlight1-y[i]) < (y[i+s.n_tot]/(s.acc_fac*dec1)) ) & \
                ((s.trafficlight1-y[i])>-1.)):
                ydot[i+s.n_tot]=acceleration(y[i+s.n_tot],0.,s.vehicle[s.i,i])
                ydot[i]=y[i+s.n_tot]
                    
    # slow vehicles for traffic light 2
    if(s.red_light2==1):
        for i in range(s.n_car):
            if(s.vehicle[s.i,i]==0):
                continue
            elif(s.vehicle[s.i,i]==1):
                dec1=s.dec_car
            elif(s.vehicle[s.i,i]==2):
                dec1=s.dec_hgv
            elif(s.vehicle[s.i,i]==3):
                dec1=s.dec_ups

            if(((s.trafficlight2-y[i]) < (y[i+s.n_tot]/(s.acc_fac*dec1)) ) & \
                ((s.trafficlight2-y[i])>-1.)):
                ydot[i+s.n_tot]=acceleration(y[i+s.n_tot],0.,s.vehicle[s.i,i])
                ydot[i]=y[i+s.n_tot]
                
    
    return ydot

def acceleration(u,v,vehicle):
    """
        acceleration and deceleration depending on vehicle type
    """
    if(vehicle==1): # car
        peak_acc = s.acc_car if (v>u) else -s.dec_car
    elif(vehicle==2): # hgv
        peak_acc = s.acc_hgv if (v>u) else -s.dec_hgv
    elif(vehicle==3): # ups
        peak_acc = s.acc_ups if (v>u) else -s.dec_ups
    
    # scale acceleration so that peak is when velocity is zero    
    return np.abs(v-u)*peak_acc*s.acc_fac
    
def greenlanetraffic01(rands01):
    """
        This Python script simulates the flow of traffic down a road with the 
        following assumptions:
        
        1. Traffic moves in one direction down the road, starting at a roundabout
        2. The vehicles accelerate coming off the roundabout to 30 mph
        3. Vehicles can turn into two industrial areas (Mitchell-Shackleton and 
                UPS - United Parcel Service)
        4. Cars travel to the end of the road - there are two sets of traffic
            lights positioned near the end of the road. 
            
        The output is as follows:
        t:          time of the simulation output
        y:          solution, which is an array that contains the positions and
                        velocities of each vehicle
        num_pass:   number of active vehicles that have passed end of green lane
                    this is just a counter, the total can be plotted as follows
                    >> plt.plot(t, np.cumsum(num_passs))
                    
        lights:     flag for lights being on or off (they can only be red or green)
                    it is an array that is [len(time),2]. Plot like this:
                    
                    >> plt.plot(t, lights[:,0],'k')
                    >> plt.plot(t, lights[:,1],'b')
                    
        vehicle:    this array stores the type of vehicle:
                    1 = car
                    2 = HGV
                    the array is bigger than it needs to be, depending on how many
                    vehicles are on the road at any one time
                    
        timer:      When a vehicle starts going down green lane a timer is started.
                    this timer is stopped when it finishes
                    you can plot the timers by doing:
                    
                    >> ntot=np.cumsum(num_pass);
                    >> ntot=ntot[-1]
                    >> plt.hist(timer[0:ntot,1]-timer[0:ntot,0])  
                    
        Running the model: base-line case, type:
        
        >> rands01=np.rand((10000,4)) # only needs to be typed once
        >> (t,y,num_pass,lights, vehicle, timer) = greenlanetraffic01(rands01)
        
        Note, if you call with a second argument, it is the hgv_turn flag
        set to 0 or 1. When set to 1, the HGVs turn into mitchell-shackleton
        when set to 0 they do not, e.g. 
        
        >> (t,y,num_pass,lights, vehicle, timer) = greenlanetraffic01(rands01,0)
        
        note, if you call with a third argument we assume that the cars_rate is 
        the total amount of traffic. The third argument is the fraction of 
        vehicles that are HGVs, e.g. (the following is 50% HGV traffic)
        
        >> (t,y,num_pass,lights, vehicle, timer) = greenlanetraffic01(rands01,1,0.5)
    """
    

    
    # start model simulation here
    for i in tqdm(range(len(s.t)-1)):
        s.i=i+1
        s.vehicle[s.i,:]=s.vehicle[s.i-1,:]
        s.flag[:,:]=0
        
        # call the solver over one time-step (10s)
        sol=solve_ivp(go1,[s.t[i],s.t[i+1]],np.squeeze(s.y[s.i-1,:]),method='RK45',\
            atol=0.001,rtol=1.e-8,max_step=10., \
            vectorized=False,t_eval=[s.t[i+1]],dense_output=False)
        
        # store the results and process output
        s.y[s.i,:]=sol.y[:,-1]
        s.tt1 = sol.t[-1]
    
        ind,=np.where((s.flag[0:s.n_car,0]==1))
        s.num_pass[s.i]=len(ind)
        # hgvs - if gone past m-s set to the end of road
        if(s.hgv_turn_flag == 1):
            ind,=np.where((s.y[s.i,0:s.n_tot]>=s.m_s) & (s.vehicle[s.i,:]==2) )
            s.y[s.i,ind]=s.length*2

        # ups vehicles        
        if(s.ups_turn_flag == 1):
            ind,=np.where((s.y[s.i,0:s.n_tot]>=s.ups) & (s.vehicle[s.i,:]==3) )
            s.y[s.i,ind]=s.length*2
    
    
    
    
        # sort the positions of the cars along the road - they shouldn't overtake, but
        # may do due to inaccuracies in solver 
#         ind,=np.where((s.y[s.i,0:s.n_car]<s.length))
        ind,=np.where((s.flag[0:s.n_car,0]==0) & (s.vehicle[s.i,0:s.n_car]!=0))
        ind2=np.argsort(s.y[s.i,ind]);ind=ind[ind2]
        s.n_car=len(ind)
        if(len(ind)>0):
            s.y[s.i,0:len(ind)]=s.y[s.i,ind]
            s.vehicle[s.i,0:len(ind)]=s.vehicle[s.i,ind]
            s.y[s.i,(s.n_tot):(len(ind)+s.n_tot)] = s.y[s.i,ind+s.n_tot]
            s.n_car=len(ind)
            
            s.y[s.i,(len(ind)):s.n_tot]=0.
            s.y[s.i,(len(ind)+s.n_tot):2*s.n_tot]=0.
            s.vehicle[s.i,(len(ind)):s.n_tot]=0.
        
        
        # this is to keep everything that has length < street length
        ind,=np.where((s.y[s.i,0:s.n_car]<s.length))
        ind2=np.argsort(s.y[s.i,ind]);ind=ind[ind2]
        s.n_car=len(ind)
        if(len(ind)>0):
            s.y[s.i,0:len(ind)]=s.y[s.i,ind]
            s.vehicle[s.i,0:len(ind)]=s.vehicle[s.i,ind]
            s.y[s.i,(s.n_tot):(len(ind)+s.n_tot)] = s.y[s.i,ind+s.n_tot]
            s.n_car=len(ind)
            
            s.y[s.i,(len(ind)):s.n_tot]=0.
            s.y[s.i,(len(ind)+s.n_tot):2*s.n_tot]=0.
            s.vehicle[s.i,(len(ind)):s.n_tot]=0.
        
        # store traffic light state    
        s.lights[s.i-1,0]=s.red_light1
        s.lights[s.i-1,1]=s.red_light2
    # end of main loop
    
    
defRand=False
if(defRand):
    rands01=np.random.rand(10000,5)   
    
# instance a state object
s = State1(runtime=3600,rands01=rands01)
s.alter_HGVs(0.,0.)
s.update()  
greenlanetraffic01(s.rands01)

#plt.hist(s.timer[0:s.j_stop,1]-s.timer[0:s.j_stop,0],alpha=0.5,density=True)
#plt.plot(s.timer[0:s.j_stop,1]-s.timer[0:s.j_stop,0])
# plt.figure()
# ind1,ind2,=np.where((s.vehicle !=0) & (s.y[:,0:s.n_tot]>0) & (s.y[:,0:s.n_tot]<s.length)) 
# ys=s.y[ind1,ind2]       
# (n,x)=np.histogram(ys[:])
# n2=np.append(np.insert(n,0,0),0) 
# x2=0.5*(x[0:-1]+x[1:])
# x2=np.insert(x2,0,x2[0]-(x2[1]-x2[0])) 
# x2=np.append(x2,x2[-1]+(x2[1]-x2[0])) 
# plt.step(x2,n2/len(s.t))

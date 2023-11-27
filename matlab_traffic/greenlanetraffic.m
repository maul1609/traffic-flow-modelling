function [t,y,num_pass,lights,vehicle,timer]=...
    greenlanetraffic(rands01,varargin)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% This MATLAB script simulates the flow of traffic down a road with the
% following assumptions:
%
% 1. Traffic moves in one direction down the road, starting at a roundabout
% 2. The vehicles accelerate coming off the roundabout to 30 mph 
% 3. Vehicles can turn into two industrial areas (Mitchell-Shackleton and
%       UPS - United Parcel Service)
% 4. Cars travel to the end of the road - there are two sets of traffic
%       lights positioned near the end of the road.
% 
% The output is as follows:
% t:        time of the simulation output
% y:        solution, which is an array that contains the positions and 
%               velocities of each vehicle
% num_pass: number of active vehicles that have passed end of green lane 
%           this is just a counter, the total can be plotted as follows:
%
%           >> plot(t,cumsum(num_pass))
%
% lights:   flag for lights being on or off (they can only be red or green)
%           it is an array that is time x 2. Plot like this:
%
%           >> plot(t,lights(:,1),'k'); hold on;
%           >> plot(t,lights(:,2),'b')
%
% vehicle:  this array stores the type of vehicle:
%           1 = car
%           2 = HGV
%           the array is bigger than it needs to be, depending on how many
%           vehicles are on the road at any one time
%
% timer:    When a vehicle starts going down green lane a timer is started.
%           this timer is stopped when it finishes
%           you can plot the timers by doing:
%           
%           >> ntot=cumsum(num_pass);
%           >> ntot=ntot(end);
%           >> hist(timer(1:ntot,2)-timer(1:ntot,1))
%
% Running the model: base-line case, type:
%
% >>rands01=rand(10000,4); % only needs to be typed once
% >>[t,y,num_pass,lights,vehicle,timer]=greenlanetraffic(rands01);
% 
% note, if you call with a second argument, it is the hgv_turn flag
% set to 0 or 1. When set to 1, the HGVs turn into mitchell-shackleton
% when set to 0 they do not, e.g.
%
% >>[t,y,num_pass,lights,vehicle,timer]=greenlanetraffic(rands01,0);
%
% note, if you call with a third argument we assume that the cars_rate is
% the total amount of traffic. The third argument is the fraction of
% vehicles that are HGVs, e.g (the following is 50% HGV traffic)
%
% >>[t,y,num_pass,lights,vehicle,timer]=greenlanetraffic(rands01,1,0.5);
%--------------------------------------------------------------------------

global s;

s.hgv_turn_flag=1;
if(nargin > 1)
    s.hgv_turn_flag=varargin{1};
end
s.n_car=0;
s.rands01=rands01;
s.count_car=1;
s.count_hgv=1;


% parameters for the model+++++++++++++++++++++++++++++++++++++++++++++++++
t=0:10:3600;            % run for one hour, outputting every 10s
s.length=1250;          % length of green lane in metres.
s.m_s=490;              % distance from roundabout to mitchell-shackleton.
s.ups=780;              % distance from roundabout to UPS 
s.cromwell=940;         % distance from roundabout to cromwell.
s.trafficlight1=1000;   % distance from roundabout to traffic light 1.
s.trafficlight2=1250;   % distance from roundabout to traffic light 2.

s.n_cars_rate=200*24.;      % number of cars per day
s.n_ups_rate=50*24.;        % number of ups vans per day
s.n_hgv_rate=80*24.*8;        % number of hgvs going in and out of 
                        %   mitchell-shackleton each day
s.n_tot=(s.n_cars_rate+s.n_ups_rate+s.n_hgv_rate)/24;

% How do cars drive?
% accelerate up to 30 mph http://physics.info/acceleration/
s.acc_car=0.4;          % acceleration of car
s.dec_car=0.8;          % deceleration of car
s.acc_hgv=0.2;          % acceleration of large truck
s.dec_hgv=0.6;          % deceleration of large truck

s.top_speed=30.*1609.344./3600; 
                        % speed limit along road  - 30 mph.

% rules 
% lights change every 1 minute, standard deviation of 5 seconds 
s.traffic_change_mean=60;
s.traffic_change_std=5;
s.traffic_inter1=norminv(rands01(:,1),...
    s.traffic_change_mean,s.traffic_change_std);
s.traffic_inter2=norminv(rands01(:,2),...
    s.traffic_change_mean,s.traffic_change_std);
s.traffic_ind1=1;
s.traffic_ind2=2;

% cars arrive so that mean interarrival time is 5000/day and 
%   standard deviation of 10% 
s.car_arrive_mean=86400./(s.n_cars_rate); 
s.car_arrive_std=s.car_arrive_mean*0.1;

% ups vans arrive so that mean interarrival time is 50/day and 
%   standard deviation of 40% 
s.ups_arrive_mean=86400./(s.n_ups_rate); 
s.ups_arrive_std=s.ups_arrive_mean*0.4;

% hgvs arrive so that mean interarrival time is 80/day 
%   (only in 8 hour period and standard deviation of 40%) 
s.hgv_arrive_mean=86400./(s.n_hgv_rate); 
if(nargin==3)
    % number of cars in one hour:
    num_car1=3600./s.car_arrive_mean;
    % if supplied, varargin{2} is the fraction of traffic that is HGVs
    % and num_car1 is the total number of vehicles
    s.car_arrive_mean=3600./((1-varargin{2}).*num_car1);
    s.hgv_arrive_mean=3600./(varargin{2}.*num_car1);
end
s.hgv_arrive_std=s.hgv_arrive_mean*0.4;


s.roundabout_waiting_time=10;   % roundabout waiting time (s)
%--------------------------------------------------------------------------




% do not edit++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
s.waiting_time=norminv(s.rands01(1,3),s.car_arrive_mean,s.car_arrive_std);
s.count_car=s.count_car+1;

s.waiting_time_hgv=norminv(s.rands01(1,4),s.hgv_arrive_mean,...
    s.hgv_arrive_std);
s.count_hgv=s.count_hgv+1;
s.traffic_time1=s.traffic_inter1(s.traffic_ind1);
s.traffic_ind1=s.traffic_ind1+1;
s.traffic_time2=s.traffic_inter2(s.traffic_ind2);
s.traffic_ind2=s.traffic_ind2+1;
%--------------------------------------------------------------------------












% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% START OF THE MODEL ++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
num_pass=zeros(size(t));
y(:,:)=zeros(length(t),(s.n_tot).*2);
s.timer=zeros((s.n_tot),2);
s.j_start=0;s.j_stop=0;
s.t=0;s.t1=0.;
s.thgv=0;
s.tr1=0.;s.red_light1=0;
s.tr2=0.;s.red_light2=0;
lights=zeros(length(t),2);
s.vehicle=zeros(length(t),(s.n_tot));

h = waitbar(0,'Please wait...');
for i=1:length(t)-1     % LOOP OVER TIME
    s.i=i+1;
    s.vehicle(i+1,:)=s.vehicle(i,:);
    s.flag=zeros(s.n_tot,1);

    
    
    % SET-UP SOLVER FOR THIS TIME-STEP
    ode1=odeset('AbsTol',...
        0.01.*ones(size(1:(s.n_tot).*2)));

    % CALL THE SOLVER OVER ONE TIME-STEP (10s):
    [t1,y1]=ode45(@go1,[0 10],y(i,:),ode1); 
    
   
    % STORE RESULTS and process output
    y(i+1,:)=y1(end,:);
    s.t1=s.t1+10;
    % test
    num_pass(i+1)=length(find(y(i+1,1:s.n_tot)>=s.length));
    % hgvs - if gone past m/s set to end.
    if(s.hgv_turn_flag==1)
        ind=find(y(i+1,1:s.n_tot)>=s.m_s &s.vehicle(i+1,:)==2);
        y(i+1,ind)=s.length.*2;
    end    
    ind=find(y(i+1,1:s.n_car)<=s.length);% & y(i+1,1:s.n_tot) >0 );
    [s1,ind2]=sort(y(i+1,ind));ind=ind(ind2);
    if(length(ind))
        %disp(['Here ',num2str(i)]);
        y(i+1,1:length(ind))=y(i+1,ind);
        s.vehicle(i+1,1:length(ind))=s.vehicle(i+1,ind);
        y(i+1,[1:length(ind)]+s.n_tot)=y(i+1,s.n_tot+ind);
        s.n_car=length(ind);
        y(i+1,length(ind)+1:s.n_tot)=0.;  y(i+1,[length(ind)+1:s.n_tot]+s.n_tot)=0.;
        s.vehicle(i+1,length(ind)+1:s.n_tot)=0.;
    end
    
    lights(i,1)=s.red_light1;
    lights(i,2)=s.red_light2;
    %[y2(i+1,1:s.n_tot),ii]=sort(y(i+1,1:s.n_tot));
    %y2(i+1,s.n_tot+1:end)=y(i+1,ii+s.n_tot);  
    waitbar(i./(length(t)-1),h);
end
close(h);
vehicle=s.vehicle;
timer=s.timer;
% -------------------------------------------------------------------------
% END OF THE MODEL --------------------------------------------------------
% -------------------------------------------------------------------------



%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%CALCULATES THE SPEEDS AND ACCELERATIONS ++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function ydot = go1(t,y)
global s;
ydot=zeros((s.n_tot).*2,1);


if(t+s.t1-s.tr1>=s.traffic_time1)
    s.traffic_time1=s.traffic_inter1(s.traffic_ind1);s.traffic_ind1=s.traffic_ind1+1;
    s.tr1=s.t1+t;
    if(s.red_light1==1)
        s.red_light1=0;
    else
        s.red_light1=1;
    end
end
if(t+s.t1-s.tr2>=s.traffic_time2)
    s.traffic_time2=s.traffic_inter2(s.traffic_ind2);s.traffic_ind2=s.traffic_ind2+1;
    %s.traffic_time2=norminv(rand(1,1),s.traffic_change_mean,s.traffic_change_std);
    s.tr2=s.t1+t;
    if(s.red_light2==1)
        s.red_light2=0;
    else
        s.red_light2=1;
    end
end
% decide what vehicle is coming from the roundabout
if(t+s.t1-s.t>=s.waiting_time)
    s.n_car=s.n_car+1;
    s.vehicle(s.i,s.n_car)=1;
    if(s.n_car>s.n_cars_rate)
        error('n_cars > n_cars_rate');
    end
    s.waiting_time=norminv(s.rands01(s.count_car,3),s.car_arrive_mean,s.car_arrive_std);
    s.count_car=s.count_car+1;
    s.t=s.t1+t;
    s.j_start=s.j_start+1;
    s.timer(s.j_start,1)=t+s.t1;
end
for i=1:s.n_car
    if(y(i)>=s.length & s.flag(i)==0)
        s.j_stop=s.j_stop+1;
        s.timer(s.j_stop,2)=t+s.t1;
        s.flag(i)=1;
    end    
end

if(t+s.t1-s.thgv>=s.waiting_time_hgv)
    s.n_car=s.n_car+1;
    s.vehicle(s.i,s.n_car)=2;
    if(s.n_car>s.n_cars_rate)
        error('n_car > n_cars_rate');
    end
    s.waiting_time_hgv=norminv(s.rands01(s.count_hgv,4),s.hgv_arrive_mean,s.hgv_arrive_std);
    s.count_hgv=s.count_hgv+1;
    s.thgv=s.t1+t;
end


for i=1:s.n_car
    ydot(i+s.n_tot)=acceleration(y(i+s.n_tot),s.top_speed,s.vehicle(s.i,i));
    ydot(i)=y(i+s.n_tot);
end

% now adjust acceleration so that cars do not collide
for i=2:s.n_car
    if (s.flag(i) == 1) continue; end
    % if the stopping distance is greater than 5 m
    if((y(i)-5-y(i-1)) < y(i-1+s.n_tot)/(5.*s.dec_car./30)  )
%     if(y(i)-5-y(i-1) < y(i-1+s.n_tot)^2/(2.*s.dec_hgv) )
        ydot(i+s.n_tot)=acceleration(y(i+s.n_tot),0,s.vehicle(s.i,i));
        ydot(i)=y(i+s.n_tot);
    end
end

if(s.hgv_turn_flag==1)
    for i=1:s.n_car
        if(s.vehicle(s.i,i)==2)
            if(s.m_s-y(i) < y(i+s.n_tot)/(5.*s.dec_hgv./30) )
                ydot(i+s.n_tot)=acceleration(y(i+s.n_tot),0,s.vehicle(s.i,i));
                ydot(i)=y(i+s.n_tot);
            end
        end
    end    
end
if(s.red_light1)
    for i=1:s.n_car
        if(s.trafficlight1-y(i) < y(i+s.n_tot)/(5.*s.dec_car./30) &&...
                (s.trafficlight1-y(i)>0))
            ydot(i+s.n_tot)=acceleration(y(i+s.n_tot),0,s.vehicle(s.i,i));
            ydot(i)=y(i+s.n_tot);
        end
    end    
end
if(s.red_light2 )
    for i=1:s.n_car
        if((s.trafficlight2-y(i)) < y(i+s.n_tot)/(5.*s.dec_car./30) &&...
                (s.trafficlight2-y(i)>0))
            ydot(i+s.n_tot)=acceleration(y(i+s.n_tot),0,s.vehicle(s.i,i));
            ydot(i)=y(i+s.n_tot);
        end
    end    
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------




%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%PROVIDES THE ACCELERATIONS +++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function a=acceleration(u,v,vehicle)
global s;
switch vehicle
    case 1
        if(v>u) peak_acc = s.acc_car;elseif(u>=v) peak_acc=-s.dec_car;end
    case 2
        if(v>u) peak_acc = s.acc_hgv;elseif(u>=v) peak_acc=-s.dec_hgv;end   
    case 3
        if(v>u) peak_acc = s.acc_ups;elseif(u>=v) peak_acc=-s.dec_ups;end        
    otherwise
end
    
a=abs(v-u)./30*peak_acc.*5; % scale acceleration so that peak is when velocity is zero
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

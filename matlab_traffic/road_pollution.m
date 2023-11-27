%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% THIS MATLAB SCRIPT JUST DOES A TRAFFIC SIMULATION AND THEN CALCULATES THE
% DISTRIBUTION OF POLLUTION USING A GAUSSIAN PLUME MODEL
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

runModel=true;

zlevel=0; % ground level
% run model if needed
if runModel
    load random_numbers
    disp('Running model')
    [t,y,num_pass,lights,vehicle,timer]=greenlanetraffic(rands01,1);
end


ypos=y(:,1:330);
ymax=max(ypos(:));
% now loop over all vehicles and their positions to calculate the map of
% pollution
[r,c]=size(vehicle);

% Loop over time-steps
ncalcs=0;
for i=1:r
    disp(['Time: ',num2str(t(i))]);
    
    % loop over each vehicle
    for j=1:c
        if vehicle(i,j)==0
            continue
        end
        
        % use gaussian plume model to calculate emission based on position
        posx=y(i,j);
        speedx=y(i,j+c); % m/s
        if vehicle(i,j)==1
            emission=10; % a car g/s - needs to be researched
%             emission=10.*speedx./1609.344; % g/s if it is 10 g/mile 
%                                            % (multiply by miles per second)
        elseif vehicle(i,j)==2
            emission=100; % a diesel g/s - needs to be researched
%             emission=100.*speedx./1609.344; % g/s if it is 100 g/mile 
%                                            % (multiply by miles per second)
        end
        
        
        [pollution,xs,ys]=gaussian_plume_model([posx,500],emission,zlevel,ymax);
        if ncalcs==0
            pollution_store=zeros(size(pollution));
        end
        pollution_store=pollution_store+pollution.*(t(2)-t(1));
        ncalcs=ncalcs+1;
    end
end
pollution_store=pollution_store./(t(end)-t(1));


% note, the -8 and 0 means 10^-8 to 10^0 in 9 steps:
% [c,h]=contour(xs,ys,pollution_store,10);
% clabel(c,h)
pcolor(xs,ys,pollution_store);shading flat
xlabel('x position along road (m)')
ylabel('y position perpendicular to road (m)')
colorbar('horiz')
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
    [t,y,num_pass,lights,vehicle,timer]=princessparkwaytraffic(rands01,1);
end

xrec=[200, 800, 1000];
yrec=[100, 300, 400];
nrec=length(yrec);

[r,c]=size(y);
ypos=y(:,1:c/2);
ymax=max(ypos(:));
% now loop over all vehicles and their positions to calculate the map of
% pollution
[r,c]=size(vehicle);

pollution_rec=zeros(r,nrec);

% call just to define xs
[pollution,xs,ys]=gaussian_plume_model([0,500],10,zlevel,ymax);

% Loop over time-steps
for i=1:r
    disp(['Time: ',num2str(t(i))]);
    
    ncalcs=0;
    pollution_store=zeros(size(pollution));
    % loop over each vehicle
    for j=1:c
        if vehicle(i,j)==0
            continue
        end
        
        % use gaussian plume model to calculate emission based on position
        posx=y(i,j);
        speedx=y(i,j+c); % m/s
        if vehicle(i,j)==1
%             emission=10; % a car g/s - needs to be researched
            emission=10.*speedx./1609.344; % g/s if it is 10 g/mile 
                                           % (multiply by miles per second)
        elseif vehicle(i,j)==2
%             emission=100; % a diesel g/s - needs to be researched
            emission=100.*speedx./1609.344; % g/s if it is 100 g/mile 
                                           % (multiply by miles per second)
        end
        
        
        [pollution,xs,ys]=gaussian_plume_model([posx,500],emission,zlevel,ymax);
        pollution_store(:,:)=pollution_store(:,:)+pollution; % add up for every car
        ncalcs=ncalcs+1;
    end
    % find the position of the receptors
    for k=1:nrec
        ipos=interp1(xs(1,:),1:length(xs(1,:)),xrec(k),'nearest',length(xs(1,:)));
        jpos=interp1(ys(:,1),1:length(ys(:,1)),yrec(k),'nearest',length(ys(:,1)));
        pollution_rec(i,k)=pollution_store(jpos,ipos);
    end
end


% note, the -8 and 0 means 10^-8 to 10^0 in 9 steps:
% [c,h]=contour(xs,ys,pollution_store,10);
% clabel(c,h)
plot(t,pollution_rec);
xlabel('time')
ylabel('Pollution level (g/m^{-3})')

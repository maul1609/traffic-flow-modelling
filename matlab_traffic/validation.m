%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% THIS MATLAB SCRIPT JUST DOES SEVERAL SIMULATIONS LOOKING AT THE 
% PROPORTION OF HGV TRAFFIC, THEN CALCULATES THE NOISE
% LEVELS. IT THEN CALCULATES L10 (the 90th percentile of NOISE IN DBS)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


load random_numbers



p=[0.01 0.05 0.1 0.5];
L10=zeros(size(p));

for i=1:length(p)
    disp(['Running model with ',num2str(p(i)),' fraction of HGV traffic'])
    [t,y,num_pass,lights,vehicle,timer]=greenlanetraffic(rands01,1,p(i));
    sound_level=noise_pollution_hgvs(t,y,vehicle);
    L10(i)=prctile(sound_level,90);
end

plot(p,L10,'+-','markersize',40)
xlabel('Proportion of traffic that is HGVs')
ylabel('L10');hold on;

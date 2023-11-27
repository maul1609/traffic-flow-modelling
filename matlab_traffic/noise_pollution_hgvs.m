function sound_level=noise_pollution_hgvs(t,y,vehicle)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% This MATLAB script simulates the noise levels from simulated vehicles
% travelling down a road
%
% in order to estimate the increase in noise that hgvs have on Green Lane
% assume that the sound always adds constructively
%
% Call as follows (note you must have run a simulation first)
%
% >> sound_level=noise_pollution_hgvs(t,y,vehicle)
% the sound_level variable contains noise level in dBs
%--------------------------------------------------------------------------



% Reference position on Green Lane - houses are here
house_location1=450;
house_location2=25; 


% http://www.google.co.uk/url?sa=t&rct=j&q=hgv%20cars%20decibels%20noise%20levels%20&source=web&cd=7&ved=0CFAQFjAG&url=http%3A%2F%2Fwww.transportresearchfoundation.co.uk%2FPDF%2FPPR202_Characteristics%2520of%2520vehicles%2520producing%2520excessive%2520noise%2520and%2520ground-borne%2520vibration-%2520phase%25201%2520.pdf&ei=kvaFT-XWNoqT8gOY77S0Bw&usg=AFQjCNGy9kIuAUEukCfD4wsAAIj3jIW8kw&cad=rja
cars_decibel=inline('65','x');%inline('(81.5-79)/(145-70)*x+79-(81.5-79)/(145-70)*70','x');65;

% 
% cars_decibel=inline('12.161*log(x)+26.308','x'); % Figure 6 of Steven et al 2005
cars_decibel=inline('13.257*log(x)+20.624','x'); % Figure 6 of Steven et al 2005
cars_reference=10; % 70 db at 15 m
cars_reference=7.5; % 70 db at 15 m
hgvs_decibel=inline('75','x');;%inline('(88.5-80)/(140-70)*x+80-(88.5-80)/(140-70)*70','x');75;
hgvs_decibel=inline('0.1233*x+75.769','x'); % Figure 13 of Steven et al 2005
hgvs_reference=10; % 84 db at 15 m.
hgvs_reference=7.5; % 84 db at 15 m.

pref=20e-6; % Pa reference

% loop through all vehicles and calculate the decibel level at the house
for i=1:length(t)
    ind1=find(vehicle(i,:)==1); % cars
    ind2=find(vehicle(i,:)==2); % hgvs
    
    % Distance of sound sources from house - Pythagoras' theorem - these
    % are root mean square values
    distance1=sqrt((y(i,ind1)-house_location1).^2+house_location2.^2);
    distance2=sqrt((y(i,ind2)-house_location1).^2+house_location2.^2);

    car_db=cars_decibel(y(i,ind1+330).*3600./1000);
    if(car_db<60) car_db=60;end;
    hgv_db=hgvs_decibel(y(i,ind2+330).*3600./1000);
    if(hgv_db<75) hgv_db=75;end;
    
    cars_pref=pref.*10.^(car_db./20);
    hgvs_pref=pref.*10.^(hgv_db./20);


    try
        % Now add up all cars and hgvs, using superposition principle.
        pressure_amplitudes=sqrt(2).*[(cars_pref.*cars_reference./distance1) ...
            (hgvs_pref.*hgvs_reference./distance2)];
        
        rand1=2.*pi.*rand(1000,length(pressure_amplitudes)); 
                    % randomize the phases of each wave arriving at house

                    % find the average of the total pressure wave
        total_pressure=...
            sqrt(mean(sum(...
                repmat(pressure_amplitudes,[1000 1]).*sin(rand1),2).^2));
    catch
        % find the average of the total pressure wave
        total_pressure=0;
    end
    pressure_level=total_pressure;
    
    sound_level(i)=20.*log10(pressure_level./pref);
end

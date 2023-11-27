%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% THIS MATLAB SCRIPT ANIMATES THE RESULTS FROM THE TRAFFIC MODEL, COMPARING
% TWO SETS OF SIMULATIONS
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
runModel=true;
printFigures=false;


% run model if needed
if runModel
    load random_numbers
    disp('Running first model - HGVs not turning')
    [t1,y1,num_pass1,lights1,vehicle1,timer1]=princessparkwaytraffic(rands01,0);
    disp('Running second model - HGVs turning')
    [t2,y2,num_pass2,lights2,vehicle2,timer2]=princessparkwaytraffic(rands01,1);
end
[r,c]=size(y1);
ypos=y1(:,1:c/2);
ymax=max(ypos(:));


for i=1:length(t1)
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % First plot: HGVs not turning
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subplot(121)
    ind=find(y1(i,1:c/2)>0 & vehicle1(i,1:c/2)==1);
    hold off;
    if(length(ind))
        plot(zeros(length(ind),1),y1(i,ind),'.',...
            'markersize',6,'color','red');
    end
    hold on;
    ind=find(y1(i,1:c/2)>0 & vehicle1(i,1:c/2)==2);
    if(length(ind))
        plot(zeros(length(ind),1),y1(i,ind),'.',...
            'markersize',20,'color','green');
    end
    ylim([0 ymax]);xlim([-10 5]);
    text(-8,0,'Start');
    text(-8,1431,'Turn off');
    text(-8,578.9,'Traffic light');
    text(-8,893.1,'Traffic light');
    set(gca,'xtick',[]);

    if(lights1(i,1)==0)
        plot(-1,578.9,'gp','markersize',8);
    else
        plot(-1,578.9,'rp','markersize',8);
    end

    if(lights1(i,2)==0)
        plot(-1,893.1,'gp','markersize',8);
    else
        plot(-1,893.1,'rp','markersize',8);
    end
    title(['without hgvs (time (s): ',num2str(t1(i)),')']);
    %----------------------------------------------------------------------
    
    
    
    
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Second plot: HGVs turning
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subplot(122)
    ind=find(y2(i,1:c/2)>0 & vehicle2(i,1:c/2)==1);
    hold off;
    if(length(ind))
    plot(zeros(length(ind),1),y2(i,ind),'.','markersize',6,'color','red');
    end
    hold on;
    ind=find(y2(i,1:c/2)>0 & vehicle2(i,1:c/2)==2);
    if(length(ind))
        plot(zeros(length(ind),1),y2(i,ind),'.',...
            'markersize',20,'color','green');
    end
    ylim([0 ymax]);xlim([-10 5]);
    text(-8,10,'Start');
    text(-8,1431,'Turn off');
    text(-8,578.9,'Traffic light');
    text(-8,893.1,'Traffic light');
    set(gca,'xtick',[]);

    if(lights2(i,1)==0)
        plot(-1,578.9,'gp','markersize',8);
    else
        plot(-1,578.9,'rp','markersize',8);
    end

    if(lights2(i,2)==0)
        plot(-1,893.1,'gp','markersize',8);
    else
        plot(-1,893.1,'rp','markersize',8);
    end
    title(['with hgvs (time (s): ',num2str(t2(i)),')']);
    %----------------------------------------------------------------------
    
    
    % CHECK IF FIGURES SHOULD BE PRINTED TO FILE
    if printFigures
        print('-dpng',['output/gla',num2str(i,'%03d'),'.png'])
    end
    
    pause(0.5); % pause briefly
end

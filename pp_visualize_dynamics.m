function pp_visualize_dynamics(filename,folder)
% This version is for cluster, where Xinyun delete the "MPEG-4" profile in VideoWriter function, otherwise an error occurs 
% Also, only remain the particle movie

    % Check if file has already been created
    if isfile(strcat(folder,"/PP_",filename))
        disp(strcat("File exists. Skipping... ",filename))
        return % Stop function to skip to next loop iteration
    else
        disp(strcat("RUNNING... ",filename))
    end



%% Visualize particle dynamics

fprintf(1, 'Now reading %s\n', filename);

load(strcat(folder,"/",filename))


% figure('visible','off')
% surf(peaks(256));
% CM = colormap(bluewhitered);
% colormap;
% vid = (CM(:,1) == 1) & (CM(:,2) == 1) & (CM(:,3) == 1)
% CM
% Ckeyboard

%% First loop in time to get Triangulation -> bin size

% Need to do the average of N/2 and then loop in time

    meandistance = zeros(size(ti,1),1);

    ri = sqrt(xi.^2+yi.^2);

    % Max radius to fix figures axes
    rmax = ceil(max(max(ri)));
    xLimits = [-1 1]*rmax;  
    yLimits = xLimits;  
   
    


% In meandistance we have the mean distance of the N/2 particles thar are
% closer to each other. The binsize is 3 times this distance so there are
% ~3x3=9 particles per bin (continuum limit)
binsize = mean(meandistance)*8;

border_mean_r = zeros(size(ti));
border_r_deformation = zeros(size(ti));



 movifilename = strcat(folder,'/movie_particles_',filename(1:end-4));
 vidfile2 = VideoWriter(movifilename);                               
 vidfile2.Quality = 75;
 vidfile2.FrameRate = 10;

 open(vidfile2);
 fig4 = figure('visible','off');

 
for j=1:size(ti,1)
    % disp(strcat('Loop 2 j = ', num2str(j)))
	
	% --------------------
    % Plot cloud of points with border
    set(0,'CurrentFigure',fig4)
    set(gcf,'color','w');

    scatter(xi(j,:),yi(j,:),"filled",'b'); 
	hold on
    scatter(0,0,100,"filled",'g'); % Origin
    scatter(mean(xi(j,:)),mean(yi(j,:)),100,"filled",'r'); % Center of Mass
    xlim(xLimits)
    ylim(yLimits)    
    axis square
    hold off

    frame = getframe(gcf);
    writeVideo(vidfile2, frame);
end

 close(vidfile2);
    
end

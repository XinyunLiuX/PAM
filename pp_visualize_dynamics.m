function pp_visualize_dynamics(filename,folder)
% This version is for cluster, where Xinyun delete the "MPEG-4" profile in VideoWriter function, otherwise an error occurs 
% Also, only remain the particle movie

    % Check if file has already been created
    if isfile(strcat(folder,'/movie_particles_',filename(1:end-4), '.avi'))
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


xmin = floor(min(min(xi)));
xmax = ceil(max(max(xi)));
ymin = floor(min(min(yi)));
ymax = ceil(max(max(yi)));  
% Max radius to fix figures axes


 movifilename = strcat(folder,'/movie_particles_',filename(1:end-4));
 vidfile = VideoWriter(movifilename);                               
 vidfile.Quality = 75;
 vidfile.FrameRate = 10;

 open(vidfile); 

% Prepare figure (offscreen)
fig = figure('Visible','off');
ax = axes(fig);
xlim(ax, [xmin xmax]*1.1);
ylim(ax, [ymin ymax]*1.1);
set(fig,'Color','w');
axis equal

scatterPts = scatter(ax, xi(1,:), yi(1,:), 40, 'b', 'filled');

scatterPts.LineWidth = 0.6;
scatterPts.MarkerEdgeColor = 'k';
scatterPts.MarkerFaceAlpha = 0.5;
scatterPts.MarkerFaceColor = [0 0.4470 0.7410];

for j=1:size(ti,1)
    scatterPts.XData = xi(j,:);
    scatterPts.YData = yi(j,:);
    frame = getframe(gcf);
    writeVideo(vidfile, frame);
end

 close(vidfile);
 close(fig);   

end

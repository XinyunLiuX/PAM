clear all
clc
close all


%% --------------
runner_PAM
myfiles = dir(strcat(folder,'/omega*.mat'));
 
 
 parfor k = 1:length(myfiles)
     filename = myfiles(k).name;
     pp_visualize_dynamics(filename,folder)
 end
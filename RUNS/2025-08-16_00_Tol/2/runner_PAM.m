% Simulations with parametric active matter in a harmonic well in 2D

clear all
clc
close all

% Directory
folder = "RES";
mkdir(folder)
BASE_FOLDER = "../../../";
addpath(BASE_FOLDER)


% PARAMETERS FOR BATCH
% We are going to do a sweep in omega & epsilon
batch_omega   = 2:0.1:4;
batch_epsilon = 0:0.1:0.9;

[batch_omega,batch_epsilon] = meshgrid(batch_omega,batch_epsilon);

numsimulations = size(batch_omega(:),1);

for i=1:numsimulations

    % Model parameters
    N       = 700;                % Number of particles
    tau     = 1.0;                % Relaxation time
    epsilon = batch_epsilon(i);   % Forcing amplitude (needs to be < 1)
    omega   = batch_omega(i);     % Forcing frequency
    T       = 2*pi/omega;         % Forcing period
    alpha   = 1;                  % Interaction force constant
    p_repul = 2;                  % Interaction force power (must be integer > 1)

    % Numerical tolerances
    RelTol = 1e-6;
    AbsTol = 1e-6;
    
    % Collect parameters in structure
    p = struct('N', {N}, 'tau', {tau}, 'epsilon', {epsilon},...
        'omega', {omega}, 'alpha', {alpha}, 'p', {p_repul}, 'RelTol', {RelTol}, 'AbsTol', {AbsTol});
    
    % Simulation time
    times_per_period  = 20;           % Number of points in time saved per period
    number_of_periods = 200;          % Number of total periods for simulation
    tspan = 0:(T/times_per_period):(T*number_of_periods);

    % Filename
    filename = strcat("omega_",num2str(omega,'%.2f'),"_epsilon_",num2str(epsilon,'%.2f'),".mat");
    
    % Check if simulation already done
    if isfile(strcat(folder,"/",filename))
        disp(strcat("Simulation exists. Skipping... ",filename))
        continue % Skip to next loop iteration
    else
        disp(strcat("RUNNING... ",filename))
    end

    % Initial conditions
    rng("shuffle")

    % Random radial coordinate and angle
    R = 10;
    r = R^2.*rand(N,1);
    theta = 2*pi*rand(N,1);

    x0 = sqrt(r).*cos(theta);
    y0 = sqrt(r).*sin(theta);
    u0 = -1 + 2*rand(N,1);  %Random velocity in interval (-1,1)
    v0 = -1 + 2*rand(N,1);  %Random velocity in interval (-1,1)

    y0 = [x0;y0;u0;v0];

    % ODE solver
    opts = odeset('RelTol',RelTol,'AbsTol',AbsTol);
    [t_sol,y] = ode45(@(t,y) active_particles_in_well(t,y,p), tspan, y0, opts);

    % Results (Rows is time, and columns is particle id)
    xi_sol = y(:,      1:N);
    yi_sol = y(:,  N+1:2*N);
    ui_sol = y(:,2*N+1:3*N);
    vi_sol = y(:,3*N+1:4*N);


    % Save data for last 6 periods
    saverange = (size(t_sol,1)-6*times_per_period):size(t_sol,1);

    ti = t_sol(saverange);
    xi = xi_sol(saverange,:);
    yi = yi_sol(saverange,:);
    ui = ui_sol(saverange,:);
    vi = vi_sol(saverange,:);

    % Function to save file within parfor loop
    parsave(strcat(folder,"/",filename), p,ti,xi,yi,ui,vi)
    
end


function parsave(fname, p,ti,xi,yi,ui,vi)
    save(fname, 'p','ti','xi','yi','ui','vi')
end
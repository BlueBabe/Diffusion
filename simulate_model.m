%% DESCRIPTION
%
% Author: Svetlana Lockwood
%
% This file is a wrapper to run the diffusion model
% (diffusion_model.m). The diffusion model describes the spread
% and persistance of antimicrobial resistance (AMR) at the population level. The
% model investigates relationship between (a) population density (high,
% low) and variable fitness cost of antimicrobial trait in bacteria (high,
% low), and (b) fitness cost and variable use levels of antimicrobial
% medications. The population is simulated on a square NxN grid. At each
% time step the model calculates the number of people who have AMR bacteria
% in them. These numbers are output into the respective files. To compute
% percentages, divide them by N^2. The model file is partitioned in 5 parts.
% The first part sets the model macro parameters. The macro parameters are 
% common to all the subsequent four parts. The next four parts
% simulate conditions described in (a) and (b), each part writes its own
% CSV file on the disk. A GIF video file is recorded and displayed for each
% of the four loops.
%
%% DEPENDENCIES
%
% Requires: diffusion_model.m
%
%% SETTING MACRO PARAMETERS
% Macro variables
time_step = 0.05; % time step
time_max = 10; % maximum time to run the model
N = 40; % side of the simulation square; total population = N^2
init_prev = 0.6; % initial level of AMR in the population
sparsity_coef = 1; % link sparsity, 1 is default (i.e. no effect)
reach_radius = 1; % link enrichment, 1 is default (i.e. no effect)

%% HIGH POPULATION DENSITY, VARIABLE FITNESS COST (Figure 1)
coef_new_cases_antibiotic_use = 0.05; % percent of new antibotic use cases

density_coef = 1.0; % high population density, fixed
fitness_coef = [0.05, 0.2]; % low and high fitness cost 

% X - stores number of AMR bacteria carriers
X = zeros(length(fitness_coef), round(time_max/time_step)+1);

% simulation loop for fixed density and variable fitness cost
for i=1:length(fitness_coef)
    i
    video_file = strcat('fitness_cost_', num2str(fitness_coef(i)), '_density_coef_', num2str(density_coef), '.gif');
    X(i, :) = diffusion_model(N, sparsity_coef, init_prev, density_coef, ...
        video_file, time_step, time_max, reach_radius, fitness_coef(i), coef_new_cases_antibiotic_use);
end
csvwrite('results1.csv', X)

%% LOW POPULATION DENSITY, VARIABLE FITNESS COST (Figure 1)
coef_new_cases_antibiotic_use = 0.05; % percent of new antibotic use cases

density_coef = 0.2; % low population density, fixed
fitness_coef = [0.05, 0.2]; % low and high fitness cost 

% X - stores number of AMR bacteria carriers
X = zeros(length(fitness_coef), round(time_max/time_step)+1);

% simulation loop for fixed density and variable fitness cost
for i=1:length(fitness_coef)
    i
    video_file = strcat('fitness_cost_', num2str(fitness_coef(i)), '_density_coef_', num2str(density_coef), '.gif');
    X(i, :) = diffusion_model(N, sparsity_coef, init_prev, density_coef, ...
        video_file, time_step, time_max, reach_radius, fitness_coef(i), coef_new_cases_antibiotic_use);
end
csvwrite('results2.csv', X)

%% LOW FITNESS COST, VARIABLE ANTIBIOTIC USE (Figure 2)
density_coef = 0.6; % medium population density

fitness_coef = 0.05; % low fitness cost, fixed
coef_new_cases_antibiotic_use = [0.1, 0.01]; % high and low percent of new antibotic use cases

% X - stores number of AMR bacteria carriers
X = zeros(length(coef_new_cases_antibiotic_use), round(time_max/time_step)+1);

% simulation loop for fixed fitness cost and variable antibiotic use
for i=1:length(coef_new_cases_antibiotic_use)
    i
    video_file = strcat('fitness_cost_', num2str(fitness_coef), '_percent_new_', num2str(coef_new_cases_antibiotic_use(i)), '.gif');
    X(i, :) = diffusion_model(N, sparsity_coef, init_prev, density_coef, ...
        video_file, time_step, time_max, reach_radius, fitness_coef, coef_new_cases_antibiotic_use(i));
end
csvwrite('results3.csv', X)

%% HIGH FITNESS COST, VARIABLE ANTIBIOTIC USE (Figure 2)
density_coef = 0.6; % medium population density

fitness_coef = 0.2; % high fitness cost, fixed
coef_new_cases_antibiotic_use = [0.1, 0.01]; % high and low percent of new antibotic use cases

% X - stores number of AMR bacteria carriers
X = zeros(length(coef_new_cases_antibiotic_use), round(time_max/time_step)+1);

% simulation loop for fixed fitness cost and variable antibiotic use
for i=1:length(coef_new_cases_antibiotic_use)
    i
    video_file = strcat('fitness_cost_', num2str(fitness_coef), '_percent_new_', num2str(coef_new_cases_antibiotic_use(i)), '.gif');
    X(i, :) = diffusion_model(N, sparsity_coef, init_prev, density_coef, ...
        video_file, time_step, time_max, reach_radius, fitness_coef, coef_new_cases_antibiotic_use(i));
end

csvwrite('results4.csv', X)


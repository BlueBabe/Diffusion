%% DIFFUSION MODEL
% 
% Atuhor: Svetlana Lockwood
% 
% This file runs the diffusion model for the given time. The diffusion is
% implemented iteratively updating system at every time step.
%
%% INPUT PARAMETERS
% N - side of simulation grid, e.g. if N=40, then total people = N^2=1600
% sparsity_coef - Enforces link sparsity: not used, default value is 1
% init_prev - Initial prevalence of antimicrobial resistance (AMR) at the
%             population level, percent [0..1]
% density_coef - population density
% video_file - filename to record video of diffusion, must be .gif, string
% time_step - time step
% time_max - maximum time to run the model
% reach_radius - Enforces link saturation; not used, default value is 1
% fitness_coef - fitness cost of antimicrobial trait in bacteria
% percent_new - percent of new antibotic use cases
%
%% OUTPUT PARAMEETRS
%
% res - number of people with AMR bacteria at each time step, vector
%
%% DEPENDENCIES
%
% Requires: create_adjacency_matrix.m
%
function [res] = diffusion_model(N, sparsity_coef, init_prev, density_coef, ...
    video_file, time_step, time_max, reach_radius, fitness_coef, percent_new)
%% DESCRIPTION
% An example of diffusion over time through a graph. 
% The graph in this example is constructed on a 2D discrete grid, with 
% points on the grid connected to their eight neighbors. Over time, 
% the exponential decay acts to distribute AMR bacteria throughout the 
% entire grid.

%%
% SYSTEM MACROPARAMETERS

% The number of pixels along a dimension of the image: 
% 1 pixel = 1 household. Ex.: if N=40, number of people = N^2=1600.
% N = 40;

% Microbial max log10 load per gram (or per specified unit weight).
max_intensity = 10;

%Percent of natural or induced microbial load clearence
% clearance_rate = 0.1;

%Frequency with which new sporadic "infection" cases occur. For example, if
%frequency = 5, that means new  cases occur after each five time steps.
frequency = 5;

% Coefficients for the main loop
delta_t = time_step; % time step
max_time = time_max; % max time to run the main loop

% Color to display image
colormap jet
% The household adjacency matrix
Adj = create_adjacency_matrix(N, reach_radius, sparsity_coef);

%% PREPARING FOR MAIN LOOP
Deg = diag(sum(Adj, 2)); % Compute the degree matrix
L = Deg - Adj; % Compute the laplacian matrix in terms of the degree and adjacency matrices
[V, D] = eig(L); % Compute the eigenvalues/vectors of the laplacian matrix
D = diag(D);

%Initial condition: set initial prevalence with max load at random
% locations.
C0 = zeros(N, N);
init_number = round(init_prev*N*N);
n = randperm(N*N, init_number);
C0(n) = max_intensity;
C0 = C0(:);

% Transform the initial condition into the coordinate system 
% of the eigenvectors
C0V = V'*C0;
Phi_current = C0V;
PhiV = V*Phi_current;

%% MAIN LOOP

% Counts loop iterations for recording each time step output and (random)
% time events such as sporadic new cases of infection, etc.
count = 0;

% Record output for each time step
res = zeros(size(0:delta_t:max_time));

% Loop through times and decay each initial component
for t = 0:delta_t:max_time
   count = count + 1;
   % Exponential decay for each component
   Phi_current = V\PhiV(:);
   Phi_current = Phi_current.*exp(-D*density_coef*delta_t);

   % Transform from eigenvector coordinate system to original coordinate system
   PhiV = V*Phi_current;
   % Clearing "infection"
   PhiV = PhiV * (1-fitness_coef);
   PhiV = reshape(PhiV, N, N);
   % Select only people with with microbial load > 10^2
   [m, ~] = size(PhiV(PhiV > 2));
   % Record the number of people whose load is grater 10^2
   res(count) = m;
   % Compute the number of newly infected nodes. Current schema: the amount
   % of newly (randomly) infected nodes is proportional to the current
   % number of people with load > 10^2 (m).
   new_nodes = round(percent_new*m);

   % Infect new nodes: infection occur only after "frequency" number of time
   % steps. Newly infected nodes are set to max microbial load 10^10.
   if mod(count, frequency) == 0
       n = randperm(N*N, new_nodes);
       PhiV(n) = max_intensity;
   end
   
   % Display the results and write to GIF file
   imagesc(PhiV);
   caxis([0, 10]);
   title(sprintf('Time: t = %3.2f; mean infection level = %.3f', t, mean(PhiV(:))));
   frame = getframe(1);
   im = frame2im(frame);
   [imind, cm] = rgb2ind(im, 256);
   if t == 0
      imwrite(imind, cm, video_file, 'gif', 'Loopcount', inf, 'DelayTime', 0.1); 
   else
      imwrite(imind, cm, video_file, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
   end
end



end

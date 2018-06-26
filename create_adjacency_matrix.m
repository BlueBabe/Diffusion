%% DESCRIPTION
% This function creates the adjacency matrix with specified parameters. 
% reach_radius and sparsity_coef currently not used, default values 1.
%
function Adj = create_adjacency_matrix(N, reach_radius, sparsity_coef)
%% INPUT PARAMETERS
% N - side of simulation grid, e.g. if N=40, then total people = N^2=1600
% sparsity_coef - Inforces sparsity in connectivity to emulate variable 
%                 density. For example, if sparsity_coef = 2, then every 
%                 2nd link removed, if sparsity_coef = 3, then every 3rd 
%                 link removed, etc.; not used, default value is 1  
% reach_radius - similarily to sparsity_coef but enforces link saturation; 
%                not used, default value is 1
%
%% OUTPUT
% Writes and displays the image of adjacency matrix on the disk.
%
%% DEPENDENCIES
%
% None.
%
%% SETTING UP ADJACENCY MATRIX
Adj = zeros(N*N, N*N); % The adjacency matrix

count = 0; % counts every cell in the adjacency matrix
for x = 1:N
    for y = 1:N
        count = count + 1;
        index = (x-1)*N + y;
        for ii = -reach_radius:reach_radius % this is y offset
            for jj = -reach_radius:reach_radius % this is x offset
                newx = x + jj;
                newy = y + ii;
                if ~(ii == 0 && jj == 0)
                   if newx > 0 && newx <= N && newy > 0 && newy <= N
                       index2 = (newx-1)*N + newy;
%                        if (mod(count, sparsity_coef) == 0) || (mod(count, sparsity_coef+1) == 0) % Inforces sparsity 
%                            Adj(index, index2) = 1;
%                        end
                       if mod(count, sparsity_coef) == 0 % Inforces sparsity 
                           Adj(index, index2) = 1;
                       end
                   end
                end
            end
        end
    end
end

sum(sum(Adj))
imagesc(Adj);
saveas(gcf, strcat('adjacency_reach_radius_', num2str(reach_radius), '.png'));
end

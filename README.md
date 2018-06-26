# Diffusion

The files in this folder are designed to run the diffusion model. All files are written in Matlab language. No special Matlab packages are required to run these files. They should also run in Octave. The main file is simulate_model.m. 

The diffusion model here is adopted to describe the spread and persistance of antimicrobial resistance (AMR) at the population level. The model investigates relationship between (a) population density (high, low) and variable fitness cost of antimicrobial trait in bacteria (high, low), and (b) fitness cost and variable use levels of antimicrobial medications. The population is simulated on a square NxN grid. At each time step the model calculates the number of people who have AMR bacteria in them. These numbers are output into the respective files. To compute percentages, divide them by N^2. All three Matlab files shoud be placed in the same folder.

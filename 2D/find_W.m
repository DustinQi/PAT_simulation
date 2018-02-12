%Signal_n_recon_2D
% =========================================================================
% SIMULATION
% =========================================================================

% load the initial pressure distribution from an image and scale the magnitude
% p0_magnitude = 3;
% p0 = p0_magnitude*loadImage('EXAMPLE_source_one.png');

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
% Nx = 420;           % number of grid points in the x (row) direction
% Ny = 561;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction  [m]
dy = 0.1e-3;        % grid point spacing in the y direction  [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% kgrid.t_array=0:2e-8:1.1980e-5;

% resize the image to match the size of the computational grid and assign
% to the source input structure
source.p0 = resize(absorb*10, [Nx, Ny]);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;

% define a centered circular sensor
sensor_radius = 5e-3;   % [m]
num_sensor_points = 50; %50;
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% =========================================================================
% VISUALISATION
% =========================================================================

% % plot the initial pressure and sensor distribution
% figure;
% imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, source.p0 + cart2grid(kgrid, sensor.mask), [-1 1]);
% colormap(getColorMap);
% ylabel('x-position [mm]');
% xlabel('y-position [mm]'); 
% axis image;
% 
% % plot the simulated sensor data
% figure;
% imagesc(sensor_data, [-1, 1]);
% colormap(getColorMap);
% ylabel('Sensor Position');
% xlabel('Time Step');
% colorbar;

% ===============================================================================================
% time-reversal reconstruction
% ===============================================================================================
% create the time array
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% %dustin 2017/8/30
% kgrid.t_array=0:2e-8:1.1980e-5;

% reset the initial pressure
source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% run the time reversal reconstruction
p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% plot the reconstructed init pressure
figure('Name','Recon');
imagesc(p0_recon,[-1,1]);
colormap(getColorMap);
title('reconstruction')
colorbar;
clearvars -except absorb p0_recon kgrid fW
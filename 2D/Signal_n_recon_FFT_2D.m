% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
PML_size = 20;     
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
% Nx = 420;           % number of grid points in the x (row) direction
% Ny = 561;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction  [m]
dy = 0.1e-3;        % grid point spacing in the y direction  [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

source.p0 = resize(absorb*10, [Nx, Ny]);
% define the properties of the propagation medium

% resize the image to match the size of the computational grid and assign
% to the source input structure
medium.sound_speed = 1500;           % [m/s]
 
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, source.p0, true);

% define a binary line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;

% create the time array
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PMLSize', PML_size, 'PlotPML', false, 'Smooth', false};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% reconstruct the initial pressure
p_xy = kspaceLineRecon(sensor_data.', dy, dt, medium.sound_speed, 'Plot', true, 'PosCond', true);

% define a second k-space grid using the dimensions of p_xy
[Nx_recon, Ny_recon] = size(p_xy);
kgrid_recon = makeGrid(Nx_recon, dt*medium.sound_speed, Ny_recon, dy);

% resample p_xy to be the same size as source.p0
p_xy_rs = interp2(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)), p_xy, kgrid.y, kgrid.x - min(kgrid.x(:)));

% =========================================================================
% VISUALISATION
% =========================================================================

% % plot the initial pressure and sensor distribution
% figure;
% imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, source.p0 + sensor.mask*disc_magnitude);
% colormap(getColorMap);
% ylabel('x-position [mm]');
% xlabel('y-position [mm]');
% axis image;
% colorbar;
% 
% % plot the simulated sensor data
% figure;
% imagesc(sensor_data, [-1, 1]);
% colormap(getColorMap);
% ylabel('Sensor Position');
% xlabel('Time Step');
% colorbar;

% plot the reconstructed initial pressure 
figure;
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, p_xy_rs);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
colorbar;

% % plot a profile for comparison
% figure;
% plot(kgrid.y_vec*1e3, source.p0(disc_x_pos, :), 'k-', kgrid.y_vec*1e3, p_xy_rs(disc_x_pos, :), 'r--');
% xlabel('y-position [mm]');
% ylabel('Pressure');
% legend('Initial Pressure', 'Reconstructed Pressure');
% axis tight;
% set(gca, 'YLim', [0 5.1]);
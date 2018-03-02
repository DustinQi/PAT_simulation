% 3D FFT Reconstruction For A Planar Sensor Example(need absorb)

Nx = 64;
PML_size = 10;
% Nx = Nx - 2*PML_size;
Ny = Nx;
Nz = Nx;               
dx = 0.2e-3;            
dy = dx;                
dz = dx;                
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
medium.sound_speed = 1500;	% [m/s]

% define initial pressure
ball_magnitude = 1; 
p0 = absorb * ball_magnitude;
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, p0, true);

% define a binary planar sensor
sensor.mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor.mask(1, :, :) = 1;

% create the time array
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% set the input arguements
input_args = {'PMLSize', PML_size, 'PMLInside', false, 'PlotPML', false, 'Smooth', false, 'DataCast', 'single'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% reshape sensor data to y, z, t
sensor_data_rs = reshape(sensor_data, Ny, Nz, length(kgrid.t_array));

% reconstruct the initial pressure
p_xyz = kspacePlaneRecon(sensor_data_rs, kgrid.dy, kgrid.dz, dt, medium.sound_speed,...
    'DataOrder', 'yzt', 'PosCond', true, 'Plot', false);

% =========================================================================
% VISUALISATION
% =========================================================================
% define a k-space grid using the dimensions of p_xyz
[Nx_recon, Ny_recon, Nz_recon] = size(p_xyz);
kgrid_recon = makeGrid(Nx_recon, dt*medium.sound_speed, Ny_recon, kgrid.dy, Nz_recon, kgrid.dz);

% define a k-space grid with the same z-spacing as p0
[Nx_p0, Ny_p0, Nz_p0] = size(source.p0);
kgrid_interp = makeGrid(Nx_p0, kgrid.dx, Ny_p0, kgrid.dy, Nz_p0, kgrid.dz);

% resample the p_xyz to be the same size as p0; for a matrix indexed as 
% [M, N, P], the axis variables passed to interp3 must be given in the 
% order N, M, P
p_xyz_rs = interp3(kgrid_recon.y - min(kgrid_recon.y(:)), kgrid_recon.x - min(kgrid_recon.x(:)), kgrid_recon.z - min(kgrid_recon.z(:)), p_xyz, kgrid_interp.y - min(kgrid_interp.y(:)), kgrid_interp.x  - min(kgrid_interp.x(:)), kgrid_interp.z - min(kgrid_interp.z(:)));

% % plot the initial pressure and sensor surface in voxel form
% voxelPlot(double(p0 | sensor.mask));
% set(gca, 'Projection', 'perspective');
% view([0, 99]);
% 
% % plot the initial pressure
% figure;
% plot_scale = [-10 10];
% subplot(2, 2, 1), imagesc(kgrid_interp.y_vec*1e3, kgrid_interp.x_vec*1e3, squeeze(source.p0(:, :, Nz/2)), plot_scale);
% title('x-y plane');
% axis image;
% subplot(2, 2, 2), imagesc(kgrid_interp.z_vec*1e3, kgrid_interp.x_vec*1e3, squeeze(source.p0(:, Ny/2, :)), plot_scale);
% title('x-z plane');
% axis image;
% xlabel('(All axes in mm)');
% subplot(2, 2, 3), imagesc(kgrid_interp.z_vec*1e3, kgrid_interp.y_vec*1e3, squeeze(source.p0(Nx/2, :, :)), plot_scale);
% title('y-z plane');
% axis image;
% colormap(getColorMap);

% plot the reconstructed initial pressure
figure;
subplot(2, 2, 1), imagesc(kgrid_interp.y_vec*1e3, kgrid_interp.x_vec*1e3, squeeze(p_xyz_rs(:, :, Nz/2)));
title('x-y plane');
axis image;
subplot(2, 2, 2), imagesc(kgrid_interp.z_vec*1e3, kgrid_interp.x_vec*1e3, squeeze(p_xyz_rs(:, Ny/2, :)));
title('x-z plane');
axis image;
xlabel('(All axes in mm)');
subplot(2, 2, 3), imagesc(kgrid_interp.z_vec*1e3, kgrid_interp.y_vec*1e3, squeeze(p_xyz_rs(Nx/2, :, :)));
title('y-z plane');
axis image;
colormap(getColorMap);

% % view reconstruction slice by slice
% flyThrough(p_xyz_rs);
% 测试随机N_d个探测器进行重建得到的结果
% 进行tt次模拟,得到out_sensors和out_sensor_mask，以及tt次的均值out_means

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
p0_binary = absorb * ball_magnitude;

% smooth the initial pressure distribution and restore the magnitude
p0 = smooth(kgrid, p0_binary,true);

% assign to the source structure
source.p0 = p0;

% define a Cartesian spherical sensor
sensor_radius = 0.95 * Nx*dx/2;       % [m]
center_pos = [0, 0, 0];     % [m]
num_sensor_points = 64;
sensor_mask_init = makeCartSphere(sensor_radius, num_sensor_points, center_pos, false);
sensor_mask=sensor_mask_init;

%-----------------------------------------------------------------------------------------------------
N_d = 6;
tt=10; 
index_d=1:num_sensor_points;

out_sensors=zeros(64,64,64,tt);
out_sensor_mask=zeros(3,N_d,tt);

for i=1:tt
    
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
p0_binary = absorb * ball_magnitude;

% smooth the initial pressure distribution and restore the magnitude
p0 = smooth(kgrid, p0_binary,true);

% assign to the source structure
source.p0 = p0;

% define a Cartesian spherical sensor
sensor_radius = 0.95 * Nx*dx/2;       % [m]
center_pos = [0, 0, 0];     % [m]
num_sensor_points = 64;
sensor_mask_init = makeCartSphere(sensor_radius, num_sensor_points, center_pos, false);
sensor_mask=sensor_mask_init;
    
%=================================================
index_d=1:num_sensor_points;
sensor_mask=sensor_mask_init;
K=randperm(length(index_d));
bb = (K(1:(num_sensor_points-N_d)));
sensor_mask(:,bb)=[];  
out_sensor_mask(:,:,i)=sensor_mask;
%=================================================

sensor.mask = sensor_mask;
% create the time array
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% set the input arguements
input_args = {'PMLSize', PML_size, 'PMLInside', false, 'PlotPML', true, ...
    'Smooth', false, 'DataCast', 'single', 'CartInterp', 'nearest'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% reset the initial pressure
source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% run the time-reversal reconstruction
p0_recon = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% create a binary sensor mask of an equivalent continuous sphere 
sensor_radius_grid_points = round(sensor_radius/kgrid.dx);
binary_sensor_mask = makeSphere(kgrid.Nx, kgrid.Ny, kgrid.Nz, sensor_radius_grid_points);

% assign to the sensor structure
sensor.mask = binary_sensor_mask;

% interpolate data to remove the gaps and assign to time reversal data
sensor.time_reversal_boundary_data = interpCartData(kgrid, sensor_data, sensor_mask, binary_sensor_mask);

% run the time-reversal reconstruction
p0_recon_interp = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
out_sensors(:,:,:,i)=p0_recon_interp;
clearvars -except absorb tt out_sensors out_sensor_mask kgrid N_d
end
clear K bb tt i
out_means=sum(out_sensors,4);

% %-----------------------------------------------------------------------------------
% ball_x_pos = 32; % y-z plane   	
% ball_y_pos = 32; % x-z plane   	
% ball_z_pos = 32; % x-y plane	
% figure;
% tr1=out_sensors(:,:,:,8);
% subplot(2, 2, 1), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, squeeze(tr1(:, :, ball_z_pos)));
% title('x-y plane');
% axis image;
% subplot(2, 2, 2), imagesc(kgrid.z_vec*1e3, kgrid.x_vec*1e3, squeeze(tr1(:, ball_y_pos, :)));
% title('x-z plane');
% axis image;
% xlabel('(All axes in mm)');
% subplot(2, 2, 3), imagesc(kgrid.z_vec*1e3, kgrid.y_vec*1e3, squeeze(tr1(ball_x_pos, :, :)));
% title('y-z plane');
% axis image;
% colormap(jet);
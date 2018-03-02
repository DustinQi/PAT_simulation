% 3D Time Reversal Reconstruction For A Spherical Sensor Example(need absorb)

% =========================================================================
% SIMULATION
% =========================================================================
% absorb2 = absorb;
% absorb2(absorb<0.55*max(max(max(absorb))))=0;
% create the computational grid
% compute dx and Nx based on a desired x_size(m) and f_max(Hz)
% f_max = 0.4e7;
% c0_min = 1500;
% x_size = 0.01;
% points_per_wavelength = 3;
% %{
% Linear simulations in homogeneous media, simulations can be run using close
% to the Nyquist limit of two points per wavelength.At least three points per 
% wavelength is recommended for the PML.For heterogeneous media,using four or
% more points per wavelength is recommended 
% %}
% dx = c0_min/(points_per_wavelength*f_max);
% Nx = round(x_size/dx);

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
sensor_radius = 0.9 * Nx*dx/2;       % [m]
center_pos = [0, 0, 0];     % [m]
num_sensor_points = 64;
sensor_mask_init = makeCartSphere(sensor_radius, num_sensor_points, center_pos, false);
sensor_mask=sensor_mask_init;

%****************************结构光**********************************

% %选取随机N_d个探测器
% N_d = 6;
% index_d=1:num_sensor_points;
% K=randperm(length(index_d));
% bb = (K(1:(num_sensor_points-N_d)));
% sensor_mask(:,bb)=[];

% %只使用x<0卦限的探测器
% pp = sensor_mask(1,:)>0;
% sensor_mask(:,pp)=[];

% % 只使用一卦限的探测器
% pp = sensor_mask(1,:)<0; 
% sensor_mask(:,pp)=[];
% pp = sensor_mask(2,:)<0;
% sensor_mask(:,pp)=[];
% pp = sensor_mask(3,:)<0;
% sensor_mask(:,pp)=[];

% %只使用一三卦限的探测器
% sensor_mask1=sensor_mask_init;
% sensor_mask2=sensor_mask_init;
% 
% pp = sensor_mask1(1,:)<0; 
% sensor_mask1(:,pp)=[];
% pp = sensor_mask1(2,:)<0;
% sensor_mask1(:,pp)=[];
% pp = sensor_mask1(3,:)<0;
% sensor_mask1(:,pp)=[];
% 
% pp = sensor_mask2(1,:)>0; 
% sensor_mask2(:,pp)=[];
% pp = sensor_mask2(2,:)>0;
% sensor_mask2(:,pp)=[];
% pp = sensor_mask2(3,:)<0;
% sensor_mask2(:,pp)=[];
% sensor_mask=[sensor_mask1,sensor_mask2];
% clear sensor_mask1 sensor_mask2

% %只使用一七卦限的探测器
% sensor_mask1=sensor_mask_init;
% sensor_mask2=sensor_mask_init;
% 
% pp = sensor_mask1(1,:)<0; 
% sensor_mask1(:,pp)=[];
% pp = sensor_mask1(2,:)<0;
% sensor_mask1(:,pp)=[];
% pp = sensor_mask1(3,:)<0;
% sensor_mask1(:,pp)=[];
% 
% pp = sensor_mask2(1,:)>0; 
% sensor_mask2(:,pp)=[];
% pp = sensor_mask2(2,:)>0;
% sensor_mask2(:,pp)=[];
% pp = sensor_mask2(3,:)>0;
% sensor_mask2(:,pp)=[];
% sensor_mask=[sensor_mask1,sensor_mask2];
% clear sensor_mask1 sensor_mask2

% %只使用一三七卦限的探测器
% sensor_mask1=sensor_mask_init;
% sensor_mask2=sensor_mask_init;
% sensor_mask3=sensor_mask_init;
% 
% pp = sensor_mask1(1,:)<0; 
% sensor_mask1(:,pp)=[];
% pp = sensor_mask1(2,:)<0;
% sensor_mask1(:,pp)=[];
% pp = sensor_mask1(3,:)<0;
% sensor_mask1(:,pp)=[];
% 
% pp = sensor_mask2(1,:)>0; 
% sensor_mask2(:,pp)=[];
% pp = sensor_mask2(2,:)>0;
% sensor_mask2(:,pp)=[];
% pp = sensor_mask2(3,:)>0;
% sensor_mask2(:,pp)=[];
% 
% pp = sensor_mask3(1,:)>0; 
% sensor_mask3(:,pp)=[];
% pp = sensor_mask3(2,:)>0;
% sensor_mask3(:,pp)=[];
% pp = sensor_mask3(3,:)<0;
% sensor_mask3(:,pp)=[];
% sensor_mask=[sensor_mask1,sensor_mask2,sensor_mask3];
% clear sensor_mask1 sensor_mask2 sensor_mask3

% %plot sensor_mask
% % select suitable axis scaling factor
% [x_sc, scale, prefix] = scaleSI(max(sensor_mask(:)));
% % create the figure
% figure;
% plot3(sensor_mask(1, :)*scale, sensor_mask(2,:)*scale, sensor_mask(3,:)*scale, '.');
% xlabel(['[' prefix 'm]']);
% ylabel(['[' prefix 'm]']);
% zlabel(['[' prefix 'm]']);
% axis equal;
% grid on;
% box on;
%******************************************************************

% assign to the sensor structure
sensor.mask = sensor_mask;

% % create the time array
% [kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed,0.7);

% set the input arguements
input_args = {'PMLSize', PML_size, 'PMLInside', false, 'PlotPML', true, ...
    'Smooth', false, 'DataCast', 'single', 'CartInterp', 'nearest'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% % plot photoacoustic signal
% subplot(241);
% plot (sensor_data(1,:));
% title('sensor1');
% subplot(242);
% plot (sensor_data(8,:));
% title('sensor8');
% subplot(243);
% plot (sensor_data(16,:));
% title('sensor16');
% subplot(244);
% plot (sensor_data(32,:));
% title('sensor32');
% subplot(245);
% plot (sensor_data(42,:));
% title('sensor42');
% subplot(246);
% plot (sensor_data(48,:));
% title('sensor48');
% subplot(247);
% plot (sensor_data(56,:));
% title('sensor56');
% subplot(248); 
% plot (sensor_data(64,:));
% title('sensor64');

%---------------------------------------------------------------------------------------------------
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

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the initial pressure and sensor surface in voxel form
% voxelPlot(double(p0_binary | cart2grid(kgrid, sensor_mask)));
% view([60, 20]);

ball_x_pos = 32; % y-z plane   	
ball_y_pos = 32; % x-z plane   	
ball_z_pos = 32; % x-y plane	

% % plot the initial pressure
% figure;
% plot_scale = [0 0.5];
% subplot(2, 2, 1), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, squeeze(p0(:, :, ball_z_pos)), plot_scale);
% title('x-y plane');
% axis image;
% subplot(2, 2, 2), imagesc(kgrid.z_vec*1e3, kgrid.x_vec*1e3, squeeze(p0(:, ball_y_pos, :)), plot_scale);
% title('x-z plane');
% axis image;
% xlabel('(All axes in mm)');
% subplot(2, 2, 3), imagesc(kgrid.z_vec*1e3, kgrid.y_vec*1e3, squeeze(p0(ball_x_pos, :, :)), plot_scale);
% title('y-z plane');
% axis image;
% colormap(getColorMap);

% % plot the reconstructed initial pressure
% figure;
% plot_scale = [0, 0.5];
% subplot(2, 2, 1), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, squeeze(p0_recon(:, :, ball_z_pos)), plot_scale);
% title('x-y plane');
% axis image;
% subplot(2, 2, 2), imagesc(kgrid.z_vec*1e3, kgrid.x_vec*1e3, squeeze(p0_recon(:, ball_y_pos, :)), plot_scale);
% title('x-z plane');
% axis image;
% xlabel('(All axes in mm)');
% subplot(2, 2, 3), imagesc(kgrid.z_vec*1e3, kgrid.y_vec*1e3, squeeze(p0_recon(ball_x_pos, :, :)), plot_scale);
% title('y-z plane');
% axis image;
% colormap(getColorMap);

% plot the reconstructed initial pressure
figure;
%plot_scale = [-10 10];
% plot_scale = [-2 2];
% plot_scale = [-1.2 1.2];
subplot(2, 2, 1), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, squeeze(p0_recon_interp(:, :, ball_z_pos)));
title('x-y plane');
axis image;
subplot(2, 2, 2), imagesc(kgrid.z_vec*1e3, kgrid.x_vec*1e3, squeeze(p0_recon_interp(:, ball_y_pos, :)));
title('x-z plane');
axis image;
xlabel('(All axes in mm)');
subplot(2, 2, 3), imagesc(kgrid.z_vec*1e3, kgrid.y_vec*1e3, squeeze(p0_recon_interp(ball_x_pos, :, :)));
title('y-z plane');
axis image;
colormap(jet);



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

% absorb=absorb_3Dtest2;
% absorb = p0_2cy_equal;
absize=size(absorb);

PML_size = 5;
Nx = 30;
Ny = 30;
Nz = 30;               
dx = 0.25e-3;            
dy = 0.25e-3;                 
dz = 0.25e-3;                
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
medium.sound_speed = 1500;	% [m/s]

% define initial pressure
ball_magnitude = 1e6; 
absorb = absorb * ball_magnitude;

p0 = resize(absorb, [Nx,Ny,Nz]);
% smooth the initial pressure distribution and restore the magnitude

p0 = smooth(kgrid, p0,true);

% assign to the source structure
source.p0 = p0;

%% define 3D Cartesian circle sensor-array demo1: cylinder surface
sensor_radius = 0.8 * Nx*dx/2;       % [m]
center_pos = [0, 0];     % [m]
num_sensor_points = 8;
sensor_mask_init = makeCartCircle(sensor_radius, num_sensor_points, center_pos, 2*pi, true);

layer1=ones(1,num_sensor_points);
layer1=layer1.*(1/6)*(dz*Nz);
layer2=ones(1,num_sensor_points);
layer2=layer2.*(1/3)*(dz*Nz);
layer3=ones(1,num_sensor_points);
layer3=layer3.*(0/4)*(dz*Nz);
layer4=ones(1,num_sensor_points);
layer4=layer4.*(-1/6)*(dz*Nz);
layer5=ones(1,num_sensor_points);
layer5=layer5.*(-1/3)*(dz*Nz);

sensor_mask=[sensor_mask_init, sensor_mask_init, sensor_mask_init, sensor_mask_init, sensor_mask_init];
layer=[layer1,layer2,layer3,layer4,layer5];
sensor_mask=[sensor_mask;layer];
clear layer1 layer2 layer3 layer4 layer5 layer

%% define 3D Cartesian circle sensor-array demo2: top&bottom
% sensor_radius = 0.8 * Nx*dx/2;       % [m]
% center_pos = [0, 0];     % [m]
% num_sensor_points1 = 4;
% sensor_mask_init1 = makeCartCircle(sensor_radius, num_sensor_points1, center_pos, 2*pi, true);
% 
% sensor_radius = 0.4 * Nx*dx/2;       % [m]
% center_pos = [0, 0];     % [m]
% num_sensor_points2 = 4;
% sensor_mask_init2 = makeCartCircle(sensor_radius, num_sensor_points2, center_pos, 2*pi, true);
% 
% layer1=ones(1,num_sensor_points1);
% layer1=layer1.*(2/5)*(dz*Nz);
% layer2=ones(1,num_sensor_points2);
% layer2=layer2.*(-2/5)*(dz*Nz);
% 
% sensor_mask=[sensor_mask_init1, sensor_mask_init2, sensor_mask_init1, sensor_mask_init2];
% layer=[layer1,layer1,layer2,layer2];
% sensor_mask=[sensor_mask;layer];
% clear layer1 layer2 layer
%% ****************************结构光**********************************

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
%% ******************************************************************

% assign to the sensor structure
sensor.mask = sensor_mask;

% % create the time array
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
% [kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed,0.7);

% set the input arguements
input_args = {'PMLSize', PML_size, 'PMLInside', true, 'PlotPML', true, ...
    'Smooth', false, 'DataCast', 'single', 'CartInterp', 'nearest'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

%% plot photoacoustic signal
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
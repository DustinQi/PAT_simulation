function sphere = makeCartSphere(radius, num_points, center_pos, plot_sphere)
%MAKECARTSPHERE     Create a 3D Cartesian sphere.
%
% DESCRIPTION:
%       makeCartSphere creates a 3 x num_points array of the Cartesian
%       coordinates of points evenly distributed over a sphere using the
%       Golden Section Spiral method. 

% USAGE:
%       sphere = makeCartSphere(radius, num_points)
%       sphere = makeCartSphere(radius, num_points, center_pos)
%       sphere = makeCartSphere(radius, num_points, center_pos, plot_sphere)

% INPUTS:
%       radius          - sphere radius [m]
%       num_points      - number of points on the sphere

% OPTIONAL INPUTS:
%       center_pos      - [x, y, z] position of the circle center [m] 
%                         (default = [0, 0, 0])
%       plot_sphere     - Boolean controlling whether the Cartesian points
%                         are plotted (default = false)
%
% OUTPUTS:
%       sphere          - 3 x num_points array of Cartesian coordinates
%
% ABOUT:
%       modified by     - Bradley Treeby
%       date            - 10th June 2009
%       last update     - 17th July 2009
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also cart2grid, makeCartCircle, makeSphere

% check for plot_sphere input
if nargin < 4
    plot_sphere = false;
end

% check for center_pos input
if nargin < 3
    cx = 0;
    cy = 0;
    cz = 0;
else
    cx = center_pos(1);
    cy = center_pos(2);
    cz = center_pos(3);    
end

% generate angle functions using the Golden Section Spiral method
inc = pi * (3-sqrt(5));
off = 2/num_points;
k = 0:num_points-1;
y = k * off - 1 + (off/2);
r = sqrt(1 - (y.^2));
phi = k * inc;

% create the sphere
sphere = radius.*[cos(phi).*r; y ;sin(phi).*r];

% offset if needed
sphere(1, :) = sphere(1, :) + cx;
sphere(2, :) = sphere(2, :) + cy;
sphere(3, :) = sphere(3, :) + cz;

% plot results
if plot_sphere
    
    % select suitable axis scaling factor
    [x_sc, scale, prefix] = scaleSI(max(sphere(:)));  %#ok<ASGLU>
    
    % create the figure
    figure('Name','position of detectors');
   
    plot3(sphere(1, :)*scale, sphere(2,:)*scale, sphere(3,:)*scale, '.');
    xlabel(['[' prefix 'm]']);
    ylabel(['[' prefix 'm]']);
    zlabel(['[' prefix 'm]']);
    axis equal;
    grid on;
    box on;
    title('position of detectors');
end
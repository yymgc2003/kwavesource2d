location = csvread(fullfile(config.location_seedfiles_path, 'location1.csv'));
config = jsondecode(fileread('../config.json'));
% Each row of 'location' corresponds to a 3D vector [x, y, z]
vec = location(1, :); % Get the first row as a 3D vector
% You can print the whole vector directly using fprintf if the number of elements is known.
fprintf('location: %f %f %f\n', vec);

% If the number of elements is variable, you can use:
fprintf('location: %s\n', num2str(vec));

Nx=config.grid.Nx;
Ny=config.grid.Ny;
Nz=config.grid.Nz;
dx=config.grid.dx;
dy=config.grid.dy;
dz=config.grid.dz;
pipe_center_x=config.pipe.center_x;
pipe_center_y=config.pipe.center_y;
pipe_center_z=config.pipe.center_z;
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask([Nx/4, Nx/2, 3*Nx/4], Ny/2, Nz/2) = 1;
cx = Nx/2+106;            % 
cy = Ny/2;            % 
cz = Nz/2;
outer_r_mm = 16;  %  [mm]
inner_r_mm = 13;  %  [mm]
outer_r = round((outer_r_mm * 1e-3) / dx);  % 
inner_r = round((inner_r_mm * 1e-3) / dx);  % 

bx=round(inner_r*vec(1))+pipe_center_x;
by=round(inner_r*vec(2))+pipe_center_y;
bz=round(inner_r*vec(3))+pipe_center_z;

%glass_mask = makeBall(Nx, Ny, Nz, bx, by, bz, radius_pts);
fprintf('bx: %d, by: %d, bz: %d\n', bx, by, bz);
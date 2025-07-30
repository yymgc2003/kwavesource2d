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
num_particles = 1;
sphere_radius = 1.25; % [mm]
sphere_volume = 4/3*pi*sphere_radius^3;
height = dz*Nz*1e3; % [mm]
fprintf('height: %f\n', height);
surface_pipe = pi*inner_r_mm^2; % [mm^2]
V_sphere = num_particles*sphere_volume; % [mm^3]
V_pipe = surface_pipe*height; % [mm^3]
fraction = V_sphere/V_pipe;
fprintf('sphere_radius: %f\n', sphere_radius);
fprintf('sphere_volume: %f\n', sphere_volume);
fprintf('V_sphere: %f\n', V_sphere);
fprintf('V_pipe: %f\n', V_pipe);
fprintf('fraction: %f\n', fraction);
medium.sound_speed = config.medium.water.sound_speed * ones(Nx, Ny, Nz);
medium.density = config.medium.water.density * ones(Nx, Ny, Nz);
medium.alpha_coeff = config.medium.water.alpha_coeff * ones(Nx, Ny, Nz);
medium.alpha_power = config.medium.water.alpha_power;
medium.BonA = 6;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
t_end = config.simulation.t_end;
cfl = config.simulation.CFL;
kgrid.makeTime(medium.sound_speed, cfl, t_end);
source_strength = config.source.source_strength;
tone_burst_freq = config.source.tone_burst_freq;
tone_burst_cycles = config.source.tone_burst_cycles;

    % Generate input signal
    % To plot input_signal, you should plot it against its time axis.
    input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
    %original_pulse = load('/home/matsubara/Scripts/kwavesource/src/init_pulse_upsampled.mat');
    %input_signal = original_pulse.pulse_upsampled;
    %input_signal = load('/home/matsubara/Scripts/kwavesource/src/init_pulse.mat');
    %input_signal = input_signal.pulse;
    input_signal = (source_strength / (config.medium.water.sound_speed * config.medium.water.density)) * input_signal;
    % Display the class type of input_signal
    fprintf('The class of input_signal is: %s\n', class(input_signal));
    t = (0:length(input_signal)-1) * kgrid.dt * 1e6; % Time in microseconds
    fprintf('input_signal.shape: %s\n', mat2str(size(input_signal)));
    figure(1);
    plot(t, input_signal , 'b-'); % Plot input_signal in MPa vs time in microseconds
    % The x-axis is time [μs], and the y-axis is pressure [MPa].
    xlabel('Time [\mus]');
    ylabel('Pressure [MPa]');
    ylim([-2 2]);
    title('tutorial signal');
    grid on;
    saveas(gcf, fullfile(config.save_full_path, 'sandbox.png'));








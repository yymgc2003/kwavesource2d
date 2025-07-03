clearvars;

% simulation settings
DATA_CAST = 'gpuArray-single';

config = jsondecode(fileread('config.json'));
save_path = config.save_path;
% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 20;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]
PML_Z_SIZE = 10;            % [grid points]

% set total number of grid points not including the PML
Nx = config.grid.Nx - 2*PML_X_SIZE;    % [grid points]
Ny = config.grid.Ny - 2*PML_Y_SIZE;    % [grid points]
Nz = config.grid.Nz - 2*PML_Z_SIZE;     % [grid points]

% set desired grid size in the x-direction not including the PML
%x = 40e-3;                  % [m]

% calculate the spacing between the grid points

dx=config.grid.dx; dy=config.grid.dy; dz=config.grid.dz;
% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
medium.sound_speed = config.medium.water.sound_speed * ones(Nx, Ny, Nz);      % [m/s]
medium.density = config.medium.water.density * ones(Nx, Ny, Nz);          % [kg/m^3]
medium.alpha_coeff = config.medium.water.alpha_coeff * ones(Nx, Ny, Nz);      % [dB/(MHz^y cm)]
medium.alpha_power = config.medium.water.alpha_power;
medium.BonA = 6;

% create the time array
t_end = config.simulation.t_end;                  % [s]
cfl= config.simulation.CFL;
kgrid.makeTime(medium.sound_speed, cfl, t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = config.source.source_strength;          % [Pa]
tone_burst_freq = config.source.tone_burst_freq;        % [Hz]
tone_burst_cycles = config.source.tone_burst_cycles;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
%input_signal = (source_strength / (medium.sound_speed .* medium.density)) .* input_signal;
input_signal = (source_strength / (config.medium.water.sound_speed * config.medium.water.density)) * input_signal;
% =========================================================================
% DEFINE THE ULTRASOUND transducer_transmit
% =========================================================================

% physical properties of the transducer_transmit
transducer_transmit.number_elements = 72;    % total number of transducer_transmit elements
transducer_transmit.element_width = 1;       % width of each element [grid points]
transducer_transmit.element_length = 12;     % length of each element [grid points]
transducer_transmit.element_spacing = 0;     % spacing (kerf width) between the elements [grid points]
transducer_transmit.radius = inf;            % radius of curvature of the transducer_transmit [m]

% calculate the width of the transducer_transmit in grid points
transducer_transmit_width = transducer_transmit.number_elements * transducer_transmit.element_width ...
    + (transducer_transmit.number_elements - 1) * transducer_transmit.element_spacing;

% use this to position the transducer_transmit in the middle of the computational grid
transducer_transmit.position = round([10, Ny/2 - transducer_transmit_width/2, Nz/2 - transducer_transmit.element_length/2]);

% properties used to derive the beamforming delays
transducer_transmit.sound_speed = 1540;                  % sound speed [m/s]
transducer_transmit.focus_distance = 20e-3;              % focus distance [m]
transducer_transmit.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
transducer_transmit.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer_transmit.transmit_apodization = 'Rectangular';    
transducer_transmit.receive_apodization = 'Rectangular';

% define the transducer_transmit elements that are currently active
transducer_transmit.active_elements = zeros(transducer_transmit.number_elements, 1);
transducer_transmit.active_elements(21:52) = 1;

% append input signal used to drive the transducer_transmit
transducer_transmit.input_signal = input_signal;

% create the transducer_transmit using the defined settings
transducer_transmit = kWaveTransducer(kgrid, transducer_transmit);

% print out transducer_transmit properties
transducer_transmit.properties;

transducer_trans.number_elements = 90;    % total number of transducer elements
transducer_trans.element_width = 1;       % width of each element [grid points/voxels]
transducer_trans.element_length = 12;     % length of each element [grid points/voxels]
transducer_trans.element_spacing = 0;     % spacing (kerf  width) between the elements [grid points/voxels]
transducer_trans.radius = inf;            % radius of curvature of the transducer [m]
transducer_width = transducer_trans.number_elements * transducer_trans.element_width ...
    + (transducer_trans.number_elements - 1) * transducer_trans.element_spacing;
transducer_trans.position = round([Nx-10, Ny/2 - transducer_width/2, Nz/2 - transducer_trans.element_length/2]);
transducer_trans.sound_speed = config.medium.water.sound_speed;
transducer_trans.focus_distance = 25e-3;
transducer_trans.elevation_focus_distance = 19e-3;
transducer_trans.steering_angle = 0;
transducer_trans.transmit_apodization = 'Rectangular';
transducer_trans.receive_apodization = 'Rectangular';
transducer_trans.active_elements = zeros(transducer_trans.number_elements, 1);
transducer_trans.active_elements(21:52) = 1;
transducer_trans = kWaveTransducer(kgrid, transducer_trans);
transducer_trans.properties;
% =========================================================================
% DEFINE SENSOR MASK
% =========================================================================

% create a binary sensor mask with four detection positions
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask([Nx/4, Nx/2, 3*Nx/4], Ny/2, Nz/2) = 1;
cx = Nx/2+106;            % 
cy = Ny/2;            % 
cz = Nz/2;
outer_r_mm = 16;  %  [mm]
inner_r_mm = 13;  %  [mm]
outer_r = round((outer_r_mm * 1e-3) / dx);  % 
inner_r = round((inner_r_mm * 1e-3) / dx);  % 

[Xg, Yg] = ndgrid(1:Nx, 1:Ny);
ring2d = sqrt((Xg-cx).^2 + (Yg-cy).^2);
%ring2d = sqrt((Xg-cx).^2 + (Yg-cy).^2*16);
ringMask = (ring2d <= outer_r) & (ring2d >= inner_r);   % Nx×Ny logical

pipe_mask = repmat(ringMask, [1 1 Nz]);   % Nx×Ny×Nz  
medium.sound_speed(pipe_mask == 1) = config.medium.vinyl.sound_speed;
medium.density(pipe_mask == 1) = config.medium.vinyl.density;

% =========================================================================
% RUN THE SIMULATION
% =========================================================================
display_mask = transducer_transmit.all_elements_mask | pipe_mask;
% set the input settings
input_args = {'DisplayMask', display_mask, ...
    'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
    'RecordMovie', true, ...
    'MovieName', fullfile(save_path, 'liquid_only.avi'), ...
    'DataCast', DATA_CAST, 'PlotScale', [-1/2, 1/2] * source_strength};

% run the simulation
sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer_transmit, transducer_transmit, input_args{:});
save(fullfile(save_path, 'liquid_only.mat'), 'sensor_data', 'kgrid', '-v7.3');
%[f_input, as_input] = spect([input_signal, zeros(1, 2 * length(input_signal))], 1/kgrid.dt);
%[~, as_1] = spect(sensor_data(1, :), 1/kgrid.dt);
%[~, as_2] = spect(sensor_data(2, :), 1/kgrid.dt);
%[f, as_3] = spect(sensor_data(3, :), 1/kgrid.dt);
scan_line = transducer_transmit.scan_line(sensor_data);
figure(1);
plot(kgrid.t_array * 1e6, scan_line * 1e-6, 'b-');
xlabel('Time [\mus]');
ylabel('Pressure [MPa]');
ylim([-2 2]);
title('Signal from Transducer transmit');
grid on;
saveas(gcf, fullfile(save_path, 'signal_liquid_only_reflector.png'));


figure(21); clf;
hold on;

% Permute mask dimensions to (x, y, z) for correct visualization
transmit_mask = permute(transducer_transmit.all_elements_mask, [2,1,3]);
sensor_mask   = permute(sensor.mask, [2,1,3]);
pipe_mask_p   = permute(pipe_mask, [2,1,3]);

% Convert to double for visualization
transmit_mask_double = double(transmit_mask);
sensor_mask_double   = double(sensor_mask);
pipe_mask_double     = double(pipe_mask_p);

% Plot transmit transducer mask (blue)
if any(transmit_mask_double(:))
    p1 = patch(isosurface(transmit_mask_double, 0.5));
    isonormals(transmit_mask_double, p1);
    set(p1, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
end

% Plot pipe mask (green, with lower threshold if needed)
if any(pipe_mask_double(:))
    p3 = patch(isosurface(pipe_mask_double, 0.3)); % try 0.3 if 0.5 does not work
    isonormals(pipe_mask_double, p3);
    set(p3, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end

% Set axis and view
axis equal;
xlim([1, Nx]);
ylim([1, Ny]);
zlim([1, Nz]);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Transducer and Pipe Mask Visualization');
view(80, 30);
camlight;
lighting gouraud;
grid on;
hold off;

saveas(gcf, fullfile(save_path, 'liquid_only.png'));


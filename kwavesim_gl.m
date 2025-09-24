function kwavesim_gl(config_file, location_csv, locnum_str)
% Main simulation logic for k-Wave, extracted for modular use.

    config = jsondecode(fileread(config_file));
    save_logs_path = fullfile(config.save_full_path, 'logs');
    save_data_path = fullfile(config.save_full_path, 'data');
    % Check if save_logs_path exists, and create it if it does not exist
    if ~exist(save_logs_path, 'dir')
        mkdir(save_logs_path);
        fprintf('Directory %s created.\n', save_logs_path);
    end
    % Check if save_data_path exists, and create it if it does not exist
    if ~exist(save_data_path, 'dir')
        mkdir(save_data_path);
        fprintf('Directory %s created.\n', save_data_path);
    end
    % Simulation settings
    DATA_CAST = 'gpuArray-single';

    % PML size settings
    PML_SIZE = config.simulation.pml_size; % [grid points]

    % Number of grid points excluding PML
    Nx = config.grid.Nx;
    Ny = config.grid.Ny;
    Nz = config.grid.Nz;

    % Grid spacing
    dx = config.grid.dx; dy = config.grid.dy; dz = config.grid.dz;

    % Create k-space grid
    kgrid = kWaveGrid(Nx, dx, Ny, dy);

    % Medium properties
    medium.sound_speed = config.medium.water.sound_speed * ones(Nx, Ny);
    medium.density = config.medium.water.density * ones(Nx, Ny);
    medium.alpha_coeff = config.medium.water.alpha_coeff * ones(Nx, Ny);
    medium.alpha_power = config.medium.water.alpha_power;
    medium.BonA = config.medium.water.BonA;

    % Create time array
    t_end = config.simulation.t_end;
    cfl = config.simulation.CFL;
    kgrid.makeTime(medium.sound_speed, cfl, t_end);

    % Input signal properties
    source_strength = config.source.source_strength;
    tone_burst_freq = config.source.tone_burst_freq;
    tone_burst_cycles = config.source.tone_burst_cycles;

    % Generate input signal
    % Read upsampled pulse from mat file
    %original_pulse = load('/home/matsubara/Scripts/kwavesource/src/pulse_4000_upsampled.mat');
    %input_signal = original_pulse.pulse_upsampled;
    
    input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

    % Scale by impedance
    input_signal = (source_strength / (config.medium.water.sound_speed * config.medium.water.density)) * input_signal(1:end-10);
    fprintf('input_signal: %f\n', size(input_signal));

    % Source mask
    source_diameter = config.source.diameter/dx;
    source_radius = round(source_diameter/2);
    dist_pipe_source = round(config.source.distance_pipe_source/dy);
    source.uy = input_signal;
    source.u_mask = zeros(Nx, Ny);
    source.u_mask(Nx/2-source_radius:Nx/2+source_radius, ...
        Ny/2-dist_pipe_source) = 1;

    % Sensor mask
    sensor.mask = zeros(Nx, Ny);
    sensor.mask(Nx/2-source_radius:Nx/2+source_radius, ...
        Ny/2-dist_pipe_source) = 1;
    sensor.mask(Nx/2-source_radius:Nx/2+source_radius, ...
        Ny/2+dist_pipe_source) = 1;
    sensor.record = {'p'};

    % Pipe mask
    cx = round(config.pipe.center_x/dx);
    cy = round(config.pipe.center_y/dy);
    outer_r = round(config.pipe.outer_radius/dx);
    inner_r = round(config.pipe.inner_radius/dy);
    [Xg, Yg] = ndgrid(1:Nx, 1:Ny);
    ring2d = sqrt((Xg-cx-Nx/2).^2 + (Yg-cy-Ny/2).^2);
    ringMask = (ring2d <= outer_r) & (ring2d >= inner_r);

    pipe_mask = repmat(ringMask, [1 1]);
    medium.sound_speed(pipe_mask == 1) = config.medium.vinyl.sound_speed;
    medium.density(pipe_mask == 1) = config.medium.vinyl.density;
    medium.alpha_coeff(pipe_mask == 1) = config.medium.vinyl.alpha_coeff;

    % Glass mask
    % Read coordinates from locationX.csv
    location = csvread(location_csv);
    glass_radius = config.simulation.glass_radius;

    % Initialize glass_mask
    gas_mask = zeros(Nx, Ny);

    % Place balls at each coordinate
    % 三次元ではz方向が0.05*128=6.4mmだった。
    % そのため、location(3)は-1:1を-3.2mm:3.2mmに変換
    % 0mm部分で断面を取り、その形状を考える
    for i = 1:size(location,1)
        loc_seed = location(i,:);
        bx = round(inner_r * loc_seed(2)) + cx + Nx/2;
        by = round(inner_r * loc_seed(1)) + cy + Ny/2;
        gas_radius = round(inner_r * loc_seed(3));
        gas_mask = gas_mask | makeDisc(Nx, Ny, bx, by, gas_radius);
    end
    %glass_mask = zeros(Nx, Ny, Nz); %check validity for liquid only setting
    medium.sound_speed(gas_mask == 1) = config.medium.air.sound_speed;
    medium.density(gas_mask == 1) = config.medium.air.density;
    medium.alpha_coeff(gas_mask == 1) = config.medium.air.alpha_coeff;

    figure;
    hold on;
    for i = 1:size(location,1)
        loc_seed = location(i,:);
        loc_seed(1) = config.pipe.inner_radius * loc_seed(1);
        loc_seed(2) = config.pipe.inner_radius * loc_seed(2);
        gas_radius = config.pipe.inner_radius * loc_seed(3);
        rectangle('Position',[loc_seed(2)-gas_radius,...
            loc_seed(1)-gas_radius,...
            2*gas_radius,2*gas_radius], ...
            'Curvature', [1 1], ...
            'FaceColor', [0.8 0.8 0.8], ...
            'EdgeColor', 'black');
        hold on;
    end
    rectangle('Position',[-16e-3,-16e-3,32e-3,32e-3], ...
        'Curvature', [1 1], ...
        'EdgeColor', 'black');
    hold on;
    rectangle('Position',[-13e-3,-13e-3,26e-3,26e-3], ...
        'Curvature', [1 1], ...
        'EdgeColor', 'black');
    hold on;
    axis equal;
    xlim([-19e-3 19e-3]);
    ylim([-19e-3 19e-3]);
    axis off;
    saveas(gcf, fullfile(save_logs_path, ['experimental_setup' locnum_str '.png']));

    % Run simulation
    input_args = {'PlotPML', false, 'PMLSize', PML_SIZE, ...
        'DataCast', DATA_CAST, 'DeviceNum', 1};

    sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
    save(fullfile(save_data_path, ['solid_liquid_reflector' locnum_str '.mat']), 'sensor_data', 'kgrid', '-v7.3');

    % Plot and save signal waveform
    sensor_len = length(sensor_data.p(:,1));
    reflector = sensor_data.p(1:sensor_len/2,:);
    transparent = sensor_data.p(sensor_len/2+1:sensor_len,:);
    reflector = mean(reflector);
    transparent = mean(transparent);
    figure;
    plot(kgrid.t_array * 1e6, reflector * 1e-6, 'b-');
    xlabel('Time [\mus]');
    ylabel('Pressure [MPa]');
    xlim([0 100]);
    ylim([-0.2 0.2]);
    title('Signal from Transducer transmit');
    grid on;
    saveas(gcf, fullfile(save_logs_path, ['signal_solid_liquid_reflector' locnum_str '.png']));
    figure;
    plot(kgrid.t_array * 1e6, transparent * 1e-6, 'b-');
    xlabel('Time [\mus]');
    ylabel('Pressure [MPa]');
    xlim([0 100]);
    ylim([-0.2 0.2]);
    title('Signal from Transducer receiver');
    grid on;
    saveas(gcf, fullfile(save_logs_path, ['signal_solid_liquid_receiver' locnum_str '.png']));
end 
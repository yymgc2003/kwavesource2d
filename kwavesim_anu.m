function kwavesim_anu(config_file, location_csv, locnum_str)
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
        medium.BonA = config.medium.water.BonA * ones(Nx, Ny);
    
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
        medium.BonA(pipe_mask == 1) = config.medium.vinyl.BonA;
    
        % Glass mask
        % Read coordinates from locationX.csv
        location = csvread(location_csv);
        glass_radius = config.simulation.glass_radius;
    
        % Initialize glass_mask
        gas_mask = zeros(Nx, Ny);

        spline_point_num = config.simulation.annular_spline_point_num;
        annular_radius = config.simulation.annular_radius_mean;
        samples = zeros(1,spline_point_num+1);
        for i=1:spline_point_num+1
            samples(1,i) = location(i)*annular_radius;
        end
        theta = linspace(0,2*pi,spline_point_num+1);
        cs = spline(theta, [0 samples 0]);
        for i=1:Ny
            for j = 1:Nx
                cur_x = (i-Ny/2)*dy;
                cur_y = (j-Nx/2)*dx;
                if cur_x > 0
                    cur_theta = atan(cur_y/cur_x);
                elseif cur_x < 0
                    cur_theta = atan(cur_y/cur_x) + pi;
                elseif cur_y > 0
                    cur_theta = pi*0.5;
                elseif cur_y < 0
                    cur_theta = pi*1.5;
                else
                    cur_theta = 0;
                end
                cur_radius = ppval(cs, cur_theta);
                if cur_x^2 + cur_y^2 < cur_radius^2
                    gas_mask(j, i) = 1;
                end
            end
        end
        %glass_mask = zeros(Nx, Ny, Nz); %check validity for liquid only setting
        medium.sound_speed(gas_mask == 1) = config.medium.air.sound_speed;
        medium.density(gas_mask == 1) = config.medium.air.density;
        medium.alpha_coeff(gas_mask == 1) = config.medium.air.alpha_coeff;
    
        fig = figure;
        ax = polaraxes(fig);
        theta_plot = linspace(0, 2*pi, 10*spline_point_num+1);
        rho_plot = ppval(cs, theta_plot);
        hold on;
        polarplot(ax, theta_plot, rho_plot, 'b');
        hold on;
        rho_i = ones(1,10*spline_point_num+1).*13e-3;
        polarplot(ax, theta_plot, rho_i, 'k');
        hold on;
        rho_o = ones(1,10*spline_point_num+1).*16e-3;
        polarplot(ax, theta_plot, rho_o, 'k');
        hold on;
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
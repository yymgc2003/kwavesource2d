function kwavesim3d_gl(config_file, location_csv, locnum_str)
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
    PML_X_SIZE = config.simulation.pml_size; % [grid points]
    PML_Y_SIZE = config.simulation.pml_size/2; % [grid points]
    PML_Z_SIZE = config.simulation.pml_size/2; % [grid points]

    % Number of grid points excluding PML
    Nx = config.grid.Nx - 2*PML_X_SIZE;
    Ny = config.grid.Ny - 2*PML_Y_SIZE;
    Nz = config.grid.Nz - 2*PML_Z_SIZE;

    % Grid spacing
    dx = config.grid.dx; dy = config.grid.dy; dz = config.grid.dz;

    % Create k-space grid
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

    % Medium properties
    medium.sound_speed = config.medium.water.sound_speed * ones(Nx, Ny, Nz);
    medium.density = config.medium.water.density * ones(Nx, Ny, Nz);
    medium.alpha_coeff = config.medium.water.alpha_coeff * ones(Nx, Ny, Nz);
    medium.alpha_power = config.medium.water.alpha_power;
    medium.BonA = config.medium.water.BonA * ones(Nx, Ny, Nz);

    % Create time array
    t_end = config.simulation.t_end;
    cfl = config.simulation.CFL*5/4.2;
    kgrid.makeTime(medium.sound_speed, cfl, t_end);

    % Input signal properties
    source_strength = config.source.source_strength;
    tone_burst_freq = config.source.tone_burst_freq;
    tone_burst_cycles = config.source.tone_burst_cycles;

    cx = round(config.pipe.center_x/dx*1e-3) + 10;
    cy = round(config.pipe.center_y/dy*1e-3) + Ny/2;
    cz = round(config.pipe.center_z/dz*1e-3) + Nz/2;

    % Generate input signal
    % Read upsampled pulse from mat file
    % original_pulse = load('/home/matsubara/Scripts/kwavesource/src/pulse_4000_upsampled.mat');
    %input_signal = original_pulse.pulse_upsampled;
    
    input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

    % Scale by impedance
    input_signal = (source_strength / (config.medium.water.sound_speed * config.medium.water.density)) * input_signal(1:end-10);
    fprintf('input_signal: %f\n', size(input_signal));
    % Define transmit transducer
    % transducer_transmit.number_elements = config.transducer.elements;
    % transducer_transmit.element_width = config.transducer.element_width;
    % transducer_transmit.element_length = config.transducer.element_length;
    % transducer_transmit.element_spacing = config.transducer.element_spacing;
    transducer_transmit.number_elements = round(9e-3/dy);
    transducer_transmit.element_width = config.transducer.element_width;
    transducer_transmit.element_length = round(0.3e-3/dz);
    transducer_transmit.element_spacing = config.transducer.element_spacing;
    transducer_transmit.radius = inf;

    % Transducer width
    transducer_transmit_width = transducer_transmit.number_elements * transducer_transmit.element_width ...
        + (transducer_transmit.number_elements - 1) * transducer_transmit.element_spacing;

    % Position
    transducer_transmit.position = round([10, cy - transducer_transmit_width/2, cz - transducer_transmit.element_length/2]);

    % Beamforming properties
    transducer_transmit.sound_speed = config.transducer.sound_speed;
    transducer_transmit.focus_distance = config.transducer.focus_distance;
    transducer_transmit.elevation_focus_distance = config.transducer.elevation_focus_distance;
    transducer_transmit.steering_angle = config.transducer.steering_angle;

    % Apodization
    transducer_transmit.transmit_apodization = 'Rectangular';
    transducer_transmit.receive_apodization = 'Rectangular';

    % Active elements
    transducer_transmit.active_elements = zeros(transducer_transmit.number_elements, 1);
    transducer_transmit.active_elements(round(transducer_transmit.number_elements*2/9):round(transducer_transmit.number_elements*7/9)) = 1;

    % Input signal
    transducer_transmit.input_signal = input_signal;

    % Create transmit transducer
    transducer_transmit = kWaveTransducer(kgrid, transducer_transmit);

    % Display properties
    transducer_transmit.properties;

    % Define receive transducer
    transducer_receive.number_elements = 90;
    transducer_receive.element_width = 1;
    transducer_receive.element_length = 12;
    transducer_receive.element_spacing = 0;
    transducer_receive.radius = inf;
    transducer_width = transducer_receive.number_elements * transducer_receive.element_width ...
        + (transducer_receive.number_elements - 1) * transducer_receive.element_spacing;
    transducer_receive.position = round([Nx-10, Ny/2 - transducer_width/2, Nz/2 - transducer_transmit.element_length/2]);
    transducer_receive.sound_speed = config.medium.water.sound_speed;
    transducer_receive.focus_distance = 25e-3;
    transducer_receive.elevation_focus_distance = 19e-3;
    transducer_receive.steering_angle = 0;
    transducer_receive.transmit_apodization = 'Rectangular';
    transducer_receive.receive_apodization = 'Rectangular';
    transducer_receive.active_elements = zeros(transducer_receive.number_elements, 1);
    transducer_receive.active_elements(21:52) = 1;
    transducer_receive = kWaveTransducer(kgrid, transducer_receive);

    % Sensor mask
    sensor.mask = zeros(Nx, Ny, Nz);
    sensor.mask([Nx/4, Nx/2, 3*Nx/4], Ny/2, Nz/2) = 1;
    outer_r_mm = config.pipe.outer_radius;
    inner_r_mm = config.pipe.inner_radius;
    outer_r = round((outer_r_mm * 1e-3) / dx);
    inner_r = round((inner_r_mm * 1e-3) / dx);
    [Xg, Yg] = ndgrid(1:Nx, 1:Ny);
    ring2d = sqrt((Xg-cx).^2 + (Yg-cy).^2);
    ringMask = (ring2d <= outer_r) & (ring2d >= inner_r);

    pipe_mask = repmat(ringMask, [1 1 Nz]);
    medium.sound_speed(pipe_mask == 1) = config.medium.vinyl.sound_speed;
    medium.density(pipe_mask == 1) = config.medium.vinyl.density;
    medium.alpha_coeff(pipe_mask == 1) = config.medium.vinyl.alpha_coeff;
    medium.BonA(pipe_mask == 1) = config.medium.vinyl.BonA;

    % Glass mask

    % Initialize glass_mask
    bubble_mask = zeros(Nx, Ny, Nz);
    % Read coordinates from locationX.csv
    location_df = readtable(location_csv);
    location = table2array(location_df);
    
    if config.simulation.flow_pattern == "slug"
        %location: スラグの中心、スラグ長さ、楕円の累乗の値
        bz = round(location(1)*inner_r) + cz;
        major_axis_length = round(location(2)*inner_r);
        minor_axis_length = round(location(3)*inner_r);
        slug_pow_num = location(4);
        for z=1:Nz
            pow_z_relative =abs((z-bz)/major_axis_length)^slug_pow_num;
            for x = cx-inner_r-1:cx+inner_r+1
                for y = cy-inner_r-1:cy+inner_r+1
                    pow_r_relative = (((x-cx)^2+(y-cy)^2)/minor_axis_length^2)^(slug_pow_num/2);
                    if pow_z_relative + pow_r_relative <= 1
                        bubble_mask(x, y, z) = 1;
                    end
                end
            end
        end
    end
    if config.simulation.flow_pattern == "bubble"
        radius_pts = round(config.simulation.glass_radius / dx);

        % Place balls at each coordinate
        for i = 1:size(location,1)
            loc_seed = location(i,:);
            bx = round(inner_r * loc_seed(1)) + cx;
            by = round(inner_r * loc_seed(2)) + cy;
            bz = round(inner_r * loc_seed(3)) + cz;
            radius_pts = round(inner_r*loc_seed(4)/2);
            if bz-radius_pts-1<Nz & bz+radius_pts+1>1
                radius_pts_short = round(inner_r*loc_seed(5)/2);
                %fprintf('bx: %d, by: %d, bz: %d\n', bx, by, bz);
                th1 = loc_seed(6); th2 = loc_seed(7); th3 = loc_seed(8);
                roll = [1,0,0;0,cos(th1),-sin(th1);0,sin(th1),cos(th1)];
                pitch =[cos(th2),0,sin(th2);0,1,0;-sin(th2),0,cos(th2)];
                yaw =  [cos(th3),-sin(th3),0;sin(th3),cos(th3),0;0,0,1];
                Q = roll*pitch*yaw;
                for x=bx-radius_pts-1:bx+radius_pts+1
                    for y=by-radius_pts-1:by+radius_pts+1
                        for z=max(1,bz-radius_pts-1):min(Nz,bz+radius_pts+1)
                            bubble_mask_relative = [x, y, z]-[bx, by, bz];
                            if bubble_mask_relative*Q*diag([1/radius_pts^2,1/radius_pts^2,1/radius_pts_short^2])*transpose(Q)*transpose(bubble_mask_relative)<=1
                                bubble_mask(x, y, z) = 1;
                            end
                        end
                    end
                end
            end
            fprintf('%d th sample generated\n', i);
        end
    end
    if config.simulation.flow_pattern == "slug-bubble"
        % Place balls at each coordinate
        i=0;
        while i < size(location,1)
            i = i+1;
            loc_seed = location(i,:);
            bx = round(inner_r * loc_seed(1)) + cx;
            by = round(inner_r * loc_seed(2)) + cy;
            bz = round(inner_r * loc_seed(3)) + cz;
            radius_pts = round(inner_r*loc_seed(4)/2);
            radius_pts_short = round(inner_r*loc_seed(5)/2);
            %fprintf('bx: %d, by: %d, bz: %d\n', bx, by, bz);
            th1 = loc_seed(6); th2 = loc_seed(7); th3 = loc_seed(8);
            roll = [1,0,0;0,cos(th1),-sin(th1);0,sin(th1),cos(th1)];
            pitch =[cos(th2),0,sin(th2);0,1,0;-sin(th2),0,cos(th2)];
            yaw =  [cos(th3),-sin(th3),0;sin(th3),cos(th3),0;0,0,1];
            Q = roll*pitch*yaw;
            for x=bx-radius_pts-1:bx+radius_pts+1
                for y=by-radius_pts-1:by+radius_pts+1
                    for z=max(1,bz-radius_pts-1):min(Nz,bz+radius_pts+1)
                        bubble_mask_relative = [x, y, z]-[bx, by, bz];
                        if bubble_mask_relative*Q*diag([1/radius_pts^2,1/radius_pts^2,1/radius_pts_short^2])*transpose(Q)*transpose(bubble_mask_relative)<=1
                            bubble_mask(x, y, z) = 1;
                        end
                    end
                end
            end
            fprintf('%d th sample generated\n', i);
        end
    end
    %glass_mask = zeros(Nx, Ny, Nz); %check validity for liquid only setting
    medium.sound_speed(bubble_mask == 1) = config.medium.air.sound_speed;
    medium.density(bubble_mask == 1) = config.medium.air.density;
    medium.alpha_coeff(bubble_mask == 1) = config.medium.air.alpha_coeff;
    medium.BonA(bubble_mask == 1) = config.medium.air.BonA;

    % Run simulation
    display_mask = transducer_transmit.all_elements_mask | pipe_mask | bubble_mask;
    input_args = {'DisplayMask', display_mask, ...
        'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
        'DataCast', DATA_CAST, 'PlotScale', [-1/2, 1/2] * source_strength};
    
    sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer_transmit, transducer_transmit, input_args{:});
    save(fullfile(save_data_path, ['solid_liquid_reflector' locnum_str '.mat']), 'sensor_data', 'kgrid', '-v7.3');

    % Visualize and save experimental setup
    figure(21); clf;
    hold on;

    % Permute mask dimensions for visualization
    transmit_mask = permute(transducer_transmit.all_elements_mask, [2,1,3]);  %NECESSARY FOR VISUALIZATION
    sensor_mask   = permute(sensor.mask, [2,1,3]);
    pipe_mask_p   = permute(pipe_mask, [2,1,3]);
    bubble_mask_p  = permute(bubble_mask, [2,1,3]);
    % Convert to double for visualization
    transmit_mask_double = double(transmit_mask);
    sensor_mask_double   = double(sensor_mask);
    pipe_mask_double     = double(pipe_mask_p);
    bubble_mask_double    = double(bubble_mask_p);

    % Transmit transducer mask (blue)
    if any(transmit_mask_double(:))
        p1 = patch(isosurface(transmit_mask_double, 0.5));
        isonormals(transmit_mask_double, p1);
        set(p1, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    end

    % Pipe mask (green)
    if any(pipe_mask_double(:))
        p3 = patch(isosurface(pipe_mask_double, 0.3));
        isonormals(pipe_mask_double, p3);
        set(p3, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end

    % Glass mask (red)
    if any(bubble_mask_double(:))
        p4 = patch(isosurface(bubble_mask_double, 0.3));
        isonormals(bubble_mask_double, p4);
        set(p4, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end

    % Axis and view settings
    axis equal;
    xlim([1, Nx]);
    ylim([1, Ny]);
    zlim([1, Nz]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Transducer, Sensor, and Pipe Mask Visualization');
    view(80, 20);
    camlight;
    lighting gouraud;
    grid on;
    hold off;

    saveas(gcf, fullfile(save_logs_path, ['slug_experimental_setup' locnum_str '.png']));

    % Plot and save signal waveform
    scan_line = transducer_transmit.scan_line(sensor_data);
    figure(1);
    plot(kgrid.t_array * 1e6, scan_line * 1e-6, 'b-');
    xlabel('Time [\mus]');
    ylabel('Pressure [MPa]');
    ylim([-2 2]);
    title('Signal from Transducer transmit');
    grid on;
    saveas(gcf, fullfile(save_logs_path, ['signal_solid_liquid_reflector' locnum_str '.png']));

    
end 
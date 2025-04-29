% =========================================================================
% k-Wave transducer tutorial
% =========================================================================
clearvars;
close all;
DATA_CAST = 'gpuArray-single';
% -------------------------------------------------------------------------
% 1) 設定ファイルの読み込み
% -------------------------------------------------------------------------
config = jsondecode(fileread('../config.json'));
USE_STATISTICS = true;
% -------------------------------------------------------------------------
% 2) シミュレーション用グリッドの定義
% -------------------------------------------------------------------------
Nx=config.grid.Nx; Ny=config.grid.Ny; Nz=config.grid.Nz;
dx=config.grid.dx; dy=config.grid.dy; dz=config.grid.dz;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
save_path = config.save_path;

% -------------------------------------------------------------------------
% 3) 媒質パラメータ
% -------------------------------------------------------------------------
%medium.sound_speed = config.medium.water.sound_speed;
%medium.density = config.medium.water.density;
%medium.alpha_coeff = config.medium.water.alpha_coeff;
medium.sound_speed = config.medium.water.sound_speed * ones(Nx,Ny,Nz,'single');
medium.density     = config.medium.water.density     * ones(Nx,Ny,Nz,'single');
medium.alpha_coeff = config.medium.water.alpha_coeff * ones(Nx,Ny,Nz,'single');
medium.alpha_power = config.medium.water.alpha_power;
kgrid.makeTime(medium.sound_speed, 0.05, config.simulation.t_end);

cx = Nx/2; cy = Ny/2; cz = Nz/2;
radius_pts = round(2.5e-3 / dx);   
glass_mask1 = makeBall(Nx, Ny, Nz, cx, cy, cz, radius_pts);
glass_mask2 = makeBall(Nx, Ny, Nz, cx, cy+40, cz-40, radius_pts);
glass_mask3 = makeBall(Nx, Ny, Nz, cx-30, cy-40, cz+40, radius_pts);
glass_mask4 = makeBall(Nx, Ny, Nz, cx+30, cy-40, cz+40, radius_pts);
glass_mask5 = makeBall(Nx, Ny, Nz, cx-20, cy-20, cz-20, radius_pts);
glass_mask = glass_mask1 | glass_mask2 | glass_mask3 | glass_mask4 | glass_mask5;

% -------------------------------------------------------------------------
% 4) トランスデューサーの設定
% -------------------------------------------------------------------------
transducer.number_elements = 90;    % total number of transducer elements
transducer.element_width = 1;       % width of each element [grid points/voxels]
transducer.element_length = 12;     % length of each element [grid points/voxels]
transducer.element_spacing = 0;     % spacing (kerf width) between the elements [grid points/voxels]
transducer.radius = inf;            % radius of curvature of the transducer [m]
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;
transducer.position = round([5, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);
transducer.sound_speed = config.medium.water.sound_speed;
transducer.focus_distance = 25e-3;
transducer.elevation_focus_distance = 19e-3;
transducer.steering_angle = 0;
transducer.transmit_apodization = 'Rectangular';
transducer.receive_apodization = 'Rectangular';
transducer.active_elements = zeros(transducer.number_elements, 1);
transducer.active_elements(21:52) = 1;

transducer_trans.number_elements = 90;    % total number of transducer elements
transducer_trans.element_width = 1;       % width of each element [grid points/voxels]
transducer_trans.element_length = 12;     % length of each element [grid points/voxels]
transducer_trans.element_spacing = 0;     % spacing (kerf width) between the elements [grid points/voxels]
transducer_trans.radius = inf;            % radius of curvature of the transducer [m]
transducer_width = transducer_trans.number_elements * transducer_trans.element_width ...
    + (transducer_trans.number_elements - 1) * transducer_trans.element_spacing;
transducer_trans.position = round([Nx-5, Ny/2 - transducer_width/2, Nz/2 - transducer_trans.element_length/2]);
transducer_trans.sound_speed = config.medium.water.sound_speed;
transducer_trans.focus_distance = 25e-3;
transducer_trans.elevation_focus_distance = 19e-3;
transducer_trans.steering_angle = 0;
transducer_trans.transmit_apodization = 'Rectangular';
transducer_trans.receive_apodization = 'Rectangular';
transducer_trans.active_elements = zeros(transducer_trans.number_elements, 1);
transducer_trans.active_elements(21:52) = 1;

cx = Nx/2;            % X 方向の中心
cy = Ny/2;            % Y 方向の中心

% 外径・内径をグリッドポイント数に変換（mmからmに変換してから計算）
outer_r_mm = 16;  % 外径 [mm]
inner_r_mm = 13;  % 内径 [mm]
outer_r = round((outer_r_mm * 1e-3) / dx);  % 外径をグリッドポイント数に変換
inner_r = round((inner_r_mm * 1e-3) / dx);  % 内径をグリッドポイント数に変換

% 2次元の円環マスクを作成
[Xg, Yg] = ndgrid(1:Nx, 1:Ny);
ring2d = sqrt((Xg-cx).^2 + (Yg-cy).^2);
ringMask = (ring2d <= outer_r) & (ring2d >= inner_r);   % Nx×Ny logical

pipe_mask = repmat(ringMask, [1 1 Nz]);   % Nx×Ny×Nz   
medium.sound_speed(pipe_mask == 1) = config.medium.vinyl.sound_speed;
medium.density(pipe_mask == 1) = config.medium.vinyl.density;
medium.alpha_coeff(pipe_mask) = config.medium.vinyl.alpha_coeff;
%medium.alpha_power(pipe_mask) = config.medium.vinyl.alpha_power;

medium.sound_speed(glass_mask == 1) = config.medium.glass.sound_speed;
medium.density(glass_mask == 1) = config.medium.glass.density;
medium.alpha_coeff(glass_mask) = config.medium.glass.alpha_coeff;
%medium.alpha_power(glass_mask) = config.medium.glass.alpha_power;
% -------------------------------------------------------------------------
% 5) ソース波形の設定
% -------------------------------------------------------------------------
source_strength = 1e6;          % [Pa]
tone_burst_freq = 4e6;        % [Hz]
tone_burst_cycles = 4;
source_signal = zeros(size(kgrid.t_array));
T_prf = 1 / config.source.prf;
t_array = kgrid.t_array;

for n = 0:config.source.max_n
    t_start = n * T_prf;
    t_end = t_start + config.source.pulse_length;
    
    if t_start > kgrid.t_array(end)
        break;
    end
    
    idx_on = (t_array >= t_start) & (t_array < t_end);
    source_signal(idx_on) = config.source.magnitude * sin(2*pi * config.source.frequency * (t_array(idx_on) - t_start));
end
transducer.input_signal = source_signal;

% -------------------------------------------------------------------------
% 6) トランスデューサーの初期化
% -------------------------------------------------------------------------
transducer = kWaveTransducer(kgrid, transducer);
transducer_trans = kWaveTransducer(kgrid, transducer_trans);
%transducer.properties;
% -------------------------------------------------------------------------
% 7) センサーの設定
% -------------------------------------------------------------------------
sensor.mask = zeros(Nx, Ny, Ny);
sensor_x = Nx/2 + config.sensor.x_offset;
sensor_y = Ny/2 + config.sensor.y_offset;

% stream the data to disk in blocks of 100 if storing the complete time
% history 
%if ~USE_STATISTICS
%    input_args = [input_args {'StreamToDisk', 100}];
%end

display_mask = transducer.active_elements_mask | transducer_trans.active_elements_mask | pipe_mask | glass_mask;
input_args = {'DisplayMask', display_mask, ...
    'RecordMovie', true, ...
    'MovieName', fullfile(save_path, 'vinyl_pipe_3d.avi'), ...
    'DataCast', DATA_CAST, 'PlotScale', [-1/4, 1/4] * source_strength};
sensor.record = {'p','p_max'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer_trans, input_args{:});

% =========================================================================
% COMPUTE THE BEAM PATTERN USING SIMULATION STATISTICS
% =========================================================================
% Convert logical masks to double for visualization
glass_mask_double = double(glass_mask);
transducer_mask_double = double(transducer.active_elements_mask);
transducer_trans_mask_double = double(transducer_trans.active_elements_mask);
pipe_mask_double = double(pipe_mask);

% Create separate masks for visualization
glass_transducer_mask = glass_mask_double + transducer_mask_double + transducer_trans_mask_double;
pipe_only_mask = pipe_mask_double;


% トランスデューサーからの信号を取得
scan_line = transducer.scan_line(sensor_data);
scan_line_trans = transducer_trans.scan_line(sensor_data);

% トランスデューサー1の信号を可視化
figure(1);
plot(kgrid.t_array * 1e6, scan_line * 1e-6, 'b-');
xlabel('Time [\mus]');
ylabel('Pressure [MPa]');
title('Signal from Transducer 1');
grid on;
saveas(gcf, fullfile(save_path, 'Transducer1_Signal.png'));

% トランスデューサー2の信号を可視化
figure(2); 
plot(kgrid.t_array * 1e6, scan_line_trans * 1e-6, 'r-');
xlabel('Time [\mus]');
ylabel('Pressure [MPa]');
title('Signal from Transducer 2');
grid on;
saveas(gcf, fullfile(save_path, 'Transducer2_Signal.png'));

figure(3);
plot(kgrid.t_array * 1e6, scan_line_trans * 1e-6, 'k-');
xlabel('Time [\mus]');
ylabel('Pressure [MPa]');
title('Scan Line After Beamforming');
saveas(gcf, fullfile(save_path, 'Beamforming.png'));

figure(4);
p1 = patch(isosurface(glass_transducer_mask, 0.5));
isonormals(glass_transducer_mask, p1);
set(p1, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.8);

p2 = patch(isosurface(pipe_only_mask, 0.5));
isonormals(pipe_only_mask, p2);
set(p2, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% Set up the figure
daspect([1 1 1]);
view(80, 30);
camlight;
lighting gouraud;
axis([1 size(glass_transducer_mask,1) 1 size(glass_transducer_mask,2) 1 size(glass_transducer_mask,3)]);
title('Combined Visualization (Transparent Pipe)');
saveas(gcf, fullfile(save_path, 'combined_visualization_transparent.png'));


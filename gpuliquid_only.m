% =========================================================================
% k-Wave liquid -only tutorial (GPU,no movie) 
% =========================================================================
clearvars;
close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));
DATA_CAST = 'gpuArray-single';
% -------------------------------------------------------------------------
% 1) 設定ファイルの読み込み
% -------------------------------------------------------------------------
config = jsondecode(fileread('config.json'));
%USE_STATISTICS = true;
% -------------------------------------------------------------------------
% 2) シミュレーション用グリッドの定義
% -------------------------------------------------------------------------
% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 40;            % [grid points]
PML_Y_SIZE = 40;            % [grid points]
PML_Z_SIZE = 20;            % [grid points]

Nx=config.grid.Nx- 2*PML_X_SIZE; Ny=config.grid.Ny - 2*PML_Y_SIZE; Nz=config.grid.Nz- 2*PML_Z_SIZE;
dx=config.grid.dx; dy=config.grid.dy; dz=config.grid.dz;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
save_path = config.save_path;

% -------------------------------------------------------------------------
% 3) 媒質パラメータ
% -------------------------------------------------------------------------

medium.sound_speed = config.medium.water.sound_speed * ones(Nx,Ny,Nz,'single');
medium.density     = config.medium.water.density     * ones(Nx,Ny,Nz,'single');
medium.alpha_coeff = config.medium.water.alpha_coeff * ones(Nx,Ny,Nz,'single');
medium.alpha_power = config.medium.water.alpha_power;
medium.BonA = 6;
kgrid.makeTime(medium.sound_speed, 0.03, config.simulation.t_end);
% -------------------------------------------------------------------------
% 4) トランスデューサーの設定
% -------------------------------------------------------------------------
transducer_transmit.number_elements = 48;    % total number of transducer_transmit elements
transducer_transmit.element_width = 1;       % width of each element [grid points/voxels]
transducer_transmit.element_length = 12;     % length of each element [grid points/voxels]
transducer_transmit.element_spacing = 0;     % spacing (kerf  width) between the elements [grid points/voxels]
transducer_transmit.radius = inf;            % radius of curvature of the transducer_transmit [m]
transducer_transmit_width = transducer_transmit.number_elements * transducer_transmit.element_width ...
    + (transducer_transmit.number_elements - 1) * transducer_transmit.element_spacing;

transducer_transmit.position = round([PML_X_SIZE+30, Ny/2 - transducer_transmit_width/2, config.grid.Nz/2 - transducer_transmit.element_length/2]);
transducer_transmit.sound_speed = config.medium.water.sound_speed;
transducer_transmit.focus_distance = 25e-3;
transducer_transmit.elevation_focus_distance = Inf;
transducer_transmit.steering_angle = 0;
transducer_transmit.transmit_apodization = 'Hanning';
transducer_transmit.active_elements = zeros(transducer_transmit.number_elements, 1);
transducer_transmit.active_elements(5:43) = 1;

transducer_receive.number_elements = 48;    % total number of transducer_transmit elements
transducer_receive.element_width = 1;       % width of each element [grid points/voxels]
transducer_receive.element_length = 12;     % length of each element [grid points/voxels]
transducer_receive.element_spacing = 0;     % spacing (kerf  width) between the elements [grid points/voxels]
transducer_receive.radius = inf;            % radius of curvature of the transducer_transmit [m]
transducer_transmit_width = transducer_transmit.number_elements * transducer_transmit.element_width ...
    + (transducer_transmit.number_elements - 1) * transducer_transmit.element_spacing;
transducer_receive.position = round([PML_X_SIZE+630, Ny/2 - transducer_transmit_width/2+1, config.grid.Nz/2 - transducer_transmit.element_length/2+1]);
transducer_receive.sound_speed = config.medium.water.sound_speed;
transducer_receive.steering_angle = 0;
transducer_receive.elevation_focus_distance = Inf;
transducer_receive.receive_apodization = 'Hanning';
transducer_receive.active_elements = zeros(transducer_receive.number_elements, 1);
transducer_receive.active_elements(5:43) = 1;

% ビニール円環のマスクを作成 - サイズを調整
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
%medium.alpha_coeff(pipe_mask) = config.medium.vinyl.alpha_coeff;
%medium.alpha_power(pipe_mask) = config.medium.vinyl.alpha_power;

% -------------------------------------------------------------------------
% 5) ソース波形の設定
% -------------------------------------------------------------------------

base_signal = toneBurst(1/kgrid.dt, config.source.tone_burst_freq, config.source.tone_burst_cycles);
base_signal = (config.source.source_strength ./ (config.medium.water.sound_speed * config.medium.water.density)) .* base_signal;
source_signal = zeros(size(kgrid.t_array));
T_prf = 1 / config.source.prf;  
burst_length = length(base_signal);
dt = kgrid.dt;
burst_time = burst_length * dt;

% 最初のパルスのみを設定
t_start = 0; % 開始時刻
idx_start = round(t_start / kgrid.dt) + 1;
idx_end = idx_start + length(base_signal) - 1;

if idx_end <= length(source_signal)
    source_signal(idx_start:idx_end) = base_signal;
end

transducer_transmit.input_signal = source_signal;


% -------------------------------------------------------------------------
% 6) トランスデューサーの初期化
% -------------------------------------------------------------------------
transducer_transmit = kWaveTransducer(kgrid, transducer_transmit);
transducer_receive = kWaveTransducer(kgrid, transducer_receive);
%transducer_transmit.properties;
%transducer_receive.properties;


% stream the data to disk in blocks of 100 if storing the complete time
% history 
%if ~USE_STATISTICS
%    input_args = [input_args {'StreamToDisk', 100}];
%end
%display_mask = transducer_transmit.active_elements_mask | transducer_receive.active_elements_mask | pipe_mask | glass_mask;
display_mask = transducer_transmit.active_elements_mask | transducer_receive.active_elements_mask | pipe_mask;
input_args = {'DisplayMask', display_mask, ...
    'RecordMovie', true, ...
    'MovieName', fullfile(save_path, 'liquid_only_d.avi'), ...
    'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], 'PMLAlpha', 2, ...
    'PlotScale', [-1/8, 1/8] * source_strength, ...
    'DataCast', DATA_CAST, ...
    };  


% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer_transmit, transducer_transmit, input_args{:});
%save(fullfile(save_path, 'liquid_only_ref.mat'), ...
%    'sensor_data', ...           % 送信用トランスデューサーで記録したデータ
%   'kgrid', ...                 % グリッド情報
%   'kgrid.t_array', ...               % 時間配列
%    'source_signal', ...         % 入力信号
%    '-v7.3');                    % 大きなデータセット用に-v7.3フォーマットを使用
% =========================================================================
% COMPUTE THE BEAM PATTERN USING SIMULATION STATISTICS
% =========================================================================

scan_line = transducer_transmit.scan_line(sensor_data);
plotsignalwaveform(kgrid, scan_line, save_path, 'signal_liquid_only_transmit.png');
% Method 3: Alternative visualization using isosurface
% % Plot the source signal
figure(2);
plot(kgrid.t_array * 1e6, source_signal * 1e-6, 'r-');
xlabel('Time [\mus]');
ylabel('Pressure [MPa]');
title('Source Signal');
grid on;
saveas(gcf, fullfile(save_path, 'source_signal.png'));
figure(3);
% Convert logical masks to double for visualization
transducer_transmit_mask_double = double(transducer_transmit.active_elements_mask);
transducer_receive_mask_double = double(transducer_receive.active_elements_mask);
pipe_mask_double = double(pipe_mask);

% Create separate masks for visualization
%glass_transducer_transmit_mask = glass_mask_double + transducer_transmit_mask_double + transducer_receive_mask_double;
transducer_transmit_mask = transducer_transmit_mask_double + transducer_receive_mask_double;
pipe_only_mask = pipe_mask_double;

% Create isosurface plots
% Glass and transducer_transmits (solid)
p1 = patch(isosurface(transducer_transmit_mask, 0.5));
isonormals(transducer_transmit_mask, p1);
set(p1, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.8);

% Pipe (transparent)
p2 = patch(isosurface(pipe_only_mask, 0.5));
isonormals(pipe_only_mask, p2);
set(p2, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% Set up the figure
daspect([1 1 1]);
view(80, 30);
camlight;
lighting gouraud;
axis([1 size(transducer_transmit_mask,1) 1 size(transducer_transmit_mask,2) 1 size(transducer_transmit_mask,3)]);
title('Liquid only experimental settings');
saveas(gcf, fullfile(save_path, 'liquid_only.png'));

% センサーデータの形状を確認
size(sensor_data)

% 必要に応じてリシェイプまたは並べ替え
if ~ismatrix(sensor_data)
    sensor_data = reshape(sensor_data, [], size(sensor_data, ndims(sensor_data)));
end

% stackedPlot を試す
try
    stackedPlot(kgrid.t_array * 1e6, sensor_data);
catch
    % エラーが発生した場合は基本的なプロットに切り替え
    figure;
    plot(kgrid.t_array * 1e6, sensor_data');
    xlabel('Time [µs]');
    ylabel('Amplitude');
end
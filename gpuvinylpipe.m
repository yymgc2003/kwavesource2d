% =========================================================================
% 2D k-Wave simulation with a focusing arc transducer and vinyl ring
% =========================================================================
clearvars;
close all;
DATA_CAST = 'gpuArray-single';
% -------------------------------------------------------------------------
% 1) シミュレーション用グリッドの定義
% -------------------------------------------------------------------------
Nx = 1024;               % x方向グリッド数 (行方向)
Ny = 1024;               % y方向グリッド数 (列方向)
dx = 0.1e-3;            % グリッド間隔 [m] (0.1 mm)
dy = 0.1e-3;            % グリッド間隔 [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);
save_path = '/mnt/sdb/matsubara/tmp'; %for dl-box
%save_path = '/mnt/matsubara/rawdata' ;% for jacob
% -------------------------------------------------------------------------
% 2) 媒質パラメータ
% -------------------------------------------------------------------------
% 水のパラメータ
medium.sound_speed = 1500;     % [m/s]
medium.density     = 1000;     % [kg/m^3]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;

% 塩ビのパラメータ
vinyl.sound_speed = 2390;      % [m/s] 塩ビの音速
vinyl.density     = 1400;      % [kg/m^3] 塩ビの密度

distance_pipe_source = 0.05; % [m] distance between glass and source
diameter = 0.009; % [m] effective transducer diameter
% -------------------------------------------------------------------------
% 3) ソースとガラス円環のマスクを作成
% -------------------------------------------------------------------------
source.p_mask = zeros(Nx, Ny);
source.p_mask(Nx/2-diameter/(2*dx):Nx/2+diameter/(2*dx), Ny/2-distance_pipe_source/dy) = 1; % locate at the center

% ガラス円環のマスクを作成
%glass_mask = zeros(Nx, Ny);
center_x = Nx/2+160;
center_y = Ny/2;
outer_radius = 160;  % 外側の半径
inner_radius = 130;   % 内側の半径
thickness = outer_radius - inner_radius;

% 円環のマスクを作成
[X, Y] = meshgrid(1:Nx, 1:Ny);
X = X - center_x;
Y = Y - center_y;
R = sqrt(X.^2 + Y.^2);
pipe_mask = (R <= outer_radius) & (R >= inner_radius);

% 媒質パラメータのマスクを作成
medium.sound_speed = medium.sound_speed * ones(Nx, Ny);
medium.density = medium.density * ones(Nx, Ny);

% 塩化ビニル円環のパラメータを設定
medium.sound_speed(pipe_mask == 1) = vinyl.sound_speed;
medium.density(pipe_mask == 1) = vinyl.density;

% -------------------------------------------------------------------------
% 4) シミュレーション時間配列の作成
% -------------------------------------------------------------------------
t_end = 1e-3; % longer than 1e-3 period simulation result in too long computational time. for movie visualization, 0.1ms simulation time is enough
kgrid.makeTime(medium.sound_speed, [], t_end);

% -------------------------------------------------------------------------
% 5) ソース波形の設定
% -------------------------------------------------------------------------
source_freq = 4e6;
source_mag = 1;
source_signal = zeros(size(kgrid.t_array));
prf = 3000;
T_prf = 1 / prf;
t_array = kgrid.t_array;
pulse_length = 1e-6;
max_n = 1000;

for n = 0:max_n
    t_start = n * T_prf;
    t_end = t_start + pulse_length;
    
    if t_start > kgrid.t_array(end)
        break;
    end
    
    idx_on = (t_array >= t_start) & (t_array < t_end);
    source_signal(idx_on) = source_mag * sin(2*pi * source_freq * (t_array(idx_on) - t_start));
end
source.p = source_signal;

% -------------------------------------------------------------------------
% 6) センサーの設定
% -------------------------------------------------------------------------
sensor.mask = zeros(Nx, Ny);
sensor_x = Nx/2;
sensor_y = Ny/2 + 400;
sensor.mask(sensor_x, sensor_y) = 1;
sensor.record = {'p'};

% -------------------------------------------------------------------------
% 7) シミュレーションのオプション設定
% -------------------------------------------------------------------------
input_args = {
    'PMLInside', false, 'PlotPML', false, ...
    'DataCast', DATA_CAST, ...
    };

% -------------------------------------------------------------------------
% 8) シミュレーション実行
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});

% -------------------------------------------------------------------------
% 9) 結果の可視化
% -------------------------------------------------------------------------
figure;
plot(kgrid.t_array*1e3, sensor_data.p(1, :));
xlabel('Time [ms]');
ylabel('Pressure [Pa]');
title('Pressure at the sensor with vinyl pipe');
saveas(gcf, fullfile(save_path, 'sensor_vinyl_pipe_bad.png')); 
% kspaceFirstOrder2D実行後にGPU配列をCPUに集約
sensor_data_cpu = structfun(@gather, sensor_data, 'UniformOutput', false);

% v7.3 形式で.matファイル保存
save(fullfile(save_path, 'sensor_data_pipe.mat'), 'sensor_data_cpu', '-v7.3');
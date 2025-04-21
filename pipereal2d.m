% =========================================================================
% 2D k-Wave simulation with a focusing arc transducer, vinyl ring and glass beads
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
save_path = '/mnt/sdb/matsubara/tmp';

% -------------------------------------------------------------------------
% 2) 媒質パラメータ
% -------------------------------------------------------------------------
% 水のパラメータ
medium.sound_speed = 1500;     % [m/s]
medium.density     = 1000;     % [kg/m^3]
%medium.alpha_coeff = 0;        % dB/(MHz^y cm)
%medium.alpha_power = 1.0;
%medium.alpha_mode  = 'no_dispersion';

% ガラスのパラメータ
vinyl.sound_speed = 2390;      % [m/s] 塩ビの音速
vinyl.density     = 1400;      % [kg/m^3] 塩ビの密度

% ガラスビーズのパラメータ
bead.sound_speed = 5500;       % [m/s] ガラス中の音速
bead.density     = 2500;       % [kg/m^3] ガラスの密度
bead.diameter    = 5e-3;       % [m] ガラスビーズの直径 (5mm)
bead.radius_grid = round(bead.diameter/(2*dx)); % グリッド単位での半径

distance_pipe_source = 0.05; % [m] distance between glass and source
diameter = 0.009; % [m] effective transducer diameter
% -------------------------------------------------------------------------
% 3) ソースとガラス円環のマスクを作成
% -------------------------------------------------------------------------
source.p_mask = zeros(Nx, Ny);
source.p_mask(Nx/2-diameter/(2*dx):Nx/2+diameter/(2*dx), Ny/2-distance_pipe_source/dy) = 1; % locate at the center

% ガラス円環のマスクを作成
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

% ガラスビーズのマスクを作成
bead_mask = zeros(Nx, Ny);
num_beads = 3; % ガラスビーズの数

% 円管内部にランダムにガラスビーズを配置
for i = 1:num_beads
    % 円管内部の範囲を定義（円環の内側）
    % 円管内のランダムな位置を生成
    theta = 2 * pi * rand(); % ランダムな角度
    r = inner_radius * sqrt(rand()); % ランダムな半径（均一分布のため平方根を取る）
    
    % 円管内の座標を計算
    bead_x = center_x + r * cos(theta);
    bead_y = center_y + r * sin(theta);
    
    % ガラスビーズのマスクを作成（円形）
    bead_distance = sqrt((X + center_x - bead_x).^2 + (Y + center_y - bead_y).^2);
    bead_mask(bead_distance <= bead.radius_grid) = 1;
end

% 媒質パラメータのマスクを作成
medium.sound_speed = medium.sound_speed * ones(Nx, Ny);
medium.density = medium.density * ones(Nx, Ny);

% ガラス円環のパラメータを設定
medium.sound_speed(pipe_mask == 1) = vinyl.sound_speed;
medium.density(pipe_mask == 1) = vinyl.density;

% ガラスビーズのパラメータを設定
medium.sound_speed(bead_mask == 1) = bead.sound_speed;
medium.density(bead_mask == 1) = bead.density;

% -------------------------------------------------------------------------
% 4) シミュレーション時間配列の作成
% -------------------------------------------------------------------------
t_end = 1e-4; % longer than 1e-3 period simulation result in too long computational time. for movie visualization, 0.1ms simulation time is enough
kgrid.makeTime(medium.sound_speed, [], t_end);

% -------------------------------------------------------------------------
% 5) ソース波形の設定
% -------------------------------------------------------------------------
source_freq = 4e6;
source_mag = 10;
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
    source_signal(idx_on) = sin(2*pi * source_freq * (t_array(idx_on) - t_start));
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
    'RecordMovie', true, ...
    'MovieName', fullfile(save_path, 'pipereal2d_with_beads.avi'), ...
    'DataCast', DATA_CAST, ...
    };

% -------------------------------------------------------------------------
% 8) シミュレーション実行
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% -------------------------------------------------------------------------
% 9) 結果の可視化
% -------------------------------------------------------------------------
figure;
plot(kgrid.t_array*1e6, sensor_data.p(1, :));
xlabel('Time [\mus]');
ylabel('Pressure [Pa]');
title('Pressure at the sensor with vinyl pipe and glass beads');
saveas(gcf, fullfile(save_path, 'sensor_pipe_with_beads.png')); 
% kspaceFirstOrder2D実行後にGPU配列をCPUに集約
sensor_data_cpu = structfun(@gather, sensor_data, 'UniformOutput', false);

% v7.3 形式で.matファイル保存
save(fullfile(save_path, 'sensor_data_pipe_with_beads.mat'), 'sensor_data_cpu', '-v7.3'); 
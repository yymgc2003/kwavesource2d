% =========================================================================
% 2D k-Wave simulation with a focusing arc transducer and glass ring
% =========================================================================
clearvars;
close all;

% -------------------------------------------------------------------------
% 1) シミュレーション用グリッドの定義
% -------------------------------------------------------------------------
Nx = 256;               % x方向グリッド数 (行方向)
Ny = 256;               % y方向グリッド数 (列方向)
dx = 0.1e-3;            % グリッド間隔 [m] (0.1 mm)
dy = 0.1e-3;            % グリッド間隔 [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);
save_path = '/home/matsubara/Scripts/tmp';

% -------------------------------------------------------------------------
% 2) 媒質パラメータ
% -------------------------------------------------------------------------
% 水のパラメータ
medium.sound_speed = 1500;     % [m/s]
medium.density     = 1000;     % [kg/m^3]
medium.alpha_coeff = 0;        % dB/(MHz^y cm)
%medium.alpha_power = 1.0;
medium.alpha_mode  = 'no_dispersion';

% ガラスのパラメータ
glass.sound_speed = 5500;      % [m/s] ガラスの音速
glass.density     = 2500;      % [kg/m^3] ガラスの密度

% -------------------------------------------------------------------------
% 3) ソースとガラス円環のマスクを作成
% -------------------------------------------------------------------------
source.p_mask = zeros(Nx, Ny);
source.p_mask(50:200, Ny/2) = 1;

% ガラス円環のマスクを作成
%glass_mask = zeros(Nx, Ny);
center_x = Nx/2;
center_y = Ny/2;
outer_radius = 100;  % 外側の半径
inner_radius = 90;   % 内側の半径
thickness = outer_radius - inner_radius;

% 円環のマスクを作成
[X, Y] = meshgrid(1:Nx, 1:Ny);
X = X - center_x;
Y = Y - center_y;
R = sqrt(X.^2 + Y.^2);
glass_mask = (R <= outer_radius) & (R >= inner_radius);

% 媒質パラメータのマスクを作成
medium.sound_speed_map = medium.sound_speed * ones(Nx, Ny);
medium.density_map = medium.density * ones(Nx, Ny);

% ガラス円環のパラメータを設定
medium.sound_speed(glass_mask == 1) = glass.sound_speed;
medium.density(glass_mask == 1) = glass.density;

% -------------------------------------------------------------------------
% 4) シミュレーション時間配列の作成
% -------------------------------------------------------------------------
t_end = 1e-3;
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
    source_signal(idx_on) = sin(2*pi * source_freq * (t_array(idx_on) - t_start));
end
source.p = source_signal;

% -------------------------------------------------------------------------
% 6) センサーの設定
% -------------------------------------------------------------------------
sensor.mask = zeros(Nx, Ny);
sensor_x = Nx/2;
sensor_y = Ny/2 + 70;
sensor.mask(sensor_x, sensor_y) = 1;
sensor.record = {'p'};

% -------------------------------------------------------------------------
% 7) シミュレーションのオプション設定
% -------------------------------------------------------------------------
input_args = {
    'DataCast', 'gpuArray-double', ...
};

% -------------------------------------------------------------------------
% 8) シミュレーション実行
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});

% -------------------------------------------------------------------------
% 9) 結果の可視化
% -------------------------------------------------------------------------
figure;
plot(kgrid.t_array*1e6, sensor_data.p(1, :));
xlabel('Time [\mus]');
ylabel('Pressure [Pa]');
title('Pressure at the sensor with glass ring');
saveas(gcf, fullfile(save_path, 'sensor_glass_ring.png')); 
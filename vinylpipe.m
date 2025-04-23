% =========================================================================
% 2D k-Wave simulation with a focusing arc transducer and vinyl pipe
% =========================================================================
clearvars;
close all;
DATA_CAST='gpuArray-single';
% -------------------------------------------------------------------------
% 1) 設定ファイルの読み込み
% -------------------------------------------------------------------------
config = jsondecode(fileread('config.json'));

% -------------------------------------------------------------------------
% 2) シミュレーション用グリッドの定義
% -------------------------------------------------------------------------
kgrid = kWaveGrid(config.grid.Nx, config.grid.dx, config.grid.Ny, config.grid.dy);
save_path = config.save_path;

% -------------------------------------------------------------------------
% 3) 媒質パラメータ
% -------------------------------------------------------------------------
% 水のパラメータ
medium.sound_speed = config.medium.water.sound_speed;
medium.density     = config.medium.water.density;
medium.alpha_coeff = config.medium.water.alpha_coeff;
medium.alpha_power = config.medium.water.alpha_power;
medium.BonA = config.medium.water.BonA;

% ガラスのパラメータ
vinyl.sound_speed = config.medium.vinyl.sound_speed;
vinyl.density     = config.medium.vinyl.density;

% -------------------------------------------------------------------------
% 4) トランスデューサーの設定
% -------------------------------------------------------------------------
% トランスデューサーの位置と向きを設定
transducer_pos = config.transducer.position;
transducer_rot = config.transducer.rotation;
transducer_focus = config.transducer.focus;

% ガラス円環のマスクを作成
%glass_mask = zeros(Nx, Ny);
center_x = config.grid.Nx/2+160;
center_y = config.grid.Ny/2;
outer_radius = 160;  % 外側の半径
inner_radius = 130;   % 内側の半径
thickness = outer_radius - inner_radius;

% -------------------------------------------------------------------------
% 5) 塩ビパイプの設定
% -------------------------------------------------------------------------
% パイプの中心位置
center_x = config.grid.Nx/2 + config.pipe.center_x_offset;
center_y = config.grid.Ny/2 + config.pipe.center_y_offset;

% パイプのマスクを作成
[X, Y] = meshgrid(1:config.grid.Nx, 1:config.grid.Ny);
X = X - center_x;
Y = Y - center_y;
R = sqrt(X.^2 + Y.^2);
pipe_mask = (R <= config.pipe.outer_radius) & (R >= config.pipe.inner_radius);

% 媒質パラメータのマスクを作成
medium.sound_speed = medium.sound_speed * ones(config.grid.Nx, config.grid.Ny);
medium.density = medium.density * ones(config.grid.Nx, config.grid.Ny);

% 塩ビパイプのパラメータを設定
medium.sound_speed(pipe_mask == 1) = vinyl.sound_speed;
medium.density(pipe_mask == 1) = vinyl.density;

% -------------------------------------------------------------------------
% 6) シミュレーション時間配列の作成
% -------------------------------------------------------------------------
kgrid.makeTime(medium.sound_speed, 0.04, config.simulation.t_end);

% -------------------------------------------------------------------------
% 7) ソース波形の設定
% -------------------------------------------------------------------------
distance_pipe_source = config.source.distance_pipe_source;
diameter = config.source.diameter;
source.p_mask = zeros(config.grid.Nx, config.grid.Ny);
source.p_mask(config.grid.Nx/2-diameter/(2*config.grid.dx):config.grid.Nx/2+diameter/(2*config.grid.dx), config.grid.Ny/2-distance_pipe_source/config.grid.dy) = 1;
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
source.p = source_signal;

% -------------------------------------------------------------------------
% 8) センサーの設定
% -------------------------------------------------------------------------
sensor.mask = zeros(config.grid.Nx, config.grid.Ny);
sensor_x = config.grid.Nx/2;
sensor_y = config.grid.Ny/2 + 400;
sensor.mask(sensor_x, sensor_y) = 1;
sensor.record = {'p'};

% -------------------------------------------------------------------------
% 9) シミュレーションのオプション設定
% -------------------------------------------------------------------------
input_args = {
    'PMLInside', false, 'PlotPML', false, ...
    'RecordMovie', true, ...
    'MovieName', fullfile(save_path, 'vinylpipe_exp.avi'), ...
    'DataCast', DATA_CAST, ...
    };

% -------------------------------------------------------------------------
% 10) シミュレーション実行
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% -------------------------------------------------------------------------
% 11) 結果の可視化
% -------------------------------------------------------------------------
figure;
plot(kgrid.t_array*1e3, sensor_data.p(1, :));
xlabel('Time [ms]');
ylabel('Pressure [Pa]');
title('Pressure at the sensor with vinyl pipe');
saveas(gcf, fullfile(save_path, 'sensor_vinyl_pipe.png')); 

% v7.3 形式で.matファイル保存
save(fullfile(save_path, 'sensor_data_pipe.mat'), 'sensor_data', '-v7.3');

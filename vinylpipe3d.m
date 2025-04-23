% =========================================================================
% 3D k-Wave simulation with a focusing arc transducer and vinyl pipe
% =========================================================================
clearvars;
close all;

% -------------------------------------------------------------------------
% 1) 設定ファイルの読み込み
% -------------------------------------------------------------------------
config = jsondecode(fileread('config.json'));

% -------------------------------------------------------------------------
% 2) シミュレーション用グリッドの定義
% -------------------------------------------------------------------------
kgrid = kWaveGrid(config.grid.Nx, config.grid.dx, config.grid.Ny, config.grid.dy, config.grid.Nz, config.grid.dz);
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

% トランスデューサーのマスクを作成
source.p_mask = zeros(config.grid.Nx, config.grid.Ny, config.grid.Nz);
% ここにトランスデューサーの形状に基づいたマスク作成コードを追加

% -------------------------------------------------------------------------
% 5) ガラスパイプの設定
% -------------------------------------------------------------------------
% パイプの中心位置
center_x = config.grid.Nx/2 + config.pipe.center_x_offset;
center_y = config.grid.Ny/2 + config.pipe.center_y_offset;
center_z = config.grid.Nz/2 + config.pipe.center_z_offset;

% パイプのマスクを作成
[X, Y, Z] = meshgrid(1:config.grid.Nx, 1:config.grid.Ny, 1:config.grid.Nz);
X = X - center_x;
Y = Y - center_y;
Z = Z - center_z;
R = sqrt(X.^2 + Y.^2);
pipe_mask = (R <= config.pipe.outer_radius) & (R >= config.pipe.inner_radius) & ...
            (abs(Z) <= config.pipe.length/(2*config.grid.dz));

% 媒質パラメータのマスクを作成
medium.sound_speed = medium.sound_speed * ones(config.grid.Nx, config.grid.Ny, config.grid.Nz);
medium.density = medium.density * ones(config.grid.Nx, config.grid.Ny, config.grid.Nz);

% ガラスパイプのパラメータを設定
medium.sound_speed(pipe_mask == 1) = vinyl.sound_speed;
medium.density(pipe_mask == 1) = vinyl.density;

% -------------------------------------------------------------------------
% 6) シミュレーション時間配列の作成
% -------------------------------------------------------------------------
kgrid.makeTime(medium.sound_speed, [], config.simulation.t_end);

% -------------------------------------------------------------------------
% 7) ソース波形の設定
% -------------------------------------------------------------------------
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
sensor.mask = zeros(config.grid.Nx, config.grid.Ny, config.grid.Nz);
sensor_x = config.grid.Nx/2 + config.sensor.x_offset;
sensor_y = config.grid.Ny/2 + config.sensor.y_offset;
sensor_z = config.grid.Nz/2 + config.sensor.z_offset;

if strcmp(config.sensor.type, 'point')
    sensor.mask(sensor_x, sensor_y, sensor_z) = 1;
elseif strcmp(config.sensor.type, 'array')
    % アレイセンサーの設定
    % ここにアレイセンサーの設定コードを追加
end
sensor.record = {'p'};

% -------------------------------------------------------------------------
% 9) シミュレーションのオプション設定
% -------------------------------------------------------------------------
input_args = {
    'PMLInside', false, ...
    'PlotPML', config.visualization.plot_pml, ...
    'DataCast', config.simulation.data_cast, ...
    'PMLSize', config.simulation.pml_size, ...
    'PMLAlpha', config.simulation.pml_alpha, ...
    'PlotSim', config.visualization.plot_sim, ...
    'PlotScale', config.visualization.plot_scale, ...
    'RecordMovie', config.visualization.record_movie, ...
    'MovieName', config.visualization.movie_name
    };

% -------------------------------------------------------------------------
% 10) シミュレーション実行
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% -------------------------------------------------------------------------
% 11) 結果の可視化
% -------------------------------------------------------------------------
if config.visualization.plot_sim
    figure;
    plot(kgrid.t_array*1e3, sensor_data.p(1, :));
    xlabel('Time [ms]');
    ylabel('Pressure [Pa]');
    title('Pressure at the sensor with vinyl pipe');
    saveas(gcf, fullfile(save_path, 'sensor_vinyl_pipe_3d.png')); 
end

% kspaceFirstOrder3D実行後にGPU配列をCPUに集約
sensor_data_cpu = structfun(@gather, sensor_data, 'UniformOutput', false);

% v7.3 形式で.matファイル保存
save(fullfile(save_path, 'sensor_data_pipe_3d.mat'), 'sensor_data_cpu', '-v7.3');
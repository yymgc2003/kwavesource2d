% =========================================================================
% 2D k-Wave simulation with a focusing arc transducer and vinyl pipe
% =========================================================================
clearvars;
close all;
DATA_CAST='gpuArray-single';
% -------------------------------------------------------------------------
% 1) 設定ファイルの読み込み
% -------------------------------------------------------------------------
config = jsondecode(fileread('/home/user01/Document/yyamaguchi/documents/kwavesource2d/config.json'));

% -------------------------------------------------------------------------
% 2) シミュレーション用グリッドの定義
% -------------------------------------------------------------------------
Nx = config.grid.Nx;
Ny = config.grid.Ny;
dx = config.grid.dx;
dy = config.grid.dy;
kgrid = kWaveGrid(Nx, dx, Ny, dy);
save_path = config.save_path;
t_end = config.simulation.t_end
CFL = config.simulation.CFL

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
vinyl.alpha_coeff = config.medium.vinyl.alpha_coeff;

% -------------------------------------------------------------------------
% 4) トランスデューサーの設定
% -------------------------------------------------------------------------
% トランスデューサーの位置と向きを設定
transducer_pos = config.transducer.position;
transducer_rot = config.transducer.rotation;
transducer_focus = config.transducer.focus;

% ガラス円環のマスクを作成
%glass_mask = zeros(Nx, Ny);
outer_radius = config.pipe.outer_radius;  % 外側の半径
inner_radius = config.pipe.inner_radius;   % 内側の半径
thickness = outer_radius - inner_radius;

% -------------------------------------------------------------------------
% 5) 塩ビパイプの設定
% -------------------------------------------------------------------------
% パイプの中心位置
center_x = Nx/2;
center_y = 30 + round(config.pipe.center_y/dy);

% 媒質パラメータのマスクを作成
medium.sound_speed = medium.sound_speed * ones(Nx, Ny);
medium.density = medium.density * ones(Nx, Ny);
medium.alpha_coeff = medium.alpha_coeff * ones(Nx, Ny);
medium.BonA = medium.BonA * ones(Nx, Ny);

% 塩ビパイプのパラメータを設定
for i = 1:Nx
    for j = 1:Ny
        if dx*sqrt((i-Nx/2)^2+(j-center_y)^2)<= outer_radius & dx*sqrt((i-Nx/2)^2+(j-center_y)^2)>= inner_radius
            medium.sound_speed(i, j) = vinyl.sound_speed;
            medium.density(i, j) = vinyl.density;
            medium.alpha_coeff(i, j) = vinyl.alpha_coeff;
        end
    end
end

bubble_r = round(0.75e-3 / dx);

for i = Nx/2-bubble_r-1:Nx/2+bubble_r+1
    for j = center_y-bubble_r-1:center_y+bubble_r+1
        if sqrt((i-Nx/2)^2+(j-center_y)^2) <= bubble_r
            medium.sound_speed(i, j) = config.medium.air.sound_speed;
            medium.density(i, j) = config.medium.air.density;
            medium.alpha_coeff(i, j) = config.medium.air.alpha_coeff;
        end
    end
end

% -------------------------------------------------------------------------
% 6) シミュレーション時間配列の作成
% -------------------------------------------------------------------------
kgrid.makeTime(config.medium.water.sound_speed, CFL, t_end);

% -------------------------------------------------------------------------
% 7) ソース波形の設定
% -------------------------------------------------------------------------
distance_pipe_source = config.source.distance_pipe_source;
distance_pipe_source = distance_pipe_source/dy;
source_diameter = config.source.diameter;
source_diameter = source_diameter/dx;
source.p_mask = zeros(Nx, Ny);
source.p_mask(round(Nx/2-source_diameter/2):round(Nx/2+source_diameter/2), 30)= 1;
source_signal = zeros(size(kgrid.t_array));
source_frequency = config.source.tone_burst_freq;
t_array = kgrid.t_array;

% Input signal properties
source_strength = config.source.source_strength;
tone_burst_freq = config.source.tone_burst_freq;
tone_burst_cycles = config.source.tone_burst_cycles;

% Generate input signal
% Read upsampled pulse from mat file
% original_pulse = load('/home/matsubara/Scripts/kwavesource/src/pulse_4000_upsampled.mat');
%input_signal = original_pulse.pulse_upsampled;

input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% Scale by impedance
input_signal = (source_strength / (config.medium.water.sound_speed * config.medium.water.density)) * input_signal(1:end-10);
source.p = input_signal;

% -------------------------------------------------------------------------
% 8) センサーの設定
% -------------------------------------------------------------------------
sensor.mask = zeros(Nx, Ny);
NNx = Nx/10;
NNy = Ny/10;
for i = 1:NNx
    for j = 1:NNy
        sensor.mask(-Nx/NNx/2+Nx/NNx*i, -Ny/NNy/2+Ny/NNy*j ) = 1;
    end
end
sensor.record = {'p'};

% -------------------------------------------------------------------------
% 9) シミュレーションのオプション設定
% -------------------------------------------------------------------------
input_args = {
    'PMLInside', true, 'PlotPML', true, ...
    'PMLSize', config.simulation.pml_size, ...
    'RecordMovie', false, ...
    'PlotFreq', 50, ...
    'MovieName', fullfile(save_path, 'vinylpipe_exp.avi'), ...
    'DataCast', DATA_CAST, ...
    };

% -------------------------------------------------------------------------
% 10) シミュレーション実行
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});

% -------------------------------------------------------------------------
% 11) 結果の可視化
% -------------------------------------------------------------------------

meshy = 0.0;

for t = 1:100
    file_id = fopen(fullfile(save_path, sprintf('vinylpipe_map%d.xy', t)), 'w');
    for i = 1:NNy
        for j = 1:NNx
            meshy = (sensor_data.p(NNx*(i-1)+j, 700*(t-1)+1));
            fprintf(file_id, '%.2f %.2f %.2e\n', 0.21+0.42*(i-1), 0.21+0.42*(j-1), meshy);
        end
        fprintf(file_id, '\n');
    end
    fclose(file_id);
end
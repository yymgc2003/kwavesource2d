% =========================================================================
% 3D k-Wave simulation with a focusing transducer and vinyl ring (GPU version)
% =========================================================================
clearvars;
close all;
DATA_CAST = 'gpuArray-single';

% -------------------------------------------------------------------------
% 1) シミュレーション用グリッドの定義
% -------------------------------------------------------------------------
Nx = 256;               % x方向グリッド数
Ny = 256;               % y方向グリッド数
Nz = 256;               % z方向グリッド数
dx = 0.2e-3;            % グリッド間隔 [m] (0.2 mm)
dy = 0.2e-3;            % グリッド間隔 [m]
dz = 0.2e-3;            % グリッド間隔 [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

save_path = '/mnt/sdb/matsubara/tmp'; %for dl-box
%save_path = '/mnt/matsubara/rawdata' ;% for jacob

% -------------------------------------------------------------------------
% 2) 媒質パラメータ
% -------------------------------------------------------------------------
% 水のパラメータ
medium.sound_speed = 1500;     % [m/s]
medium.density     = 1000;     % [kg/m^3]

% ビニールのパラメータ
vinyl.sound_speed = 2390;      % [m/s]
vinyl.density     = 1400;      % [kg/m^3]

% トランスデューサーのパラメータ
transducer_freq = 4e6;         % トランスデューサー周波数 [Hz]
focus_distance = 0.05;         % 焦点距離 [m]
diameter = 0.009;             % トランスデューサー直径 [m]

% -------------------------------------------------------------------------
% 3) トランスデューサーとビニール円環の設定
% -------------------------------------------------------------------------
% トランスデューサーの作成
bowl_pos = [kgrid.x_vec(end/2), kgrid.y_vec(end/2), kgrid.z_vec(end/2) - focus_distance];
focus_pos = [kgrid.x_vec(end/2), kgrid.y_vec(end/2), kgrid.z_vec(end/2)];
source = kWaveTransducer(kgrid, medium, bowl_pos, focus_pos, diameter);
source.frequency = transducer_freq;

% ビニール円環のマスクを作成
center = [Nx/2, Ny/2, Nz/2];
outer_radius = 80;  % 外側の半径（グリッドポイント）
inner_radius = 65;  % 内側の半径（グリッドポイント）

[X, Y, Z] = meshgrid(1:Nx, 1:Ny, 1:Nz);
X = X - center(1);
Y = Y - center(2);
Z = Z - center(3);
R = sqrt(X.^2 + Y.^2);
pipe_mask = (R <= outer_radius) & (R >= inner_radius);

% 媒質パラメータの設定
medium.sound_speed = medium.sound_speed * ones(Nx, Ny, Nz, DATA_CAST);
medium.density = medium.density * ones(Nx, Ny, Nz, DATA_CAST);

% ビニール円環のパラメータを設定
medium.sound_speed(pipe_mask) = vinyl.sound_speed;
medium.density(pipe_mask) = vinyl.density;

% -------------------------------------------------------------------------
% 4) シミュレーション時間の設定
% -------------------------------------------------------------------------
t_end = 40e-6;
kgrid.makeTime(medium.sound_speed);

% -------------------------------------------------------------------------
% 5) センサーの設定
% -------------------------------------------------------------------------
sensor.mask = zeros(Nx, Ny, Nz);
sensor_plane = zeros(Nx, Ny);
sensor_plane(Nx/2-50:Nx/2+50, Ny/2-50:Ny/2+50) = 1;
sensor.mask(:, :, Nz/2) = sensor_plane;
sensor.record = {'p'};

% -------------------------------------------------------------------------
% 6) シミュレーションのオプション設定
% -------------------------------------------------------------------------
input_args = {...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DataCast', DATA_CAST ...
    };

% -------------------------------------------------------------------------
% 7) トランスデューサーの可視化
% -------------------------------------------------------------------------
figure;
voxelPlot(pipe_mask);
hold on;
source.plotTransducer();
title('Transducer and Vinyl Pipe Configuration');
xlabel('X [grid points]');
ylabel('Y [grid points]');
zlabel('Z [grid points]');
view(45, 30);
saveas(gcf, fullfile(save_path, 'transducer_config_3d_gpu.png'));

% -------------------------------------------------------------------------
% 8) シミュレーション実行
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% -------------------------------------------------------------------------
% 9) 結果の可視化と保存
% -------------------------------------------------------------------------
% GPU配列をCPUに移動
sensor_data_cpu = structfun(@gather, sensor_data, 'UniformOutput', false);

figure;
imagesc(sensor_data_cpu.p(:, :, end/2));
colorbar;
title('Pressure Field at Central Plane');
xlabel('X [grid points]');
ylabel('Y [grid points]');
saveas(gcf, fullfile(save_path, 'pressure_field_3d_gpu.png'));

% データの保存
save(fullfile(save_path, 'sensor_data_3d_gpu.mat'), 'sensor_data_cpu', '-v7.3'); 
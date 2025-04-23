% =========================================================================
% 3D k-Wave simulation with a focusing arc transducer and vinyl pipe (GPU version)
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
kgrid.makeTime(medium.sound_speed, 0.04, config.simulation.t_end);

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
    
% using transducer1 as source 
transducer1.input_signal = input_signal;
% calculate the width of the transducer1 in grid points
transducer_width = transducer2.number_elements * transducer2.element_width ...
    + (transducer2.number_elements - 1) * transducer2.element_spacing;
transducer1.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer1.element_length/2]);
transducer2.position = round([Nx-10, Ny/2 - transducer_width/2, Nz/2 - transducer1.element_length/2]);

transducer1.active_elements = zeros(transducer1.number_elements, 1);
transducer1.active_elements(4:30) = 1;
transducer2.active_elements = zeros(transducer2.number_elements, 1);
transducer2.active_elements(4:30) = 1;

transducer1 = kWaveTransducer(kgrid, transducer1);          
transducer2 = kWaveTransducer(kgrid, transducer2);
% print out transducer1 properties
transducer1.properties;
transducer2.properties;


% ビニール円環のマスクを作成 - サイズを調整
cx = Nx/2;            % X 方向の中心
cy = Ny/2 + 40;       % Y 方向の中心
outer_r = 160; inner_r = 116;

[Xg, Yg] = ndgrid(1:Nx, 1:Ny);
ring2d = sqrt( (Xg-cx).^2 + (Yg-cy).^2 );
ringMask = (ring2d <= outer_r) & (ring2d >= inner_r);   % Nx×Ny logical

pipe_mask = repmat(ringMask, [1 1 Nz]);   % Nx×Ny×Nz     % Nx×Ny×Nz

% Initialize medium parameters
medium.sound_speed = medium.sound_speed .* ones(Nx, Ny, Nz);
medium.density = medium.density .* ones(Nx, Ny, Nz);

medium.sound_speed(pipe_mask) = vinyl.sound_speed;
medium.density(pipe_mask) = vinyl.density;


% -------------------------------------------------------------------------
% 8) センサーの設定
% -------------------------------------------------------------------------
sensor.mask = zeros(Nx, Ny, Nz);
sensor_plane = zeros(Nx, Ny);
sensor_plane(Nx/2-50:Nx/2+50, Ny/2-50:Ny/2+50) = 1;
sensor.mask(:, :, Nz/2) = sensor_plane;
sensor.record = {'p'};

% -------------------------------------------------------------------------
% 9) シミュレーションのオプション設定
% -------------------------------------------------------------------------
display_mask = transducer1.active_elements_mask | transducer2.active_elements_mask | pipe_mask;
input_args = {...
    'DisplayMask', display_mask ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'RecordMovie', true, ...
    'MovieName', fullfile(save_path, 'vinyl_pipe_3d.avi'), ...
    'DataCast', DATA_CAST
    };

% -------------------------------------------------------------------------
% 7) トランスデューサーの可視化
% -------------------------------------------------------------------------


% create voxel plot of transducer1 mask and 
voxelPlot(single(transducer1.active_elements_mask | transducer2.active_elements_mask));
view(126, 25);
saveas(gcf, fullfile(save_path, 'trans_config_3d.png'));


% -------------------------------------------------------------------------
% 11) 結果の可視化
% -------------------------------------------------------------------------
figure;
plot(kgrid.t_array*1e3, sensor_data.p(1, :));
xlabel('Time [ms]');
ylabel('Pressure [Pa]');
title('Pressure at the sensor with vinyl pipe (GPU)');
saveas(gcf, fullfile(save_path, 'sensor_vinyl_pipe_3d_gpu.png')); 

% kspaceFirstOrder3DG実行後にGPU配列をCPUに集約
sensor_data_cpu = structfun(@gather, sensor_data, 'UniformOutput', false);

% v7.3 形式で.matファイル保存
save(fullfile(save_path, 'sensor_data_pipe_3d_gpu.mat'), 'sensor_data_cpu', '-v7.3');
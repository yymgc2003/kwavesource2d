% =========================================================================
% 3D k-Wave simulation with a focusing arc transducer and vinyl pipe
% =========================================================================
clearvars;
close all;
DATA_CAST = 'gpuArray-single';

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

% 塩ビ管のパラメータ
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

% 塩ビ管のパラメータを設定
medium.sound_speed(pipe_mask == 1) = vinyl.sound_speed;
medium.density(pipe_mask == 1) = vinyl.density;

% -------------------------------------------------------------------------
% 6) シミュレーション時間配列の作成
% -------------------------------------------------------------------------
kgrid.makeTime(medium.sound_speed, [], config.simulation.t_end);
% create the input signal using toneBurst 
source_strength = 1e2;          % [Pa]
tone_burst_freq = 4e6;    	% [Hz]
tone_burst_cycles = 4;
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength ./ (medium.sound_speed(1) * medium.density(1))) .* input_signal;
transducer1.number_elements = config.transducer.elements;
transducer2.number_elements = config.transducer.elements;
transducer2.element_width = config.transducer.element_width;
transducer2.element_length = config.transducer.element_height;
transducer2.element_spacing = config.transducer.element_spacing;
transducer1.radius = inf;
transducer2.radius = inf;
% using transducer1 as source 
transducer1.input_signal = input_signal;
% calculate the width of the transducer1 in grid points
transducer_width = transducer2.number_elements * transducer2.element_width ...
    + (transducer2.number_elements - 1) * transducer2.element_spacing;
transducer1.position = round([1, config.grid.Ny/2 - transducer_width/2, config.grid.Nz/2 - transducer2.element_length/2]);
transducer2.position = round([config.grid.Nx-10, config.grid.Ny/2 - transducer_width/2, config.grid.Nz/2 - transducer2.element_length/2]);

transducer1.active_elements = zeros(transducer1.number_elements, 1);
transducer1.active_elements(4:30) = 1;
transducer2.active_elements = zeros(transducer2.number_elements, 1);
transducer2.active_elements(4:30) = 1;

transducer1 = kWaveTransducer(kgrid, transducer1);          
transducer2 = kWaveTransducer(kgrid, transducer2);
% print out transducer1 properties
%transducer1.properties;
%transducer2.properties;


% ビニール円環のマスクを作成 - サイズを調整
cx = config.grid.Nx/2;            % X 方向の中心
cy = config.grid.Ny/2;       % Y 方向の中心
outer_r = 160; inner_r = 116;

[Xg, Yg] = ndgrid(1:config.grid.Nx, 1:config.grid.Ny);
ring2d = sqrt( (Xg-cx).^2 + (Yg-cy).^2 );
ringMask = (ring2d <= outer_r) & (ring2d >= inner_r);   % Nx×Ny logical

pipe_mask = repmat(ringMask, [1 1 config.grid.Nz]);   % Nx×Ny×Nz   

% Initialize medium parameters
medium.sound_speed = medium.sound_speed .* ones(config.grid.Nx, config.grid.Ny, config.grid.Nz);
medium.density = medium.density .* ones(config.grid.Nx, config.grid.Ny, config.grid.Nz);

medium.sound_speed(pipe_mask) = vinyl.sound_speed;
medium.density(pipe_mask) = vinyl.density;


% -------------------------------------------------------------------------
% 8) センサーの設定
% -------------------------------------------------------------------------
sensor.mask = zeros(config.grid.Nx, config.grid.Ny, config.grid.Nz);
sensor_plane = zeros(config.grid.Nx, config.grid.Ny);
sensor_plane(config.grid.Nx/2-50:config.grid.Nx/2+50, config.grid.Ny/2-50:config.grid.Ny/2+50) = 1;
sensor.mask(:, :, config.grid.Nz/2) = sensor_plane;
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
% 10) シミュレーション実行
% -------------------------------------------------------------------------


% create voxel plot of transducer1 mask and 
voxelPlot(single(transducer1.active_elements_mask | transducer2.active_elements_mask | pipe_mask));
view(126, 25);
saveas(gcf, fullfile(save_path, 'trans_config_3d.png'));


% -------------------------------------------------------------------------
% 11) 結果の可視化
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer1, transducer2, input_args{:});

% -------------------------------------------------------------------------
% 9) 結果の可視化
% -------------------------------------------------------------------------
field=gather(sensor_data);
figure;
imagesc(field(:,:));
colorbar;
title('Pressure Field at Central Plane');
xlabel('X [grid points]');
ylabel('Y [grid points]');
saveas(gcf, fullfile(save_path, 'pressure_field_3d.png'));

% データの保存
save(fullfile(save_path, 'sensor_data_3d.mat'), 'sensor_data', '-v7.3');
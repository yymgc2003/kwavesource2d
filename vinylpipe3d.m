% =========================================================================
% 3D k-Wave simulation with a focusing transducer and vinyl ring
% =========================================================================
clearvars;
close all;
DATA_CAST = 'gpuArray-single';

% -------------------------------------------------------------------------
% 1) Grids for simulation
% -------------------------------------------------------------------------
Nx = 1024;               % minimum and maximum considering the limitation of memory
Ny = 1024;               % minimum considering the size of the experimental settings
Nz = 128;                % minimum considering the size of the transducers
dx = 0.1e-3;            % maximum considering the source ultrasonic frequency
dy = 0.1e-3;            
dz = 0.1e-3;            
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

%save_path = '/mnt/sdb/matsubara/tmp'; %for dl-box
save_path = '/mnt/matsubara/rawdata' ;% for jacob

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

t_end = 40e-6;
kgrid.makeTime(medium.sound_speed, [], t_end);

% -------------------------------------------------------------------------
% 3) トランスデューサーとビニール円環の設定
% -------------------------------------------------------------------------
% ソースマスクを簡略化
source.p_mask = makeBall(Nx, Ny, Nz, round(12e-3/dx), Ny/2, Nz/2, 3);

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 4e6;        % [Hz]
tone_burst_cycles = 4;

% create the input signal using toneBurst 
source.p = source_strength .* toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% トランスデューサーの物理的特性
transducer.number_elements = 36;    % 要素数を減らす (72 -> 36)
transducer.element_width = 1;       
transducer.element_length = 6;      % 長さを短く (12 -> 6)
transducer.element_spacing = 0;     
transducer.radius = inf;            

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = 1540;                  % sound speed [m/s]
transducer.focus_distance = 25e-3;              % focus distance [m]
transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = zeros(transducer.number_elements, 1);
transducer.active_elements(11:26) = 1;          % アクティブ要素数を調整

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

% ビニール円環のマスクを作成 - サイズを調整
center = [Nx/2, Ny/2+40];
outer_radius = 40;  % 外側の半径を小さく (160 -> 40)
inner_radius = 32;  % 内側の半径を小さく (130 -> 32)

[X, Y] = meshgrid(1:Nx, 1:Ny);
X = X - center(1);
Y = Y - center(2);
R = sqrt(X.^2 + Y.^2);

% Create mask for vinyl pipe aligned with Z-axis
pipe_mask_2d = (R <= outer_radius) & (R >= inner_radius);
pipe_mask = repmat(permute(pipe_mask_2d,[2 1]), [1 1 Nz]);      % Nx×Ny×Nz

% Initialize medium parameters
medium.sound_speed = medium.sound_speed .* ones(Nx, Ny, Nz);
medium.density = medium.density .* ones(Nx, Ny, Nz);

medium.sound_speed(pipe_mask) = vinyl.sound_speed;
medium.density(pipe_mask) = vinyl.density;


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
    'DisplayMask', transducer.active_elements_mask, ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'RecordMovie', true, ...
    'MovieName', fullfile(save_path, 'vinyl_pipe_3d.avi'), ...
    'DataCast', DATA_CAST
    };

% -------------------------------------------------------------------------
% 7) トランスデューサーの可視化
% -------------------------------------------------------------------------
figure;
voxelPlot(single(pipe_mask));
hold on;
title('Transducer and Vinyl Pipe Configuration');
xlabel('X [grid points]');
ylabel('Y [grid points]');
zlabel('Z [grid points]');
view(45, 30);
saveas(gcf, fullfile(save_path, 'transducer_config_3d.png'));

% -------------------------------------------------------------------------
% 8) シミュレーション実行
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, transducer, input_args{:});

% -------------------------------------------------------------------------
% 9) 結果の可視化
% -------------------------------------------------------------------------
field=gather(sensor_data);
figure;
imagesc(field(:,:,end/2));
colorbar;
title('Pressure Field at Central Plane');
xlabel('X [grid points]');
ylabel('Y [grid points]');
saveas(gcf, fullfile(save_path, 'pressure_field_3d.png'));

% データの保存
save(fullfile(save_path, 'sensor_data_3d.mat'), 'sensor_data', '-v7.3');
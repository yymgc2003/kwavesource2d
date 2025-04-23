% =========================================================================
% 3D k-Wave simulation with a focusing transducer1 and vinyl ring
% =========================================================================
clearvars;
close all;

% Load configuration from JSON file
try
    config_file = fullfile(fileparts(mfilename('fullpath')), 'config.json');
    if ~exist(config_file, 'file')
        error('設定ファイルが見つかりません: %s', config_file);
    end
    config = jsondecode(fileread(config_file));
catch ME
    error('設定ファイルの読み込みに失敗しました: %s', ME.message);
end

DATA_CAST = config.simulation.data_cast;

% -------------------------------------------------------------------------
% 1) Grids for simulation
% -------------------------------------------------------------------------
Nx = config.grid.Nx;
Ny = config.grid.Ny;
Nz = config.grid.Nz;
dx = config.grid.dx;
dy = config.grid.dy;
dz = config.grid.dz;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

save_path = '/mnt/sdb/matsubara/tmp'; %for dl-box
%save_path = '/mnt/matsubara/rawdata' ;% for jacob
number_scan_lines = 4;

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

t_end = 100e-6;
kgrid.makeTime(medium.sound_speed, [], t_end);

% -------------------------------------------------------------------------
% 3) トランスデューサーとビニール円環の設定
% -------------------------------------------------------------------------

% define properties of the input signal
source_strength = 1e5;          % [Pa]
tone_burst_freq = 4e6;        % [Hz]
tone_burst_cycles = 4;
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
input_signal = (source_strength ./ (medium.sound_speed * medium.density)) .* input_signal;


% トランスデューサーのハードウェア情報
transducer1.number_elements = 90;    
transducer1.element_width = 1;       
transducer1.element_length = 6;      
transducer1.element_spacing = 0;     
transducer1.radius = inf;            
transducer2.number_elements = 90;    
transducer2.element_width = 1;       
transducer2.element_length = 6;      
transducer2.element_spacing = 0;     
transducer2.radius = inf;        

% calculate the width of the transducer1 in grid points
transducer_width = transducer1.number_elements * transducer1.element_width ...
    + (transducer1.number_elements - 1) * transducer1.element_spacing;

% properties used to derive the beamforming delays
transducer1.sound_speed = 1540;                  % sound speed [m/s]
transducer1.focus_distance = 25e-3;              % focus distance [m]
transducer1.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
transducer1.steering_angle = 0;                  % steering angle [degrees]
transducer1.transmit_apodization = 'Rectangular';    
transducer1.receive_apodization = 'Rectangular';
transducer2.sound_speed = 1540;                  % sound speed [m/s]
transducer2.focus_distance = 25e-3;              % focus distance [m]
transducer2.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
transducer2.steering_angle = 0;                  % steering angle [degrees]
transducer2.transmit_apodization = 'Rectangular';    
transducer2.receive_apodization = 'Rectangular';



    
% using transducer1 as source 
transducer1.input_signal = input_signal;
% calculate the width of the transducer1 in grid points
transducer_width = transducer2.number_elements * transducer2.element_width ...
    + (transducer2.number_elements - 1) * transducer2.element_spacing;
transducer1.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer1.element_length/2]);
transducer2.position = round([Nx-10, Ny/2 - transducer_width/2, Nz/2 - transducer1.element_length/2]);

transducer1.active_elements = zeros(transducer1.number_elements, 1);
transducer1.active_elements(config.transducer.config.active_elements_range(1):config.transducer.config.active_elements_range(2)) = 1;
transducer2.active_elements = zeros(transducer2.number_elements, 1);
transducer2.active_elements(config.transducer.config.active_elements_range(1):config.transducer.config.active_elements_range(2)) = 1;

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
% 8) シミュレーション実行
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer1, transducer2, input_args{:});

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
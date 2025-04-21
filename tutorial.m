% =========================================================================
% 3D k‑Wave simulation with a focusing transducer and vinyl ring (GPU)
% =========================================================================
clearvars; close all;
DATA_CAST = 'gpuArray-single';          % GPU を使う

% -------------------------------------------------------------------------
% 1) 計算グリッド
% -------------------------------------------------------------------------
Nx = 256; Ny = 256; Nz = 64;            % グリッド数
dx = 0.2e-3; dy = dx; dz = dx;          % 格子間隔 0.2 mm
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

save_path = '/mnt/sdb/matsubara/tmp';   % 出力先

% -------------------------------------------------------------------------
% 2) 媒質パラメータ（背景は水）
% -------------------------------------------------------------------------
medium.sound_speed = 1500 * ones(Nx,Ny,Nz,'single');
medium.density     = 1000 * ones(Nx,Ny,Nz,'single');

vinyl.sound_speed  = 2390;
vinyl.density      = 1400;

% -------------------------------------------------------------------------
% 3) トランスデューサ（数を 36 要素に削減）
% -------------------------------------------------------------------------
transducer.number_elements           = 36;
transducer.element_width             = 1;     % [grid pts]
transducer.element_length            = 6;     % [grid pts]
transducer.element_spacing           = 0;
transducer.position                  = ...
    round([1, Ny/2 - (transducer.number_elements*transducer.element_width)/2, Nz/2 - transducer.element_length/2]);

transducer.sound_speed               = 1540;
transducer.focus_distance            = 25e-3;
transducer.elevation_focus_distance  = 19e-3;
transducer.steering_angle            = 0;
transducer.transmit_apodization      = 'Rectangular';
transducer.receive_apodization       = 'Rectangular';
transducer.active_elements           = zeros(transducer.number_elements,1);
transducer.active_elements(11:26)    = 1;   % 16 個を発振

% k‑Wave オブジェクト化
transducer = kWaveTransducer(kgrid, transducer);
transducer.properties;                % コンソールに表示

% -------------------------------------------------------------------------
% 4) 送波信号
% -------------------------------------------------------------------------
source_strength   = 1e6;              % [Pa]
tone_burst_freq   = 4e6;              % [Hz]
tone_burst_cycles = 4;
source.p          = source_strength .* ...
                    toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% 送信面はトランスデューサ要素（kWaveTransducer が内部で設定）
source.p_mask     = transducer.source_mask;

% -------------------------------------------------------------------------
% 5) ビニール円環のマスク
% -------------------------------------------------------------------------
center = [Nx/2, Ny/2];                % x‑y 中心（z は全層同じ）
outer_radius = 40; inner_radius = 32; % [grid pts]

% x‑y 平面の円環
[x2d, y2d] = meshgrid(1:Nx, 1:Ny);    % Ny×Nx
x2d = x2d - center(1);
y2d = y2d - center(2);
R2d = sqrt(x2d.^2 + y2d.^2);
ring2d = (R2d <= outer_radius) & (R2d >= inner_radius);   % Ny×Nx

% Nx×Ny×Nz へ展開（Z 方向に複製）
pipe_mask = repmat(permute(ring2d,[2 1]), [1 1 Nz]);      % Nx×Ny×Nz

% 媒質を書き換え
medium.sound_speed(pipe_mask) = vinyl.sound_speed;
medium.density(pipe_mask)     = vinyl.density;

% -------------------------------------------------------------------------
% 6) センサー設定（中央 z スライスに 100×100 の平面）
% -------------------------------------------------------------------------
sensor.mask = zeros(Nx,Ny,Nz);
sensor_plane = zeros(Nx,Ny);
sensor_plane(Nx/2-50:Nx/2+50, Ny/2-50:Ny/2+50) = 1;
sensor.mask(:,:,Nz/2) = sensor_plane;
sensor.record = {'p'};

% -------------------------------------------------------------------------
% 7) シミュレーション設定
% -------------------------------------------------------------------------
t_end = 40e-6; kgrid.makeTime(1500, [], t_end);   % 40 µs まで
input_args = { ...
    'DataCast',            DATA_CAST, ...
    'PMLInside',           false, ...
    'PlotPML',             false, ...
    'PlotSim',             false, ...             % 描画オフで速度アップ
    'RecordMovie',         true, ...
    'MovieName',           fullfile(save_path,'vinyl_pipe_3d.mp4'), ...
    'MovieProfile',        'MPEG-4', ...
    'MovieArgs',           {'FrameRate',30} ...
    };

% -------------------------------------------------------------------------
% 8) トランスデューサとパイプの可視化
% -------------------------------------------------------------------------
figure;
voxelPlot(pipe_mask);
hold on; transducer.plotTransducer();
title('Transducer and Vinyl Pipe');
view(45,30); axis tight; xlabel('x'); ylabel('y'); zlabel('z');
saveas(gcf, fullfile(save_path,'transducer_config_3d.png'));

% -------------------------------------------------------------------------
% 9) シミュレーション実行（GPU）
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% -------------------------------------------------------------------------
% 10) 結果の可視化
% -------------------------------------------------------------------------
figure;
imagesc(squeeze(sensor_data.p(:, :, end/2)));
axis image; colorbar;
title('Pressure Field (Central z Slice)');
xlabel('x grid'); ylabel('y grid');
saveas(gcf, fullfile(save_path,'pressure_field_3d.png'));

% データ保存
save(fullfile(save_path,'sensor_data_3d.mat'), 'sensor_data','-v7.3');

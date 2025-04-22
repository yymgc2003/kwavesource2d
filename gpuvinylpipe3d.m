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
Nz = 64;                % z方向グリッド数
dx = 0.2e-3;            % グリッド間隔 [m]
dy = 0.2e-3;            
dz = 0.2e-3;            
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

save_path = '/mnt/matsubara/rawdata' ;% for jacob
if ~exist(save_path,'dir'); mkdir(save_path); end

% -------------------------------------------------------------------------
% 2) 媒質パラメータ
% -------------------------------------------------------------------------
% 水のパラメータ
medium.sound_speed = 1500 * ones(Nx,Ny,Nz,'single');   % [m/s]
medium.density     = 1000 * ones(Nx,Ny,Nz,'single');   % [kg/m^3]

% ビニールのパラメータ
vinyl.sound_speed = 2390;      % [m/s]
vinyl.density     = 1400;      % [kg/m^3]

% -------------------------------------------------------------------------
% 3) ビニール円環のマスク作成
% -------------------------------------------------------------------------
center = [Nx/2, Ny/2+40];         % グリッド上の中心 (x,y)
outer_radius = 40; 
inner_radius = 32;
[X2d,Y2d] = meshgrid(1:Nx,1:Ny);                 % Ny×Nx
ring2d = sqrt((X2d-center(1)).^2 + (Y2d-center(2)).^2);
ringMask = (ring2d<=outer_radius) & (ring2d>=inner_radius);
pipe_mask = repmat(permute(ringMask,[2 1]),[1 1 Nz]);  % Nx×Ny×Nz

% ビニール領域の物性値を設定
medium.sound_speed = gpuArray(medium.sound_speed);
medium.density = gpuArray(medium.density);
medium.sound_speed(pipe_mask) = vinyl.sound_speed;
medium.density(pipe_mask) = vinyl.density;

% -------------------------------------------------------------------------
% 4) トランスデューサ定義
% -------------------------------------------------------------------------
% 幾何情報
tr.number_elements = 36;
tr.element_width   = 1;    % [grid pt]
tr.element_length  = 6;
tr.element_spacing = 0;
tr.radius          = inf;

% 位置：計算領域 x=1 面に貼り付け，y‑z 中央付近
tr_width = tr.number_elements*tr.element_width + ...
           (tr.number_elements-1)*tr.element_spacing;
tr.position = round([1, Ny/2 - tr_width/2, Nz/2 - tr.element_length/2]);

% フォーカス・アポダイゼーション
tr.sound_speed               = 1540;
tr.focus_distance           = 25e-3;
tr.elevation_focus_distance = 19e-3;
tr.steering_angle          = 0;
tr.transmit_apodization    = 'Rectangular';
tr.receive_apodization     = 'Rectangular';
tr.active_elements         = zeros(tr.number_elements,1);
tr.active_elements(11:26)  = 1;     % 16エレメント励振

% k‑Wave オブジェクト化
transducer = kWaveTransducer(kgrid, tr);

% -------------------------------------------------------------------------
% 5) 送波信号
% -------------------------------------------------------------------------
tone_burst_freq   = 4e6;      % [Hz]
tone_burst_cycles = 4;
transducer.input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% -------------------------------------------------------------------------
% 6) センサー設定
% -------------------------------------------------------------------------
sensor.mask = zeros(Nx,Ny,Nz);
sensor.mask(round(Nx/2), :, round(Nz/2)) = 1;    % x,z 中央の y ライン
sensor.record = {'p'};

% -------------------------------------------------------------------------
% 7) シミュレーションオプション
% -------------------------------------------------------------------------
t_end = 40e-6; 
kgrid.makeTime(medium.sound_speed,[],t_end);

input_args = {...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DataCast', DATA_CAST ...
    };

% -------------------------------------------------------------------------
% 8) トランスデューサーとビニール管の可視化
% -------------------------------------------------------------------------
% トランスデューサー要素の中心座標を取得
elem_centres = transducer.element_pos;           % 3×Ne
elem_centres = elem_centres.' * 1e-3;            % [m] に換算

% cart2grid でグリッドマスク化
circle3D = cart2grid(kgrid, elem_centres);

% 3Dプロット
figure('Name','Masks','Position',[100 100 800 600]);
hp1 = voxelPlot(gather(single(pipe_mask)), kgrid);        % ビニール円環
set(hp1,'FaceColor',[1 0 0],'EdgeColor','none','FaceAlpha',0.3); 
hold on;
hp2 = voxelPlot(single(circle3D), kgrid);         % エレメント中心 (黄色)
set(hp2,'FaceColor',[1 1 0],'EdgeColor','none','FaceAlpha',0.8);

% プロット設定
axis tight; 
daspect([1 1 1]);
camlight; 
lighting gouraud;
xlabel('x [m]'); 
ylabel('y [m]'); 
zlabel('z [m]');
title('Transducer Elements and Vinyl Pipe Configuration (GPU)');
view(45,25);

% アクティブ要素のハイライト表示
active_mask = zeros(size(circle3D));
active_idx = find(tr.active_elements);
for idx = active_idx'
    elem_pos = transducer.element_pos(:,idx);
    active_mask = active_mask | (cart2grid(kgrid, elem_pos'*1e-3) > 0);
end
hp3 = voxelPlot(single(active_mask), kgrid);
set(hp3,'FaceColor',[0 1 0],'EdgeColor','none','FaceAlpha',0.8);

% 凡例追加
legend([hp1,hp2,hp3], 'Vinyl Pipe', 'All Elements', 'Active Elements', 'Location', 'northeastoutside');

% 画像保存
saveas(gcf, fullfile(save_path,'transducer_config_3d_gpu.png'));
saveas(gcf, fullfile(save_path,'transducer_config_3d_gpu.fig'));

% -------------------------------------------------------------------------
% 9) シミュレーション実行
% -------------------------------------------------------------------------
disp('=== k-Wave simulation start (GPU) ===');
sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, sensor, input_args{:});
disp('=== simulation finished ===');

% -------------------------------------------------------------------------
% 10) 結果の可視化
% -------------------------------------------------------------------------
% GPU配列をCPUに移動
pressure = gather(sensor_data);              % GPU→CPU
t = kgrid.t_array * 1e6;                    % µs

% y–time ヒートマップ
figure('Position',[100 100 600 400]);
imagesc(t, (0:Ny-1)*dy*1e3, squeeze(pressure(:)));  % Ny×Nt
axis xy;
xlabel('Time [µs]'); 
ylabel('Y [mm]');
title('Pressure along Y-axis at x,z center (GPU)');
colorbar;
saveas(gcf, fullfile(save_path,'pressure_vs_y_time_gpu.png'));

% 時刻 t=中央 の Y 断面
tIdx = round(numel(t)/2);
figure('Position',[100 100 600 400]);
plot((0:Ny-1)*dy*1e3, pressure(:,tIdx),'-o');
xlabel('Y [mm]'); 
ylabel('Pressure [Pa]');
title(sprintf('Pressure along Y (t = %.2f µs) - GPU', t(tIdx)));
grid on;
saveas(gcf, fullfile(save_path,'pressure_profile_midT_gpu.png'));

% データの保存
save(fullfile(save_path,'sensor_data_3d_gpu.mat'), 'pressure','t','-v7.3');
disp(['All results saved in ', save_path]); 
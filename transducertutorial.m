% =========================================================================
% k-Wave transducer tutorial
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
kgrid = kWaveGrid(config.grid.Nx, config.grid.dx, config.grid.Ny, config.grid.dy);
save_path = config.save_path;

% -------------------------------------------------------------------------
% 3) 媒質パラメータ
% -------------------------------------------------------------------------
medium.sound_speed = config.medium.water.sound_speed;
medium.density = config.medium.water.density;

% -------------------------------------------------------------------------
% 4) トランスデューサーの設定
% -------------------------------------------------------------------------
transducer.number_elements = config.transducer.elements;
transducer.element_width = config.transducer.element_width;
transducer.element_length = config.transducer.element_height;
transducer.element_spacing = config.transducer.element_spacing;
transducer.radius = config.transducer.radius;

% トランスデューサーの位置と向きを設定
transducer.position = config.transducer.position;
transducer.rotation = config.transducer.rotation;
transducer.focus = config.transducer.focus;

% トランスデューサーのプロパティを設定
transducer.sound_speed = config.medium.water.sound_speed;
transducer.focus_distance = norm(config.transducer.focus - config.transducer.position);
transducer.elevation_focus_distance = transducer.focus_distance;
transducer.steering_angle = 0;
transducer.transmit_apodization = 'Rectangular';
transducer.receive_apodization = 'Rectangular';

% -------------------------------------------------------------------------
% 5) ソース波形の設定
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
transducer.input_signal = source_signal;

% -------------------------------------------------------------------------
% 6) トランスデューサーの初期化
% -------------------------------------------------------------------------
transducer = kWaveTransducer(kgrid, transducer);

% -------------------------------------------------------------------------
% 7) センサーの設定
% -------------------------------------------------------------------------
sensor.mask = zeros(config.grid.Nx, config.grid.Ny);
sensor_x = config.grid.Nx/2 + config.sensor.x_offset;
sensor_y = config.grid.Ny/2 + config.sensor.y_offset;

% stream the data to disk in blocks of 100 if storing the complete time
% history 
if ~USE_STATISTICS
    input_args = [input_args {'StreamToDisk', 100}];
end
sensor.record = {'p'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, sensor, input_args{:});

% =========================================================================
% COMPUTE THE BEAM PATTERN USING SIMULATION STATISTICS
% =========================================================================
voxelPlot(double(transducer | cart2grid(kgrid, sensor.mask)));
view(127, 18);
saveas(gcf, fullfile(save_path, 'trans_config_3d_tut.png'));
if USE_STATISTICS
    
    % reshape the returned rms and max fields to their original position
    sensor_data.p_rms = reshape(sensor_data.p_rms, [Nx, Nj]);
    sensor_data.p_max = reshape(sensor_data.p_max, [Nx, Nj]);
    
    % plot the beam pattern using the pressure maximum
    figure;
    imagesc(j_vec * 1e3, (kgrid.x_vec - min(kgrid.x_vec(:))) * 1e3, sensor_data.p_max * 1e-6);
    xlabel([j_label '-position [mm]']);
    ylabel('x-position [mm]');
    title('Total Beam Pattern Using Maximum Of Recorded Pressure');
    colormap(jet(256));
    c = colorbar;
    ylabel(c, 'Pressure [MPa]');
    axis image;
    
    % plot the beam pattern using the pressure rms
    figure;
    plot(kgrid.t_array*1e3, sensor_data.p(1, :));
    xlabel('Time [ms]');
    ylabel('Pressure [Pa]');
    title('Pressure at the sensor');
    saveas(gcf, fullfile(save_path, 'sensor_transducer_tutorial.png')); 
end

% kspaceFirstOrder2D実行後にGPU配列をCPUに集約
sensor_data_cpu = structfun(@gather, sensor_data, 'UniformOutput', false);

% v7.3 形式で.matファイル保存
save(fullfile(save_path, 'sensor_data_transducer_tutorial.mat'), 'sensor_data_cpu', '-v7.3'); 
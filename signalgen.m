function signalgen_module_all()
% この関数はconfig.jsonを読み込み、指定ディレクトリ内のlocation1.csv, location2.csv, ... すべてに対して
% signalgen_moduleを順次実行します。

    % configファイルのパス
    config_file = 'config.json';

    % config.jsonからlocation_seedfiles_pathを取得
    config = jsondecode(fileread(config_file));
    location_dir = config.location_seedfiles_path;

    % location*.csvファイルをすべて取得
    files = dir(fullfile(location_dir, 'location*.csv'));
    fprintf('--- %d location*.csv files found. ---\n', length(files));
    if isempty(files)
        error('No location*.csv files found in %s.', location_dir);
    end

    % 各CSVファイルに対してsignalgen_moduleを実行
    for i = 1:length(files)
        location_csv = fullfile(location_dir, files(i).name);
        fprintf('--- Processing %s ---\n', location_csv);
        signalgen_module(config_file, location_csv);
    end
end

function signalgen_module(config_file, location_csv)
% この関数は指定されたconfig.jsonとlocationX.csvを使ってk-Waveシミュレーションを実行し、
% location番号付きの信号データファイルと画像を出力します。
% 例: signalgen_module('config.json', 'location3.csv')

    % location_csvファイル名から番号を抽出
    [~, locname, ~] = fileparts(location_csv);
    locnum = regexp(locname, '\d+', 'match');
    if isempty(locnum)
        locnum_str = '';
    else
        locnum_str = locnum{1};
    end

    % 設定ファイルの読み込み
    config = jsondecode(fileread(config_file));
    save_path = config.save_path;

    % シミュレーション設定
    DATA_CAST = 'gpuArray-single';

    % PMLサイズ設定
    PML_X_SIZE = 20; % [grid points]
    PML_Y_SIZE = 10; % [grid points]
    PML_Z_SIZE = 10; % [grid points]

    % PMLを除いたグリッドポイント数
    Nx = config.grid.Nx - 2*PML_X_SIZE;
    Ny = config.grid.Ny - 2*PML_Y_SIZE;
    Nz = config.grid.Nz - 2*PML_Z_SIZE;

    % グリッド間隔
    dx = config.grid.dx; dy = config.grid.dy; dz = config.grid.dz;

    % k-spaceグリッド作成
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

    % 伝播媒体のプロパティ
    medium.sound_speed = config.medium.water.sound_speed * ones(Nx, Ny, Nz);
    medium.density = config.medium.water.density * ones(Nx, Ny, Nz);
    medium.alpha_coeff = config.medium.water.alpha_coeff * ones(Nx, Ny, Nz);
    medium.alpha_power = config.medium.water.alpha_power;
    medium.BonA = 6;

    % 時間配列作成
    t_end = config.simulation.t_end;
    cfl = config.simulation.CFL;
    kgrid.makeTime(medium.sound_speed, cfl, t_end);

    % 入力信号のプロパティ
    source_strength = config.source.source_strength;
    tone_burst_freq = config.source.tone_burst_freq;
    tone_burst_cycles = config.source.tone_burst_cycles;

    % 入力信号生成
    input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

    % インピーダンスでスケーリング
    input_signal = (source_strength / (config.medium.water.sound_speed * config.medium.water.density)) * input_signal;

    % 送信トランスデューサ定義
    transducer_transmit.number_elements = 180;
    transducer_transmit.element_width = 1;
    transducer_transmit.element_length = 12;
    transducer_transmit.element_spacing = 0;
    transducer_transmit.radius = inf;

    % トランスデューサ幅
    transducer_transmit_width = transducer_transmit.number_elements * transducer_transmit.element_width ...
        + (transducer_transmit.number_elements - 1) * transducer_transmit.element_spacing;

    % 配置
    transducer_transmit.position = round([10, Ny/2 - transducer_transmit_width/2, Nz/2 - transducer_transmit.element_length/2]);

    % ビームフォーミング用プロパティ
    transducer_transmit.sound_speed = 1540;
    transducer_transmit.focus_distance = 20e-3;
    transducer_transmit.elevation_focus_distance = 19e-3;
    transducer_transmit.steering_angle = 0;

    % アポダイゼーション
    transducer_transmit.transmit_apodization = 'Rectangular';
    transducer_transmit.receive_apodization = 'Rectangular';

    % アクティブエレメント
    transducer_transmit.active_elements = zeros(transducer_transmit.number_elements, 1);
    transducer_transmit.active_elements(41:141) = 1;

    % 入力信号
    transducer_transmit.input_signal = input_signal;

    % トランスデューサ作成
    transducer_transmit = kWaveTransducer(kgrid, transducer_transmit);

    % プロパティ表示
    transducer_transmit.properties;

    % 受信トランスデューサ定義
    transducer_receive.number_elements = 90;
    transducer_receive.element_width = 1;
    transducer_receive.element_length = 12;
    transducer_receive.element_spacing = 0;
    transducer_receive.radius = inf;
    transducer_width = transducer_receive.number_elements * transducer_receive.element_width ...
        + (transducer_receive.number_elements - 1) * transducer_receive.element_spacing;
    transducer_receive.position = round([Nx-10, Ny/2 - transducer_width/2, Nz/2 - transducer_transmit.element_length/2]);
    transducer_receive.sound_speed = config.medium.water.sound_speed;
    transducer_receive.focus_distance = 25e-3;
    transducer_receive.elevation_focus_distance = 19e-3;
    transducer_receive.steering_angle = 0;
    transducer_receive.transmit_apodization = 'Rectangular';
    transducer_receive.receive_apodization = 'Rectangular';
    transducer_receive.active_elements = zeros(transducer_receive.number_elements, 1);
    transducer_receive.active_elements(21:52) = 1;
    transducer_receive = kWaveTransducer(kgrid, transducer_receive);

    % センサーマスク
    sensor.mask = zeros(Nx, Ny, Nz);
    sensor.mask([Nx/4, Nx/2, 3*Nx/4], Ny/2, Nz/2) = 1;
    cx = config.pipe.center_x;
    cy = config.pipe.center_y;
    cz = config.pipe.center_z;
    outer_r_mm = 16;
    inner_r_mm = 13;
    outer_r = round((outer_r_mm * 1e-3) / dx);
    inner_r = round((inner_r_mm * 1e-3) / dx);
    [Xg, Yg] = ndgrid(1:Nx, 1:Ny);
    ring2d = sqrt((Xg-cx).^2 + (Yg-cy).^2);
    ringMask = (ring2d <= outer_r) & (ring2d >= inner_r);

    pipe_mask = repmat(ringMask, [1 1 Nz]);
    medium.sound_speed(pipe_mask == 1) = config.medium.vinyl.sound_speed;
    medium.density(pipe_mask == 1) = config.medium.vinyl.density;

    % ガラスマスク
    % locationX.csvから座標を読み込む
    location = csvread(location_csv);
    radius_pts = round(1.25e-3 / dx);

    % glass_mask初期化
    glass_mask = zeros(Nx, Ny, Nz);

    % 各座標にボールを配置
    for i = 1:size(location,1)
        loc_seed = location(i,:);
        bx = round(inner_r * loc_seed(1)) + cx;
        by = round(inner_r * loc_seed(2)) + cy;
        bz = round(Nz * loc_seed(3));
        fprintf('bx: %d, by: %d, bz: %d\n', bx, by, bz);
        glass_mask = glass_mask | makeBall(Nx, Ny, Nz, bx, by, bz, radius_pts);
    end
    medium.sound_speed(glass_mask == 1) = config.medium.glass.sound_speed;
    medium.density(glass_mask == 1) = config.medium.glass.density;
    medium.alpha_coeff(glass_mask == 1) = config.medium.glass.alpha_coeff;

    % シミュレーション実行
    display_mask = transducer_transmit.all_elements_mask | pipe_mask | glass_mask;
    input_args = {'DisplayMask', display_mask, ...
        'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
        'RecordMovie', true, ...
        'MovieName', fullfile(save_path, ['solid_liquid' locnum_str '.avi']), ...
        'DataCast', DATA_CAST, 'PlotScale', [-1/2, 1/2] * source_strength};

    sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer_transmit, transducer_transmit, input_args{:});
    save(fullfile(save_path, ['solid_liquid' locnum_str '.mat']), 'sensor_data', 'kgrid', '-v7.3');

    % 信号波形のプロットと保存
    scan_line = transducer_transmit.scan_line(sensor_data);
    figure(1);
    plot(kgrid.t_array * 1e6, scan_line * 1e-6, 'b-');
    xlabel('Time [\mus]');
    ylabel('Pressure [MPa]');
    ylim([-2 2]);
    title('Signal from Transducer transmit');
    grid on;
    saveas(gcf, fullfile(save_path, ['signal_solid_liquid_reflector' locnum_str '.png']));

    % 実験セットアップの可視化と保存
    figure(21); clf;
    hold on;

    % マスクの次元を可視化用に並べ替え
    transmit_mask = permute(transducer_transmit.all_elements_mask, [2,1,3]);
    sensor_mask   = permute(sensor.mask, [2,1,3]);
    pipe_mask_p   = permute(pipe_mask, [2,1,3]);
    glass_mask_p  = permute(glass_mask, [2,1,3]);
    % 可視化用にdouble型へ
    transmit_mask_double = double(transmit_mask);
    sensor_mask_double   = double(sensor_mask);
    pipe_mask_double     = double(pipe_mask_p);
    glass_mask_double    = double(glass_mask_p);

    % 送信トランスデューサマスク（青）
    if any(transmit_mask_double(:))
        p1 = patch(isosurface(transmit_mask_double, 0.5));
        isonormals(transmit_mask_double, p1);
        set(p1, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    end

    % パイプマスク（緑）
    if any(pipe_mask_double(:))
        p3 = patch(isosurface(pipe_mask_double, 0.3));
        isonormals(pipe_mask_double, p3);
        set(p3, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end

    % ガラスマスク（赤）
    if any(glass_mask_double(:))
        p4 = patch(isosurface(glass_mask_double, 0.3));
        isonormals(glass_mask_double, p4);
        set(p4, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end

    % 軸とビュー設定
    axis equal;
    xlim([1, Nx]);
    ylim([1, Ny]);
    zlim([1, Nz]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Transducer, Sensor, and Pipe Mask Visualization');
    view(80, 20);
    camlight;
    lighting gouraud;
    grid on;
    hold off;

    saveas(gcf, fullfile(save_path, ['experimental_setup' locnum_str '.png']));
end
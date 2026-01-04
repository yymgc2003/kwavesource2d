function kwavesim(config_path, save_path, loc_path, loc_path_solid, locnum)
    DATA_CAST='gpuArray-single';

    config=jsondecode(fileread(config_path));
    save_loc_path=fullfile(save_path,"location_seed");
    save_data_path=fullfile(save_path,"data");
    save_pic_path=fullfile(save_path,"logs");

    num_gnu_pic = 200;

    % (48mm+50grid)×(32mm+20grid)必要
    % dx=0.012mmの場合、Nx=4096, Ny=2916
    % dx=0.02mmの場合、Nx=2592, Ny=1728
    % dx=0.04mmの場合、Nx=1296, Ny=864
    % dx=0.05mmの場合、Nx=1024, Ny=768

    Nx = config.grid.Nx;
    Ny = config.grid.Ny;
    dx = config.grid.dx;
    dy = config.grid.dy;
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    t_end = config.simulation.t_end;
    cfl = config.simulation.CFL;

    water.sound_speed = config.medium.water.sound_speed;
    water.density     = config.medium.water.density;
    water.alpha_coeff = config.medium.water.alpha_coeff;
    water.alpha_power = config.medium.water.alpha_power;
    water.BonA = config.medium.water.BonA;

    medium.sound_speed = water.sound_speed * ones(Nx, Ny);
    medium.density = water.density * ones(Nx, Ny);
    medium.alpha_coeff = water.alpha_coeff * ones(Nx, Ny);
    medium.alpha_power = water.alpha_power;
    medium.BonA = water.BonA * ones(Nx, Ny);

    air.sound_speed = config.medium.air.sound_speed;
    air.density = config.medium.air.density;
    air.alpha_coeff = config.medium.air.alpha_coeff;
    air.BonA = config.medium.air.BonA;

    glass.sound_speed = config.medium.glass.sound_speed;
    glass.density = config.medium.glass.density;
    glass.alpha_coeff = config.medium.glass.alpha_coeff;
    glass.BonA = config.medium.glass.BonA;

    % config.json内のパイプの中心はtdx1からの相対座標！！
    % tdx1の位置は(30, Ny/2)で固定！！
    distance_tdx = config.pipe.distance_tdx;
    inner_radius = config.pipe.inner_radius;
    outer_radius = config.pipe.outer_radius;
    pcx = 30 + round(distance_tdx/dx); % pipe center x in grid
    pcy = Ny/2; % pipe center y in grid
    outer_radius_sqrt2_grid = round(outer_radius/dx/sqrt(2));

    pipe_mask = zeros(Nx, Ny);
    pipe_mask = pipe_mask | makeDisc(Nx, Ny, pcx, pcy, ...
                            round(inner_radius/dx));

    bubble_diameter = config.bubble.diameter;
    glass_diameter = config.glass.diameter;

    sample_table=readtable(loc_path);
    samples=table2array(sample_table);
    samples=transpose(samples);

    bubble_mask = zeros(Nx, Ny);
    glass_mask = zeros(Nx, Ny);
    sz = size(samples);
    sz = sz(2);
    disp(sz);

    for count=1:sz
        candidate = samples(1:2, count);
        cur_diameter_bubble = samples(3, count);
        minor_axis_length = samples(4, count);
        candidate=(candidate/dx);
        euler_angles = samples(5, count);
        bc=candidate+[pcx;pcy];
        bx=bc(1);by=bc(2);
        r_maj=(cur_diameter_bubble/dx/2);
        r_min=(minor_axis_length/dx/2);
        Q = [cos(euler_angles), -sin(euler_angles);sin(euler_angles),cos(euler_angles)];
        for xx=max(1,floor(bx-r_maj)-1):min(Nx,floor(1+bx+r_maj)+1)
            for yy=max(1,floor(by-r_maj)-1):min(Ny,floor(1+by+r_maj)+1)
                bubble_mask_relative = [xx, yy]-[bx, by];
                if bubble_mask_relative*Q*diag([1/r_maj^2,1/r_min^2])*transpose(Q)*transpose(bubble_mask_relative)<=1
                    bubble_mask(xx, yy) = 1;
                end
            end
        end
    end

    sample_table=readtable(loc_path_solid);
    samples=table2array(sample_table);
    samples=transpose(samples);
    sz = size(samples);
    sz = sz(2);
    disp(sz);
    gr=round(glass_diameter/2/dx);
    for count=1:sz
        bx=samples(1,count);
        by=samples(2,count);
        bx=round(bx/dx)+pcx;
        by=round(by/dx)+pcy;
        glass_mask=glass_mask | makeDisc(Nx,Ny,bx,by,gr);
    end

    % bubble_mask = bubble_mask | makeDisc(Nx, Ny, pcx, pcy-y_offset,...
    %                     round(bubble_diameter/2/dx));
    % bubble_mask = bubble_mask | makeDisc(Nx, Ny, pcx, pcy+y_offset,...
    %                     round(bubble_diameter/2/dx));   

    medium.sound_speed(glass_mask==1) = glass.sound_speed;
    medium.density(glass_mask==1) = glass.density;
    medium.alpha_coeff(glass_mask==1) = glass.alpha_coeff;

    medium.sound_speed(bubble_mask==1)=air.sound_speed;
    medium.density(bubble_mask==1)=air.density;
    medium.alpha_coeff(bubble_mask==1)=air.alpha_coeff;

    % medium.sound_speed_ref = water.sound_speed;
    kgrid.makeTime(water.sound_speed, cfl, t_end);

    % tdx1の位置は(30, Ny/2)で固定！！
    % config.jsonのtdx2以上はパイプ中心からの相対座標！！

    tone_burst_cycles = config.transducer1.tone_burst_cycles;
    tone_burst_freq = config.transducer1.tone_burst_freq;

    tdx1_pressure = config.transducer1.pressure;

    tdx1_diameter = config.transducer1.diameter;
    tdx1_radius_grid = round(tdx1_diameter/2/dx);
    % tdx1_pos = [30, Ny/2];
    % tdx1_focus_distance = config.transducer1.focus_distance;
    % tdx1_steering_angle = config.transducer1.steering_angle;
    % tdx1_mask = makeArc([Nx, Ny], tdx1_pos, ...
    %                     round(2*tdx1_focus_distance/dx), ...
    %                     round(tdx1_diameter/dx), ...
    %                     [tdx1_pos(1)+round(tdx1_focus_distance/dx), tdx1_pos(2)]);
    tdx1_mask = zeros(Nx, Ny);
    tdx1_mask(30, Ny/2-tdx1_radius_grid:Ny/2+tdx1_radius_grid) = 1;
    source.u_mask= tdx1_mask;

    input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

    input_signal = tdx1_pressure/water.sound_speed/water.density*input_signal;
    source.ux = input_signal;

    sensor_mask = zeros(Nx, Ny);
    sensor_mask(30, pcy) = 1;
    % distance_tdx_sqrt2_grid = round(distance_tdx/sqrt(2)/dx);
    % tdx2_mask = makeArc([Nx,Ny],...
    %                     [pcx-outer_radius_sqrt2_grid, pcy-outer_radius_sqrt2_grid],...
    %                     inf,2*tdx1_radius_grid+1,[pcx, pcy]);
    % sensor_mask(tdx2_mask==1) = 1;
    % tdx3_mask = makeArc([Nx,Ny],...
    %                     [pcx+outer_radius_sqrt2_grid, pcy-outer_radius_sqrt2_grid],...
    %                     inf,2*tdx1_radius_grid+1,[pcx, pcy]);
    % sensor_mask(tdx3_mask==1) = 1;
    % tdx4_mask = zeros(Nx, Ny);
    % tdx4_mask(pcx+round(outer_radius/dx), Ny/2-tdx1_radius_grid:Ny/2+tdx1_radius_grid) = 1;
    % sensor_mask(tdx4_mask==1) = 1;
    sensor_mask(pcx-outer_radius_sqrt2_grid,pcy-outer_radius_sqrt2_grid) = 1;
    sensor_mask(pcx+outer_radius_sqrt2_grid,pcy-outer_radius_sqrt2_grid) = 1;
    sensor_mask(pcx+round((outer_radius)/dx),pcy) = 1;
    sensor_mask(pcx,pcy-round((outer_radius)/dx)) = 1;
    sensor.mask = zeros(Nx, Ny);
    sensor.mask(sensor_mask==1)=1;
    sensor.record = {'p'};

    if locnum~=1
        input_args = {
            'PMLInside', true, 'PlotPML', true, ...
            'PMLSize', config.simulation.pml_size, ...
            'RecordMovie', false, ...
            'PlotFreq', 50, ...
            'DataCast', DATA_CAST, ...
            'DeviceNum', 2, ...
            };
    else
        input_args = {
            'PMLInside', true, 'PlotPML', true, ...
            'PMLSize', config.simulation.pml_size, ...
            'RecordMovie', false, ...
            'PlotFreq', 50, ...
            'DataCast', DATA_CAST, ...
            'UseGnuplot', true, ...
            'GnuplotFPS', num_gnu_pic/t_end, ...
            'GnuplotOutpath', fullfile(save_path, "gnu_csv"), ...
            'GnuplotTimeRange', [0e-6, 40e-6], ...
            };
    end

    transmit_mask_double = double(tdx1_mask);
    pipe_mask_double   = double(pipe_mask);
    sensor_mask_double     = double(sensor_mask);
    bubble_mask_double    = double(bubble_mask);
    glass_mask_double     = double(glass_mask);

    pipe_mask_double = pipe_mask_double - glass_mask_double - bubble_mask_double;

    composite_mask = pipe_mask_double * 2 + ...
                    glass_mask_double * 1 + ...
                    bubble_mask_double * 3;

    my_colormap = [
    0.0 0.0 1.0;
    1.0 0.0 0.0;
    0.5 0.5 1.0;  % 0: 背景 (青)
    1.0 1.0 1.0;  % 3: bubble_mask (白)
    ];

    glass_fraction=sum(glass_mask,"all")/sum(pipe_mask,"all");
    bubble_fraction=sum(bubble_mask,"all")/sum(pipe_mask,"all");
    liquid_fraction=1-glass_fraction-bubble_fraction;
    frac_path=fullfile(save_loc_path, sprintf('fraction%d.csv',locnum));
    frac_id=fopen(frac_path, 'w');
    fprintf(frac_id, "%f,%f,%f",glass_fraction,bubble_fraction,liquid_fraction);
    fclose(frac_id);

    % プロット
    figure;
    composite_mask = flip(composite_mask, 2);
    composite_mask = transpose(composite_mask);
    imagesc(composite_mask); 
    colormap(my_colormap); % 定義したカラーマップを適用
    axis off;
    fontsize(24, "points");
    title(sprintf("sgl=%.2g, %.2g, %.2g", glass_fraction, bubble_fraction,liquid_fraction))
    axis equal tight;
                
    saveas(gcf, fullfile(save_pic_path, sprintf('experimental_setup%d.png',locnum)));

    if locnum~=1
        sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
    else
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    end

    p = sensor_data.p;
    p=gather(p);
    [filepath, filename, ext] = fileparts(mfilename('fullpath'));
    save(fullfile(save_data_path,sprintf('sim_data%d.mat',locnum)), ...
        'p','kgrid', '-v7.3');
    
    figure;
    plot(1e6*kgrid.t_array, p(1,:));
    xlim([0,t_end*1e6]);
    saveas(gcf, fullfile(save_pic_path,sprintf('tdx1_signal%d.png',locnum)));

    figure;
    plot(1e6*kgrid.t_array, p(2,:));
    xlim([0,t_end*1e6]);
    saveas(gcf, fullfile(save_pic_path,sprintf('tdx2_signal%d.png',locnum)));

    figure;
    plot(1e6*kgrid.t_array, p(3,:));
    xlim([0,t_end*1e6]);
    saveas(gcf, fullfile(save_pic_path,sprintf('tdx3_signal%d.png',locnum)));

    figure;
    plot(1e6*kgrid.t_array, p(4,:));
    xlim([0,t_end*1e6]);
    saveas(gcf, fullfile(save_pic_path,sprintf('tdx4_signal%d.png',locnum)));

    figure;
    plot(1e6*kgrid.t_array, p(5,:));
    xlim([0,t_end*1e6]);
    saveas(gcf, fullfile(save_pic_path,sprintf('tdx5_signal%d.png',locnum)));
end
function plot_gen3d_gl(config_file, location_csv, locnum_str, ... 
    save_full_path)

    config = jsondecode(fileread(config_file));
    save_logs_path = fullfile(save_full_path, 'logs');
    save_data_path = fullfile(save_full_path, 'data');

    glass_radius = config.simulation.glass_radius;

    Nx = config.grid.Nx;
    Ny = config.grid.Ny; 
    Nz_sim = config.grid.Nz;
    dx = config.grid.dx; dy = config.grid.dy; dz = config.grid.dz;
    cx = config.pipe.center_x;
    cy = config.pipe.center_y;
    cz = config.pipe.center_z;
    outer_r_mm = config.pipe.outer_radius;
    inner_r_mm = config.pipe.inner_radius;
    outer_r = round((outer_r_mm * 1e-3) / dx);
    inner_r = round((inner_r_mm * 1e-3) / dx);
    Nz = inner_r * 2;
    [Xg, Yg] = ndgrid(1:Nx, 1:Ny);
    ring2d = sqrt((Xg-cx).^2 + (Yg-cy).^2);
    ringMask = (ring2d <= outer_r) & (ring2d >= inner_r);

    pipe_mask = repmat(ringMask, [1 1 Nz]);

    bubble_mask = zeros(Nx, Ny, Nz);

    flow_pattern = config.simulation.flow_pattern;

    figure;
    hold on;
    if flow_pattern == "bubble"
        location_df = readtable(location_csv);
        location = table2array(location_df);
        for i = 1:size(location,1)
            loc_seed = location(i,:);
            bx = round(inner_r * loc_seed(1)) + cx;
            by = round(inner_r * loc_seed(2)) + cy;
            bz = round(inner_r * loc_seed(3)) + Nz/2;
            radius_pts = round(inner_r*loc_seed(4)/2);
            if bz-radius_pts-1<Nz & bz+radius_pts+1>1
                radius_pts_short = round(inner_r*loc_seed(5)/2);
                %fprintf('bx: %d, by: %d, bz: %d\n', bx, by, bz);
                th1 = loc_seed(6); th2 = loc_seed(7); th3 = loc_seed(8);
                roll = [1,0,0;0,cos(th1),-sin(th1);0,sin(th1),cos(th1)];
                pitch =[cos(th2),0,sin(th2);0,1,0;-sin(th2),0,cos(th2)];
                yaw =  [cos(th3),-sin(th3),0;sin(th3),cos(th3),0;0,0,1];
                Q = roll*pitch*yaw;
                for x=bx-radius_pts-1:bx+radius_pts+1
                    for y=by-radius_pts-1:by+radius_pts+1
                        for z=max(1,bz-radius_pts-1):min(Nz,bz+radius_pts+1)
                            bubble_mask_relative = [x, y, z]-[bx, by, bz];
                            if bubble_mask_relative*Q*diag([1/radius_pts^2,1/radius_pts^2,1/radius_pts_short^2])*transpose(Q)*transpose(bubble_mask_relative)<=1
                                bubble_mask(x, y, z) = 1;
                            end
                        end
                    end
                end
            end
        end
        fprintf('Mask Finished.\n')
    end
    if flow_pattern == "slug"
        location_df = readtable(location_csv);
        location = table2array(location_df);
        %location: スラグの中心、スラグ長さ、楕円の累乗の値
        bz = round(location(1)*inner_r) + Nz/2;
        major_axis_length = round(location(2)*inner_r);
        minor_axis_length = round(location(3)*inner_r);
        slug_pow_num = location(4);
        for z=1:Nz
            pow_z_relative =abs((z-bz)/major_axis_length)^slug_pow_num;
            for x = cx-inner_r-1:cx+inner_r+1
                for y = cy-inner_r-1:cy+inner_r+1
                    pow_r_relative = (((x-cx)^2+(y-cy)^2)/minor_axis_length^2)^(slug_pow_num/2);
                    if pow_z_relative + pow_r_relative <= 1 & z > bz
                        bubble_mask(x, y, z) = 1;
                    end
                end
            end
        end
    end
    if flow_pattern == "slug-bubble"
        location_csv(1)
        location_df = readtable(location_csv{1});
        location = table2array(location_df);
        for i = 1:size(location,1)
            loc_seed = location(i,:);
            bx = round(inner_r * loc_seed(1)) + cx;
            by = round(inner_r * loc_seed(2)) + cy;
            bz = round(inner_r * loc_seed(3)) + Nz/2;
            radius_pts = round(inner_r*loc_seed(4)/2);
            if bz-radius_pts-1<Nz & bz+radius_pts+1>1
                radius_pts_short = round(inner_r*loc_seed(5)/2);
                %fprintf('bx: %d, by: %d, bz: %d\n', bx, by, bz);
                th1 = loc_seed(6); th2 = loc_seed(7); th3 = loc_seed(8);
                roll = [1,0,0;0,cos(th1),-sin(th1);0,sin(th1),cos(th1)];
                pitch =[cos(th2),0,sin(th2);0,1,0;-sin(th2),0,cos(th2)];
                yaw =  [cos(th3),-sin(th3),0;sin(th3),cos(th3),0;0,0,1];
                Q = roll*pitch*yaw;
                for x=bx-radius_pts-1:bx+radius_pts+1
                    for y=by-radius_pts-1:by+radius_pts+1
                        for z=max(1,bz-radius_pts-1):min(Nz,bz+radius_pts+1)
                            bubble_mask_relative = [x, y, z]-[bx, by, bz];
                            if bubble_mask_relative*Q*diag([1/radius_pts^2,1/radius_pts^2,1/radius_pts_short^2])*transpose(Q)*transpose(bubble_mask_relative)<=1
                                bubble_mask(x, y, z) = 1;
                            end
                        end
                    end
                end
            end
        end
        location_df = readtable(location_csv(2));
        location = table2array(location_df);
        %location: スラグの中心、スラグ長さ、楕円の累乗の値
        bz = round(location(1)*inner_r);
        major_axis_length = round(location(2)*inner_r);
        minor_axis_length = round(location(3)*inner_r);
        slug_pow_num = location(4);
        for z=1:Nz
            pow_z_relative =abs((z-bz)/major_axis_length)^slug_pow_num;
            for x = cx-inner_r-1:cx+inner_r+1
                for y = cy-inner_r-1:cy+inner_r+1
                    pow_r_relative = (((x-cx)^2+(y-cy)^2)/minor_axis_length^2)^(slug_pow_num/2);
                    if pow_z_relative + pow_r_relative <= 1 & z > bz
                        bubble_mask(x, y, z) = 1;
                    end
                end
            end
        end
    end
    fprintf('Mask Finished.\n')
    s_idx = cx - Ny/2 + 1;
    e_idx = cx + Ny/2;
    rgb_img = ones(Ny, Ny, 3);
    rgb_img(:,:,:) = rgb_img(:,:,:) - 0.5 * pipe_mask(s_idx:e_idx, :, Nz/2);
    for i=1:2
        rgb_img(:,:,i) = rgb_img(:,:,i) - permute(bubble_mask(s_idx:e_idx, :, Nz/2), [2 1 3]);
    end
    figure;
    imshow(rgb_img);
    axis tight;
    axis equal;
    axis xy;
    axis image;
    set(gca,'Units','normalized');
    set(gca,'Position',[0 0 1 1]);
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, fullfile(save_logs_path, ['experimental_setup' locnum_str '.png']));


    s_idx = cx - Ny/2 + 1;
    e_idx = cx + Ny/2;
    rgb_img = ones(Ny, Ny, 3);
    s_idx_z = (Ny - Nz)/2 + 1;
    e_idx_z = (Ny + Nz)/2;
    rgb_img(s_idx_z:e_idx_z,:,:) = rgb_img(s_idx_z:e_idx_z,:,:) - 0.5 * permute(pipe_mask(s_idx:e_idx, cy, :),[3 1 2]);
    for i=1:2
        rgb_img(s_idx_z:e_idx_z,:,i) = rgb_img(s_idx_z:e_idx_z,:,i) - permute(bubble_mask(s_idx:e_idx, cy, :), [3 1 2]);
    rgb_img((Ny-Nz_sim)/2+1:(Ny+Nz_sim)/2,1:6,:) = 0;
    rgb_img((Ny-Nz_sim)/2+1:(Ny+Nz_sim)/2,Ny-5:Ny,:) = 0;
    rgb_img((Ny-Nz_sim)/2-4:(Ny-Nz_sim)/2+1,:,:) = 0;
    rgb_img((Ny+Nz_sim)/2:(Ny+Nz_sim)/2+5,:,:) = 0;

    end
    figure;
    imshow(rgb_img);
    axis tight;
    axis equal;
    axis xy;
    axis image;
    set(gca,'Units','normalized');
    set(gca,'Position',[0 0 1 1]);
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, fullfile(save_logs_path, ['experimental_setup_vertical' locnum_str '.png']));
end
function C=pipe_location_gen(m, bubble_num_max, config_file)
    config = jsondecode(fileread(config_file));
    mu = [0, 0];
    sigma = eye(2);
    k_gamma=config.bubble.k_gamma;
    theta_gamma=config.bubble.theta_gamma;

    inner_radius=config.pipe.inner_radius;
    min_diameter_bubble=config.bubble.min_diameter/inner_radius;
    max_diameter_bubble=config.bubble.max_diameter/inner_radius;
    min_dist=config.bubble.distance/inner_radius;
    min_dist_wall=config.bubble.distance_wall/inner_radius;
    slug_min_d = config.bubble.slug_min_diameter/inner_radius;
    slug_max_d = config.bubble.slug_max_diameter/inner_radius;
    annular_max_d = config.bubble.annular_max_diameter/inner_radius;
    annular_min_d = config.bubble.annular_min_diameter/inner_radius;

    glass_radius = config.glass.diameter/2/inner_radius;

    flow_pattern=config.bubble.flow_pattern;

    oblate_border=1.2e-3/inner_radius;
    prolate_border=6e-3/inner_radius;
    Nx=300;
    num_bubble_max = 100;
    num_bubble_max_annular = 5;
    pd = makedist('Gamma','a',k_gamma,'b',theta_gamma);
    max_attempts_radius=100;
    max_attempts = 4000;
    max_attempts_solid=4000;

    attempts_radius=0;
    while attempts_radius < max_attempts_radius
        num_bubble = 0;
        cur_gas_fraction = 0;
        diameter_bubble = zeros(num_bubble_max,2);

        if flow_pattern=="slug"
            num_bubble = num_bubble+1;
            cur_diameter_bubble = (slug_max_d-slug_min_d)*rand+slug_min_d;
            minor_axis_length=cur_diameter_bubble;
            diameter_bubble(num_bubble,:)=[cur_diameter_bubble,minor_axis_length];
            cur_gas_fraction = cur_gas_fraction + minor_axis_length*cur_diameter_bubble*0.25;
        elseif flow_pattern=="annular"
            num_bubble = num_bubble+1;
            cur_diameter_bubble = (annular_max_d-annular_min_d)*rand+annular_min_d;
            minor_axis_length=cur_diameter_bubble;
            diameter_bubble(num_bubble,:)=[cur_diameter_bubble,minor_axis_length];
            cur_gas_fraction = cur_gas_fraction + minor_axis_length*cur_diameter_bubble*0.25;
        end

        while num_bubble < bubble_num_max
            % cur_diameter_bubble = random(pd)*1e-3/inner_radius;
            cur_diameter_bubble = rand*(max_diameter_bubble-min_diameter_bubble)+min_diameter_bubble;
            if cur_diameter_bubble < max_diameter_bubble && cur_diameter_bubble> min_diameter_bubble
                num_bubble = num_bubble+1;
                if cur_diameter_bubble>oblate_border && cur_diameter_bubble<prolate_border
                    minor_axis_length = cur_diameter_bubble/(0.3*rand^2+1);
                else
                    minor_axis_length=cur_diameter_bubble;
                end
                diameter_bubble(num_bubble,:)=[cur_diameter_bubble,minor_axis_length];
                cur_gas_fraction = cur_gas_fraction + minor_axis_length*cur_diameter_bubble*0.25;
            end
        end
        diameter_bubble = sortrows(diameter_bubble, 1, 'descend');
        samples = zeros(num_bubble,5);
        fprintf("Number of bubbles: %d\n", num_bubble);
        count = 0;
        attempts = 0;

        c_mask=zeros(2*Nx,2*Nx);

        if flow_pattern=="annular" || flow_pattern=="slug"

            major_axis_length = diameter_bubble(count+1,1);
            minor_axis_length = diameter_bubble(count+1,2);
            cur_min_dist_wall = min_dist_wall + major_axis_length/2;

            while attempts<max_attempts && count<1
                xy = (mvnrnd([0, 0], 0.01*eye(2), 1)); 
                candidate = xy;
                euler_angles = 0;

                if (candidate(1))^2 + (candidate(2))^2 <= (1-cur_min_dist_wall)^2
                    count = count + 1;
                    samples(count, 1:2) = candidate;
                    samples(count, 3) = major_axis_length;
                    samples(count, 4) = minor_axis_length;
                    samples(count, 5) = euler_angles;
                    attempts = 0;
                    bc=round(candidate.*Nx)+[Nx,Nx];
                    bx=bc(1);by=bc(2);
                    r_maj=round(major_axis_length*Nx/2);
                    r_min=round(minor_axis_length*Nx/2);
                    Q = [cos(euler_angles), -sin(euler_angles);sin(euler_angles),cos(euler_angles)];
                    for xx=max(1,bx-r_maj-1):min(2*Nx,bx+r_maj+1)
                        for yy=max(1,by-r_maj-1):min(2*Nx,by+r_maj+1)
                            bubble_mask_relative = [xx, yy]-[bx, by];
                            if bubble_mask_relative*Q*diag([1/r_maj^2,1/r_min^2])*transpose(Q)*transpose(bubble_mask_relative)<=1
                                c_mask(xx, yy) = 1;
                            end
                        end
                    end
                else
                    attempts = attempts+1;    
                end
            end
        end

        samples_solid=zeros(m,2);
        count_s=0;
        attempts_solid=0;
        while count_s<m && attempts_solid<max_attempts_solid
            xy = [2*rand-1, 2*rand-1];
            if flow_pattern=="annular" || flow_pattern=="slug"
                if xy(1)^2+xy(2)^2 <= (1-glass_radius-min_dist_wall)^2
                    if (xy(1)-samples(1,1))^2+(xy(2)-samples(1,2))^2>=min_dist+samples(1,3)/2+glass_radius
                        if count_s==0
                            count_s = count_s+1;
                            samples_solid(count_s,:)=xy;
                        else
                            dists = sqrt(sum((samples_solid(1:count_s,:)-xy).^2, 2));
                            if all(dists>=min_dist+2*glass_radius)
                                count_s=count_s+1;
                                samples_solid(count_s,:)=xy;
                                attempts_solid=0;
                            end
                        end
                    end
                end
            else
                if xy(1)^2+xy(2)^2 <= (1-glass_radius-min_dist_wall)^2
                    if count_s==0
                        count_s = count_s+1;
                        samples_solid(count_s,:)=xy;
                    else
                        dists = sqrt(sum((samples_solid(1:count_s,:)-xy).^2, 2));
                        if all(dists>=min_dist+2*glass_radius)
                            count_s=count_s+1;
                            samples_solid(count_s,:)=xy;
                            attempts_solid=0;
                        end
                    end
                end
            end
            attempts_solid=attempts_solid+1;
        end
        if count_s<m
            samples_solid = -1;
            fprintf("Number of particles: %d\n", count_s);
        else
            samples_solid = samples_solid*inner_radius;
            fprintf("Number of particles: %d\n", count_s);
        end

        for mm=1:count_s
            c_mask=c_mask | makeDisc(2*Nx,2*Nx,...
            round(samples_solid(mm,1)/inner_radius*Nx)+Nx,...
            round(samples_solid(mm,2)/inner_radius*Nx)+Nx,...
            round(glass_radius*Nx));
        end

        while count < num_bubble && attempts < max_attempts

            major_axis_length = diameter_bubble(count+1,1);
            minor_axis_length = diameter_bubble(count+1,2);
            cur_min_dist_wall = min_dist_wall + major_axis_length/2;
            % x, yはガウス分布、zは[-1,1]の一様分布からサンプリング
            % if major_axis_length >= 5.5e-3/inner_radius
            %     xy = (mvnrnd([0, 0], 0.02*eye(2), 1)); % 2x1ベクトル
            % elseif major_axis_length < 5.5e-3/inner_radius && major_axis_length >= 3e-3/inner_radius && flow_pattern=="bubble"
            %     xy = (mvnrnd([0, 0], 0.1*eye(2), 1));
            % else
            %     xy = [2*rand-1, 2*rand-1];
            % end
            xy = [2*rand-1, 2*rand-1];
            candidate = xy;
            euler_angles = 2*pi*rand;
            if (candidate(1))^2 + (candidate(2))^2 <= (1-cur_min_dist_wall)^2
                overlapping = 0;
                bc=round(candidate.*Nx)+[Nx,Nx];
                bx=bc(1);by=bc(2);
                r_maj=round((major_axis_length/2+min_dist)*Nx);
                r_min=round((minor_axis_length/2+min_dist)*Nx);
                Q = [cos(euler_angles), -sin(euler_angles);sin(euler_angles),cos(euler_angles)];
                xx=max(1,bx-r_maj-1);
                yy=max(1,by-r_maj-1);
                while overlapping==0 & xx<=min(2*Nx,bx+r_maj+1)
                    while overlapping==0 & yy<=min(2*Nx,by+r_maj+1)
                        bubble_mask_relative = [xx, yy]-[bx, by];
                        if bubble_mask_relative*Q*diag([1/r_maj^2,1/r_min^2])*transpose(Q)*transpose(bubble_mask_relative)<=1
                            if c_mask(xx,yy)==1
                                overlapping = 1;
                            end
                        end
                        yy = yy+1;
                    end
                    xx=xx+1;
                    yy=max(1,by-r_maj-1);
                end
                if ~overlapping
                    count = count + 1;
                    samples(count, 1:2) = candidate;
                    samples(count, 3) = major_axis_length;
                    samples(count, 4) = minor_axis_length;
                    samples(count, 5) = euler_angles;
                    attempts = 0;
                    bc=round(candidate.*Nx)+[Nx,Nx];
                    bx=bc(1);by=bc(2);
                    r_maj=round(major_axis_length*Nx/2);
                    r_min=round(minor_axis_length*Nx/2);
                    Q = [cos(euler_angles), -sin(euler_angles);sin(euler_angles),cos(euler_angles)];
                    for xx=max(1,bx-r_maj-1):min(2*Nx,bx+r_maj+1)
                        for yy=max(1,by-r_maj-1):min(2*Nx,by+r_maj+1)
                            bubble_mask_relative = [xx, yy]-[bx, by];
                            if bubble_mask_relative*Q*diag([1/r_maj^2,1/r_min^2])*transpose(Q)*transpose(bubble_mask_relative)<=1
                                c_mask(xx, yy) = 1;
                            end
                        end
                    end
                end       
            end
            attempts = attempts+1;
        end
        if count == num_bubble
            samples_ans = samples * inner_radius;
            samples_ans(:,5) = samples_ans(:,5)/inner_radius;
            break;
        end
        if count < num_bubble
            attempts_radius = attempts_radius + 1;
        end
        if attempts_radius == max_attempts_radius
            samples_ans = -1;
            break;
        end
    end
    C = {samples_ans,samples_solid};
end
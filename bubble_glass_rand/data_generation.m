function signalgen_module_all()
    config_file = 'config.json';

    config = jsondecode(fileread(config_file));

    save_dir = config.save_path;
    [filepath, filename, ext] = fileparts(mfilename('fullpath'));

    dataset_name = strsplit(filepath, '/');
    dataset_name = string(dataset_name(end));
    % dataset_name = glass_rand
    save_dir = fullfile(save_dir, dataset_name);
    disp(save_dir);

    if ~exist(save_dir)
        mkdir(save_dir);
    end

    num_repeat = config.simulation.num_dataset;
    for case_num = 1:16
        save_path = sprintf("%s/dataset20_case%d",save_dir,case_num);
        if ~exist(save_path)
            mkdir(save_path);
        end
        disp(save_path);
        loc_dir=fullfile(save_path,"location_seed");
        data_dir=fullfile(save_path,"data");
        logs_dir=fullfile(save_path,"logs");
        gnu_dir =fullfile(save_path,"gnu_csv");
        if ~exist(loc_dir)
            mkdir(loc_dir);
        end
        if ~exist(data_dir)
            mkdir(data_dir);
        end
        if ~exist(logs_dir)
            mkdir(logs_dir);
        end
        if ~exist(gnu_dir)
            mkdir(gnu_dir)
        end
        gd=config.glass.diameter*1e3/2;
        bubble_num = case_num;
        encoded=jsonencode(config);
        fid=fopen(fullfile(save_path,config_file), 'w');
        fprintf(fid,'%s',encoded);
        fclose(fid);
        i=2;
        while i <=11
            glass_num = randi([1,6]);
            C = pipe_location_gen(glass_num, bubble_num, config_file);
            samples=C{1};samples_solid=C{2};
            if samples_solid ~= -1
                if samples ~= -1
                    loc_path=fullfile(loc_dir,sprintf("location%d.csv",i));
                    sample_table=array2table(samples,'VariableNames', {'X', 'Y', 'RM', 'Rm', 'Theta'});
                    writetable(sample_table, loc_path);
                    loc_path_solid=fullfile(loc_dir,sprintf("location_solid%d.csv",i));
                    sample_table=array2table(samples_solid,'VariableNames', {'X', 'Y'});
                    writetable(sample_table, loc_path_solid);
                    if i==1
                        loctable = readtable(sprintf("%s/location%d.csv", loc_dir, i));
                        loc = table2array(loctable);
                        loc(:,1:4)=loc(:,1:4)*1e3;
                        loc(:,3:4)=loc(:,3:4)/2;
                        sz=size(loc);

                        fid = fopen(sprintf("%s/gnu_text%d.txt", loc_dir, i), "w");
                        fprintf(fid,'plot ');
                        for ii=1:sz(1)
                            x0=loc(ii,1);y0=loc(ii,2);aa=loc(ii,3);bb=loc(ii,4);th=loc(ii,5);
                            fprintf(fid,'(%.5g)*((%.5g)*cos(t))-(%.5g)*((%.5g)*sin(t))+(%.5g),\\\n',...
                                cos(th),aa,sin(th),bb,x0);
                            fprintf(fid,'(%.5g)*((%.5g)*cos(t))+(%.5g)*((%.5g)*sin(t))+(%.5g)\\\n',...
                                sin(th),aa,cos(th),bb,y0);
                            fprintf(fid,'lw 3 lc "white",\\\n');
                        end
                        loctable = readtable(sprintf("%s/location%d.csv", loc_dir, i));
                        loc = table2array(loctable);
                        loc=loc*1e3;
                        sz=size(loc);
                        for ii=1:sz(1)
                            x0=loc(ii,1);y0=loc(ii,2);
                            fprintf(fid,'(%.5g)*cos(t)+(%.5g),(%.5g)*sin(t)+(%.5g)\\\n',...
                                gd,x0,gd,y0);
                            fprintf(fid,'lw 3 lc "black",\\\n');
                        end
                        fclose(fid);
                    end
                    kwavesim(config_file,save_path,loc_path,loc_path_solid,i);
                    i = i+1;
                end
            end
        end
    end
function signalgen_module_all1()
    % This function reads config.json and sequentially runs signalgen_module for all location1.csv, location2.csv, ... in the specified directory.
    
        % Path to config file
        config_file = 'config.json';
    
        % Get location_seedfiles_path from config.json
        config = jsondecode(fileread(config_file));
        location_dir = config.location_seedfiles_path1;
    
        % Get all location*.csv files
        files = dir(fullfile(location_dir, 'location*.csv'));
        fprintf('--- %d location*.csv files found. ---\n', length(files));
        if isempty(files)
            error('No location*.csv files found in %s.', location_dir);
        end
        % Delete all files in the save_data_path directory before data generation
        save_data_path = config.save_data_path1;
        save_logs_path = config.save_logs_path1;
    
        % Delete all files in the save_data_path directory before data generation
        if exist(save_data_path, 'dir')
            files_to_delete = dir(fullfile(save_data_path, '*'));
            for k = 1:length(files_to_delete)
                fname = files_to_delete(k).name;
                % Skip '.' and '..'
                if ~strcmp(fname, '.') && ~strcmp(fname, '..')
                    fpath = fullfile(save_data_path, fname);
                    if isfile(fpath)
                        delete(fpath); % Delete file
                    elseif isfolder(fpath)
                        % If there are subfolders, remove them recursively
                        rmdir(fpath, 's');
                    end
                end
            end
            fprintf('All files in %s have been deleted before data generation.\n', save_data_path);
        else
            % If the directory does not exist, create it
            mkdir(save_data_path);
            fprintf('Directory %s did not exist and was created.\n', save_data_path);
        end
    
        % Delete all files in the save_logs_path directory before data generation
        if exist(save_logs_path, 'dir')
            files_to_delete = dir(fullfile(save_logs_path, '*'));
            for k = 1:length(files_to_delete)
                fname = files_to_delete(k).name;
                % Skip '.' and '..'
                if ~strcmp(fname, '.') && ~strcmp(fname, '..')
                    fpath = fullfile(save_logs_path, fname);
                    if isfile(fpath)
                        delete(fpath); % Delete file
                    elseif isfolder(fpath)
                        % If there are subfolders, remove them recursively
                        rmdir(fpath, 's');
                    end
                end
            end
            fprintf('All files in %s have been deleted before data generation.\n', save_logs_path);
        else
            % If the directory does not exist, create it
            mkdir(save_logs_path);
            fprintf('Directory %s did not exist and was created.\n', save_logs_path);
        end
        % Run signalgen_module for each CSV file
        % Initialize (delete all files and folders) in save_full_path directory before data generation
        save_full_path = config.save_full_path1;
        if exist(save_full_path, 'dir')
            files_to_delete = dir(fullfile(save_full_path, '*'));
            for k = 1:length(files_to_delete)
                fname = files_to_delete(k).name;
                % Skip '.' and '..'
                if ~strcmp(fname, '.') && ~strcmp(fname, '..')
                    fpath = fullfile(save_full_path, fname);
                    if isfile(fpath)
                        delete(fpath); % Delete file
                    elseif isfolder(fpath)
                        % Remove subfolders recursively
                        rmdir(fpath, 's');
                    end
                end
            end
            fprintf('All files and folders in %s have been deleted before data generation.\n', save_full_path);
        else
            % If the directory does not exist, create it
            mkdir(save_full_path);
            fprintf('Directory %s did not exist and was created.\n', save_full_path);
        end
    
        % Copy config2d.json file to save_full_path
        config_file_src = config_file; % assuming config_file is the path to config.json
        config_file_dst = fullfile(save_full_path, 'config.json');
        copyfile(config_file_src, config_file_dst);
        fprintf('config.json has been copied to %s.\n', save_full_path);
    
        % Copy location_seed folder to save_full_path
        location_seed_src = config.location_seedfiles_path1;
        [~, location_seed_folder] = fileparts(location_seed_src);
        location_seed_dst = fullfile(save_full_path, location_seed_folder);
        if exist(location_seed_dst, 'dir')
            rmdir(location_seed_dst, 's');
        end
        device_num = 1;
        copyfile(location_seed_src, location_seed_dst);
        fprintf('location_seed folder has been copied to %s.\n', save_full_path);
        for i = 1:length(files)
            location_csv = fullfile(location_dir, files(i).name);
            fprintf('--- Processing %s ---\n', location_csv);
            simulation_execution(config_file, location_csv, ...
            device_num, save_full_path);
        end
    end
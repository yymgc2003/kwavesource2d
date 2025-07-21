function signalgen_module_all()
% This function reads config.json and sequentially runs signalgen_module for all location1.csv, location2.csv, ... in the specified directory.

    % Path to config file
    config_file = 'config.json';

    % Get location_seedfiles_path from config.json
    config = jsondecode(fileread(config_file));
    location_dir = config.location_seedfiles_path;

    % Get all location*.csv files
    files = dir(fullfile(location_dir, 'location*.csv'));
    fprintf('--- %d location*.csv files found. ---\n', length(files));
    if isempty(files)
        error('No location*.csv files found in %s.', location_dir);
    end
    % Delete all files in the save_data_path directory before data generation
    save_data_path = config.save_data_path;
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
    % Run signalgen_module for each CSV file
    for i = 1:length(files)
        location_csv = fullfile(location_dir, files(i).name);
        fprintf('--- Processing %s ---\n', location_csv);
        simulation_execution(config_file, location_csv);
    end
end
function generate_spaced_samples(number_samples, save_path)
    % Generate n samples of 3D points with minimum distance constraint
    % Input: number_samples - number of samples to generate
    %        save_path - path where to save the CSV file
    % Output: CSV file with coordinates (x, y, z) for each sample
    
    % Set minimum distance threshold between any two points
    min_distance = 2.0;
    
    % Initialize storage for coordinates
    coordinates = zeros(number_samples, 3);
    
    % Generate first point randomly
    coordinates(1, :) = randn(1, 3);
    
    % Generate remaining points with distance constraint
    for i = 2:number_samples
        max_attempts = 1000;
        attempts = 0;
        valid_point_found = false;
        
        while attempts < max_attempts && ~valid_point_found
            % Generate candidate point from normal distribution
            candidate = randn(1, 3);
            
            % Check distance to all existing points
            valid_point_found = true;
            for j = 1:(i-1)
                distance = norm(candidate - coordinates(j, :));
                if distance < min_distance
                    valid_point_found = false;
                    break;
                end
            end
            
            attempts = attempts + 1;
        end
        
        if valid_point_found
            coordinates(i, :) = candidate;
        else
            error('Could not find valid point for sample %d after %d attempts', i, max_attempts);
        end
    end
    
    % Create table for CSV output
    sample_table = table(coordinates(:, 1), coordinates(:, 2), coordinates(:, 3), ...
        'VariableNames', {'x', 'y', 'z'});
    
    % Save to CSV file
    filename = sprintf('spaced_samples_%d.csv', number_samples);
    full_filename = fullfile(save_path, filename);
    writetable(sample_table, full_filename);
    
    fprintf('Generated %d samples with minimum distance %.2f\n', number_samples, min_distance);
    fprintf('Saved to file: %s\n', full_filename);
end

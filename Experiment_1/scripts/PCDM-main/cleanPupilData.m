function in_cleaned = cleanPupilData(in)
    % Loop through each trial/run in the cell array
    for i = 1:length(in.pupilArea)
        data = in.pupilArea{i};
        
        % Check if data is all-NaN, which can crash fillmissing
        if all(isnan(data))
            % Handle this case, e.g., set to zero or keep as-is if later code handles it
            fprintf('Warning: Trial %d is all NaNs.\n', i);
        else
            % Use linear interpolation to fill missing values
            % 'fillmissing' is available in recent MATLAB versions
            % For older versions, you might need a custom interpolation function
            data_filled = fillmissing(data, 'linear', 'EndValues', 'nearest');
            in.pupilArea{i} = data_filled;
        end
    end
    in_cleaned = in;
end
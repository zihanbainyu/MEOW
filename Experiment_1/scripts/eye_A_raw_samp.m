clear;
clc;
close all;
dbstop if error;

%==========================================================================
% SCRIPT 1: RAW ASC SAMPLE PARSER
%==========================================================================
fprintf('--- starting script 1: raw asc parser ---\n');

subj_id = 501;
base_dir = '..';
subj_folder = sprintf('sub%03d', subj_id); 
results_dir = fullfile(base_dir, 'data', subj_folder);
asc_file = fullfile(results_dir, sprintf('%03d_2_1.asc', subj_id));
output_file = fullfile(results_dir, sprintf('%03d_2_1_raw_samples.mat', subj_id));

% --- 1. open the file ---
fid = fopen(asc_file);
if fid == -1, error('asc file not found: %s', asc_file); end
fprintf('file opened. parsing... \n');

% --- 2. pre-allocate (this is the key to not crashing) ---
% a 10-minute 1000hz file has 600,000 samples. 
% we'll set a 5-million-sample buffer (overkill is safe)
n_buffer = 5000000;
raw_timestamps = zeros(n_buffer, 1);
raw_xpos = zeros(n_buffer, 1);
raw_ypos = zeros(n_buffer, 1);
raw_pupil = zeros(n_buffer, 1);

% for events
n_event_buffer = 5000;
event_timestamps = zeros(n_event_buffer, 1);
event_messages = cell(n_event_buffer, 1);

sample_count = 0;
event_count = 0;
tline = fgetl(fid);

% --- 3. read the file line by line ---
while ischar(tline)
    
    if isempty(tline)
        tline = fgetl(fid);
        continue;
    end
    
    % if the first character is a number, it's a sample
    if tline(1) >= '0' && tline(1) <= '9'
        sample_count = sample_count + 1;
        % use sscanf for speed. it's built for this.
        % format: timestamp, x, y, pupil, (skip rest)
        try
            sample_data = sscanf(tline, '%d %f %f %d %*s');
            raw_timestamps(sample_count) = sample_data(1);
            raw_xpos(sample_count) = sample_data(2);
            raw_ypos(sample_count) = sample_data(3);
            raw_pupil(sample_count) = sample_data(4);
        catch
            % a sample line was malformed, or a blink ('.')
            % we set to nan
            sample_count = sample_count + 1;
            raw_timestamps(sample_count) = sscanf(tline, '%d %*s');
            raw_xpos(sample_count) = nan;
            raw_ypos(sample_count) = nan;
            raw_pupil(sample_count) = nan;
        end
        
    % if it starts with 'MSG', it's an event
    elseif startsWith(tline, 'MSG')
        event_count = event_count + 1;
        parts = split(tline); % split is fine for slow events
        event_timestamps(event_count) = str2double(parts{2});
        event_messages{event_count} = tline;
    end
    
    tline = fgetl(fid);
end
fclose(fid);

% --- 4. truncate and save ---
fprintf('parsing complete. read %d samples and %d events.\n', sample_count, event_count);
raw_samples.timestamps = raw_timestamps(1:sample_count);
raw_samples.xpos_px = raw_xpos(1:sample_count);
raw_samples.ypos_px = raw_ypos(1:sample_count);
raw_samples.pupil_raw = raw_pupil(1:sample_count);

raw_samples.events.timestamps = event_timestamps(1:event_count);
raw_samples.events.messages = event_messages(1:event_count);

fprintf('saving raw .mat file to: %s\n', output_file);
save(output_file, 'raw_samples', '-v7.3'); % -v7.3 is for >2gb files
fprintf('--- script 1 complete ---\n');
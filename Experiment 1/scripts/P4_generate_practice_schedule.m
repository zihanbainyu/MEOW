
%==========================================================================
%                 Generate Practice Schedule
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================


%==================================================================
% Practice encoding schedule
%==================================================================
clear; clc; rng('shuffle');

practice_dir = fullfile('..','stimulus','stim_selected','practice');
out_dir = fullfile('..','subj_setup');

% Load practice targets, lures, and foils
targ_files = dir(fullfile(practice_dir,'prac_*_targ_l2.png'));
lure_files = dir(fullfile(practice_dir,'prac_*_lure_l2.png'));
foil_files = dir(fullfile(practice_dir,'prac_foil_*.png'));

% Convert names to strings
targ_names = string({targ_files.name});
lure_names = string({lure_files.name});
foil_names = string({foil_files.name});


seq = { ...
    targ_names(1), "repeat", "target", "none"; ...
    targ_names(1), "repeat", "target", "j"; ...
    foil_names(3), "new", "foil", "none"; ...
    foil_names(4), "new", "foil", "none"; ...
    targ_names(2), "similar", "target", "none"; ...
    lure_names(2), "similar", "lure", "k"; ...
};

col_names = {'stimulus_id','condition','role','correct_response','block','nback_target_id','fix_duration'};
variable_types = {'string','string','string','string','string','string','double'}; 
practice_schedule = table('Size',[0,7], 'VariableTypes', variable_types, ...
                          'VariableNames', col_names);

block_num = "1";

for i = 1:size(seq,1)
    stim_id = seq{i,1};
    cond = seq{i,2};
    role = seq{i,3};
    correct_resp = seq{i,4};
    
    if i == 1
        nback = "";
    else
        nback = seq{i-1,1};
    end
    
    % jittered fixation 0.5-1.0s
    fix_dur = 0.75 + (rand*0.5 - 0.25);
    
    % Append row
    new_row = table(stim_id, cond, role, correct_resp, ...
                block_num, nback, fix_dur, ...
                'VariableNames', col_names);
    practice_schedule = [practice_schedule; new_row];
end

% Save
save(fullfile(out_dir,'practice_encoding_schedule.mat'),'practice_schedule');

%==================================================================
% Practice memory schedule
%==================================================================
clear; clc; rng('shuffle');

practice_dir = fullfile('..','stimulus','stim_selected','practice');
out_dir = fullfile('..','subj_setup');

% Load practice targets, lures, and foils
targ_files = dir(fullfile(practice_dir,'prac_*_targ_l2.png'));
lure_files = dir(fullfile(practice_dir,'prac_*_lure_l2.png'));
foil_files = dir(fullfile(practice_dir,'prac_foil_*.png'));

% Convert names to strings
targ_names = string({targ_files.name});
lure_names = string({lure_files.name});
foil_names = string({foil_files.name});

seq = { ...
    % repeat
    targ_names(3), "repeat", "target", "none"; ...
    foil_names(3), "repeat", "foil", "none"; ...
    targ_names(3), "repeat", "target", "j"; ...
    % new
    foil_names(4), "new", "foil", "none"; ...
    foil_names(5), "new", "foil", "none"; ...
    foil_names(6), "new", "foil", "none"; ...
    % similar
    targ_names(4), "similar", "target", "none"; ...
    foil_names(7), "similar", "foil", "none"; ...
    lure_names(4), "similar", "lure", "k"; ...
    % new
    foil_names(8), "new", "foil", "none"; ...
    foil_names(9), "new", "foil", "none"; ...
    foil_names(10), "new", "foil", "none"; ...
    % similar
    targ_names(5), "similar", "target", "none"; ...
    foil_names(11), "similar", "foil", "none"; ...
    lure_names(5), "similar", "lure", "k"; ...
    % new
    foil_names(12), "new", "foil", "none"; ...
    foil_names(13), "new", "foil", "none"; ...
    foil_names(14), "new", "foil", "none"; ...
    % repeat
    targ_names(6), "repeat", "target", "none"; ...
    foil_names(15), "repeat", "foil", "none"; ...
    targ_names(6), "repeat", "target", "j"; ...
};

col_names = {'stimulus_id','condition','role','correct_response','block','nback_target_id','fix_duration'};
variable_types = {'string','string','string','string','string','string','double'}; 
practice_schedule = table('Size',[0,7], 'VariableTypes', variable_types, ...
                          'VariableNames', col_names);

block_num = "1";

for i = 1:size(seq,1)
    stim_id = seq{i,1};
    cond = seq{i,2};
    role = seq{i,3};
    correct_resp = seq{i,4};
    
    if i <= 2
        nback = "";
    else
        nback = seq{i-2,1};  % two trials ago
    end
    
    % jittered fixation 0.5–1.0 s
    fix_dur = 0.75 + (rand * 0.5 - 0.25);
    
    % append row to schedule
    new_row = table(stim_id, cond, role, correct_resp, ...
                    block_num, nback, fix_dur, ...
                    'VariableNames', col_names);
    
    practice_schedule = [practice_schedule; new_row];
end

% Save
save(fullfile(out_dir,'practice_memory_schedule.mat'),'practice_schedule');

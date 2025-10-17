% --- Keyboard Sanity Check Script ---
sca; clear;
fprintf('STARTING KEYBOARD TEST...\n');
fprintf('Listening ONLY to device 8.\n');
fprintf('Press the SPACEBAR to continue or ESCAPE to quit.\n\n');

% Basic setup
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1);
screenNumber = max(Screen('Screens'));
[window, ~] = PsychImaging('OpenWindow', screenNumber, 128); % Open gray window
p.keys.device = 8; % Hard-code the device index

% Define keys
KbName('UnifyKeyNames');
spaceKey = KbName('space');
escapeKey = KbName('escape');

% The Test Loop
while true
    [keyIsDown, ~, keyCode] = KbCheck(p.keys.device); % IMPORTANT: Check device 8
    
    if keyIsDown
        if keyCode(spaceKey)
            fprintf('SUCCESS! Spacebar press detected from device %d.\n', p.keys.device);
            break; % Exit loop
        elseif keyCode(escapeKey)
            fprintf('ESCAPE detected. Quitting.\n');
            break; % Exit loop
        end
    end
end

sca; % Close screen
fprintf('TEST COMPLETE.\n');
function [pressed, firstPress] = resp_queue(action, p, keyList)
% KbQueue on every device in p.keys.queue_devices. KbQueueCreate does not
% accept the merged index -3 (it falls back to the default keyboard), so
% one queue per keyboard/keypad is run and their presses are merged.
% Devices whose queue cannot be created are skipped. Escape aborts on any
% flush or check.
persistent devs
pressed = false;
firstPress = zeros(1, 256);
esc = KbName(p.keys.quit);
switch action
    case 'create'
        devs = [];
        for d = p.keys.queue_devices
            try
                KbQueueCreate(d, keyList);
                devs(end+1) = d; %#ok<AGROW>
            catch
                fprintf(2, 'resp_queue: skipped device %d\n', d);
            end
        end
        if isempty(devs), error('resp_queue: no response device could be opened.'); end
    case 'start'
        for d = devs, KbQueueStart(d); end
    case 'flush'
        for d = devs
            [pr, fp] = KbQueueCheck(d);
            if pr && fp(esc) > 0
                error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
            end
            KbQueueFlush(d);
        end
    case 'check'
        for d = devs
            [pr, fp] = KbQueueCheck(d);
            if pr
                pressed = true;
                hit = fp > 0 & (firstPress == 0 | fp < firstPress);
                firstPress(hit) = fp(hit);
            end
        end
        if firstPress(esc) > 0
            error('USER_ABORT:ExperimentAborted', 'Experiment aborted by user.');
        end
    case 'release'
        for d = devs, KbQueueRelease(d); end
        devs = [];
end
end

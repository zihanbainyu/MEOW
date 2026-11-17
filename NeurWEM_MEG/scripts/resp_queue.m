function [pressed, firstPress] = resp_queue(action, p, keyList)
% KbQueue on every device in p.keys.queue_devices. KbQueueCreate does not
% accept the merged index -3 (it falls back to the default keyboard), so
% one queue per keyboard/keypad is run and their presses are merged.
pressed = false;
firstPress = zeros(1, 256);
devs = p.keys.queue_devices;
switch action
    case 'create'
        for d = devs, KbQueueCreate(d, keyList); end
    case 'start'
        for d = devs, KbQueueStart(d); end
    case 'flush'
        esc = KbName(p.keys.quit);
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
    case 'release'
        for d = devs, KbQueueRelease(d); end
end
end

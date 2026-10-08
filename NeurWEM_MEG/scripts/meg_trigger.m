function meg_trigger(p, code)
% Send an 8-bit event code to the MEG trigger channels (pulse of p.trig.pulse
% seconds, then reset to 0) and mirror it to the EyeLink as 'TRIG <code>' so
% both recordings share the same event stream. The port write itself is the
% function handle p.trig.write, configured in main.m.
if p.trig.on
    p.trig.write(code);
    WaitSecs(p.trig.pulse);
    p.trig.write(0);
end
if p.eyetracking == 1
    Eyelink('Message', 'TRIG %d', code);
end
end

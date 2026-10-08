function draw_photodiode(p, is_on)
% Photodiode patch in a screen corner: white on image frames, black otherwise,
% so the true on-screen stimulus onset is recorded in an MEG analog channel.
% Drawn before the Flip it belongs to. Disabled when p.pd.on is 0.
if ~p.pd.on, return; end
if is_on, col = p.colors.white; else, col = p.colors.black; end
Screen('FillRect', p.window, col, p.pd.rect);
end

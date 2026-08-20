function s = pstr(p)
% PSTR  Journal-style p-value string.
if p < .001, s = 'p < .001'; else, s = sprintf('p = %.3f', p); end
end

function Dc = nm_load_long(s)
% Load one subject's MST-back behavior as a struct of recoded per-trial tables
% with item number + lure level parsed, correctness scored with the house
% conventions (resp '1'=same/old, '2'=similar, none/NA=new; correct = corr_resp==resp_key).
% Returns [] if the subject has no concat.
D=sprintf('/Users/bai/Documents/GitHub/MEOW/NeurWEM_CN/data/sub%d',s);
f=sprintf('%s/sub%d_concat.mat',D,s);
if ~isfile(f), Dc=[]; return; end
L=load(f,'final_data_output'); o=L.final_data_output;
Dc=struct('subj',s);
Dc.one = rc(o.results_1_back_all);
Dc.two = rc(o.results_2_back_all);
if isfield(o,'results_mst'),         Dc.mst = rc(o.results_mst); else, Dc.mst=[]; end
if isfield(o,'results_recognition'), Dc.rec = rc(o.results_recognition); else, Dc.rec=[]; end
end

function T=rc(T)
% recode a results table: resp_key NA->none; correct flag; parse item + lure level
T.resp_key = cellstr(string(T.resp_key)); T.resp_key(strcmp(T.resp_key,'NA'))={'none'};
T.corr_resp= cellstr(string(T.corr_resp));
T.correct  = strcmp(T.corr_resp, T.resp_key);
sid=string(T.stim_id);
T.item = nan(height(T),1);
for i=1:height(T)
    tok=regexp(sid(i),'mst_(\d+)_','tokens','once');
    if ~isempty(tok), T.item(i)=str2double(tok{1}); end
end
lv=regexp(sid,'_(l\d|foil)\.png$','tokens','once');
T.lure=strings(height(T),1);
for i=1:height(T)
    if ~isempty(lv(i)) && ~isempty(lv{i}), T.lure(i)=string(lv{i}); end
end
if ~isnumeric(T.rt), T.rt=str2double(string(T.rt)); end
end

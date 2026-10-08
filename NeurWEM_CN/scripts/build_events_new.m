function build_events_new(subj)
TR=1.5;
base='/Users/bai/Documents/GitHub/MEOW/NeurWEM_CN/data';
out=fullfile('/Users/bai/Documents/GitHub/MEOW/NeurWEM_CN/events_staging',sprintf('sub-%d',subj));
if ~exist(out,'dir'), mkdir(out); end
S=load(fullfile(base,sprintf('sub%d',subj),sprintf('sub%d_concat.mat',subj)));
fdo=S.final_data_output;
nb1=write_nback(fdo.results_1_back_all,'1back',false,subj,out,TR);
nb2=write_nback(fdo.results_2_back_all,'2back',true,subj,out,TR);
nm=0;
if isfield(fdo,'results_mst')
    try
        tm=fdo.results_mst;
        if istable(tm) && height(tm)>0, nm=write_mst(tm,subj,out,TR); end
    catch, nm=0; end
end
fprintf('sub-%d: 1back %d trials, 2back %d trials, mst %d trials -> %s\n',subj,nb1,nb2,nm,out);
end

function n=write_nback(t,task,hasgoal,subj,out,TR)
blocks=unique(t.block); fprintf('  %s blocks: %s\n',task,mat2str(blocks(:)'));
n=0;
for bi=1:numel(blocks)
    b=blocks(bi); r=t(t.block==b,:);
    E=table();
    E.onset=r.stim_onset_tr*TR;
    E.duration=repmat(TR,height(r),1);
    E.trial_type=r.condition;
    E.stim_onset_tr=r.stim_onset_tr;
    E.fix_onset_tr=r.fix_onset_tr;
    E.response_time=r.rt;
    E.stim_file=r.stim_id;
    E.condition=r.condition;
    E.identity=r.identity;
    if hasgoal, E.goal=r.goal; end
    E.corr_resp=r.corr_resp;
    E.resp_key=r.resp_key;
    fn=fullfile(out,sprintf('sub-%d_task-%s_run-%02d_events.tsv',subj,task,b));
    writetable(E,fn,'FileType','text','Delimiter','\t');
    n=n+height(r);
end
end

function n=write_mst(t,subj,out,TR)
E=table();
E.onset=t.stim_onset_tr*TR;
E.duration=repmat(TR,height(t),1);
E.trial_type=t.trial_type;
E.stim_onset_tr=t.stim_onset_tr;
E.fix_onset_tr=t.fix_onset_tr;
E.response_time=t.rt;
E.stim_file=t.stim_id;
E.corr_resp=t.corr_resp;
E.resp_key=t.resp_key;
fn=fullfile(out,sprintf('sub-%d_task-mst_events.tsv',subj));
writetable(E,fn,'FileType','text','Delimiter','\t');
n=height(t);
end

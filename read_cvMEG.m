function [meg_data, trl_meg] = read_cvMEG(meg_file, save_path, params)
%prprocess_cvMEG Read conventional MEG data for benchmarking
% recordings. 
% Requires the following arguments:
% Path: containing save_path and meg_file
% Params: containing pre, post (pre- and poststim).

%%
pre = params.pre;
post = params.post;
filter = params.filter;

%% --- Read triggers ---
trl_meg = [];
cfg             = [];
cfg.datafile    = meg_file;
meg_raw         = ft_preprocessing(cfg);
meg_trig = find(contains(meg_raw.label,'STI101'));
trig = meg_raw.trial{1}(meg_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
smpl = 1:meg_raw.sampleinfo(2);
trl_meg(:,1) = smpl(trig)-pre*meg_raw.fsample;
trl_meg(:,2) = smpl(trig)+post*meg_raw.fsample;
trl_meg(:,3) = -pre*meg_raw.fsample;
trl_meg(:,4) = meg_raw.trial{1}(meg_trig,smpl(trig));

%% --- Read data ---
cfg = [];
cfg.dataset         = meg_file;
cfg.trl             = trl_meg;
%cfg.checkmaxfiler = 'no';
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = filter.hp_freq;
cfg.padding         = 2;
cfg.padtype         = 'data';
meg_data = ft_preprocessing(cfg);

%% Notch filter
cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = filter.notch;
meg_data_filt = ft_preprocessing(cfg,meg_data);

%% Save 
save(fullfile(save_path, 'meg_data'), 'meg_data_filt', 'trl_meg',"-v7.3"); disp('done');

end
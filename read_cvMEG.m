function [squid_cleaned, squideeg_cleaned] = read_cvMEG(squid_file, save_path, params)
%prprocess_cvMEG Read conventional MEG data for benchmarking
% recordings. 
% Requires the following arguments:
% Path: containing save_path and squid_file
% Params: containing pre, post (pre- and poststim).

%% --- Read triggers ---
trl_meg = [];
cfg             = [];
cfg.datafile    = squid_file;
squid_raw         = ft_preprocessing(cfg);
squid_trig = find(contains(squid_raw.label,'STI101'));
trig = squid_raw.trial{1}(squid_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
smpl = 1:squid_raw.sampleinfo(2);
trl_meg(:,1) = smpl(trig)-params.pre*squid_raw.fsample;
trl_meg(:,2) = smpl(trig)+params.post*squid_raw.fsample;
trl_meg(:,3) = -params.pre*squid_raw.fsample;
trl_meg(:,4) = squid_raw.trial{1}(squid_trig,smpl(trig));

%% MEG data filter & epoch
cfg = [];
cfg.datafile        = squid_file;
cfg.trl             = trl_meg;
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
cfg.hpinstabilityfix  = 'reduce';
cfg.padding         = params.pre + params.post + 2;
cfg.paddingtype     = 'data';
cfg.demean          = 'yes';
cfg.baselinewindow  = [-inf 0];
squid_epo = ft_preprocessing(cfg);

%% MEG 
cfg = [];
cfg.channel = squid_epo.label(find(~contains(squid_epo.label,'eeg')));
squid_cleaned = ft_selectdata(cfg, squid_epo);

% no bad channel detection since maxfilter already does that

% Reject jump trials
cfg = [];
cfg.channel = 'meg';
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,badtrl_squid_jump] = ft_badsegment(cfg, squid_cleaned);
squid_cleaned = ft_rejectartifact(cfg,squid_cleaned);

% Reject noisy trials
cfg = [];
cfg.channel = 'squidmag';
cfg.metric = 'std';
cfg.threshold = params.squidmag_std_threshold;
[cfg,badtrl_squidmag_std] = ft_badsegment(cfg, squid_cleaned);
squid_cleaned = ft_rejectartifact(cfg,squid_cleaned);

cfg = [];
cfg.channel = 'squidgrad';
cfg.metric = 'std';
cfg.threshold = params.squidgrad_std_threshold;
[cfg,badtrl_squidgrad_std] = ft_badsegment(cfg, squid_cleaned);
squid_cleaned = ft_rejectartifact(cfg,squid_cleaned);

%% EEG
cfg = [];
cfg.channel = squid_epo.label(find(~contains(squid_epo.label,'MEG')));
squideeg_cleaned = ft_selectdata(cfg, squid_epo);

% Reject bad channels
cfg = [];
cfg.trl = trl_meg;
squid_raw_epo = ft_redefinetrial(cfg,squid_raw);
cfg = [];
cfg.z_threshold = params.z_threshold;
cfg.corr_threshold = params.corr_threshold;
[badchs_opmeeg, badchs_squideeg_flat, badchs_squideeg_neighbors, badchs_squideeg_zmax, badtrl_squideeg_zmax] = eeg_badchannels(cfg,squid_raw_epo);
cfg = [];
cfg.channel = setdiff(squideeg_cleaned.label,badchs_opmeeg);
cfg.trials  = setdiff(1:length(squideeg_cleaned.trial),badtrl_squideeg_zmax); % remove bad trials
squideeg_cleaned = ft_selectdata(cfg, squideeg_cleaned);

% Reject jump trials
cfg = [];
cfg.channel = {'EEG*'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg, badtrl_squideeg_jump] = ft_badsegment(cfg, squideeg_cleaned);
squideeg_cleaned = ft_rejectartifact(cfg,squideeg_cleaned);

% Reject noisy trials
cfg = [];
cfg.channel = {'EEG*'};
cfg.metric = 'std';
cfg.threshold = params.eeg_std_threshold;
[cfg, badtrl_squideeg_std] = ft_badsegment(cfg, squideeg_cleaned);
squideeg_cleaned = ft_rejectartifact(cfg,squideeg_cleaned);

%% Save 
save(fullfile(save_path, [params.sub '_squideeg_badchs']), ...
    'badchs_squideeg_flat', ...
    'badchs_squideeg_neighbors', ...
    'badchs_squideeg_zmax' ,"-v7.3"); 

[~,idx]=ismember(squid_cleaned.sampleinfo,badtrl_squid_jump,'rows');
badtrl_squid_jump = find(idx);
[~,idx]=ismember(squid_cleaned.sampleinfo,badtrl_squidmag_std,'rows');
badtrl_squid_std = find(idx);
[~,idx]=ismember(squid_cleaned.sampleinfo,badtrl_squidgrad_std,'rows');
badtrl_squid_std = unique([badtrl_squid_std; find(idx)]);
save(fullfile(save_path, [params.sub '_squid_badtrls']), ...
    'badtrl_squid_jump', ...
    'badtrl_squid_std',"-v7.3"); 

[~,idx]=ismember(squideeg_cleaned.sampleinfo,badtrl_squideeg_jump,'rows');
badtrl_squideeg_jump = find(idx);
[~,idx]=ismember(squideeg_cleaned.sampleinfo,badtrl_squideeg_std,'rows');
badtrl_squideeg_std = find(idx);
save(fullfile(save_path, [params.sub '_squideeg_badtrls']), ...
    'badtrl_squideeg_jump', ...
    'badtrl_squideeg_std', ...
    'badtrl_squideeg_zmax',"-v7.3"); 

%save(fullfile(save_path, [params.sub '_squid_cleaned']), 'squid_cleaned',"-v7.3");
%save(fullfile(save_path, [params.sub '_squideeg_cleaned']), 'squideeg_cleaned',"-v7.3"); disp('done');

end
function [meg_cleaned, megeeg_cleaned] = read_cvMEG(meg_file, save_path, params)
%prprocess_cvMEG Read conventional MEG data for benchmarking
% recordings. 
% Requires the following arguments:
% Path: containing save_path and meg_file
% Params: containing pre, post (pre- and poststim).

%% --- Read triggers ---
trl_meg = [];
cfg             = [];
cfg.datafile    = meg_file;
meg_raw         = ft_preprocessing(cfg);
meg_trig = find(contains(meg_raw.label,'STI101'));
trig = meg_raw.trial{1}(meg_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
smpl = 1:meg_raw.sampleinfo(2);
trl_meg(:,1) = smpl(trig)-params.pre*meg_raw.fsample;
trl_meg(:,2) = smpl(trig)+params.post*meg_raw.fsample;
trl_meg(:,3) = -params.pre*meg_raw.fsample;
trl_meg(:,4) = meg_raw.trial{1}(meg_trig,smpl(trig));

%% MEG data filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.hpinstabilityfix  = 'reduce';
cfg.padding         = 1;
cfg.padtype         = 'data';
meg_epo = ft_preprocessing(cfg,meg_raw);

% Epoch
cfg = [];
cfg.trl             = trl_meg;
meg_epo = ft_redefinetrial(cfg,meg_epo);

% Notch filter
cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
meg_epo = ft_preprocessing(cfg,meg_epo);

%% MEG 
cfg = [];
cfg.channel = meg_epo.label(find(~contains(meg_epo.label,'eeg')));
meg_cleaned = ft_selectdata(cfg, meg_epo);

% no bad channel detection since maxfilter already does that

% Reject jump trials
cfg = [];
cfg.channel = 'meg';
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,badtrl_meg_jump] = ft_badsegment(cfg, meg_cleaned);
meg_cleaned = ft_rejectartifact(cfg,meg_cleaned);

% Reject noisy trials
cfg = [];
cfg.channel = 'megmag';
cfg.metric = 'std';
cfg.threshold = params.megmag_std_threshold;
[cfg,badtrl_megmag_std] = ft_badsegment(cfg, meg_cleaned);
meg_cleaned = ft_rejectartifact(cfg,meg_cleaned);

cfg = [];
cfg.channel = 'megplanar';
cfg.metric = 'std';
cfg.threshold = params.megplanar_std_threshold;
[cfg,badtrl_megplanar_std] = ft_badsegment(cfg, meg_cleaned);
meg_cleaned = ft_rejectartifact(cfg,meg_cleaned);

%% EEG
cfg = [];
cfg.channel = meg_epo.label(find(~contains(meg_epo.label,'MEG')));
megeeg_cleaned = ft_selectdata(cfg, meg_epo);

% Reject bad channels
cfg = [];
cfg.trl = trl_meg;
meg_raw_epo = ft_redefinetrial(cfg,meg_raw);
cfg = [];
cfg.z_threshold = params.z_threshold;
cfg.corr_threshold = params.corr_threshold;
[badchs_opmeeg, badchs_megeeg_flat, badchs_megeeg_neighbors, badchs_megeeg_zmax] = eeg_badchannels(cfg,meg_raw_epo);
cfg = [];
cfg.channel = setdiff(megeeg_cleaned.label,badchs_opmeeg);
megeeg_cleaned = ft_selectdata(cfg, megeeg_cleaned);

% Reject jump trials
cfg = [];
cfg.channel = {'EEG*'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg, badtrl_megeeg_jump] = ft_badsegment(cfg, megeeg_cleaned);
megeeg_cleaned = ft_rejectartifact(cfg,megeeg_cleaned);

% Reject noisy trials
cfg = [];
cfg.channel = {'EEG*'};
cfg.metric = 'std';
cfg.threshold = params.eeg_std_threshold;
[cfg, badtrl_megeeg_std] = ft_badsegment(cfg, megeeg_cleaned);
megeeg_cleaned = ft_rejectartifact(cfg,megeeg_cleaned);

%% Save 
save(fullfile(save_path, [params.sub '_megeeg_badchs']), ...
    'badchs_megeeg_flat', ...
    'badchs_megeeg_neighbors', ...
    'badchs_megeeg_zmax' ,"-v7.3"); 

save(fullfile(save_path, [params.sub '_meg_badtrls']), ...
    'badtrl_meg_jump', ...
    'badtrl_megmag_std', ...
    'badtrl_megplanar_std' ,"-v7.3"); 
save(fullfile(save_path, [params.sub '_megeeg_badtrls']), ...
    'badtrl_megeeg_jump', ...
    'badtrl_megeeg_std',"-v7.3"); 

%save(fullfile(save_path, [params.sub '_meg_cleaned']), 'meg_cleaned',"-v7.3");
%save(fullfile(save_path, [params.sub '_megeeg_cleaned']), 'megeeg_cleaned',"-v7.3"); disp('done');

end
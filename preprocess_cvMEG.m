function [meg_cleaned] = preprocess_cvMEG(data,save_path)
%prprocess_osMEG Preprocessing on-scalp MEG data for benchmarking
% recordings. Requires arguments:
% Path: containing save_path and meg_file
% Params: containing pre, post (pre- and poststim), lp_freq, hp_freq,
% bp_freq and notch filter frequencies (corresponding filters are only
% applied if the frequency is defined), n_comp and coh_cutoff (for 
% automated ICA), and ds_freq (downsampling frequency).

%% Reject bad channels/trials
% Jumps
cfg = [];
cfg.method = 'summary';
cfg.keepchannel = 'no'; 
cfg.keeptrials = 'no'; 
cfg.channel = {'megmag'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
jump_rej = ft_rejectvisual(cfg, data);

jump_badchs = jump_rej.cfg.badchannel(contains(jump_rej.cfg.badchannel,'bz'));
include_chs = setdiff(data.label,jump_badchs);
cfg = [];
cfg.channel = include_chs;
cfg.trials  = jump_rej.cfg.trials; % remove bad trials
data_badchs = ft_selectdata(cfg, data);

%% Visual rejection
cfg = [];
cfg.channel = 'megmag';
mag_data = ft_selectdata(cfg, data_badchs);
cfg.channel = 'meggrad';
grad_data = ft_selectdata(cfg, data_badchs);
cfg.channel = 'EEG';
eeg_data = ft_selectdata(cfg, data_badchs);

cfg = [];
cfg.method = 'summary';
cfg.keepchannel = 'no'; 
cfg.keeptrials = 'no'; 

mag_data_rej = ft_rejectvisual(cfg, mag_data);
mag_badchs = setdiff(mag_data.label,mag_data_rej.label);
grad_data_rej = ft_rejectvisual(cfg, grad_data);
grad_badchs = setdiff(grad_data.label,grad_data_rej.label);
eeg_data_rej = ft_rejectvisual(cfg, eeg_data);
eeg_badchs = setdiff(eeg_data.label,eeg_data_rej.label);

cfg = [];
cfg.channel = setdiff(data_badchs.label,[mag_badchs; eeg_badchs; grad_badchs]);
cfg.trials  = intersect(mag_data_rej.cfg.trials,intersect(grad_data_rej.cfg.trials, eeg_data_rej.cfg.trials));
meg_cleaned = ft_selectdata(cfg, data_badchs);

save(fullfile(save_path, 'meg_badchannels'), 'jump_badchs', 'mag_badchs','eeg_badchs'); disp('done');
clear temp

%% Browse raw
% cfg = [];
% cfg.continuous          = 'no';
% cfg.allowoverlap        = 'yes';
% cfg.viewmode            = 'vertical';
% cfg.preproc.demean      = 'yes';
% cfg.channel             = {'EEG', 'megmag'};
% cfg.ylim                = [-1e-12 1e-12];
% cfg.eogscale            = 1e-7;
% cfg.eegscale            = 1e-7;
% meg_artif = ft_databrowser(cfg, meg_cleaned);
% save(fullfile(save_path, 'meg_artifacts'), 'meg_artif', '-v7.3'); disp('done');
% 
% % Remove artifacts
% if isfield(meg_artif,'artfctdef')
%     cfg             = [];
%     cfg.artfctdef   = meg_artif.artfctdef;
%     meg_cleaned       = ft_rejectartifact(cfg, meg_cleaned);
% end

%% Save
save(fullfile(save_path, 'meg_cleaned'), 'meg_cleaned',"-v7.3"); disp('done');

end
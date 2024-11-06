function [opm_cleaned] = preprocess_osMEG(data,save_path)
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
cfg.channel = {'*bz'};
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


cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = 20;
jump_rej2 = ft_badsegment(cfg, data);
jump_badchs2 = jump_rej2.cfg.badchannel(contains(jump_rej2.cfg.badchannel,'bz'));
include_chs = setdiff(data.label,jump_badchs2);


%% Visual rejection
cfg = [];
cfg.channel = '*bz';
mag_data = ft_selectdata(cfg, data_badchs);
cfg.channel = 'EEG';
eeg_data = ft_selectdata(cfg, data_badchs);

cfg = [];
cfg.method = 'summary';
cfg.keepchannel = 'no'; 
cfg.keeptrials = 'no'; 

mag_data_rej = ft_rejectvisual(cfg, mag_data);
mag_badchs = setdiff(mag_data.label,mag_data_rej.label);
eeg_data_rej = ft_rejectvisual(cfg, eeg_data);
eeg_badchs = setdiff(eeg_data.label,eeg_data_rej.label);

cfg = [];
cfg.channel = setdiff(data_badchs.label,[mag_badchs; eeg_badchs]);
cfg.trials  = intersect(mag_data_rej.cfg.trials, eeg_data_rej.cfg.trials);
opm_cleaned = ft_selectdata(cfg, data_badchs);

save(fullfile(save_path, 'opm_badchannels'), 'jump_badchs', 'mag_badchs','eeg_badchs'); disp('done');
clear temp

%% Browse raw
% cfg = [];
% cfg.continuous          = 'no';
% cfg.allowoverlap        = 'yes';
% cfg.viewmode            = 'vertical';
% cfg.preproc.demean      = 'yes';
% cfg.channel             = {'EEG', '*bz'};
% cfg.ylim                = [-1e-12 1e-12];
% cfg.eogscale            = 1e-7;
% cfg.eegscale            = 1e-7;
% opm_artif = ft_databrowser(cfg, opm_cleaned);
% save(fullfile(save_path, 'opm_artifacts'), 'opm_artif', '-v7.3'); disp('done');
% 
% % Remove artifacts
% if isfield(opm_artif,'artfctdef')
%     cfg             = [];
%     cfg.artfctdef   = opm_artif.artfctdef;
%     %cfg.artfctdef.reject  = 'nan';
%     opm_cleaned       = ft_rejectartifact(cfg, opm_cleaned);
% end

%% Save
save(fullfile(save_path, 'opm_cleaned'), 'opm_cleaned',"-v7.3"); disp('done');

end
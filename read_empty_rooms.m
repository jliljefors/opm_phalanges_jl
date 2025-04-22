function read_empty_rooms(opm_file, squid_file, opm_chs, squid_chs, save_path, params)
%prprocess_osMEG Read on-scalp MEG data for benchmarking
% recordings and combine with auxiliary TRIUX data/EEG. 
% Requires the following arguments:
% Path: containing save_path and meg_file
% Params: containing pre, post (pre- and poststim), and ds_freq 
% (downsampling frequency).

%% OPM
ft_hastoolbox('mne', 1);

% Read triggers
cfg = [];
cfg.datafile        = opm_file;
cfg.coordsys        = 'dewar';
cfg.coilaccuracy    = 0;
cfg.channel         = opm_chs;
data_raw = ft_preprocessing(cfg);
n_smpl = round((params.pre+params.post)*data_raw.fsample);
n_trl = floor(data_raw.sampleinfo(2)/n_smpl-2);
if n_trl > 120
    n_trl = 120;
end
trl = zeros(n_trl,4);
trl(:,1) = data_raw.fsample + n_smpl*(0:(n_trl-1))';
trl(:,2) = data_raw.fsample + n_smpl*(1:n_trl)' - 1;
trl(:,3) = -params.pre*data_raw.fsample;
trl(:,4) = ones(length(trl(:,1)),1);

% Filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.hpinstabilityfix  = 'reduce';
data_epo = ft_preprocessing(cfg, data_raw);

cfg = [];
cfg.trl = trl;
data_epo = ft_redefinetrial(cfg,data_epo);

cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
data_epo = ft_preprocessing(cfg,data_epo);

% Resample 
cfg                 = [];
cfg.method          = 'resample';
cfg.resamplefs      = 1000;
empty_room_cleaned = ft_resampledata(cfg, data_epo);
empty_room_cleaned.sampleinfo = round(trl(:,1:2)/5);

cfg = [];
cfg.trl = trl;
cfg.z_threshold = params.z_threshold;
cfg.corr_threshold = params.corr_threshold;
[~, ~, ~, ~, ~, ~, badtrl_zmax] = opm_badchannels(cfg, data_raw);
cfg = [];
cfg.trials  = setdiff(1:length(empty_room_cleaned.trial),badtrl_zmax); % remove bad trials
empty_room_cleaned = ft_selectdata(cfg, empty_room_cleaned);
empty_room_cleaned.sampleinfo = round(trl(:,1:2)/5);
empty_room_cleaned.sampleinfo(:,2) = empty_room_cleaned.sampleinfo(:,2)-1;
clear trl;

% HFC
i_chs = find(contains(empty_room_cleaned.label,'bz'));
chs = empty_room_cleaned.label(i_chs);
ori = zeros(length(chs),3);
for i = 1:length(chs)
    i_chs_grad(i) = find(strcmp(empty_room_cleaned.grad.label,chs{i}));
    ori(i,:) = empty_room_cleaned.grad.chanori(i_chs_grad(i),:);
end
empty_room_cleaned.grad.M = eye(size(ori,1)) - ori*pinv(ori);
for i = 1:length(empty_room_cleaned.trial)
    empty_room_cleaned.trial{i}(i_chs,:) = params.sign * empty_room_cleaned.grad.M*empty_room_cleaned.trial{i}(i_chs,:);
end
% update grad
empty_room_cleaned.grad.tra(i_chs_grad,i_chs_grad) = empty_room_cleaned.grad.M * empty_room_cleaned.grad.tra(i_chs_grad,i_chs_grad); 

% Reject jump trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,badtrl_jump] = ft_badsegment(cfg, empty_room_cleaned);
empty_room_cleaned = ft_rejectartifact(cfg,empty_room_cleaned);

% Reject noisy trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'std';
cfg.threshold = params.opm_std_threshold;
[cfg,badtrl_std] = ft_badsegment(cfg, empty_room_cleaned);
empty_room_cleaned = ft_rejectartifact(cfg,empty_room_cleaned);

[~,idx]=ismember(empty_room_cleaned.sampleinfo,badtrl_jump,'rows');
badtrl_jump = find(idx);
[~,idx]=ismember(empty_room_cleaned.sampleinfo,badtrl_std,'rows');
badtrl_std = find(idx);
save(fullfile(save_path, [params.sub '_opm_badtrls']), ...
    'badtrl_jump', ...
    'badtrl_std', ...
    'badtrl_zmax' ,"-v7.3"); 

% Flip? 
if params.flip
    for i = 1:length(empty_room_cleaned.trial)
        empty_room_cleaned.trial{i}(i_chs,:) = -empty_room_cleaned.trial{i}(i_chs,:);
    end
end

% Save
save(fullfile(save_path, [params.sub '_opm_empty_room_cleaned']), 'empty_room_cleaned',"-v7.3");

clear empty_room_cleaned data_raw data_epo

%% SQUID
ft_hastoolbox('mne', 1);

% Read triggers
cfg             = [];
cfg.datafile    = squid_file;
cfg.channel         = squid_chs;
data_raw         = ft_preprocessing(cfg);
n_smpl = round((params.pre+params.post)*data_raw.fsample);
n_trl = floor(data_raw.sampleinfo(2)/n_smpl -2);
if n_trl > 120
    n_trl = 120;
end
trl = zeros(n_trl,4);
trl(:,1) = data_raw.fsample + n_smpl*(0:(n_trl-1))';
trl(:,2) = data_raw.fsample + n_smpl*(1:n_trl)' - 1;
trl(:,3) = -params.pre*data_raw.fsample;
trl(:,4) = ones(length(trl(:,1)),1);

% Filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.hpinstabilityfix  = 'reduce';
data_epo = ft_preprocessing(cfg, data_raw);

cfg = [];
cfg.trl = trl;
data_epo = ft_redefinetrial(cfg,data_epo);
clear trl 

cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
data_epo = ft_preprocessing(cfg,data_epo);

% Reject jump trials
cfg = [];
cfg.channel = 'meg';
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,badtrl_squid_jump] = ft_badsegment(cfg, data_epo);
empty_room_cleaned = ft_rejectartifact(cfg,data_epo);

% Reject noisy trials
cfg = [];
cfg.channel = 'squidmag';
cfg.metric = 'std';
cfg.threshold = params.squidmag_std_threshold;
[cfg,badtrl_squidmag_std] = ft_badsegment(cfg, empty_room_cleaned);
empty_room_cleaned = ft_rejectartifact(cfg,empty_room_cleaned);

cfg = [];
cfg.channel = 'squidgrad';
cfg.metric = 'std';
cfg.threshold = params.squidgrad_std_threshold;
[cfg,badtrl_squidgrad_std] = ft_badsegment(cfg, empty_room_cleaned);
empty_room_cleaned = ft_rejectartifact(cfg,empty_room_cleaned);

[~,idx]=ismember(empty_room_cleaned.sampleinfo,badtrl_squid_jump,'rows');
badtrl_squid_jump = find(idx);
[~,idx]=ismember(empty_room_cleaned.sampleinfo,badtrl_squidmag_std,'rows');
badtrl_squid_std = find(idx);
[~,idx]=ismember(empty_room_cleaned.sampleinfo,badtrl_squidgrad_std,'rows');
badtrl_squid_std = unique([badtrl_squid_std; find(idx)]);
save(fullfile(save_path, [params.sub '_squid_badtrls']), ...
    'badtrl_squid_jump', ...
    'badtrl_squid_std',"-v7.3"); 

% Save
save(fullfile(save_path, [params.sub '_squid_empty_room_cleaned']), 'empty_room_cleaned',"-v7.3");

clear empty_room_cleaned data_raw data_epo
end
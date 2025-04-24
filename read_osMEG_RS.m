function [] = read_osMEG_RS(opm_file, aux_file, opm_chs, save_path, params)
%prprocess_osMEG Read on-scalp MEG data for benchmarking
% recordings and combine with auxiliary TRIUX data/EEG. 
% Requires the following arguments:
% Path: containing save_path and meg_file
% Params: containing pre, post (pre- and poststim), and ds_freq 
% (downsampling frequency).

%% --- Read triggers ---
% AUX
cfg = [];
cfg.datafile        = aux_file;
aux_raw = ft_preprocessing(cfg);
aux_trig = contains(aux_raw.label,'STI101');
trig = aux_raw.trial{1}(aux_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trig = find(trig);
aux_period = aux_raw.time{1}(trig(end))-aux_raw.time{1}(trig(1));
n_smpl = round((params.pre+params.post)*aux_raw.fsample);
n_trl = floor((trig(end)-trig(1))/n_smpl-1);
trl_aux = zeros(n_trl,4);
trl_aux(:,1) = trig(1) + n_smpl*(0:(n_trl-1))';
trl_aux(:,2) = trig(1) + n_smpl*(1:n_trl)' - 1;
trl_aux(:,3) = -params.pre*aux_raw.fsample;
trl_aux(:,4) = ones(length(trl_aux(:,1)),1);

% OPM
cfg = [];
cfg.datafile        = opm_file;
cfg.coordsys        = 'dewar';
cfg.coilaccuracy    = 0;
cfg.channel         = opm_chs;
opm_raw = ft_preprocessing(cfg);
opm_trig = contains(opm_raw.label,'di');
trig = opm_raw.trial{1}(opm_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trig = find(trig);
opm_period = opm_raw.time{1}(trig(end))-opm_raw.time{1}(trig(1));
n_smpl = round((opm_period/aux_period)*(params.pre+params.post)*opm_raw.fsample);
n_trl = floor((trig(end)-trig(1))/n_smpl-1);
trl_opm = zeros(n_trl,4);
trl_opm(:,1) = trig(1) + n_smpl*(0:(n_trl-1))';
trl_opm(:,2) = trig(1) + n_smpl*(1:n_trl)' - 1;
trl_opm(:,3) = -params.pre*opm_raw.fsample;
trl_opm(:,4) = ones(length(trl_opm(:,1)),1);

% Check if uneven amount of trial. If so assume error in beginning.
if size(trl_aux,1) > size(trl_opm,1)
    trl_aux = trl_aux((end-size(trl_opm,1)+1):end,:);
elseif size(trl_aux,1) < size(trl_opm,1)
    trl_opm = trl_opm((end-size(trl_aux,1)+1):end,:);
end
if trl_aux(:,4) ~= trl_opm(:,4) % Throw error if trials don't match.
    disp()
    error('events do not match')
end

%% AUX data filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.hpinstabilityfix  = 'reduce';
aux_epo = ft_preprocessing(cfg, aux_raw);

cfg = [];
cfg.trl = trl_aux;
aux_epo = ft_redefinetrial(cfg,aux_epo);

cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
aux_epo = ft_preprocessing(cfg,aux_epo);

%% OPM data filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.hpinstabilityfix  = 'reduce';
opm_epo = ft_preprocessing(cfg,opm_raw);

cfg = [];
cfg.trl             = trl_opm;
opm_epo = ft_redefinetrial(cfg,opm_epo);

cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
opm_epo = ft_preprocessing(cfg,opm_epo);

%% --- Resample --- 
cfg            = [];
cfg.time = aux_epo.time;
cfg.detrend    = 'no';
opm_epo_ds = ft_resampledata(cfg, opm_epo);

%% Combine data
EOG_channels = find(contains(aux_epo.label,'EOG'));
ECG_channels = find(contains(aux_epo.label,'ECG'));
EEG_channels = find(contains(aux_epo.label,'EEG'));
MISC_channels = find(contains(aux_epo.label,'MISC'));
TRIG_channels = find(contains(aux_epo.label,'STI101'));
include_channels = [EOG_channels; ECG_channels; EEG_channels; MISC_channels; TRIG_channels];

comb = opm_epo_ds; 
comb.elec = aux_epo.elec;
comb.time = aux_epo.time;
comb.label = [comb.label; aux_epo.label(include_channels)];
comb.hdr.label = comb.label;
comb.hdr.nChans = comb.hdr.nChans + length(include_channels);
comb.hdr.chantype = [comb.hdr.chantype; aux_epo.hdr.chantype(include_channels)];
comb.hdr.chanunit = [comb.hdr.chanunit; aux_epo.hdr.chanunit(include_channels)];
comb.sampleinfo = aux_epo.sampleinfo;
n_smpl = size(aux_epo.trial{1},2);
for i = 1:length(comb.trial)
    comb.trial{i} = [comb.trial{i}(:,1:n_smpl); aux_epo.trial{i}(include_channels,:)]; 
end

%% OPM 
cfg = [];
cfg.channel = comb.label(find(~contains(comb.label,'eeg')));
opm_cleaned = ft_selectdata(cfg, comb);

% Reject bad channels
cfg = [];
cfg.trl = trl_opm;
cfg.z_threshold = params.z_threshold;
cfg.corr_threshold = params.corr_threshold;
[~, ~, ~, ~, ~, ~, badtrl_opm_zmax] = opm_badchannels_restingstate(cfg, opm_raw);
cfg = [];
%cfg.channel = setdiff(opm_cleaned.label,badchs_opm); % No channels excluded since we already selected
cfg.trials  = setdiff(1:length(opm_cleaned.trial),badtrl_opm_zmax); % remove bad trials
opm_cleaned = ft_selectdata(cfg, opm_cleaned);

% HFC
i_chs = find(contains(opm_cleaned.label,'bz'));
chs = opm_cleaned.label(i_chs);
ori = zeros(length(chs),3);
for i = 1:length(chs)
    i_chs_grad(i) = find(strcmp(empty_room_cleaned.grad.label,chs{i}));
    ori(i,:) = empty_room_cleaned.grad.chanori(i_chs_grad(i),:);
end
opm_cleaned.grad.M = eye(size(ori,1)) - ori*pinv(ori);
for i = 1:length(opm_cleaned.trial)
    opm_cleaned.trial{i}(i_chs,:) = params.sign*opm_cleaned.grad.M*opm_cleaned.trial{i}(i_chs,:);
end
% update grad
opm_cleaned.grad.tra(i_chs_grad,i_chs_grad) = opm_cleaned.grad.M * opm_cleaned.grad.tra(i_chs_grad,i_chs_grad); 

% Reject jump trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,~] = ft_badsegment(cfg, opm_cleaned);
opm_cleaned = ft_rejectartifact(cfg,opm_cleaned);

% Reject noisy trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'std';
cfg.threshold = params.opm_std_threshold;
[cfg,~] = ft_badsegment(cfg, opm_cleaned);
opm_cleaned = ft_rejectartifact(cfg,opm_cleaned);

%% ICA
params.modality = 'opm';
params.layout = 'fieldlinebeta2bz_helmet.mat';
params.chs = '*bz';
opm_RS_ica = ica_MEG(opm_cleaned, save_path, params, 0);

cfg = [];
cfg.channel = '*bz';
opm_RS_ica = ft_selectdata(cfg,opm_RS_ica);

%% Save
save(fullfile(save_path, [params.sub '_resting_state_opm_ica']), 'opm_RS_ica',"-v7.3");

end
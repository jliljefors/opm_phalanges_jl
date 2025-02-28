function [squid_cleaned, squideeg_cleaned] = read_cvMEG(squid_file, save_path, params)
%prprocess_cvMEG Read conventional MEG data for benchmarking
% recordings. 
% Requires the following arguments:
% Path: containing save_path and squid_file
% Params: containing pre, post (pre- and poststim).

%% --- Read triggers ---
trl_squid = [];
cfg             = [];
cfg.datafile    = squid_file;
squid_raw         = ft_preprocessing(cfg);
squid_trig = find(contains(squid_raw.label,'STI101'));
trig = squid_raw.trial{1}(squid_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trl_squid(:,1) = find(trig)-(params.pre+params.pad)*squid_raw.fsample;
trl_squid(:,2) = find(trig)+(params.post+params.pad)*squid_raw.fsample;
trl_squid(:,3) = -(params.pre+params.pad)*squid_raw.fsample;
trl_squid(:,4) = squid_raw.trial{1}(squid_trig,trig);
trl_squid(:,1:2) = trl_squid(:,1:2) + floor(0.041*squid_raw.fsample); % adjust for stim delay
trl_squid = round(trl_squid);

%% MEG data filter & epoch
cfg = [];
%cfg.datafile        = squid_file;
%cfg.trl             = trl_meg;
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.hpinstabilityfix  = 'reduce';
%cfg.padding         = params.pre + params.post + 3;
%cfg.paddingtype     = 'data';
squid_epo = ft_preprocessing(cfg,squid_raw);

cfg = [];
cfg.trl             = trl_squid;
squid_epo = ft_redefinetrial(cfg,squid_epo);

cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
squid_epo = ft_preprocessing(cfg,squid_epo);

% Find bad channels
cfg = [];
cfg.trl = trl_squid;
cfg.trl(:,2) = cfg.trl(:,2) + squid_raw.fsample; % use 1sec longer trials for better neighborscorr
squid_raw = ft_redefinetrial(cfg,squid_raw);
cfg = [];
cfg.z_threshold = params.z_threshold;
cfg.corr_threshold = params.corr_threshold;
[badchs_opmeeg, badchs_squideeg_flat, badchs_squideeg_neighbors] = eeg_badchannels(cfg,squid_raw);
clear squid_raw

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

cfg = [];
cfg.channel = {'EOG', 'ECG'};
exg = ft_selectdata(cfg, squideeg_cleaned);

% Interpolate bad chs
cfg = [];
cfg.method = 'triangulation';
cfg.senstype = 'EEG';
neighbors = ft_prepare_neighbours(cfg,squideeg_cleaned);
cfg = [];
cfg.method = 'spline';
cfg.neighbors = neighbors;
cfg.badchannel = badchs_opmeeg;
cfg.senstype = 'EEG';
squideeg_cleaned = ft_channelrepair(cfg, squideeg_cleaned);

% Re-reference
% cfg = [];
% cfg.refef = 'yes';
% cfg.reffchannel = 'EEG023';
% squideeg_cleaned = ft_preprocessing(cfg,squideeg_cleaned);

cfg = [];
squideeg_cleaned = ft_appenddata(cfg,squideeg_cleaned,exg);

%cfg = [];
%cfg.channel = setdiff(squideeg_cleaned.label,badchs_opmeeg);
%squideeg_cleaned = ft_selectdata(cfg, squideeg_cleaned);

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

%% Spectra
cfg = [];
cfg.channel = 'EEG*';
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foilim = [1 100];
freq = ft_freqanalysis(cfg, squideeg_cleaned);
h = figure;
semilogy(freq.freq,freq.powspctrm)
xlabel('Frequency (Hz)')
ylabel('Power (T^2)')
title('SQUID-EEG spectrum - preICA')
saveas(h, fullfile(save_path, 'figs', [params.sub '_squideeg_spectrum.jpg']))

cfg = [];
cfg.channel = 'megmag';
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foilim = [1 100];
freq = ft_freqanalysis(cfg, squid_cleaned);
h = figure;
semilogy(freq.freq,freq.powspctrm)
xlabel('Frequency (Hz)')
ylabel('Power (T^2)')
title('SQUID-MAG spectrum - preICA')
saveas(h, fullfile(save_path, 'figs', [params.sub '_squidmag_spectrum1.jpg']))

cfg = [];
cfg.channel = 'meggrad';
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foilim = [1 100];
freq = ft_freqanalysis(cfg, squid_cleaned);
h = figure;
semilogy(freq.freq,freq.powspctrm)
xlabel('Frequency (Hz)')
ylabel('Power (T^2)')
title('SQUID-MAG spectrum - preICA')
saveas(h, fullfile(save_path, 'figs', [params.sub '_squidgrad_spectrum1.jpg']))


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
function [meg_ica] = ica_cvMEG(data,save_path,params)
%prprocess_osMEG Preprocessing on-scalp MEG data for benchmarking
% recordings. Requires arguments:
% Path: containing save_path and meg_file
% Params: containing pre, post (pre- and poststim), lp_freq, hp_freq,
% bp_freq and notch filter frequencies (corresponding filters are only
% applied if the frequency is defined), n_comp and coh_cutoff (for 
% automated ICA), and ds_freq (downsampling frequency).

%%
n_comp = params.n_comp;
ica_threshold = params.ica_threshold;

%% Downsample
cfg             = [];
cfg.resamplefs  = 200;
cfg.detrend     = 'no';
data_ds = ft_resampledata(cfg, data);

%% Calculate ICA
cfg            = [];
cfg.method     = 'runica';
cfg.numcomponent = n_comp;
cfg.channel    = 'meg';                     
comp           = ft_componentanalysis(cfg, data_ds);

%% Plot
cfg           = [];
cfg.component = 1:n_comp;       
cfg.layout    = 'neuromag306mag.lay'; 
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)

savefig(fullfile(save_path, 'meg_comp.fig')) 

%% --- ECG ----
% Find ECG artifacts
cfg                       = [];
cfg.continuous            = 'no';
cfg.artfctdef.ecg.pretim  = 0.25;
cfg.artfctdef.ecg.psttim  = 0.40-1/1200;
cfg.channel               = 'ECG';
cfg.artfctdef.ecg.inspect = 'ECG';
[cfg, ecg_artifact]       = ft_artifact_ecg(cfg,data);

%% Create ECG-locked
cfg = [];
cfg.dftfilter  = 'yes';
cfg.demean     = 'yes';
cfg.trl        = [ecg_artifact zeros(size(ecg_artifact,1), 2)];
temp = ft_redefinetrial(cfg, data);

% Separate MEG and ECG data
cfg.channel    = 'meg';
data_ecg = ft_selectdata(cfg, temp);
cfg.channel    = 'ECG';
ecg = ft_selectdata(cfg, temp);
ecg.channel{:} = 'ECG';

% Filter power line noise
cfg = [];
cfg.dftfilter       = 'yes';
cfg.dftfreq         = [50, 100, 150];
ecg = ft_preprocessing(cfg, ecg);

% Decompose the ECG-locked data
cfg = [];
cfg.unmixing  = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp_ecg = ft_componentanalysis(cfg, data_ecg);

% Combine ECG and ECG-locked data
comp_ecg = ft_appenddata([], ecg, comp_ecg);

% Plot averaged comp
cfg = [];
timelock = ft_timelockanalysis(cfg, comp_ecg);
figure
subplot(2,1,1); plot(timelock.time, timelock.avg(1,:)); title('ECG')
subplot(2,1,2); plot(timelock.time, timelock.avg(2:end,:));  title('ICA comp')

%% Correlation
ecg_comp_idx = [];
for i = 2:size(timelock.avg,1)
    tmp = corrcoef(timelock.avg(1,:), timelock.avg(i,:));
    R(i-1,1) = tmp(2,1);
end

idx = 1:n_comp;
ecg_comp_idx=idx(abs(R)>ica_threshold);
h = figure;
for i = 1:length(ecg_comp_idx)
    subplot(length(ecg_comp_idx),1,i)
    yyaxis left
    plot(timelock.time, timelock.avg(1,:));
    yyaxis right
    plot(timelock.time, timelock.avg(ecg_comp_idx(i)+1,:));  
    title(['Comp: ' num2str(ecg_comp_idx(i)) '; R_{ecg} = ' num2str(R(ecg_comp_idx(i),1))])
end
savefig(h,fullfile(save_path, 'figs', 'meg_ECG_comps.fig'))

%% --- EOG ---
% Find EOG artifacts
cfg = [];
cfg.continuous            = 'no';
cfg.channel               = 'EOG';
[~, eog_artifact] = ft_artifact_eog(cfg, data);

% Make artifact epochs
cfg = [];
cfg.dftfilter  = 'yes';
cfg.demean     = 'yes';
cfg.trl        = [eog_artifact zeros(size(eog_artifact,1), 1)];
temp = ft_redefinetrial(cfg, data);
    
% Separate MEG and EOG data
cfg.channel    = 'meg';
data_eog = ft_selectdata(cfg, temp);
cfg.channel    = 'EOG';
eog = ft_selectdata(cfg, temp);
eog.channel{:} = 'EOG';         % renaming for bookkeeping
    
% Filter power line noise
cfg = [];
cfg.dftfilter  = 'yes';
cfg.dftfreq    = [50, 100, 150];
eog = ft_preprocessing(cfg, eog);

% Decompose EOG-locked data
cfg = [];
cfg.unmixing  = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp_eog = ft_componentanalysis(cfg, data_eog);

% Combine EOG and EOG-locked data
comp_eog = ft_appenddata([], eog, comp_eog);

%% Correlation
cfg = [];
timelock = ft_timelockanalysis(cfg, comp_eog);
figure
subplot(2,1,1); plot(timelock.time, timelock.avg(1:2,:)); title('EOG')
subplot(2,1,2); plot(timelock.time, timelock.avg(3:end,:));  title('ICA comp')

eog1_comp_idx = [];
eog2_comp_idx = [];
for i = 3:size(timelock.avg,1)
    tmp = corrcoef(timelock.avg(1,:), timelock.avg(i,:));
    R(i-2,1) = tmp(2,1);
    tmp = corrcoef(timelock.avg(2,:), timelock.avg(i,:));
    R(i-2,2) = tmp(2,1);
end

eog1_comp_idx=idx(abs(R(:,1))>ica_threshold);
eog2_comp_idx=idx(abs(R(:,2))>ica_threshold);

h = figure;
for i = 1:length(eog1_comp_idx)
    subplot(length(eog1_comp_idx),1,i)
    yyaxis left
    plot(timelock.time, timelock.avg(1,:));
    yyaxis right
    plot(timelock.time, timelock.avg(eog1_comp_idx(i)+2,:));  
    title(['Comp: ' num2str(eog1_comp_idx(i)) '; R_{eog1} = ' num2str(R(eog1_comp_idx(i),1))])
end
savefig(h,fullfile(save_path, 'figs', 'meg_EOG1_comps.fig'))

h = figure;
for i = 1:length(eog2_comp_idx)
    subplot(length(eog2_comp_idx),1,i)
    yyaxis left
    plot(timelock.time, timelock.avg(2,:));
    yyaxis right
    plot(timelock.time, timelock.avg(eog2_comp_idx(i)+2,:));  
    title(['Comp: ' num2str(eog2_comp_idx(i)) '; R_{eog2} = ' num2str(R(eog2_comp_idx(i),2))])
end
savefig(h,fullfile(save_path, 'figs', 'meg_EOG2_comps.fig'))

%% Remove components
% Make a list of all "bad" components
reject_comp = unique([ecg_comp_idx eog1_comp_idx eog2_comp_idx]);

% Remove components
cfg = [];
cfg.component   = reject_comp;
cfg.channel     = 'meg';
cfg.updatesens  = 'no';
meg_ica = ft_rejectcomponent(cfg, comp, data);

save(fullfile(save_path, 'meg_ica_comp'), 'comp', 'ecg_comp_idx', 'eog1_comp_idx', 'eog2_comp_idx'); disp('done');

%% Save
save(fullfile(save_path, 'meg_ica_cleaned'), 'meg_ica',"-v7.3"); disp('done');

cfg = [];
cfg.latency = [-0.03 0.3];
meg_ica_ds = ft_selectdata(cfg, meg_ica);
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = [-0.03 0];
meg_ica_ds = ft_preprocessing(cfg,meg_ica_ds);
cfg = [];
cfg.resamplefs = 512;
meg_ica_ds = ft_resampledata(cfg, meg_ica_ds);
save(fullfile(save_path, 'meg_ica_ds.mat'), 'meg_ica_ds','-v7.3');   

end
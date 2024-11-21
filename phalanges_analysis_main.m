%% Set up
clear all
close all
restoredefaultpath
addpath('/home/chrpfe/Documents/MATLAB/fieldtrip-20231220/') % Fieldtrip path
addpath('/home/chrpfe/Documents/MATLAB/fieldtrip_private') % Fieldtrip private functions
addpath('/home/chrpfe/Documents/MATLAB/21099_opm/phalanges')
ft_defaults

global ft_default
ft_default.showcallinfo = 'no';

%% Params
overwrite = [];
overwrite.preproc = true;
overwrite.coreg = true;
overwrite.dip = true;
overwrite.mne = true;

params = [];
params.pre = 0.1; %sec
params.post = 0.4; %sec
params.filter = [];
params.filter.hp_freq = 3;
params.filter.lp_freq = 100;
params.filter.bp_freq = [];
params.filter.notch = sort([50:50:150 60:60:120]);
params.n_comp = 40;
params.ica_threshold = 0.8; % cutoff for EOG/ECG coherence
params.z_threshold = 20;
params.corr_threshold = 0.7; % correlation threshold for badchannel neighbors
params.opm_std_threshold = 2.5e-12;
params.eeg_std_threshold = 1e-4;
params.megmag_std_threshold = 5e-12;
params.megplanar_std_threshold = 5e-11;
params.hpi_freq = 33;

params.trigger_code = [2 4 8 16 32];
params.phalange_labels = {'I3' 'I2' 'I1' 'T1' 'I2b'};

%% Subjects + dates
subses = {'0005' '240208';
    '0905' '240229';
    '0916' '240320';
    '0953' '241104';
    '1096' '241022';
    '1153' '240321';
    '1167' '240425';
    '1186' '240925';
    '1190' '241023';
    '1191' '241024';
    '1193' '241029';
    '1194' '241029';
    '1195' '241030'};
mri_files = {'00000001.dcm' '/nifti/anat/sub-15931_T1w.nii.gz'  '/nifti/anat/sub-15985_T1w.nii.gz'};

%% Base paths
base_data_path = '/archive/21099_opm/';
base_save_path = '/home/chrpfe/Documents/21099_opm/Phalanges';

%% Loop over subjects
for i_sub = 1:size(subses,1)
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);

    %% Paths
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    if ~exist(save_path, 'dir')
       mkdir(save_path)
    end
    if ~exist(fullfile(save_path,'figs'), 'dir')
       mkdir(fullfile(save_path,'figs'))
    end
    meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
    if i_sub == 9
        meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98.fif');
    end
    opm_file = fullfile(raw_path, 'osmeg', 'PhalangesOPM_raw.fif');
    aux_file = fullfile(raw_path, 'meg', 'PhalangesEEG.fif');
    hpi_file = fullfile(raw_path, 'osmeg', 'HPIpre_raw.fif');

    %% --- OPM-MEG --------------------------------------------------------
    %% Preprocess
    % Read data, filter and reject bad channels/trials
    if exist(fullfile(save_path, [params.sub '_opm_cleaned.mat']),'file') && overwrite.preproc==false
        load(fullfile(save_path, [params.sub '_opm_cleaned.mat']));
        load(fullfile(save_path, [params.sub '_opmeeg_cleaned.mat']));
    else
        [opm_cleaned, opmeeg_cleaned] = read_osMEG(opm_file, aux_file, save_path, params); % Read data
    end

    %% Remove eye and heart artifacts with ICA
    if exist(fullfile(save_path, [params.sub '_opmeeg_ica.mat']),'file') && overwrite.preproc==false
        load(fullfile(save_path, [params.sub '_opm_ica.mat']));
        opm_ica = data_ica;
        load(fullfile(save_path, [params.sub '_opmeeg_ica.mat']));
        opmeeg_ica = data_ica;
        clear data_ica
        load(fullfile(save_path, [params.sub '_opmeeg_layout.mat']));
    else
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        opm_ica = ica_MEG(opm_cleaned, save_path, params);

        cfg = [];
        cfg.elec = opmeeg_cleaned.elec;
        cfg.output = fullfile(save_path, [params.sub '_opmeeg_layout.mat']);
        opmeeg_layout = ft_prepare_layout(cfg);
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        params.modality = 'opmeeg';
        opmeeg_ica = ica_MEG(opmeeg_cleaned, save_path, params);
        close all
    end
    %% Average
    if exist(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat']),'file') && overwrite.preproc==false
        fullfile(save_path, [params.sub '_opm_timelocked.mat'])
        fullfile(save_path, [params.sub '_opmeeg_timelocked.mat'])
    else
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        opm_timelocked = timelock_MEG(opm_ica, save_path, params);
        close all

        params.modality = 'opmeeg';
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        opmeeg_timelocked = timelock_MEG(opmeeg_ica, save_path, params);
        close all
    end

    %% --- SQUID-MEG ------------------------------------------------------
    ft_hastoolbox('mne', 1);
    %% Preprocess
    % Read data, filter and reject bad channels/trials
    if exist(fullfile(save_path, [params.sub '_megeeg_cleaned.mat']),'file') && overwrite.preproc==false
        load(fullfile(save_path, [params.sub '_meg_cleaned.mat']));
        load(fullfile(save_path, [params.sub '_megeeg_cleaned.mat']));
    else
        [meg_cleaned, megeeg_cleaned] = read_cvMEG(meg_file, save_path, params); % Read data
    end

    %% Remove eye and heart artifacts with ICA
    if exist(fullfile(save_path, [params.sub '_megeeg_ica.mat']),'file') && overwrite.preproc==false
        load(fullfile(save_path, [params.sub '_meg_ica.mat']));
        meg_ica = data_ica;
        load(fullfile(save_path, [params.sub '_megeeg_ica.mat']));
        megeeg_ica = data_ica;
        clear data_ica
        load(fullfile(save_path, [params.sub '_megeeg_layout.mat']));
    else
        params.modality = 'meg';
        params.layout = 'neuromag306all.lay';
        params.chs = 'MEG*';
        meg_ica = ica_MEG(meg_cleaned, save_path, params);

        cfg = [];
        cfg.elec = megeeg_cleaned.elec;
        cfg.output = fullfile(save_path, [params.sub '_megeeg_layout.mat']);
        megeeg_layout = ft_prepare_layout(cfg);
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        params.modality = 'megeeg';
        megeeg_ica = ica_MEG(megeeg_cleaned, save_path, params);
        close all
    end
    %% Average
    if exist(fullfile(save_path, [params.sub '_megeeg_timelocked.mat']),'file') && overwrite.preproc==false
        fullfile(save_path, [params.sub '_meg_timelocked.mat'])
        fullfile(save_path, [params.sub '_megeeg_timelocked.mat'])
    else
        params.modality = 'meg';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'megmag';
        opm_timelocked = timelock_MEG(meg_ica, save_path, params);
        close all

        params.modality = 'megeeg';
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        opmeeg_timelocked = timelock_MEG(megeeg_ica, save_path, params);
        close all
    end
end

%% --- Group sensor level -------------------------------------------------
%% Peak amplitude ratio
%peak_ratio_meg = zeros(size(subses,1),length(params.phalange_labels));
%peak_ratio_eeg = zeros(size(subses,1),length(params.phalange_labels));
for i_sub = 1:size(subses,1)
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);
    save_path = fullfile(base_save_path,params.sub);
    load(fullfile(save_path, [params.sub '_opm_M100'])); 
    M100_opm{i_sub} = M100;
    load(fullfile(save_path, [params.sub '_opmeeg_M100'])); 
    M100_opmeeg{i_sub} = M100;
    load(fullfile(save_path, [params.sub '_meg_M100'])); 
    M100_meg{i_sub} = M100;
    load(fullfile(save_path, [params.sub '_megeeg_M100']));
    M100_megeeg{i_sub} = M100;
    
    load(fullfile(save_path, [params.sub '_meg_timelocked'])); 
    meg_timelocked = timelocked;
    load(fullfile(save_path, [params.sub '_opm_timelocked'])); 
    meg_chs = find(contains(meg_timelocked{1}.label,'MEG'));
    opm_timelocked = timelocked;
    opm_chs = find(contains(opm_timelocked{1}.label,'bz'));

    for i_phalange = 1:length(params.phalange_labels)
        peak_ratio_meg(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.max_amplitude/M100_meg{i_sub}{i_phalange}.max_amplitude;
        peak_ratio_eeg(i_sub,i_phalange) = M100_opmeeg{i_sub}{i_phalange}.max_amplitude/M100_megeeg{i_sub}{i_phalange}.max_amplitude;
        peak_ratio_meg(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.max_amplitude/M100_meg{i_sub}{i_phalange}.max_amplitude;
        peak_ratio_eeg(i_sub,i_phalange) = M100_opmeeg{i_sub}{i_phalange}.max_amplitude/M100_megeeg{i_sub}{i_phalange}.max_amplitude;
        snr_error_opm(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.max_amplitude/M100_opm{i_sub}{i_phalange}.std_error;
        snr_error_meg(i_sub,i_phalange) = M100_meg{i_sub}{i_phalange}.max_amplitude/M100_meg{i_sub}{i_phalange}.std_error;
        snr_prestim_opm(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.max_amplitude/M100_opm{i_sub}{i_phalange}.prestim_std;
        snr_prestim_meg(i_sub,i_phalange) = M100_meg{i_sub}{i_phalange}.max_amplitude/M100_meg{i_sub}{i_phalange}.prestim_std;
        snr_ratio_error(i_sub,i_phalange) = snr_error_opm(i_sub,i_phalange)/snr_error_meg(i_sub,i_phalange);
        snr_ratio_prestim(i_sub,i_phalange) = snr_prestim_opm(i_sub,i_phalange)/snr_prestim_meg(i_sub,i_phalange);
        latency_opm(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.peak_latency;
        latency_meg(i_sub,i_phalange) = M100_meg{i_sub}{i_phalange}.peak_latency;
        latency_opmeeg(i_sub,i_phalange) = M100_opmeeg{i_sub}{i_phalange}.peak_latency;
        latency_megeeg(i_sub,i_phalange) = M100_megeeg{i_sub}{i_phalange}.peak_latency;

        h = figure;
        subplot(2,1,1)
        plot(meg_timelocked{i_phalange}.time*1e3,meg_timelocked{i_phalange}.avg(meg_chs(1:3:end),:)*1e15)
        xlabel('t [msec]')
        ylabel('B [fT]')
        title(['Evoked SQUID MAG - phalange ' params.phalange_labels{i_phalange} ' (n_trls=' num2str(length(meg_timelocked{i_phalange}.cfg.trials)) ')'])
        subplot(2,1,2)
        plot(opm_timelocked{i_phalange}.time*1e3,meg_timelocked{i_phalange}.avg(opm_chs,:)*1e15)
        xlabel('t [msec]')
        ylabel('B [fT]')
        title(['Evoked OPM MAG - phalange ' params.phalange_labels{i_phalange} ' (n_trls=' num2str(length(opm_timelocked{i_phalange}.cfg.trials)) ')'])
        saveas(h, fullfile(save_path, 'figs', [params.sub '_megopm_butterfly_ph-' params.phalange_labels{i_phalange} '.jpg']))
    end
    close all
end
%% Plot group figures
save_path = base_save_path;
if ~exist(fullfile(save_path,'figs'), 'dir')
       mkdir(fullfile(save_path,'figs'))
end

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(peak_ratio_meg,1));
hold on
er = errorbar(1:5,mean(peak_ratio_meg,1), mean(peak_ratio_meg,1)-min(peak_ratio_meg,[],1), mean(peak_ratio_meg,1)-max(peak_ratio_meg,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 peak amplitude ratio (mean = ' num2str(mean(mean(peak_ratio_meg))) ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'Peak_amplitude_ratios_meg.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(peak_ratio_eeg,1));
hold on
er = errorbar(1:5,mean(peak_ratio_eeg,1), mean(peak_ratio_eeg,1)-min(peak_ratio_eeg,[],1), mean(peak_ratio_eeg,1)-max(peak_ratio_eeg,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 peak amplitude ratio (mean = ' num2str(mean(mean(peak_ratio_eeg))) ')'])
ylabel('OPMEEG/SQUIDEEG')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'Peak_amplitude_ratios_eeg.jpg'))

% SNR
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(snr_ratio_error,1));
hold on
er = errorbar(1:5,mean(snr_ratio_error,1), mean(snr_ratio_error,1)-min(snr_ratio_error,[],1), mean(snr_ratio_error,1)-max(snr_ratio_error,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 SNR_{stderror} ratio (mean = ' num2str(mean(mean(snr_ratio_error))) ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'SNR_ratios_error.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(snr_ratio_prestim,1));
hold on
er = errorbar(1:5,mean(snr_ratio_prestim,1), mean(snr_ratio_prestim,1)-min(snr_ratio_prestim,[],1), mean(snr_ratio_prestim,1)-max(snr_ratio_prestim,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 SNR_{prestim} ratio (mean = ' num2str(mean(mean(snr_ratio_prestim))) ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'SNR_ratios_prestim.jpg'))

%% Peak latency
tmp = latency_opm-latency_meg;
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(tmp,1));
hold on
er = errorbar(1:5,mean(tmp,1), mean(tmp,1)-min(tmp,[],1), mean(tmp,1)-max(tmp,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 latency diff (opm mean = ' num2str(mean(mean(latency_opm))) '; meg mean = ' num2str(mean(mean(latency_meg))) ')'])
ylabel('t_{OPM}-t_{SQUID}')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'Latency.jpg'))


%% HPI localization
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    aux_file = fullfile(raw_path, 'meg', 'PhalangesEEG.fif');
    hpi_file = fullfile(raw_path, 'osmeg', 'HPIpre_raw.fif');

    if exist(fullfile(save_path, 'opm_trans.mat'),'file') && overwrite.coreg==false
        load(fullfile(save_path, 'hpi_fit.mat'));
        load(fullfile(save_path, 'opm_trans.mat'));
    else
        ft_hastoolbox('mne', 1);
        load(fullfile(save_path, [params.sub '_opm_ica_ds']));
        params.include_chs = data_ica_ds.label(find(contains(data_ica_ds.label,'bz')));
        params.include_chs = params.include_chs([1:75 77:end]);
        [hpi_fit, opm_trans, hpi_fit_tf] = fit_hpi(hpi_file, aux_file, save_path, params);
        close all
    end
end

%% Prepare MRIs
for i_sub = 1:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    if exist(fullfile(save_path, 'headmodels.mat'),'file') && overwrite.coreg==false
        load(fullfile(save_path, 'headmodels.mat'));
        load(fullfile(save_path, 'meshes.mat'));
        load(fullfile(save_path, 'mri_resliced_cm.mat'));
    else
        meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_tsss.fif');
        mri_file = fullfile(mri_path, mri_files{i_sub});
        [headmodels, meshes, mri_resliced_cm] = prepare_mri(mri_file,meg_file,aux_file,save_path);
        close all
    end
    
    %% Transform for OPM
    opm_timelockedT = opm_timelocked;
    for i = 1:5
        opm_timelockedT{i}.grad.chanpos = opm_trans.transformPointsForward(opm_timelocked{i}.grad.chanpos*1e2)*1e-2;
        opm_timelockedT{i}.grad.coilpos = opm_trans.transformPointsForward(opm_timelocked{i}.grad.coilpos*1e2)*1e-2;
        opm_timelockedT{i}.elec.chanpos = meg_timelocked{i}.elec.chanpos;
        opm_timelockedT{i}.elec.elecpos = meg_timelocked{i}.elec.elecpos;
    end
    
    h = figure; 
    hold on; 
    ft_plot_sens(opm_timelockedT{1}.grad,'unit','cm')
    ft_plot_sens(opm_timelockedT{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
    ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.7,'FaceColor',[229 194 152]/256,'unit','cm')
    ft_plot_headmodel(headmodels.headmodel_meg)
    hold off;
    title('OPM-MEG')
    view([-140 10])
    savefig(h, fullfile(save_path, 'figs', 'opm_layout.fig'))
    saveas(h, fullfile(save_path, 'figs', 'opm_layout.jpg'))

    h = figure; 
    hold on
    ft_plot_sens(meg_timelocked{1}.grad,'unit','cm')
    ft_plot_sens(meg_timelocked{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
    ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.7,'FaceColor',[229 194 152]/256,'unit','cm')
    ft_plot_headmodel(headmodels.headmodel_meg)
    hold off;
    title('SQUID-MEG')
    view([-140 10])
    savefig(h, fullfile(save_path, 'figs', 'meg_layout.fig'))
    saveas(h, fullfile(save_path, 'figs', 'meg_layout.jpg'))

    %% Dipole fits
    if exist(fullfile(save_path, 'dipoles.mat'),'file') && overwrite.dip==false
        load(fullfile(save_path, 'dipoles.mat'));
    else
        [megmag_dipole, megplanar_dipole, opm_dipole, eeg_dipole] = fit_dipoles(save_path,meg_timelocked,opm_timelockedT,headmodels,mri_resliced_cm,params);
    end

    %% Read cotrical restrained source model
    subjectname = ['NatMEG_' sub{i_sub}];
    filename = fullfile(mri_path,'workbench',[subjectname,'.L.midthickness.8k_fs_LR.surf.gii']);
    sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});

    T = mri_resliced_cm.transform/mri_resliced_cm.hdr.vox2ras;
    sourcemodelT = ft_transform_geometry(T, sourcemodel);
    sourcemodelT.inside = sourcemodelT.atlasroi>0;
    sourcemodelT = rmfield(sourcemodelT, 'atlasroi');

    %sourcemodelOPM = sourcemodelT;
    %sourcemodelOPM.pos = opm_trans.transformPointsForward(sourcemodelOPM.pos);

    filename =fullfile(mri_path,'workbench',[subjectname,'.L.midthickness.32k_fs_LR.surf.gii']);
    cortex = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});
    cortex = ft_transform_geometry(T, cortex);
    cortex = rmfield(cortex, 'sulc');
    cortex = rmfield(cortex, 'curv');
    cortex = rmfield(cortex, 'thickness');
    cortex = rmfield(cortex, 'brainstructure');
    cortex = rmfield(cortex, 'atlasroi');

    % Plot source and head models
    h=figure; 
    ft_plot_mesh(sourcemodelT, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
    hold on; 
    ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.2,'FaceColor',[229 194 152]/256,'unit','cm')
    ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
    ft_plot_sens(opm_timelockedT{1}.grad,'unit','cm')
    ft_plot_sens(opm_timelockedT{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
    hold off;
    title('OPM-MEG')
    view([-140 10])
    savefig(h, fullfile(save_path, 'figs', 'opm_layout2.fig'))
    saveas(h, fullfile(save_path, 'figs', 'opm_layout2.jpg'))

    h=figure; 
    ft_plot_mesh(sourcemodelT, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
    hold on; 
    ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.2,'FaceColor',[229 194 152]/256,'unit','cm')
    ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
    ft_plot_sens(meg_timelocked{1}.grad,'unit','cm')
    ft_plot_sens(meg_timelocked{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
    hold off;
    title('SQUID-MEG')
    view([-140 10])
    savefig(h, fullfile(save_path, 'figs', 'meg_layout2.fig'))
    saveas(h, fullfile(save_path, 'figs', 'meg_layout2.jpg'))

    %% MNE fit
    if exist(fullfile(save_path, 'mne_fits.mat'),'file') && overwrite.mne==false
        load(fullfile(save_path, 'mne_fits.mat'));
    else
        [megmag_mne, megplanaer_mne, opm_mne, eeg_mne, FAHM] = fit_mne(save_path,meg_timelocked,opm_timelockedT,headmodels,sourcemodelT,params);
    end
end
close all

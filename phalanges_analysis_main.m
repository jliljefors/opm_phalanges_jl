%% Reset all
clear all
close all
restoredefaultpath

%% Base paths
if contains(pwd,'/home/chrpfe')
    % Server:
    base_data_path = '/archive/21099_opm/';
    base_save_path = '/home/chrpfe/Documents/21099_opm/Phalanges';
    base_matlab_path = '/home/chrpfe/Documents/MATLAB/';
    project_scripts_path = '/home/chrpfe/Documents/MATLAB/21099_opm/phalanges';
else
    % Laptop:
    base_data_path = '/Users/christophpfeiffer/data_archive/21099_opm';
    base_save_path = '/Users/christophpfeiffer/data_local/Benchmarking_phalanges';
    base_matlab_path = '/Users/christophpfeiffer/Dropbox/Mac/Documents/MATLAB';
    project_scripts_path = '/Users/christophpfeiffer/opm_phalanges';
end

%% Set up fieldtrip
addpath(fullfile(base_matlab_path,'fieldtrip-20231220/')) % Fieldtrip path
addpath(fullfile(base_matlab_path,'fieldtrip_private')) % Fieldtrip private functions
addpath(project_scripts_path)
ft_defaults

global ft_default
ft_default.showcallinfo = 'no';

%% Params
overwrite = [];
overwrite.preproc = true;
overwrite.coreg = true;
overwrite.mri = false;
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
mri_files = {'00000001.dcm' 
    '/mri/sub-15931_T1w.nii.gz'  
    '/nifti/anat/sub-15985_T1w.nii.gz'};

%% Loop over subjects
for i_sub = 1:size(subses,1)
    params.sub = ['sub_' num2str(i_sub,'%02d')];

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

    %% OPM-MEG 
    if exist(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat']),'file') && overwrite.preproc==false
        load(fullfile(save_path, [params.sub '_opm_timelocked.mat']))
        opm_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat']))
        opmeeg_timelocked = timelocked;
        clear timelocked
    else
        ft_hastoolbox('mne', 1);

        % Read data
        [opm_cleaned, opmeeg_cleaned] = read_osMEG(opm_file, aux_file, save_path, params); % Read data
        
        % ICA
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

        % Average
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

    %% SQUID-MEG 
    if exist(fullfile(save_path, [params.sub '_megeeg_timelocked.mat']),'file') && overwrite.preproc==false
        load(fullfile(save_path, [params.sub '_megmag_timelocked.mat']))
        megmag_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_megplanar_timelocked.mat']))
        megplanar_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_megeeg_timelocked.mat']))
        megeeg_timelocked = timelocked;
        clear timelocked
    else
        ft_hastoolbox('mne', 1);

        % Read data
        [meg_cleaned, megeeg_cleaned] = read_cvMEG(meg_file, save_path, params); % Read data
        
        % ICA
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

        % Average
        params.modality = 'megmag';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'megmag';
        meg_timelocked = timelock_MEG(meg_ica, save_path, params);
        close all

        params.modality = 'megplanar';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'megplanarr';
        megplanar_timelocked = timelock_MEG(meg_ica, save_path, params);
        close all

        params.modality = 'megeeg';
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        meegeeg_timelocked = timelock_MEG(megeeg_ica, save_path, params);
        close all
    end
end

%% --- Group sensor level -------------------------------------------------
if ~exist(fullfile(base_save_path,'figs'), 'dir')
       mkdir(fullfile(base_save_path,'figs'))
end
subs = 1:13;
sensor_results_goup(base_save_path,subs, params)

%% Prepare MRIs
for i_sub = 6%7:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}],'mri');
    if exist(fullfile(save_path, 'headmodels.mat'),'file') && overwrite.mri==false
        load(fullfile(save_path, 'headmodels.mat'));
        load(fullfile(save_path, 'meshes.mat'));
        load(fullfile(save_path, 'mri_resliced.mat'));
    else
        meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98.fif');
        end
        mri_file = fullfile(mri_path, 'orig','001.mgz');
        [headmodels, meshes] = prepare_mri(mri_file,meg_file,save_path);
        close all
    end
end

%% HPI localization
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    hpi_file = fullfile(raw_path, 'osmeg', 'HPIpre_raw.fif');

    if exist(fullfile(save_path, 'opm_trans.mat'),'file') && overwrite.coreg==false
        load(fullfile(save_path, 'hpi_fit.mat'));
        load(fullfile(save_path, 'opm_trans.mat'));
    else
        meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98.fif');
        end
        ft_hastoolbox('mne', 1);
        load(fullfile(save_path, [params.sub '_opm_ica_ds']));
        params.include_chs = data_ica_ds.label(find(contains(data_ica_ds.label,'bz')));
        params.include_chs = params.include_chs([1:75 77:end]);
        [hpi_fit, opm_trans, hpi_fit_tf] = fit_hpi(hpi_file, meg_file, save_path, params);
        close all
    end
end

% Transform for OPM data
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    if exist(fullfile(save_path, 'headmodels.mat'),'file') && exist(fullfile(save_path, 'opm_trans.mat'),'file') && overwrite.coreg==true
        load(fullfile(save_path, 'opm_trans.mat'));
        load(fullfile(save_path, 'headmodels.mat'));
        load(fullfile(save_path, 'meshes.mat'));
        load(fullfile(save_path, [params.sub '_opm_timelocked.mat']))
        opm_timelocked = timelocked;
        opm_timelockedT = opm_timelocked;
        load(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat']))
        opmeeg_timelocked = timelocked;
        opmeeg_timelockedT = opmeeg_timelocked;
        load(fullfile(save_path, [params.sub '_meg_timelocked.mat']))
        meg_timelocked = timelocked;
        clear timelocked;
        meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98.fif');
        end
        headshape = ft_read_headshape(meg_file);
        for i = 1:5
            opm_timelockedT{i}.grad.chanpos = opm_trans.transformPointsForward(opm_timelocked{i}.grad.chanpos*1e2)*1e-2;
            opm_timelockedT{i}.grad.coilpos = opm_trans.transformPointsForward(opm_timelocked{i}.grad.coilpos*1e2)*1e-2;
            opmeeg_timelockedT{i}.elec.chanpos = meg_timelocked{i}.elec.chanpos;
            opmeeg_timelockedT{i}.elec.elecpos = meg_timelocked{i}.elec.elecpos;
        end
        
        h = figure; 
        hold on; 
        ft_plot_sens(opm_timelockedT{1}.grad,'unit','cm')
        ft_plot_sens(opmeeg_timelockedT{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
        ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.7,'FaceColor',[229 194 152]/256,'unit','cm')
        ft_plot_headmodel(headmodels.headmodel_meg)
        ft_plot_headshape(headshape)
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
        ft_plot_headshape(headshape)
        hold off;
        title('SQUID-MEG')
        view([-140 10])
        savefig(h, fullfile(save_path, 'figs', 'meg_layout.fig'))
        saveas(h, fullfile(save_path, 'figs', 'meg_layout.jpg'))
        close all

        %% Save
        save(fullfile(save_path, [params.sub '_opm_timelockedT']), 'opm_timelockedT', '-v7.3');
        save(fullfile(save_path, [params.sub '_opmeeg_timelockedT']), 'opm_timelockedT', '-v7.3');

    else
        disp('Required files not found. No transformed OPM data was saved.')
    end
end

%% Dipole fits
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);

    if exist(fullfile(save_path, 'dipoles.mat'),'file') && overwrite.dip==false
        load(fullfile(save_path, 'dipoles.mat'));
    else
        clear headmodels mri_resliced
        clear megmag_timelocked magplanar_timelocked megeeg_timelocked
        clear opm_timelockedT opmeeg_timelockedT
        load(fullfile(save_path, [params.sub '_opm_timelockedT.mat']))
        load(fullfile(save_path, [params.sub '_opmeeg_timelockedT.mat']))
        load(fullfile(save_path, [params.sub '_megmag_timelocked.mat']))
        megmag_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_megplanar_timelocked.mat']))
        megplanar_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_megeeg_timelocked.mat']))
        megeeg_timelocked = timelocked;
        clear timelocked
        m100_latency = cell(5,1);
        for i_ph = 1:5
            m100_latency{i_ph} = [];
            load(fullfile(save_path, [params.sub '_opm_M100'])); 
            m100_latency{i_ph}.opm = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_opmeeg_M100'])); 
            m100_latency{i_ph}.opmeeg = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_megmag_M100'])); 
            m100_latency{i_ph}.megmag = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_megplanar_M100'])); 
            m100_latency{i_ph}.megplanar = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_megeeg_M100']));
            m100_latency{i_ph}.megeeg = M100{i_ph}.peak_latency;
        end
        load(fullfile(save_path, 'headmodels.mat'));
        load(fullfile(save_path, 'mri_resliced.mat'));
        [megmag_dipole, megplanar_dipole, opm_dipole, eeg_dipole] = fit_dipoles(save_path,megmag_timelocked,megplanar_timelocked,megeeg_timelocked,opm_timelockedT,opmeg_timelockedT,headmodels,mri_resliced,m100_latency,params);
    end
end

%% --- Group source level -------------------------------------------------
if ~exist(fullfile(base_save_path,'figs'), 'dir')
       mkdir(fullfile(base_save_path,'figs'))
end
subs = 2:13;
dipole_results_goup(base_save_path,subs, params)

%% Prepare MNE sourcemodel 
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);

    clear headmodels megmag_timelocked opm_timelockedT meshes filename
    clear headmodels meshes filename
    load(fullfile(save_path,'headmodels.mat'));
    load(fullfile(save_path,'meshes.mat'));    
    load(fullfile(save_path, 'mri_resliced.mat'));
    load(fullfile(save_path, [params.sub '_opm_timelockedT.mat']))
    load(fullfile(save_path, [params.sub '_opmeeg_timelockedT.mat']))
    load(fullfile(save_path, [params.sub '_megmag_timelocked.mat']))
    megmag_timelocked = timelocked;
    load(fullfile(save_path, [params.sub '_megeeg_timelocked.mat']))
    megeeg_timelocked = timelocked;
    clear timelocked

    % Read and transform cortical restrained source model
    files = dir(fullfile(mri_path,'workbench'));
    for i = 1:length(files)
        if endsWith(files(i).name,'.L.midthickness.8k_fs_LR.surf.gii')
            filename = fullfile(mri_path,'workbench',files(i).name);
            %return;
        end
    end
    sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});

    T = mri_resliced.transform/mri_resliced.hdr.vox2ras;
    sourcemodelT = ft_transform_geometry(T, sourcemodel);
    sourcemodelT.inside = true(size(sourcemodelT.pos,1),1);
    %sourcemodelOPM = sourcemodelT;
    %sourcemodelOPM.pos = opm_trans.transformPointsForward(sourcemodelOPM.pos);

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
    ft_plot_sens(megmag_timelocked{1}.grad,'unit','cm')
    ft_plot_sens(megmag_timelocked{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
    hold off;
    title('SQUID-MEG')
    view([-140 10])
    savefig(h, fullfile(save_path, 'figs', 'meg_layout2.fig'))
    saveas(h, fullfile(save_path, 'figs', 'meg_layout2.jpg'))
    
    close all
    % Save
    save(fullfile(save_path, [params.sub '_sourcemodelT']), 'sourcemodelT', '-v7.3');
end

%% MNE
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);

    %% MNE fit
    if exist(fullfile(save_path, 'mne_fits.mat'),'file') && overwrite.mne==false
        load(fullfile(save_path, 'mne_fits.mat'));
    else
        clear headmodels sourcemodelT
        clear megmag_timelocked magplanar_timelocked megeeg_timelocked
        clear opm_timelockedT opmeeg_timelockedT
        load(fullfile(save_path, [params.sub '_opm_timelockedT.mat']))
        load(fullfile(save_path, [params.sub '_opmeeg_timelockedT.mat']))
        load(fullfile(save_path, [params.sub '_megmag_timelocked.mat']))
        megmag_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_megplanar_timelocked.mat']))
        megplanar_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_megeeg_timelocked.mat']))
        megeeg_timelocked = timelocked;
        clear timelocked
        m100_latency = cell(5,1);
        for i_ph = 1:5
            m100_latency{i_ph} = [];
            load(fullfile(save_path, [params.sub '_opm_M100'])); 
            m100_latency{i_ph}.opm = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_opmeeg_M100'])); 
            m100_latency{i_ph}.opmeeg = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_megmag_M100'])); 
            m100_latency{i_ph}.megmag = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_megplanar_M100'])); 
            m100_latency{i_ph}.megplanar = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_megeeg_M100']));
            m100_latency{i_ph}.megeeg = M100{i_ph}.peak_latency;
        end
        [megmag_mne, megplanaer_mne, opm_mne, eeg_mne, FAHM] = fit_mne(save_path,megmag_timelocked,megplanar_timelocked,opm_timelockedT,headmodels,sourcemodelT,m100_latency,params);
    end
end
close all

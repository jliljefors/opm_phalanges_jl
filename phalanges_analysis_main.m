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

%% Overwrite
overwrite = [];
overwrite.preproc = false;
overwrite.coreg = false;
overwrite.mri = false;
overwrite.dip = false;
overwrite.empty_room = true;
overwrite.mne = true;

%% Params
params = [];
params.pre = 0.05; %sec
params.post = 0.3; %sec
params.pad = 0.2; %sec
params.filter = [];
params.filter.hp_freq = 1;
params.filter.lp_freq = 70;
params.filter.bp_freq = [];
params.filter.notch = [50 60 80]; %[50 60 100 120 150];
params.n_comp = 40;
params.ica_cor = 0.8; % cutoff for EOG/ECG coherence
params.ica_coh = 0.95; % cutoff for EOG/ECG coherence
params.z_threshold = 20;
params.corr_threshold = 0.6; % correlation threshold for badchannel neighbors
params.opm_std_threshold = 2.5e-12;
params.eeg_std_threshold = 1e-4;
params.squidmag_std_threshold = 5e-12;
params.squidgrad_std_threshold = 5e-11;
params.hpi_freq = 33;
params.hpi_gof = 0.9;
params.M60 = [-0.02 0.02]+0.06;

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
    if exist(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat']),'file') && exist(fullfile(save_path, [params.sub '_squideeg_timelocked.mat']),'file') && overwrite.preproc==false
        disp(['Not overwriting preproc for ' params.sub]);
    else
        ft_hastoolbox('mne', 1);

        if i_sub <=3 % Flip amplitudes in old recordings
            params.sign  = -1;
        else
            params.sign = 1;
        end
        % Read data
        [opm_cleaned, opmeeg_cleaned] = read_osMEG(opm_file, aux_file, save_path, params); % Read data
        opm_cleaned.grad = ft_convert_units(opm_cleaned.grad,'cm');

        % OPM ICA
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        data_ica = ica_MEG(opm_cleaned, save_path, params, 1);
        clear opm_cleaned

        % OPM Average
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        timelock_MEG(data_ica, save_path, params);
        close all
        clear data_ica

        % OPM-EEG ICA 
        cfg = [];
        cfg.elec = opmeeg_cleaned.elec;
        cfg.output = fullfile(save_path, [params.sub '_opmeeg_layout.mat']);
        opmeeg_layout = ft_prepare_layout(cfg);
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        params.modality = 'opmeeg';
        data_ica = ica_MEG(opmeeg_cleaned, save_path, params, 1);
        close all
        clear opmeeg_cleaned

        params.modality = 'opmeeg';
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        timelock_MEG(data_ica, save_path, params);
        close all
        clear data_ica

        params = rmfield(params,{'modality', 'layout', 'chs', 'amp_scaler', 'amp_label'}); % remove fields used for picking modality
        
        %% SQUID-MEG 
        ft_hastoolbox('mne', 1);

        % Read data
        [squid_cleaned, squideeg_cleaned] = read_cvMEG(meg_file, save_path, params); % Read data
        
        % SQUID-MAG ICA
        params.modality = 'squid';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'meg';
        data_ica = ica_MEG(squid_cleaned, save_path, params, 1);      

        % SQUID-MAG timelock
        params.modality = 'squidmag';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'megmag';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        timelock_MEG(data_ica, save_path, params);
        close all

        % SQUID-GRAD timelock
        params.modality = 'squidgrad';
        params.layout = 'neuromag306planar.lay';
        params.chs = 'megplanar';
        params.amp_scaler = 1e15/100;
        params.amp_label = 'B [fT/cm]';
        timelock_MEG(data_ica, save_path, params);
        close all
        clear data_ica

        % EEG ICA
        cfg = [];
        cfg.elec = squideeg_cleaned.elec;
        cfg.output = fullfile(save_path, [params.sub '_megeeg_layout.mat']);
        megeeg_layout = ft_prepare_layout(cfg);
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        params.modality = 'squideeg';
        data_ica = ica_MEG(squideeg_cleaned, save_path, params, 1);
        close all
        clear squideeg_cleaned

        % EEG timelock
        params.modality = 'squideeg';
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        timelock_MEG(data_ica, save_path, params);
        close all
        clear data_íca

        params = rmfield(params,{'modality', 'layout', 'chs', 'amp_scaler', 'amp_label'}); % remove fields used for picking modality    
    
        %% Save results in report
        create_bads_reports(base_save_path, i_sub, params);
    end
end

%% Sensor level group analysis
if ~exist(fullfile(base_save_path,'figs'), 'dir')
       mkdir(fullfile(base_save_path,'figs'))
end
subs = 1:13;
sensor_results_goup(base_save_path,subs, params)
close all

%% Prepare MRIs
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}],'mri');
    if exist(fullfile(save_path, 'headmodels.mat'),'file') && overwrite.mri==false
        disp(['Not overwriting MRI for ' params.sub]);
    else
        meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98.fif');
        end
        mri_file = fullfile(mri_path, 'orig','001.mgz');
        prepare_mri(mri_file,meg_file,save_path);
        close all
    end
end

%% HPI localization
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    hpi_path = fullfile(raw_path, 'osmeg'); %hpi_file = fullfile(raw_path, 'osmeg', 'HPIpre_raw.fif');

    if exist(fullfile(save_path, 'opm_trans.mat'),'file') && overwrite.coreg==false
        disp(['Not overwriting OPM transform for ' params.sub]);
    else
        meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98.fif');
        end
        ft_hastoolbox('mne', 1);
        load(fullfile(save_path, [params.sub '_opm_ica_ds']));
        params.include_chs = data_ica_ds.label(find(contains(data_ica_ds.label,'bz')));
        fit_hpi(hpi_path, meg_file, save_path, params);
        close all
    end
end

%% Transform for OPM data
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    if exist(fullfile(save_path, 'headmodels.mat'),'file') && exist(fullfile(save_path, 'opm_trans.mat'),'file') && or(overwrite.preproc,or(overwrite.coreg,overwrite.mri))
        clear headmodels meshes filename headshape 
        headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;
        meshes = load(fullfile(save_path,'meshes.mat')).meshes;    
        mri_resliced = load(fullfile(save_path, 'mri_resliced.mat')).mri_resliced;
        opm_trans = load(fullfile(save_path, 'opm_trans.mat')).opm_trans;
        
        clear opm_timelockedT opmeeg_timelcokedT 
        clear squideeg_timelocked squidmag_timelocked
        opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelocked.mat'])).timelocked;
        opmeeg_timelockedT = load(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat'])).timelocked;
        squideeg_timelocked = load(fullfile(save_path, [params.sub '_squideeg_timelocked.mat'])).timelocked;
        squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat'])).timelocked;

        % Transform opm & opmeeg data 
        meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98.fif');
        end
        headshape = ft_read_headshape(meg_file);

        for i = 1:5
            opm_timelockedT{i}.grad.chanpos = opm_trans.transformPointsForward(opm_timelockedT{i}.grad.chanpos);
            opm_timelockedT{i}.grad.coilpos = opm_trans.transformPointsForward(opm_timelockedT{i}.grad.coilpos);
            opm_timelockedT{i}.grad.chanori = (opm_trans.Rotation'*opm_timelockedT{i}.grad.chanori')';
            opm_timelockedT{i}.grad.coilori = (opm_trans.Rotation'*opm_timelockedT{i}.grad.coilori')';
            opmeeg_timelockedT{i}.elec.chanpos = squideeg_timelocked{i}.elec.chanpos;
            opmeeg_timelockedT{i}.elec.elecpos = squideeg_timelocked{i}.elec.elecpos;
        end

        % Read and transform cortical restrained source model
        clear sourcemodel sourcemodel_inflated
        files = dir(fullfile(mri_path,'workbench'));
        if i_sub ==5
            files = dir(fullfile(save_path,'wb'));
        end
        for i = 1:length(files)
            if endsWith(files(i).name,'.L.midthickness.8k_fs_LR.surf.gii')
                filename = fullfile(mri_path,'workbench',files(i).name);
            elseif endsWith(files(i).name,'.L.aparc.8k_fs_LR.label.gii')
                filename2 = fullfile(mri_path,'workbench',files(i).name);
            end
        end
        sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});

        aparc_L = ft_read_atlas({filename2,filename});
        aparc_R = ft_read_atlas({strrep(filename2,'.L.','.R.'),strrep(filename,'.L.','.R.')});
        tmp = ft_read_atlas(strrep(filename2, '.L.', '.R.'),'format','caret_label');
        n_labels = length(aparc_L.parcellationlabel);
        atlas = [];
        atlas.parcellationlabel = [aparc_L.parcellationlabel; aparc_R.parcellationlabel];
        atlas.parcellation = [aparc_L.parcellation; aparc_R.parcellation + n_labels];
        atlas.rgba = [aparc_L.rgba; aparc_R.rgba; [0 0 0 1]];
        n_labels = length(atlas.parcellationlabel);
        atlas.parcellation(isnan(atlas.parcellation))=n_labels+1;
        sourcemodel.brainstructure = atlas.parcellation;
        sourcemodel.brainstructurelabel = atlas.parcellationlabel;
        sourcemodel.brainstructurecolor = atlas.rgba;

        T = mri_resliced.transform/mri_resliced.hdr.vox2ras;
        sourcemodel = ft_transform_geometry(T, sourcemodel);
        sourcemodel.inside = true(size(sourcemodel.pos,1),1);

        for i = 1:length(files)
            if endsWith(files(i).name,'.L.inflated.8k_fs_LR.surf.gii')
                filename = fullfile(mri_path,'workbench',files(i).name);
            end
        end
        sourcemodel_inflated = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});
        sourcemodel_inflated = ft_transform_geometry(T, sourcemodel_inflated);
        sourcemodel_inflated.inside = true(size(sourcemodel_inflated.pos,1),1);
        sourcemodel_inflated.brainstructure = sourcemodel.brainstructure;
        sourcemodel_inflated.brainstructurelabel = sourcemodel.brainstructurelabel;
        sourcemodel_inflated.brainstructurecolor = sourcemodel.brainstructurecolor;

        % Plot source and head models
        h=figure; 
        ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
        hold on; 
        ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.2,'FaceColor',[229 194 152]/256,'unit','cm')
        ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
        ft_plot_sens(opm_timelockedT{1}.grad,'unit','cm')
        ft_plot_sens(opmeeg_timelockedT{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
        hold off;
        title('OPM-MEG')
        view([-140 10])
        savefig(h, fullfile(save_path, 'figs', 'opm_layout2.fig'))
        saveas(h, fullfile(save_path, 'figs', 'opm_layout2.jpg'))
    
        h=figure; 
        ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
        hold on; 
        ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.2,'FaceColor',[229 194 152]/256,'unit','cm')
        ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
        ft_plot_sens(squidmag_timelocked{1}.grad,'unit','cm')
        ft_plot_sens(squideeg_timelocked{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
        hold off;
        title('SQUID-MEG')
        view([-140 10])
        savefig(h, fullfile(save_path, 'figs', 'meg_layout2.fig'))
        saveas(h, fullfile(save_path, 'figs', 'meg_layout2.jpg'))

        close all

        %% Save
        save(fullfile(save_path, [params.sub '_opm_timelockedT']), 'opm_timelockedT', '-v7.3');
        save(fullfile(save_path, [params.sub '_opmeeg_timelockedT']), 'opmeeg_timelockedT', '-v7.3');
        save(fullfile(save_path, [params.sub '_sourcemodel']), 'sourcemodel', '-v7.3');
        save(fullfile(save_path, [params.sub '_sourcemodel_inflated']), 'sourcemodel_inflated', '-v7.3');

    else
        disp(['Required files not found. No transformed OPM/sourcemodel data was saved for ' params.sub])
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
        disp(['Not overwriting dipole source reconstruction for ' params.sub]);
    else
        clear headmodels mri_resliced
        headmodels = load(fullfile(save_path, 'headmodels.mat')).headmodels;
        mri_resliced = load(fullfile(save_path, 'mri_resliced.mat')).mri_resliced;
        
        clear squimdag_timelocked squidgrad_timelocked opm_timelockedT
        opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
        squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat'])).timelocked;
        squidgrad_timelocked = load(fullfile(save_path, [params.sub '_squidgrad_timelocked.mat'])).timelocked;
        clear M60_opm M60_squidmag M60_squidgrad
        M60_opm = load(fullfile(save_path, [params.sub '_opm_M60'])).M60; 
        M60_squidmag = load(fullfile(save_path, [params.sub '_squidmag_M60'])).M60; 
        M60_squidgrad = load(fullfile(save_path, [params.sub '_squidgrad_M60'])).M60; 
        [squidmag_dipole, squidgrad_dipole, opm_dipole] = fit_dipoles(save_path, squidmag_timelocked, squidgrad_timelocked, opm_timelockedT, headmodels, mri_resliced, M60_squidmag, M60_squidgrad, M60_opm, params);
    end
end

%% Dipole group analysis
subs = [2:13];
dipole_results_goup(base_save_path,subs, params)

%% Empty room & resting state for noise covariances
for i_sub = 2%2:size(subses,1)
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    save_path = fullfile(base_save_path,params.sub);
    if i_sub <=3 % Flip amplitudes in old recordings
        params.sign = -1;
    else
        params.sign = 1;
    end
    if exist(fullfile(save_path, 'opm_resting_state_ica.mat'),'file') && overwrite.empty_room == false
        disp(['Not overwriting MNE source reconstruction for ' params.sub]);
    else
        data_ica = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT{1};
        opm_chs = data_ica.label(contains(data_ica.label,'bz'));
        clear data_ica
        data_ica = load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat'])).timelocked{1};
        squidmag_chs = data_ica.label(contains(data_ica.label,'MEG'));
        clear data_ica
        data_ica = load(fullfile(save_path, [params.sub '_squidgrad_timelocked.mat'])).timelocked{1};
        squidgrad_chs = data_ica.label(contains(data_ica.label,'MEG'));
        clear data_ica
        squid_chs = [squidmag_chs; squidgrad_chs];
    
        % Empty room
        opm_file = fullfile(raw_path, 'osmeg', 'EmptyRoomOPM_raw.fif');
        squid_file = fullfile(raw_path, 'meg', 'EmptyRoomMEG_tsss.fif');
        if exist(opm_file,'file') && exist(squid_file,'file')
            read_empty_rooms(opm_file, squid_file, opm_chs, squid_chs, save_path, params);
        end
    
        % RESO
        opm_file = fullfile(raw_path, 'osmeg', 'RSEOOPM_raw.fif');
        aux_file = fullfile(raw_path, 'meg', 'RSEOEEG.fif');
        read_osMEG_RS(opm_file, aux_file, opm_chs, save_path, params)
        squid_file = fullfile(raw_path, 'meg', 'RSEOMEG_tsss.fif');
        read_cvMEG_RS(squid_file, squid_chs, save_path, params)
    end
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
        disp(['Not overwriting MNE source reconstruction for ' params.sub]);
    else
        clear headmodels sourcemodel sourcemodel_inflated
        sourcemodel = load(fullfile(save_path, [params.sub '_sourcemodel'])).sourcemodel;
        sourcemodel_inflated = load(fullfile(save_path, [params.sub '_sourcemodel_inflated'])).sourcemodel_inflated;
        headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;
        
        clear squimdag_timelocked squidgrad_timelocked opm_timelockedT
        opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
        squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat'])).timelocked;
        squidgrad_timelocked = load(fullfile(save_path, [params.sub '_squidgrad_timelocked.mat'])).timelocked;

        %params.use_cov_all = true;

        fit_mne(save_path, squidmag_timelocked, squidgrad_timelocked, opm_timelockedT, headmodels, sourcemodel, sourcemodel_inflated, params);
    end
end
close all

%% MNE group analysis
subs = [2:10 12:13];
mne_results_goup(base_save_path, subs, params);

%%
close all
clear all
exit
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
overwrite.preproc = true;
overwrite.coreg = false;
overwrite.mri = false;
overwrite.dip = false;
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
    if exist(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat']),'file') && overwrite.preproc==false
        disp(['Not overwriting OPM preproc for ' params.sub]);
    else
        ft_hastoolbox('mne', 1);

        % Read data
        [opm_cleaned, opmeeg_cleaned] = read_osMEG(opm_file, aux_file, save_path, params); % Read data
        
        if i_sub <=3 % Flip amplitudes in old recordings
            chs = find(contains(opm_cleaned.label,'bz'));
            for i_trl = 1:length(opm_cleaned.trial)
                opm_cleaned.trial{i_trl}(chs,:) = -opm_cleaned.trial{i_trl}(chs,:);
            end
        end
        opm_cleaned.grad = ft_convert_units(opm_cleaned.grad,'cm');

        % ICA
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        opm_ica = ica_MEG(opm_cleaned, save_path, params);
        clear opm_cleaned

        cfg = [];
        cfg.elec = opmeeg_cleaned.elec;
        cfg.output = fullfile(save_path, [params.sub '_opmeeg_layout.mat']);
        opmeeg_layout = ft_prepare_layout(cfg);
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        params.modality = 'opmeeg';
        opmeeg_ica = ica_MEG(opmeeg_cleaned, save_path, params);
        close all
        clear opmeeg_cleaned

        % Average
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        opm_timelocked = timelock_MEG(opm_ica, save_path, params);
        close all
        clear opm_ica

        params.modality = 'opmeeg';
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        opmeeg_timelocked = timelock_MEG(opmeeg_ica, save_path, params);
        close all
        clear opmeeg_ica
        clear opm_timelocked opmeeg_timelocked

        params = rmfield(params,{'modality', 'layout', 'chs', 'amp_scaler', 'amp_label'}); % remove fields used for picking modality
    end

    %% SQUID-MEG 
    if exist(fullfile(save_path, [params.sub '_squideeg_timelocked.mat']),'file') && overwrite.preproc==false
        disp(['Not overwriting SQUID preproc for ' params.sub]);
    else
        ft_hastoolbox('mne', 1);

        % Read data
        [squid_cleaned, squideeg_cleaned] = read_cvMEG(meg_file, save_path, params); % Read data
        
        % SQUID-MAG ICA
        params.modality = 'squidmag';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'megmag';
        squidmag_ica = ica_MEG(squid_cleaned, save_path, params);      

        % SQUID-MAG timelock
        params.modality = 'squidmag';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'megmag';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        squid_timelocked = timelock_MEG(squidmag_ica, save_path, params);
        close all
        clear squidmag_ica squid_timelocked

        % SQUID-GRAD timelock
        params.modality = 'squidgrad';
        params.layout = 'neuromag306planar.lay';
        params.chs = 'meggrad';
        squidgrad_ica = ica_MEG(squid_cleaned, save_path, params);
        clear squid_cleaned 

        % SQUID-GRAD timelock
        params.modality = 'squidgrad';
        params.layout = 'neuromag306planar.lay';
        params.chs = 'meggrad';
        params.amp_scaler = 1e15/100;
        params.amp_label = 'B [fT/cm]';
        squidgrad_timelocked = timelock_MEG(squidgrad_ica, save_path, params);
        close all
        clear squidgrad_ica squidgrad_timelocked

        % EEG ICA
        cfg = [];
        cfg.elec = squideeg_cleaned.elec;
        cfg.output = fullfile(save_path, [params.sub '_megeeg_layout.mat']);
        megeeg_layout = ft_prepare_layout(cfg);
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        params.modality = 'squideeg';
        squideeg_ica = ica_MEG(squideeg_cleaned, save_path, params);
        close all
        clear squideeg_cleaned

        % EEG timelock
        params.modality = 'squideeg';
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        squideeg_timelocked = timelock_MEG(squideeg_ica, save_path, params);
        close all
        clear squideeg_íca
        clear squideeg_timelocked

        params = rmfield(params,{'modality', 'layout', 'chs', 'amp_scaler', 'amp_label'}); % remove fields used for picking modality    
    end

    create_bads_reports(base_save_path, i_sub, params);
    close all

end

% Sensor level group analysis
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
        clear squidmag_timelocked opm_timelocked opmeeg_timelocked
        clear headmodels meshes filename
        load(fullfile(save_path,'headmodels.mat'));
        load(fullfile(save_path,'meshes.mat'));    
        load(fullfile(save_path, 'mri_resliced.mat'));
        load(fullfile(save_path, 'opm_trans.mat'));
        load(fullfile(save_path, [params.sub '_opm_timelocked.mat']))
        opm_timelockedT = timelocked;
        clear timelocked
        load(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat']))
        opmeeg_timelockedT = timelocked;
        clear timelocked
        load(fullfile(save_path, [params.sub '_squideeg_timelocked.mat']))
        squideeg_timelocked = timelocked;
        clear timelocked
        load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat']))
        squidmag_timelocked = timelocked;
        clear timelocked;

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
            opmeeg_timelockedT{i}.elec.chanpos = squideeg_timelockedT{i}.elec.chanpos;
            opmeeg_timelockedT{i}.elec.elecpos = squideeg_timelockedT{i}.elec.elecpos;
        end

        % Read and transform cortical restrained source model
        files = dir(fullfile(mri_path,'workbench'));
        if i_sub ==5
            files = dir(fullfile(save_path,'wb'));
        end
        for i = 1:length(files)
            if endsWith(files(i).name,'.L.midthickness.4k_fs_LR.surf.gii')
                filename = fullfile(mri_path,'workbench',files(i).name);
            end
        end
        sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});
        T = mri_resliced.transform/mri_resliced.hdr.vox2ras;
        sourcemodel = ft_transform_geometry(T, sourcemodel);
        sourcemodel.inside = true(size(sourcemodel.pos,1),1);

        for i = 1:length(files)
            if endsWith(files(i).name,'.L.inflated.4k_fs_LR.surf.gii')
                filename = fullfile(mri_path,'workbench',files(i).name);
            end
        end
        sourcemodel_inflated = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});
        sourcemodel_inflated = ft_transform_geometry(T, sourcemodel_inflated);
        sourcemodel_inflated.inside = true(size(sourcemodel_inflated.pos,1),1);
        
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
        clear headmodels mri_resliced timelocked M60
        load(fullfile(save_path, 'headmodels.mat'));
        load(fullfile(save_path, 'mri_resliced.mat'));
        clear megmag_timelocked magplanar_timelocked squideeg_timelocked
        clear opm_timelockedT opmeeg_timelockedT
        load(fullfile(save_path, [params.sub '_opm_timelockedT.mat']));
        load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat']));
        squidmag_timelocked = timelocked;
        clear timelocked
        load(fullfile(save_path, [params.sub '_squidgrad_timelocked.mat']));
        squidgrad_timelocked = timelocked;
        clear timelocked
        load(fullfile(save_path, [params.sub '_opm_M60'])); 
        M60_opm = M60;
        clear M60
        load(fullfile(save_path, [params.sub '_squidmag_M60'])); 
        M60_squidmag = M60;
        clear M60
        load(fullfile(save_path, [params.sub '_squidgrad_M60'])); 
        M60_squidgrad = M60;
        clear M60
        [squidmag_dipole, squidgrad_dipole, opm_dipole] = fit_dipoles(save_path, squidmag_timelocked, squidgrad_timelocked, opm_timelockedT, headmodels, mri_resliced, M60_squidmag, M60_squidgrad, M60_opm, params);
    end
end

%% Dipole group analysis
subs = [2:13];
dipole_results_goup(base_save_path,subs, params)

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
        load(fullfile(save_path, [params.sub '_sourcemodel']));
        load(fullfile(save_path, [params.sub '_sourcemodel_inflated']));
        load(fullfile(save_path,'headmodels.mat'));
        clear squidmag_timelocked squidgrad_timelocked OPM_timelockedT
        clear timelocked
        load(fullfile(save_path, [params.sub '_opm_timelockedT.mat']));
        load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat']));
        squidmag_timelocked = timelocked;
        clear timelocked
        load(fullfile(save_path, [params.sub '_squidgrad_timelocked.mat']));
        squidgrad_timelocked = timelocked;
        clear timelocked
        clear M60
        load(fullfile(save_path, [params.sub '_opm_M60'])); 
        M60_opm = M60;
        clear M60
        load(fullfile(save_path, [params.sub '_squidmag_M60'])); 
        M60_squidmag = M60;
        clear M60
        load(fullfile(save_path, [params.sub '_squidgrad_M60'])); 
        M60_squidgrad = M60;
        clear M60
        fit_mne(save_path, squidmag_timelocked, squidgrad_timelocked, opm_timelockedT, headmodels, sourcemodel, sourcemodel_inflated, params);
    end
end
close all

%% MNE group analysis
subs = [2:10 12:13];
mne_results_goup(base_save_path, subs, params);


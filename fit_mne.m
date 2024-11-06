function [megmag_mne, megplanar_mne, opm_mne, megeeg_mne, opmeeg_mne, FAHM] = fit_mne(save_path,meg_timelocked,opm_timelocked,headmodels,sourcemodel,params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

colors = [[0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; [0.6350 0.0780 0.1840]];
latency_m100 = 0.105;

%% Prepare leadfields
cfg = [];
cfg.grad             = meg_timelocked{1}.grad;              % sensor positions
cfg.channel          = 'megmag';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodels.headmodel_meg;          % volume conduction model
leadfield_megmag = ft_prepare_leadfield(cfg,meg_timelocked{1});

cfg = [];
cfg.grad             = meg_timelocked{1}.grad;              % sensor positions
cfg.channel          = 'megplanar';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodels.headmodel_meg;          % volume conduction model
leadfield_megplanar = ft_prepare_leadfield(cfg,meg_timelocked{1});

cfg = [];
cfg.grad             = opm_timelocked{1}.grad;              % sensor positions
cfg.channel          = '*bz';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodels.headmodel_meg;          % volume conduction model
leadfield_opm = ft_prepare_leadfield(cfg,opm_timelocked{1});

cfg = [];
cfg.elec             = meg_timelocked{1}.elec;              % sensor positions
cfg.channel          = 'eeg';                  % the used channels
cfg.senstype         = 'eeg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodels.headmodel_eeg;          % volume conduction model
leadfield_megeeg = ft_prepare_leadfield(cfg,meg_timelocked{1});

cfg = [];
cfg.elec             = opm_timelocked{1}.elec;              % sensor positions
cfg.channel          = 'eeg';                  % the used channels
cfg.senstype         = 'eeg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodels.headmodel_eeg;          % volume conduction model
leadfield_opmeeg = ft_prepare_leadfield(cfg,opm_timelocked{1});

%% Fit dipoles
for i_phalange = 1:5
    % MEG
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodels.headmodel_meg;    % supply the headmodel
    %cfg.sourcemodel         = sourcemodel;
    %cfg.sourcemodel.leadfield = leadfield_megmag.leadfield;
    cfg.sourcemodel         = leadfield_megmag;
    cfg.senstype            = 'meg';            % sensor type
    cfg.channel             = 'megmag';         % which channels to use
    megmag_mne{i_phalange} = ft_sourceanalysis(cfg, meg_timelocked{i_phalange});
    megmag_mne{i_phalange}.tri = sourcemodel.tri;
    cfg.channel             = 'meggrad';            % which channels to use
    %cfg.sourcemodel.leadfield = leadfield_megplanar.leadfield;
    cfg.sourcemodel         = leadfield_megplanar;
    megplanar_mne{i_phalange} = ft_sourceanalysis(cfg, meg_timelocked{i_phalange});
    megplanar_mne{i_phalange}.tri = sourcemodel.tri;
    cfg.headmodel           = headmodels.headmodel_eeg;    % supply the headmodel
    cfg.senstype            = 'eeg';            % sensor type
    cfg.channel             = 'eeg';         % which channels to use
    %cfg.sourcemodel.leaddield = leadfield_megeeg.leadfield;
    cfg.sourcemodel         = leadfield_megeeg;
    megeeg_mne{i_phalange} = ft_sourceanalysis(cfg, meg_timelocked{i_phalange});
    megeeg_mne{i_phalange}.tri = sourcemodel.tri;
    
    % OPM
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodels.headmodel_meg;    % supply the headmodel
    %cfg.sourcemodel         = sourcemodel;
    %cfg.sourcemodel.leadfield = leadfield_opm.leadfield;
    cfg.sourcemodel         = leadfield_opm;
    cfg.senstype            = 'meg';            % sensor type
    cfg.channel             = '*bz';         % which channels to use
    opm_mne{i_phalange} = ft_sourceanalysis(cfg, opm_timelocked{i_phalange});
    opm_mne{i_phalange}.tri = sourcemodel.tri;
    cfg.headmodel           = headmodels.headmodel_eeg;    % supply the headmodel
    cfg.senstype            = 'eeg';            % sensor type
    cfg.channel             = 'eeg';         % which channels to use
    %cfg.sourcemodel.leadfield = leadfield_opmeeg.leadfield;
    cfg.sourcemodel         = leadfield_opmeeg;
    opmeeg_mne{i_phalange} = ft_sourceanalysis(cfg, opm_timelocked{i_phalange});
    opmeeg_mne{i_phalange}.tri = sourcemodel.tri;
    

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    % Change for better color options
    cfg.latency         = latency_m100;     % The time-point to plot (s)
    cfg.colorbar        = 'no';

    %% Activation size (Full Area - Half Max)
    % Find max activation at latency_m100, find all vertices >= half max,
    % find all triangles that include those vertices, calculate their area
    % and divide by 3 (if all three vertices are connected to a triangle it
    % will sum to the whole area).
    FAHM{i_phalange} = [];

    tmp = megmag_mne{i_phalange};
    [~,i_latency] = min(abs(tmp.time-latency_m100));
    half_max = max(tmp.avg.pow(:,i_latency))/2;
    i_vertices = find(tmp.avg.pow(:,i_latency)>=half_max);
    [triangles,~] = find(ismember(tmp.tri,i_vertices)); 
    triangles = tmp.tri(triangles,:);
    FAHM{i_phalange}.megmag = sum(calculateTriangleAreas(tmp.pos, triangles))/3;

    tmp = megplanar_mne{i_phalange};
    [~,i_latency] = min(abs(tmp.time-latency_m100));
    half_max = max(tmp.avg.pow(:,i_latency))/2;
    i_vertices = find(tmp.avg.pow(:,i_latency)>=half_max);
    [triangles,~] = find(ismember(tmp.tri,i_vertices)); 
    triangles = tmp.tri(triangles,:);
    FAHM{i_phalange}.megplanar = sum(calculateTriangleAreas(tmp.pos, triangles))/3;

    tmp = opm_mne{i_phalange};
    [~,i_latency] = min(abs(tmp.time-latency_m100));
    half_max = max(tmp.avg.pow(:,i_latency))/2;
    i_vertices = find(tmp.avg.pow(:,i_latency)>=half_max);
    [triangles,~] = find(ismember(tmp.tri,i_vertices)); 
    triangles = tmp.tri(triangles,:);
    FAHM{i_phalange}.opm = sum(calculateTriangleAreas(tmp.pos, triangles))/3;

    tmp = opmeeg_mne{i_phalange};
    [~,i_latency] = min(abs(tmp.time-latency_m100));
    half_max = max(tmp.avg.pow(:,i_latency))/2;
    i_vertices = find(tmp.avg.pow(:,i_latency)>=half_max);
    [triangles,~] = find(ismember(tmp.tri,i_vertices)); 
    triangles = tmp.tri(triangles,:);
    FAHM{i_phalange}.opmeeg = sum(calculateTriangleAreas(tmp.pos, triangles))/3;

    tmp = megeeg_mne{i_phalange};
    [~,i_latency] = min(abs(tmp.time-latency_m100));
    half_max = max(tmp.avg.pow(:,i_latency))/2;
    i_vertices = find(tmp.avg.pow(:,i_latency)>=half_max);
    [triangles,~] = find(ismember(tmp.tri,i_vertices)); 
    triangles = tmp.tri(triangles,:);
    FAHM{i_phalange}.megeeg = sum(calculateTriangleAreas(tmp.pos, triangles))/3;

    %% Plots
    h = figure;
    ft_sourceplot(cfg, opm_mne{i_phalange})
    title(['OPM (FAHM=' num2str(FAHM{i_phalange}.opm,3) ')'])
    savefig(h, fullfile(save_path,'figs', ['opm_mne_ph' params.phalange_labels{i_phalange} '.fig']))
    saveas(h, fullfile(save_path,'figs', ['opm_mne_ph' params.phalange_labels{i_phalange} '.jpg']))

    h = figure;
    ft_sourceplot(cfg, megmag_mne{i_phalange})
    title(['MEG-MAG (FAHM=' num2str(FAHM{i_phalange}.megmag,3) ')'])
    savefig(h, fullfile(save_path,'figs', ['megmag_mne_ph' params.phalange_labels{i_phalange} '.fig']))
    saveas(h, fullfile(save_path,'figs', ['megmag_mne_ph' params.phalange_labels{i_phalange} '.jpg']))

    h = figure;
    ft_sourceplot(cfg, megplanar_mne{i_phalange})
    title(['MEG-PLANAR (FAHM=' num2str(FAHM{i_phalange}.megplanar,3) ')'])
    savefig(h, fullfile(save_path,'figs', ['megplanar_mne_ph' params.phalange_labels{i_phalange} '.fig']))
    saveas(h, fullfile(save_path,'figs', ['megplanar_mne_ph' params.phalange_labels{i_phalange} '.jpg']))

    h = figure;
    ft_sourceplot(cfg, opmeeg_mne{i_phalange})
    title(['OPM-EEG (FAHM=' num2str(FAHM{i_phalange}.opmeeg,3) ')'])
    savefig(h, fullfile(save_path,'figs', ['opmeeg_mne_ph' params.phalange_labels{i_phalange} '.fig']))
    saveas(h, fullfile(save_path,'figs', ['opmeeg_mne_ph' params.phalange_labels{i_phalange} '.jpg']))

    h = figure;
    ft_sourceplot(cfg, megeeg_mne{i_phalange})
    title(['MEG-EEG (FAHM=' num2str(FAHM{i_phalange}.megeeg,3) ')'])
    savefig(h, fullfile(save_path,'figs', ['megeeg_mne_ph' params.phalange_labels{i_phalange} '.fig']))
    saveas(h, fullfile(save_path,'figs', ['megeeg_mne_ph' params.phalange_labels{i_phalange} '.jpg']))

    close all

end

%% Save
save(fullfile(save_path, 'mne_fits'), 'megmag_mne', 'megplanar_mne', 'opm_mne', 'megeeg_mne', 'opmeeg_mne', 'FAHM'); disp('done');

end
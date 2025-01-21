function [megmag_mne, megplanar_mne, opm_mne, megeeg_mne, opmeeg_mne, FAHM] = fit_mne(save_path,megmag_timelocked,megplanar_timelocked,megeeg_timelocked,opm_timelocked, opmeeg_timelocked,headmodels,sourcemodel,latency,params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

colors = [[0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; [0.6350 0.0780 0.1840]];

%% Prepare leadfields
cfg = [];
cfg.grad             = megmag_timelocked{1}.grad;              % sensor positions
cfg.channel          = 'megmag';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodels.headmodel_meg;          % volume conduction model
leadfield_megmag = ft_prepare_leadfield(cfg,megmag_timelocked{1});

cfg = [];
cfg.grad             = megplanar_timelocked{1}.grad;              % sensor positions
cfg.channel          = 'megplanar';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodels.headmodel_meg;          % volume conduction model
leadfield_megplanar = ft_prepare_leadfield(cfg,megplanar_timelocked{1});

cfg = [];
cfg.grad             = opm_timelocked{1}.grad;              % sensor positions
cfg.channel          = '*bz';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodels.headmodel_meg;          % volume conduction model
leadfield_opm = ft_prepare_leadfield(cfg,opm_timelocked{1});

if ~isempty(headmodels.headmodel_eeg)
    cfg = [];
    cfg.elec             = megeeg_timelocked{1}.elec;              % sensor positions
    cfg.channel          = 'eeg';                  % the used channels
    cfg.senstype         = 'eeg';            % sensor type
    cfg.grid.pos         = sourcemodel.pos;           % source points
    cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
    cfg.headmodel        = headmodels.headmodel_eeg;          % volume conduction model
    leadfield_megeeg = ft_prepare_leadfield(cfg,megeeg_timelocked{1});
    
    cfg = [];
    cfg.elec             = opmeeg_timelocked{1}.elec;              % sensor positions
    cfg.channel          = 'eeg';                  % the used channels
    cfg.senstype         = 'eeg';            % sensor type
    cfg.grid.pos         = sourcemodel.pos;           % source points
    cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
    cfg.headmodel        = headmodels.headmodel_eeg;          % volume conduction model
    leadfield_opmeeg = ft_prepare_leadfield(cfg,opmeeg_timelocked{1});
end

%% Fit dipoles
megmag_mne = [];
megmag_mne.avg = cell(5,1);
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
    tmp{i_phalange} = ft_sourceanalysis(cfg, megmag_timelocked{i_phalange});
    megmag_mne.avg{i_phalange} = [];
    megmag_mne.avg{i_phalange}.pow = tmp.avg.pow;
    megmag_mne.avg{i_phalange}.mom = tmp.avg.mom;
end
megmag_mne.time = tmp.time;
megmag_mne.cfg = tmp.cfg;
megmag_mne.method = tmp.method;
megmag_mne.pos = tmp.pos;
megmag_mne.tri = sourcemodel.tri;
clear tmp

megplanar_mne = [];
megplanar_mne.avg = cell(5,1);
for i_phalange = 1:5
    cfg.channel             = 'meggrad';            % which channels to use
    %cfg.sourcemodel.leadfield = leadfield_megplanar.leadfield;
    cfg.sourcemodel         = leadfield_megplanar;
    tmp = ft_sourceanalysis(cfg, megplanar_timelocked{i_phalange});
    megplanar_mne.avg{i_phalange} = [];
    megplanar_mne.avg{i_phalange}.pow = tmp.avg.pow;
    megplanar_mne.avg{i_phalange}.mom = tmp.avg.mom;
end
megplanar_mne.time = tmp.time;
megplanar_mne.cfg = tmp.cfg;
megplanar_mne.method = tmp.method;
megplanar_mne.pos = tmp.pos;
megplanar_mne.tri = sourcemodel.tri;
clear tmp

megeeg_mne = [];
megeeg_mne.avg = cell(5,1);
for i_phalange = 1:5
    megeeg_mne.avg{i_phalange} = [];
    if ~isempty(headmodels.headmodel_eeg)
        cfg.headmodel           = headmodels.headmodel_eeg;    % supply the headmodel
        cfg.senstype            = 'eeg';            % sensor type
        cfg.channel             = 'eeg';         % which channels to use
        %cfg.sourcemodel.leaddield = leadfield_megeeg.leadfield;
        cfg.sourcemodel         = leadfield_megeeg;
        tmp = ft_sourceanalysis(cfg, megeeg_timelocked{i_phalange});
        megeeg_mne.avg{i_phalange}.pow = tmp.avg.pow;
        megeeg_mne.avg{i_phalange}.mom = tmp.avg.mom;
        megeeg_mne.time = tmp.time;
        megeeg_mne.cfg = tmp.cfg;
        megeeg_mne.method = tmp.method;
        megeeg_mne.pos = tmp.pos;
        megeeg_mne.tri = sourcemodel.tri;
    end
end
clear tmp

% OPM
opm_mne = [];
opm_mne.avg = cell(5,1);
for i_phalange = 1:5
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
    tmp = ft_sourceanalysis(cfg, opm_timelocked{i_phalange});
    opm_mne.avg{i_phalange} = [];
    opm_mne.avg{i_phalange}.pow = tmp.avg.pow;
    opm_mne.avg{i_phalange}.mom = tmp.avg.mom;
end
opm_mne.time = tmp.time;
opm_mne.cfg = tmp.cfg;
opm_mne.method = tmp.method;
opm_mne.pos = tmp.pos;
opm_mne.tri = sourcemodel.tri;
clear tmp


opmeeg_mne = [];
opmeeg_mne.avg = cell(5,1);
for i_phalange = 1:5
    opmeeg_mne.avg{i_phalange} = [];
    if ~isempty(headmodels.headmodel_eeg)
        cfg.headmodel           = headmodels.headmodel_eeg;    % supply the headmodel
        cfg.senstype            = 'eeg';            % sensor type
        cfg.channel             = 'eeg';         % which channels to use
        %cfg.sourcemodel.leadfield = leadfield_opmeeg.leadfield;
        cfg.sourcemodel         = leadfield_opmeeg;
        tmp = ft_sourceanalysis(cfg, opmeeg_timelocked{i_phalange});
        opmeeg_mne.avg{i_phalange}.pow = tmp.avg.pow;
        opmeeg_mne.avg{i_phalange}.mom = tmp.avg.mom;
        opmeeg_mne.time = tmp.time;
        opmeeg_mne.cfg = tmp.cfg;
        opmeeg_mne.method = tmp.method;
        opmeeg_mne.pos = tmp.pos;
        opmeeg_mne.tri = sourcemodel.tri;
    end
end
clear tmp

for i_phalange = 1:5
    %% Activation size (Full Area - Half Max)
    % Find max activation at latency_m100, find all vertices >= half max,
    % find all triangles tat include those vertices, calculate their area
    % and divide by 3 (if all three vertices are connected to a triangle it
    % will sum to the whole area).
    FAHM{i_phalange} = [];

    tmp = megmag_mne;
    [~,i_latency] = min(abs(tmp.time-latency{i_phalange}.megmag));
    half_max = max(tmp.avg{i_phalange}.pow(:,i_latency))/2;
    i_vertices = find(tmp.avg{i_phalange}.pow(:,i_latency)>=half_max);
    [triangles,~] = find(ismember(tmp.tri,i_vertices)); 
    triangles = tmp.tri(triangles,:);
    FAHM{i_phalange}.megmag = sum(calculateTriangleAreas(tmp.pos, triangles))/3;

    tmp = megplanar_mne;
    [~,i_latency] = min(abs(tmp.time-latency{i_phalange}.megplanar));
    half_max = max(tmp.avg{i_phalange}.pow(:,i_latency))/2;
    i_vertices = find(tmp.avg{i_phalange}.pow(:,i_latency)>=half_max);
    [triangles,~] = find(ismember(tmp.tri,i_vertices)); 
    triangles = tmp.tri(triangles,:);
    FAHM{i_phalange}.megplanar = sum(calculateTriangleAreas(tmp.pos, triangles))/3;

    tmp = opm_mne;
    [~,i_latency] = min(abs(tmp.time-latency{i_phalange}.opm));
    half_max = max(tmp.avg{i_phalange}.pow(:,i_latency))/2;
    i_vertices = find(tmp.avg{i_phalange}.pow(:,i_latency)>=half_max);
    [triangles,~] = find(ismember(tmp.tri,i_vertices)); 
    triangles = tmp.tri(triangles,:);
    FAHM{i_phalange}.opm = sum(calculateTriangleAreas(tmp.pos, triangles))/3;

    if ~isempty(headmodels.headmodel_eeg)
        tmp = opmeeg_mne;
        [~,i_latency] = min(abs(tmp.time-latency{i_phalange}.opmeeg));
        half_max = max(tmp.avg{i_phalange}.pow(:,i_latency))/2;
        i_vertices = find(tmp.avg{i_phalange}.pow(:,i_latency)>=half_max);
        [triangles,~] = find(ismember(tmp.tri,i_vertices)); 
        triangles = tmp.tri(triangles,:);
        FAHM{i_phalange}.opmeeg = sum(calculateTriangleAreas(tmp.pos, triangles))/3;
    
        tmp = megeeg_mne;
        [~,i_latency] = min(abs(tmp.time-latency{i_phalange}.megeeg));
        half_max = max(tmp.avg{i_phalange}.pow(:,i_latency))/2;
        i_vertices = find(tmp.avg{i_phalange}.pow(:,i_latency)>=half_max);
        [triangles,~] = find(ismember(tmp.tri,i_vertices)); 
        triangles = tmp.tri(triangles,:);
        FAHM{i_phalange}.megeeg = sum(calculateTriangleAreas(tmp.pos, triangles))/3;
    else
        FAHM{i_phalange}.opmeeg = NaN;
        FAHM{i_phalange}.megeeg = NaN;
    end

    %% Plots
    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    % Change for better color options
    cfg.colorbar        = 'no';

    tmp = opm_mne;
    tmp.avg = tmp.avg{i_phalange};
    cfg.latency         = latency{i_phalange}.opm;
    h = figure;
    ft_sourceplot(cfg, tmp)
    title(['OPM (FAHM=' num2str(FAHM{i_phalange}.opm,3) ')'])
    savefig(h, fullfile(save_path,'figs', ['opm_mne_ph' params.phalange_labels{i_phalange} '.fig']))
    saveas(h, fullfile(save_path,'figs', ['opm_mne_ph' params.phalange_labels{i_phalange} '.jpg']))

    tmp = megmag_mne;
    tmp.avg = tmp.avg{i_phalange};
    cfg.latency         = latency{i_phalange}.megmag;
    h = figure;
    ft_sourceplot(cfg, tmp)
    title(['MEG-MAG (FAHM=' num2str(FAHM{i_phalange}.megmag,3) ')'])
    savefig(h, fullfile(save_path,'figs', ['megmag_mne_ph' params.phalange_labels{i_phalange} '.fig']))
    saveas(h, fullfile(save_path,'figs', ['megmag_mne_ph' params.phalange_labels{i_phalange} '.jpg']))

    tmp = megplanar_mne;
    tmp.avg = tmp.avg{i_phalange};
    cfg.latency         = latency{i_phalange}.megplanar;
    h = figure;
    ft_sourceplot(cfg, tmp)
    title(['MEG-PLANAR (FAHM=' num2str(FAHM{i_phalange}.megplanar,3) ')'])
    savefig(h, fullfile(save_path,'figs', ['megplanar_mne_ph' params.phalange_labels{i_phalange} '.fig']))
    saveas(h, fullfile(save_path,'figs', ['megplanar_mne_ph' params.phalange_labels{i_phalange} '.jpg']))

    if ~isempty(headmodels.headmodel_eeg)
        tmp = opmeeg_mne;
        tmp.avg = tmp.avg{i_phalange};
        cfg.latency         = latency{i_phalange}.opmeeg;
        h = figure;
        ft_sourceplot(cfg, tmp)
        title(['OPM-EEG (FAHM=' num2str(FAHM{i_phalange}.opmeeg,3) ')'])
        savefig(h, fullfile(save_path,'figs', ['opmeeg_mne_ph' params.phalange_labels{i_phalange} '.fig']))
        saveas(h, fullfile(save_path,'figs', ['opmeeg_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    
        tmp = megeeg_mne;
        tmp.avg = tmp.avg{i_phalange};
        cfg.latency         = latency{i_phalange}.megeeg;
        h = figure;
        ft_sourceplot(cfg, tmp)
        title(['MEG-EEG (FAHM=' num2str(FAHM{i_phalange}.megeeg,3) ')'])
        savefig(h, fullfile(save_path,'figs', ['megeeg_mne_ph' params.phalange_labels{i_phalange} '.fig']))
        saveas(h, fullfile(save_path,'figs', ['megeeg_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    end

    close all
end

%% Save
save(fullfile(save_path, 'mne_fits'), 'megmag_mne', 'megplanar_mne', 'opm_mne', 'megeeg_mne', 'opmeeg_mne'); 
save(fullfile(save_path, 'mne_fahm'), 'FAHM'); disp('done');

end
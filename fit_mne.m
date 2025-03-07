function fit_mne(save_path, squidmag_timelocked, squidgrad_timelocked, opm_timelocked, headmodels, sourcemodel, sourcemodel_inflated, M60_squidmag, M60_squidgrad, M60_opm, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Prepare leadfields
headmodel = headmodels.headmodel_meg;
%headmodel.order = 20;
sourcemodel.mom = surface_normals(sourcemodel.pos, sourcemodel.tri, 'vertex')';

cfg = [];
cfg.grad             = squidmag_timelocked{1}.grad; % sensor positions
cfg.channel          = 'megmag';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.sourcemodel      = sourcemodel;           % source points
cfg.headmodel        = headmodel;          % volume conduction model
leadfield_squidmag = ft_prepare_leadfield(cfg,squidmag_timelocked{1});

cfg = [];
cfg.grad             = squidgrad_timelocked{1}.grad;              % sensor positions
cfg.channel          = 'megplanar';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.sourcemodel      = sourcemodel;           % source points
cfg.headmodel        = headmodel;          % volume conduction model
leadfield_squidgrad = ft_prepare_leadfield(cfg,squidgrad_timelocked{1});

cfg = [];
cfg.grad             = opm_timelocked{1}.grad;              % sensor positions
cfg.channel          = '*bz';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.sourcemodel      = sourcemodel;           % source points
cfg.headmodel        = headmodel;          % volume conduction model
leadfield_opm = ft_prepare_leadfield(cfg,opm_timelocked{1});

%% MNE invserse
% MEG-MAG
squidmag_mne = [];
squidmag_mne.avg = cell(5,1);
squidmag_mne_M60 = cell(5,1);
for i_phalange = 1:5
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield_squidmag;
    cfg.senstype            = 'meg';            % sensor type
    cfg.channel             = 'megmag';         % which channels to use
    tmp = ft_sourceanalysis(cfg, squidmag_timelocked{i_phalange});
    tmp.tri = sourcemodel.tri;
    %cfg = [];
    %cfg.projectmom = 'yes';
    %tmp = ft_sourcedescriptives(cfg,tmp);
    squidmag_mne.avg{i_phalange} = [];
    squidmag_mne.avg{i_phalange}.pow = tmp.avg.pow;
    squidmag_mne.avg{i_phalange}.mom = tmp.avg.mom;
    
    squidmag_mne_M60{i_phalange} = FullAreaHalfMax(tmp,sourcemodel);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = squidmag_mne_M60{i_phalange}.peak_latency;
    %tmp.pos = sourcemodel_inflated.pos;
    %tmp.tri = sourcemodel_inflated.tri;
    h = figure;
    ft_sourceplot(cfg, tmp)
    lighting gouraud
    material dull
    title(['SQUID-MAG (FAHM=' num2str(squidmag_mne_M60{i_phalange}.fahm,3) ')'])
    saveas(h, fullfile(save_path,'figs', [params.sub '_squidmag_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    close all
    %%
end
squidmag_mne.time = tmp.time;
squidmag_mne.cfg = tmp.cfg;
squidmag_mne.method = tmp.method;
squidmag_mne.pos = tmp.pos;
squidmag_mne.tri = sourcemodel.tri;

save(fullfile(save_path, 'squidmag_mne'), 'squidmag_mne'); 
save(fullfile(save_path, 'squidmag_mne_M60'), 'squidmag_mne_M60'); 
clear tmp squidmag_mne leadfield_squidmag

% MEG-GRAD
squidgrad_mne = [];
squidgrad_mne.avg = cell(5,1);
squidgrad_mne_M60 = cell(5,1);
for i_phalange = 1:5
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.senstype            = 'meg';
    cfg.channel             = 'megplanar';            % which channels to use
    cfg.sourcemodel         = leadfield_squidgrad;
    tmp = ft_sourceanalysis(cfg, squidgrad_timelocked{i_phalange});
    tmp.tri = sourcemodel.tri;
    %cfg = [];
    %cfg.projectmom = 'yes';
    %tmp = ft_sourcedescriptives(cfg,tmp);
    squidgrad_mne.avg{i_phalange} = [];
    squidgrad_mne.avg{i_phalange}.pow = tmp.avg.pow;
    squidgrad_mne.avg{i_phalange}.mom = tmp.avg.mom;
    
    squidgrad_mne_M60{i_phalange} = FullAreaHalfMax(tmp,sourcemodel);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = squidgrad_mne_M60{i_phalange}.peak_latency;
    %tmp.pos = sourcemodel_inflated.pos;
    %tmp.tri = sourcemodel_inflated.tri;
    h = figure;
    ft_sourceplot(cfg, tmp)
    lighting gouraud
    material dull
    title(['SQUID-GRAD (FAHM=' num2str(squidgrad_mne_M60{i_phalange}.fahm,3) ')'])
    saveas(h, fullfile(save_path,'figs', [params.sub '_squidgrad_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    close all
end
squidgrad_mne.time = tmp.time;
squidgrad_mne.cfg = tmp.cfg;
squidgrad_mne.method = tmp.method;
squidgrad_mne.pos = tmp.pos;
squidgrad_mne.tri = sourcemodel.tri;

save(fullfile(save_path, 'squidgrad_mne'), 'squidgrad_mne'); 
save(fullfile(save_path, 'squidgrad_mne_M60'), 'squidgrad_mne_M60'); 
clear tmp squidgrad_mne leadfield_squidgrad

% OPM
opm_mne = [];
opm_mne.avg = cell(5,1);
opm_mne_M60 = cell(5,1);
for i_phalange = 1:5
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield_opm;
    cfg.senstype            = 'meg';            % sensor type
    cfg.channel             = '*bz';         % which channels to use
    tmp = ft_sourceanalysis(cfg, opm_timelocked{i_phalange});
    tmp.tri = sourcemodel.tri;
    %cfg = [];
    %cfg.projectmom = 'yes';
    %tmp = ft_sourcedescriptives(cfg,tmp);
    opm_mne.avg{i_phalange} = [];
    opm_mne.avg{i_phalange}.pow = tmp.avg.pow;
    opm_mne.avg{i_phalange}.mom = tmp.avg.mom;

    opm_mne_M60{i_phalange} = FullAreaHalfMax(tmp,sourcemodel);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = opm_mne_M60{i_phalange}.peak_latency;
    %tmp.pos = sourcemodel_inflated.pos;
    %tmp.tri = sourcemodel_inflated.tri;
    h = figure;
    ft_sourceplot(cfg, tmp)
    lighting gouraud
    material dull
    title(['OPM (FAHM=' num2str(opm_mne_M60{i_phalange}.fahm,3) ')'])
    saveas(h, fullfile(save_path,'figs', [params.sub '_opm_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    close all
end
opm_mne.time = tmp.time;
opm_mne.cfg = tmp.cfg;
opm_mne.method = tmp.method;
opm_mne.pos = tmp.pos;
opm_mne.tri = sourcemodel.tri;

save(fullfile(save_path, 'opm_mne'), 'opm_mne'); 
save(fullfile(save_path, 'opm_mne_M60'), 'opm_mne_M60'); 
clear tmp opm_mne leadfield_opm

end
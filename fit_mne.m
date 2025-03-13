function fit_mne(save_path, squidmag_timelocked, squidgrad_timelocked, opm_timelocked, headmodels, sourcemodel, sourcemodel_inflated, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Prepare leadfields
headmodel = headmodels.headmodel_meg;
sourcemodel.mom = surface_normals(sourcemodel.pos, sourcemodel.tri, 'vertex')';
sourcemodel.unit = 'cm';

%% MNE invserse
squidmag_mne_M60 = cell(5,1);
squidgrad_mne_M60 = cell(5,1);
opm_mne_M60 = cell(5,1);

for i_phalange = 1:5
    params.i_phalange = i_phalange;

    if isfield(params,'use_cov_all') 
        if params.use_cov_all
            squidmag_timelocked{i_phalange}.cov = squidmag_timelocked{i_phalange}.cov_all;
            squidgrad_timelocked{i_phalange}.cov = squidgrad_timelocked{i_phalange}.cov_all;
            opm_timelocked{i_phalange}.cov = opm_timelocked{i_phalange}.cov_all;
        end
    end
    
    %% MEG-MAG
    cfg = [];
    cfg.grad             = squidmag_timelocked{i_phalange}.grad; % sensor positions
    cfg.senstype         = 'meg';            % sensor type
    cfg.sourcemodel      = sourcemodel;           % source points
    cfg.headmodel        = headmodel;          % volume conduction model
    leadfield = ft_prepare_leadfield(cfg,squidmag_timelocked{i_phalange});

    cfg = [];
    cfg.method              = 'mne';
    %cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield;
    cfg.senstype            = 'meg';            % sensor type
    tmp = ft_sourceanalysis(cfg, squidmag_timelocked{i_phalange});
    tmp.tri = sourcemodel.tri;

    params.modality = 'squidmag';
    squidmag_mne_M60{i_phalange} = FullAreaHalfMax(tmp,sourcemodel,params, save_path);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = squidmag_mne_M60{i_phalange}.peak_latency;
    tmp.pos = sourcemodel_inflated.pos;
    tmp.tri = sourcemodel_inflated.tri;
    h = figure;
    h.Position(3) = round(h.Position(3)*1.2);
    ft_sourceplot(cfg, tmp)
    lighting gouraud
    material dull
    title(['SQUID-MAG (FAHM=' num2str(squidmag_mne_M60{i_phalange}.fahm,3) 'cm^2; t=' num2str(round(squidmag_mne_M60{i_phalange}.peak_latency*1e3)) 'ms)'])
    saveas(h, fullfile(save_path,'figs', [params.sub '_squidmag_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    close all
    
    clear tmp leadfield

    %% MEG-GRAD
    cfg = [];
    cfg.grad             = squidgrad_timelocked{i_phalange}.grad;              % sensor positions
    cfg.senstype         = 'meg';            % sensor type
    cfg.sourcemodel      = sourcemodel;           % source points
    cfg.headmodel        = headmodel;          % volume conduction model
    leadfield = ft_prepare_leadfield(cfg,squidgrad_timelocked{i_phalange});

    cfg = [];
    cfg.method              = 'mne';
    %cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.senstype            = 'meg';
    cfg.sourcemodel         = leadfield;
    tmp = ft_sourceanalysis(cfg, squidgrad_timelocked{i_phalange});
    tmp.tri = sourcemodel.tri;
    
    params.modality = 'squidgrad';
    squidgrad_mne_M60{i_phalange} = FullAreaHalfMax(tmp,sourcemodel,params, save_path);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = squidgrad_mne_M60{i_phalange}.peak_latency;
    tmp.pos = sourcemodel_inflated.pos;
    tmp.tri = sourcemodel_inflated.tri;
    h = figure;
    ft_sourceplot(cfg, tmp)
    lighting gouraud
    material dull
    title(['SQUID-GRAD (FAHM=' num2str(squidgrad_mne_M60{i_phalange}.fahm,3) 'cm^2; t=' num2str(round(squidgrad_mne_M60{i_phalange}.peak_latency*1e3)) 'ms)'])
    saveas(h, fullfile(save_path,'figs', [params.sub '_squidgrad_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    close all

    clear tmp leadfield

    %% OPM
    cfg = [];
    cfg.grad             = opm_timelocked{i_phalange}.grad;              % sensor positions
    cfg.senstype         = 'meg';            % sensor type
    cfg.sourcemodel      = sourcemodel;           % source points
    cfg.headmodel        = headmodel;          % volume conduction model
    leadfield = ft_prepare_leadfield(cfg,opm_timelocked{i_phalange});

    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield;
    cfg.senstype            = 'meg';            % sensor type
    tmp = ft_sourceanalysis(cfg, opm_timelocked{i_phalange});
    tmp.tri = sourcemodel.tri;
    
    params.modality = 'opm';
    opm_mne_M60{i_phalange} = FullAreaHalfMax(tmp,sourcemodel,params, save_path);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = opm_mne_M60{i_phalange}.peak_latency;
    tmp.pos = sourcemodel_inflated.pos;
    tmp.tri = sourcemodel_inflated.tri;
    h = figure;
    ft_sourceplot(cfg, tmp)
    lighting gouraud
    material dull
    title(['OPM (FAHM=' num2str(opm_mne_M60{i_phalange}.fahm,3) 'cm^2; t=' num2str(round(opm_mne_M60{i_phalange}.peak_latency*1e3)) 'ms)'])
    saveas(h, fullfile(save_path,'figs', [params.sub '_opm_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    close all

    tmp = [];
    clear tmp leadfield

    %% Overlaps
    i_vertices = opm_mne_M60{i_phalange}.halfmax_distribution & squidmag_mne_M60{i_phalange}.halfmax_distribution;
    [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
    triangles = sourcemodel.tri(triangles,:);
    opm_mne_M60{i_phalange}.overlap_squidmag = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
    squidmag_mne_M60{i_phalange}.overlap_opm = opm_mne_M60{i_phalange}.overlap_squidmag;

    i_vertices = opm_mne_M60{i_phalange}.halfmax_distribution & squidgrad_mne_M60{i_phalange}.halfmax_distribution;
    [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
    triangles = sourcemodel.tri(triangles,:);
    opm_mne_M60{i_phalange}.overlap_squidgrad = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
    squidgrad_mne_M60{i_phalange}.overlap_opm = opm_mne_M60{i_phalange}.overlap_squidgrad;

    i_vertices = squidmag_mne_M60{i_phalange}.halfmax_distribution & squidgrad_mne_M60{i_phalange}.halfmax_distribution;
    [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
    triangles = sourcemodel.tri(triangles,:);
    squidmag_mne_M60{i_phalange}.overlap_squidgrad = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
    squidgrad_mne_M60{i_phalange}.overlap_squidmag = squidmag_mne_M60{i_phalange}.overlap_squidgrad;

end
save(fullfile(save_path, 'squidmag_mne_M60'), 'squidmag_mne_M60'); 
save(fullfile(save_path, 'squidgrad_mne_M60'), 'squidgrad_mne_M60'); 
save(fullfile(save_path, 'opm_mne_M60'), 'opm_mne_M60'); 

end
function [hpi_fit, opm_trans, dip_pos_tf] = fit_hpi(file, aux_file, save_path, params)
%prprocess_osMEG Read on-scalp MEG data for benchmarking
% recordings and combine with auxiliary TRIUX data/EEG. 
% Requires the following arguments:
% Path: containing save_path and meg_file
% Params: hpi_freq.

%% --- Read triggers ---
% OPM
cfg = [];
cfg.datafile        = file;
cfg.coordsys        = 'dewar';
cfg.coilaccuracy    = 0;
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.hpi_freq+10;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.hpi_freq-10;
cfg.padding         = 1;
cfg.padtype         = 'data';
raw = ft_preprocessing(cfg);

%% Epoch
cfg = [];
cfg.length = 0.25;
cfg.overlap = 0;
epo = ft_redefinetrial(cfg,raw);

cfg = [];
timelocked = ft_timelockanalysis(cfg,epo);
timelocked.avg = zeros(size(timelocked.avg,1),1);
timelocked.time = zeros(size(timelocked.time,1),1);

hpi_chs = find(contains(raw.label,'hpiin'));
hpi_labels = raw.label(hpi_chs);
hpi_trials = false(length(hpi_chs),length(epo.trial));
for trl = 1:length(epo.trial)
    hpi_trials(:,trl) = (max(epo.trial{trl}(hpi_chs,:),[],2)-min(epo.trial{trl}(hpi_chs,:),[],2))>1e-3;
end

for i = 1:length(hpi_chs)
    hpi_trl{i} = find(hpi_trials(i,:));
    hpi_trl{i} = hpi_trl{i}(3:end-2);
    if ~isempty(hpi_trl{i})
        hpi_on(i) = true;
    else
        hpi_on(i) = false;
    end
end

hpi_chs = hpi_chs(hpi_on);

%% Prepare for dipole grid search
[X,Y,Z] = meshgrid((min(epo.grad.chanpos(:,1))-0.02):0.01:(max(epo.grad.chanpos(:,1))+0.02),(min(epo.grad.chanpos(:,2))-0.02):0.01:(max(epo.grad.chanpos(:,2))+0.02),(min(epo.grad.chanpos(:,3))-0.02):0.01:(max(epo.grad.chanpos(:,3))+0.02));
pos = [X(:) Y(:) Z(:)];
addpath('/Users/christophpfeiffer/Dropbox/Mac/Documents/MATLAB/myFunctions/')
inside = insidePointcloud(pos,epo.grad.chanpos+epo.grad.chanori*5e-3);

%% Freq analysis
amp = zeros(size(epo.trial{1},1),length(hpi_chs));
dip_pos = zeros(length(hpi_chs),3);
for coil = 1:length(hpi_chs)
    % Lock-in
    R = zeros(size(epo.trial{hpi_trl{coil}(1)},1),length(hpi_trl{coil}));
    Theta = zeros(size(epo.trial{hpi_trl{coil}(1)},1),length(hpi_trl{coil}));
    for i_trl = 1:length(hpi_trl{coil})
        X = mean(cos(2*pi*params.hpi_freq*epo.time{hpi_trl{coil}(i_trl)}).*epo.trial{hpi_trl{coil}(i_trl)},2);
        Y = mean(sin(2*pi*params.hpi_freq*epo.time{hpi_trl{coil}(i_trl)}).*epo.trial{hpi_trl{coil}(i_trl)},2);
        tmp = complex(X,Y);
        R(:,i_trl) = abs(tmp);
        Theta(:,i_trl) = angle(tmp./repmat(tmp(hpi_chs(coil)),[length(tmp) 1]));
    end
    amp(:,coil) = mean(R,2);
    amp(abs(mean(Theta,2))>pi/2,coil) = -amp(abs(mean(Theta,2))>pi/2,coil);
    
    timelocked.avg = amp(:,coil);

    cfg = [];
    cfg.layout = 'fieldlinebeta2bz_helmet.mat'; 
    cfg.parameter = 'avg';
    cfg.channel = '*bz';
    h = figure; ft_topoplotER(cfg,timelocked); colorbar
    saveas(h, fullfile(save_path, 'figs', ['hpi_topo_coil-' num2str(coil) '.jpg']))

    %% Dipole fit
    cfg = [];
    cfg.method = 'infinite';
    headmodel = ft_prepare_headmodel(cfg);
    
    %[~, i_maxchan] = max(abs(mean(opm_fft{coil}.fourierspctrm(:,:,1),1)));
    cfg = [];
    %cfg.frequency       = 33;
    cfg.numdipoles      = 1;
    %cfg.unit            = 'cm';
    cfg.gridsearch      = 'yes';
    cfg.channel = '*bz';
    cfg.sourcemodel     = [];
    cfg.sourcemodel.pos = pos;
    cfg.sourcemodel.inside = inside;
    %cfg.xgrid           = (min(opm_fft{coil}.grad.chanpos(:,1))-0.01):0.01:(max(opm_fft{coil}.grad.chanpos(:,1))+0.01);
    %cfg.ygrid           = (min(opm_fft{coil}.grad.chanpos(:,2))-0.01):0.01:(max(opm_fft{coil}.grad.chanpos(:,2))+0.01);
    %cfg.zgrid           = (min(opm_fft{coil}.grad.chanpos(:,3))-0.01):0.01:(max(opm_fft{coil}.grad.chanpos(:,3))+0.01);
    cfg.nonlinear       = 'yes';
    cfg.headmodel       = headmodel;
%     hpi_fit{coil} = ft_dipolefitting(cfg,opm_fft{coil});
%     hpi_fit{coil}.dip.ori = mean(hpi_fit{coil}.dip.mom(:,1:(size(hpi_fit{coil}.dip.mom,2)/2)),2);
%     hpi_fit{coil}.dip.ori = hpi_fit{coil}.dip.ori/norm(hpi_fit{coil}.dip.ori);
%     hpi_fit{coil}.dip.gof = 1-mean(hpi_fit{coil}.dip.rv(1:(size(hpi_fit{coil}.dip.mom,2)/2)));
    hpi_fit{coil} = ft_dipolefitting(cfg,timelocked);
    hpi_fit{coil}.dip.ori = hpi_fit{coil}.dip.mom/norm(hpi_fit{coil}.dip.mom);
    hpi_fit{coil}.dip.gof = 1-hpi_fit{coil}.dip.rv;

    dip_pos(coil,:) = hpi_fit{coil}.dip.pos*1e2;
    dip_ori(coil,:) = hpi_fit{coil}.dip.ori;
end

%%
headshape = ft_read_headshape(aux_file);
hpi_polhemus = headshape.pos(find(contains(headshape.label,'hpi')),:);
fixed = pointCloud(hpi_polhemus);
moving = pointCloud(dip_pos);
[opm_trans, dip_pos_tf, dist] = pcregistericp(moving, fixed);

n_dipoles = size(dip_pos,1);
i_min = zeros(n_dipoles,1);
d_min = zeros(n_dipoles,1);
for i = 1:n_dipoles
    [d_min(i),i_min(i)] = min(vecnorm(hpi_polhemus-repmat(dip_pos_tf.Location(i,:),[n_dipoles 1]),2,2));
end

dip_ori_tf = dip_ori*opm_trans.Rotation;

epoT = epo;
epoT.grad.chanpos = opm_trans.transformPointsForward(epo.grad.chanpos*1e2)*1e-2;
epoT.grad.coilpos = opm_trans.transformPointsForward(epo.grad.coilpos*1e2)*1e-2;

%%
colors = [[0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; [0.6350 0.0780 0.1840]];

h = figure;
ft_plot_sens(epoT.grad,'unit','cm','DisplayName','senspos'); 
hold on 
for coil = 1:length(hpi_chs)
    quiver3(dip_pos_tf.Location(coil,1),dip_pos_tf.Location(coil,2),dip_pos_tf.Location(coil,3),dip_ori_tf(coil,1),dip_ori_tf(coil,2),dip_ori_tf(coil,3),'*','Color',colors(coil,:),'DisplayName',[hpi_labels{coil} ' (GOF=' num2str((hpi_fit{coil}.dip.gof)*100,'%.2f') '%)'],'LineWidth',2);
end
scatter3(hpi_polhemus(:,1),hpi_polhemus(:,2),hpi_polhemus(:,3),'r','DisplayName','polhemus'); 
hold off
title(['HPI fits (mean dist = ' num2str(dist*10) ' mm)'])
legend
saveas(h, fullfile(save_path, 'figs', 'hpi_fits.jpg'))

%% Save 
save(fullfile(save_path, 'hpi_fit'), 'hpi_fit'); disp('done');
save(fullfile(save_path, 'opm_trans'), 'opm_trans'); disp('done');

end
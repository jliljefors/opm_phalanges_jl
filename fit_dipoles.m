function [megmag_dipole, megplanar_dipole, opm_dipole, eeg_dipole, opmeeg_dipole] = fit_dipoles(save_path,meg_timelocked,opm_timelocked,headmodels,mri,params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

colors = [[0 0.4470 0.7410]; % blue
    [0.8500 0.3250 0.0980]; % red
    [0.9290 0.6940 0.1250]; % yellow
    [0.4940 0.1840 0.5560]; % purple
    [0.4660 0.6740 0.1880]; % green
    [0.6350 0.0780 0.1840]]; % light blue

latency_m100 = [0.08 0.13];

cfg              = [];
cfg.resolution   = 1;
cfg.tight        = 'yes';
cfg.inwardshift  = 0;
cfg.headmodel    = headmodels.headmodel_meg;
sourcemodel    = ft_prepare_sourcemodel(cfg);

%% Fit dipoles
for i_phalange = 1:5
    % MEG
    cfg = [];
    cfg.gridsearch      = 'yes';           
    cfg.numdipoles      = 1;                
    cfg.sourcemodel     = sourcemodel;           
    cfg.headmodel       = headmodels.headmodel_meg;    
    cfg.senstype        = 'meg';            
    cfg.channel         = 'megmag';         
    cfg.nonlinear       = 'yes';           
    cfg.latency         = latency_m100;   
    cfg.dipfit.checkinside = 'yes';
    %cfg.dipfit.noisecov = meg_timelocked{i_phalange}.cov;
    megmag_dipole{i_phalange} = ft_dipolefitting(cfg, meg_timelocked{i_phalange});
    cfg.channel         = 'meglpanar';           
    megplanar_dipole{i_phalange} = ft_dipolefitting(cfg, meg_timelocked{i_phalange});
    cfg.headmodel       = headmodels.headmodel_eeg;    
    cfg.senstype        = 'eeg';            
    cfg.channel         = 'eeg';        
    eeg_dipole{i_phalange} = ft_dipolefitting(cfg, meg_timelocked{i_phalange});

    % OPM
    cfg = [];
    cfg.gridsearch      = 'yes';           
    cfg.numdipoles      = 1;                
    cfg.sourcemodel     = sourcemodel;            
    cfg.headmodel       = headmodels.headmodel_meg;    
    cfg.senstype        = 'meg';            
    cfg.channel         = '*bz';        
    cfg.nonlinear       = 'yes';            
    cfg.latency         = latency_m100;   
    cfg.dipfit.checkinside = 'yes';
    %cfg.dipfit.noisecov = opm_timelocked{i_phalange}.cov;
    opm_dipole{i_phalange} = ft_dipolefitting(cfg, opm_timelocked{i_phalange});
    cfg.headmodel       = headmodels.headmodel_eeg;  
    cfg.senstype        = 'eeg';           
    cfg.channel         = 'eeg';        
    opmeeg_dipole{i_phalange} = ft_dipolefitting(cfg, opm_timelocked{i_phalange});

    % Plot OPM vs SQUID
    pos_mag = megmag_dipole{i_phalange}.dip.pos;
    [~,idx] = max(vecnorm(megmag_dipole{i_phalange}.dip.mom,2,1));
    ori_mag = megmag_dipole{i_phalange}.dip.mom(:,idx);

    pos_planar = megplanar_dipole{i_phalange}.dip.pos;
    [~,idx] = max(vecnorm(megplanar_dipole{i_phalange}.dip.mom,2,1));
    ori_planar = megplanar_dipole{i_phalange}.dip.mom(:,idx);
    
    pos_opm = opm_dipole{i_phalange}.dip.pos;
    %pos_opm = opm_trans.transformPointsInverse(pos_opm);
    [~,idx] = max(vecnorm(opm_dipole{i_phalange}.dip.mom,2,1));
    ori_opm = -opm_dipole{i_phalange}.dip.mom(:,idx);

    pos_eeg = eeg_dipole{i_phalange}.dip.pos;
    [~,idx] = max(vecnorm(eeg_dipole{i_phalange}.dip.mom,2,1));
    ori_eeg = eeg_dipole{i_phalange}.dip.mom(:,idx);

    pos_opmeeg = opmeeg_dipole{i_phalange}.dip.pos;
    %pos_opmeeg = opm_trans.transformPointsInverse(pos_opmeeg);
    [~,idx] = max(vecnorm(opmeeg_dipole{i_phalange}.dip.mom,2,1));
    ori_opmeeg = -opmeeg_dipole{i_phalange}.dip.mom(:,idx);

    h = figure;
    ft_plot_dipole(pos_mag,ori_mag,'color',colors(1,:))
    hold on;
    ft_plot_dipole(pos_opm,ori_opm,'color',colors(2,:))
    ft_plot_dipole(pos_planar,ori_planar,'color',colors(3,:))
    ft_plot_headmodel(headmodels.headmodel_meg,'EdgeAlpha',0,'FaceAlpha',0.3,'FaceColor',[229 194 152]/256,'unit','cm') 
    hold off
    title([params.phalange_labels{i_phalange} ' (SqMAG-OPM = ' num2str(norm(pos_mag-pos_opm)*10,'%.1f') 'mm / SqGRAD-OPM = ' num2str(norm(pos_planar-pos_opm)*10,'%.1f') 'mm)'])
    legend('SQUIDMAG','OPM','SQUIDPLANAR','brain')
    saveas(h, fullfile(save_path, 'figs', ['dipfit_SQUIDvOPM_ph-' params.phalange_labels{i_phalange} '.jpg']))
    close

    h = figure;
    ft_plot_dipole(pos_eeg,ori_eeg,'color',colors(1,:))
    hold on;
    ft_plot_dipole(pos_opmeeg,ori_opmeeg,'color',colors(2,:))
    ft_plot_headmodel(headmodels.headmodel_meg,'EdgeAlpha',0,'FaceAlpha',0.3,'FaceColor',[229 194 152]/256,'unit','cm') 
    hold off
    title([params.phalange_labels{i_phalange} ' (SqEEG-OpmEEG = ' num2str(norm((pos_eeg-pos_opmeeg))*10,'%.1f') 'mm)'])
    legend('SQUIDEEG','OPMEEG','brain')
    saveas(h, fullfile(save_path, 'figs', ['dipfit_EEGvOPMEEG_ph-' params.phalange_labels{i_phalange} '.jpg']))
    close
end
close all

%% Plot phalanges jointly
% SQUID
pos_mag = zeros(5,3);
ori_mag = zeros(5,3);
for i = 1:5
    pos_mag(i,:) = megmag_dipole{i}.dip.pos;
    [~,idx] = max(vecnorm(megmag_dipole{i}.dip.mom,2,1));
    ori_mag(i,:) = megmag_dipole{i}.dip.mom(:,idx);
end
mean_pos = mean(pos_mag,1);
tmp = pos_mag;
tmp(:,1) = mean_pos(1)-1;
%100
h=figure;
hold on
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [1 0 0], 'resolution', 0.1)
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_mag(i,:),'color',colors(i,:))
end
hold off
view(-90,0)
title('SQUID-MAG')
saveas(h, fullfile(save_path, 'figs', 'dipfit-SQUIDMAG_100.jpg'))
%010
h=figure;
hold on
tmp = pos_mag;
tmp(:,2) = mean_pos(2)-1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 1 0], 'resolution', 0.1)
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_mag(i,:),'color',colors(i,:))
end
hold off
view(0,0)
title('SQUID-MAG')
saveas(h, fullfile(save_path, 'figs', 'dipfit-SQUIDMAG_010.jpg'))
%001
h=figure;
hold on
tmp = pos_mag;
tmp(:,3) = mean_pos(3)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 0 1], 'resolution', 0.1)       
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_mag(i,:),'color',colors(i,:))
end
hold off
title('SQUID-MAG')
saveas(h, fullfile(save_path, 'figs', 'dipfit-SQUIDMAG_001.jpg'))
close all

% OPM
pos_opm = zeros(5,3);
ori_opm = zeros(5,3);
for i = 1:5
    pos_opm(i,:) = opm_dipole{i}.dip.pos;
    [~,idx] = max(vecnorm(opm_dipole{i}.dip.mom,2,1));
    ori_opm(i,:) = -opm_dipole{i}.dip.mom(:,idx);
end
mean_pos = mean(pos_opm,1);
tmp = pos_opm;
tmp(:,1) = mean_pos(1)-1;
%100
h=figure;
hold on
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [1 0 0], 'resolution', 0.1)
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_opm(i,:),'color',colors(i,:))
end
hold off
view(-90,0)
title('OPM')
saveas(h, fullfile(save_path, 'figs', 'dipfit-OPM_100.jpg'))
%010
h=figure;
hold on
tmp = pos_opm;
tmp(:,2) = mean_pos(2)-1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 1 0], 'resolution', 0.1)
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_opm(i,:),'color',colors(i,:))
end
hold off
view(0,0)
title('OPM')
saveas(h, fullfile(save_path, 'figs', 'dipfit-OPM_010.jpg'))
%001
h=figure;
hold on
tmp = pos_opm;
tmp(:,3) = mean_pos(3)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 0 1], 'resolution', 0.1)       
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_opm(i,:),'color',colors(i,:))
end
hold off
title('OPM')
saveas(h, fullfile(save_path, 'figs', 'dipfit-OPM_001.jpg'))
close all
%% Save 
save(fullfile(save_path, 'dipoles'), 'megmag_dipole', 'megplanar_dipole', 'opm_dipole', 'eeg_dipole', 'opmeeg_dipole'); disp('done');

end
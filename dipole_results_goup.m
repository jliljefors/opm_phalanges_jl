function dipole_results_goup(base_save_path, subs, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i_sub = subs
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);
    save_path = fullfile(base_save_path,params.sub);
    load(fullfile(save_path, 'dipoles')); 
    dipole_squidmag{i_sub} = megmag_dipole;
    dipole_squidgrad{i_sub} = megplanar_dipole;
    dipole_opm{i_sub} = opm_dipole;
    dipole_squideeg{i_sub} = eeg_dipole;
    dipole_opmeeg{i_sub} = opmeeg_dipole;

    n_ph = length(params.phalange_labels);
    % Metrics: 
    % - distance between dipoles for same phalange different systems
    % - over phalanges: average distance from mean location within distance
    pos_squidmag = zeros(n_ph,3);
    pos_squidgrad = zeros(n_ph,3);
    pos_opm = zeros(n_ph,3);
    pos_squideeg = zeros(n_ph,3);
    pos_opmeeg = zeros(n_ph,3);
    for i_phalange = 1:n_ph
        pos_squidmag(i_phalange,:) = dipole_squidmag{i_sub}.pos;
        pos_squidgrad(i_phalange,:) = dipole_squidgrad{i_sub}.pos;
        pos_opm(i_phalange,:) = dipole_opm{i_sub}.pos;
        pos_squideeg(i_phalange,:) = dipole_squideeg{i_sub}.pos;
        pos_opmeeg(i_phalange,:) = dipole_opmeeg{i_sub}.pos;

        dist_sqmag_opm(i_sub,i_phalange) = norm(pos_squidmag(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqgrad_opm(i_sub,i_phalange) = norm(pos_squidgrad(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqmag_sqgrad(i_sub,i_phalange) = norm(pos_squidmag(i_phalange,:)-pos_squidgrad(i_phalange,:));
        dist_sqeeg_opmeeg(i_sub,i_phalange) = norm(pos_squideeg(i_phalange,:)-pos_opmeeg(i_phalange,:));
    end
    spread_opm(i_sub,:) = vecnorm(pos_opm-repmat(mean(pos_opm,1),[1 3]),2,2)'; % mean distance from center of phalanges
    spread_squidmag(i_sub,:) = vecnorm(pos_squidmag-repmat(mean(pos_squidmag,1),[1 3]),2,2)'; % mean distance from center of phalanges
    spread_squidgrad(i_sub,:) = vecnorm(pos_squidgrad-repmat(mean(pos_squidgrad,1),[1 3]),2,2)'; % mean distance from center of phalanges
    spread_squideeg(i_sub,:) = vecnorm(pos_squideeg-repmat(mean(pos_squideeg,1),[1 3]),2,2)'; % mean distance from center of phalanges
    spread_opmeeg(i_sub,:) = vecnorm(pos_opmeeg-repmat(mean(pos_opmeeg,1),[1 3]),2,2)'; % mean distance from center of phalanges
end

%% Plot distances
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqmag_opm,1));
hold on
er = errorbar(1:5,mean(dist_sqmag_opm,1), mean(dist_sqmag_opm,1)-min(pdist_sqmag_opm,[],1), mean(dist_sqmag_opm,1)-max(dist_sqmag_opm,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Distance SQUID-MAG to OPM (mean = ' num2str(mean(mean(dist_sqmag_opm))) ')'])
ylabel('Distance [cm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'dipole_squidmag_to_opm_dist.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqgrad_opm,1));
hold on
er = errorbar(1:5,mean(dist_sqgrad_opm,1), mean(dist_sqgrad_opm,1)-min(pdist_sqgrad_opm,[],1), mean(dist_sqgard_opm,1)-max(dist_sqgrad_opm,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Distance SQUID-GRAD to OPM (mean = ' num2str(mean(mean(dist_sqgrad_opm))) ')'])
ylabel('Distance [cm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'dipole_squidgrad_to_opm_dist.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqmag_sqgrad,1));
hold on
er = errorbar(1:5,mean(dist_sqmag_sqgrad,1), mean(dist_sqmag_sqgrad,1)-min(dist_sqmag_sqgrad,[],1), mean(dist_sqmag_sqgrad,1)-max(dist_sqmag_sqgrad,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Distance SQUID-MAG to SQUID-GRAD (mean = ' num2str(mean(mean(dist_sqmag_sqgrad))) ')'])
ylabel('Distance [cm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'dipole_squidmag_to_squidgrad_dist.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqeeg_opmeeg,1));
hold on
er = errorbar(1:5,mean(dist_sqeeg_opmeeg,1), mean(dist_sqeeg_opmeeg,1)-min(dist_sqeeg_opmeeg,[],1), mean(dist_sqeeg_opmeeg,1)-max(dist_sqeeg_opmeeg,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Distance SQUID-EEG to OPM-EEG (mean = ' num2str(mean(mean(dist_sqeeg_opmeeg))) ')'])
ylabel('Distance [cm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'dipole_squideeg_to_opmeeg_dist.jpg'))

%% Plot spread
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(spread_opm,1));
hold on
er = errorbar(1:5,mean(spread_opm,1), mean(spread_opm,1)-min(spread_opm,[],1), mean(spread_opm,1)-max(spread_opm,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['OPM spread (mean = ' num2str(mean(mean(spread_opm))) ')'])
ylabel('Distance from center of phalanges [cm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'dipole_opm_spread.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(spread_squidmag,1));
hold on
er = errorbar(1:5,mean(spread_squidmag,1), mean(spread_squidmag,1)-min(spread_squidmag,[],1), mean(spread_squidmag,1)-max(spread_squidmag,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['SQUID-MAG spread (mean = ' num2str(mean(mean(spread_squidmag))) ')'])
ylabel('Distance from center of phalanges [cm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'dipole_squidmag_spread.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(spread_squidgrad,1));
hold on
er = errorbar(1:5,mean(spread_squidgrad,1), mean(spread_squidgrad,1)-min(spread_squidgrad,[],1), mean(spread_squidgrad,1)-max(spread_squidgrad,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['SQUID-GRAD spread (mean = ' num2str(mean(mean(spread_squidgrad))) ')'])
ylabel('Distance from center of phalanges [cm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'dipole_squidgrad_spread.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(spread_opmeeg,1));
hold on
er = errorbar(1:5,mean(spread_opmeeg,1), mean(spread_opmeeg,1)-min(spread_opmeeg,[],1), mean(spread_opmeeg,1)-max(spread_opmeeg,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['OPM-EEG spread (mean = ' num2str(mean(mean(spread_opmeeg))) ')'])
ylabel('Distance from center of phalanges [cm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'dipole_opmeeg_spread.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(spread_squideeg,1));
hold on
er = errorbar(1:5,mean(spread_squideeg,1), mean(spread_squideeg,1)-min(spread_squideeg,[],1), mean(spread_squideeg,1)-max(spread_squideeg,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['SQUID-EEG spread (mean = ' num2str(mean(mean(spread_squideeg))) ')'])
ylabel('Distance from center of phalanges [cm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'dipole_squideeg_spread.jpg'))
end
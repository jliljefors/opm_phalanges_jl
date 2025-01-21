function dipole_results_goup(base_save_path, subs, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n_subs = max(subs);
n_ph = length(params.phalange_labels);
dist_sqmag_opm = nan(n_subs,n_ph);
dist_sqgrad_opm = nan(n_subs,n_ph);
dist_sqmag_sqgrad = nan(n_subs,n_ph);
dist_sqeeg_opmeeg = nan(n_subs,n_ph);
spread_opm = nan(n_subs,n_ph);
spread_squidmag = nan(n_subs,n_ph);
spread_squidgrad = nan(n_subs,n_ph);
spread_squideeg = nan(n_subs,n_ph);
spread_opmeeg = nan(n_subs,n_ph);
for i_sub = subs
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);
    save_path = fullfile(base_save_path,params.sub);
    load(fullfile(save_path, 'dipoles')); 
    dipole_squidmag{i_sub} = megmag_dipole;
    dipole_squidgrad{i_sub} = megplanar_dipole;
    dipole_opm{i_sub} = opm_dipole;
    dipole_squideeg{i_sub} = megeeg_dipole;
    dipole_opmeeg{i_sub} = opmeeg_dipole;
  
    % Metrics: 
    % - distance between dipoles for same phalange different systems
    % - over phalanges: average distance from mean location within distance
    pos_squidmag = zeros(n_ph,3);
    pos_squidgrad = zeros(n_ph,3);
    pos_opm = zeros(n_ph,3);
    pos_squideeg = zeros(n_ph,3);
    pos_opmeeg = zeros(n_ph,3);
    for i_phalange = 1:n_ph
        pos_squidmag(i_phalange,:) = dipole_squidmag{i_sub}{i_phalange}.dip.pos;
        pos_squidgrad(i_phalange,:) = dipole_squidgrad{i_sub}{i_phalange}.dip.pos;
        pos_opm(i_phalange,:) = dipole_opm{i_sub}{i_phalange}.dip.pos;

        dist_sqmag_opm(i_sub,i_phalange) = 1e1*norm(pos_squidmag(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqgrad_opm(i_sub,i_phalange) = 1e1*norm(pos_squidgrad(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqmag_sqgrad(i_sub,i_phalange) = 1e1*norm(pos_squidmag(i_phalange,:)-pos_squidgrad(i_phalange,:));
        if ~isempty(dipole_squideeg{i_sub})
            pos_squideeg(i_phalange,:) = 1e1*dipole_squideeg{i_sub}{i_phalange}.dip.pos;
            pos_opmeeg(i_phalange,:) = 1e1*dipole_opmeeg{i_sub}{i_phalange}.dip.pos;
            dist_sqeeg_opmeeg(i_sub,i_phalange) = 1e1*norm(pos_squideeg(i_phalange,:)-pos_opmeeg(i_phalange,:));
        end
    end
    spread_opm(i_sub,:) = 1e1*vecnorm(pos_opm-repmat(mean(pos_opm,1),[n_ph 1]),2,2)'; % mean distance from center of phalanges
    spread_squidmag(i_sub,:) = 1e1*vecnorm(pos_squidmag-repmat(mean(pos_squidmag,1),[n_ph 1]),2,2)'; % mean distance from center of phalanges
    spread_squidgrad(i_sub,:) = 1e1*vecnorm(pos_squidgrad-repmat(mean(pos_squidgrad,1),[n_ph 1]),2,2)'; % mean distance from center of phalanges
    spread_squideeg(i_sub,:) = 1e1*vecnorm(pos_squideeg-repmat(mean(pos_squideeg,1),[n_ph 1]),2,2)'; % mean distance from center of phalanges
    spread_opmeeg(i_sub,:) = 1e1*vecnorm(pos_opmeeg-repmat(mean(pos_opmeeg,1),[n_ph 1]),2,2)'; % mean distance from center of phalanges
end

%% Plot distances
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqmag_opm,1,'omitnan'));
hold on
er = errorbar(1:5,mean(dist_sqmag_opm,1,'omitnan'), mean(dist_sqmag_opm,1,'omitnan')-min(dist_sqmag_opm,[],1,'omitnan'), mean(dist_sqmag_opm,1,'omitnan')-max(dist_sqmag_opm,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Dist SQ-MAG to OPM (mean = ' num2str(mean(mean(dist_sqmag_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'dipole_squidmag_to_opm_dist.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqgrad_opm,1,'omitnan'));
hold on
er = errorbar(1:5,mean(dist_sqgrad_opm,1,'omitnan'), mean(dist_sqgrad_opm,1,'omitnan')-min(dist_sqgrad_opm,[],1,'omitnan'), mean(dist_sqgrad_opm,1,'omitnan')-max(dist_sqgrad_opm,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Dist SQ-GRAD to OPM (mean = ' num2str(mean(mean(dist_sqgrad_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'dipole_squidgrad_to_opm_dist.jpg'))
close

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqmag_sqgrad,1,'omitnan'));
hold on
er = errorbar(1:5,mean(dist_sqmag_sqgrad,1,'omitnan'), mean(dist_sqmag_sqgrad,1,'omitnan')-min(dist_sqmag_sqgrad,[],1,'omitnan'), mean(dist_sqmag_sqgrad,1,'omitnan')-max(dist_sqmag_sqgrad,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Dist SQ-MAG to SQ-GRAD (mean = ' num2str(mean(mean(dist_sqmag_sqgrad,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'dipole_squidmag_to_squidgrad_dist.jpg'))
close

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqeeg_opmeeg,1,'omitnan'));
hold on
er = errorbar(1:5,mean(dist_sqeeg_opmeeg,1,'omitnan'), mean(dist_sqeeg_opmeeg,1,'omitnan')-min(dist_sqeeg_opmeeg,[],1,'omitnan'), mean(dist_sqeeg_opmeeg,1,'omitnan')-max(dist_sqeeg_opmeeg,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Dist SQ-EEG to OPM-EEG (mean = ' num2str(mean(mean(dist_sqeeg_opmeeg,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'dipole_squideeg_to_opmeeg_dist.jpg'))
close all

%% Plot spread
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(spread_opm,1,'omitnan'));
hold on
er = errorbar(1:5,mean(spread_opm,1,'omitnan'), mean(spread_opm,1,'omitnan')-min(spread_opm,[],1,'omitnan'), mean(spread_opm,1,'omitnan')-max(spread_opm,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['OPM spread (mean = ' num2str(mean(mean(spread_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Dipoles spread [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'dipole_opm_spread.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(spread_squidmag,1,'omitnan'));
hold on
er = errorbar(1:5,mean(spread_squidmag,1,'omitnan'), mean(spread_squidmag,1,'omitnan')-min(spread_squidmag,[],1,'omitnan'), mean(spread_squidmag,1,'omitnan')-max(spread_squidmag,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['SQ-MAG spread (mean = ' num2str(mean(mean(spread_squidmag,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Dipoles spread [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'dipole_squidmag_spread.jpg'))
close

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(spread_squidgrad,1,'omitnan'));
hold on
er = errorbar(1:5,mean(spread_squidgrad,1,'omitnan'), mean(spread_squidgrad,1,'omitnan')-min(spread_squidgrad,[],1,'omitnan'), mean(spread_squidgrad,1,'omitnan')-max(spread_squidgrad,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['SQ-GRAD spread (mean = ' num2str(mean(mean(spread_squidgrad,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Dipoles spread [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'dipole_squidgrad_spread.jpg'))
close

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(spread_opmeeg,1,'omitnan'));
hold on
er = errorbar(1:5,mean(spread_opmeeg,1,'omitnan'), mean(spread_opmeeg,1,'omitnan')-min(spread_opmeeg,[],1,'omitnan'), mean(spread_opmeeg,1,'omitnan')-max(spread_opmeeg,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['OPM-EEG spread (mean = ' num2str(mean(mean(spread_opmeeg,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Dipoles spread [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'dipole_opmeeg_spread.jpg'))
close

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(spread_squideeg,1,'omitnan'));
hold on
er = errorbar(1:5,mean(spread_squideeg,1,'omitnan'), mean(spread_squideeg,1,'omitnan')-min(spread_squideeg,[],1,'omitnan'), mean(spread_squideeg,1,'omitnan')-max(spread_squideeg,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['SQ-EEG spread (mean = ' num2str(mean(mean(spread_squideeg,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Dipoles spread [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'dipole_squideeg_spread.jpg'))
close all

end
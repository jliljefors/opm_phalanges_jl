function mne_results_goup(base_save_path, subs, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n_subs = max(subs);
n_ph = length(params.phalange_labels);
dist_sqmag_opm = nan(n_subs,n_ph);
dist_sqgrad_opm = nan(n_subs,n_ph);
dist_sqmag_sqgrad = nan(n_subs,n_ph);
fahm_opm = nan(n_subs,n_ph);
fahm_squidmag = nan(n_subs,n_ph);
fahm_squidgrad = nan(n_subs,n_ph);
for i_sub = subs
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);
    save_path = fullfile(base_save_path,params.sub);
    load(fullfile(save_path, 'squidmag_mne_M60')); 
    mne_squidmag{i_sub} = squidmag_mne_M60;
    load(fullfile(save_path, 'squidgrad_mne_M60'));
    mne_squidgrad{i_sub} = squidgrad_mne_M60;
    load(fullfile(save_path, 'opm_mne_M60'));
    mne_opm{i_sub} = opm_mne_M60;
  
    % Metrics: 
    % - distance between mnes for same phalange different systems
    % - over phalanges: average distance from mean location within distance
    pos_squidmag = zeros(n_ph,3);
    pos_squidgrad = zeros(n_ph,3);
    pos_opm = zeros(n_ph,3);
    for i_phalange = 1:n_ph
        pos_squidmag(i_phalange,:) = mne_squidmag{i_sub}{i_phalange}.peakloc;
        pos_squidgrad(i_phalange,:) = mne_squidgrad{i_sub}{i_phalange}.peakloc;
        pos_opm(i_phalange,:) = mne_opm{i_sub}{i_phalange}.peakloc;

        dist_sqmag_opm(i_sub,i_phalange) = 1e1*norm(pos_squidmag(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqgrad_opm(i_sub,i_phalange) = 1e1*norm(pos_squidgrad(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqmag_sqgrad(i_sub,i_phalange) = 1e1*norm(pos_squidmag(i_phalange,:)-pos_squidgrad(i_phalange,:));

        fahm_opm(i_sub,i_phalange) = mne_opm{i_sub}{i_phalange}.fahm; % mean distance from center of phalanges
        fahm_squidmag(i_sub,i_phalange) = mne_squidmag{i_sub}{i_phalange}.fahm; % mean distance from center of phalanges
        fahm_squidgrad(i_sub,i_phalange) = mne_squidgrad{i_sub}{i_phalange}.fahm; % mean distance from center of phalanges
    end
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
saveas(h, fullfile(base_save_path, 'figs', 'mne_squidmag_to_opm_dist.jpg'))

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
saveas(h, fullfile(base_save_path, 'figs', 'mne_squidgrad_to_opm_dist.jpg'))
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
saveas(h, fullfile(base_save_path, 'figs', 'mne_squidmag_to_squidgrad_dist.jpg'))
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
saveas(h, fullfile(base_save_path, 'figs', 'mne_squideeg_to_opmeeg_dist.jpg'))
close all

for i_ph = 1:5
    h = figure('DefaultAxesFontSize',16);
    plot(subs,dist_sqmag_opm(subs,i_ph));
    hold on
    plot(subs,dist_sqgrad_opm(subs,i_ph));
    plot(subs,dist_sqmag_sqgrad(subs,i_ph));
    hold off
    title([params.phalange_labels{i_ph} ' - mne distances over subjects'])
    ylabel('Distance [mm]')
    xlabel('Subjects')
    legend(['SQMAG-OPM   '; 'SQGRAD-OPM  '; 'SQMAG-SQGRAD'])
    saveas(h, fullfile(base_save_path, 'figs', ['mne_dist_vs_sub-' params.phalange_labels{i_ph} '.jpg']))
    close
end

%% Plot FAHM squidmag vs opm
data1 = fahm_squidmag;
data2 = fahm_opm;
data3 = fahm_squidgrad;
mean1 = mean(data1,1,'omitnan');
mean2 = mean(data2,1,'omitnan');
mean3 = mean(data3,1,'omitnan');
min1 = min(data1,[],1,'omitnan');
min2 = min(data2,[],1,'omitnan');
min3 = min(data3,[],1,'omitnan');
max1 = max(data1,[],1,'omitnan');
max2 = max(data2,[],1,'omitnan');
max3 = max(data3,[],1,'omitnan');
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];
err3 = [mean3-min3; max3-mean3];

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
    errorbar(k+0.22,mean3(k),err3(1,k),err3(2,k),'k','linestyle','none');
end

p_values = zeros(5, 3);
for i = 1:5
    [~, p_values(i, 1)] = ttest(data1(:,i), data2(:,i));
    [~, p_values(i, 2)] = ttest(data2(:,i), data3(:,i));
    [~, p_values(i, 3)] = ttest(data1(:,i), data3(:,i));
end
for i = 1:5
    sigstar({[i-0.22, i]}, p_values(i, 1));
    sigstar({[i, i+0.22]}, p_values(i, 2));
    sigstar({[i-0.22, i+0.22]}, p_values(i, 3));
end

hold off
title('MNE: Group level M100 FAHM')
ylabel('M60 FAHM [cm^2]')
xlabel('Phalange')
legend({'squidmag','opm','squidgrad'});
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'mne_fahm.jpg'))

%% Plot FAHM squidgrad vs opm
% data1 = fahm_squidgrad;
% data2 = fahm_opm;
% mean1 = mean(data1,1,'omitnan');
% mean2 = mean(data2,1,'omitnan');
% min1 = min(data1,[],1,'omitnan');
% min2 = min(data2,[],1,'omitnan');
% max1 = max(data1,[],1,'omitnan');
% max2 = max(data2,[],1,'omitnan');
% err1 = [mean1-min1; max1-mean1];
% err2 = [mean2-min2; max2-mean2];
% 
% h = figure('DefaultAxesFontSize',16);
% bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
% hold on
% for k=1:length(params.phalange_labels)
%     errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
%     errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
% end
% p_values = zeros(1, 5);
% for i = 1:5
%     [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
% end
% sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
% hold off
% title('MNE: Group level M100 FAHM')
% ylabel('M60 FAHM [cm^2]')
% xlabel('Phalange')
% legend({'squidgrad','opm'});
% xticklabels(params.phalange_labels)
% saveas(h, fullfile(base_save_path, 'figs', 'mne_fahm_squidgrad_opm.jpg'))

%% Plot FAHM squidgrad vs opm
% data1 = fahm_squidgrad;
% data2 = fahm_squidmag;
% mean1 = mean(data1,1,'omitnan');
% mean2 = mean(data2,1,'omitnan');
% min1 = min(data1,[],1,'omitnan');
% min2 = min(data2,[],1,'omitnan');
% max1 = max(data1,[],1,'omitnan');
% max2 = max(data2,[],1,'omitnan');
% err1 = [mean1-min1; max1-mean1];
% err2 = [mean2-min2; max2-mean2];
% 
% h = figure('DefaultAxesFontSize',16);
% bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
% hold on
% for k=1:length(params.phalange_labels)
%     errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
%     errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
% end
% p_values = zeros(1, 5);
% for i = 1:5
%     [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
% end
% sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
% hold off
% title('MNE: Group level M100 FAHM')
% ylabel('M60 FAHM [cm^2]')
% xlabel('Phalange')
% legend({'squidgrad','squidmag'});
% xticklabels(params.phalange_labels)
% saveas(h, fullfile(base_save_path, 'figs', 'mne_fahm_squidgrad_squidmag.jpg'))

close all

%% FAHM vs sub
for i_ph = 1:5
    h = figure('DefaultAxesFontSize',16);
    plot(subs,fahm_squidmag(subs,i_ph));
    hold on
    plot(subs,fahm_opm(subs,i_ph));
    plot(subs,fahm_squidgrad(subs,i_ph));
    hold off
    title([params.phalange_labels{i_ph} ' - mne FAHM over subjects'])
    ylabel('FAHM [cm^2]')
    xlabel('Subjects')
    legend(['SQMAG '; 'OPM   '; 'SQGRAD'])
    saveas(h, fullfile(base_save_path, 'figs', ['mne_fahm_vs_sub-' params.phalange_labels{i_ph} '.jpg']))
    close
end

end
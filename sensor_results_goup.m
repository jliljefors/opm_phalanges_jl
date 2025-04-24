function sensor_results_goup(base_save_path, subs, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

peak_ratio = [];
snr = [];
latency = [];
for i_sub = subs
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);
    save_path = fullfile(base_save_path,params.sub);
    clear M60_opm M60_opmeeg
    clear M60_squidmag M60_squidgrad M60_squideeg
    M60_opm = load(fullfile(save_path, [params.sub '_opm_M60.mat'])).M60; 
    M60_opmeeg = load(fullfile(save_path, [params.sub '_opmeeg_M60.mat'])).M60; 
    M60_squidmag = load(fullfile(save_path, [params.sub '_squidmag_M60.mat'])).M60; 
    M60_squidgrad = load(fullfile(save_path, [params.sub '_squidgrad_M60.mat'])).M60; 
    M60_squideeg = load(fullfile(save_path, [params.sub '_squideeg_M60.mat'])).M60;

    clear squidmag_timelocked opm_timelocked
    squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked'])).timelocked; 
    opm_timelocked = load(fullfile(save_path, [params.sub '_opm_timelocked'])).timelocked; 
    meg_chs = find(contains(squidmag_timelocked{1}.label,'MEG'));
    opm_chs = find(contains(opm_timelocked{1}.label,'bz'));

    for i_phalange = 1:length(params.phalange_labels)
        peak_ratio.meg(i_sub,i_phalange) = M60_opm{i_phalange}.peak_amplitude/M60_squidmag{i_phalange}.peak_amplitude;
        peak_ratio.eeg(i_sub,i_phalange) = M60_opmeeg{i_phalange}.peak_amplitude/M60_squideeg{i_phalange}.peak_amplitude;
        snr.error_opm(i_sub,i_phalange) = M60_opm{i_phalange}.peak_amplitude/M60_opm{i_phalange}.std_error;
        snr.error_squidmag(i_sub,i_phalange) = M60_squidmag{i_phalange}.peak_amplitude/M60_squidmag{i_phalange}.std_error;
        snr.error_squideeg(i_sub,i_phalange) = M60_squideeg{i_phalange}.peak_amplitude/M60_squideeg{i_phalange}.std_error;
        snr.error_opmeeg(i_sub,i_phalange) = M60_opmeeg{i_phalange}.peak_amplitude/M60_opmeeg{i_phalange}.std_error;
        snr.prestim_opm(i_sub,i_phalange) = M60_opm{i_phalange}.peak_amplitude/M60_opm{i_phalange}.prestim_std;
        snr.prestim_squidmag(i_sub,i_phalange) = M60_squidmag{i_phalange}.peak_amplitude/M60_squidmag{i_phalange}.prestim_std;
        snr.prestim_opmeeg(i_sub,i_phalange) = M60_opmeeg{i_phalange}.peak_amplitude/M60_opmeeg{i_phalange}.prestim_std;
        snr.prestim_squideeg(i_sub,i_phalange) = M60_squideeg{i_phalange}.peak_amplitude/M60_squideeg{i_phalange}.prestim_std;
        snr.ratio_error(i_sub,i_phalange) = snr.error_opm(i_sub,i_phalange)/snr.error_squidmag(i_sub,i_phalange);
        snr.ratio_prestim(i_sub,i_phalange) = snr.prestim_opm(i_sub,i_phalange)/snr.prestim_squidmag(i_sub,i_phalange);
        latency.opm(i_sub,i_phalange) = M60_opm{i_phalange}.peak_latency;
        latency.squidmag(i_sub,i_phalange) = M60_squidmag{i_phalange}.peak_latency;
        latency.squidgrad(i_sub,i_phalange) = M60_squidgrad{i_phalange}.peak_latency;
        latency.opmeeg(i_sub,i_phalange) = M60_opmeeg{i_phalange}.peak_latency;
        latency.squideeg(i_sub,i_phalange) = M60_squideeg{i_phalange}.peak_latency;
        amp.opm(i_sub,i_phalange) = M60_opm{i_phalange}.peak_amplitude;
        amp.squidmag(i_sub,i_phalange) = M60_squidmag{i_phalange}.peak_amplitude;
        amp.squidgrad(i_sub,i_phalange) = M60_squidgrad{i_phalange}.peak_amplitude;
        amp.opmeeg(i_sub,i_phalange) = M60_opmeeg{i_phalange}.peak_amplitude;
        amp.squideeg(i_sub,i_phalange) = M60_squideeg{i_phalange}.peak_amplitude;

        h = figure;
        subplot(2,1,1)
        plot(squidmag_timelocked{i_phalange}.time*1e3,squidmag_timelocked{i_phalange}.avg(meg_chs(1:3:end),:)*1e15)
        xlabel('t [msec]')
        ylabel('B [fT]')
        title(['Evoked SQUID MAG - phalange ' params.phalange_labels{i_phalange} ' (n_{trls}=' num2str(length(squidmag_timelocked{i_phalange}.cfg.trials)) ')'])
        subplot(2,1,2)
        plot(opm_timelocked{i_phalange}.time*1e3,opm_timelocked{i_phalange}.avg(opm_chs,:)*1e15)
        xlabel('t [msec]')
        ylabel('B [fT]')
        title(['Evoked OPM MAG - phalange ' params.phalange_labels{i_phalange} ' (n_{trls}=' num2str(length(opm_timelocked{i_phalange}.cfg.trials)) ')'])
        saveas(h, fullfile(save_path, 'figs', [params.sub '_squidopm_butterfly_ph-' params.phalange_labels{i_phalange} '.jpg']))
        close all
    end
end

%% Save
save(fullfile(base_save_path, 'group_sensor'),"peak_ratio","snr","latency","amp","-v7.3");

%% Plot SNR vs subs
for i_ph = 1:5
    h = figure('DefaultAxesFontSize',16);
    plot(subs,snr.error_opm(:,i_ph))
    hold on
    plot(subs,snr.error_squidmag(:,i_ph),':')
    hold off
    title("SNR_{error} over subjects (:=squidmag)")
    ylabel("SNR")
    xlabel("sub")
    saveas(h, fullfile(base_save_path, 'figs', ['Peak_SNR_vs_subs' params.phalange_labels{i_ph} '.jpg']))
end

%% Plot amp vs subs
for i_ph = 1:5
    h = figure('DefaultAxesFontSize',16);
    plot(subs,amp.opm(:,i_ph)*1e15)
    hold on
    plot(subs,amp.squidmag(:,i_ph)*1e15,':')
    hold off
    title("M60 amp over subjects (:=squidmag)")
    ylabel("fT")
    xlabel("sub")
    saveas(h, fullfile(base_save_path, 'figs', ['Peak_amp_vs_subs' params.phalange_labels{i_ph} '.jpg']))
end

%% Plot ratio
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(peak_ratio.meg,1,'omitnan'));
hold on
er = errorbar(1:5,mean(peak_ratio.meg,1,'omitnan'), mean(peak_ratio.meg,1,'omitnan')-min(peak_ratio.meg,[],1,'omitnan'), mean(peak_ratio.meg,1,'omitnan')-max(peak_ratio.meg,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 peak amp ratio (mean = ' num2str(mean(mean(peak_ratio.meg,'omitnan'),'omitnan'),'%.2f') ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'Peak_amplitude_ratios_meg.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(peak_ratio.eeg,1,'omitnan'));
hold on
er = errorbar(1:5,mean(peak_ratio.eeg,1,'omitnan'), mean(peak_ratio.eeg,1,'omitnan')-min(peak_ratio.eeg,[],1,'omitnan'), mean(peak_ratio.eeg,1,'omitnan')-max(peak_ratio.eeg,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 peak amp ratio (mean = ' num2str(mean(mean(peak_ratio.eeg,'omitnan'),'omitnan'),'%.2f') ')'])
ylabel('OPMEEG/SQUIDEEG')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'Peak_amplitude_ratios_eeg.jpg'))

%% Plot SNR - error
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(snr.ratio_error,1,'omitnan'));
hold on
er = errorbar(1:5,mean(snr.ratio_error,1,'omitnan'), mean(snr.ratio_error,1,'omitnan')-min(snr.ratio_error,[],1,'omitnan'), mean(snr.ratio_error,1,'omitnan')-max(snr.ratio_error,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 SNR_{stderror} ratio (mean = ' num2str(mean(mean(snr.ratio_error,'omitnan'),'omitnan'),'%.2f') ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'SNR_ratios_error.jpg'))

%% Plot SNR - prestim
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(snr.ratio_prestim,1,'omitnan'));
hold on
er = errorbar(1:5,mean(snr.ratio_prestim,1,'omitnan'), mean(snr.ratio_prestim,1,'omitnan')-min(snr.ratio_prestim,[],1,'omitnan'), mean(snr.ratio_prestim,1,'omitnan')-max(snr.ratio_prestim,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 SNR_{prestim} ratio (mean = ' num2str(mean(mean(snr.ratio_prestim,'omitnan'),'omitnan'),'%.2f') ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'SNR_ratios_prestim.jpg'))

%% Plot peak amp
% MEG
data1 = 1e15*amp.squidmag;
data2 = 1e15*amp.opm;
mean1 = mean(data1,1,'omitnan');
mean2 = mean(data2,1,'omitnan');
min1 = min(data1,[],1,'omitnan');
min2 = min(data2,[],1,'omitnan');
max1 = max(data1,[],1,'omitnan');
max2 = max(data2,[],1,'omitnan');
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
p_values = zeros(1, 5);
for i = 1:5
    [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
end
sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
hold off
title('Group level M100 amplitude')
ylabel('Peak amplitude [fT]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
legend({'squidmag','opm'});
saveas(h, fullfile(base_save_path, 'figs', 'Amplitude_meg.jpg'))

% EEG
data1 = 1e6*amp.squideeg;
data2 = 1e6*amp.opmeeg;
mean1 = mean(data1,1,'omitnan');
mean2 = mean(data2,1,'omitnan');
min1 = min(data1,[],1,'omitnan');
min2 = min(data2,[],1,'omitnan');
max1 = max(data1,[],1,'omitnan');
max2 = max(data2,[],1,'omitnan');
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
p_values = zeros(1, 5);
for i = 1:5
    [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
end
sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
hold off
title('Group level M100 amplitude')
ylabel('Peak amplitude [uV]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
legend({'squideeg','opmeeg'});
saveas(h, fullfile(base_save_path, 'figs', 'Amplitude_eeg.jpg'))

%% Plot peak latency
data1 = 1e3*latency.squidmag;
data2 = 1e3*latency.opm;
mean1 = mean(data1,1,'omitnan');
mean2 = mean(data2,1,'omitnan');
min1 = min(data1,[],1,'omitnan');
min2 = min(data2,[],1,'omitnan');
max1 = max(data1,[],1,'omitnan');
max2 = max(data2,[],1,'omitnan');
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
p_values = zeros(1, 5);
for i = 1:5
    [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
end
sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
hold off
title('Group level M100 latency')
ylabel('Latency [ms]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
legend({'squidmag','opm'},'Location','southeast');
saveas(h, fullfile(base_save_path, 'figs', 'Latency.jpg'))

%% Plot SNR - error
data1 = snr.error_squidmag;
data2 = snr.error_opm;
mean1 = mean(data1,1,'omitnan');
mean2 = mean(data2,1,'omitnan');
min1 = min(data1,[],1,'omitnan');
min2 = min(data2,[],1,'omitnan');
max1 = max(data1,[],1,'omitnan');
max2 = max(data2,[],1,'omitnan');
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
p_values = zeros(1, 5);
for i = 1:5
    [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
end
sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
hold off
title('Group level SNR_{m100,stderror}')
ylabel('SNR')
xlabel('Phalange')
legend({'squidmag','opm'});
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'SNR_error.jpg'))

%% Plot SNR - prestim
data1 = snr.prestim_squidmag;
data2 = snr.prestim_opm;
mean1 = mean(data1,1,'omitnan');
mean2 = mean(data2,1,'omitnan');
min1 = min(data1,[],1,'omitnan');
min2 = min(data2,[],1,'omitnan');
max1 = max(data1,[],1,'omitnan');
max2 = max(data2,[],1,'omitnan');
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
p_values = zeros(1, 5);
for i = 1:5
    [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
end
sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
hold off
title('Group level SNR_{m100,prestim}')
ylabel('SNR')
xlabel('Phalange')
legend({'squidmag','opm'});
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'SNR_prestim.jpg'))

close all
end
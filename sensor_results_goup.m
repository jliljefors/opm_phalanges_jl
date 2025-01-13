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
    load(fullfile(save_path, [params.sub '_opm_M100'])); 
    M100_opm{i_sub} = M100;
    load(fullfile(save_path, [params.sub '_opmeeg_M100'])); 
    M100_opmeeg{i_sub} = M100;
    load(fullfile(save_path, [params.sub '_megmag_M100'])); 
    M100_megmag{i_sub} = M100;
    load(fullfile(save_path, [params.sub '_megplanar_M100'])); 
    M100_megplanar{i_sub} = M100;
    load(fullfile(save_path, [params.sub '_megeeg_M100']));
    M100_megeeg{i_sub} = M100;
    
    load(fullfile(save_path, [params.sub '_meg_timelocked'])); 
    meg_timelocked = timelocked;
    load(fullfile(save_path, [params.sub '_opm_timelocked'])); 
    meg_chs = find(contains(meg_timelocked{1}.label,'MEG'));
    opm_timelocked = timelocked;
    opm_chs = find(contains(opm_timelocked{1}.label,'bz'));

    for i_phalange = 1:length(params.phalange_labels)
        peak_ratio.meg(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.max_amplitude/M100_megmag{i_sub}{i_phalange}.max_amplitude;
        peak_ratio.eeg(i_sub,i_phalange) = M100_opmeeg{i_sub}{i_phalange}.max_amplitude/M100_megeeg{i_sub}{i_phalange}.max_amplitude;
        snr.error_opm(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.max_amplitude/M100_opm{i_sub}{i_phalange}.std_error;
        snr.error_meg(i_sub,i_phalange) = M100_megmag{i_sub}{i_phalange}.max_amplitude/M100_megmag{i_sub}{i_phalange}.std_error;
        snr.error_megeeg(i_sub,i_phalange) = M100_megeeg{i_sub}{i_phalange}.max_amplitude/M100_megeeg{i_sub}{i_phalange}.std_error;
        snr.error_opmeeg(i_sub,i_phalange) = M100_opmeeg{i_sub}{i_phalange}.max_amplitude/M100_opmeeg{i_sub}{i_phalange}.std_error;
        snr.prestim_opm(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.max_amplitude/M100_opm{i_sub}{i_phalange}.prestim_std;
        snr.prestim_meg(i_sub,i_phalange) = M100_megmag{i_sub}{i_phalange}.max_amplitude/M100_megmag{i_sub}{i_phalange}.prestim_std;
        snr.prestim_opmeeg(i_sub,i_phalange) = M100_opmeeg{i_sub}{i_phalange}.max_amplitude/M100_opmeeg{i_sub}{i_phalange}.prestim_std;
        snr.prestim_megeeg(i_sub,i_phalange) = M100_megeeg{i_sub}{i_phalange}.max_amplitude/M100_megeeg{i_sub}{i_phalange}.prestim_std;
        snr.ratio_error(i_sub,i_phalange) = snr.error_opm(i_sub,i_phalange)/snr.error_meg(i_sub,i_phalange);
        snr.ratio_prestim(i_sub,i_phalange) = snr.prestim_opm(i_sub,i_phalange)/snr.prestim_meg(i_sub,i_phalange);
        latency.opm(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.peak_latency;
        latency.megmag(i_sub,i_phalange) = M100_megmag{i_sub}{i_phalange}.peak_latency;
        latency.megplanar(i_sub,i_phalange) = M100_megplanar{i_sub}{i_phalange}.peak_latency;
        latency.opmeeg(i_sub,i_phalange) = M100_opmeeg{i_sub}{i_phalange}.peak_latency;
        latency.megeeg(i_sub,i_phalange) = M100_megeeg{i_sub}{i_phalange}.peak_latency;

        h = figure;
        subplot(2,1,1)
        plot(meg_timelocked{i_phalange}.time*1e3,meg_timelocked{i_phalange}.avg(meg_chs(1:3:end),:)*1e15)
        xlabel('t [msec]')
        ylabel('B [fT]')
        title(['Evoked SQUID MAG - phalange ' params.phalange_labels{i_phalange} ' (n_{trls}=' num2str(length(meg_timelocked{i_phalange}.cfg.trials)) ')'])
        subplot(2,1,2)
        plot(opm_timelocked{i_phalange}.time*1e3,meg_timelocked{i_phalange}.avg(opm_chs,:)*1e15)
        xlabel('t [msec]')
        ylabel('B [fT]')
        title(['Evoked OPM MAG - phalange ' params.phalange_labels{i_phalange} ' (n_{trls}=' num2str(length(opm_timelocked{i_phalange}.cfg.trials)) ')'])
        saveas(h, fullfile(save_path, 'figs', [params.sub '_megopm_butterfly_ph-' params.phalange_labels{i_phalange} '.jpg']))
    end
    close all
end

%% Save
save_path = base_save_path;
save(fullfile(save_path, 'group_sensor'),"peak_ratio","snr","latency","-v7.3");

%% Plot ratio
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(peak_ratio.meg,1));
hold on
er = errorbar(1:5,mean(peak_ratio.meg,1), mean(peak_ratio.meg,1)-min(peak_ratio.meg,[],1), mean(peak_ratio.meg,1)-max(peak_ratio.meg,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 peak amp ratio (mean = ' num2str(mean(mean(peak_ratio.meg)),'%.2f') ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'Peak_amplitude_ratios_meg.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(peak_ratio.eeg,1));
hold on
er = errorbar(1:5,mean(peak_ratio.eeg,1), mean(peak_ratio.eeg,1)-min(peak_ratio.eeg,[],1), mean(peak_ratio.eeg,1)-max(peak_ratio.eeg,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 peak amp ratio (mean = ' num2str(mean(mean(peak_ratio.eeg)),'%.2f') ')'])
ylabel('OPMEEG/SQUIDEEG')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'Peak_amplitude_ratios_eeg.jpg'))

%% Plot SNR - error
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(snr.ratio_error,1));
hold on
er = errorbar(1:5,mean(snr.ratio_error,1), mean(snr.ratio_error,1)-min(snr.ratio_error,[],1), mean(snr.ratio_error,1)-max(snr.ratio_error,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 SNR_{stderror} ratio (mean = ' num2str(mean(mean(snr.ratio_error)),'%.2f') ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'SNR_ratios_error.jpg'))

%% Plot SNR - prestim
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(snr.ratio_prestim,1));
hold on
er = errorbar(1:5,mean(snr.ratio_prestim,1), mean(snr.ratio_prestim,1)-min(snr.ratio_prestim,[],1), mean(snr.ratio_prestim,1)-max(snr.ratio_prestim,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 SNR_{prestim} ratio (mean = ' num2str(mean(mean(snr.ratio_prestim)),'%.2f') ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'SNR_ratios_prestim.jpg'))

%% Plot peak latency
data1 = 1e3*latency.megmag;
data2 = 1e3*latency.opm;
mean1 = mean(data1,1);
mean2 = mean(data2,1);
min1 = min(data1,[],1);
min2 = min(data2,[],1);
max1 = max(data1,[],1);
max2 = max(data2,[],1);
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
hold off
title('Group level M100 latency')
ylabel('Latency [ms]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
legend({'megmag','opm'});
saveas(h, fullfile(save_path, 'figs', 'Latency.jpg'))

%% Plot SNR - error
data1 = snr.error_opm;
data2 = snr_error_meg;
mean1 = mean(data1,1);
mean2 = mean(data2,1);
min1 = mins(data1,[],1);
min2 = min(data2,[],1);
max1 = max(data1,[],1);
max2 = max(data2,[],1);
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
hold off
title('Group level SNR_{m100,stderror}')
ylabel('SNR')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'SNR_error.jpg'))

%% Plot SNR - prestim
data1 = snr.prestim_opm;
data2 = snr_prestim_meg;
mean1 = mean(data1,1);
mean2 = mean(data2,1);
min1 = mins(data1,[],1);
min2 = min(data2,[],1);
max1 = max(data1,[],1);
max2 = max(data2,[],1);
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
hold off
title('Group level SNR_{m100,prestim}')
ylabel('SNR')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'SNR_prestim.jpg'))

end
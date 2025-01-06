function [timelocked] = timelock_MEG(data, save_path, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
timelocked = cell(length(params.trigger_code),1);
M100 = cell(5,1);
h = figure; 
hold on
leg = [];

cfg = [];
cfg.channel = params.chs;
data = ft_selectdata(cfg, data);
for i_phalange = 1:length(params.trigger_code)
    cfg = [];
    cfg.covariance          = 'yes';
    cfg.covariancewindow    = 'prestim';
    cfg.trials = find(data.trialinfo==params.trigger_code(i_phalange));
    timelocked{i_phalange} = ft_timelockanalysis(cfg, data);

    [~, interval_M100(1)] = min(abs(dat.time-0.08));
    [~, interval_M100(2)] = min(abs(dat.time-0.125));
    [~, interval_M100(3)] = min(abs(dat.time-0));
    tmp = [];
    [tmp.max_amplitude, i_maxch] = max(max(dat.avg(:,(interval_M100(1):interval_M100(2))),[],2));
    tmp.max_channel = dat.label{i_maxch};
    [~,i_peak_latency] = max(dat.avg(i_maxch,interval_M100(1):interval_M100(2)));
    tmp.peak_latency = dat.time(interval_M100(1)-1+i_peak_latency);
    tmp.prestim_std = std(dat.avg(i_maxch,1:interval_M100(3)));
    tmp.std_error = sqrt(dat.var(i_maxch,tmp.peak_latency));
    %for i_trl = find(data.trialinfo==params.trigger_code(i_phalange))'
    %    tmp.std_error = tmp.std_error + abs(tmp.max_amplitude - data.trial{i_trl}(i_maxch,interval_M100(1)-1+i_peak_latency));
    %end
    %tmp.std_error = tmp.std_error/length(find(data.trialinfo==params.trigger_code(i_phalange)));
    M100{i_phalange} = tmp;
    
    plot(dat.time.*1e3, dat.avg(i_maxch,:).*1e15)
    leg = [leg; [num2str(i_phalange) ': ' strrep(tmp.max_channel,'_','-')]];

end
hold off
title([params.modality ' - Max channel'])
ylabel('Field [fT]')
xlabel('time [ms]')
legend(leg)
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_evoked_maxchannels.jpg']))

save(fullfile(save_path, [params.sub '_' params.modality '_timelocked']), 'timelocked', '-v7.3'); 
save(fullfile(save_path, [params.sub '_' params.modality '_M100']), 'M100', '-v7.3'); 

%% Plot max channel with variation and peak time
for i_phalange = 1:length(params.trigger_code)
    h = figure;
    hold on
    for i_trl = find(data.trialinfo==params.trigger_code(i_phalange))
        plot(data.time{i_trl}.*1e3, data.trial{i_trl}(i_maxch,:).*1e15,[211,211,211]/255)
    end
    plot(dat.time.*1e3, dat.avg(i_maxch,:).*1e15,[0,0,0]/255)
    ylimits = ylim;
    latency = 1e3*M100{i_phalange}.peakl_latency;
    plot([latency latency],ylimits,'r--')
    hold off
    title(['Peak channel ' params.modality ': ' M100{i_phalange}.max_channel])
    ylabel('Field [fT]')
    xlabel('time [ms]')
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_evoked_maxchannel_ph' num2str(i_phalange) '.jpg']))
end

%% Topoplots
for i_phalange = 1:length(params.trigger_code)
    chs = find(contains(timelocked{i_phalange}.label,ft_channelselection(params.chs,timelocked{i_phalange}.label)));
    h = figure;
    plot(timelocked{i_phalange}.time*1e3,timelocked{i_phalange}.avg(chs,:)*1e15)
    xlabel('t [msec]')
    ylabel('B [fT]')
    title(['Evoked ' params.modality ' - phalange ' params.phalange_labels{i_phalange} ' (n_trls=' num2str(length(timelocked{i_phalange}.cfg.trials)) ')'])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_butterfly_ph-' params.phalange_labels{i_phalange} '.jpg']))

    cfg = [];
    cfg.xlim = [M100{i_phalange}.peak_latency-0.01 M100{i_phalange}.peak_latency+0.01];
    cfg.layout = params.layout; 
    cfg.parameter = 'avg';
    h = figure;
    ft_topoplotER(cfg, timelocked{i_phalange});
    axis on
    colorbar
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_M100_topo_ph-' params.phalange_labels{i_phalange} '.jpg']))
end
end
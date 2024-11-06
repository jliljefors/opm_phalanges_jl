function [] = plot_sub_evoked(save_path, i_sub, phalange_labels,opm_timelocked, meg_timelocked)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
%% OPM multiplot
cfg = [];
cfg.layout = 'fieldlinebeta2bz_helmet.mat'; 
cfg.parameter = 'avg';
cfg.channel = '*bz';
h = figure;
ft_multiplotER(cfg, opm_timelocked{1},opm_timelocked{2},opm_timelocked{3},opm_timelocked{4},opm_timelocked{5});
savefig(h, fullfile(save_path, 'figs', 'opm_evoked_multiplot.fig'))

%% MEG multiplot
cfg = [];
cfg.layout = 'neuromag306mag.lay'; 
cfg.parameter = 'avg';
cfg.channel = 'megmag';
h = figure;
ft_multiplotER(cfg, meg_timelocked{1},meg_timelocked{2},meg_timelocked{3},meg_timelocked{4},meg_timelocked{5});
savefig(h, fullfile(save_path, 'figs', 'meg_evoked_multiplot.fig'))

%% MEG&OPM butterfly (separated phalanges)
h = figure('Position',[89 89 912*2 560]); 
for i_phalange = 1:5
    cfg = [];
    cfg.channel = 'megmag';
    data = ft_selectdata(cfg, meg_timelocked{i_phalange});
    subplot(2,5,i_phalange)
    plot(data.time.*1e3, data.avg(:,:).*1e15)
    title(phalange_labels(i_phalange));
    ylabel('Field [fT]')
    xlabel('time [ms]')
    xlim([-30 300])
    ylim([-300 300])

    cfg = [];
    cfg.channel = '*bz';
    data = ft_selectdata(cfg, opm_timelocked{i_phalange});
    subplot(2,5,5+i_phalange)
    plot(data.time.*1e3, data.avg(:,:).*1e15)
    title(phalange_labels(i_phalange));
    ylabel('Field [fT]')
    xlabel('time [ms]')
    xlim([-30 300])
    ylim([-600 600])
end
savefig(h, fullfile(save_path, 'figs', 'evoked_butterfly_sub.fig'))

%% Max channels
leg = [];
M100 = cell(5,1);
h = figure; 
subplot(2,1,1)
hold on
for i_phalange = 1:5
    M100{i_phalange} = [];
    % meg
    cfg = [];
    cfg.channel = 'megmag';
    data = ft_selectdata(cfg, meg_timelocked{i_phalange});
    [~, interval_M100(1)] = min(abs(data.time-0.08));
    [~, interval_M100(2)] = min(abs(data.time-0.125));
    tmp = [];
    [tmp.maxamp, tmp.i_maxch] = max(max(data.avg(:,(interval_M100(1):interval_M100(2))),[],2));
    tmp.label = data.label{tmp.i_maxch};
    [~,tmp.i_tmax] = max(data.avg(tmp.i_maxch,interval_M100(1):interval_M100(2)));
    tmp.tmax = data.time(interval_M100(1)-1+tmp.i_tmax);
    t_max_meg(i_sub,i_phalange) = tmp.tmax;
    M100{i_phalange}.meg = tmp;
    plot(data.time.*1e3, data.avg(tmp.i_maxch,:).*1e15)
    leg = [leg; [num2str(i_phalange) ': ' tmp.label]];
end
hold off
title('SQUID - max channel')
ylabel('Field [fT]')
xlabel('time [ms]')
legend(leg)

leg = [];
subplot(2,1,2)
hold on
for i_phalange = 1:5
    % opm
    cfg = [];
    cfg.channel = '*bz';
    data = ft_selectdata(cfg, opm_timelocked{i_phalange});
    [~, interval_M100(1)] = min(abs(data.time-0.08));
    [~, interval_M100(2)] = min(abs(data.time-0.125));
    tmp = [];
    [tmp.maxamp, tmp.i_maxch] = max(max(data.avg(:,interval_M100(1):interval_M100(2)),[],2));
    tmp.label = data.label{tmp.i_maxch};
    [~,tmp.i_tmax] = max(data.avg(tmp.i_maxch,interval_M100(1):interval_M100(2)));
    tmp.tmax = data.time(interval_M100(1)-1+tmp.i_tmax);
    t_max_opm(i_sub,i_phalange) = tmp.tmax;
    M100{i_phalange}.opm = tmp;
    plot(data.time.*1e3, data.avg(tmp.i_maxch,:).*1e15)
    leg = [leg; [num2str(i_phalange) ': ' tmp.label]];

    M100{i_phalange}.peak_ratio = M100{i_phalange}.opm.maxamp/M100{i_phalange}.meg.maxamp;
    peak_ratio(i_sub,i_phalange) = M100{i_phalange}.peak_ratio;
end
hold off
title('OPM - max channel')
ylabel('Field [fT]')
xlabel('time [ms]')
legend(leg)
savefig(h, fullfile(save_path, 'figs', 'evoked_maxchannel.fig'))

%% Topoplots
for i_phalange = 1:5
    cfg = [];
    cfg.xlim = [t_max_opm(i_sub,i_phalange)-0.01 t_max_opm(i_sub,i_phalange)+0.01];
    cfg.layout = 'fieldlinebeta2bz_helmet.mat'; 
    cfg.parameter = 'avg';
    h = figure;
    ft_topoplotER(cfg, opm_timelocked{i_phalange});
    axis on
    colorbar
    savefig(h, fullfile(save_path, 'figs', ['M100_opm_topo_ph-' phalange_labels{i_phalange} '.fig']))

    cfg = [];
    cfg.xlim = [t_max_meg(i_sub,i_phalange)-0.01 t_max_meg(i_sub,i_phalange)+0.01];
    cfg.layout = 'neuromag306mag.lay'; 
    cfg.parameter = 'avg';
    h = figure;
    ft_topoplotER(cfg, meg_timelocked{i_phalange});
    axis on
    colorbar
    savefig(h, fullfile(save_path, 'figs', ['M100_meg_topo_ph-' phalange_labels{i_phalange} '.fig']))
end
end
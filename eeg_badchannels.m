function [badchs, badchs_flat, badchs_neighbors] = eeg_badchannels(cfg, data)
%opm_badchannels Detects channels that are flat, have low correlation with
%their neighbors or show a lot of jumping artifacts.
%   cfg.z_threshold
%   cfg.corr_threshold
%   cfg.n_neighbors
%   cfg.njump_threshold

n_neighbors     = ft_getopt(cfg, 'n_neighbors', 4);
corr_threshold  = ft_getopt(cfg, 'corr_threshold', 0.6);

cfg = [];
cfg.channel = 'EEG';
data = ft_selectdata(cfg, data);

chs = find(contains(data.label,'EEG'));
n_chs = length(chs);

%% Find channels with flat segments or high std
cfg = [];
cfg.length = 1;
data_seg = ft_redefinetrial(cfg,data);
n_trls = length(data_seg.trial);
trl_std = zeros(n_chs,n_trls);
for i_trl = 1:n_trls
    trl_std(:,i_trl) = std(data_seg.trial{i_trl},0,2);
end
badchs_flat = find(any(trl_std<1e-9,2));

%% Neighbors
goodchs = setdiff(chs,badchs_flat);

cfg = [];
cfg.resamplefs = 200;
cfg.lpfilter = 'yes';
cfg.lpfreq = 30;
cfg.lpinstabilityfix  = 'reduce';
data_lp = ft_resampledata(cfg,data);

cfg = [];
cfg.length = 1;
data_lp = ft_redefinetrial(cfg,data_lp);

% Create neighbor structure
n_chs = length(data_lp.label);
chanpos = data_lp.elec.chanpos;
neighbors = zeros(n_chs,n_neighbors);
for i = 1:size(chanpos,1)
        [~,tmp]= sort(vecnorm(chanpos(goodchs,:)-repmat(chanpos(i,:),[length(goodchs) 1]),2,2));
        neighbors(i,:) = goodchs(tmp(2:(n_neighbors+1)));
end
neighborscorr = zeros(n_chs,n_neighbors,length(data_lp.trial));
for i_trl = 1:length(data_lp.trial)
    dat = data_lp.trial{i_trl};
    for i = 1:n_chs
        for j = 1:n_neighbors
                tmp2 = corrcoef(dat(i,:),dat(int32(neighbors(i,j)),:));
                neighborscorr(i,j,i_trl) = abs(tmp2(1,2));
        end
    end 
end
badchs_neighbors = find(max(mean(neighborscorr,3),[],2)<corr_threshold,2); % bad if no neighbors exceed correlation threshold

badchs = [badchs_flat; badchs_neighbors];

% Convert to channel labels
badchs = data.label(chs(badchs));
badchs_flat = data.label(chs(badchs_flat));
badchs_neighbors = data.label(chs(badchs_neighbors));
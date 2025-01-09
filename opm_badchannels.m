function [badchs, badchs_flat, badchs_neighbors, badchs_zmax, badtrls_zmax] = opm_badchannels(cfg, data)
%opm_badchannels Detects channels that are flat, have low correlation with
%their neighbors or show a lot of jumping artifacts.
%   cfg.z_threshold
%   cfg.corr_threshold
%   cfg.n_neighbors
%   cfg.njump_threshold

cfg.n_neighbors     = ft_getopt(cfg, 'channel', 4);
cfg.corr_threshold  = ft_getopt(cfg, 'channel', 0.6);
cfg.z_threshold     = ft_getopt(cfg, 'channel', 20);
cfg.njump_threshold = ft_getopt(cfg, 'channel', 0.1);

chs = find(contains(data.label,'_bz'));
n_chs = length(chs);
n_trls = length(data.trial);

% Create neighbor structure
neighbors = zeros(n_chs,cfg.n_neighbors);
for i = 1:size(data.grad.chanpos,1)
        [~,tmp]= sort(vecnorm(data.grad.chanpos-repmat(data.grad.chanpos(i,:),[n_chs 1]),2,2));
        neighbors(i,:) = tmp(2:(cfg.n_neighbors+1));
end

neighborscorr = zeros(n_chs,cfg.n_neighbors,n_trls);
trial_std = zeros(n_chs,length(data.trial));
trial_mean = zeros(n_chs,length(data.trial));
z_max = zeros(n_chs,length(data.trial));
z_mean = zeros(n_chs,length(data.trial));
for trial = 1:n_trls
    dat = data.trial{trial}(chs,:);
    for i = 1:n_chs
        for j = 1:cfg.n_neighbors
                tmp2 = corrcoef(dat(i,:),dat(int32(neighbors(i,j)),:));
                neighborscorr(i,j,trial) = tmp2(1,2);
        end
    end
    dat = diff(movmedian(dat,9*5,2),1,2);
    trial_std(:,trial) = std(dat,0,2);
    trial_mean = repmat(mean(dat,2),[1 size(dat,2)]);
    z_max(:,trial) = max(abs(dat-trial_mean),[],2);
end    
z_max = z_max./repmat(mean(trial_std,2),[1 n_trls]);

% Bad channels (flat, low neighbor correlation, large number of jumps)
badchs_flat = find(any(trial_std<1e-15,2));
badchs_neighbors = find(~any(trimmean(neighborscorr,0.1,3)>cfg.corr_threshold,2)); % bad if no neighbors exceed correlation threshold
badchs_zmax = find(sum(z_max(setdiff(1:end,[badchs_flat; badchs_neighbors]),:)>cfg.z_threshold,2)>(n_trls*cfg.njump_threshold));
badchs = [badchs_flat; badchs_neighbors; badchs_zmax];

% Bad trials (jumps)
badtrls_zmax = find(sum(z_max(setdiff(1:end,badchs),:)>cfg.z_threshold,1)>0);

% Convert to channel labels
badchs = data.label(chs(badchs));
badchs_flat = data.label(chs(badchs_flat));
badchs_neighbors = data.label(chs(badchs_neighbors));
badchs_zmax = data.label(chs(badchs_zmax));
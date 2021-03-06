% verification for scattering events
%%
plotON = 0;
nPhoton = [];
i = 3; % which wavelength
[traj_group, max_depth] = XINTRINSIC_SIM_PlotTraj(S.traj2{i}, cfg.bounds, [1 300 0 60], nPhoton, plotON);

%% plot visitation depth
figure, 
hist(max_depth, 300)
length(find(max_depth > cfg.bounds(1)))/length(max_depth)
%% all photons 
CosTheta = cell(1, length(traj_group));
Dist = CosTheta;
for iPhoton = 1:length(traj_group)
    for iEvent = 2:length(traj_group{iPhoton}) - 1
        u = traj_group{iPhoton}(iEvent,:) - traj_group{iPhoton}(iEvent - 1,:);
        v = traj_group{iPhoton}(iEvent+1,:) - traj_group{iPhoton}(iEvent,:);

%         CosTheta{iPhoton}(iEvent) = (dot(u,v) / (norm(u)*norm(v)));
        CosTheta{iPhoton}(iEvent) = dot(u / norm(u), v / norm(v));
        Dist{iPhoton}(iEvent) = norm(v).*cfg.unitinmm;
    end
    
end
%% only those traveled in skull
ind = find(max_depth <= cfg.bounds(1));
CosTheta = cell(1, length(ind));
Dist = CosTheta;
for iPhoton = 1:length(ind)
    for iEvent = 2:length(traj_group{ind(iPhoton)}) - 1
        u = traj_group{ind(iPhoton)}(iEvent,:) - traj_group{ind(iPhoton)}(iEvent - 1,:);
        v = traj_group{ind(iPhoton)}(iEvent+1,:) - traj_group{ind(iPhoton)}(iEvent,:);

%         CosTheta{iPhoton}(iEvent) = (dot(u,v) / (norm(u)*norm(v)));
        CosTheta{iPhoton}(iEvent - 1) = dot(u / norm(u), v / norm(v));
        Dist{iPhoton}(iEvent - 1) = norm(v).*cfg.unitinmm;
    end
end
CosTheta_mat = cell2mat(CosTheta);
CosTheta_mat(isnan(CosTheta_mat)) = [];
Dist_mat = cell2mat(Dist);

figure,
subplot(1,2,1), hist(CosTheta_mat, 360), xlabel('scattering angles (cosine)')
title(['mean = ', num2str(mean(CosTheta_mat))])
subplot(1,2,2), hist(Dist_mat, 360), xlabel('distance between events (mm)')
title(['mean = ', num2str( mean(Dist_mat))])




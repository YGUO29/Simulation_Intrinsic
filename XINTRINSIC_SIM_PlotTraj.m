% plot trajectories 
% input traj is the output structure "traj" from MCX
% AxLim = [1 300 0 60]
% [traj_group, max_depth] = XINTRINSIC_SIM_PlotTraj(S.traj2{1}, cfg.bounds, [1 300 0 60], 10, 1)
% output: traj_group -- trajectories grouped by each detected photon
% output: max_depth -- visition depth for each photon
function [traj_group, max_depth] = XINTRINSIC_SIM_PlotTraj(traj, bounds1, AxLim, nPhoton, plotON)

if plotON 
    figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
    hold on
end
% title([' photon trajectories, no skull \newline wavelength = ', num2str(WLs(i)),', #DetPhoton = ',num2str(length(photon_id)) ])
photon_id = unique(traj.id);
if isempty(nPhoton)
    nPhoton = length(photon_id);
end
traj_group = cell(1,nPhoton);
max_depth = zeros(1,nPhoton);
fwait = waitbar(0,'Progress...');

% for k = 1:length(photon_id)
for k = 1:nPhoton
    waitbar(k/nPhoton,fwait,['Photon ID: ',num2str(k),'/',num2str(nPhoton)]);

    ind_traj = find(traj.id == photon_id(k));
%     plot3(traj.pos(ind_traj,1), traj.pos(ind_traj,2), traj.pos(ind_traj,3))
%     view(-37.5, 30)
%     axis([0 300 0 300 0 300])
%     grid on
    if plotON
        scatter(traj.pos(ind_traj(1),1), traj.pos(ind_traj(1),3),50); % the start point
        plot(traj.pos(ind_traj,1), traj.pos(ind_traj,3)); % plot on x-z plane
        axis(AxLim)
        for b = 1:length(bounds1)
            plot(AxLim(1:2),[bounds1(b) bounds1(b)],'--');
        end
        pause
        cla
    end
    max_depth(k) = max(traj.pos(ind_traj,3));

%     xticklabels({'5' '6' '7' '8' '9' '10'})
%     yticks([0 20 40 60]), yticklabels({'0' '1' '2' '3'})
%     xlabel('lateral distance, mm')
%     ylabel('depth, mm')
    traj_group{k} = traj.pos(ind_traj, :);

end

end
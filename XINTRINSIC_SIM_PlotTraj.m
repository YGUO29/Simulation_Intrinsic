% plot trajectories 
% input traj is the output structure "traj" from MCX
% AxLim = [1 300 0 60]
function traj_group = XINTRINSIC_SIM_PlotTraj(traj, bounds1, AxLim, nPhoton, plotON)

if plotON 
    figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
    hold on
end
photon_id = unique(traj.id);

% title([' photon trajectories, no skull \newline wavelength = ', num2str(WLs(i)),', #DetPhoton = ',num2str(length(photon_id)) ])
traj_group = cell(1,length(photon_id));
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

%     xticklabels({'5' '6' '7' '8' '9' '10'})
%     yticks([0 20 40 60]), yticklabels({'0' '1' '2' '3'})
%     xlabel('lateral distance, mm')
%     ylabel('depth, mm')
    traj_group{k} = traj.pos(ind_traj, :);

end

end
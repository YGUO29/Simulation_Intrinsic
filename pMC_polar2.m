function [profile, x, S]  = pMC_polar2(cfg, det1, seeds, method, para)
% this version track individual detected photons, rotate it according to
% the scattering angle, and preserve the polarization angle for each step 
%%
% newcfg              = cfg;
% newcfg.seed         = seeds{i}.data(:,det1{i}.DetInd);
% newcfg.detphotons   = det1{i}.data(:,det1{i}.DetInd);
% 
% tic, [~, det_temp1, ~, ~,traj1] = mcxlab(newcfg); toc
% nPhoton = 100;
% traj_group = XINTRINSIC_SIM_PlotTraj(traj1, cfg.bounds, [1 400 1 500], nPhoton, 0);
%%
s0 = [1,0];
for j = 1
    p = traj_group{j}; % positions
    for k = 1:length(p)
        figure, hold on
        if k == 1
            v1 = [p(k,1:2),-1] - p(k,:); v1 = v1./vecnorm(v1);
        else
            v1 = p(k,:) - p(k-1,:); v1 = v1./vecnorm(v1);
        end
        v = v1;
        plot3([0, v(1)],[0, v(2)],[0, v(3)])
        
        v2 = p(k+1,:) - p(k,:); v2 = v2./vecnorm(v2);
        v = v2;
        plot3([0, v(1)],[0, v(2)],[0, v(3)])
        
        norm = cross(v1,v2); % d = 0, because this plane passes [0 0 0]
        v = norm;
        plot3([0, v(1)],[0, v(2)],[0, v(3)])
%         view(-37.5,30)
       figure,
        Xlim = get(gca,'xlim'); Ylim = get(gca,'xlim');
        [X,Y] = meshgrid(Xlim(1):0.0001:Xlim(2), Ylim(1):0.01:Ylim(2));
        Z = (- (norm(1)*X) - (norm(2)*Y))/norm(3);
        mesh(X,Y,Z)
        grid on
    end
end
    
%%
S.f2 = cell(1,length(para.ind)); 
S.det2 = S.f2;
S.traj2 = S.f2;
% S.NewDetWeight_depth = S.f2;
S.NewDetWeight = S.f2;
S.DetWeight = zeros(1,6);
x = 0;
profile = S.f2;

SimVol = para.SimVol; % "replay" simulation volume, does not have to be the entire volume (to save time)
SimCenter = cfg.size(1)/2 - 1; % the illumination center is cfg.size(1)/2-1

%% get initial total photon weights
for i = 2
% change properties
% switch cfg.config
%     case {'without skull', 'without skull off focus'}
%         if isnan(cfg.perturb)
%         cfg.prop = [0 0 1 1            % medium 0: the environment
%         para.brain_mua(para.ind(i)) para.gray_matter_mus(para.ind(i)) para.gray_matter_g(para.ind(i)), 1.37   % medium 1: gray matter
%         para.brain_mua(para.ind(i)) para.white_matter_mus(para.ind(i)) para.white_matter_g(para.ind(i)), 1.37   % medium 2: white matter
%         para.brain_mua(para.ind(i))*1.1 para.gray_matter_mus(para.ind(i)) para.gray_matter_g(para.ind(i)), 1.37   % medium 3: gray matter, perturbed
%         para.brain_mua(para.ind(i))*1.1 para.white_matter_mus(para.ind(i)) para.white_matter_g(para.ind(i)), 1.37];   % medium 4: white matter, pertubed
%         else
%         cfg.prop = [0 0 1 1            % medium 0: the environment
%         para.brain_mua(para.ind(i)) para.gray_matter_mus(para.ind(i)) para.gray_matter_g(para.ind(i)), 1.37   % medium 1: gray matter
%         para.brain_mua(para.ind(i)) para.white_matter_mus(para.ind(i)) para.white_matter_g(para.ind(i)), 1.37   % medium 2: white matter
%         para.brain_mua(para.ind(i))+cfg.perturb para.gray_matter_mus(para.ind(i)) para.gray_matter_g(para.ind(i)), 1.37   % medium 3: gray matter, perturbed
%         para.brain_mua(para.ind(i))+cfg.perturb para.white_matter_mus(para.ind(i)) para.white_matter_g(para.ind(i)), 1.37];   % medium 4: white matter, perturbed
%         end
%     case {'with skull', 'with skull off focus'}
        if isnan(cfg.perturb)
        cfg.prop=[0 0 1 1            % medium 0: the environment
        para.skull_mua(para.ind(i)) para.skull_mus(para.ind(i)) 0.9337, 1.56 % medium 1: skull
        para.brain_mua(para.ind(i)) para.gray_matter_mus(para.ind(i)) para.gray_matter_g(para.ind(i)), 1.37   % medium 2: gray matter
        para.brain_mua(para.ind(i)) para.white_matter_mus(para.ind(i)) para.white_matter_g(para.ind(i)), 1.37   % medium 3: white matter
        para.skull_mua(para.ind(i))*1.1 para.skull_mus(para.ind(i)) 0.9337, 1.56 % medium 4: skull, perturbed
        para.brain_mua(para.ind(i))*1.1 para.gray_matter_mus(para.ind(i)) para.gray_matter_g(para.ind(i)), 1.37   % medium 5: gray matter, perturbed
        para.brain_mua(para.ind(i))*1.1 para.white_matter_mus(para.ind(i)) para.white_matter_g(para.ind(i)), 1.37];   % medium 6: white matter, perturbed
        else
        cfg.prop=[0 0 1 1            % medium 0: the environment
        para.skull_mua(para.ind(i)) para.skull_mus(para.ind(i)) 0.9337, 1.56 % medium 1: skull
        para.brain_mua(para.ind(i)) para.gray_matter_mus(para.ind(i)) para.gray_matter_g(para.ind(i)), 1.37   % medium 2: gray matter
        para.brain_mua(para.ind(i)) para.white_matter_mus(para.ind(i)) para.white_matter_g(para.ind(i)), 1.37   % medium 3: white matter
        para.skull_mua(para.ind(i))+cfg.perturb para.skull_mus(para.ind(i)) 0.9337, 1.56 % medium 4: skull, perturbed
        para.brain_mua(para.ind(i))+cfg.perturb para.gray_matter_mus(para.ind(i)) para.gray_matter_g(para.ind(i)), 1.37   % medium 5: gray matter, perturbed
        para.brain_mua(para.ind(i))+cfg.perturb para.white_matter_mus(para.ind(i)) para.white_matter_g(para.ind(i)), 1.37];   % medium 6: white matter, perturbed
        end
%     otherwise 
%         disp('Check property setting (cfg.prop)')
% end
           
% calculate total photon weights after perturbation   
newcfg              = cfg;
newcfg.seed         = seeds{i}.data(:,det1{i}.DetInd);
newcfg.detphotons   = det1{i}.data(:,det1{i}.DetInd);
% newcfg.outputtype   = 'jacobian';

switch method
    case 'z'
        tic, [~, det_temp1, ~, ~,traj1] = mcxlab(newcfg); toc
        
        traj_group = XINTRINSIC_SIM_PlotTraj(traj1, cfg.bounds, [1 400 1 500], 0);
        max_temp = cellfun(@max,traj_group,'UniformOutput',false); % max along 3 axis
        det_temp1.MaxDepth = cellfun(@(x) x(3), max_temp); % take only z axis
        %  ========== change properties in Z direction  ==========
        S.NewDetWeight{i} = zeros(1,SimVol);
        fwait = waitbar(0,'Perturbation simulation started ...');
    %     SimStart = min(find(squeeze(cfg.vol(1,1,:)) == 2)); % the first gray matter layer as the layer to start
    
        
        SimStart = 1;
        for zz = SimStart : SimStart + SimVol - 1
            % ========== define tissue structure ==========
            waitbar((zz-SimStart+1)/SimVol,fwait,['Perturbation simulation processing ',num2str((zz-SimStart+1)),'/',num2str(SimVol)]);
            newcfg.vol = cfg.vol;
            switch cfg.config
                case {'without skull', 'without skull off focus'}
                    if zz <= cfg.bounds
                        newcfg.vol(:,:,zz) = 5; % gray matter, perturbed
                    else
                        newcfg.vol(:,:,zz) = 6; % white matter, perturbed
                    end
                    newcfg.vol=uint8(newcfg.vol);
                case {'with skull', 'with skull off focus'}
                    if zz <= cfg.bounds(1)
                        newcfg.vol(:,:,zz) = 4; % skull, perturbed
                    elseif zz > cfg.bounds(1) && zz <= cfg.bounds(2)
                        newcfg.vol(:,:,zz) = 5; % gray matter, perturbed           
                    else
                        newcfg.vol(:,:,zz) = 6; % white matter, perturbed
                    end
                    newcfg.vol=uint8(newcfg.vol);
                otherwise
            end
            tic, [S.f2{i}, S.det2{i}, S.vol2, S.seeds2, S.traj2{i}] = mcxlab(newcfg); toc
                    % calculate baseline measurement (sum of weights)
%             det_temp1 = det1{i};
%             det_temp1.ppath = det_temp1.ppath(det_temp1.DetInd,:); % select photons by its exit angle
%             det_temp1.nscat = det_temp1.nscat(det_temp1.DetInd,:); % select photons by its exit angle
%             det_temp1.p = det_temp1.p(det_temp1.DetInd,:); % select photons by its exit angle
            det_temp2 = S.det2{i};
            det_temp2.MaxDepth = det_temp1.MaxDepth;
%             traj_group = XINTRINSIC_SIM_PlotTraj(S.traj2{2}, cfg.bounds, [1 400 1 500], 0);
%             max_temp = cellfun(@max,traj_group,'UniformOutput',false); % max along 3 axis
%             det_temp2.MaxDepth = cellfun(@(x) x(3), max_temp); % take only z axis
            
            
            [C,ia,ib] = intersect(det_temp1.p, det_temp2.p, 'rows'); % make sure the 2nd replay simulation uses only the detected photons in 1st simulation
            S.nDet = size(C,1);
            det_temp1.ppath = det_temp1.ppath(ia,:); % select photons by its exit angle
            
            det_temp1.chi = para.chi;
            S.DetWeight(i) = sum(mcxdetweight_polar(det_temp1,cfg.prop));
            % calculate perturbed measurement (sum of weights)
            det_temp2.ppath = det_temp2.ppath(ib,:);
            det_temp2.nscat = det_temp2.nscat(ib,:);
            
            det_temp2.chi = para.chi;
            S.NewDetWeight{i}(zz-SimStart+1) = sum(mcxdetweight_polar(det_temp2,newcfg.prop));
            profile{i}(zz-SimStart+1) = -(S.NewDetWeight{i}(zz-SimStart+1) - S.DetWeight(i)).*100./S.DetWeight(i);
        end
%     profile{i} = -(S.NewDetWeight{i} - S.DetWeight(i)).*100./S.DetWeight(i);
    x = (1:SimVol).*cfg.unitinmm;

    case 'x'
        SimCenter = cfg.size(1)/2-1; 
        % ========== change properties in X direction  ==========
        S.NewDetWeight{i} = zeros(1,SimVol);
        fwait = waitbar(0,'Perturbation simulation started ...');

        for xx = SimCenter-SimVol/2: SimCenter+SimVol/2 
            waitbar((xx-(SimCenter-SimVol/2))/SimVol,fwait,['Perturbation simulation processing ',num2str((xx-(SimCenter-SimVol/2))),'/',num2str(SimVol)]);

            % ========== define tissue structure ==========
            newcfg.vol = cfg.vol;
            switch cfg.config
                case {'without skull', 'without skull off focus'} % perturbate along x, at the y-value of SimCenter+1 to avoid the dip artifact 
                    newcfg.vol(xx,SimCenter+1,1:cfg.bounds) = 5; % gray matter, perturbed
    %                 newcfg.vol(xx,SimCenter+1,cfg.bounds+1:end) = 6; % white matter, perturbed
                    newcfg.vol=uint8(newcfg.vol);
                case {'with skull', 'with skull off focus'}
                    newcfg.vol(xx,SimCenter+1,cfg.bounds(1)+1:cfg.bounds(2)) = 5; % gray matter, perturbed
    %                 newcfg.vol(xx,SimCenter+1,cfg.bounds(2)+1:end) = 6; % white matter, perturbed
                    newcfg.vol=uint8(newcfg.vol);
                otherwise
            end
            tic, [S.f2{i}, S.det2{i}, S.vol2, S.seeds2, S.traj2{i}] = mcxlab(newcfg); toc
            % calculate baseline measurement (sum of weights)
            det_temp1 = det1{i};
            det_temp1.ppath = det_temp1.ppath(det_temp1.DetInd,:); % select photons by its exit angle
            det_temp1.p = det_temp1.p(det_temp1.DetInd,:); % select photons by its exit angle
            det_temp2 = S.det2{i};

            [C,ia,ib] = intersect(det_temp1.p, det_temp2.p, 'rows');
            S.nDet = size(C,1);
            det_temp1.ppath = det_temp1.ppath(ia,:); % select photons by its exit angle
            det_temp1.chi = para.chi;
            S.DetWeight(i) = sum(mcxdetweight_polar(det_temp1,cfg.prop));
            % calculate perturbed measurement (sum of weights)
            det_temp2.ppath = det_temp2.ppath(ib,:);
            det_temp2.chi = para.chi;
            S.NewDetWeight{i}(xx-(SimCenter-SimVol/2)+1) = sum(mcxdetweight_polar(det_temp2,newcfg.prop));
            profile{i}(xx-(SimCenter-SimVol/2)+1) = -(S.NewDetWeight{i}(xx-(SimCenter-SimVol/2)+1) - S.DetWeight(i)).*100./S.DetWeight(i);
        end
    %     profile{i} = -(S.NewDetWeight{i} - S.DetWeight(i)).*100./S.DetWeight(i);
        x = [-SimVol/2:SimVol/2].*cfg.unitinmm;
    
    case 'r'
        SimCenter = cfg.size(1)/2-1; % the illumination center is cfg.size(1)/2-1
        % ========== change properties in R direction  ==========
        S.NewDetWeight{i} = zeros(1,SimVol/2);
        fwait3 = waitbar(0,'Horizontal perturbation simulation started ...');
        profile{i} = zeros(1,SimVol/2+1);
        % ========== assign r values to x and y ==========
        R = zeros(SimVol/2, SimVol/2);
        for xx = 1:SimVol/2
            for yy = 1:SimVol/2
                R(xx,yy) = round(vecnorm([xx,yy])); % the integer diameter value for each pixel, 1st quadrant
            end
        end
        temp = 1:SimVol/2;
        RR = [rot90(R,2), flipud(temp'), flipud(R); fliplr(temp), 0, temp; fliplr(R), temp', R ];
        r_real = zeros(1, SimVol/2+1);
        for rr = 0 : SimVol/2
            waitbar(rr/(SimVol/2),fwait3,['Horizontal perturbation simulation processing ',num2str(rr),'/',num2str(SimVol/2)]);
            % ========== define tissue structure ==========
            newcfg.vol = cfg.vol;
            [xx, yy, ~] = find(RR == rr); 
            xx = xx - SimVol/2 -1 + SimCenter; 
            yy = yy - SimVol/2 -1 + SimCenter; 
            switch cfg.config
                case {'without skull', 'without skull off focus'}
                    for k = 1:length(xx)
                        newcfg.vol(xx(k),yy(k),1:cfg.bounds) = 5; % gray matter, perturbed
    %                     newcfg.vol(xx(k),yy(k),cfg.bounds+1:end) = 6; % white matter, perturbed
                    end     
                case {'with skull', 'with skull off focus'}
                    for k = 1:length(xx)
        %                 newcfg.vol(xx(k),yy(k),cfg.bounds(1)+1:cfg.bounds(2)) = 4; % skull matter, perturbed
                        newcfg.vol(xx(k),yy(k),cfg.bounds(1)+1:cfg.bounds(2)) = 5; % gray matter, perturbed
    %                     newcfg.vol(xx(k),yy(k),cfg.bounds(2)+1:end) = 6; % white matter, perturbed
                    end            
                otherwise
            end
            newcfg.vol=uint8(newcfg.vol);
            xx = xx - SimCenter; yy = yy - SimCenter;
            r_real(rr+1) = mean(vecnorm([xx,yy]'));

        %     VolDiff = newcfg.vol(:,:,cfg.bounds(1)+1) - cfg.vol(:,:,cfg.bounds(1)+1);
        %     nChange(rr+1) = length(find(VolDiff(:)));
            nChange = length(xx); % how many columes were perturbed
            tic, [S.f2{i}, S.det2{i}, S.vol2, S.seeds2, S.traj2{i}] = mcxlab(newcfg); toc

            % calculate baseline measurement (sum of weights)
            det_temp1 = det1{i};
            det_temp1.ppath = det_temp1.ppath(det_temp1.DetInd,:); % select photons by its exit angle
            det_temp1.p = det_temp1.p(det_temp1.DetInd,:); % select photons by its exit angle
            det_temp2 = S.det2{i};

            [C,ia,ib] = intersect(det_temp1.p, det_temp2.p, 'rows');
            S.nDet = size(C,1);
            det_temp1.ppath = det_temp1.ppath(ia,:); % select photons by its exit angle
            det_temp1.chi = para.chi;
            S.DetWeight(i) = sum(mcxdetweight_polar(det_temp1,cfg.prop));
            % calculate perturbed measurement (sum of weights)
            det_temp2.ppath = det_temp2.ppath(ib,:);
            det_temp2.chi = para.chi;
            S.NewDetWeight{i}(rr+1) = sum(mcxdetweight_polar(det_temp2,newcfg.prop));
            profile{i}(rr+1) = -(S.NewDetWeight{i}(rr+1) - S.DetWeight(i))./nChange.*100./S.DetWeight(i);
        end

        x = [-fliplr(r_real(2:end)), r_real].*cfg.unitinmm;
    %     profile{i} = -(S.NewDetWeight{i} - S.DetWeight(i))./nChange.*100./S.DetWeight(i);
        profile{i} = [fliplr(profile{i}(2:end)) profile{i}];
    otherwise
end

end

end
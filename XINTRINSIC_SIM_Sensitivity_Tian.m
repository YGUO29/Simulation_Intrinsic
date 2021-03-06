%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCXLAB - Monte Carlo eXtreme for MATLAB/Octave by Qianqina Fang
% 
% Depth sensitivity based on 2011 Tian paper
% 
% This file is part of Monte Carlo eXtreme (MCX) URL:http://mcx.sf.net
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear cfg;
%% set tissue parameters
load('C:\Users\guoyu\Documents\MATLAB\Simulation_Intrinsic\spectra.mat')
lambda = spectra(:,1); e_HbO = spectra(:,2); e_HbR = spectra(:,3); 
c_Hb   = 129; % hemoglobin (g/L) used in Hilman 2016
perc_blood = 0.03; % blood volume / tissue volume
sat_O2 = 0.75;
c_HbO = c_Hb*sat_O2; c_HbR = c_Hb*(1-sat_O2);
mua_HbR = perc_blood.*(e_HbR.*c_HbR.*2.303./64500)./10; % unitsL 1/mm
mua_HbO = perc_blood.*(e_HbO.*c_HbO.*2.303./64500)./10;
mua_HbT = mua_HbO+mua_HbR;


figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2),
subplot(121), hold on, set(gca, 'YScale', 'log')
semilogy(lambda, mua_HbR), semilogy(spectra(:,1),mua_HbO),
semilogy(lambda,mua_HbT,'k')
scatter(lambda(ind),mua_HbT(ind),32)
for i = 1:length(ind)
text(lambda(ind(i))+5,mua_HbT(ind(i)),num2str(mua_HbT(ind(i)), 2))
end
title('Absorption coefficient')
legend('mua-HbR', 'mua-HbO2', 'mua-HbT')
subplot(122), hold on, set(gca, 'YScale', 'log')
semilogy(spectra(:,1),e_HbR), hold on, semilogy(spectra(:,1),e_HbO)
for i = 1:length(ind); xline(lambda(ind(i)),'--'); end
title('Moler extinction coefficient')
legend('\epsilon-HbR', '\epsilon-HbO2')


load('C:\Users\guoyu\Documents\MATLAB\Simulation_Intrinsic\Skull_mus.mat')
g_skull = 0.89;
lambda_skull = Skull_mus(:,1);
musr_skull = Skull_mus(:,2);
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2),
plot(lambda_skull,musr_skull)
mus_skull = interp1(lambda_skull,musr_skull,WLs,'linear','extrap')./(1-g_skull);
%% preparing the input data

    % set seed to make the simulation repeatible
    cfg.seed=hex2dec('623F9A9E'); 

    cfg.nphoton=5e7;
    cfg.unitinmm = 0.05; % define units, 10 um
    % time-domain simulation parameters
    cfg.tstart=0;
    cfg.tend=5e-9;
    cfg.tstep=5e-10;
    
% ========== define source and detector ==========
%     cfg.srctype = 'pencil';
    NA = 0.05;
    cfg.srctype = 'cone';
    cfg.srcparam1=[asin(NA) 0 0 0];
    % cfg.srcparam1 = 100; %radius for the disk
    cfg.size = [300 300 300]; % 15mm cube
    cfg.srcpos=[cfg.size(1)/2, cfg.size(2)/2, 1];
    cfg.srcdir=[0 0 1];
    cfg.detpos=[cfg.size(1)/2, cfg.size(2)/2, 0, min(cfg.size(1)/2, cfg.size(2)/2)]; % [x y z radius]
    cfg.maxdetphoton = cfg.nphoton;
    cfg.maxjumpdebug = cfg.nphoton;
    cfg.maxexitangle = 0.05;
% ========== define a 3 layer structure ==========
    cfg.vol = ones(cfg.size); % 15 mm volume
    cfg.layersize = 8; % 0.4mm per layer
    bounds1 = cfg.layersize + 0:cfg.layersize:cfg.size(3); bounds1 = bounds1(1:4);
    for l = 1:length(bounds1)-1
        cfg.vol(:,:,bounds1(l)+1:bounds1(l+1))=l+1;
    end       
    cfg.vol(:,:,bounds1(l+1):end) = l+2;
%     bounds1 = [floor(0.7/cfg.unitinmm) floor(0.7/cfg.unitinmm)+floor(0.1/cfg.unitinmm)]
%     cfg.vol(:,:,bounds1(1)+1:bounds1(2))=2; % 0.7mm: skin & skull, 0.1mm CSF
%     cfg.vol(:,:,bounds1(2)+1:end)=3; % the rest: gray matter
%     cfg.vol=uint8(cfg.vol);
    cfg.isreflect=0; % disable reflection at exterior boundary
    % GPU thread configuration
    cfg.autopilot=1;
    cfg.gpuid=1;
% ========== output control ==========
    cfg.issaveexit = 1;
%%    
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w')
% for i = 1:length(ind)
for i = 1
    prop_temp = [mua_HbT(ind(i)) 21 0.82, 1.37];
    cfg.prop=[0 0 1 1            % medium 0: the environment
%        0.019 7.8   0.89 1.55     % medium 1: skin & skull (n changed from 1.37 to 1.55)
%        0.004 0.009 0.89 1.37     % medium 2: CSF
       repmat(prop_temp,l+2,1)   ];   % medium 3: gray matter
   
    fprintf('running simulation ... this takes about 35 seconds on a GTX 470\n');
    tic;
    % f1=mcxlab(cfg);
    [f1,det1,vol,seeds,traj]=mcxlab(cfg);
    toc;
% ========== detected photons ==========
    ExitAngle   = sin(vecnorm(det1.v(:,1:2)')./vecnorm(det1.v')); % sin_theta (NA)
    DetInd      = find(ExitAngle < NA);
    DetNumber   = length(DetInd);
    DetExitPosition = det1.p(DetInd,1:2);

% ========== plot figure ==========
%     subplot(2,3,i);
%     I1 = log10(squeeze(sum(f1.data(:,size(cfg.vol,1)./2+1,:,:),4))');
%     I1(I1<0) = 0;
%         if i == 1
%         cmax1 = max(I1(:)); cmin1 = min(I1(:));
%         end
%     contourf(I1, cmin1:0.5:cmax1);
%     hold on
%     plot([0 size(cfg.vol,1)],[bounds1(1) bounds1(1)],'--',...
%         [0 size(cfg.vol,1)],[bounds1(2) bounds1(2)],'--');
%     title(['flux with no reflection, no skull, wavelength ',num2str(WLs(i)),'nm']);
%     set(gca,'clim',[cmin1 cmax1]);
%     colorbar
%     drawnow
% ========== plot depth sensitivity ==========
%     subplot(3,6,6+i);
%     h1 = histogram2(DetExitPosition(:,1),DetExitPosition(:,2),100,...
%         'XBinLimits',[cfg.detpos(1)-10, cfg.detpos(1)+10],'YBinLimits',[cfg.detpos(2)-10, cfg.detpos(2)+10],...
%         'DisplayStyle','tile','ShowEmptyBins','on');
%     DetExitProfile = sum(h1.Values,2);
%         [mm, ii] = max(DetExitProfile); cutoff = mm./2;
%         [~,xx1] = min(abs(DetExitProfile(1:ii)-cutoff));
%         [~,xx2] = min(abs(DetExitProfile(ii:end)-cutoff));    
%     DetExitSize = (xx2 + ii - xx1).*cfg.unitinmm;
%     title(['Exit position, FWHM: ',num2str(DetExitSize,3),'mm'])
%     
%     subplot(3,6,12+i);
%     DetWeight = mcxdetweight(det1,cfg.prop);
%     DetPathlength = det1.ppath.*repmat(DetWeight(:),1,size(det1.ppath,2));
%     DetMeanPathlength = cfg.unitinmm*sum(DetPathlength(:,2)) / sum(DetWeight(:));
%     h2 = histogram(DetPathlength);
%     title(['Pathlength distribution, mean: ',num2str(DetMeanPathlength,3),'mm']) 
end
%%
clear MeanWeight
det2 = det1;
det2.ppath = det1.ppath(DetInd,:);
for i = 1:length(ind)
for i_layer = 1:length(bounds1)+1
    cfg.prop(2:end,1) = 0;
    cfg.prop(i_layer+1,1) = mua_HbT(ind(i));   
    DetWeight = mcxdetweight(det2,cfg.prop);
    MeanWeight(i_layer) = 1 - sum(DetWeight)/DetNumber;
%     DetWeight = mcxdetweight(det1,cfg.prop);
%     DetPathlength = det1.ppath.*repmat(DetWeight(:),1,size(det1.ppath,2));
%     DetMeanPathlength(i) = cfg.unitinmm*sum(DetPathlength(:,2)) / sum(DetWeight(:));
%     h2 = histogram(DetPathlength);
end
hold on, plot(MeanWeight./sum(MeanWeight))
drawnow
end
set(gca,'XTickLabels',{'0-400um', '400-800um', '800-1200um', '1200-1600um', '1600-15000um'})


[a, bg] = LoadNorm_tMC( 'semi-inf-td' );

g = 0.01;
mus = 1;
musp = mus * (1-g);
mua = 0.005;

v = 3e11;
D = v/(3*musp);

dt = 0.1e-9;
t = dt:dt:5e-9;

rho = 15;
zdet = 10;

zo = 1.0/musp;


% Plot TPSF theory and MC

% From Patterson 1989
% flux out
Jmeas = squeeze(-a(30+rho,30,1,:));

Jtheory = (4*pi*D).^(-3/2) * zo * t.^(-5/2) .* exp(-mua*v*t) .* ...
	  exp( -((rho-1)^2+zo^2)./(4*D*t) );

% fluence in
r2 = ((rho-1).^2 + (zdet-zo)^2);
ri2 = ((rho-1).^2 + (zdet+zo+4/(3*musp))^2);
P2Dtheory = v*(4*pi*D*t).^(-3/2) .* exp(-v*mua*t) .* ...
      ( exp(- r2./(4*D*t) ) - ...
	exp(- ri2./(4*D*t) )  );

hf=figure(3);
hold off
h=semilogy((t-dt/2)*1e9, Jmeas, 'k.', t*1e9, Jtheory, 'k-',  ...
	 (t-dt/2)*1e9, squeeze(a(30+rho,30,zdet,:)), 'r.', t*1e9, P2Dtheory, 'r-' );
legend( 'J_{out} Measured', 'J_{out} Theory' );

set(h(1),'MarkerSize',15);
set(h(2),'Linewidth',3);
set(h(3),'MarkerSize',15);
set(h(4),'Linewidth',3);
set(get(hf,'CurrentAxes'),'FontSize',20); 
xlabel('Time (ns)');
ylabel('Flux and Fluence');
legend( 'Monte Carlo', 'Diffusion Theory' );
print -djpeg90 semi_inf_td1.jpg



% Plot contours of Half Max

tt = t - dt/2;;
x = -29:1:30;
z = 1:1:60;
[xx, zz] = meshgrid( x, z);
rr2 = (xx-0.5).^2 + (zz-zo).^2;
rri2= (xx-0.5).^2 + (zz+(zo+4/(3*musp))).^2;
for tidx=1:length(t)
  P3Dtheory(:,:,tidx) = v*(4*pi*D*tt(tidx))^(-3/2) * exp(-v*mua*tt(tidx)) * ...
      ( exp(- rr2/(4*D*tt(tidx)) ) - ...
	exp(- rri2/(4*D*tt(tidx)) )  );
end

hf=figure(4)
hold off
for tidx = [1 5 10 15 20]
  HM = max(max(a(:,30,:,tidx))) / 2;
  [c h1]=contour( max(squeeze(a(:,30,:,tidx))',1e-12), [HM HM], 'k:' );
  set(h1(:,1),'LineWidth',3);
  hold on
  HM = max(max(P3Dtheory(:,30,tidx))) / 2;
  [c h2]=contour( max(squeeze(P3Dtheory(:,:,tidx)),1e-12), [HM HM], 'k-' );
  set(h2(:,1),'LineWidth',3);
end
hold off
set(get(hf,'CurrentAxes'),'FontSize',20); 
xlabel('Lateral Position (mm)');
ylabel('Depth (mm)');
legend([h1(1); h2(1)], 'Fluence Monte Carlo','Diffusion Theory');
print -djpeg90 semi_inf_td2.jpg






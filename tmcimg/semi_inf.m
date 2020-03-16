[a, bg] = LoadNorm_tMC( 'semi-inf' );

g = 0.01;
mus = 1;
musp = mus * (1-g);
mua = 0.005;

k = sqrt(3*musp*mua);

Jmeas = -a(:,30,1);
Jcalc = a(:,30,2)/4 + 1/(6*musp) * (a(:,30,3) - a(:,30,2));

x = -29:30;
z = 0:59;

zo = 1.0/musp;

r = sqrt((x-0.5).^2 + zo^2);
ri = sqrt((x-0.5).^2 + (zo+4/(3*musp))^2);
Ptheory = 3*musp* (exp( -k*r ) ./ (4*pi*r) - exp( -k*ri ) ./ (4*pi*ri));
Jtheory = Ptheory/2;

[xx, zz] = meshgrid( x, z);
rr = ((xx-0.5).^2 + (zz-zo).^2).^0.5;
rri= ((xx-0.5).^2 + (zz+(zo+4/(3*musp))).^2).^0.5;
P3Dtheory = 3*musp* (exp( -k*rr ) ./ (4*pi*rr) - exp( -k*rri ) ./ (4*pi*rri));

P3Dtheory(1,:) = zeros(1,size(P3Dtheory,2));

hf=figure(1);
hold off
h=semilogy(x(31:60), Jmeas(31:60), 'k.', x(31:60), Jtheory(31:60), 'k-' );
set(h(1),'MarkerSize',15);
set(h(2),'Linewidth',3);
set(get(hf,'CurrentAxes'),'FontSize',20); 
xlabel('Separation (mm)');
ylabel('Remitted Flux');
legend( 'J_{out} Monte Carlo', 'J_{out} Diffusion Theory' );
print -djpeg90 semi_inf1.jpg


hf=figure(2)
hold off
clines = -1.5:-0.5:-4;
[c h1]=contour( log10(max(squeeze(a(:,30,:))',1e-8)), clines, 'k:' );
set(h1(:,1),'LineWidth',3);
hold on
[c h2]=contour( log10(max(squeeze(P3Dtheory),1e-8)), clines, 'k-' );
set(h2(:,1),'LineWidth',3);
hold off
set(get(hf,'CurrentAxes'),'FontSize',20); 
xlabel('Lateral Position (mm)');
ylabel('Depth (mm)');
legend([h1(1); h2(1)], 'Fluence Monte Carlo','Diffusion Theory');
print -djpeg90 semi_inf2.jpg

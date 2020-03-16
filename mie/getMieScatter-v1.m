function musgp = getMieScatter(lambda, dia, fv, npar,nmed)
% function musgp = getMieScatter(lambda, dia, fv)
%  fv            = volume fraction of spheres in medium (eg., fv = 0.05)
%  lambda        = wavelength in um (eg., lambda = 0.633)
%  dia           = sphere diameter in um (eg., dia_um = 0.0500)
%  npar          = particle refractive index (eg. polystyrene = 1.57)
%  nmed          = medium refractive index (eg., water = 1.33)
%                  Note: npar and nmed can be imaginary numbers.
%  returns musgp = [mus g musp]  
%       mus      = scattering coefficient [cm^-1]
%       g        = anisotropy of scattering [dimensionless]
%       musp     = reduced scattering coefficient [cm^-1]
%  Uses
%       Mie.m, which uses mie_abcd.m, from Maetzler 2002
%       
% - Steven Jacques, 2009

% Corrected some errors, added table format printing of selected
% parameters. -Rob Brown, Rockwell Collins, June 2017

Vsphere = 4/3*pi*(dia/2)^3;     % volume of sphere
rho     = fv/Vsphere;           % #/um^3, concentration of spheres

m = npar/nmed;                  % ratio of refractive indices
x = pi*dia/(lambda/nmed);       % ratio circumference/wavelength in medium

u = mie(m, x)';                 % <----- Matzler's subroutine
% u = [real(m) imag(m) x qext qsca qabs qb asy qratio];

qsca = u(5);                    % scattering efficiency, Qsca
g    = u(8);                    % anisotropy, g

A       = pi*dia^2/4;           % geometrical cross-sectional area, um^2
sigma_s = qsca*A;               % scattering cross-section, um^2
mus     = sigma_s*rho*1e4;      % scattering coeff. cm^-1
musp    = mus*(1-g);            % reduced scattering coeff. cm^-1

if 0 % 1 = print full report, 0 = disable
    disp('----- choice:')
    fprintf('lambda  \t= %0.3f um\n', lambda)
    fprintf('dia     \t= %0.3f um\n', dia)
    fprintf('rho     \t= %0.3f #/um^3\n', rho)
    fprintf('npar    \t= %0.3f\n', npar)
    fprintf('nmed    \t= %0.3f\n', nmed)
    disp('----- result:')
    fprintf('real(m) \t= %0.3f\n', u(1))
    fprintf('imag(m) \t= %0.3e\n', u(2))
    fprintf('x       \t= %0.3e\n', u(3))
    fprintf('qext    \t= %0.3e\n', u(4))
    fprintf('qsca    \t= %0.3e\n', u(5))
    fprintf('qabs    \t= %0.3e\n', u(6))
    fprintf('qb      \t= %0.3e\n', u(7))
    fprintf('asy     \t= %0.4f\n', u(8))
    fprintf('qratio  \t= %0.3e\n', u(9))
    disp('----- optical properties:')
    fprintf('mus     \t= %0.3f cm^-1\n', mus)
    fprintf('g       \t= %0.4f\n', g)
    fprintf('musp    \t= %0.3f cm^-1\n', musp)
end
%Print the just parameters you want, in table format. Mod header from
%demoMie (line 15) to properly label columns
%fprintf('%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f  \n', lambda, dia, rho, npar, nmed, u(1),u(2),u(3),u(4),u(5),u(6),u(7),u(8),u(9), mus, g, musp)
fprintf('%0.3f %0.3f %0.6f %0.3f %0.3f \n', lambda, dia, rho, mus, g)

musgp= real([mus g musp]);


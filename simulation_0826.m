N = 1000; % number of photons each run
Threshold = 0.01; % below which weight is determined using Roulette
Chance = 0.1; % used in Roulette
One_Minus_Coszero = 1e-12; 
mus = 10.0; % Scattering coefficient [mm^-1]
mua = 0.3; % Absorption coefficient [mm^-1]
g = 0.9; % [dimensionless]
w = 1.0; % initialize photon weight
albedo = mus./ (mus + mua); 
mut = mus + mua;
depth = 0.0; % [mm]
LL=[];
CC=[];
dd=[];

NR = 100; % set number of bins
nt = 1.33; % refractive index of tissue [dimensionless]
radial_size = 3.0; % total range over which bin extends [mm]
dr = radial_size / NR; % [mm]
Csph = zeros(NR, 1);
Ccyl = zeros(NR, 1);
Cpla = zeros(NR, 1);


hold on

for i = 1:1
    Photon_Num = 1;
    sum = 0;
    Count_Num = 0;
    
    Light_Density={};
    while Photon_Num <= N
        w = 1.0; % initialize photon weight
        x = 0;
        y = 0;
        z = 0;
        
        path = [x, y, z, w]; % intialize path
        L = 0;  % Pathlength [mm]
        depth = 0;
        Photon_Status = true;
        
        %{ assume theta = zero for first move %}
        costheta = 1;
        sintheta = 0;
        psi = 2 * pi * rand;
    
        ux = sintheta * cos(psi);
        uy = sintheta * sin(psi);
        uz = costheta;

        Photon_Num = Photon_Num + 1;
        
        while Photon_Status && z >= 0
            xi = rand;
            s = -log(xi) / mut;
            x  = x + s * ux;
            y  = y + s * uy;
            z  = z + s * uz;
            path = [path;x, y, z, w]; % record of path
            L = L + s;
            
            if depth <= z
                depth = z;
            end 
            
            if z < 0
                Count_Num = Count_Num + 1;
                sum = sum + L;
            end
            
            absorb = 1 - w .* (1-albedo);
            w = w - absorb;
            
            % spherical bin
            r = sqrt(x.^2 + y.^2 + z.^2);
            ir = int8(r ./ dr) + 1; %array index starting from 1 not 0
            if (ir >= NR)
                ir = NR;
            end
            Csph(ir) = Csph(ir) + absorb;
            
            % cylindrical bin
            r = sqrt(x.^2 + y.^2);
            ir = int8(r ./ dr) + 1;
            if (ir >= NR)
                ir = NR;
            end
            Ccyl(ir) = Ccyl(ir) + absorb;
            
            % planar bin
            r = abs(z);
            ir = int8(r ./ dr) + 1;
            if (ir >= NR)
                ir = NR;
            end
            Cpla(ir) = Cpla(ir) + absorb;
            
            if w > 0
                Photon_Status = true;
            else 
                Photon_Status = false;
            end 
            
            %{ fix g = 0.9 %}
            costheta=((1+g.^2-((1-g.^2)./(1-g+2*g.*rand)).^2)./(2.*g));
            sintheta = sqrt(1 - costheta.^2);

            %{ sample psi %}
            psi=2*pi*rand;
            cospsi = cos(psi);
            
            if psi < pi
                sinpsi = sqrt(1.0 - cospsi.^ 2);
            else 
                sinpsi = - sqrt(1.0 - cospsi.^ 2);
            end

            if 1 - abs(uz) <= One_Minus_Coszero
                uxx = sintheta * cospsi;
                uyy = sintheta * sinpsi;
                uzz = costheta * sign(uz);
            else
                temp = sqrt(1.0 - uz.^ 2);
                uxx = sintheta.* (ux.* uz.* cospsi - uy.* sinpsi) / temp + ux.* costheta;
                uyy = sintheta.* (uy.* uz.* cospsi + ux.* sinpsi) / temp + uy.* costheta;
                uzz = - sintheta.* cospsi.* temp + uz.* costheta;
            end

            %{ update trajectory %}
            ux = uxx;
            uy = uyy;
            uz = uzz;
            
            if w < Threshold
                if rand <= Chance
                    w = w./Chance;
                else 
                    Photon_Status = false;
                end 
            end
        end
        dd = [dd, depth];
        plot(path(:,1),path(:,3))
        Light_Density{Photon_Num - 1} = path;
    end
end

hold off

avg = sum / Count_Num;
if Count_Num > 0
    LL = [LL, avg];
    CC = [CC, Count_Num];
end


    
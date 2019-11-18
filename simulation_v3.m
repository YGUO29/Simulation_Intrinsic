N = 1e2; % number of photons each run
Threshold = 0.01; % below which weight is determined using Roulette
Chance = 0.1; % used in Roulette
One_Minus_Coszero = 1e-12; 
mus = 10.0; % Scattering coefficient [mm^-1]
mua = 0.3; % Absorption coefficient [mm^-1]
mut = mus + mua;
g = 0.9; % [dimensionless]
w = 1.0; % initialize photon weight
factor = 1 - mua / mut;
depth = 0.0; % [mm]
LL=[];
CC=[];
dd=[];

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
            
            w = factor * w;
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

% histogram(dd)
%%
c = 0.1; % the density of the contour
range = c * [0.9, 1.1];
hold on
for i = 1:N
    index = and(Light_Density{i}(:,4) >= range(1), Light_Density{i}(:,4) <= range(2)); % find the points in the range
    scatter(Light_Density{i}(index, 1), Light_Density{i}(index, 3), 'b', 'filled') % scatter (x, z)
end
hold off  


    
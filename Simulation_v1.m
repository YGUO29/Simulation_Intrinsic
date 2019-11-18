Threshold = 0.01; % below which weight is determined using Roulette
Chance = 0.1; % used in Roulette
One_Minus_Coszero = 1e-12; 
Photon_Num = 1;
mus = 1.0; % Scattering coefficient
mua = 0.3; % Absorption coefficient
mut = mus + mua;
g = rand; % anisotropy coefficient
% g = 0.5;
w = 1.0; % initialize photon weight
factor = 1 - mua / mut;
Count_Num = 0;
depth = 0.0;
sum = 0.0;

while Photon_Num <= 1000000
    x = 0;
    y = 0;
    z = 0;
    L = 0.0;  % Pathlength
    depth = 0;
    Photon_Status = true;
    
    costheta = 2.0 * rand - 1.0;
    sintheta = sqrt(1.0 - costheta ^ 2);
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
        L = L + s;
        if depth <= z
            depth = z;
        end
        
        if z < 0.0
            Count_Num = Count_Num + 1;
            sum = sum + depth;
        end
        
        w = factor * w;
        if w > 0.0
            Photon_Status = true;
        else 
            Photon_Status = false;
        end
        
        if g == 0.0
            costheta = 2.0 * rand - 1.0;
        else 
            temp = (1.0 - g ^ 2) / (1.0 - g + 2 * g * rand);
            costheta = (1.0 + g ^ 2 - temp ^ 2) / (2.0 * g);
        end
        sintheta = sqrt(1.0 - costheta ^ 2);
        
        %{ sample psi %}
        psi = 2.0 * pi * rand;
        cospsi = cos(pi);
        
        if psi < pi
            sinpsi = sqrt(1.0 - cospsi ^ 2);
        else 
            sinpsi = - sqrt(1.0 - cospsi ^ 2);
        end
        
        if 1 - abs(uz) <= One_Minus_Coszero
            uxx = sintheta * cospsi;
            uyy = sintheta * sinpsi;
            uzz = costheta * sign(uz);
        else
            temp = sqrt(1.0 - uz ^ 2);
            uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta;
            uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta;
            uzz = - sintheta * cospsi * temp + uz * costheta;
        end
        
        %{ update trajectory %}
        ux = uxx;
        uy = uyy;
        uz = uzz;
        
        if w < Threshold
            if rand <= Chance
                w = w / Chance;
            else 
                Photon_Status = false;
            end 
        end
     
    end
end

avg = sum / Count_Num
Count_Num

        
        
        

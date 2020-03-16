% addpath(genpath('/Users/jramella/MPL Dropbox/jr/A_RESEARCH/MONTECARLO_MATLAB_WORKS/MATLABMONTECARLO_BACKUP/MatScat'))
addpath(genpath('D:\=code=\Simulation_Intrinsic\JRcode\MatScat'))

format long

THRESHOLD   = 0.01
CHANCE      = 0.1
NN          = 100;
ALIVE      	= 1;
DEAD       	= 0;
MM          = NN-1

% CHOOSE MIE SCATTERING parameters */
rad         = 2/2; % microns */
lambda 		= 0.6328; % microns */
rho 		= 1.152e-4;%Dilution 1*/
Nphotons	= 1e4;
mua 		= 0.0; %µa  */

% ------------------------*/
nre_p   	= 1.59;
nim_p  	 	= 0;
nre_med		= 1.33;
nim_med 	= 0.0;
nangles 	= 1000;

m.re = nre_p/nre_med;
m.im = 0.0;
x    = 2*pi*rad/(lambda/nre_med);
vol  = 4.0/3*pi*rad*rad*rad;
A    = pi*rad*rad;
conv=1;

%% Calculate amplitude scattering matrix
[S, C, ang] = calcmie(rad, nre_p, nre_med, lambda, nangles, 'ConvergenceFactor', conv);
Q = getEfficiencies(C, rad(end), 3);

%% Calculate and plot Mueller matrix
M = getMuellerMatrix(S);
s11=squeeze(M(1,1,:));
s12=squeeze(M(1,2,:));
s33=squeeze(M(3,3,:));
s43=squeeze(M(4,3,:));

g       = 0.90809; %<<-- NEED TO SET THIS UP MANUALLY FOR NOW
mus 	= Q.ext*A*rho*1e4; % Mus is in cm^-1 */
musp 	= mus*(1-g);% (cm^-1) */
albedo 	= mus/(mus + mua);

hw 			= 7/mus; % (cm) , maximum range in x and y for output. */
dx 			= 2.0*hw/NN;
dy 			= 2.0*hw/NN;

slabsize	= 2;

IT=0;
QT=0;
UT=0;
VT=0;
IR_1=0;
QR_1=0;
UR_1=0;
VR_1=0;


for jjj = 1:4
    
    if (jjj == 1)
        
        S0=[1 1 0 0];
        display("launch H\n");
        
        TH=S0*0;
        RH=S0*0;
        
        IR=zeros(NN,NN);%W*/
        QR=zeros(NN,NN);
        UR=zeros(NN,NN);
        VR=zeros(NN,NN);
        
        
    end
    
    
    if jjj == 2
        
        S0=[1 -1 0 0];
        S=S0*0;
        RV=S0*0;
        TV=S0*0;
        IR=zeros(NN,NN);%W*/
        QR=zeros(NN,NN);
        UR=zeros(NN,NN);
        VR=zeros(NN,NN);
        
        display("launch V\n");
    end
    
    if jjj == 3
        
        S0=[1 0 1 0];
        TP=S0*0;
        RP=S0*0;
        IR=zeros(NN,NN);%W*/
        QR=zeros(NN,NN);
        UR=zeros(NN,NN);
        VR=zeros(NN,NN);
        
        display("launch P\n");
    end
    
    if jjj == 4
        
        S0=[1 0 0 1];
        TR=S0*0;
        RR=S0*0;
        IR=zeros(NN,NN);%W*/
        QR=zeros(NN,NN);
        UR=zeros(NN,NN);
        VR=zeros(NN,NN);
        
        display("launch R\n");
    end
    
    
    % LAUNCH photon */
    for i_photon = 1:Nphotons
        
        
        %pencil beam	*/
        x = 0.0;
        y = 0.0;
        z = 0.0;
        
        
        % photon direction cosines */
        
        U  =[0 0 1];
        S  = S0;
        S2 = S0*0;
        
        photon_status = ALIVE;
        W	= 1; % photon weight */
        cnt=0;
        %******** ALIVE cycle *****************/
        while (photon_status == ALIVE)
            
            %*** HOP */
            
            s = -log(rand)/(mus+mua);
            x = x+U(1)*s;
            y = y+U(2)*s;
            z = z+U(3)*s;
            
            %*** ABSORB */
            
            absorb = W*(1-albedo);
            W= W-absorb;
            cnt=cnt+1;
            
            
            if (z<=0)
                
                %return to detector reference frame*/
                phi=atan2(U(2),U(1));
                S2=rotSphi(S, phi);
                      
                if (x >= -hw)
                    ix = 1+round(abs(x + hw)/dx);
                else
                    ix=1;
                end
                
                
                if (y >= -hw)
                    iy = 1+round(abs(y + hw)/dy);
                else
                    iy=1;
                end
                
                if (ix > MM)
                    ix = MM;
                end
                
                if (iy > MM)
                    iy = MM;
                end
                
                
                if jjj==1
                    
                    RH=RH+S2*W;
                    IR(iy,ix) = IR(iy,ix) + W*S2(1);
                    
                    QR(iy,ix) = QR(iy,ix) + W*S2(2);
                    
                    UR(iy,ix) = UR(iy,ix) + W*S2(3);
                    
                    VR(iy,ix) = VR(iy,ix) + W*S2(4);
                    
                end
                
                if jjj==2
                    
                    RV=RV+S2*W;
                    
                    IR(iy,ix) = IR(iy,ix) + W*S2(1);
                    
                    QR(iy,ix) = QR(iy,ix) + W*S2(2);
                    
                    UR(iy,ix) = UR(iy,ix) + W*S2(3);
                    
                    VR(iy,ix) = VR(iy,ix) + W*S2(4);
                end
                
                if jjj==3
                    
                    RP=RP+S2*W;
                    IR(iy,ix) = IR(iy,ix) + W*S2(1);
                    
                    QR(iy,ix) = QR(iy,ix) + W*S2(2);
                    
                    UR(iy,ix) = UR(iy,ix) + W*S2(3);
                    
                    VR(iy,ix) = VR(iy,ix) + W*S2(4);
                end
                
                
                if jjj==4
                    
                    RR=RR+S2*W;
                    IR(iy,ix) = IR(iy,ix) + W*S2(1);
                    
                    QR(iy,ix) = QR(iy,ix) + W*S2(2);
                    
                    UR(iy,ix) = UR(iy,ix) + W*S2(3);
                    
                    VR(iy,ix) = VR(iy,ix) + W*S2(4);
                end
                photon_status = DEAD;
                
                
                
                
                
            elseif (z>=slabsize)
                phi = -atan2(U(2),U(1));
                
                S2  = rotSphi(S, phi);
             
                
                if jjj==1
                    
                    TH=TH+S2*W;
                end
                
                if jjj==2
                    
                    TV=TV+S2*W;
                end
                
                if jjj==3
                    
                    TP=TP+S2*W;
                    
                end
                
                
                if jjj==4
                    
                    TR=TR+S2*W;
                    
                end
                
                photon_status = DEAD;
                
            end%z>slab size*/
            
            
            % SPIN */
            
            % REJECTION METHOD to choose azimuthal angle phi and deflection angle theta */
            condition=1;
            
            while condition
                
                theta 	= acos(2*rand-1);
                
                phi = rand*2.0*pi;
                
                I0=s11(1)*S(1)+s12(1)*(S(2)*cos(2*phi)+S(3)*sin(2*phi));
                
                ithedeg = 1+floor(theta*nangles/pi);
                
                I=s11(ithedeg)*S(1)+s12(ithedeg)*(S(2)*cos(2*phi)+S(3)*sin(2*phi));
                
                xra=rand*I0;
                if (xra>=I)
                    condition=1;
                else
                    condition=0;
                    
                end
                                
            end
            
            
            %------------------------------------------------------------------------------
            %  Scattering : rotate to meridian plane	then scatter
            %------------------------------------------------------------------------------*/
            
            U2=updateU(U, phi, theta);  % update photon trajectory vector */
            
            costheta=cos(theta);
            
            S2=rotSphi(S, phi);
                       
            S(1)= s11(ithedeg)*S2(1)+s12(ithedeg)*S2(2);
            
            S(2)= s12(ithedeg)*S2(1)+s11(ithedeg)*S2(2);
            
            S(3)= s33(ithedeg)*S2(3)+s43(ithedeg)*S2(4);
            
            S(4)= -s43(ithedeg)*S2(3)+s33(ithedeg)*S2(4);
            
            
            temp=(sqrt(1-costheta*costheta)*sqrt(1-U2(3)*U2(3)));
            
            if ( temp==0)
                cosi=0;
            else
                
                if ((phi>pi) && (phi<2*pi))
                    cosi=(U2(3)*costheta-U(3))/temp;
                    
                else
                    cosi=-(U2(3)*costheta-U(3))/temp;
                    
                    if (cosi<=-1)
                        cosi=-1;
                    end
                    
                    if (cosi>=1)
                        cosi=1;
                    end
                end
            end
            %
            sini = real(sqrt(1-cosi*cosi));
            
            cos22=2*cosi*cosi-1;
            
            sin22=2*sini*cosi;
            
            S2(1)=S(1);
            S2(2)=S(2)*cos22-S(3)*sin22;
            S2(3)=S(2)*sin22+S(3)*cos22;
            S2(4)=S(4);
            
            
            S=S2/S2(1);
            
            
            U = U2; % update U */
            
            %ROULETTE*/
            
            if (W<THRESHOLD)
                if (rand<=CHANCE)
                    W=W/CHANCE;
                else
                    photon_status=DEAD;
                end
            end
            
            
        end % end of single photon launching */
        
    end% slab size*/
    
    if jjj==1
        
        display(RH)
        display(TH)
        
        outHI=IR;
        outHQ=QR;
        outHU=UR;
        outHV=VR;
        
        save('outHI','outHI');
        save('outHQ','outHQ');
        save('outHU','outHU');
        save('outHV','outHV');
        
        clear outHI outHQ outHV outHU
        
        
    end
    
    if jjj==2
        
        display(RV)
        display(TV)
        
        outVI=IR;
        outVQ=QR;
        outVU=UR;
        outVV=VR;
        
        save('outVI','outVI');
        save('outVQ','outVQ');
        save('outVU','outVU');
        save('outVV','outVV');
        
        clear outVI outVQ outVV outVU
        
        
    end
    
    if jjj==3
        
        display(RP)
        display(TP)
        
        
        outPI=IR;
        outPQ=QR;
        outPU=UR;
        outPV=VR;
        
        save('outPI','outPI');
        save('outPQ','outPQ');
        save('outPU','outPU');
        save('outPV','outPV');
        
        clear outPI outPQ outPV outPU
        
        
    end
    
    
    if jjj==4
        
        display(RR/Nphotons)
        display(TR/Nphotons)
        
        outRI=IR;
        outRQ=QR;
        outRU=UR;
        outRV=VR;
        
        save('outRI','outRI');
        save('outRQ','outRQ');
        save('outRU','outRU');
        save('outRV','outRV');
        
        clear outRI outRQ outRV outRU
        
        
    end
    
end

toc

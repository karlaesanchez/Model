
function [rv_A,r_cap,rv_V,radnr_A,radnr_V,hv_A,hv_V,hnr_A,hnr_V,Betav_A,Betav_V,Betavnr_A,Betavnr_V,N,Axa_A,Axg_C,Axa_V,Axv_A,Axv_V,Axg_A,Axg_V,AxgnrA,AxgnrV,Asg_A,Asg_V,Asnr_A,twopowerarray,Asnr_V,lv_A,lv_V,Vg_A,Vg_V,Rgnr_A,R0_C,Rgnr_V,Rg_A,Rg_V,pmmHgA,Poutnr_C,pmmHgV,Ping_A,Pout_C,Poutg_V,Sgnr_A,Sgnr_C,Sgnr_V,Sg_A,Sg_C,Sg_V]=GenM_L_101018(gamma, lambda, Pe, phbrain, Pica, Qvar)


%% Fixed input
Qdev     = 750*(1e-6/60);                     % [m^3/s]   Developmental Input Cerebral blood flow i.e. for estimating generation number. Value = 750 ml/min from Oktar et al. (2007). Based on L and R ICAs.
Peref    = 11*133-phbrain;                    % [Pa]      Reference external pressure.
vcutoff  = 4;                                 % [non-dim] Vein generation that goes to sinuses.
gammadev = 0.5547;%  0.5507;  %                        % [non-dim] Scaling factor of cross-sectional area (CSA). 
xi       = lambda/(2*gamma^2);                % [non-dim] Scaling factor of generational resistance. 
r_cap    = 3e-6;                              % [m]       Capillary radius.
u_cap    = 0.001;                             % [m/s]     Velocity in capillary.
Qcap     = u_cap*pi*r_cap^2;
nocap    = Qdev/Qcap;                         % [non-dim] Number of capillaries.
nogen    = floor(log2(nocap));                % [non-dim] Rounding. 
N        = nogen;                             % [non-dim] Number of generations.
ica_fac  = 0.75;% 0.73; %                            % [non-dim] Root vessel ratio.
va_ratio = 1.75;                              % [non-dim] Vein-artery ratio.
r0_A     = ica_fac*(r_cap*gammadev^(-(N)/2)); % [m]       Input arterial radius. Scaled by assuming two ICAs + vertebral arteries.
h0_A     = (r0_A*2)*0.1;                      % [m]       Input arterial wall thickness. Caro&Parker (taken from h/d ratio).

evar = 1;
kvar = 1;
lvar = 1;
R0var= 1;

Nac      = 8;                                 % Number of capillaries per arteriole
l0_A     = 0.085;%0.075*lvar;                        % [m]       Input arterial length. Caro&Parker. Estimated from ~half of CA (.15)
E_A      = 9e5*evar;                          % [N/m^2]   Young's Modulus. Pending ref. One order of magnitude bigger than Caro's. 
mu       = 0.0033;                            % [Pa.s]    Viscosity of blood. Guyton (Kestin et al., & Guyton), three times viscosity of water.
mu_csf   = 0.00085;                        % [Pa.s]    Viscosity of cerebrospinal fluid (CSF). Bloomfield (1998): Effects of Proteins, Blood Cells and Glucose on the Viscosity of Cerebrospinal Fluid.Similar to water. This since the only possible filtration solute would be CSF-like fluid
k_A      = kvar*1e-21;%3.375e-22*kvar;%1.5e-21*kvar;                      % [m^2]     Arterial compartmental permeability. Used in Lakin's et al., not clear how they estimated it. %N.B. Matlab's machine accuracy goes to e-7 digits. Check permeabilities.
R0_C     = (10/750)*(133*60*1e6);%3.5467e07;%(10/750)*(133*60*1e6);%3.5467e07*R0var;%%R0var*(15/750)*(133*60*1e6);%0.0150*(133*60*1e6);%(15/750)*(133*60*1e6);%106400000;%%3*1.596e7;%13*(133*60/1e-6);%3.5467e07*R0var;%                         % [Pa/m^3/s]Baseline resistance of capillaries =16mmHg/ml/min taken from Linninger et al. 2009.
r0_V     = va_ratio*r0_A;                     % [m]       Output venous radious. Changed from *0.70 as 70% being the vein of Galen.
h0_V     = h0_A/2;                            % [m]       Output venous wall thickness. Caro&Parker.
l0_V     = l0_A;                              % [m]       Output venous length. Caro&Parker. Estimated from ~half of CA (.15)
E_V      = 0.7e5*evar;                        % [N/m^2]   Young's Modulus. Pending ref. One order of magnitude bigger than Caro's. 
k_V      = k_A/2;                             % [m^2]     Venous compartmental permeability. Lakin's %1.d-12; %m^2 Arbitrary at the moment
kc       = kvar*1.8797e-11; %Guyton
lc       = 6e-4;%lv_A(N);%

%% Generation 0 and other non-loop dependent parameters

Ax0_A = pi*(r0_A^2);
As0_A = 2*pi*r0_A*l0_A;
R0_A  = (8*mu*l0_A)/(pi*r0_A^4);
V0_A  = l0_A*Ax0_A;
R0_V  = (8*mu*l0_V)/(pi*r0_V^4);

% Useful arrays
twopowerarray=zeros(N,1);
twolambdapowerarray=zeros(N,1);
    for i=0:N-1
        j=i+1;
        twopowerarray(j)=2.^i;
        twolambdapowerarray(j)=(2*lambda).^i;
    end

    
%% Arteries

%Preallocating arrays
lv_A     = zeros(N,1);    % Arterial length per vessel (and generation).
rv_A     = zeros(N,1);    % Arterial radius.
Rv_A     = zeros(N,1);    % Arterial resistance per vessel.
Paveg_A  = zeros(N,1);    % Arterial average pressure per generation.
Ping_A   = zeros(N,1);    % Arterial input pressure per generation.
Poutg_A  = zeros(N,1);    % Arterial output pressure per generation.

for n = 0:N-1
    g = n+1;
% Individual vessel
    rv_A(g)    = r0_A*sqrt(gamma^n);
    lv_A(g)    = l0_A*(lambda^n);
    Rv_A(g)    = R0_A*(lambda/gamma^2)^n;
    
% Generational
    Ping_A(g)  = Pica-Qdev*R0_A*(xi^n-1)/(xi-1);
    Poutg_A(g) = Ping_A(g)-Qdev*R0_A*xi^n;
    Paveg_A(g) = (Ping_A(g)+Poutg_A(g))/2;
end

Axv_A    = pi*rv_A.^2;
Asv_A    = 2*pi*rv_A.*lv_A;
Vv_A     = Axv_A.*lv_A;
Axg_A    = twopowerarray.*Axv_A;
Asg_A    = twopowerarray.*Asv_A;
Vg_A     = twopowerarray.*Vv_A;
Rg_A     = Rv_A./twopowerarray;
hv_A     = (h0_A/(2*r0_A))*2.*rv_A; 
Betav_A  = (4/3)*sqrt(pi)*E_A*hv_A;
Ptrans_A = Paveg_A-Peref;
Sg_A     = (k_A.*Asg_A.*Ptrans_A)./(mu_csf.*hv_A);
    

%% Capillaries

Pin_C    = Poutg_A(end);       % Capillaries input pressure.
Pout_C   = Pin_C-(R0_C*Qdev);  % Capillaries output pressure.
Paveg_C  = (Pin_C+Pout_C)/2;
Sg_C     = kc*(Paveg_C-Peref);


%% Veins

% Preallocating arrays
lv_V     = zeros(N,1);    % Venous length per vessel (and generation).
rv_V     = zeros(N,1);    % Venous radius.
Rv_V     = zeros(N,1);    % Venous resistance per vessel.
Paveg_V  = zeros(N,1);    % Venous average pressure per generation.
Ping_V   = zeros(N,1);    % Venous input pressure per generation.
Poutg_V  = zeros(N,1);    % Venous output pressure per generation.

for m = 0:N-1
    h = m+1;   
% Individual vessel
    rv_V(h)    = r0_V*sqrt(gamma^(N-1-m)); 
    lv_V(h)    = l0_V*lambda^(N-1-m);
    Rv_V(h)    = R0_V*(lambda/gamma^2)^(N-1-m);
    
% Generational
    Ping_V(h)  = Pout_C-Qdev*R0_V*xi^(N-1)*(xi^(-m)-1)/(xi^(-1)-1);
    Poutg_V(h) = Pout_C-Qdev*R0_V*xi^(N-1)*(xi^(-m-1)-1)/(xi^(-1)-1);       
    Paveg_V(h) = (Ping_V(h)+Poutg_V(h))/2;
end

Axv_V    = pi*rv_V.^2;
Asv_V    = 2*pi*rv_V.*lv_V;
Vv_V     = Axv_V.*lv_V;
Axg_V    = flipud(twopowerarray).*Axv_V;
Asg_V    = flipud(twopowerarray).*Asv_V;
Vg_V     = flipud(twopowerarray).*Vv_V;
Rg_V     = Rv_V./flipud(twopowerarray);
hv_V     = (h0_V/(2*r0_V))*2.*rv_V;
Betav_V  = (4/3)*sqrt(pi)*E_V*hv_V;
Ptrans_V = Paveg_V-Peref;
Sg_V     = (k_V.*Asg_V.*Ptrans_V)./(mu_csf.*hv_V);


%% Compliant model: Arteries & Veins

Omeg_A   = (4.*sqrt(Axv_A).*Rg_A.*Qvar)./Betav_A;
Omeg_V   = (4.*sqrt(Axv_V).*Rg_V.*Qvar)./Betav_V;

% M matrix (x4)       
dm       = ones(N,1);
diag_M1A = diag(-1*dm, 0) + diag(dm(1:N-1), -1);
diag_M2A = diag(dm.*Omeg_A, 0);   
diag_M3A = diag(-(1/2)*dm, 0) + diag(-(1/2)*dm(1:N-1), -1);
diag_M4A = diag(1*dm, 0);
MA       = [diag_M1A diag_M2A; diag_M3A diag_M4A];

diag_M1V = diag(-1*dm, 0) + diag(dm(1:N-1), -1);
diag_M2V = diag(dm.*Omeg_V, 0); 
diag_M3V = diag(-(1/2)*dm, 0) + diag(-(1/2)*dm(1:N-1), -1);
diag_M4V = diag(1*dm, 0);
MV       = [diag_M1V diag_M2V; diag_M3V diag_M4V];


% Vectors
% b1A      = [Qvar.*(Rg_A.*(1+4*(sqrt(Axv_A)./Betav_A).*(Paveg_A-Peref+Pe)));zeros(N,1)];
% b1A      = [Qvar.*(Rg_A.*(1+4*(sqrt(Axv_A)./Betav_A).*(Pe)));zeros(N,1)];

% b1A      = [Qvar.*Rg_A.*(1+4.*(sqrt(Axv_A).*(Ptrans_A+Pe)./Betav_A));zeros(N,1)];

% b1A      = [Qvar.*Rg_A.*(1+4.*(sqrt(Axv_A).*Pe./Betav_A));zeros(N,1)];


Ohm_A    = 1+(4.*sqrt(Axv_A).*Pe)./Betav_A;
b1A      = [Rg_A.*Qvar.*Ohm_A;zeros(N,1)];

b1A(1)   = b1A(1)-Pica;
b1A(N+1) = b1A(N+1)+Pica/2;

pA       = MA\b1A; 
pmmHgA   = pA./133;

Pinnr_C  = pA(N);   
Poutnr_C = Pinnr_C-(R0_C*Qvar);

% b1V      = [Qvar.*(Rg_V.*(1+4*(sqrt(Axv_V)./Betav_V).*(Paveg_V-Peref+Pe)));zeros(N,1)];
% b1V      = [Qvar.*(Rg_V.*(1+4*(sqrt(Axv_V)./Betav_V).*(Pe)));zeros(N,1)];

% b1V      = [Qvar.*Rg_V.*(1+4.*(sqrt(Axv_V).*(-Ptrans_V+Pe)./Betav_V));zeros(N,1)];
% b1V      = [Qvar.*Rg_V.*(1+4.*(sqrt(Axv_V).*(Ptrans_V+Pe)./Betav_V));zeros(N,1)];
% b1V      = [Qvar.*Rg_V.*(1+4.*(sqrt(Axv_V).*Pe./Betav_V));zeros(N,1)];

Ohm_V    = 1+(4.*sqrt(Axv_V).*Pe)./Betav_V;
b1V      = [Rg_V.*Qvar.*Ohm_V;zeros(N,1)];


b1V(1)   = b1V(1)-Poutnr_C;
b1V(N+1) = b1V(N+1)+Poutnr_C/2;

pV       = MV\b1V; 
pmmHgV   = pV./133;

% % New transmural flux with actual CSA (non-rigid) and resistances

Ptransnr_A = pA(N+1:2*N)-Pe;
% Axa_A      = Axv_A.*(1+((pA(N+1:2*N)-Pe-Paveg_A+Peref).*sqrt(Axv_A))./Betav_A).^2;
% Axa_A      = Axv_A.*(1+(Ptransnr_A-Ptrans_A).*sqrt(Axv_A)./Betav_A).^2;

% Axa_A      = Axv_A.*(1+((pA(N+1:2*N)-Pe).*sqrt(Axv_A))./Betav_A).^2;

Omeptr_A  = Ptransnr_A.*sqrt(Axv_A)./Betav_A;
Axa_A      = Axv_A.*(1+Omeptr_A).^2;


Asnr_A     = 2.*sqrt(pi.*Axa_A).*lv_A;
rnr_A      = (Axa_A./pi).^(1/2);
hnr_A      = (2.*rnr_A).*0.1;
Sgnr_A     = (k_A.*Asnr_A.*twopowerarray.*Ptransnr_A)./(mu_csf.*hnr_A); 


Pavegnr_C  = (Pinnr_C+Poutnr_C)/2;
Ptransnr_C = Pavegnr_C-Pe;
Sgnr_C     = kc*(Pavegnr_C-Pe);

Ptransnr_V = pV(N+1:2*N)-Pe;
% Axa_V      = Axv_V.*(1+((pV(N+1:2*N)-Pe-Paveg_V+Peref).*sqrt(Axv_V))./Betav_V).^2;
% Axa_V      = Axv_V.*(1+(Ptransnr_V-Ptrans_V).*sqrt(Axv_V)./Betav_V).^2;
% Axa_V      = Axv_V.*(1+((pV(N+1:2*N)-Pe).*sqrt(Axv_V))./Betav_V).^2;

Omeptr_V   = Ptransnr_V.*sqrt(Axv_V)./Betav_V;
Axa_V      = Axv_V.*(1+Omeptr_V).^2;

Axg_C = 1;
rnr_V      = (Axa_V./pi).^(1/2);
hnr_V      = (2.*rnr_V).*0.1;
Asnr_V     = 2.*sqrt(pi.*Axa_V).*lv_V; 
Sgnr_V     = (k_V.*Asnr_V.*flipud(twopowerarray).*Ptransnr_V)./(mu_csf.*hnr_V);

Rgnr_A = (8*pi*mu*lv_A)./(twopowerarray.*Axa_A.^2);
Rgnr_V = (8*pi*mu*lv_V)./(flipud(twopowerarray).*Axa_V.^2);

% Other parameters for SA purposes
radnr_A    = sqrt(Axa_A/pi);
Betavnr_A  = (4/3)*sqrt(pi)*E_A*hnr_A;

radnr_V    = sqrt(Axa_V/pi);
Betavnr_V  = (4/3)*sqrt(pi)*E_V*hnr_V;


%% Gen function output

Pout=pV(N+1-vcutoff); 
Va=sum([Axa_A.*lv_A.*twopowerarray;Axa_V.*lv_V.*flipud(twopowerarray)]);%+Vc;
Var=sum([Axv_A.*lv_A.*twopowerarray;Axv_V.*lv_V.*flipud(twopowerarray)]);%+Vc;
Sa=sum([Sgnr_A;Sgnr_V])+Sgnr_C;

% Other outputs 
VnrA= lv_A.*Axa_A.*twopowerarray;
VnrV= lv_V.*Axa_V.*flipud(twopowerarray);
AxgnrA= Axa_A.*twopowerarray;
AxgnrV= Axa_V.*flipud(twopowerarray);

V_al = sum(Axa_A.*lv_A.*twopowerarray);
V_vl = sum(Axa_A.*lv_A.*twopowerarray);
As_av= sum([Asnr_A.*twopowerarray;Asnr_V.*flipud(twopowerarray)]);
As_al= sum(Asnr_A.*twopowerarray);
As_vl= sum(Asnr_A.*twopowerarray);

Ptdiff_A = Ptransnr_A-Ptrans_A;
Ptdiff_V = Ptransnr_V-Ptrans_V;



%% Exporting values

% assignin('base','Asg_C',Asg_C)
assignin('base','Ptdiff_A',Ptdiff_A/133)
assignin('base','Ptdiff_V',Ptdiff_V/133)
assignin('base','V_al',V_al)
assignin('base','V_vl',V_vl)
assignin('base','As_av',As_av)
assignin('base','As_al',As_al)
assignin('base','As_vl',As_vl)
assignin('base','Ptransnr_A',Ptransnr_A)
assignin('base','Ptransnr_C',Ptransnr_C)
assignin('base','Ptransnr_V',Ptransnr_V)
assignin('base','Ptnrmmhg_A',Ptransnr_A/133)
assignin('base','Ptnrmmhg_V',Ptransnr_V/133)
assignin('base','r_cap',r_cap)
assignin('base','twopowerarray',twopowerarray)
assignin('base','Rgnr_A',Rgnr_A)
assignin('base','Rgnr_V',Rgnr_V)
assignin('base','R0_C',R0_C)
assignin('base','Pica',Pica/133)
assignin('base','pmmHgA',pmmHgA)
assignin('base','pmmHgV',pmmHgV)
assignin('base','Pout_C',Pout_C)
assignin('base','Sgnr_C',Sgnr_C)
assignin('base','VnrA',VnrA)
assignin('base','VnrV',VnrV)
assignin('base','AxgnrA',AxgnrA)
assignin('base','AxgnrV',AxgnrV)
assignin('base','radnr_A',radnr_A)
assignin('base','radnr_V',radnr_V)
assignin('base','hnr_A',hnr_A)
assignin('base','hnr_V',hnr_V)
assignin('base','Betavnr_A',Betavnr_A)
assignin('base','Betavnr_V',Betavnr_V)
assignin('base','Asg_A',Asg_A)
assignin('base','Asg_V',Asg_V)
assignin('base','Asnr_A',Asnr_A)
assignin('base','Asnr_V',Asnr_V)
assignin('base','Asv_A',Asv_A)
assignin('base','Asv_V',Asv_V)
assignin('base','Axa_A',Axa_A)
assignin('base','Axa_V',Axa_V)
assignin('base','Axagen_A',Axa_A.*twopowerarray)
assignin('base','Axagen_V',Axa_V.*flipud(twopowerarray))
assignin('base','Axg_A',Axg_A)
assignin('base','Axg_V',Axg_V)
assignin('base','Axg_C',Axg_C)
assignin('base','Axv_A',Axv_A)
assignin('base','Axv_V',Axv_V)
assignin('base','b1A',b1A)
assignin('base','b1V',b1V)
assignin('base','Betav_A',Betav_A)
assignin('base','Betav_V',Betav_V)
assignin('base','h0_A',h0_A)
assignin('base','h0_V',h0_V)
assignin('base','hv_A',hv_A)
assignin('base','hv_V',hv_V)
assignin('base','l0_A',l0_A)
assignin('base','l0_V',l0_V)
assignin('base','lv_A',lv_A)
assignin('base','lv_V',lv_V)
assignin('base','MV',MV)
assignin('base','N',N)
assignin('base','nocap',nocap)
assignin('base','pA',pA/133)
assignin('base','Paveg_A',Paveg_A/133)
assignin('base','Paveg_V',Paveg_V/133)
assignin('base','Ping_A',Ping_A/133)
assignin('base','Ping_V',Ping_V/133)
assignin('base','Pout',Pout/133)
assignin('base','Poutnr_C',Poutnr_C/133)
assignin('base','Poutg_A',Poutg_A/133)
assignin('base','Poutg_V',Poutg_V/133)
assignin('base','Ptrans_A',Ptrans_A/133)
assignin('base','Ptrans_V',Ptrans_V/133)
assignin('base','pV',pV/133)
assignin('base','r0_A',r0_A)
assignin('base','r0_V',r0_V)
assignin('base','R0_A',R0_A)
assignin('base','R0_V',R0_V)
% assignin('base','Rcon_A',Rcon_A)
% assignin('base','Rcon_V',Rcon_V)
assignin('base','Rg_A',Rg_A*(1e-6/60)/133)
assignin('base','Rg_V',Rg_V)
assignin('base','rv_A',rv_A)
assignin('base','rv_V',rv_V)
assignin('base','S_tot',Sa/(1e-6/60))
assignin('base','Sg_A',Sg_A)
assignin('base','Sg_V',Sg_V)
assignin('base','Sgnr_A',Sgnr_A/(1e-6/60))
assignin('base','Sgnr_V',Sgnr_V/(1e-6/60))
assignin('base','Va',Va)
assignin('base','Var',Var)
assignin('base','Vg_A',Vg_A)
assignin('base','Vg_V',Vg_V)
assignin('base','Vv_A',Vv_A)
assignin('base','Vv_V',Vv_V)


end



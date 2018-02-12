% clear all
close all
format long
clc

global gv

tic

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
disp('Pre-processing')

%-------------------------------------------------------------------------
% beginning of input section

% Geometrical input del mantello siliconico
Es          =40e3;                          % [Pa] modulo di young
nus         =500;                           % [Pa*s] viscosità
Pois        =0;                             % [-] modulo di Poisson (incomprimibile 0.5)
Gs          =Es/(2*(1+Pois));               % [Pa] modulo di taglio
h           =1e-3;                          % [m] half thickness
Rn          =12.5e-3;                       % [m] Raggio Nuzzle
Rb          =10e-3;                         % [m] Raggio bocche di entrata 
Area_out    =pi*Rn.^2;                      % [m^2] area di uscita
Area_in     =Area_out+3*pi*Rb.^2;           % [m^2] area di entrata
Js          =h.^2/3;                        % [m^2]
Is          =eps;
% posizione relativa a corpo rigido
gsr         =[0 0 1 0;0 1 0 0; -1 0 0 0; 0 0 0 1];

% Geometrical input del braccio siliconico
Eb          =110e3;                         % [Pa] modulo di young 110e3
etab        =300;                           % [Pa*s] viscosità
Poib        =0;                             % [-] modulo di Poisson
Gb          =Eb/(2*(1+Poib));               % [Pa] modulo di taglio
Rmax        =10e-3;                         % [m] Raggio di base
Rmin        =2e-3;                          % [m] Raggio finale
Lb          =245e-3;                        % [m] Lunghezza del braccio
Lr          =112e-3;                        % [m] lunghezza corpo rigido
nsezb       =floor(Lb*2e2+1);               % una sezione per mezzo centimetro
Xb          =linspace(0,Lb,nsezb);          % [m] curvilinear abscissa
Rb          =((Rmin-Rmax)/Lb)*Xb + Rmax;    % [m] Raggio della sezione attuale
Rb_primo    =(Rmin-Rmax)/Lb;                % [-] derivata del raggio
A           =pi*Rb.^2;                      % [m^2]
A_primo     =2*pi*Rb_primo*Rb;              % [m]
Jb          =pi*Rb.^4/4;                    % [m^4]
Jb_primo    =pi*Rb_primo*Rb.^3;             % [m^3]
Ib          =pi*Rb.^4/2;                    % [m^4]
Ib_primo    =2*pi*Rb_primo*Rb.^3;           % [m^3]
% posizione relativa a corpo rigido
gbr         =gsr*[1 0 0 0; 0 1 0 0; 0 0 1 -Lr; 0 0 0 1];

% Geometrical input del corpo rigido (cono)
hr          =20e-3;                         % [m] raggio base
Volr        =pi*hr^2*Lr/3;                  % [m^3] volume corpo rigido
Jr          =(3/5)*(Rmax^2/4+Lr^2);         % [m^2]
Ir          =(3/10)*Rmax^2;                 % [m^2]
% posizione del cono relativa alla punta
grr         =gsr*[1 0 0 0; 0 1 0 0; 0 0 1 -Lr*3/4; 0 0 0 1];

%-------------------------------------------------------------------------
% initial configuration t=0

% shell
Rs                       =31e-3;                % [m] Raggio parte sferica 31e-3
alpha0                   =pi/4;                 % [rad] buco anteriore pi/4 - pi/6 - pi/9 - pi/18 - pi/36
LH                       =(pi/2-alpha0)*Rs;     % [m] lunghezza parte sferica
sezH                     =20;                   % [-] sezioni parte sferica 20 - 15 - 17 - 19 - 20
dXs                      =LH/sezH;              % [m] delta X
sezS                     =floor(3.9*sezH);      % [-] sezioni parte cilindrica
LS                       =sezS*dXs;             % [m] lunghezza parte cilindrica
nsezs                    =sezS+sezH+1;          % [-] sezioni totali
Ls                       =LS+LH;                % [m] Lunghezza del profilo
Xs                       =linspace(0,Ls,nsezs); % [m] curvilinear abscissa
erre0                    =zeros(1,nsezs);
erre0_primo              =zeros(1,nsezs);
erre0_sec                =zeros(1,nsezs);
erre0_ter                =zeros(1,nsezs);
alpha                    =Xs(1:sezH+1)/Rs+alpha0;
erre0(1:sezH+1)          =Rs*sin(alpha);
erre0(sezH+2:nsezs)      =Rs*ones(1,sezS);      % [m] Raggio iniziale
erre0_primo(1:sezH+1)    =cos(alpha);
erre0_primo(sezH+2:nsezs)=zeros(1,sezS);
erre0_sec(1:sezH)        =-sin(alpha(1:sezH))/Rs;
erre0_sec(sezH+1)        =diff([erre0_primo(sezH) erre0_primo(sezH+2)],1,2)./(2*dXs);
erre0_sec(sezH+2:nsezs)  =zeros(1,sezS);
erre0_ter(1:sezH)        =-cos(alpha(1:sezH))/(Rs^2);
erre0_ter(sezH+1)        =diff([erre0_sec(sezH) erre0_sec(sezH+2)],1,2)./(2*dXs);
erre0_ter(sezH+2:nsezs)  =zeros(1,sezS);
zeta0_primo              =-sqrt(ones(1,nsezs)-erre0_primo.^2);
zeta0_sec                =-erre0_primo.*erre0_sec./zeta0_primo;
zeta0                    =dXs*cumtrapz(zeta0_primo);        % [m] altitudine iniziale
senth0                   =zeta0_primo;                      % seno theta iniziale
costh0                   =erre0_primo;                      % coseno theta iniziale
mu0                      =-erre0_sec./zeta0_primo;          % theta primo iniziale
varsigma0                =senth0./erre0;
senth0_primo             =mu0.*costh0;
varsigma0_pr             =senth0_primo./erre0-varsigma0.*(erre0_primo./erre0);
% derivate rispetto a X del 2° ordine centrata dentro e upwind ai bordi 1°
mu0_primo                =zeros(1,nsezs);
mu0_primo(:,1)           =diff([mu0(:,1) mu0(:,2)],1,2)./dXs;
mu0_primo(2:2:nsezs-1)   =diff(mu0(1:2:nsezs),1,2)./(2*dXs);
mu0_primo(3:2:nsezs-1)   =diff(mu0(2:2:nsezs),1,2)./(2*dXs);
mu0_primo(nsezs)         =diff([mu0(nsezs-1) mu0(nsezs)],1,2)./dXs;

theta0                   =zeros(1,nsezs);
for ii=1:nsezs
    if (costh0(ii) >= 0)
        theta0(ii)             =asin(senth0(ii));
    else
        theta0(ii)             =(sign(senth0(ii)))*pi-asin(senth0(ii));
    end
end
VolX                     =pi*dXs*cumtrapz(erre0.^2);
Vol                      =VolX(nsezs)-VolX(1);
AreamX                   =2*pi*dXs*cumtrapz(erre0);
Aream                    =AreamX(nsezs)-AreamX(1);
VolmX                    =2*h*AreamX;
Volm                     =VolmX(nsezs)-VolmX(1);

% beam
xcib_star    =[0 0 0 1 0 0]'*ones(1,nsezb);

%-------------------------------------------------------------------------
% actuation load

% shell
per         =0.66;                   % [s] scegli fra 0.53; 0.66; 0.79
ten         =-15*triangularPulse(0,1/2,1,[0 1/2 1]); % [N/m] internal actuation force on the circuferential direction -10

% beam
tactb       =1.0;                    % [s] tempo torque in dir z o y
trelb       =1.5;                    % [s] tempo rilassamento
ntrelb      =trelb*10^2;
deltam      =0.005;                  % [Nm] max internal actuation torque 0.005, 0.007

%-------------------------------------------------------------------------
% dinamic parameters
% Shell

ro_water    =1022;                        % [Kg/m^3] densità acqua
massas      =336e-3;                      % [Kg] massa shell piu` motori
ro_shel     =massas/Volm;                 % [Kg/m^3] densità shell 1022
B           =1.1;                         % [-] coeff. massa aggiunta 1.1
Cds         =1.7;                         % [-] coeff. viscosità 1.7
Cf          =1;                           % [-] coeff. perdita nuzzle 1
Blbs        =0;                           % [-] coeff. massa aggiunta perpendicolare 0
Epsi        =(2*Es*h)/(1-Pois^2);         % [N/m]
Epsi2       =(6*nus*h)/(1-Pois^2);
Gammas      =ro_shel*diag([Js Is Js 2*h 2*h 2*h]);
Gammas_ag   =Gammas+ro_shel*B*diag([0 0 0 2*h 2*h 2*h]);
Gammas_ag2  =Gammas_ag+diag([0 0 0 0 ro_water*2*h*Blbs 0]);

% Beam
massab      =34e-3;                       % [Kg] massa braccio
ro_arm      =1080;                        % [Kg/m^3] densità nominale 1080
gra         =[0;0;0];                 % [m/s^2] vettore gravità [0;0;-9.81]
Clt         =0.01;                        % [-] coeff. viscosità longitudinale
Cln         =2.5;                         % [-] coeff. viscosità trasversale
Clb         =2.5;                         % [-] coeff. viscosità trasversale
Bln         =1.5;                         % [-] coeff. massa aggiunta trasversale
Blbb        =1.5;                         % [-] coeff. massa aggiunta trasversale

Db          =zeros(3,3*nsezb);
Eps         =zeros(6,6*nsezb);
Eps_primo   =zeros(6,6*nsezb);
Ipsi        =zeros(6,6*nsezb);
Ipsi_primo  =zeros(6,6*nsezb);
Gammab      =zeros(6,6*nsezb);
Ag          =zeros(6,6*nsezb);
Gammab_ag   =zeros(6,6*nsezb);
for ii=1:nsezb
    Db(:,3*(ii-1)+1:3*(ii-1)+3)          =diag([0.5*pi*Rb(ii)*Clt Rb(ii)*Cln Rb(ii)*Clb]);                               % drag coef matrix
    Eps(:,6*(ii-1)+1:6*(ii-1)+6)         =diag([Gb*Ib(ii) Eb*Jb(ii) Eb*Jb(ii) Eb*A(ii) Gb*A(ii) Gb*A(ii)]);              % stifness matrix
    Ipsi(:,6*(ii-1)+1:6*(ii-1)+6)        =etab*diag([Ib(ii) 3*Jb(ii) 3*Jb(ii) 3*A(ii) A(ii) A(ii)]);                     % viscosity matrix
    Eps_primo(:,6*(ii-1)+1:6*(ii-1)+6)   =diag([Gb*Ib_primo(ii) Eb*Jb_primo(ii) Eb*Jb_primo(ii) Eb*A_primo(ii) Gb*A_primo(ii) Gb*A_primo(ii)]);
    Ipsi_primo(:,6*(ii-1)+1:6*(ii-1)+6)  =etab*diag([Ib_primo(ii) 3*Jb_primo(ii) 3*Jb_primo(ii) 3*A_primo(ii) A_primo(ii) A_primo(ii)]);
    Gammab(:,6*(ii-1)+1:6*(ii-1)+6)      =ro_arm*diag([Ib(ii) Jb(ii) Jb(ii) A(ii) A(ii) A(ii)]);                         % inertia matrix
    Ag(:,6*(ii-1)+1:6*(ii-1)+6)          =ro_water*diag([0 0 0 0 A(ii)*Bln A(ii)*Blbb]);                                 % mass to add
    Gammab_ag(:,6*(ii-1)+1:6*(ii-1)+6)   =Gammab(:,6*(ii-1)+1:6*(ii-1)+6)+Ag(:,6*(ii-1)+1:6*(ii-1)+6);                   % added mass matrix
end

% Rigid Body
Crx                 =21e-3;                % [m^2] coeff. viscosità
Brx                 =0;                    % [-] coeff. massa aggiunta
massar              =204e-3;               % [kg] massa corpo rigido
ro_root             =massar/Volr;          % [Kg/m^3] densità corpo rigido
Gammarl             =massar*diag([Jr Jr Ir 1 1 1]);
Gammarl_ag          =Gammarl+massar*Brx*diag([0 0 0 1 1 1]);
Gammar              =dinamico_coAdjoint(grr)*Gammarl*dinamico_Adjoint(grr^-1);
Gammar_ag           =dinamico_coAdjoint(grr)*Gammarl_ag*dinamico_Adjoint(grr^-1);

%-------------------------------------------------------------------------
% numerical setting
time        =5;                      % [s]
nsol        =time*10^2+1;            % una soluzione ogni centisecondo
tspan       =linspace(0,time,nsol);  % [s] time
dXb         =Lb/(nsezb-1);           % delta X beam
nspi        =2;                      % ATTENZIONE deve essere pari
dphi        =2*pi/nspi;              % delta phi
phi         =linspace(0,(2*pi-dphi),nspi);  % [rad] curvilinear abscissa 2
cosphi      =cos(phi);
senphi      =sin(phi);

%-------------------------------------------------------------------------
% osservabili

F_bire2     =zeros(1,nsol);
F_bire3     =zeros(1,nsol);
F_bere2     =zeros(1,nsol);
F_bere3     =zeros(1,nsol);
F_biare2    =zeros(1,nsol);
F_biare3    =zeros(1,nsol);
pitch       =zeros(1,nsol);
axes        =zeros(3,nsol);
magn        =zeros(1,nsol);
Volu        =zeros(1,nsol);
nstep       =1;


% global variable shell
gv.Es          =Es;
gv.nus         =nus;
gv.Gs          =Gs;
gv.Pois        =Pois;
gv.Ls          =Ls;
gv.gsr         =gsr;
gv.Area_out    =Area_out;
gv.Area_in     =Area_in;
gv.Xs          =Xs;
gv.phi         =phi;
gv.cosphi      =cosphi;
gv.senphi      =senphi;
gv.tspan       =tspan;
gv.h           =h;
gv.per         =per;
gv.ten         =ten;
gv.erre0       =erre0;
gv.erre0_primo =erre0_primo;
gv.senth0      =senth0;
gv.costh0      =costh0;
gv.theta0      =theta0;
gv.mu0         =mu0;
gv.varsigma0   =varsigma0;
gv.varsigma0_pr=varsigma0_pr;
gv.mu0_primo   =mu0_primo;
gv.Aream       =Aream;
gv.Js          =Js;
gv.Epsi        =Epsi;
gv.Epsi2       =Epsi2;
gv.ro_water    =ro_water;
gv.massas      =massas;
gv.ro_shel     =ro_shel;
gv.Cf          =Cf;
gv.Cds         =Cds;
gv.Gammas      =Gammas;
gv.Gammas_ag2  =Gammas_ag2;
gv.dXs         =dXs;
gv.dphi        =dphi;
gv.nsol        =nsol;
gv.nsezs       =nsezs;
gv.nspi        =nspi;
gv.sezH        =sezH;

% global variable beam
gv.ro_water    =ro_water;
gv.massab      =massab;
gv.ro_arm      =ro_arm;
gv.gra         =gra;
gv.Lb          =Lb;
gv.Xb          =Xb;
gv.Rmax        =Rmax;
gv.Rmin        =Rmin;
gv.gbr         =gbr;
gv.xcib_star   =xcib_star;
gv.A           =A;
gv.Db          =Db;
gv.Eps         =Eps;
gv.Ipsi        =Ipsi;
gv.Eps_primo   =Eps_primo;
gv.Ipsi_primo  =Ipsi_primo;
gv.Gammab      =Gammab;
gv.Ag          =Ag;
gv.Gammab_ag   =Gammab_ag;
gv.nsezb       =nsezb;
gv.dXb         =dXb;
gv.tactb       =tactb;
gv.trelb       =trelb;
gv.deltam      =deltam;

% global variable rigid body
gv.F_bire2      =F_bire2;
gv.F_bire3      =F_bire3;
gv.F_bere2      =F_bere2;
gv.F_bere3      =F_bere3;
gv.F_biare2     =F_biare2;
gv.F_biare3     =F_biare3;
gv.pitch        =pitch;
gv.axes         =axes;
gv.magn         =magn;
gv.Volu         =Volu;
gv.nstep        =nstep;
gv.nstep       =nstep;
gv.Lr          =Lr;
gv.massar      =massar;
gv.ro_root     =ro_root;
gv.Crx         =Crx;
gv.Gammar      =Gammar;
gv.Gammar_ag   =Gammar_ag;

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% solution initialization 
% sol=[sigmas*(nsezs-1) xcis1*(nsezs-1) r*nsezs zeta*nsezs theta*nsezs
%      sigmab*(nsezb-1) xcib*(nsezb-1) gb*nsezb sigmar gr]

disp('Time-advancing')
myopt           =odeset('RelTol',1e-4);

%-------------------------------------------------------------------------
% condizioni iniziali temporali
sigmas_0        =zeros(6,nsezs-1);
xcis1_0         =[zeros(2,nsezs-1);mu0(1:nsezs-1);ones(1,nsezs-1);zeros(2,nsezs-1)];
sigmab_0        =zeros(6,nsezb-1);
xcib_0          =[zeros(3,nsezb);ones(1,nsezb);zeros(2,nsezb)];
gb_0            =zeros(4,nsezb*4);
for ii=1:nsezb
    gb_0(:,1+4*(ii-1):4+4*(ii-1))          =[0 -1 0 0; 0 0 1 0; -1 0 0 -Xb(ii); 0 0 0 1];
end
sigmar_0        =zeros(6,1);
gr_0            =[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
            
ini_cond        =[sigmas_0(1,:) sigmas_0(2,:) sigmas_0(3,:) sigmas_0(4,:) sigmas_0(5,:) sigmas_0(6,:)...    % sigmas
                  xcis1_0(1,:) xcis1_0(2,:) xcis1_0(3,:) xcis1_0(4,:) xcis1_0(5,:) xcis1_0(6,:)...                % xcis
                  erre0 zeta0 theta0...                                                                     % r z theta 
                  sigmab_0(1,:) sigmab_0(2,:) sigmab_0(3,:) sigmab_0(4,:) sigmab_0(5,:) sigmab_0(6,:)...    % sigmab1
                  xcib_0(1,:) xcib_0(2,:) xcib_0(3,:) xcib_0(4,:) xcib_0(5,:) xcib_0(6,:)...                % xcib1
                  gb_0(1,1:4:end) gb_0(2,1:4:end) gb_0(3,1:4:end)...                                        % tb1 
                  gb_0(1,2:4:end) gb_0(2,2:4:end) gb_0(3,2:4:end)...                                        % nb1
                  gb_0(1,3:4:end) gb_0(2,3:4:end) gb_0(3,3:4:end)...                                        % bb1
                  gb_0(1,4:4:end) gb_0(2,4:4:end) gb_0(3,4:4:end)...                                        % ub1
                  sigmar_0' gr_0(1:3,1)' gr_0(1:3,2)' gr_0(1:3,3)' gr_0(1:3,4)'];                           % sigmar tr nr br ur                                                   

% integrate
[t,z]           =ode23(@octopus_derivatives,tspan,ini_cond);

toc
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% postproc

disp('Post-processing')
nsol=size(z,1);

%-------------------------------------------------------------------------
% ricostruisco la BC

% condizioni su X=0
sigmas_o        =zeros(6,nsol);
sigmab_o        =zeros(6,nsol);

% get solution
r               =z(:,12*nsezs-11:13*nsezs-12)';
zeta            =z(:,13*nsezs-11:14*nsezs-12)';
theta           =z(:,14*nsezs-11:15*nsezs-12)';
senth           =sin(theta);
costh           =cos(theta);

% strain in direzione phi
tau             =zeros(nsezs,nsol);
varsigma        =zeros(nsezs,nsol);
for ii=1:nsol
    tau(:,ii)         =r(:,ii)./erre0';
    varsigma(:,ii)    =senth(:,ii)./erre0';
end

% condizioni su X=L ... xcib(L) zero è e zero rimarrà
eta_L           =zeros(1,nsol);
lambda_L        =zeros(1,nsol);
mu_L            =zeros(1,nsol);
for ii=1:nsol
    lamqu             =1-Pois*(tau(nsezs,ii)^2-1);
    lambda_L(ii)      =sqrt(lamqu);
    mu_L(ii)          =mu0(nsezs)/lambda_L(ii)+(Pois/lambda_L(ii))*(varsigma0(nsezs)-tau(nsezs,ii)*varsigma(nsezs,ii)); 
end
xcis_L           =[zeros(2,nsol);mu_L;lambda_L;eta_L;zeros(1,nsol)];
%-------------------------------------------------------------------------

Ns      =15*nsezs-12;
Z       =[sigmas_o(1,:)' z(:,1:nsezs-1) sigmas_o(2,:)' z(:,nsezs:2*nsezs-2) sigmas_o(3,:)' z(:,2*nsezs-1:3*nsezs-3)...
          sigmas_o(4,:)' z(:,3*nsezs-2:4*nsezs-4) sigmas_o(5,:)' z(:,4*nsezs-3:5*nsezs-5) sigmas_o(6,:)' z(:,5*nsezs-4:6*nsezs-6)...
          z(:,6*nsezs-5:7*nsezs-7) xcis_L(1,:)' z(:,7*nsezs-6:8*nsezs-8) xcis_L(2,:)' z(:,8*nsezs-7:9*nsezs-9) xcis_L(3,:)'...
          z(:,9*nsezs-8:10*nsezs-10) xcis_L(4,:)' z(:,10*nsezs-9:11*nsezs-11) xcis_L(5,:)' z(:,11*nsezs-10:12*nsezs-12) xcis_L(6,:)'...
          z(:,12*nsezs-11:15*nsezs-12)...                        % completo xcis con xcis_L e sigmas sigmas_o
          sigmab_o(1,:)' z(:,Ns+1:Ns+1*nsezb-1) sigmab_o(2,:)' z(:,Ns+1*(nsezb-1)+1:Ns+2*(nsezb-1)) sigmab_o(3,:)' z(:,Ns+2*(nsezb-1)+1:Ns+3*(nsezb-1))...
          sigmab_o(4,:)' z(:,Ns+3*(nsezb-1)+1:Ns+4*(nsezb-1)) sigmab_o(5,:)' z(:,Ns+4*(nsezb-1)+1:Ns+5*(nsezb-1)) sigmab_o(6,:)' z(:,Ns+5*(nsezb-1)+1:Ns+6*(nsezb-1))...
          z(:,Ns+6*(nsezb-1)+1:end)];                           % completo sigmab4 con sigmab_o

% octopus_postproc(t,Z)               % one plot
octopus_postproc2(t,Z)               % two subplots

toc

% end
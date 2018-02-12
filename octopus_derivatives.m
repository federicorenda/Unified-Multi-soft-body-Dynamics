function z_punto = octopus_derivatives(t,z)
global gv

t

% shell
gsr         =gv.gsr;
Area_out    =gv.Area_out;
Area_in     =gv.Area_in;
nus         =gv.nus;
Gs          =gv.Gs;
Pois        =gv.Pois;
tspan       =gv.tspan;
h           =gv.h;
per         =gv.per;
ten         =gv.ten;
erre0       =gv.erre0;
erre0_primo =gv.erre0_primo;
mu0         =gv.mu0;
varsigma0   =gv.varsigma0;
varsigma0_pr=gv.varsigma0_pr;
mu0_primo   =gv.mu0_primo;
Aream       =gv.Aream;
Js          =gv.Js;
Epsi        =gv.Epsi;
Epsi2       =gv.Epsi2;
ro_water    =gv.ro_water;
massas      =gv.massas;
ro_shel     =gv.ro_shel;
Cf          =gv.Cf;
Cds         =gv.Cds;
Gammas      =gv.Gammas;
Gammas_ag2  =gv.Gammas_ag2;
nsezs       =gv.nsezs;
nspi        =gv.nspi;
dXs         =gv.dXs;
dphi        =gv.dphi;
cosphi      =gv.cosphi;
senphi      =gv.senphi;

% beam
gbr         =gv.gbr;
massab      =gv.massab;
ro_arm      =gv.ro_arm;
gra         =gv.gra;
Xb          =gv.Xb;
A           =gv.A;
Db          =gv.Db;
Eps         =gv.Eps;
Ipsi        =gv.Ipsi;
Eps_primo   =gv.Eps_primo;
Ipsi_primo  =gv.Ipsi_primo;
Gammab      =gv.Gammab;
Ag          =gv.Ag;
Gammab_ag   =gv.Gammab_ag;
nsezb       =gv.nsezb;
dXb         =gv.dXb;
xcib_star   =gv.xcib_star;
tactb       =gv.tactb;
trelb       =gv.trelb;
deltam      =gv.deltam;

% rigid body
massar      =gv.massar;
ro_root     =gv.ro_root;
Crx         =gv.Crx;
Gammar      =gv.Gammar;
Gammar_ag   =gv.Gammar_ag;

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% actual solution sigma xci g

%-------------------------------------------------------------------------
% Shell
% condizioni su X=0
sigmas_o        =zeros(6,1);

% cinematica differenziale mantello
sigmas          =[sigmas_o';...
                 z(1:nsezs-1) z(nsezs:2*nsezs-2) z(2*nsezs-1:3*nsezs-3) z(3*nsezs-2:4*nsezs-4) z(4*nsezs-3:5*nsezs-5) z(5*nsezs-4:6*nsezs-6)]';
we              =sigmas(3,:);
va              =sigmas(4,:);
vb              =sigmas(5,:);

% cinematica mantello
erre            =z(12*nsezs-11:13*nsezs-12)';
zeta            =z(13*nsezs-11:14*nsezs-12)';
theta           =z(14*nsezs-11:15*nsezs-12)';
costh           =cos(theta);
senth           =sin(theta);

% strain in direzione phi
tau             =erre./erre0;
varsigma        =senth./erre0;

% condizioni su X=L
eta_L           =0;
lamqu           =1-Pois*(tau(nsezs)^2-1);
lambda_L        =sqrt(lamqu);
mu_L            =mu0(nsezs)/lambda_L+(Pois/lambda_L)*(varsigma0(nsezs)-tau(nsezs)*varsigma(nsezs));        
xcis1_L          =[0;0;mu_L;lambda_L;eta_L;0];

% deformazioni mantello
xcis1            =[z(6*nsezs-5:7*nsezs-7) z(7*nsezs-6:8*nsezs-8) z(8*nsezs-7:9*nsezs-9) z(9*nsezs-8:10*nsezs-10) z(10*nsezs-9:11*nsezs-11) z(11*nsezs-10:12*nsezs-12);...
                  xcis1_L']';
ks1              =xcis1(1:3,:);
q1               =xcis1(4:6,:);
mu               =ks1(3,:);
lambda           =q1(1,:);
eta              =q1(2,:);

%-------------------------------------------------------------------------
% Beam
% screw
Ns              =15*nsezs-12;               % grandezza vettore di stato shell
sigmab_o        =zeros(6,1);
sigmab          =[sigmab_o';...
                 [z(Ns+1:Ns+1*nsezb-1) z(Ns+1*(nsezb-1)+1:Ns+2*(nsezb-1)) z(Ns+2*(nsezb-1)+1:Ns+3*(nsezb-1)) z(Ns+3*(nsezb-1)+1:Ns+4*(nsezb-1)) z(Ns+4*(nsezb-1)+1:Ns+5*(nsezb-1)) z(Ns+5*(nsezb-1)+1:Ns+6*(nsezb-1))]]';
xcib            =[z(Ns+6*(nsezb-1)+1:Ns+7*nsezb-6) z(Ns+7*nsezb-5:Ns+8*nsezb-6) z(Ns+8*nsezb-5:Ns+9*nsezb-6) z(Ns+9*nsezb-5:Ns+10*nsezb-6) z(Ns+10*nsezb-5:Ns+11*nsezb-6) z(Ns+11*nsezb-5:Ns+12*nsezb-6)]';
gb              =zeros(nsezb*4,4);
gb(1:4:end,1:3) =[z(Ns+12*nsezb-5:Ns+13*nsezb-6) z(Ns+13*nsezb-5:Ns+14*nsezb-6) z(Ns+14*nsezb-5:Ns+15*nsezb-6)];   % tb
gb(2:4:end,1:3) =[z(Ns+15*nsezb-5:Ns+16*nsezb-6) z(Ns+16*nsezb-5:Ns+17*nsezb-6) z(Ns+17*nsezb-5:Ns+18*nsezb-6)];   % nb
gb(3:4:end,1:3) =[z(Ns+18*nsezb-5:Ns+19*nsezb-6) z(Ns+19*nsezb-5:Ns+20*nsezb-6) z(Ns+20*nsezb-5:Ns+21*nsezb-6)];   % bb
gb(4:4:end,1:3) =[z(Ns+21*nsezb-5:Ns+22*nsezb-6) z(Ns+22*nsezb-5:Ns+23*nsezb-6) z(Ns+23*nsezb-5:Ns+24*nsezb-6)];   % ub
gb(4:4:end,4)   =ones(nsezb,1);
gb              =gb';              

% cinematica e cinematica differenziale corpo rigido
sigmar          =z(Ns+24*nsezb-5:Ns+24*nsezb);
gr              =[z(Ns+24*nsezb+1:Ns+24*nsezb+3) z(Ns+24*nsezb+4:Ns+24*nsezb+6)...                         % tr nr
                  z(Ns+24*nsezb+7:Ns+24*nsezb+9) z(Ns+24*nsezb+10:Ns+24*nsezb+12)];                        % br ur
gr              =[gr; [0 0 0 1]];

%-------------------------------------------------------------------------
% velocity screws assolute per shell e beam1
% Shell
Sigmas          =zeros(6,nsezs);
for ii=1:nsezs
    gs                =[cosphi(1)*costh(ii) -cosphi(1)*senth(ii) senphi(1) cosphi(1)*erre(ii);...
                        senphi(1)*costh(ii) -senphi(1)*senth(ii) -cosphi(1) senphi(1)*erre(ii);...
                        senth(ii) costh(ii) 0 zeta(ii); 0 0 0 1];
    Adjoint_invgs     =dinamico_Adjoint((gsr*gs)^-1);
    Sigmas(:,ii)      =sigmas(:,ii)+Adjoint_invgs*[0;0;0;sigmar(4);0;0];        
end

% Beam
Sigmab          =zeros(6,nsezb);
for ii=1:nsezb
    Adjoint_invgb     =dinamico_Adjoint((gbr*gb(:,4*(ii-1)+1:4*(ii-1)+4))^-1);
    Sigmab(:,ii)      =sigmab(:,ii)+Adjoint_invgb*sigmar;
end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% derivate rispetto a X

% Shell
% derivata spaziale del 2° ordine centrata dentro e upwind ai bordi 1°
sigmas_primo                    =zeros(6,nsezs);
xcis1_primo                     =zeros(6,nsezs);
sigmas_primo(:,1)               =diff([sigmas(:,1) sigmas(:,2)],1,2)./dXs;
sigmas_primo(:,2:2:nsezs-1)     =diff(sigmas(:,1:2:nsezs),1,2)./(2*dXs);
sigmas_primo(:,3:2:nsezs-1)     =diff(sigmas(:,2:2:nsezs),1,2)./(2*dXs);
sigmas_primo(:,nsezs)           =diff([sigmas(:,nsezs-1) sigmas(:,nsezs)],1,2)./dXs;
xcis1_primo(:,1)                =diff([xcis1(:,1) xcis1(:,2)],1,2)./dXs;
xcis1_primo(:,2:2:nsezs-1)      =diff(xcis1(:,1:2:nsezs),1,2)./(2*dXs);
xcis1_primo(:,3:2:nsezs-1)      =diff(xcis1(:,2:2:nsezs),1,2)./(2*dXs);
xcis1_primo(:,nsezs)            =diff([xcis1(:,nsezs-1) xcis1(:,nsezs)],1,2)./dXs;

% Beam
% derivate rispetto a X upwinding
sigmab_primo                    =diff([sigmab sigmab(:,nsezb)],1,2)./dXb;
xcib_primo                      =diff([xcib(:,1) xcib],1,2)./dXb;

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% calcolo le derivate

%-------------------------------------------------------------------------
% kinematics

% Rigid Body
gr_punto           =gr*dinamico_hat(sigmar);

% Shell
r_punto            =costh.*va-senth.*vb;
zeta_punto         =costh.*vb+senth.*va;
theta_punto        =we;

% Beam
gb_punto           =zeros(4,4*nsezb);
for ii=1:nsezb
    gb_punto(:,4*(ii-1)+1:4*(ii-1)+4)    =gb(:,4*(ii-1)+1:4*(ii-1)+4)*dinamico_hat(sigmab(:,ii));
end

%-------------------------------------------------------------------------
% compatibility
% Shell
xcis1_punto                   =zeros(6,nsezs);
for ii=1:nsezs-1
    xcis1_punto(:,ii)         =sigmas_primo(:,ii)+dinamico_adj(xcis1(:,ii))*sigmas(:,ii);
end

tau_punto           =va.*costh./erre0-vb.*varsigma;
varsigma_punto      =we.*costh./erre0;

% condizioni su X=L
etaL_punto          =0;
lambdaL_punto       =-Pois*tau(nsezs)*tau_punto(nsezs)/lambda_L;
muL_punto           =-(lambdaL_punto/lambda_L^2)*(mu0(nsezs)+Pois*(varsigma0(nsezs)-tau(nsezs)*varsigma(nsezs)))...
                     -(Pois/lambda_L)*(tau(nsezs)*varsigma_punto(nsezs)+tau_punto(nsezs)*varsigma(nsezs));
xcis1L_punto        =[0;0;muL_punto;lambdaL_punto;etaL_punto;0];
xcis1_punto(:,nsezs)=xcis1L_punto;

% Beam
xcib_punto          =zeros(6,nsezb);
for ii=1:nsezb-1
    xcib_punto(:,ii)          =sigmab_primo(:,ii)+dinamico_adj(xcib(:,ii))*sigmab(:,ii);
end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% derivate rispetto a X

% Shell
% derivata spaziale del 2° ordine centrata dentro e upwind ai bordi 1°
xcis1_punto_primo                =zeros(6,nsezs);
xcis1_punto_primo(:,1)           =diff([xcis1_punto(:,1) xcis1_punto(:,2)],1,2)./dXs;
xcis1_punto_primo(:,2:2:nsezs-1) =diff(xcis1_punto(:,1:2:nsezs),1,2)./(2*dXs);
xcis1_punto_primo(:,3:2:nsezs-1) =diff(xcis1_punto(:,2:2:nsezs),1,2)./(2*dXs);
xcis1_punto_primo(:,nsezs)       =diff([xcis1_punto(:,nsezs-1) xcis1_punto(:,nsezs)],1,2)./dXs;

% Beam
xcib_punto_primo                 =diff([xcib_punto(:,1) xcib_punto],1,2)./dXb;

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% dinamica

%-------------------------------------------------------------------------
% Internal Wrenches shell
mu_punto          =xcis1_punto(3,:);
lambda_punto      =xcis1_punto(4,:);
eta_punto         =xcis1_punto(5,:);

e_11              =0.5*(lambda.^2+eta.^2-ones(1,nsezs));
e_22              =0.5*(tau.^2-ones(1,nsezs));
d_11              =mu0-mu.*lambda;
d_22              =varsigma0-tau.*varsigma;

ep_11             =lambda.*lambda_punto+eta.*eta_punto;
ep_22             =tau.*tau_punto;
dp_11             =-mu_punto.*lambda-mu.*lambda_punto;
dp_22             =-tau_punto.*varsigma-tau.*varsigma_punto;

E_11              =e_11+Pois.*e_22;
E_22              =e_22+Pois.*e_11;
D_11              =d_11+Pois.*d_22;
D_22              =d_22+Pois.*d_11;

Ep_11             =ep_11+Pois.*ep_22;
Ep_22             =ep_22+Pois.*ep_11;
Dp_11             =dp_11+Pois.*dp_22;
Dp_22             =dp_22+Pois.*dp_11;

%-------------------------------------------------------------------------
w_primo           =sigmas_primo(3,:);
va_primo          =sigmas_primo(4,:);
vb_primo          =sigmas_primo(5,:);

mu_primo          =xcis1_primo(3,:);
lambda_primo      =xcis1_primo(4,:);
eta_primo         =xcis1_primo(5,:);

tau_primo         =lambda.*costh./erre0-eta.*varsigma-tau.*(erre0_primo./erre0);
varsigma_pr       =mu.*costh./erre0-varsigma.*(erre0_primo./erre0);

mu_punto_primo    =xcis1_punto_primo(3,:);
lambda_punto_primo=xcis1_punto_primo(4,:);
eta_punto_primo   =xcis1_punto_primo(5,:);

tau_punto_primo   =va_primo.*costh./erre0-vb_primo.*varsigma...
                   -va.*mu.*varsigma-vb.*mu.*costh./erre0...
                   -tau_punto.*erre0_primo./erre0;
varsigma_punto_pr =w_primo.*costh./erre0-we.*mu.*varsigma...
                   -varsigma_punto.*erre0_primo./erre0;

e_11_primo        =lambda.*lambda_primo+eta.*eta_primo;
e_22_primo        =tau.*tau_primo;
d_11_primo        =mu0_primo-mu.*lambda_primo-mu_primo.*lambda;
d_22_primo        =varsigma0_pr-tau_primo.*varsigma-tau.*varsigma_pr;

ep_11_primo       =lambda_primo.*lambda_punto+lambda.*lambda_punto_primo...
                   +eta_primo.*eta_punto+eta.*eta_punto_primo;
ep_22_primo       =tau_primo.*tau_punto+tau.*tau_punto_primo;
dp_11_primo       =-mu_punto_primo.*lambda-mu_punto.*lambda_primo...
                   -mu_primo.*lambda_punto-mu.*lambda_punto_primo;
dp_22_primo       =-tau_punto_primo.*varsigma-tau_punto.*varsigma_pr...
                   -tau_primo.*varsigma_punto-tau.*varsigma_punto_pr;
          
E_11_primo        =e_11_primo+Pois.*e_22_primo;
D_11_primo        =d_11_primo+Pois.*d_22_primo;

Ep_11_primo       =ep_11_primo+Pois.*ep_22_primo;
Dp_11_primo       =dp_11_primo+Pois.*dp_22_primo;

%-------------------------------------------------------------------------
ntil11            =Epsi*E_11+Epsi2*Ep_11;
ntil22            =Epsi*E_22+Epsi2*Ep_22;
mtil11            =Epsi*Js*D_11+Epsi2*Js*Dp_11;
mtil22            =Epsi*Js*D_22+Epsi2*Js*Dp_22;
q1                =2*h*Gs*eta+2*h*nus*eta_punto;

ntil11_primo      =Epsi*E_11_primo+Epsi2*Ep_11_primo;
mtil11_primo      =Epsi*Js*D_11_primo+Epsi2*Js*Dp_11_primo;
q1_primo          =2*h*Gs*eta_primo+2*h*nus*eta_punto_primo;

%-------------------------------------------------------------------------
N_X               =lambda.*ntil11 - mu.*mtil11;
N_phi             =tau.*ntil22-varsigma.*mtil22;
H                 =q1+eta.*ntil11;
M_X               =-lambda.*mtil11;
M_phi             =-tau.*mtil22;
N_X_primo         =lambda_primo.*ntil11+lambda.*ntil11_primo...
                   -mu_primo.*mtil11-mu.*mtil11_primo;
H_primo           =q1_primo+eta_primo.*ntil11+eta.*ntil11_primo;
M_X_primo         =-lambda_primo.*mtil11-lambda.*mtil11_primo;
zetas1_i          =[zeros(2,nsezs);M_X;N_X;H;zeros(1,nsezs)];
zetas2_i          =[M_phi;zeros(4,nsezs);-N_phi];
zetas1_i_primo    =[zeros(2,nsezs);M_X_primo;N_X_primo;H_primo;zeros(1,nsezs)];

%-------------------------------------------------------------------------
% external and actuation wrenches

% shell
xcis2             =[varsigma;costh./erre0;zeros(3,nsezs);-tau];
r0pr_r0           =erre0_primo./erre0;

% calcolo la velocità di nuotata
sigmar_bar    =dinamico_Adjoint(gsr^-1)*sigmar;
vel           =sigmar_bar(6);

% calcolo geometrie
Aref          =pi*max(erre).^2;
VolX_punto    =pi*dXs*cumtrapz((2.*erre.*r_punto+erre.^2.*2.*ep_11./(2.*e_11+ones(1,nsezs))).*sqrt(2.*e_11+ones(1,nsezs)));
Vol_punto     =VolX_punto(nsezs)-VolX_punto(1);
Ten           =interp1([0 1/2 1],ten,mod(t,per)/per);

% calcolo moduli delle forze fluidiche
mod_thrust    =-sign(Vol_punto)*(ro_water*Vol_punto^2)*...
               ((Vol_punto<0)*Cf/(Area_out*Aream)+(Vol_punto>0)*Cf/(Area_in*Aream));
mod_drag      =-(1/2)*Cds*Aref*ro_water*abs(vel)*vel/Aream;
mod_fluid     =mod_thrust+mod_drag;

% beam
act_beam      =zeros(1,nsezb);
act_beam_primo=zeros(1,nsezb);
if t<=tactb                                         % prima virata
    force_beam    =linspace(deltam,0,floor(nsezb*t/tactb));
    act_beam      =[-force_beam zeros(1,nsezb-floor(nsezb*t/tactb))];
    if floor(nsezb*t/tactb)>=2
        act_beam_primo=[deltam/Xb(floor(nsezb*t/tactb))*ones(1,floor(nsezb*t/tactb)) zeros(1,nsezb-floor(nsezb*t/tactb))];
    end
end
if (t>tactb && t<tactb+trelb)                      % primo rilascio
    act_beam      =zeros(1,nsezb);
    act_beam_primo=zeros(1,nsezb);
end
if (t>=tactb+trelb && t<=2*tactb+trelb)              % seconda virata
    force_beam    =linspace(deltam,0,floor(nsezb*(t-(tactb+trelb))/tactb));
    act_beam      =[-force_beam zeros(1,nsezb-floor(nsezb*(t-(tactb+trelb))/tactb))];
    if floor(nsezb*(t-(tactb+trelb))/tactb)>=2
        act_beam_primo=[deltam/Xb(floor(nsezb*(t-(tactb+trelb))/tactb))*ones(1,floor(nsezb*(t-(tactb+trelb))/tactb))...
                        zeros(1,nsezb-floor(nsezb*(t-(tactb+trelb))/tactb))];
    end
end
if t>2*tactb+trelb                                  % secondo rilascio
    act_beam      =zeros(1,nsezb);
    act_beam_primo=zeros(1,nsezb);
end

%----------------------------------------------------------------------
% equilibrium

% Rigid Body
% intrisic loads
F_gal           =-(massab*(1-ro_water/ro_arm)+massas*(1-ro_water/ro_shel)...
                   +massar*(1-ro_water/ro_root))*gr^-1*[gra;0];
F_peso          =(1-ro_water/ro_root)*Gammar*dinamico_Adjoint(gr^-1)*[0;0;0;gra];
F_drag          =-((1/2)*Crx*ro_water*sigmar(4:6).^2).*sign(sigmar(4:6));
F_ir            =[0;0;0;F_gal(1:3)]+[0;0;0;F_drag]+F_peso;
% beam boundary load
F_ib            =Eps(:,1:6)*(xcib(:,1)-xcib_star(:,1)) + Ipsi(:,1:6)*(xcib_punto(:,1));
F_ab            =[0;0;act_beam(1);0;0;0]*(t<tactb+trelb)+...
                 [0;act_beam(1);0;0;0;0]*(t>=tactb+trelb);
F_ibr           =dinamico_coAdjoint(gbr*gb(:,1:4))*F_ib;
F_abr           =dinamico_coAdjoint(gbr*gb(:,1:4))*F_ab;
% shell boundary load
F_is1r          =zeros(6,1);
F_is2r          =zeros(6,1);
for jj=1:nspi
    gs          =[cosphi(jj)*costh(1) -cosphi(jj)*senth(1) senphi(jj) cosphi(jj)*erre(1);...
                  senphi(jj)*costh(1) -senphi(jj)*senth(1) -cosphi(jj) senphi(jj)*erre(1);...
                  senth(1) costh(1) 0 zeta(1); 0 0 0 1];
    F_is1r       =F_is1r+erre0(1)*dinamico_coAdjoint(gsr*gs)*zetas1_i(:,1);
end
F_is1r          =dphi*F_is1r;
for ii=1:nsezs
    for jj=1:nspi
        gs      =[cosphi(jj)*costh(ii) -cosphi(jj)*senth(ii) senphi(jj) cosphi(jj)*erre(ii);...
                  senphi(jj)*costh(ii) -senphi(jj)*senth(ii) -cosphi(jj) senphi(jj)*erre(ii);...
                  senth(ii) costh(ii) 0 zeta(ii); 0 0 0 1];
        F_is2r  =F_is2r+erre0(ii)*(1-ro_water/ro_shel)*dinamico_coAdjoint(gsr*gs)*Gammas*dinamico_Adjoint((gr*gsr*gs)^-1)*[0;0;0;gra];
    end
end
F_is2r          =dXs*dphi*F_is2r;

F_r             =F_ibr-F_abr+F_is1r+F_is2r+F_ir;
sigmar_punto    =Gammar_ag^-1*(F_r+dinamico_coadj(sigmar)*(Gammar_ag*sigmar));

% Shell
% initialization
sigmas_punto        =zeros(6,nsezs);

for ii=2:nsezs
    
    % external and actuation wrench
    fluid                 =mod_fluid*[senth(ii);costh(ii);0];
    zetas2_a              =[0;0;0;0;0;-Ten];
    zetas_e               =[0;0;0;fluid];
    
    % accelerazione di trascinamento
    gs                    =[cosphi(1)*costh(ii) -cosphi(1)*senth(ii) senphi(1) cosphi(1)*erre(ii);...
                            senphi(1)*costh(ii) -senphi(1)*senth(ii) -cosphi(1) senphi(1)*erre(ii);...
                            senth(ii) costh(ii) 0 zeta(ii); 0 0 0 1];
    Adjoint_invgs         =dinamico_Adjoint((gsr*gs)^-1);
    sigmar_bar_punto      =Adjoint_invgs*[0;0;0;sigmar_punto(4);0;0];
    sigmar_bar            =Adjoint_invgs*[0;0;0;sigmar(4);0;0];
    
    % internal dynamic
    sigmas_punto(:,ii)    =Gammas_ag2^-1*(r0pr_r0(ii)*zetas1_i(:,ii) + zetas1_i_primo(:,ii)...
                           - dinamico_coadj(xcis1(:,ii))*zetas1_i(:,ii) - dinamico_coadj(xcis2(:,ii))*(zetas2_i(:,ii)-zetas2_a)...
                           + zetas_e + dinamico_coadj(Sigmas(:,ii))*(Gammas*Sigmas(:,ii))) - sigmar_bar_punto + dinamico_adj(sigmas(:,ii))*sigmar_bar;
end

% Beam
% initialization
sigmab_punto    =zeros(6,nsezb);

for ii=2:nsezb
    
    % internal actuation
    act                  =[0;0;act_beam(ii);0;0;0]*(t<=tactb+trelb)+...
                          [0;act_beam(ii);0;0;0;0]*(t>tactb+trelb);
    act_primo            =[0;0;act_beam_primo(ii);0;0;0]*(t<=tactb+trelb)+...
                          [0;act_beam_primo(ii);0;0;0;0]*(t>tactb+trelb);
    
    % internal wrench
    intw                 =Eps(:,6*(ii-1)+1:6*(ii-1)+6)*(xcib(:,ii)-xcib_star(:,ii)) + Ipsi(:,6*(ii-1)+1:6*(ii-1)+6)*(xcib_punto(:,ii)) - act;
    intw_primo           =Eps(:,6*(ii-1)+1:6*(ii-1)+6)*(xcib_primo(:,ii)) + Ipsi(:,6*(ii-1)+1:6*(ii-1)+6)*(xcib_punto_primo(:,ii))...
                          +Eps_primo(:,6*(ii-1)+1:6*(ii-1)+6)*(xcib(:,ii)-xcib_star(:,ii)) + Ipsi_primo(:,6*(ii-1)+1:6*(ii-1)+6)*(xcib_punto(:,ii))...
                          -act_primo;
    
    % gravity & buoyancy
    g_b                  =gr*gbr*gb(:,4*(ii-1)+1:4*(ii-1)+4);
    gra_local            =g_b^-1*[gra;0];
    grbu                 =(ro_arm-ro_water)*A(ii)*gra_local(1:3);
    
    % drag load
    drag                 =-(ro_water*Db(:,3*(ii-1)+1:3*(ii-1)+3)*((Sigmab(4:6,ii)).^2)).*sign(Sigmab(4:6,ii));
    
    % external wrench
    extw_primo           =[0;0;0;grbu]+[0;0;0;drag];                        
    
    % accelerazione di trascinamento
    Adjoint_invgb        =dinamico_Adjoint((gbr*gb(:,4*(ii-1)+1:4*(ii-1)+4))^-1);
    sigmar_bar_punto     =Adjoint_invgb*sigmar_punto;
    sigmar_bar           =Adjoint_invgb*sigmar;
    
    % beam equilibrium
    sigmab_punto(:,ii)  =Gammab_ag(:,6*(ii-1)+1:6*(ii-1)+6)^-1*(intw_primo + extw_primo - dinamico_coadj(xcib(:,ii))*intw...
                          + dinamico_coadj(Sigmab(:,ii))*(Gammab_ag(:,6*(ii-1)+1:6*(ii-1)+6)*Sigmab(:,ii)))...
                          - sigmar_bar_punto + dinamico_adj(sigmab(:,ii))*sigmar_bar;
end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% osservabili, SRA reaction loads, rigid motion screw, volume interno

nstep           =gv.nstep;
if t>tspan(nstep)
    F_bire2         =gv.F_bire2;
    F_bire3         =gv.F_bire3;
    F_bere2         =gv.F_bere2;
    F_bere3         =gv.F_bere3;
    F_biare2        =gv.F_biare2;
    F_biare3        =gv.F_biare3;
    pitch           =gv.pitch;
    axes            =gv.axes;
    magn            =gv.magn;
    Volu            =gv.Volu;
    
    % shell
    VolX            =pi*dXs*cumtrapz(erre.^2.*sqrt(2.*e_11+ones(1,nsezs)));
    Volu(nstep)     =VolX(nsezs)-VolX(1);
    
    % beam
    F_bir_dis       =zeros(6,nsezb);
    F_ber_dis       =zeros(6,nsezb);
    for ii=1:nsezb
        Adjoint_invgb        =dinamico_Adjoint((gbr*gb(:,4*(ii-1)+1:4*(ii-1)+4))^-1);
        sigmar_bar_punto     =Adjoint_invgb*sigmar_punto;
        sigmar_bar           =Adjoint_invgb*sigmar;
        F_bir_dis(:,ii)      =dinamico_coAdjoint(gbr*gb(:,4*(ii-1)+1:4*(ii-1)+4))*...
                              (-Gammab(:,6*(ii-1)+1:6*(ii-1)+6)*(sigmab_punto(:,ii)+sigmar_bar_punto - dinamico_adj(sigmab(:,ii))*sigmar_bar)+...
                              dinamico_coadj(Sigmab(:,ii))*(Gammab(:,6*(ii-1)+1:6*(ii-1)+6)*Sigmab(:,ii)));
        g_b                  =gr*gbr*gb(:,4*(ii-1)+1:4*(ii-1)+4);
        gra_local            =g_b^-1*[gra;0];
        grbu                 =(ro_arm-ro_water)*A(ii)*gra_local(1:3);
        drag                 =-(ro_water*Db(:,3*(ii-1)+1:3*(ii-1)+3)*((Sigmab(4:6,ii)).^2)).*sign(Sigmab(4:6,ii));
        addmas               =-Ag(:,6*(ii-1)+1:6*(ii-1)+6)*(sigmab_punto(:,ii)+sigmar_bar_punto - dinamico_adj(sigmab(:,ii))*sigmar_bar)+...
                              dinamico_coadj(Sigmab(:,ii))*(Ag(:,6*(ii-1)+1:6*(ii-1)+6)*Sigmab(:,ii));
        F_ber_dis(:,ii)      =dinamico_coAdjoint(gbr*gb(:,4*(ii-1)+1:4*(ii-1)+4))*([0;0;0;grbu]+[0;0;0;drag]+addmas);
    end
    F_birX          =dXb*cumtrapz(F_bir_dis,2);
    F_berX          =dXb*cumtrapz(F_ber_dis,2);
    F_bir           =F_birX(:,nsezb)-F_birX(:,1);
    F_ber           =F_berX(:,nsezb)-F_berX(:,1);
    F_biar          =F_ibr-F_abr;
    F_bire2(nstep)  =F_bir(2);
    F_bire3(nstep)  =F_bir(3);
    F_bere2(nstep)  =F_ber(2);
    F_bere3(nstep)  =F_ber(3);
    F_biare2(nstep) =F_biar(2);
    F_biare3(nstep) =F_biar(3);

    % rigid body
    pitch(nstep)    =sigmar(1:3)'*sigmar(4:6)/(sigmar(1:3)'*sigmar(1:3));
    axes(:,nstep)   =cross(sigmar(1:3),sigmar(4:6))/(sigmar(1:3)'*sigmar(1:3));
    magn(nstep)     =sqrt(sigmar(1:3)'*sigmar(1:3));
    
    % somma
    nstep           =nstep+1;
    
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
end

%----------------------------------------------------------------------
z_punto         =[sigmas_punto(1,2:nsezs)';sigmas_punto(2,2:nsezs)';sigmas_punto(3,2:nsezs)';sigmas_punto(4,2:nsezs)';sigmas_punto(5,2:nsezs)';sigmas_punto(6,2:nsezs)';...
                  xcis1_punto(1,1:nsezs-1)';xcis1_punto(2,1:nsezs-1)';xcis1_punto(3,1:nsezs-1)';xcis1_punto(4,1:nsezs-1)';xcis1_punto(5,1:nsezs-1)';xcis1_punto(6,1:nsezs-1)';...
                  r_punto'; zeta_punto'; theta_punto';...
                  sigmab_punto(1,2:nsezb)';sigmab_punto(2,2:nsezb)';sigmab_punto(3,2:nsezb)';sigmab_punto(4,2:nsezb)';sigmab_punto(5,2:nsezb)';sigmab_punto(6,2:nsezb)';...
                  xcib_punto(1,:)';xcib_punto(2,:)';xcib_punto(3,:)';xcib_punto(4,:)';xcib_punto(5,:)';xcib_punto(6,:)';...
                  gb_punto(1,1:4:end)';gb_punto(2,1:4:end)';gb_punto(3,1:4:end)';...            % tb
                  gb_punto(1,2:4:end)';gb_punto(2,2:4:end)';gb_punto(3,2:4:end)';...            % nb
                  gb_punto(1,3:4:end)';gb_punto(2,3:4:end)';gb_punto(3,3:4:end)';...            % bb
                  gb_punto(1,4:4:end)';gb_punto(2,4:4:end)';gb_punto(3,4:4:end)';...            % ub
                  sigmar_punto;gr_punto(1:3,1);gr_punto(1:3,2);gr_punto(1:3,3);gr_punto(1:3,4)];   % nur tr nr br ur
              
% eof
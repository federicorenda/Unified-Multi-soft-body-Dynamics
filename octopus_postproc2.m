function octopus_postproc2(t,z)
global gv

% shell
gsr     =gv.gsr;
Xs      =gv.Xs;
nsol    =gv.nsol;
nsezs   =gv.nsezs;
F_bire2 =gv.F_bire2;
F_bire3 =gv.F_bire3;
F_bere2 =gv.F_bere2;
F_bere3 =gv.F_bere3;
F_biare2=gv.F_biare2;
F_biare3=gv.F_biare3;
pitch   =gv.pitch;
axes    =gv.axes;
magn    =gv.magn;
Volu    =gv.Volu;

% beam
gbr     =gv.gbr;
Lb      =gv.Lb;
Rmax    =gv.Rmax;
Rmin    =gv.Rmin;
nsezb   =gv.nsezb;

% rigid body
Lr      =gv.Lr;

%-------------------------------------------------------------------------
% pre-processing

Xb       =linspace(0,Lb,nsezb);
[Ts,Ss]  =meshgrid(t,Xs);
[Tb,Sb]  =meshgrid(t,Xb);
mkdir('.\LAST RUN\');
sradius  =zeros(1,nsol);
for ii=1:nsol
    sradius(ii) =norm(axes(:,ii));
end

% get the solution
% Shell
we       =z(:,2*nsezs+1:3*nsezs);
va       =z(:,3*nsezs+1:4*nsezs);
vb       =z(:,4*nsezs+1:5*nsezs);
mu       =z(:,8*nsezs+1:9*nsezs);
lambda   =z(:,9*nsezs+1:10*nsezs);
eta      =z(:,10*nsezs+1:11*nsezs);
erre     =z(:,12*nsezs+1:13*nsezs);
zeta     =z(:,13*nsezs+1:14*nsezs);
theta    =z(:,14*nsezs+1:15*nsezs);
senth    =sin(theta);
costh    =cos(theta);

% Beam1
Ns       =15*nsezs;
wbt      =z(:,Ns+1:Ns+1*nsezb);
wbn      =z(:,Ns+1*nsezb+1:Ns+2*nsezb);
wbb      =z(:,Ns+2*nsezb+1:Ns+3*nsezb);
vbt      =z(:,Ns+3*nsezb+1:Ns+4*nsezb);
vbn      =z(:,Ns+4*nsezb+1:Ns+5*nsezb);
vbb      =z(:,Ns+5*nsezb+1:Ns+6*nsezb);
tau      =z(:,Ns+6*nsezb+1:Ns+7*nsezb);
xci      =z(:,Ns+7*nsezb+1:Ns+8*nsezb);
k        =z(:,Ns+8*nsezb+1:Ns+9*nsezb);
q        =z(:,Ns+9*nsezb+1:Ns+10*nsezb);
p        =z(:,Ns+10*nsezb+1:Ns+11*nsezb);
r        =z(:,Ns+11*nsezb+1:Ns+12*nsezb);
tbx      =z(:,Ns+12*nsezb+1:Ns+13*nsezb);
tby      =z(:,Ns+13*nsezb+1:Ns+14*nsezb);
tbz      =z(:,Ns+14*nsezb+1:Ns+15*nsezb);
nbx      =z(:,Ns+15*nsezb+1:Ns+16*nsezb);
nby      =z(:,Ns+16*nsezb+1:Ns+17*nsezb);
nbz      =z(:,Ns+17*nsezb+1:Ns+18*nsezb);
bbx      =z(:,Ns+18*nsezb+1:Ns+19*nsezb);
bby      =z(:,Ns+19*nsezb+1:Ns+20*nsezb);
bbz      =z(:,Ns+20*nsezb+1:Ns+21*nsezb);
ubx      =z(:,Ns+21*nsezb+1:Ns+22*nsezb);
uby      =z(:,Ns+22*nsezb+1:Ns+23*nsezb);
ubz      =z(:,Ns+23*nsezb+1:Ns+24*nsezb);

% rigid body
wr_x     =z(:,Ns+24*nsezb+1);
wr_y     =z(:,Ns+24*nsezb+2);
wr_z     =z(:,Ns+24*nsezb+3);
vr_x     =z(:,Ns+24*nsezb+4);
vr_y     =z(:,Ns+24*nsezb+5);
vr_z     =z(:,Ns+24*nsezb+6);
tr_x     =z(:,Ns+24*nsezb+7);
tr_y     =z(:,Ns+24*nsezb+8);
tr_z     =z(:,Ns+24*nsezb+9);
nr_x     =z(:,Ns+24*nsezb+10);
nr_y     =z(:,Ns+24*nsezb+11);
nr_z     =z(:,Ns+24*nsezb+12);
br_x     =z(:,Ns+24*nsezb+13);
br_y     =z(:,Ns+24*nsezb+14);
br_z     =z(:,Ns+24*nsezb+15);
ur_x     =z(:,Ns+24*nsezb+16);
ur_y     =z(:,Ns+24*nsezb+17);
ur_z     =z(:,Ns+24*nsezb+18);

% Screw
w        =[wr_x';wr_y';wr_z'];

%-------------------------------------------------------------------------
% save risults
% Shell
save('.\LAST RUN\postproc','t','z')
save('.\LAST RUN\long vel','t','va');
save('.\LAST RUN\tras b vel','t','vb');
save('.\LAST RUN\radius','t','erre');
save('.\LAST RUN\altitude','t','zeta');
save('.\LAST RUN\angle','t','theta');
save('.\LAST RUN\vel rotazione su -e_phi','t','we');
save('.\LAST RUN\long strain','t','lambda');
save('.\LAST RUN\tras b strain','t','eta');
save('.\LAST RUN\curvatura su -e_phi','t','mu');

% Beam
save('.\LAST RUN\long vel','t','vbt');
save('.\LAST RUN\tras n vel','t','vbn');
save('.\LAST RUN\tras b vel','t','vbb');
save('.\LAST RUN\vel rotazione su t','t','wbt');
save('.\LAST RUN\vel rotazione su n','t','wbn');
save('.\LAST RUN\vel rotazione su b','t','wbb');
save('.\LAST RUN\long strain','t','q');
save('.\LAST RUN\tras n strain','t','p');
save('.\LAST RUN\tras b strain','t','r');
save('.\LAST RUN\torsione','t','tau');
save('.\LAST RUN\curavtura su n','t','xci');
save('.\LAST RUN\curvatura su b','t','k');
save('.\LAST RUN\backbone','t','ubx','uby','ubz');

% Rigid Body
save('.\LAST RUN\rigid long vel','t','vr_x');
save('.\LAST RUN\rigid tras n vel','t','vr_y');
save('.\LAST RUN\rigid tras b vel','t','vr_z');
save('.\LAST RUN\rigid vel rotazione su t','t','wr_x');
save('.\LAST RUN\rigid vel rotazione su n','t','wr_y');
save('.\LAST RUN\rigid vel rotazione su b','t','wr_z');
save('.\LAST RUN\rigid body position','t','ur_x','ur_y','ur_z');
save('.\LAST RUN\Torque beam inertial e2','F_bire2');
save('.\LAST RUN\Torque beam inertial e3','F_bire3');
save('.\LAST RUN\Torque beam external e2','F_bere2');
save('.\LAST RUN\Torque beam external e3','F_bere3');
save('.\LAST RUN\Torque beam e2','F_biare2');
save('.\LAST RUN\Torque beam e3','F_biare3');
save('.\LAST RUN\Pitch rigid','pitch');
save('.\LAST RUN\Axis rigid','axes');
save('.\LAST RUN\Magnitude rigid','magn');
save('.\LAST RUN\Volume','Volu');

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% plots

% Shell

figure
surf(Ts,Ss,va')
grid on
title('vel long')
ylabel('X [m]')
xlabel('t [s]')
zlabel('va [m/s]')
auxstr  =strcat('.\LAST RUN\','velocità longitudinale shell.png');
print('-dpng',auxstr)

figure
surf(Ts,Ss,vb')
grid on
title('vel tras b')
ylabel('X [m]')
xlabel('t [s]')
zlabel('vb [m/s]')
auxstr  =strcat('.\LAST RUN\','velocità tras b shell.png');
print('-dpng',auxstr)

figure
surf(Ts,Ss,we')
grid on
title('rot on -e_{phi}')
ylabel('X [m]')
xlabel('t [s]')
zlabel('w [1/s]')
auxstr  =strcat('.\LAST RUN\','rotazione su -e_{phi} shell.png');
print('-dpng',auxstr)

figure
surf(Ts,Ss,lambda')
grid on
title('longitudinal strain')
ylabel('X [m]')
xlabel('t [s]')
zlabel('lambda [-]')
auxstr  =strcat('.\LAST RUN\','longitudinal_strain shell.png');
print('-dpng',auxstr)

figure
surf(Ts,Ss,eta')
grid on
title('tras b strain')
ylabel('X [m]')
xlabel('t [s]')
zlabel('eta [-]')
auxstr  =strcat('.\LAST RUN\','tras_n_strain shell.png');
print('-dpng',auxstr)

figure
surf(Ts,Ss,mu')
grid on
title('curvature on -e_{phi}')
ylabel('X [m]')
xlabel('t [s]')
zlabel('mu [1/m]')
auxstr  =strcat('.\LAST RUN\','curvature_module_on_menoe_phi shell.png');
print('-dpng',auxstr)

figure
surf(Ts,Ss,erre')
grid on
title('raggio')
ylabel('X [m]')
xlabel('t [s]')
zlabel('r [m]')
auxstr  =strcat('.\LAST RUN\','raggio shell.png');
print('-dpng',auxstr)

figure
surf(Ts,Ss,zeta')
grid on
title('altitudine')
ylabel('X [m]')
xlabel('t [s]')
zlabel('z [m]')
auxstr  =strcat('.\LAST RUN\','altitudine shell.png');
print('-dpng',auxstr)

figure
surf(Ts,Ss,theta')
grid on
title('rotation')
ylabel('X [m]')
xlabel('t [s]')
zlabel('\theta [-]')
auxstr  =strcat('.\LAST RUN\','theta shell.png');
print('-dpng',auxstr)

figure
plot(t,Volu)
grid on
title('Volume')
ylabel('V [m^3]')
xlabel('t [s]')
auxstr  =strcat('.\LAST RUN\','volume.png');
print('-dpng',auxstr)

% Beam

figure
surf(Tb,Sb,vbt')
grid on
title('vel long')
ylabel('X [m]')
xlabel('t [s]')
zlabel('vt [m/s]')
auxstr  =strcat('.\LAST RUN\','velocità longitudinale Beam.png');
print('-dpng',auxstr)

figure
surf(Tb,Sb,vbn')
grid on
title('vel tras n')
ylabel('X [m]')
xlabel('t [s]')
zlabel('vn [m/s]')
auxstr  =strcat('.\LAST RUN\','velocità tras n Beam.png');
print('-dpng',auxstr)

figure
surf(Tb,Sb,vbb')
grid on
title('vel tras b')
ylabel('X [m]')
xlabel('t [s]')
zlabel('vb [m/s]')
auxstr  =strcat('.\LAST RUN\','velocità tras b Beam.png');
print('-dpng',auxstr)

figure
surf(Tb,Sb,wbt')
grid on
title('rot on t')
ylabel('X [m]')
xlabel('t [s]')
zlabel('wt [1/s]')
auxstr  =strcat('.\LAST RUN\','rotazione su t Beam.png');
print('-dpng',auxstr)

figure
surf(Tb,Sb,wbn')
grid on
title('rot on n')
ylabel('X [m]')
xlabel('t [s]')
zlabel('wn [1/s]')
auxstr  =strcat('.\LAST RUN\','rotazione su n Beam.png');
print('-dpng',auxstr)

figure
surf(Tb,Sb,wbb')
grid on
title('rot on b')
ylabel('X [m]')
xlabel('t [s]')
zlabel('wb [1/s]')
auxstr  =strcat('.\LAST RUN\','rotazione su b Beam.png');
print('-dpng',auxstr)

figure
surf(Tb,Sb,q')
grid on
title('longitudinal strain')
ylabel('X [m]')
xlabel('t [s]')
zlabel('q [-]')
auxstr  =strcat('.\LAST RUN\','longitudinal_strain Beam.png');
print('-dpng',auxstr)

figure
surf(Tb,Sb,p')
grid on
title('tras n strain')
ylabel('X [m]')
xlabel('t [s]')
zlabel('p [-]')
auxstr  =strcat('.\LAST RUN\','tras_n_strain Beam.png');
print('-dpng',auxstr)

figure
surf(Tb,Sb,r')
grid on
title('tras b strain')
ylabel('X [m]')
xlabel('t [s]')
zlabel('r [-]')
auxstr  =strcat('.\LAST RUN\','tras_b_strain Beam.png');
print('-dpng',auxstr)

figure
surf(Tb,Sb,tau')
grid on
title('torsion')
ylabel('X [m]')
xlabel('t [s]')
zlabel('tau [1/m]')
auxstr  =strcat('.\LAST RUN\','torsione Beam.png');
print('-dpng',auxstr)

figure
surf(Tb,Sb,xci')
grid on
title('curvature on n')
ylabel('X [m]')
xlabel('t [s]')
zlabel('xci [1/m]')
auxstr  =strcat('.\LAST RUN\','curvature_module_on_n Beam.png');
print('-dpng',auxstr)

figure
surf(Tb,Sb,k')
grid on
title('curvature on b')
ylabel('X [m]')
xlabel('t [s]')
zlabel('k [1/m]')
auxstr  =strcat('.\LAST RUN\','curvature_module_on_b Beam.png');
print('-dpng',auxstr)

% Rigid Body

figure
plot(t,vr_x)
grid on
title('rigid vel x')
xlabel('t [s]')
ylabel('vr_x [m/s]')
auxstr  =strcat('.\LAST RUN\',' rigid velocity x Rigid Body.png');
print('-dpng',auxstr)

figure
plot(t,vr_y)
grid on
title('rigid vel y')
xlabel('t [s]')
ylabel('vr_y [m/s]')
auxstr  =strcat('.\LAST RUN\',' rigid velocity y Rigid Body.png');
print('-dpng',auxstr)

figure
plot(t,vr_z)
grid on
title('rigid vel z')
xlabel('t [s]')
ylabel('vr_z [m/s]')
auxstr  =strcat('.\LAST RUN\',' rigid velocity z Rigid Body.png');
print('-dpng',auxstr)

figure
plot(t,wr_x)
grid on
title('rigid ang vel on x')
xlabel('t [s]')
ylabel('wr_x [1/s]')
auxstr  =strcat('.\LAST RUN\',' rigid angular velocity x Rigid Body.png');
print('-dpng',auxstr)

figure
plot(t,wr_y)
grid on
title('rigid ang vel on y')
xlabel('t [s]')
ylabel('wr_y [1/s]')
auxstr  =strcat('.\LAST RUN\',' rigid angular velocity y Rigid Body.png');
print('-dpng',auxstr)

figure
plot(t,wr_z)
grid on
title('rigid ang vel on z')
xlabel('t [s]')
ylabel('wr_z [1/s]')
auxstr  =strcat('.\LAST RUN\',' rigid angular velocity z Rigid Body.png');
print('-dpng',auxstr)

figure
plot(t,ur_x)
grid on
title('rigid pos x')
xlabel('t [s]')
ylabel('ur_x [m]')
auxstr  =strcat('.\LAST RUN\',' rigid position x Rigid Body.png');
print('-dpng',auxstr)

figure
plot(t,ur_y)
grid on
title('rigid pos y')
xlabel('t [s]')
ylabel('ur_y [m]')
auxstr  =strcat('.\LAST RUN\',' rigid position y Rigid Body.png');
print('-dpng',auxstr)

figure
plot(t,ur_z)
grid on
title('rigid pos z')
xlabel('t [s]')
ylabel('ur_z [m]')
auxstr  =strcat('.\LAST RUN\',' rigid position z Rigid Body.png');
print('-dpng',auxstr)

% Osservabili

figure
plot(t,F_bire2)
grid on
title('Torque beam inertial E2')
ylabel('F_{bi_{E2}} [Nm]')
xlabel('t [s]')
auxstr  =strcat('.\LAST RUN\','Torque inertial e2.png');
print('-dpng',auxstr)

figure
plot(t,F_bire3)
grid on
title('Torque beam inertial E3')
ylabel('F_{bi_{E3}} [Nm]')
xlabel('t [s]')
auxstr  =strcat('.\LAST RUN\','Torque inertial e3.png');
print('-dpng',auxstr)

figure
plot(t,F_bere2)
grid on
title('Torque beam external E2')
ylabel('F_{be_{E2}} [Nm]')
xlabel('t [s]')
auxstr  =strcat('.\LAST RUN\','Torque external e2.png');
print('-dpng',auxstr)

figure
plot(t,F_bere3)
grid on
title('Torque beam external E3')
ylabel('F_{be_{E3}} [Nm]')
xlabel('t [s]')
auxstr  =strcat('.\LAST RUN\','Torque external e3.png');
print('-dpng',auxstr)

figure
plot(t,F_biare2)
grid on
title('Torque beam E2')
ylabel('F_{b_{E2}} [Nm]')
xlabel('t [s]')
auxstr  =strcat('.\LAST RUN\','Torque e2.png');
print('-dpng',auxstr)

figure
plot(t,F_biare3)
grid on
title('Torque beam E3')
ylabel('F_{b_{E3}} [Nm]')
xlabel('t [s]')
auxstr  =strcat('.\LAST RUN\','Torque e3.png');
print('-dpng',auxstr)

figure
plot(t,sradius)
grid on
title('Steering Radius')
ylabel('R_c [m]')
xlabel('t [s]')
auxstr  =strcat('.\LAST RUN\','Steering radius.png');
print('-dpng',auxstr)

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Video Poseidrone

mov               =VideoWriter(strcat('.\LAST RUN\Monopus dynamics'),'MPEG-4');
mov.FrameRate     =10^2;                                    % fps
open(mov)
scrsz             =get(0,'ScreenSize');
figure('color','w','Position',[scrsz(3)/24 2*scrsz(4)/48 11*scrsz(3)/12 9*scrsz(4)/10])

% figure invariants
ang               =linspace(0,2*pi,180);

for ii=1:nsol                               % per ogni istante
    
    % get the solution
    % Shell
    usx_now  =erre(ii,:);
    usy_now  =zeros(1,nsezs);
    usz_now  =zeta(ii,:);
    tsx_now  =costh(ii,:);
    tsy_now  =zeros(1,nsezs);
    tsz_now  =senth(ii,:);
    nsx_now  =-senth(ii,:);
    nsy_now  =zeros(1,nsezs);
    nsz_now  =costh(ii,:);
    bsx_now  =senth(ii,:);
    bsy_now  =-costh(ii,:);
    bsz_now  =zeros(1,nsezs);
    
    % Beam
    ubx_now  =ubx(ii,:);
    uby_now  =uby(ii,:);
    ubz_now  =ubz(ii,:);
    tbx_now  =tbx(ii,:);
    tby_now  =tby(ii,:);
    tbz_now  =tbz(ii,:);
    nbx_now  =nbx(ii,:);
    nby_now  =nby(ii,:);
    nbz_now  =nbz(ii,:);
    bbx_now  =bbx(ii,:);
    bby_now  =bby(ii,:);
    bbz_now  =bbz(ii,:);
    
    % Rigid Body
    urx_now  =ur_x(ii);
    ury_now  =ur_y(ii);
    urz_now  =ur_z(ii);
    trx_now  =tr_x(ii);
    try_now  =tr_y(ii);
    trz_now  =tr_z(ii);
    nrx_now  =nr_x(ii);
    nry_now  =nr_y(ii);
    nrz_now  =nr_z(ii);
    brx_now  =br_x(ii);
    bry_now  =br_y(ii);
    brz_now  =br_z(ii);
    
    % Screw
    uniw_now        =w(:,ii)/norm(w(:,ii));
    magn_now        =magn(ii);
    pitch_now       =pitch(ii);
    
    % get rototraslation
    % Shell
    gs              =zeros(4,nsezs*4);
    gs(1:3,1:4:end) =[tsx_now;tsy_now;tsz_now];   % tb
    gs(1:3,2:4:end) =[nsx_now;nsy_now;nsz_now];   % nb
    gs(1:3,3:4:end) =[bsx_now;bsy_now;bsz_now];   % bb
    gs(1:3,4:4:end) =[usx_now;usy_now;usz_now];   % ub
    gs(4,4:4:end)   =ones(nsezs,1);
    
    % Beam
    gb              =zeros(4,nsezb*4);
    gb(1:3,1:4:end) =[tbx_now;tby_now;tbz_now];   % tb1
    gb(1:3,2:4:end) =[nbx_now;nby_now;nbz_now];   % nb1
    gb(1:3,3:4:end) =[bbx_now;bby_now;bbz_now];   % bb1
    gb(1:3,4:4:end) =[ubx_now;uby_now;ubz_now];   % ub1
    gb(4,4:4:end)   =ones(nsezb,1);
    
    % Rigid Body
    gr              =[trx_now nrx_now brx_now urx_now;...
                      try_now nry_now bry_now ury_now;...
                      trz_now nrz_now brz_now urz_now;...
                      0 0 0 1];
    
    % rototraslo
    % Shell
    g_s             =gr*gsr*gs;
    usx_now         =g_s(1,4:4:end);
    usy_now         =g_s(2,4:4:end);
    usz_now         =g_s(3,4:4:end);
    
    % Beam1234
    g_b             =gr*gbr*gb;
    ubx_now         =g_b(1,4:4:end);
    uby_now         =g_b(2,4:4:end);
    ubz_now         =g_b(3,4:4:end);
    
    % Screw
    axesE_now       =gr*[axes(:,ii);0];
    axesE_now       =axesE_now(1:3);
    uniwE_now       =gr*[uniw_now;0];
    uniwE_now       =uniwE_now(1:3);
   
    clf
    
    % comera fixed to ambient (E1,E2,E3)
    subplot(1,2,1)
    
    % set the graph options
    set(gca,'CameraPosition',[urx_now+(Lr+Lb) ury_now+(Lr+Lb) urz_now+(Lr+Lb)],...           
        'CameraTarget',[urx_now ury_now urz_now],...
        'CameraUpVector',[0 0 1])
    axis equal
    grid on
    hold on
    xlabel('e1 [m]')
    ylabel('e2 [m]')
    zlabel('e3 [m]')
    title(strcat('t= ',num2str(t(ii))))
    axis ([min(ur_x)-(Lr+Lb) max(ur_x)+(Lr+Lb) min(ur_y)-(Lr+Lb) max(ur_y)+(Lr+Lb) min(ur_z)-(Lr+Lb) max(ur_z)+(Lr+Lb)])       

    % disegno il poseidrone
    % Shell
    for jj=1:nsezs
       erre_now=erre(ii,jj);
       zeta_now=zeta(ii,jj);
       sez     =[erre_now*cos(ang);erre_now*sin(ang);zeta_now*ones(1,180);ones(1,180)];
       sez     =gr*gsr*sez;
       plot3(sez(1,:),sez(2,:),sez(3,:),'Color',[0,0,1])
    end
    % profile
    plot3(usx_now,usy_now,usz_now,'LineWidth',2,'Color',[0,0.5,0])
    
     % Beam
    for jj=1:nsezb
       Rnow    =((Rmin-Rmax)/Lb)*Xb(jj) + Rmax;
       sez     =[zeros(1,180);Rnow*sin(ang);Rnow*cos(ang);ones(1,180)];
       sez     =g_b(1:4,4*(jj-1)+1:4*(jj-1)+4)*sez;
       plot3(sez(1,:),sez(2,:),sez(3,:),'Color',[1,0,0])
    end
    % backbone
    plot3(ubx_now,uby_now,ubz_now,'LineWidth',2,'Color',[0,0.5,0])
    
    % force drawing
    drawnow
    
    % comera fixed to reference (e1,e2,e3)
    subplot(1,2,2)
  
    cam_pos         =gr*gsr*[-Lr-Lb;Lr+Lb;Lr+Lb;1];
    set(gca,'CameraPosition',cam_pos(1:3)',...
        'CameraTarget',[urx_now ury_now urz_now],...
        'CameraUpVector',gr(1:3,3)')
    axis equal
    grid on
    hold on
    xlabel('e1 [m]')
    ylabel('e2 [m]')
    zlabel('e3 [m]')
    title(strcat('t= ',num2str(t(ii))))
    axis ([urx_now-(Lr+Lb) urx_now+(Lr+Lb) ury_now-(Lr+Lb) ury_now+(Lr+Lb) urz_now-(Lr+Lb) urz_now+(Lr+Lb)])

    % disegno il poseidrone
    % Shell
    for jj=1:nsezs
       erre_now=erre(ii,jj);
       zeta_now=zeta(ii,jj);
       sez     =[erre_now*cos(ang);erre_now*sin(ang);zeta_now*ones(1,180);ones(1,180)];
       sez     =gr*gsr*sez;
       plot3(sez(1,:),sez(2,:),sez(3,:),'Color',[0,0,1])
    end
    % profile
    plot3(usx_now,usy_now,usz_now,'LineWidth',2,'Color',[0,0.5,0])
    
     % Beam
    for jj=1:nsezb
       Rnow    =((Rmin-Rmax)/Lb)*Xb(jj) + Rmax;
       sez     =[zeros(1,180);Rnow*sin(ang);Rnow*cos(ang);ones(1,180)];
       sez     =g_b(1:4,4*(jj-1)+1:4*(jj-1)+4)*sez;
       plot3(sez(1,:),sez(2,:),sez(3,:),'Color',[1,0,0])
    end
    % backbone
    plot3(ubx_now,uby_now,ubz_now,'LineWidth',2,'Color',[0,0.5,0])
    
    % force drawing
    drawnow
    
    % for movie
    F   = getframe(gcf);
    writeVideo(mov,F);
end

close(mov);

% eof
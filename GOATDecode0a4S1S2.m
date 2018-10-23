dim=20;
Iq=qeye(2);
Is=qeye(dim);
sigmaz=[0,0;0,1]; 
sigmax=[0,1;1,0];
sigmay=[0,-1i;1i,0];
sz1=tensor(sigmaz,Is,Is,Iq);
sx1=tensor(sigmax,Is,Is,Iq);
sy1=tensor(sigmay,Is,Is,Iq);
sz2=tensor(Iq,Is,Is,sigmaz);
sx2=tensor(Iq,Is,Is,sigmax);
sy2=tensor(Iq,Is,Is,sigmay);
tempaD=create(dim);
tempa=destroy(dim);
aD1=tensor(Iq,tempaD,Is,Iq);
a1=tensor(Iq,tempa,Is,Iq);
aD2=tensor(Iq,Is,tempaD,Iq);
a2=tensor(Iq,Is,tempa,Iq);

chisq1=1.759e6*2*pi;
Krs1 = 8.7e3*2*pi;
chisq2=2.777e6*2*pi;
Krs2 = 24.0e3*2*pi;
chiS1S2=15.5e3*2*pi;

amp180=200*2^4;
mysigma=5e-9;
NBAR1=2420;
NBAR2=1210;
CavityPulseWidth=20e-9;
AQ=pi/2/(amp180)/sqrt(2*pi)/mysigma;
AS1=1/(NBAR1*CavityPulseWidth);
AS2=1/(NBAR2*CavityPulseWidth);
H0=-chisq1*sz1*aD1*a1-Krs1/2*aD1*aD1*a1*a1-chisq2*sz2*aD2*a2-Krs2/2*aD2*aD2*a2*a2-chiS1S2*aD1*a1*aD2*a2;
H1=AQ*sx1;
H2=AQ*sy1;
H3=AS1*(aD1+a1);
H4=-AS1*(1i*a1-1i*aD1);
H5=AQ*sx2;
H6=AQ*sy2;
H7=AS2*(aD2+a2);
H8=-AS2*(1i*a2-1i*aD2);
dt=1e-9;
tstart=1e-9;
tstop=500e-9;
num=round((tstop-tstart)/dt)+1;
tspan=linspace(tstart, tstop, 500);

C0 = (tensor(basis(0),fock(dim,0),fock(dim,0),basis(0)));
r0 = (tensor(basis(0),unit(fock(dim,0)+fock(dim,4)),unit(fock(dim,0)+fock(dim,4)),basis(0)));

% r0 = state2dm(tensor(basis(0),fock(dim,0),fock(dim,0),basis(0)));
% C0 = state2dm(tensor(basis(0),unit(fock(dim,0)+fock(dim,4)),fock(dim,0),basis(0)));

C1 = (tensor(unit(basis(0)+basis(1)),fock(dim,0),fock(dim,0),basis(0)));
r1 = (tensor(basis(0),unit(unit(fock(dim,0)+fock(dim,4))+fock(dim,2)),unit(fock(dim,0)+fock(dim,4)),basis(0)));

% r1 = state2dm(tensor(unit(basis(0)+basis(1)),fock(dim,0),fock(dim,0),basis(0)));
% C1 = state2dm(tensor(basis(0),unit(unit(fock(dim,0)+fock(dim,4))+fock(dim,2)),fock(dim,0),basis(0)));

C2 = (tensor(basis(0),fock(dim,0),fock(dim,0),unit(basis(0)+basis(1))));
r2 = (tensor(basis(0),unit(fock(dim,0)+fock(dim,4)),unit(unit(fock(dim,0)+fock(dim,4))+fock(dim,2)),basis(0)));
% 
C3 = (tensor(unit(basis(0)+basis(1)),fock(dim,0),fock(dim,0),unit(basis(0)+basis(1))));
r3 = (tensor(basis(0),unit(unit(fock(dim,0)+fock(dim,4))+fock(dim,2)),unit(unit(fock(dim,0)+fock(dim,4))+fock(dim,2)),basis(0)));
% r2=r0;
% C2=C0;
% 
% r3=r1;
% C3=C1;

% load('uEncode0a4S1Q1.at');
% load('uEncode0a4S2Q2.mat');
% u0=[u,x];
% load('uDecode0a4S1S21013.mat');
% u0=u;

rs={r0,r1,r2,r3};
Cs={C0,C1,C2,C3};%check
Hs={H1,H2,H3,H4,H5,H6,H7,H8};
n_har=35;

% u0=cat( 1, 8000/n_har*(rand(n_har*length(Hs),1)-0.5), pi/(tstop-tstart)*ones(length(Hs),1) );
load('uGOATDecode0a4S1S21013.mat');
u0=u;

penaltyl=1e-7;
h=cat( 1, 8000/n_har*ones(n_har*length(Hs),1), pi/dt/n_har/10*ones(length(Hs),1) );

clear GOATPulse;
Myfun = @(u) GOATPulse(u,H0,Hs,rs,Cs,tspan,n_har,h,penaltyl,'uGOATDecode0a4S1S21013');

options = optimoptions('fminunc');
%% Modify options setting
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'MaxFunctionEvaluations', 1e100);
options = optimoptions(options,'MaxIterations', 1e100);
options = optimoptions(options,'OptimalityTolerance', 1e-100);
options = optimoptions(options,'FunctionTolerance', 1e-100);
options = optimoptions(options,'StepTolerance', 1e-100);
options = optimoptions(options,'PlotFcn', {  @optimplotfval @optimplotstepsize });
options = optimoptions(options,'SpecifyObjectiveGradient', true);
options = optimoptions(options,'Algorithm', 'quasi-newton');
options = optimoptions(options,'HessUpdate', 'steepdesc');

[x,fval] = fminunc(Myfun,u0,options);

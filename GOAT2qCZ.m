dim=2;

dw1 = 0;
dw2 = 0.5e9;
Ec=0.2e9;
g=0.017e9;

Iq=qeye(dim);
tempaD=create(dim);
tempa=destroy(dim);
aD1=tensor(Iq,tempaD);
a1=tensor(Iq,tempa);
aD2=tensor(tempaD,Iq);
a2=tensor(tempa,Iq);

% amp180=200*2^4;
% mysigma=5e-1;
% NBAR1=2420;
% NBAR2=1210;
% CavityPulseWidth=20e-1;
% AQ=pi/2/(amp180)/sqrt(2*pi)/mysigma;
% AS1=1/(NBAR1*CavityPulseWidth);
% AS2=1/(NBAR2*CavityPulseWidth);

H0 = 2*pi*( dw1*aD1*a1 - Ec/2*aD1*aD1*a1*a1 + Ec*aD2*a2 - Ec/2*aD2*aD2*a2*a2 + g*(aD1*a2+aD2*a1) );
% H0=-chisq1*sz1*aD1*a1-Krs1/2*aD1*aD1*a1*a1-chisq2*sz2*aD2*a2-Krs2/2*aD2*aD2*a2*a2-chiS1S2*aD1*a1*aD2*a2;

H1=2*pi*1e9*aD1*a1;
H2=2*pi*1e9*aD2*a2;

dt=1e-9;
tstart=1e-9;
tstop=50e-9;
num=round((tstop-tstart)/dt)+1;
tspan=linspace(tstart, tstop, num);

C0 = (tensor(basis(dim,0),basis(dim,0)));
r0 = (tensor(basis(dim,0),basis(dim,0)));

C1 = (tensor(basis(dim,1),basis(dim,0)));
r1 = (tensor(basis(dim,1),basis(dim,0)));

C2 = (tensor(basis(dim,0),basis(dim,1)));
r2 = (tensor(basis(dim,0),basis(dim,1)));

C3 = (tensor(unit(basis(dim,0)+basis(dim,1)),unit(basis(dim,0)+basis(dim,1))));
r3 = unit(tensor(basis(dim,0),basis(dim,0))+tensor(basis(dim,1),basis(dim,0))+tensor(basis(dim,0),basis(dim,1))-tensor(basis(dim,1),basis(dim,1)));
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
Hs={H1,H2};
n_har=15;

filename = 'uTGOAT2qCZ.mat';
u=cat( 1, 0.1/n_har*(rand(n_har*length(Hs),1)-0.5), pi/(tstop-tstart)/10*ones(length(Hs),1), zeros(length(Hs),1) );
save( filename, 'u');
% load(filename, 'u');
u0=u;

penaltyl=1e-5;
h=cat( 1, 0.1/n_har*ones(n_har*length(Hs),1), pi/dt/n_har*ones(length(Hs),1) );

clear GOATPulse;
Myfun = @(u) GOATPulse(u,H0,Hs,rs,Cs,tspan,n_har,h,penaltyl,filename);

options = optimoptions('fminunc');
%% Modify options setting
% options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'MaxFunctionEvaluations', 1e100);
options = optimoptions(options,'DerivativeCheck', 'off');
options = optimoptions(options,'MaxIterations', 1e100);
options = optimoptions(options,'OptimalityTolerance', 1e-100);
options = optimoptions(options,'FunctionTolerance', 1e-100);
options = optimoptions(options,'StepTolerance', 1e-100);
options = optimoptions(options,'PlotFcn', {  @optimplotfval @optimplotstepsize });
options = optimoptions(options,'SpecifyObjectiveGradient', true);
options = optimoptions(options,'Algorithm', 'quasi-newton');
% options = optimoptions(options,'HessUpdate', 'steepdesc');

[x,fval] = fminunc(Myfun,u0,options);

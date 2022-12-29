%% Levi Dekker 4224175
% Homework set 2 for Multibody Dynamics B

%% Constants etc.
clear all
format short e

% constants
g = 9.81;
l1 = 0.55;
l2 = 0.55;
rho = 1180;
wi = 0.05;
th = 0.004;

m1 = rho*l1*wi*th;
I1 = (1/12)*m1*(l1^2 + wi^2);
m2 = rho*l2*wi*th;
I2 = (1/12)*m2*(l2^2 + wi^2);

% gravity forces
Fg = [m1*g 0 0 m2*g 0 0].'

% mass matrix
M = diag([m1 m1 I1 m2 m2 I2]);

%% This section is copied from hwsym.m by A.L. Schwab, TUDelft 2018
% ADJUSTMENTS: 
% Symbolic variable declaration "real"
% l1 l2 constants

% use the symbolic toolbox to get the derivatives of the constraints
%
% by A.L. Schwab, TUDelft 2018

% define some symbolic variables
syms x1 y1 phi1 x2 y2 phi2 real
syms xd1 yd1 phid1 xd2 yd2 phid2 real % for dx/dt I use xd etc

% put the cm coordinates in a column vector x
x = [x1; y1; phi1; x2; y2; phi2]
% and the time derivaties
xd = [xd1; yd1; phid1; xd2; yd2; phid2]

% the constraints 
dxA = x1-l1/2*cos(phi1);
dyA = y1-l1/2*sin(phi1);
dxB = (x2-l2/2*cos(phi2))-(x1+l1/2*cos(phi1));
dyB = (y2-l2/2*sin(phi2))-(y1+l1/2*sin(phi1));
% put all constraints in one vector C
C = [dxA; dyA; dxB; dyB];

% Cx is the jacobian dC/dx
Cx = jacobian(C,x);
Cx = simplify(Cx)

% and C2 are the convective terms C,xx*xd*xd
% which are by definition d(dC/dt)/dx*xd
% first determine dC/dt
Cd = Cx*xd; % this is dC/dt=dC/dx*xd
% and next the convective terms d(dC/dt)/dx*xd
C2 = jacobian(Cd,x)*xd;
C2 = simplify(C2)

MC = [M Cx'; Cx zeros(4,4)];
FC = [Fg; -C2];


%% Initial conditions 1 (both bars vertical up and zero speeds)
MC1 = subs(MC, phi1, pi/2);
MC1 = subs(MC1, phi2, pi/2);

FC1 = subs(FC, phi1, pi/2);
FC1 = subs(FC1, phi2, pi/2);
FC1 = subs(FC1, phid1, 0);
FC1 = subs(FC1, phid2, 0);

result1 = inv(MC1) * FC1;
result1 = vpa(result1);

%% Initial conditions 2 (both bars horizontal to the right and zero speeds)
MC2 = subs(MC, phi1, 0);
MC2 = subs(MC2, phi2, 0);

FC2 = subs(FC, phi1, 0);
FC2 = subs(FC2, phi2, 0);
FC2 = subs(FC2, phid1, 0);
FC2 = subs(FC2, phid2, 0);

result2 = inv(MC2) * FC2;
result2 = vpa(result2);

%% Initial conditions 3 (both bars horizontal with initial angular speed)
MC3 = subs(MC, phi1, 0);
MC3 = subs(MC3, phi2, 0);

FC3 = subs(FC, phi1, 0);
FC3 = subs(FC3, phi2, 0);
FC3 = subs(FC3, phid1, 2*pi);
FC3 = subs(FC3, phid2, 2*pi);

result3 = inv(MC3) * FC3;
result3 = vpa(result3);


%% Extra constraint
c5 = x2 + 0.5*l2*cos(phi2); % vertical slider constraint through origin
D = [dxA; dyA; dxB; dyB; c5]; % new C-matrix

% Dx is the jacobian 
Dx = jacobian(D,x);
Dx = simplify(Dx)

% D2 convective terms 
Dd = Dx*xd; 
D2 = jacobian(Dd,x)*xd;
D2 = simplify(D2)

MD = [M Dx'; Dx zeros(5,5)];
FD = [Fg; -D2];

%% Question e: both bars vertical up and zero speeds
MDe = subs(MD, phi1, pi/2);
MDe = subs(MDe, phi2, pi/2);

FDe = subs(FD, phi1, pi/2);
FDe = subs(FDe, phi2, pi/2);
FDe = subs(FDe, phid1, 0);
FDe = subs(FDe, phid2, 0);

resulte = inv(MDe) * FDe;
resulte = vpa(resulte);


%% Question f: both bars vertical up and initial angular speed on bar 1
MDf = subs(MD, phi1, pi/2);
MDf = subs(MDf, phi2, pi/2);

FDf = subs(FD, phi1, pi/2);
FDf = subs(FDf, phi2, pi/2);
FDf = subs(FDf, phid1, -2*pi);
FDf = subs(FDf, phid2, 0);

resultf = inv(MDf) * FDf;
resultf = vpa(resultf);


%% Question g: extra force etc.

% Extra force
F = [m1*g 10 (0.5*l1*10) m2*g 0 0].';
FD2 = [F; -D2];

MDg = subs(MD, phi1, 0);
MDg = subs(MDg, phi2, pi);
 
FDg = subs(FD2, phi1, 0);
FDg = subs(FDg, phi2, pi);
FDg = subs(FDg, phid1, 0);
FDg = subs(FDg, phid2, 0);
 
resultg = inv(MDg) * FDg;
resultg = vpa(resultg);

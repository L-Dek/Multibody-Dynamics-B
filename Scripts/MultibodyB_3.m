%% Levi Dekker 4224175
% Homework set 3 for Multibody Dynamics B

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
% THIS SECTION OCCURS MULTIPLE TIMES IN THE SCRIPT
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


%% Spring
l0 = 2*l1/3;
k = (15/2)*m1*g/l1;

Cs = (     (  x1+(1/6)*l1*cos(phi1) + (l1/2)  )^2    +    (  y1 + (1/6)*l1*sin(phi1)  )^2      )^0.5 - l0;   % spring length
Csi = jacobian(Cs, x);
Csi = simplify(Csi);

% Initial conditions
Cs_ev = subs(Cs, x1, 0);
Cs_ev = subs(Cs_ev, y1, 0.5*l1);
Cs_ev = subs(Cs_ev, phi1, pi/2);

Csi_ev = subs(Csi, x1, 0);
Csi_ev = subs(Csi_ev, y1, 0.5*l1);
Csi_ev = subs(Csi_ev, phi1, pi/2);

sigmas = k*Cs_ev;

FC_new = [(Fg - sigmas*Csi_ev'); -C2];

MC1 = subs(MC, phi1, pi/2);
MC1 = subs(MC1, phi2, pi/2);

FC1 = subs(FC_new, phi1, pi/2);
FC1 = subs(FC1, phi2, pi/2);
FC1 = subs(FC1, phid1, 0);
FC1 = subs(FC1, phid2, 0);

result1 = inv(MC1) * FC1;
result1 = vpa(result1);


%testlength = ( (0.55/2)^2+(2*0.55/3)^2 )^0.5
%uitslag = k * (testlength - l0)


%% Question B: Motor
clear all

% constants
omega = (-1)*120*2*pi / 60; % -4pi rad/s

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

% define some symbolic variables
syms x1 y1 phi1 x2 y2 phi2 real
syms xd1 yd1 phid1 xd2 yd2 phid2 real % for dx/dt I use xd etc
syms t real

% put the cm coordinates in a column vector x
x = [x1; y1; phi1; x2; y2; phi2]
% and the time derivaties
xd = [xd1; yd1; phid1; xd2; yd2; phid2]

% the constraints 
dxA = x1-l1/2*cos(phi1);
dyA = y1-l1/2*sin(phi1);
dxB = (x2-l2/2*cos(phi2))-(x1+l1/2*cos(phi1));
dyB = (y2-l2/2*sin(phi2))-(y1+l1/2*sin(phi1));

% electric motor constraint
Cl = phi1 - omega*t; % motor kinematic constraint

% put all constraints in one vector C
C = [dxA; dyA; dxB; dyB; Cl];

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

MC = [M Cx'; Cx zeros(5,5)];
FC = [Fg; -C2];


%Cli = jacobian(Cl,x);


%% Initial conditions: both bars vertical up and omega on both bars
MC2 = subs(MC, phi1, pi/2);
MC2 = subs(MC2, phi2, pi/2);

FC2 = subs(FC, phi1, pi/2);
FC2 = subs(FC2, phi2, pi/2);
FC2 = subs(FC2, phid1, omega);
FC2 = subs(FC2, phid2, omega);

result2 = inv(MC2) * FC2;
result2 = vpa(result2);


%% Question C: Impulses
% Repeat 
clear all

% constants
g = 9.81;
l1 = 0.55;
l2 = 0.55;
rho = 1180;
wi = 0.05;
th = 0.004;
omega = 120 * 2*pi/60;

m1 = rho*l1*wi*th;
I1 = (1/12)*m1*(l1^2 + wi^2);
m2 = rho*l2*wi*th;
I2 = (1/12)*m2*(l2^2 + wi^2);

% gravity forces
Fg = [m1*g 0 0 m2*g 0 0].'

% mass matrix
M = diag([m1 m1 I1 m2 m2 I2]);

syms x1 y1 phi1 x2 y2 phi2 real
syms xd1 yd1 phid1 xd2 yd2 phid2 real % for dx/dt I use xd etc

% put the cm coordinates in a column vector x
x = [x1; y1; phi1; x2; y2; phi2]
% and the time derivaties
xd = [xd1; yd1; phid1; xd2; yd2; phid2]

velocities = [-omega*0.5*l1; 0; omega; -omega*1.5*l1; 0; omega];
xdmin = subs(xd,xd1,velocities(1));
xdmin = subs(xdmin,yd1,velocities(2));
xdmin = subs(xdmin,phid1,velocities(3));
xdmin = subs(xdmin,xd2,velocities(4));
xdmin = subs(xdmin,yd2,velocities(5));
xdmin = subs(xdmin,phid2,velocities(6));

% the constraints 
dxA = x1-l1/2*cos(phi1);
dyA = y1-l1/2*sin(phi1);
dxB = (x2-l2/2*cos(phi2))-(x1+l1/2*cos(phi1));
dyB = (y2-l2/2*sin(phi2))-(y1+l1/2*sin(phi1));
C5 = x2 + 0.5*l2*cos(phi2);
% put all constraints in one vector C
C = [dxA; dyA; dxB; dyB; C5];

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

MC = [M Cx'; Cx zeros(5,5)];

kinmin = 0.5*(xdmin(1)^2 + xdmin(2)^2)*m1 + 0.5*I1*xdmin(3)^2 + 0.5*(xdmin(4)^2 + xdmin(5)^2)*m2 + 0.5*I2*xdmin(6)^2;

%% e = 1
e = 1;
Imp = [M*xdmin; -e*Cx*xdmin];
MC3 = subs(MC, phi1, pi/2);
MC3 = subs(MC3, phi2, pi/2);
Imp = subs(Imp, phi1, pi/2);
Imp = subs(Imp, phi2, pi/2);
Imp = subs(Imp, phid1, omega);
Imp = subs(Imp, phid2, omega);

result3 = inv(MC3) * Imp;
result3 = vpa(result3);

kin1 = 0.5*(result3(1)^2 + result3(2)^2)*m1 + 0.5*I1*result3(3)^2 + 0.5*(result3(4)^2 + result3(5)^2)*m2 + 0.5*I2*result3(6)^2;
diff1 = kinmin - kin1;


%% e = 0.5
e = 0.5;
Imp = [M*xdmin; -e*Cx*xdmin];
MC3 = subs(MC, phi1, pi/2);
MC3 = subs(MC3, phi2, pi/2);
Imp = subs(Imp, phi1, pi/2);
Imp = subs(Imp, phi2, pi/2);
Imp = subs(Imp, phid1, omega);
Imp = subs(Imp, phid2, omega);

result4 = inv(MC3) * Imp;
result4 = vpa(result4);

kin2 = 0.5*(result4(1)^2 + result4(2)^2)*m1 + 0.5*I1*result4(3)^2 + 0.5*(result4(4)^2 + result4(5)^2)*m2 + 0.5*I2*result4(6)^2;
diff2 = kinmin - kin2;

%% e = 0
e = 0;
Imp = [M*xdmin; -e*Cx*xdmin];
MC3 = subs(MC, phi1, pi/2);
MC3 = subs(MC3, phi2, pi/2);
Imp = subs(Imp, phi1, pi/2);
Imp = subs(Imp, phi2, pi/2);
Imp = subs(Imp, phid1, omega);
Imp = subs(Imp, phid2, omega);

result5 = inv(MC3) * Imp;
result5 = vpa(result5);

kin3 = 0.5*(result5(1)^2 + result5(2)^2)*m1 + 0.5*I1*result5(3)^2 + 0.5*(result5(4)^2 + result5(5)^2)*m2 + 0.5*I2*result5(6)^2;
diff3 = kinmin - kin3;

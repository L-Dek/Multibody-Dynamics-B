%% Levi Dekker 4224175
% Homework set 4 for Multibody Dynamics B
% 20-03-2018

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
syms xdd1 ydd1 phidd1 xdd2 ydd2 phidd2 real % 
syms t real;

% put the cm coordinates in a column vector x
x = [x1; y1; phi1; x2; y2; phi2];
% and the time derivaties
xd = [xd1; yd1; phid1; xd2; yd2; phid2];

%% Bar 1 coordinates
x1 = 0.5*l1*cos(phi1);
y1 = 0.5*l1*sin(phi1);
cm1 = [x1; y1];

xd1 = simplify(jacobian(x1,phi1)*phid1);
yd1 = simplify(jacobian(y1,phi1)*phid1);
cmd1 = [xd1; yd1];

xdd1 = simplify(jacobian(xd1,phi1)*phid1 + jacobian(xd1,phid1)*phidd1);
ydd1 = simplify(jacobian(yd1,phi1)*phid1 + jacobian(yd1,phid1)*phidd1);
cmdd1 = [xdd1; ydd1];

%% Bar 2 coordinates
x2 = 2*x1 + 0.5*l1*cos(phi2);
y2 = 2*y1 + 0.5*l1*sin(phi2);
cm2 = [x2; y2];

xd2 = simplify(jacobian(x2,phi1)*phid1 + jacobian(x2,phi2)*phid2);
yd2 = simplify(jacobian(y2,phi1)*phid1 + jacobian(y2,phi2)*phid2);
cmd2 = [xd2; yd2];

xdd2 = simplify(jacobian(xd2,phi1)*phid1 + jacobian(xd2,phid1)*phidd1 + jacobian(xd2,phi2)*phid2 + jacobian(xd2,phid2)*phidd2);
ydd2 = simplify(jacobian(yd2,phi1)*phid1 + jacobian(yd2,phid1)*phidd1 + jacobian(yd2,phi2)*phid2 + jacobian(yd2,phid2)*phidd2);
cmdd2 = [xdd2; ydd2];

%% Kinetic energies
T1 = simplify((1/2) * m1 * (cmd1' * cmd1) + (1/2) * I1 * phid1^2);
T2 = simplify((1/2) * m2 * (cmd2' * cmd2) + (1/2) * I1 * phid2^2);
T = T1 + T2;
T = vpa(simplify(T));

%% Potential energies
V1 = -m1*g*x1;
V2 = -m2*g*x2;
V = V1 + V2;
V = vpa(simplify(V));

%% Generalized coordinates
q = [phi1; phi2];
qd = [phid1; phid2];
qdd = [phidd1; phidd2];

%% Equations of Motion (method from Farbod Alijani from Engineering Dynamics course)

%%%%=======================================================================
%%%% II. Build Lagrange equations
%%%%     D(dT/dqd)/Dt - dT/dq + dV/dq - Q = 0 
%%%%=======================================================================

% - compute dT/dqd
dT_dqd     = simplify(jacobian(T,qd))';  % note: transposed used to adjust dimension

% - compute D(dT/dqd)/Dt
% --> total D/Dt of dT/dqd =  partial d(dT/dqd)/dt 
%                            + d(dT/dqd)/dq*qd
%                            + d(dT/dqd)/dqd*qdd
DdT_Dtdqd =  simplify(jacobian(dT_dqd,t)) + simplify(jacobian(dT_dqd,q)) * qd + simplify(jacobian(dT_dqd,qd)) * qdd;

% - compute dT/dq
dT_dq      = simplify(jacobian(T,q))' ;    
    
% - compute dV/dq
dV_dq      = simplify(jacobian(V,q))' ; 
       
%- Set up equations
Equations = vpa(simplify(DdT_Dtdqd - dT_dq + dV_dq));

T_1a = T
V_1a = V
Equations_1a = Equations

%% HW1 b. Initial Conditions 1: both bars vertical up, zero speeds
q0_1 = [pi/2; pi/2];
qd0_1 = [0; 0];

eqsubs1 = subs(Equations, q, q0_1);
eqsubs1 = subs(eqsubs1, qd, qd0_1);

ic1eq1 = eqsubs1(1) == 0;
ic1eq2 = eqsubs1(2) == 0;

sol1 = solve([ic1eq1, ic1eq2], [phidd1, phidd2]);
%sol1phidd1 = sol1.phidd1;
%sol1phidd2 = sol1.phidd2;
qdd0_1 = [sol1.phidd1; sol1.phidd2];

sol1xdd1 = subs(xdd1, q, q0_1);
sol1xdd1 = subs(sol1xdd1, qd, qd0_1);
sol1xdd1 = subs(sol1xdd1, qdd, qdd0_1);

sol1ydd1 = subs(ydd1, q, q0_1);
sol1ydd1 = subs(sol1ydd1, qd, qd0_1);
sol1ydd1 = subs(sol1ydd1, qdd, qdd0_1);

sol1xdd2 = subs(xdd2, q, q0_1);
sol1xdd2 = subs(sol1xdd2, qd, qd0_1);
sol1xdd2 = subs(sol1xdd2, qdd, qdd0_1);

sol1ydd2 = subs(ydd2, q, q0_1);
sol1ydd2 = subs(sol1ydd2, qd, qd0_1);
sol1ydd2 = subs(sol1ydd2, qdd, qdd0_1);

result1 = vpa([sol1xdd1; sol1ydd1; qdd0_1(1); sol1xdd2; sol1ydd2; qdd0_1(2)])

%% HW1 c. Initial Conditions 2: both bars horizontal to the right and zero speeds
q0_2 = [0; 0];
qd0_2 = [0; 0];

eqsubs2 = subs(Equations, q, q0_2);
eqsubs2 = subs(eqsubs2, qd, qd0_2);

ic2eq1 = eqsubs2(1) == 0;
ic2eq2 = eqsubs2(2) == 0;

sol2 = solve([ic2eq1, ic2eq2], [phidd1, phidd2]);
qdd0_2 = [sol2.phidd1; sol2.phidd2];

sol2xdd1 = subs(xdd1, q, q0_2);
sol2xdd1 = subs(sol2xdd1, qd, qd0_2);
sol2xdd1 = subs(sol2xdd1, qdd, qdd0_2);

sol2ydd1 = subs(ydd1, q, q0_2);
sol2ydd1 = subs(sol2ydd1, qd, qd0_2);
sol2ydd1 = subs(sol2ydd1, qdd, qdd0_2);

sol2xdd2 = subs(xdd2, q, q0_2);
sol2xdd2 = subs(sol2xdd2, qd, qd0_2);
sol2xdd2 = subs(sol2xdd2, qdd, qdd0_2);

sol2ydd2 = subs(ydd2, q, q0_2);
sol2ydd2 = subs(sol2ydd2, qd, qd0_2);
sol2ydd2 = subs(sol2ydd2, qdd, qdd0_2);

result2 = vpa([sol2xdd1; sol2ydd1; qdd0_2(1); sol2xdd2; sol2ydd2; qdd0_2(2)])

%% HW1 d. Initial Conditions 3: both bars horizontal, omega = 60 rpm
omega3 = 60 * 2*pi / 60;

q0_3 = [0; 0];
qd0_3 = [omega3; omega3];

eqsubs3 = subs(Equations, q, q0_3);
eqsubs3 = subs(eqsubs3, qd, qd0_3);

ic3eq1 = eqsubs3(1) == 0;
ic3eq2 = eqsubs3(2) == 0;

sol3 = solve([ic3eq1, ic3eq2], [phidd1, phidd2]);
qdd0_3 = [sol3.phidd1; sol3.phidd2];

sol3xdd1 = subs(xdd1, q, q0_3);
sol3xdd1 = subs(sol3xdd1, qd, qd0_3);
sol3xdd1 = subs(sol3xdd1, qdd, qdd0_3);

sol3ydd1 = subs(ydd1, q, q0_3);
sol3ydd1 = subs(sol1ydd1, qd, qd0_3);
sol3ydd1 = subs(sol1ydd1, qdd, qdd0_3);

sol3xdd2 = subs(xdd2, q, q0_3);
sol3xdd2 = subs(sol3xdd2, qd, qd0_3);
sol3xdd2 = subs(sol3xdd2, qdd, qdd0_3);

sol3ydd2 = subs(ydd2, q, q0_3);
sol3ydd2 = subs(sol3ydd2, qd, qd0_3);
sol3ydd2 = subs(sol3ydd2, qdd, qdd0_3);

result3 = vpa([sol3xdd1; sol3ydd1; qdd0_3(1); sol3xdd2; sol3ydd2; qdd0_3(2)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  SCRIPT CLEARED SECOND TIME IN NEXT SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% New slider constraint through origin is added

%% HW2 e. Extra vertical slide constraint through origin
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

%% This section is copied from hwsym.m by A.L. Schwab, TUDelft 2018

% define some symbolic variables
syms x1 y1 phi1 x2 y2 phi2 real
syms xd1 yd1 phid1 xd2 yd2 phid2 real % for dx/dt I use xd etc
syms xdd1 ydd1 phidd1 xdd2 ydd2 phidd2 real % 
syms t real;

% put the cm coordinates in a column vector x
x = [x1; y1; phi1; x2; y2; phi2];
% and the time derivaties
xd = [xd1; yd1; phid1; xd2; yd2; phid2];

%% Bar 1 coordinates
x1 = 0.5*l1*cos(phi1);
y1 = 0.5*l1*sin(phi1);
cm1 = [x1; y1];

xd1 = simplify(jacobian(x1,phi1)*phid1);
yd1 = simplify(jacobian(y1,phi1)*phid1);
cmd1 = [xd1; yd1];

xdd1 = simplify(jacobian(xd1,phi1)*phid1 + jacobian(xd1,phid1)*phidd1);
ydd1 = simplify(jacobian(yd1,phi1)*phid1 + jacobian(yd1,phid1)*phidd1);
cmdd1 = [xdd1; ydd1];

%% Bar 2 coordinates
x2 = x1; % Slider constraint incorporated
y2 = 3*y1; %2*y1 + 0.5*l1*sin(phi1);
cm2 = [x2; y2];

xd2 = simplify(jacobian(x2,phi1)*phid1);
yd2 = simplify(jacobian(y2,phi1)*phid1);
cmd2 = [xd2; yd2];

xdd2 = simplify(jacobian(xd2,phi1)*phid1 + jacobian(xd2,phid1)*phidd1);
ydd2 = simplify(jacobian(yd2,phi1)*phid1 + jacobian(yd2,phid1)*phidd1);
cmdd2 = [xdd2; ydd2];

%% Kinetic energies
T1 = simplify((1/2) * m1 * (cmd1' * cmd1) + (1/2) * I1 * phid1^2);
T2 = simplify((1/2) * m2 * (cmd2' * cmd2) + (1/2) * I1 * phid1^2); %phid2^2 = phid1^2 ?
T = T1 + T2;
T = vpa(simplify(T));

%% Potential energies
V1 = -m1*g*x1;
V2 = -m2*g*x2;
V = V1 + V2;
V = vpa(simplify(V));

%% Generalized coordinates
q = phi1;
qd = phid1;
qdd = phidd1;

%% Equations of Motion (method from Farbod Alijani from Engineering Dynamics course)

% - compute dT/dqd
dT_dqd     = simplify(jacobian(T,qd))';  % note: transposed used to adjust dimension

% - compute D(dT/dqd)/Dt
DdT_Dtdqd =  simplify(jacobian(dT_dqd,t)) + simplify(jacobian(dT_dqd,q)) * qd + simplify(jacobian(dT_dqd,qd)) * qdd;

% - compute dT/dq
dT_dq      = simplify(jacobian(T,q))' ;    
    
% - compute dV/dq
dV_dq      = simplify(jacobian(V,q))' ; 
       
%- Set up equations
Equations = vpa(simplify(DdT_Dtdqd - dT_dq + dV_dq));

T_2e = T
V_2e = V
Equations_2e = Equations

%% HW2 e. Initial condition 4: Both bars vertical up and zero speeds
q0_4 = pi/2;
qd0_4 = 0;

eqsubs4 = subs(Equations, q, q0_4);
eqsubs4 = subs(eqsubs4, qd, qd0_4);

ic4eq1 = eqsubs4(1) == 0;

sol4 = solve(ic4eq1, phidd1);
qdd0_4 = sol4;

sol4xdd1 = subs(xdd1, q, q0_4);
sol4xdd1 = subs(sol4xdd1, qd, qd0_4);
sol4xdd1 = subs(sol4xdd1, qdd, qdd0_4);

sol4ydd1 = subs(ydd1, q, q0_4);
sol4ydd1 = subs(sol4ydd1, qd, qd0_4);
sol4ydd1 = subs(sol4ydd1, qdd, qdd0_4);

sol4xdd2 = subs(xdd2, q, q0_4);
sol4xdd2 = subs(sol4xdd2, qd, qd0_4);
sol4xdd2 = subs(sol4xdd2, qdd, qdd0_4);

sol4ydd2 = subs(ydd2, q, q0_4);
sol4ydd2 = subs(sol4ydd2, qd, qd0_4);
sol4ydd2 = subs(sol4ydd2, qdd, qdd0_4);

result4 = vpa([sol4xdd1; sol4ydd1; qdd0_4; sol4xdd2; sol4ydd2; -qdd0_4])


%% HW2 f. Initial condition 4: Both bars vertical up with initial omega = 60 rpm on bar 1
q0_5 = pi/2;
qd0_5 = 60 * 2*pi / 60;

eqsubs5 = subs(Equations, q, q0_5);
eqsubs5 = subs(eqsubs5, qd, qd0_5);

ic5eq1 = eqsubs5(1) == 0;

sol5 = solve(ic5eq1, phidd1);
qdd0_5 = sol5;

sol5xdd1 = subs(xdd1, q, q0_5);
sol5xdd1 = subs(sol5xdd1, qd, qd0_5);
sol5xdd1 = subs(sol5xdd1, qdd, qdd0_5);

sol5ydd1 = subs(ydd1, q, q0_5);
sol5ydd1 = subs(sol5ydd1, qd, qd0_5);
sol5ydd1 = subs(sol5ydd1, qdd, qdd0_5);

sol5xdd2 = subs(xdd2, q, q0_5);
sol5xdd2 = subs(sol5xdd2, qd, qd0_5);
sol5xdd2 = subs(sol5xdd2, qdd, qdd0_5);

sol5ydd2 = subs(ydd2, q, q0_5);
sol5ydd2 = subs(sol5ydd2, qd, qd0_5);
sol5ydd2 = subs(sol5ydd2, qdd, qdd0_5);

result5 = vpa([sol5xdd1; sol5ydd1; qdd0_5; sol5xdd2; sol5ydd2; -qdd0_5])


%% CHANGED Potential energies  (FORCE ADDED)
V1 = -m1*g*x1 - 10*2*y1;
V2 = -m2*g*x2;
V = V1 + V2;
V = vpa(simplify(V));

%% Equations of Motion  (New Equations of Motion with added force integrated into potential energy)

% - compute dT/dqd
dT_dqd     = simplify(jacobian(T,qd))';  % note: transposed used to adjust dimension

% - compute D(dT/dqd)/Dt
DdT_Dtdqd =  simplify(jacobian(dT_dqd,t)) + simplify(jacobian(dT_dqd,q)) * qd + simplify(jacobian(dT_dqd,qd)) * qdd;

% - compute dT/dq
dT_dq      = simplify(jacobian(T,q))' ;    
    
% - compute dV/dq
dV_dq      = simplify(jacobian(V,q))' ; 
       
%- Set up equations
Equations = vpa(simplify(DdT_Dtdqd - dT_dq + dV_dq));

T_2g = T
V_2g = V
Equations_2g = Equations


%% HW2 g. Initial conditions 6: Both bars horizontal, zero speeds, vertical force of 10 N at point B
q0_6 = 0;
qd0_6 = 0;

eqsubs6 = subs(Equations, q, q0_6);
eqsubs6 = subs(eqsubs6, qd, qd0_6);

ic6eq1 = eqsubs6(1) == 0;

sol6 = solve(ic6eq1, phidd1);
qdd0_6 = sol6;

sol6xdd1 = subs(xdd1, q, q0_6);
sol6xdd1 = subs(sol6xdd1, qd, qd0_6);
sol6xdd1 = subs(sol6xdd1, qdd, qdd0_6);

sol6ydd1 = subs(ydd1, q, q0_6);
sol6ydd1 = subs(sol6ydd1, qd, qd0_6);
sol6ydd1 = subs(sol6ydd1, qdd, qdd0_6);

sol6xdd2 = subs(xdd2, q, q0_6);
sol6xdd2 = subs(sol6xdd2, qd, qd0_6);
sol6xdd2 = subs(sol6xdd2, qdd, qdd0_6);

sol6ydd2 = subs(ydd2, q, q0_6);
sol6ydd2 = subs(sol6ydd2, qd, qd0_6);
sol6ydd2 = subs(sol6ydd2, qdd, qdd0_6);

result6 = vpa([sol6xdd1; sol6ydd1; qdd0_6; sol6xdd2; sol6ydd2; -qdd0_6])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  SCRIPT CLEARED SECOND TIME IN NEXT SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% HW3 a. Spring added to the model


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

kspring = (15/2)*(m1*g/l1);
l0 = 2*l1/3;


%% This section is copied from hwsym.m by A.L. Schwab, TUDelft 2018
% define some symbolic variables
syms x1 y1 phi1 x2 y2 phi2 real
syms xd1 yd1 phid1 xd2 yd2 phid2 real % for dx/dt I use xd etc
syms xdd1 ydd1 phidd1 xdd2 ydd2 phidd2 real % 
syms t real;

% put the cm coordinates in a column vector x
x = [x1; y1; phi1; x2; y2; phi2];
% and the time derivaties
xd = [xd1; yd1; phid1; xd2; yd2; phid2];

%% Bar 1 coordinates
x1 = 0.5*l1*cos(phi1);
y1 = 0.5*l1*sin(phi1);
cm1 = [x1; y1];

xd1 = simplify(jacobian(x1,phi1)*phid1);
yd1 = simplify(jacobian(y1,phi1)*phid1);
cmd1 = [xd1; yd1];

xdd1 = simplify(jacobian(xd1,phi1)*phid1 + jacobian(xd1,phid1)*phidd1);
ydd1 = simplify(jacobian(yd1,phi1)*phid1 + jacobian(yd1,phid1)*phidd1);
cmdd1 = [xdd1; ydd1];

%% Bar 2 coordinates
x2 = 2*x1 + 0.5*l1*cos(phi2);
y2 = 2*y1 + 0.5*l1*sin(phi2);
cm2 = [x2; y2];

xd2 = simplify(jacobian(x2,phi1)*phid1 + jacobian(x2,phi2)*phid2);
yd2 = simplify(jacobian(y2,phi1)*phid1 + jacobian(y2,phi2)*phid2);
cmd2 = [xd2; yd2];

xdd2 = simplify(jacobian(xd2,phi1)*phid1 + jacobian(xd2,phid1)*phidd1 + jacobian(xd2,phi2)*phid2 + jacobian(xd2,phid2)*phidd2);
ydd2 = simplify(jacobian(yd2,phi1)*phid1 + jacobian(yd2,phid1)*phidd1 + jacobian(yd2,phi2)*phid2 + jacobian(yd2,phid2)*phidd2);
cmdd2 = [xdd2; ydd2];

%% Spring length
pointD = [-l1/2; 0];
pointE = [(2/3)*2*x1; (2/3)*2*y1];
springvec = pointE - pointD;
length = (springvec' * springvec)^0.5;
s = length - l0;

%% Kinetic energies
T1 = simplify((1/2) * m1 * (cmd1' * cmd1) + (1/2) * I1 * phid1^2);
T2 = simplify((1/2) * m2 * (cmd2' * cmd2) + (1/2) * I1 * phid2^2);
T = T1 + T2;
T = vpa(simplify(T));

%% Potential energies
V1 = -m1*g*x1;
V2 = -m2*g*x2;
Vspring = 0.5 * kspring * s^2;
V = V1 + V2 + Vspring;
V = vpa(simplify(V));

%% Generalized coordinates
q = [phi1; phi2];
qd = [phid1; phid2];
qdd = [phidd1; phidd2];

%% Equations of Motion (method from Farbod Alijani from Engineering Dynamics course)

% - compute dT/dqd
dT_dqd     = simplify(jacobian(T,qd))';  % note: transposed used to adjust dimension

% - compute D(dT/dqd)/Dt
DdT_Dtdqd =  simplify(jacobian(dT_dqd,t)) + simplify(jacobian(dT_dqd,q)) * qd + simplify(jacobian(dT_dqd,qd)) * qdd;

% - compute dT/dq
dT_dq      = simplify(jacobian(T,q))' ;    
    
% - compute dV/dq
dV_dq      = simplify(jacobian(V,q))' ; 
       
%- Set up equations
Equations = vpa(simplify(DdT_Dtdqd - dT_dq + dV_dq));

T_3a = T
V_3a = V
Equations_3a = Equations


%% HW3 a. Initial condition 7: (vertical up, zero speeds)
q0_7 = [pi/2; pi/2];
qd0_7 = [0; 0];

eqsubs7 = subs(Equations, q, q0_7);
eqsubs7 = subs(eqsubs7, qd, qd0_7);

ic7eq1 = eqsubs7(1) == 0;
ic7eq2 = eqsubs7(2) == 0;

sol7 = solve([ic7eq1, ic7eq2], [phidd1, phidd2]);
qdd0_7 = [sol7.phidd1; sol7.phidd2];

sol7xdd1 = subs(xdd1, q, q0_7);
sol7xdd1 = subs(sol7xdd1, qd, qd0_7);
sol7xdd1 = subs(sol7xdd1, qdd, qdd0_7);

sol7ydd1 = subs(ydd1, q, q0_7);
sol7ydd1 = subs(sol7ydd1, qd, qd0_7);
sol7ydd1 = subs(sol7ydd1, qdd, qdd0_7);

sol7xdd2 = subs(xdd2, q, q0_7);
sol7xdd2 = subs(sol7xdd2, qd, qd0_7);
sol7xdd2 = subs(sol7xdd2, qdd, qdd0_7);

sol7ydd2 = subs(ydd2, q, q0_7);
sol7ydd2 = subs(sol7ydd2, qd, qd0_7);
sol7ydd2 = subs(sol7ydd2, qdd, qdd0_7);

result7 = vpa([sol7xdd1; sol7ydd1; qdd0_7(1); sol7xdd2; sol7ydd2; qdd0_7(2)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  SCRIPT CLEARED SECOND TIME IN NEXT SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% HW3 b. Motor added to model

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

omega = -120 * 2*pi / 60; % motor speed

%% This section is copied from hwsym.m by A.L. Schwab, TUDelft 2018

% define some symbolic variables
syms x1 y1 phi1 x2 y2 phi2 real
syms xd1 yd1 phid1 xd2 yd2 phid2 real % for dx/dt I use xd etc
syms xdd1 ydd1 phidd1 xdd2 ydd2 phidd2 real % 
syms t real;

% put the cm coordinates in a column vector x
x = [x1; y1; phi1; x2; y2; phi2];
% and the time derivaties
xd = [xd1; yd1; phid1; xd2; yd2; phid2];

%% Bar 1 coordinates
x1 = 0.5*l1*cos(phi1);
y1 = 0.5*l1*sin(phi1);
cm1 = [x1; y1];

xd1 = simplify(jacobian(x1,phi1)*phid1);
yd1 = simplify(jacobian(y1,phi1)*phid1);
cmd1 = [xd1; yd1];

xdd1 = simplify(jacobian(xd1,phi1)*phid1 + jacobian(xd1,phid1)*phidd1);
ydd1 = simplify(jacobian(yd1,phi1)*phid1 + jacobian(yd1,phid1)*phidd1);
cmdd1 = [xdd1; ydd1];

%% Bar 2 coordinates
x2 = 2*x1 + 0.5*l1*cos(phi2);
y2 = 2*y1 + 0.5*l1*sin(phi2);
cm2 = [x2; y2];

xd2 = simplify(jacobian(x2,phi1)*phid1 + jacobian(x2,phi2)*phid2);
yd2 = simplify(jacobian(y2,phi1)*phid1 + jacobian(y2,phi2)*phid2);
cmd2 = [xd2; yd2];

xdd2 = simplify(jacobian(xd2,phi1)*phid1 + jacobian(xd2,phid1)*phidd1 + jacobian(xd2,phi2)*phid2 + jacobian(xd2,phid2)*phidd2);
ydd2 = simplify(jacobian(yd2,phi1)*phid1 + jacobian(yd2,phid1)*phidd1 + jacobian(yd2,phi2)*phid2 + jacobian(yd2,phid2)*phidd2);
cmdd2 = [xdd2; ydd2];


%% Kinetic energies
T1 = simplify((1/2) * m1 * (cmd1' * cmd1) + (1/2) * I1 * phid1^2);
T2 = simplify((1/2) * m2 * (cmd2' * cmd2) + (1/2) * I1 * phid2^2);
T = T1 + T2;
T = vpa(simplify(T));

%% Potential energies
V1 = -m1*g*x1;
V2 = -m2*g*x2;
V = V1 + V2;
V = vpa(simplify(V));

%% Generalized coordinates
q = [phi1; phi2];
qd = [phid1; phid2];
qdd = [phidd1; phidd2];

%% Equations of Motion (method from Farbod Alijani from Engineering Dynamics course)

% MOTOR CONSTRAINT
Cm = phi1 - omega*t;
Cmj = jacobian(Cm,phi1)
syms lambdam real % Lagrange multiplier


% - compute dT/dqd
dT_dqd     = simplify(jacobian(T,qd))';  % note: transposed used to adjust dimension

% - compute D(dT/dqd)/Dt
DdT_Dtdqd =  simplify(jacobian(dT_dqd,t)) + simplify(jacobian(dT_dqd,q)) * qd + simplify(jacobian(dT_dqd,qd)) * qdd;

% - compute dT/dq
dT_dq      = simplify(jacobian(T,q))' ;    
    
% - compute dV/dq
dV_dq      = simplify(jacobian(V,q))' ; 
       
%- Set up equations
Equations = vpa(simplify(DdT_Dtdqd - dT_dq + dV_dq));

T_3b = T
V_3b = V
Equations_3b = Equations

Equations(1) = Equations(1) + Cmj*lambdam;


%% HW3 b. Initial conditions: (Both bars vertical up, initial speeds omega)
q0_8 = [pi/2; pi/2];
qd0_8 = [omega; omega];

eqsubs8 = subs(Equations, q, q0_8);
eqsubs8 = subs(eqsubs8, qd, qd0_8);
eqsubs8 = subs(eqsubs8, phidd1, 0);

ic8eq1 = eqsubs8(1) == 0;
ic8eq2 = eqsubs8(2) == 0;

sol8 = solve([ic8eq1, ic8eq2], [phidd1, phidd2, lambdam]);
qdd0_8 = [sol8.phidd1; sol8.phidd2];

sol8xdd1 = subs(xdd1, q, q0_8);
sol8xdd1 = subs(sol8xdd1, qd, qd0_8);
sol8xdd1 = subs(sol8xdd1, qdd, qdd0_8);

sol8ydd1 = subs(ydd1, q, q0_8);
sol8ydd1 = subs(sol8ydd1, qd, qd0_8);
sol8ydd1 = subs(sol8ydd1, qdd, qdd0_8);

sol8xdd2 = subs(xdd2, q, q0_8);
sol8xdd2 = subs(sol8xdd2, qd, qd0_8);
sol8xdd2 = subs(sol8xdd2, qdd, qdd0_8);

sol8ydd2 = subs(ydd2, q, q0_8);
sol8ydd2 = subs(sol8ydd2, qd, qd0_8);
sol8ydd2 = subs(sol8ydd2, qdd, qdd0_8);

result8 = vpa([sol8xdd1; sol8ydd1; qdd0_8(1); sol8xdd2; sol8ydd2; qdd0_8(2)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  SCRIPT CLEARED SECOND TIME IN NEXT SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% HW 3 c. Impulses
clear all
format short e

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

syms x1 y1 phi1 x2 y2 phi2 real
syms xd1 yd1 phid1 xd2 yd2 phid2 real % for dx/dt I use xd etc
syms xdd1 ydd1 phidd1 xdd2 ydd2 phidd2 real % 
syms t real;

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

%% Bar 1 coordinates
x1 = 0.5*l1*cos(phi1);
y1 = 0.5*l1*sin(phi1);
cm1 = [x1; y1];

xd1 = simplify(jacobian(x1,phi1)*phid1);
yd1 = simplify(jacobian(y1,phi1)*phid1);
cmd1 = [xd1; yd1];

xdd1 = simplify(jacobian(xd1,phi1)*phid1 + jacobian(xd1,phid1)*phidd1);
ydd1 = simplify(jacobian(yd1,phi1)*phid1 + jacobian(yd1,phid1)*phidd1);
cmdd1 = [xdd1; ydd1];

%% Bar 2 coordinates
x2 = 2*x1 + 0.5*l1*cos(phi2);
y2 = 2*y1 + 0.5*l1*sin(phi2);
cm2 = [x2; y2];

xd2 = simplify(jacobian(x2,phi1)*phid1 + jacobian(x2,phi2)*phid2);
yd2 = simplify(jacobian(y2,phi1)*phid1 + jacobian(y2,phi2)*phid2);
cmd2 = [xd2; yd2];

xdd2 = simplify(jacobian(xd2,phi1)*phid1 + jacobian(xd2,phid1)*phidd1 + jacobian(xd2,phi2)*phid2 + jacobian(xd2,phid2)*phidd2);
ydd2 = simplify(jacobian(yd2,phi1)*phid1 + jacobian(yd2,phid1)*phidd1 + jacobian(yd2,phi2)*phid2 + jacobian(yd2,phid2)*phidd2);
cmdd2 = [xdd2; ydd2];

%% Kinetic energies
T1 = simplify((1/2) * m1 * (cmd1' * cmd1) + (1/2) * I1 * phid1^2);
T2 = simplify((1/2) * m2 * (cmd2' * cmd2) + (1/2) * I1 * phid2^2);
T = T1 + T2;
T = vpa(simplify(T));

%% Potential energies
V1 = -m1*g*x1;
V2 = -m2*g*x2;
V = V1 + V2;
V = vpa(simplify(V));

%% Generalized coordinates
q = [phi1; phi2];
qd = [phid1; phid2];
qdd = [phidd1; phidd2];

%% Equations of Motion (method from Farbod Alijani from Engineering Dynamics course)

% Wall constraint
C5 = x2 + 0.5*l2*cos(phi2);
C5j = jacobian(C5,q);


% - compute dT/dqd
dT_dqd     = simplify(jacobian(T,qd))';  % note: transposed used to adjust dimension

% - compute D(dT/dqd)/Dt
DdT_Dtdqd =  simplify(jacobian(dT_dqd,t)) + simplify(jacobian(dT_dqd,q)) * qd + simplify(jacobian(dT_dqd,qd)) * qdd;

% - compute dT/dq
dT_dq      = simplify(jacobian(T,q))' ;    
    
% - compute dV/dq
dV_dq      = simplify(jacobian(V,q))' ; 
       
%- Set up equations
Equations = vpa(simplify(DdT_Dtdqd - dT_dq + dV_dq));

T_3c = T
V_3c = V
Equations_3c = Equations

%% Levi Dekker 4224175
% Homework set 1 for Multibody Dynamics B
%%
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

% symbolic variables
syms x1 y1 phi1 x2 y2 phi2 real; % coordinates
syms x1d y1d phi1d x2d y2d phi2d real; % velocities
syms x1dd y1dd phi1dd x2dd y2dd phi2dd real; % accelerations

% gravity forces
Fg = [m1*g 0 0 m2*g 0 0].'

% matrices
a = [0.5*l1*phi1d^2 * cos(phi1);
    0.5*l1*phi1d^2 * sin(phi1);
    0.5*l1*phi1d^2 * cos(phi1) + 0.5*l2*phi2d^2 * cos(phi2);
    0.5*l1*phi1d^2 * sin(phi1) + 0.5*l2*phi2d^2 * sin(phi2)]

M = diag([m1 m1 I1 m2 m2 I2]);

A = [-1 0 1 0; 
    0 -1 0 1; 
    -0.5*l1*sin(phi1) 0.5*l1*cos(phi1) -0.5*l1*sin(phi1) 0.5*l1*cos(phi1);
    0 0 -1 0;
    0 0 0 -1;
    0 0 -0.5*l2*sin(phi2) 0.5*l2*cos(phi2)]

B = [-1 0 -0.5*l1*sin(phi1) 0 0 0;
    0 -1 0.5*l1*cos(phi1) 0 0 0;
    1 0 -0.5*l1*sin(phi1) -1 0 -0.5*l2*sin(phi2);
    0 1 0.5*l1*cos(phi1) 0 -1 0.5*l2*cos(phi2)]

MAB = [M A; B zeros(4,4)];
Fanda = [Fg; a];

%% Initial conditions 1
MAB1 = subs(MAB, phi1, pi/2);
MAB1 = subs(MAB1, phi2, pi/2);

Fanda1 = subs(Fanda, phi1, pi/2);
Fanda1 = subs(Fanda1, phi2, pi/2);
Fanda1 = subs(Fanda1, phi1d, 0);
Fanda1 = subs(Fanda1, phi2d, 0);

result1 = inv(MAB1) * Fanda1;
result1 = vpa(result1);

%% Initial conditions 2
MAB2 = subs(MAB, phi1, 0);
MAB2 = subs(MAB2, phi2, 0);

Fanda2 = subs(Fanda, phi1, 0);
Fanda2 = subs(Fanda2, phi2, 0);
Fanda2 = subs(Fanda2, phi1d, 0);
Fanda2 = subs(Fanda2, phi2d, 0);

result2 = inv(MAB2) * Fanda2;
result2 = vpa(result2);

%% Initial conditions 3
% 60 rpm = 60 * 2pi rad / 60s = 2pi rad/s
MAB3 = subs(MAB, phi1, 0);
MAB3 = subs(MAB3, phi2, 0);

Fanda3 = subs(Fanda, phi1, 0);
Fanda3 = subs(Fanda3, phi2, 0);
Fanda3 = subs(Fanda3, phi1d, 2*pi);
Fanda3 = subs(Fanda3, phi2d, 2*pi);

result3 = inv(MAB2) * Fanda3;
result3 = vpa(result3);

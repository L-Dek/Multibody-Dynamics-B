% Levi Dekker 4224175
% Homework set 9 for Multibody Dynamics B 07-06-2018

clear all

%% Constants and symbolics
syms alpha beta gamma real
syms alphad betad gammad real
syms alphadd betadd gammadd real

angles = [alpha; beta; gamma];
q = angles;
qd = [alphad; betad; gammad];
qdd = [alphadd; betadd; gammadd];

d = 0.3;
e = 0.4;

md = 3; % Upper arm mass. CM at distance d/3 from shoulder.
me = 3; % Lower arm mass. CM at distance e/2 from elbow.
g = 9.81;

%% Initial position (and other initial conditions)
elbow0 = [0; 0; -d];
hand0 = [0; e; -d] - elbow0; % Relative to elbow!  


q0 = (pi/180)*[110; -20; -20];     % first initial position      
%q0 = (pi/180)*[30; -20; -20];       % second initial position   
%q0 = (pi/180)*[0; 0; -90];
%q0 = [0; 0; 0];
%CHANGE EOM FUNCTION BELOW WHEN USING ONE OF THE Q0'S TO APPLY CORRECT TORQUES!

qd0 = [0; 0; 0];


%% Rotation matrices

% Alpha rotation about x
Ralpha = [1,        0,           0;
    0, cos(alpha), -sin(alpha);
    0, sin(alpha),  cos(alpha)];

% Beta rotation about y
Rbeta = [cos(beta), 0, sin(beta);
    0, 1, 0;
    -sin(beta), 0, cos(beta)];

% Gamma rotation about x
Rgamma = [1,        0,           0;
    0, cos(gamma), -sin(gamma);
    0, sin(gamma),  cos(gamma)];

elbowrot = Ralpha*Rbeta;
handrot = elbowrot*Rgamma;


%% Equations of motion
M = diag([md md md me me me]);

% Mappings
pos1 = elbowrot*elbow0/3;                       %elbow0 = [0; 0; -d];
pos2 = elbowrot*elbow0 + handrot*hand0/2;       %hand0 = [0; e; 0]; relative to elbow

x1 = pos1(1);
y1 = pos1(2);
z1 = pos1(3);

x2 = pos2(1);
y2 = pos2(2);
z2 = pos2(3);

x = [x1; y1; z1; x2; y2; z2];
xd = simplify(jacobian(x,q)*qd);
F = [0; 0; -md*g; 0; 0; -me*g];

% TMT
Tmatrix = jacobian(x,q);
A = simplify(Tmatrix' * M * Tmatrix);
Z = simplify(Tmatrix' * (F - M * jacobian(xd,q)*qd));

eom = A*qdd == Z;


%%%% EOM 2:  Torques included

syms Talpha Tbeta Tgamma
torques = [Talpha; Tbeta; Tgamma];


M_2 = diag([md md md me me me 0 0 0]);
x_2 = [x1; y1; z1; x2; y2; z2; alpha; beta; gamma];
xd_2 = simplify(jacobian(x_2,q)*qd);
F_2 = [0; 0; -md*g; 0; 0; -me*g; Talpha; Tbeta; Tgamma];

% TMT
Tmatrix_2 = jacobian(x_2,q);
A_2 = simplify(Tmatrix_2' * M_2 * Tmatrix_2);
Z_2 = simplify(Tmatrix_2' * (F_2 - M_2 * jacobian(xd_2,q)*qd));

eom2 = A_2*qdd == Z_2;


%% Calculating torques

% First set. Question d.
eom2_0 = subs(eom2,q,q0);
eom2_0 = subs(eom2_0,qd,qd0);

% Initial accelerations. We are interested in finding the torques when the
% arm is held in a static position by the torques. So all accelerations 0.
qdd_0 = [0; 0; 0];

eom2_0 = subs(eom2_0,qdd,qdd_0);
Tvals = solve(eom2_0,torques);

Talpha0 = double(Tvals.Talpha);
Tbeta0 = double(Tvals.Tbeta);
Tgamma0 = double(Tvals.Tgamma);


% Second set. Question e.
q0_2 = (pi/180)*[30; -20; -20];
eom3_0 = subs(eom2,q,q0_2);
eom3_0 = subs(eom3_0,qd,qd0);
eom3_0 = subs(eom3_0,qdd,qdd_0);

Tvals2 = solve(eom3_0,torques);

Talpha0_2 = double(Tvals2.Talpha);
Tbeta0_2 = double(Tvals2.Tbeta);
Tgamma0_2 = double(Tvals2.Tgamma);


%% State space
u = qd;

%q0 already done in previous sections
u0 = qd0;

y = [q; u];
y0 = subs(y,q,q0);
y0 = subs(y0,u,u0);
y0 = double(y0);


%% Numerical integration
time = 5;
h = 0.001;
ydata = zeros(time/h,6);
yold = y0;
time_now = 0;

accel = zeros(time/h,3);

% RK4 scheme 
for i = 1:(time/h) 
    ydata(i,:) = yold(1:6);
    
    tempaccel = eomeval2(yold(1:3),yold(4:6));
    accel(i,:) = tempaccel(4:6);
    
    k1 = eomeval2(yold(1:3),yold(4:6));
    yk1 = yold + 0.5*h*k1;
    k2 = eomeval2(yk1(1:3),yk1(4:6));
    yk2 = yold + 0.5*h*k2;
    k3 = eomeval2(yk2(1:3),yk2(4:6));
    yk3 = yold + h*k3;
    k4 = eomeval2(yk3(1:3),yk3(4:6));
    ynew = yold + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);

    yold = ynew;
end


%% Animation
%{
elbowpos = elbowrot*elbow0;
handpos = elbowpos + handrot*hand0;

for i = 1:length(ydata)
    clf
    
    elbowpos_new = subs(elbowpos,q,ydata(i,1:3)');
    handpos_new = subs(handpos,q,ydata(i,1:3)');
    
    xplot = [0; elbowpos_new(1); handpos_new(1)];
    yplot = [0; elbowpos_new(2); handpos_new(2)];
    zplot = [0; elbowpos_new(3); handpos_new(3)];
    
    plot3(xplot,yplot,zplot,'LineWidth',1.5);
    %plot(xplot,yplot,'LineWidth',1.5);
    %plot(yplot,zplot,'LineWidth',1.5);

    %axis([-0.5 0.5 -1 0.5])
    axis([-1 1 -1 1 -1 1])
    daspect([1 1 1])
    
    
    %hold on
    %cm1 = subs(pos1,q,ydata(i,1:3)');
    %cm2 = subs(pos2,q,ydata(i,1:3)');
    %plot(cm1(2),cm1(3),'-o');
    %plot(cm2(2),cm2(3),'-o');
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    grid on
    
    drawnow
end
%}


%%  Plots
hand_plot = elbow0+hand0;

xplot0 = [0; elbow0(1); hand_plot(1)];
yplot0 = [0; elbow0(2); hand_plot(2)];
zplot0 = [0; elbow0(3); hand_plot(3)];

figure
plot3(xplot0,yplot0,zplot0,'LineWidth',1.5);
hold on
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
daspect([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')

% rotation to ball catch posture
angles1 = (pi/180)*[110; -20; -20];
elbowrot1 = double(subs(elbowrot,angles,angles1));
handrot1 = double(subs(handrot,angles,angles1));

elbow1 = elbowrot1*elbow0;
hand1 = elbow1 + handrot1*hand0;

xplot1 = [0; elbow1(1); hand1(1)];
yplot1 = [0; elbow1(2); hand1(2)];
zplot1 = [0; elbow1(3); hand1(3)];
plot3(xplot1,yplot1,zplot1,'LineWidth',1.5);


%%%% Plotting angular velocities
tijd = linspace(0,time,time/h);
tijd = tijd';

figure
subplot(3,1,1);
plot(tijd,ydata(:,4),'LineWidth',1.5);
title("omega_a_l_p_h_a");
ylabel('(rad/s)')

subplot(3,1,2);
plot(tijd,ydata(:,5),'LineWidth',1.5);
title("omega_b_e_t_a");
ylabel('(rad/s)')

subplot(3,1,3);
plot(tijd,ydata(:,6),'LineWidth',1.5);
title("omega_g_a_m_m_a");
xlabel('time (s)')
ylabel('(rad/s)')


%% Functions
function y = eomeval(q,qd)
alpha = q(1);
beta = q(2);
gamma = q(3);

alphad = qd(1);
betad = qd(2);
gammad = qd(3);

 
A = [ (21*cos(2*beta))/100 - (9*sin(gamma))/50 + (3*cos(gamma)^2)/50 - (9*cos(2*beta)*sin(gamma))/50 - (3*cos(2*beta)*cos(gamma)^2)/50 + 21/100, -(3*cos(gamma)*sin(beta)*(2*sin(gamma) - 3))/50, -(3*cos(beta)*(3*sin(gamma) - 2))/50;
                                                                                            -(3*cos(gamma)*sin(beta)*(2*sin(gamma) - 3))/50,  (3*sin(gamma)^2)/25 - (9*sin(gamma))/25 + 3/10,                                    0;
                                                                                                        -(3*cos(beta)*(3*sin(gamma) - 2))/50,                                               0,                                 3/25];


Z = [ (21*alphad*betad*sin(2*beta))/50 - (2943*cos(beta)*sin(alpha))/250 - (2943*cos(alpha)*cos(gamma))/500 + (3*alphad*gammad*sin(2*gamma))/25 - (9*betad^2*cos(beta)*cos(gamma))/50 + (9*gammad^2*cos(beta)*cos(gamma))/50 + (2943*cos(beta)*sin(alpha)*sin(gamma))/500 + (9*alphad*gammad*cos(beta)^2*cos(gamma))/25 + (6*betad*gammad*cos(gamma)^2*sin(beta))/25 + (3*betad^2*cos(beta)*cos(gamma)*sin(gamma))/25 - (18*alphad*betad*cos(beta)*sin(beta)*sin(gamma))/25 - (6*alphad*betad*cos(beta)*cos(gamma)^2*sin(beta))/25 - (6*alphad*gammad*cos(beta)^2*cos(gamma)*sin(gamma))/25;
                                                                                                                                                                       (2943*cos(alpha)*sin(beta)*sin(gamma))/500 - (2943*cos(alpha)*sin(beta))/250 - (3*betad*gammad*sin(2*gamma))/25 - (21*alphad^2*sin(2*beta))/100 + (9*betad*gammad*cos(gamma))/25 - (6*alphad*gammad*sin(beta))/25 + (6*alphad*gammad*cos(gamma)^2*sin(beta))/25 + (9*alphad^2*cos(beta)*sin(beta)*sin(gamma))/25 + (9*alphad*gammad*sin(beta)*sin(gamma))/25 + (3*alphad^2*cos(beta)*cos(gamma)^2*sin(beta))/25;
                                                                                                                                                                                           (2943*sin(alpha)*sin(gamma))/500 - (3*alphad^2*sin(2*gamma))/50 + (3*betad^2*sin(2*gamma))/50 - (9*betad^2*cos(gamma))/50 - (9*alphad^2*cos(beta)^2*cos(gamma))/50 - (2943*cos(alpha)*cos(beta)*cos(gamma))/500 + (6*alphad*betad*sin(beta))/25 - (6*alphad*betad*cos(gamma)^2*sin(beta))/25 - (9*alphad*betad*sin(beta)*sin(gamma))/25 + (3*alphad^2*cos(beta)^2*cos(gamma)*sin(gamma))/25];
 

B = A\Z;

y = [qd; B]
end


function y = eomeval2(q,qd)
alpha = q(1);
beta = q(2);
gamma = q(3);

alphad = qd(1);
betad = qd(2);
gammad = qd(3);

% Zero set
%Talpha = 0;
%Tbeta = 0;
%Tgamma = 0;

% First set
Talpha = 10.2809;
Tbeta = 1.6126;
Tgamma = 0.1141;

% Second set
%Talpha = 11.2669;
%Tbeta = -4.0831;
%Tgamma = 5.5077;

A = [ (21*cos(2*beta))/100 - (9*sin(gamma))/50 + (3*cos(gamma)^2)/50 - (9*cos(2*beta)*sin(gamma))/50 - (3*cos(2*beta)*cos(gamma)^2)/50 + 21/100, -(3*cos(gamma)*sin(beta)*(2*sin(gamma) - 3))/50, -(3*cos(beta)*(3*sin(gamma) - 2))/50;
                                                                                            -(3*cos(gamma)*sin(beta)*(2*sin(gamma) - 3))/50,  (3*sin(gamma)^2)/25 - (9*sin(gamma))/25 + 3/10,                                    0;
                                                                                                       -(3*cos(beta)*(3*sin(gamma) - 2))/50,                                               0,                                 3/25];
 
Z = [ Talpha - (2943*cos(alpha)*cos(gamma))/500 - (2943*cos(beta)*sin(alpha))/250 + (21*alphad*betad*sin(2*beta))/50 + (3*alphad*gammad*sin(2*gamma))/25 - (9*betad^2*cos(beta)*cos(gamma))/50 + (9*gammad^2*cos(beta)*cos(gamma))/50 + (2943*cos(beta)*sin(alpha)*sin(gamma))/500 + (9*alphad*gammad*cos(beta)^2*cos(gamma))/25 + (6*betad*gammad*cos(gamma)^2*sin(beta))/25 + (3*betad^2*cos(beta)*cos(gamma)*sin(gamma))/25 - (18*alphad*betad*cos(beta)*sin(beta)*sin(gamma))/25 - (6*alphad*betad*cos(beta)*cos(gamma)^2*sin(beta))/25 - (6*alphad*gammad*cos(beta)^2*cos(gamma)*sin(gamma))/25;
                                                                                                                                                                        Tbeta - (21*alphad^2*sin(2*beta))/100 - (2943*cos(alpha)*sin(beta))/250 - (3*betad*gammad*sin(2*gamma))/25 + (2943*cos(alpha)*sin(beta)*sin(gamma))/500 + (9*betad*gammad*cos(gamma))/25 - (6*alphad*gammad*sin(beta))/25 + (6*alphad*gammad*cos(gamma)^2*sin(beta))/25 + (9*alphad^2*cos(beta)*sin(beta)*sin(gamma))/25 + (9*alphad*gammad*sin(beta)*sin(gamma))/25 + (3*alphad^2*cos(beta)*cos(gamma)^2*sin(beta))/25;
                                                                                                                                                                                           Tgamma + (2943*sin(alpha)*sin(gamma))/500 - (3*alphad^2*sin(2*gamma))/50 + (3*betad^2*sin(2*gamma))/50 - (9*betad^2*cos(gamma))/50 - (9*alphad^2*cos(beta)^2*cos(gamma))/50 - (2943*cos(alpha)*cos(beta)*cos(gamma))/500 + (6*alphad*betad*sin(beta))/25 - (6*alphad*betad*cos(gamma)^2*sin(beta))/25 - (9*alphad*betad*sin(beta)*sin(gamma))/25 + (3*alphad^2*cos(beta)^2*cos(gamma)*sin(gamma))/25];
                                                                                                                                                                                       
B = A\Z;

y = [qd; B];
end

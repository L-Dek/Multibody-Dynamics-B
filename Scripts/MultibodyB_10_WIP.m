% Levi Dekker 4224175
% Homework set 10 for Multibody Dynamics B

% Accelerations and lambda (multiplier) vector
% B = [xdd; ydd; zdd; q0dd; q1dd; q2dd; q3dd; lambda];
% coordinates = [x; y; z; q0; q1; q2; q3];

clear all

%% Constants and symbolics
r = 1;                  % radius of sphere
m = 420;                % total mass of vehicle
J = [170,   0,   0;     % mass moment of inertia tensor
       0, 120,   0;
       0,   0, 140];
   
cm = [-0.01; 0.01; -0.1]; % CM pos. relative to sphere center

rho = 1.25;     % air density
A = pi*r^2;     % frontal area of sphere
cd = 0.5;       % drag coefficient
height = 25;    % height of structure
w = 18;         % width of structure
g = 9.81;       % gravity
zeta = 0.1;     % 10% damping ratio

% attachment of bungee cords at [0; r; 0] left and right
attach = [0; r; 0];

syms x y z lam0 lam1 lam2 lam3 vx vy vz Bomegax Bomegay Bomegaz real

X = [x; y; z; lam0; lam1; lam2; lam3; ...
    vx; vy; vz; Bomegax; Bomegay; Bomegaz]; % State vector

% All these variables describe the state and are known initially.
% The equations of motion are evaluated to find the accelerations.
% The next state is then calculated by numerically integrating the
% acclerations.

% Derivatives of state variables (used for the constraint virtual power)
syms lam0d lam1d lam2d lam3d vxd vyd vzd Bomegaxd Bomegayd Bomegazd real
Xd = [vx; vy; vz; lam0d; lam1d; lam2d; lam3d; ... 
    vxd; vyd; vzd; Bomegaxd; Bomegayd; Bomegazd];


Bomega = [Bomegax; Bomegay; Bomegaz];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start with initial accelerations

% integrate one time step
% find state vector

% x_n+1         = x_n + vx_n * h
% lam0_n+1      = lam0_n + lam0d_n * h
% vx_n+1        = vx_n + xdd_n * h;
% Bomegax_n+1   = Bomegax_n + Bomegaxd_n * h;

% all state variables are not found by EOM eval
% only accelerations are found by eomeval

% state space representation for numerical integration
% the state vector X already contains the speeds 
% in the same we did it in previous assignments


% evaluate EoM to find accelerations:


% start with inital config
% state known  body omega's   and  lambda's
% evaluate EOM
% find euler parameter accelerations and xdd ydd zdd etc.

% the in the numerical integration you first calculate
% omegas (from the lambda dots)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Determining k, c_damp and L0  (model simplified to mass-spring system)
amp = 0.5*(height+w/2);                  % amplitude of oscillation
k_single = (m/amp)*4.8*9.81;        % k for theoretical single spring
L0 = m*g/k_single;                  % rest length of bungee cords
c_single = zeta*2*sqrt(k_single*m); % damping coefficient

k = 0.5*k_single;                   % stiffness of one bungee cord
c_damp = 0.5*c_single;              % damping coefficient for one cord


%% Rotation matrices
q0 = lam0;
q1 = lam1;
q2 = lam2;
q3 = lam3;


% USED NOWHERE:
% Body to inertial map
R11 = q0^2 + q1^2 - q2^2 - q3^2;
R12 = 2*(q1*q2 - q0*q3);
R13 = 2*(q1*q3 - q0*q2);
R21 = 2*(q2*q1 - q0*q3);
R22 = q0^2 - q1^2 + q2^2 - q3^2;
R23 = 2*(q2*q3 - q0*q1);
R31 = 2*(q3*q1 - q0*q2);
R32 = 2*(q3*q2 - q0*q1);
R33 = q0^2 - q1^2 - q2^2 + q3^2;

R = [R11 R12 R13;
    R21 R22 R23;
    R31 R32 R33];


% Inertial to body map = R'


%% Euler parameters
q = [q1; q2; q3];
qtilde = [0 -q3 q2;
          q3 0 -q1;
         -q2 q1 0];

Q = ... 
    [q0, -q';
    q, q0*eye(3) + qtilde];

q_quat_dot = simplify((1/2)*Q*[0; Bomega]);
qd = [q_quat_dot(2); q_quat_dot(3); q_quat_dot(4)];


Qbar = ... 
    [q0, -q';
    q, q0*eye(3) - qtilde];
          
quat_inertia = ...
    [0, [0 0 0];
    [0; 0; 0], J];


% Qdot, time derivative of each entry
q0dot = q_quat_dot(1);
q1dot = q_quat_dot(2);
q2dot = q_quat_dot(3);
q3dot = q_quat_dot(4);

qtilded = ... 
[0 -q3dot q2dot;
q3dot 0 -q1dot;
-q2dot q1dot 0];

Qdot = ...
    [q0dot, -qd';
    qd, q0dot*eye(3) + qtilded];


qdot_norm = sqrt(q0dot^2 + qd'*qd);


%% Coordinates 
coordinates = [x; y; z; q0; q1; q2; q3];


%% Constraints (springs)
spherecenter = [x; y; z] - R*[-0.01; 0.01; -0.1];
springleft  =   [0; -w/2; height] - (spherecenter - R*attach);
springright =   [0; w/2; height] - (spherecenter + R*attach);

c1 = norm(springleft) - L0;
cdampjac1 = simplify(jacobian(c1,coordinates));
sigma1 = k*c1 + c_damp*cdampjac1*[vx; vy; vz; q_quat_dot];

c2 = norm(springright) - L0;
cdampjac2 = simplify(jacobian(c2,coordinates));
sigma2 = k*c2 + c_damp*cdampjac2*[vx; vy; vz; q_quat_dot]; 

C = [c1; c2];
sigmas = [sigma1; sigma2];
sigmas = simplify(sigmas);


%% Drag forces
F_drag = -(1/2)*rho*A*cd*sqrt(vx^2+vy^2+vz^2)*[vx;vy;vz];

% Act in x, y, and z direction


%% Equations of motion
M = diag([m m m]);
F = [0; 0; -m*g] + F_drag;

Z1 = F;

M_torques = [0; 0; 0]; % applied torques

A2 = ...
    [4*Q*quat_inertia*Q', 2*[q0; q];
    2*[q0, q'], 0];

Z2 = ...
    [2*Q*[0; M_torques] + 8*Qdot*quat_inertia*Qdot' * [q0; q];
    -2*(qdot_norm)^2];


% Constraint power
C1 = simplify(jacobian(C,coordinates));
Cpower = C1'*sigmas;
Cpower = [Cpower; 0];

Afull = ...
    [M, zeros(3,5);
    zeros(5,3), A2]; 

Zfull = [Z1; Z2];

Afull = simplify(Afull);   
Zfull = simplify(Zfull) - Cpower;   % SIMPLIFYING Cpower FREEZES!

A_eval = matlabFunction(Afull); % @(lam0,lam1,lam2,lam3)
Z_eval = matlabFunction(Zfull); % @(Bomegax,Bomegay,Bomegaz,lam0,lam1,lam2,lam3,vx,vy,vz,x,y,z)

%% Initial conditions
y_0 = [0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0];  % Center of mass at [0; 0; 0];

%% Numerical integration
time = 10;
h = 0.01;
ydata = zeros(time/h,13);
yddata = zeros(time/h,7);
yold = y_0;

% RK4 scheme   
for i = 1:(time/h)
    ydata(i,:) = yold(1:13);
    yddata(i,:) = eomeval(A_eval,Z_eval,yold);
        
    % k_i = [xd yd zd lam0d lam1d lam2d lam3d vxd vyd vzd Bomegaxd Bomegayd Bomegazd]
    
    % k1
    eval1 = eomeval(A_eval,Z_eval,yold);
    lamd1 = lamd(yold(4:7),yold(11:13));
    Bomegad1 = Bomegad(yold(4:7),yold(11:13),eval1(4:7));
    k1 = [yold(8:10); lamd1; eval1(1:3); Bomegad1];
    yk1 = yold + 0.5*h*k1;
    
    % k2        
    eval2 = eomeval(A_eval,Z_eval,yk1);
    lamd2 = lamd(yk1(4:7),yk1(11:13));
    Bomegad2 = Bomegad(yk1(4:7),yk1(11:13),eval2(4:7));
    k2 = [yk1(8:10); lamd2; eval2(1:3); Bomegad2];
    yk2 = yold + 0.5*h*k2;
    
    % k3
    eval3 = eomeval(A_eval,Z_eval,yk2);
    lamd3 = lamd(yk2(4:7),yk2(11:13));
    Bomegad3 = Bomegad(yk2(4:7),yk2(11:13),eval3(4:7));
    k3 = [yk2(8:10); lamd3; eval3(1:3); Bomegad3];
    yk3 = yold + h*k3;   
    
    % k4
    eval4 = eomeval(A_eval,Z_eval,yk3);
    lamd4 = lamd(yk3(4:7),yk3(11:13));
    Bomegad4 = Bomegad(yk3(4:7),yk3(11:13),eval4(4:7));
    k4 = [yk3(8:10); lamd4; eval4(1:3); Bomegad4];
    ynew = yold + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
    
    % Normalize lambda
    % lam_i = lam_i/(sqrt(lam_j * lam_j))
    ynew(4:7) = ynew(4:7)/sqrt(sum(ynew(4:7).*ynew(4:7)));    
    yold = ynew;
end


%% 3D Animation
[sx, sy, sz] = sphere;
figure

for i = 1:length(ydata)
    clf
    
    sSurf = surf(sx+ydata(i,1),sy+ydata(i,2),sz+ydata(i,3));
    
    origin = ydata(i,1:3);
    direction = ydata(i,5:7);
    angle_rad = 2*acos(ydata(i,4));
    angle_deg = 360*angle_rad/(2*pi);
    rotate(sSurf,direction,angle_deg,origin);   % rotates sphere
    
    hold on
    xp = ydata(i,1);
    yp = ydata(i,2);
    zp = ydata(i,3);
    quiver3(xp,yp,zp, 0-xp, w/2-yp, height-zp)
    
    axis([-2 2 -w/2 w/2 0 50])
    daspect([1 1 1])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    drawnow
   
end


%direction = [1 1 1];
%rotate(hSurf,direction,25)

%% Plots
tijd = linspace(0,time,time/h);
tijd = tijd';

plot(tijd,ydata(:,3)); % z-position in time


%% Functions
function yd = eomeval(func1,func2,yvec)
% State X = [x; y; z; lam0; lam1; lam2; lam3; vx; vy; vz; Bomegax; Bomegay; Bomegaz]

% eomeval should give back Xd?

% make it so eomeval gives back accelerations only
% that is: xdd ydd zdd  lamdd's 
x = yvec(1);
y = yvec(2);
z = yvec(3);
lam0 = yvec(4);
lam1 = yvec(5);
lam2 = yvec(6);
lam3 = yvec(7);
vx = yvec(8);      
vy = yvec(9);      
vz = yvec(10);     
Bomegax = yvec(11);                 
Bomegay = yvec(12);
Bomegaz = yvec(13);

%q0 = lam0;
%q1 = lam1;
%q2 = lam2;
%q3 = lam3;
%Bomega = [Bomegax; Bomegay; Bomegaz];

%A_eval = matlabFunction(Afull); % @(lam0,lam1,lam2,lam3)
%Z_eval = matlabFunction(Zfull); % @(Bomegax,Bomegay,Bomegaz,lam0,lam1,lam2,lam3,vx,vy,vz,x,y,z)

A = func1(lam0,lam1,lam2,lam3);
Z = func2(Bomegax,Bomegay,Bomegaz,lam0,lam1,lam2,lam3,vx,vy,vz,x,y,z);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
B = A\Z; % column vector of 8 components
yd = [B(1:7)];
end


% generate Bomegad's from lam0dd lam1dd lam2dd lam3dd
function y = Bomegad(x1,x2,x3)
q0 = x1(1);
q1 = x1(2);
q2 = x1(3);
q3 = x1(4);

Bomega = x2;
qdd = x3;

q = [q1; q2; q3];
qtilde = [0 -q3 q2;
          q3 0 -q1;
         -q2 q1 0];

Q = ... 
    [q0, -q';
    q, q0*eye(3) + qtilde];

q_quat_dot = (1/2)*Q*[0; Bomega];
q0dot = q_quat_dot(1);
qd = [q_quat_dot(2); q_quat_dot(3); q_quat_dot(4)];

qdot_norm_squared = q0dot^2 + qd'*qd;

Bomegad1 = 2*Q'*qdd + 2*[qdot_norm_squared; zeros(3,1)];
Bomegad1 = Bomegad1(2:4);

y = Bomegad1;
end

% generate [lam0d; ... ; lam3d] from omega's
function y = lamd(x1,x2)
q0 = x1(1);
q1 = x1(2);
q2 = x1(3);
q3 = x1(4);

Bomega = x2;

q = [q1; q2; q3];
qtilde = [0 -q3 q2;
          q3 0 -q1;
         -q2 q1 0];

Q = ... 
    [q0, -q';
    q, q0*eye(3) + qtilde];

y = (1/2) * Q * [0; Bomega];
end

%{
% express lam0d ... lam3d  in terms of Bomegax ... Bomegaz
function y = convertlam2ome(x)

end
%}

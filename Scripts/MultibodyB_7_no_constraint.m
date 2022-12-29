%% Levi Dekker 4224175
% Homework set 7 for Multibody Dynamics B
% 23-05-2018

%%
clear all
format short e

%% Constants
global O2A O4B BC O4O2 O4G4 BG5 yc m3 m4 m5 m6 J2 J4 J5

O2A = 0.2;
O4B = 0.7;
BC = 0.6;
O4O2 = 0.3;
O4G4 = 0.4;
BG5 = 0.3;
yc = 0.9;
m3 = 0.5;
m4 = 6;
m5 = 4;
m6 = 2;
J2 = 100; 
J4 = 10;
J5 = 6;
Force = 0;
T = 0;


%% Given initial conditions
omega2_0 = 75 * 2*pi / 60; %theta2dd
theta2_0 = 0;


%% Generalized coordinates
syms theta2 theta4 theta5 real
syms theta2d theta4d theta5d real
syms theta2dd theta4dd theta5dd real

q = [theta2; theta4; theta5];
qd = [theta2d; theta4d; theta5d];
qdd = [theta2dd; theta4dd; theta5dd];


%% Mappings
x3 = O2A * cos(theta2);
y3 = O4O2 + O2A * sin(theta2);
x4 = O4G4 * cos(theta4);
y4 = O4G4 * sin(theta4);
x5 = O4B * cos(theta4) + BG5 * cos(theta5);
y5 = O4B * sin(theta4) + BG5 * sin(theta5);
%x6 = O4B * cos(theta4) + BC * cos(theta5);

x = [theta2; x3; y3; x4; y4; theta4; x5; y5; theta5];
xd = simplify(jacobian(x,q)*qd);


%% Constraints
%c1 = O4B * sin(theta4) + BC * sin(theta5) - yc;
c2 = x3 * sin(theta4) - y3 * cos(theta4);   % y3/x3  =  sin4/cos4  =  tan4
%c2 = y3/x3 - tan(theta4);

Cc = c2; %[c1; c2];
C1 = simplify(jacobian(Cc,q));

C_dd = C1*qd;
C2 = simplify(jacobian(C_dd,q));

C3 = simplify(jacobian(C_dd,qd));


%% TMT-method for equations of motion

% Mass matrix
M = diag([J2,m3,m3,m4,m4,J4,m5,m5,J5]);

% Forces
F = zeros(9,1);
F(3) = -1000;
F(5) = -1000;
F(8) = -1000;

Tmatrix = jacobian(x,q);

%testmatrix = Tmatrix*qd;              % called xd below
%testmatrix2 = jacobian(testmatrix,q);

%A = [Tmatrix' * M * Tmatrix, C1';
%    C1, 0]; %klopt dit wel?

A = [Tmatrix' * M * Tmatrix, C1';
    C3, 0]; % test

Z = [Tmatrix' * (F - M * jacobian(xd,q)*qd);
    -C2*qd];

syms lambda1 lambda2 real
lambda = lambda2; %[lambda1; lambda2];
B = [qdd; lambda];
equations = A * B == Z
%solution = solve(equations,B)                   % To solve!

%solution.theta2dd
%solution.theta4dd
%solution.theta5dd
%solution.lambda1
%solution.lambda2

% B = inv(A) * Z                                 % Or like this!
% evaluate A and Z -> calculate inv(A) * Z      


%% Initial positions and velocities
% Constraints
% c1 = O4B * sin(theta4) + BC * sin(theta5) - yc;
% c2 = x3 * sin(theta4) - y3 * cos(theta4);

% theta4 and theta5 made explicit from constraints
theta4_expl = atan(  (O2A*sin(theta2) + O4O2) / (O2A*cos(theta2))  );
theta5_expl = pi - asin(  (yc - O4B*sin(theta4))/BC   ); % changed to pi - asin

% Initial positions
theta4_0 = vpa(subs(theta4_expl, theta2, theta2_0));
theta5_0 = vpa(subs(theta5_expl, theta4, theta4_0)); 

% Initial velocities
theta4d_0 = jacobian(theta4_expl, theta2) * omega2_0;
theta4d_0 = vpa( subs(theta4d_0, theta2, theta2_0) ); 

theta5d_0 = jacobian(theta5_expl, theta4) * theta4d_0;
theta5d_0 = vpa( subs(theta5d_0, theta4, theta4_0) );

q0 = [theta2_0; theta4_0; theta5_0]
qd0 = [omega2_0; theta4d_0; theta5d_0]
u0 = qd0;


%% Convert equations to state space form and substitute parameters
u = qd;

y = [q; u];
y0 = subs(y,q,q0);
y0 = subs(y0,u,u0);
y0 = double(y0);


%% Numerical integration
time = 10;
h = 0.01;
ydata = zeros(time/h,6);
lamdata = zeros(time/h,2); % Lambda forces
accdata = zeros(time/h,3); % Accelerations
yold = y0;

% RK4 scheme   (with Euler the speed increases because of instability)
for i = 1:(time/h)
    ydata(i,:) = yold(1:6);
        
    [k1,lam] = eomeval5(yold(1),yold(2),yold(3),yold(4),yold(5),yold(6));
    yk1 = yold + 0.5*h*k1;
    accdata(i,:) = k1(4:6); % Extract accelerations for plots
    
    k2 = eomeval5(yk1(1),yk1(2),yk1(3),yk1(4),yk1(5),yk1(6));
    yk2 = yold + 0.5*h*k2;
    
    k3 = eomeval5(yk2(1),yk2(2),yk2(3),yk2(4),yk2(5),yk2(6));
    yk3 = yold + h*k3;
    
    k4 = eomeval5(yk3(1),yk3(2),yk3(3),yk3(4),yk3(5),yk3(6));
    ynew = yold + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
    
    % CORRECT for constraint drift each step
    %ycorrected1 = [projection(ynew(1),ynew(2),ynew(3)); ynew(4); ynew(5); ynew(6)];
    
    % CORRECT velocities
    %ycorrected2 = [ycorrected1(1);ycorrected1(2);ycorrected1(3);velocity(ycorrected1)];
    
    yold = ynew% ycorrected2; %yold = ynew;
    
    %lamdata(i,:) = lam(1);
end


%% Plots and assignment questions
tijd = linspace(0,time,time/h);
tijd = tijd';


% Question d1. The normal force exerted by slider 3 on rocker 4.
%figure;
%plot(tijd,lamdata(:,1),'LineWidth',1.5,'Color','b');
%title("Normal force of slider 3 on rocker 4");
%xlabel('t')
%ylabel('Force (N)')
grid on
hold on


%% Animation
figure;

x3p = O2A * cos(ydata(:,1));
y3p = O4O2 + O2A * sin(ydata(:,1));

xB = O4B * cos(ydata(:,2));
yB = O4B * sin(ydata(:,2));

xC = xB + BC*cos(ydata(:,3));
yC = yB + BC*sin(ydata(:,3));

for i = 1:length(x3p)
    %plot([xstart xend],[ystart yend])
    clf
    hold on
    plot([0, x3p(i)],[O4O2, y3p(i)],[0, xB(i)],[0, yB(i)],[xB(i),xC(i)],[yB(i),yC(i)],'LineWidth',2);
    plot([-1.1, 1.1],[yc, yc]);
    axis([-1.1 1.1 -0.1 1.1])
    daspect([1 1 1])
    drawnow
    %pause(0.01)
end


%% Functions
%{
% EoM evaluation with matrix inversal. Global variables
function [x,y] = eomeval4(theta2,theta4,theta5,theta2d,theta4d,theta5d)
global O2A O4B BC O4O2 O4G4 BG5 m3 m4 m5 m6 J2 J4 J5
                                               
A = [        m3*O2A^2 + J2,                                                                                                         0,                                                                                                         0,               0,                   -O2A*cos(theta2 - theta4);
                         0,                                                        m6*O4B^2*sin(theta4)^2 + m5*O4B^2 + m4*O4G4^2 + J4, (BC*O4B*m6*cos(theta4 - theta5))/2 - (BC*O4B*m6*cos(theta4 + theta5))/2 + BG5*O4B*m5*cos(theta4 - theta5), O4B*cos(theta4), O4O2*sin(theta4) + O2A*cos(theta2 - theta4);
                         0, (BC*O4B*m6*cos(theta4 - theta5))/2 - (BC*O4B*m6*cos(theta4 + theta5))/2 + BG5*O4B*m5*cos(theta4 - theta5),                                                                     m6*BC^2*sin(theta5)^2 + m5*BG5^2 + J5,  BC*cos(theta5),                                           0;
                         0,                                                                                           O4B*cos(theta4),                                                                                            BC*cos(theta5),               0,                                           0;
 -O2A*cos(theta2 - theta4),                                                               O4O2*sin(theta4) + O2A*cos(theta2 - theta4),                                                                                                         0,               0,                                           0];
 
 
 
Z = [                                                                                                                                                                                                      0;
          -O4B*(1000*sin(theta4) + BC*m6*theta5d^2*cos(theta5)*sin(theta4) - BG5*m5*theta5d^2*cos(theta4)*sin(theta5) + BG5*m5*theta5d^2*cos(theta5)*sin(theta4) + O4B*m6*theta4d^2*cos(theta4)*sin(theta4));
 BG5*O4B*m5*theta4d^2*cos(theta5)*sin(theta4) - BC^2*m6*theta5d^2*cos(theta5)*sin(theta5) - BC*O4B*m6*theta4d^2*cos(theta4)*sin(theta5) - BG5*O4B*m5*theta4d^2*cos(theta4)*sin(theta5) - 1000*BC*sin(theta5);
                                                                                                                                                        O4B*sin(theta4)*theta4d^2 + BC*sin(theta5)*theta5d^2;
                                                 - theta4d*(theta4d*(O4O2*cos(theta4) + O2A*sin(theta2 - theta4)) - O2A*theta2d*sin(theta2 - theta4)) - O2A*theta2d*sin(theta2 - theta4)*(theta2d - theta4d)];
 
B = A\Z;
x = [theta2d;theta4d;theta5d;B(1:3)];
y = B(4:5);
end
%}


function [x,y] = eomeval5(theta2,theta4,theta5,theta2d,theta4d,theta5d)



%GRAVITY TEST
 
A = [ cos(theta2)^2/50 + sin(theta2)^2/50 + 100,                                                                 0,                                                                 0,                     -cos(theta2 - theta4)/5;
                                         0,                (73*cos(theta4)^2)/25 + (73*sin(theta4)^2)/25 + 10, (21*cos(theta4)*cos(theta5))/25 + (21*sin(theta4)*sin(theta5))/25, cos(theta2 - theta4)/5 + (3*sin(theta4))/10;
                                         0, (21*cos(theta4)*cos(theta5))/25 + (21*sin(theta4)*sin(theta5))/25,                   (9*cos(theta5)^2)/25 + (9*sin(theta5)^2)/25 + 6,                                           0;
                   -cos(theta2 - theta4)/5,                       cos(theta2 - theta4)/5 + (3*sin(theta4))/10,                                                                 0,                                           0];
 

 
Z = [                                                                                                                                                                                   (cos(theta2)*((sin(theta2)*theta2d^2)/10 - 1000))/5 - (theta2d^2*cos(theta2)*sin(theta2))/50;
 (2*cos(theta4)*((12*sin(theta4)*theta4d^2)/5 - 1000))/5 - (7*sin(theta4)*((14*cos(theta4)*theta4d^2)/5 + (6*cos(theta5)*theta5d^2)/5))/10 + (7*cos(theta4)*((14*sin(theta4)*theta4d^2)/5 + (6*sin(theta5)*theta5d^2)/5 - 1000))/10 - (24*theta4d^2*cos(theta4)*sin(theta4))/25;
                                                                                                       (3*cos(theta5)*((14*sin(theta4)*theta4d^2)/5 + (6*sin(theta5)*theta5d^2)/5 - 1000))/10 - (3*sin(theta5)*((14*cos(theta4)*theta4d^2)/5 + (6*cos(theta5)*theta5d^2)/5))/10;
                                                                                                                    - theta4d*(theta4d*(sin(theta2 - theta4)/5 + (3*cos(theta4))/10) - (theta2d*sin(theta2 - theta4))/5) - (theta2d*sin(theta2 - theta4)*(theta2d - theta4d))/5];
 

                                                                                                            
                                                                                                            
%{
A = [ 
 cos(theta2)^2/50 + sin(theta2)^2/50 + 100,                                                                 0,                                                                 0,                     -cos(theta2 - theta4)/5;
                                         0,                (73*cos(theta4)^2)/25 + (73*sin(theta4)^2)/25 + 10, (21*cos(theta4)*cos(theta5))/25 + (21*sin(theta4)*sin(theta5))/25, cos(theta2 - theta4)/5 + (3*sin(theta4))/10;
                                         0, (21*cos(theta4)*cos(theta5))/25 + (21*sin(theta4)*sin(theta5))/25,                   (9*cos(theta5)^2)/25 + (9*sin(theta5)^2)/25 + 6,                                           0;
                   -cos(theta2 - theta4)/5,                       cos(theta2 - theta4)/5 + (3*sin(theta4))/10,                                                                 0,                                           0];
 


 
                                                                                                                                Z = [                            0;
 (7*cos(theta4)*((14*sin(theta4)*theta4d^2)/5 + (6*sin(theta5)*theta5d^2)/5))/10 - (7*sin(theta4)*((14*cos(theta4)*theta4d^2)/5 + (6*cos(theta5)*theta5d^2)/5))/10;
 (3*cos(theta5)*((14*sin(theta4)*theta4d^2)/5 + (6*sin(theta5)*theta5d^2)/5))/10 - (3*sin(theta5)*((14*cos(theta4)*theta4d^2)/5 + (6*cos(theta5)*theta5d^2)/5))/10;
       - theta4d*(theta4d*(sin(theta2 - theta4)/5 + (3*cos(theta4))/10) - (theta2d*sin(theta2 - theta4))/5) - (theta2d*sin(theta2 - theta4)*(theta2d - theta4d))/5];
%}

B = A\Z;
x = [theta2d;theta4d;theta5d;B(1:3)];
y = B(4);
end


% Coordinate projection
function y = projection(theta2,theta4,theta5)
iter = 0;
tol = 1e-12;
qb_n1 = [theta2;theta4;theta5];
q_n1 = qb_n1;

Cc = [ cos(qb_n1(2))*(sin(qb_n1(1))/5 + 3/10) - (cos(qb_n1(1))*sin(qb_n1(2)))/5];

while (max(abs(Cc))>tol && iter<11)
    Cc = [ cos(qb_n1(2))*(sin(qb_n1(1))/5 + 3/10) - (cos(qb_n1(1))*sin(qb_n1(2)))/5];
    
    C1 = [        cos(qb_n1(1) - qb_n1(2))/5, - cos(qb_n1(1) - qb_n1(2))/5 - (3*sin(qb_n1(2)))/10,                   0];
    
    %C1 = [ -cos(theta2 - theta4)/5, cos(theta2 - theta4)/5 + (3*sin(theta4))/10, 0];
    e = -Cc;
    
    %pseudo = C1' * CC_trans_inv;   
    pseudo = pinv(C1);    
    delta = pseudo*e;
    q_n1 = qb_n1 + delta;
    qb_n1 = q_n1; 
    iter = iter + 1;
    
end

y = q_n1;
end


% Velocity correction
function y = velocity(x) %(theta2,theta4,theta5,theta2d,theta4d,theta5d)
theta2 = x(1);
theta4 = x(2);
theta5 = x(3);
theta2d = x(4);
theta4d = x(5);
theta5d = x(6);

qb_n1 = [theta2;theta4;theta5];
qd_n1 = [theta2d;theta4d;theta5d];
       
%Cd = C1 * qd_n1;
%Cd = [                                 (7*theta4d*cos(theta4))/10 + (3*theta5d*cos(theta5))/5;
%     (theta2d*cos(theta2 - theta4))/5 - theta4d*(cos(theta2 - theta4)/5 + (3*sin(theta4))/10)];

%C1 = [        cos(qb_n1(1) - qb_n1(2))/5, - cos(qb_n1(1) - qb_n1(2))/5 - (3*sin(qb_n1(2)))/10,                   0];
%Cd = C1*qd_n1;

%Cdqd = jacobian(Cd,qd) 
%Cdqd = [                      0,                            (7*cos(theta4))/10, (3*cos(theta5))/5;
%         cos(theta2 - theta4)/5, - cos(theta2 - theta4)/5 - (3*sin(theta4))/10,                 0];

%syms t2 t4 t5
%syms t2d t4d t5d
%qd = [t2d;t4d;t5d];

%C1_syms = [        cos(t2 - t4q)/5, - cos(t2 - t4)/5 - (3*sin(t4))/10,                   0];
%Cd_sysm = C1_syms*qd;
%Cdqd_syms = jacobian(Cd,qd_n1);

Cd = [theta4d*(cos(theta2 - theta4)/5 + (3*sin(theta4))/10) - (theta2d*cos(theta2 - theta4))/5];

Cdqd = [ -cos(theta2 - theta4)/5, cos(theta2 - theta4)/5 + (3*sin(theta4))/10, 0];
G = Cdqd;
e = -Cd;

deltaqd = pinv(G)*e;          % Pseudo-inverse
%deltaqd = G'*inv(G*G')*e;    % slower
%deltaqd = G' * ((G*G')\e);   % slower

y = qd_n1 + deltaqd; % Corrected velocities!


% Tests. test1 and test2 should be 0. Succes!
%test = Cd + G*deltaqd

%C1 = [                      0,                            (7*cos(theta4))/10, (3*cos(theta5))/5;
%       cos(theta2 - theta4)/5, - cos(theta2 - theta4)/5 - (3*sin(theta4))/10,                 0];
%test2 = C1 * y


%{
C1 = [ -cos(theta2 - theta4)/5, cos(theta2 - theta4)/5 + (3*sin(theta4))/10, 0]

Cd:
C1*qd
 
ans = theta4d*(cos(theta2 - theta4)/5 + (3*sin(theta4))/10) - (theta2d*cos(theta2 - theta4))/5
%}


end

%% Levi Dekker 4224175
% Homework set 8 for Multibody Dynamics B
% 29-05-2018

%% Constants and symbolic variables etc.
clear all

syms m1 I1 m2 I2 real
syms a b c d real
syms x1 y1 phi1 x2 y2 phi2 real
syms x1d y1d phi1d x2d y2d phi2d real
syms x1dd y1dd phi1dd x2dd y2dd phi2dd real
syms torque real

constants = [m1; I1; m2; I2; a; b; c; d];
vals = [1; 0.1; 0; 0; 0.5; 0.5; 0.125; 0.125];
X = [x1;y1;phi1;x2;y2;phi2];
Xd = [x1d;y1d;phi1d;x2d;y2d;phi2d];
Xdd = [x1dd;y1dd;phi1dd;x2dd;y2dd;phi2dd];

% Ma = F
M = diag([m1,m1,I1,m2,m2,I2]);
%F = [0;0;0;0;0;0];
F = [0;0;torque;0;0;-torque];

%% Constraints
% Holonomic constraints
c1 = x1 + b*cos(phi1) - x2 + d*cos(phi2);
c2 = y1 + b*sin(phi1) - y2 + d*sin(phi2);

C = [c1;c2];

% Cdot = jacobian(C,X)*Xd = C1*Xd = Cd
% Cdotdot = jacobian(Cd,X)*Xd + jacobian(Cd,Xd)*Xdd
% Cdotdot = C2*Xd + C3*Xdd
% Cdotdot = C2*Xd + C1*Xdd    (C3 = C1, because logic)
% Cdotdot = C2m + C1*Xdd;

C1 = jacobian(C,X);
Cd = C1*Xd;
C2 = jacobian(Cd,X);
C2m = C2*Xd;


% Non-holonomic constraints
% direction of velocity: [cos(phi1); sin(phi1)]
% perpendicular to this: [-sin(phi1); cos(phi1)]
xA = x1 - a*cos(phi1);
xAd = x1d + a*sin(phi1)*phi1d;

yA = y1 - a*sin(phi1);
yAd = y1d - a*cos(phi1)*phi1d;

xC = x2 + c*cos(phi2);
xCd = x2d - c*sin(phi2)*phi2d;

yC = y2 + c*sin(phi2);
yCd = y2d + c*cos(phi2)*phi2d;

s1 = xAd*(-sin(phi1)) + yAd*cos(phi1);
s2 = xCd*(-sin(phi2)) + yCd*cos(phi2);
S = [s1;s2];

% Sdot = jacobian(S,X)*Xd + jacobian(S,Xd)*Xdd  =  0
% Sdot = S1*Xd + S2*Xdd
% Sdot = S1m + S2*Xdd

S1 = jacobian(S,X);
S2 = jacobian(S,Xd);
S1m = S1*Xd;

zeroCC = zeros(2,2);
zeroSC = zeros(2,2);
zeroSS = zeros(2,2);


%% DAE
A = [M, C1', S2';  
    C1, zeroCC, zeroSC;
    S2, zeroSC', zeroSS];

Z = [F; -C2m; -S1m];                          

syms lambda1 lambda2 lambda3 lambda4 real
lambda = [lambda1;lambda2;lambda3;lambda4];
B = [Xdd; lambda];

DAE = A*B == Z;


%% Initial positions and velocities
eqn_c1 = x1 + b*cos(phi1) - x2 + d*cos(phi2) == 0;
eqn_c2 = y1 + b*sin(phi1) - y2 + d*sin(phi2) == 0;
eqn_c1d = x1d - x2d - b*phi1d*sin(phi1) - d*phi2d*sin(phi2) == 0;
eqn_c2d = y1d - y2d + b*phi1d*cos(phi1) + d*phi2d*cos(phi2) == 0;
eqn_s1 = xAd*(-sin(phi1)) + yAd*cos(phi1) == 0;
eqn_s2 = xCd*(-sin(phi2)) + yCd*cos(phi2) == 0;

x1_0 = 0;
y1_0 = 0;
phi1_0 = 0;
x2_0 = x1_0 + 0.5 - 0.125; 
y2_0 = 0; 
phi2_0 = pi;
x1d_0 = 0;
y1d_0 = 0;
phi1d_0 = 0;
x2d_0 = 0;
y2d_0 = 0;
phi2d_0 = 0; 

y_0 = [x1_0; y1_0; phi1_0; x2_0; y2_0; phi2_0; x1d_0; y1d_0; phi1d_0; x2d_0; y2d_0; phi2d_0];

%% State space
u = Xd;
y = [X; u];

%% Numerical integration
time = 100;
h = 0.001;
ydata = zeros(time/h,12);
yold = y_0;
time_now = 0;
torque_data = zeros(time/h,1);

% RK4 scheme 
for i = 1:(time/h)
    torque_now = -0.1*cos(pi*time_now);
    torque_data(i) = torque_now;
    
    ydata(i,:) = yold(1:12);
    k1 = eomeval(yold,torque_now);
    yk1 = yold + 0.5*h*k1;
    k2 = eomeval(yk1,torque_now);
    yk2 = yold + 0.5*h*k2;
    k3 = eomeval(yk2,torque_now);
    yk3 = yold + h*k3;
    k4 = eomeval(yk3,torque_now);
    ynew = yold + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
    
    % CORRECT for constraint drift each step
    ycorrected1 = [projection(ynew(1:6)); ynew(7:12)];
    
    % CORRECT velocities
    ycorrected2 = [ycorrected1(1:6); velocity(ycorrected1(1:6),ycorrected1(7:12))];
    
    yold = ycorrected2;
    
    time_now = time_now + h;
    progress = time_now/time
end


%% Animation

ap = 0.5;
bp = 0.5;
cp = 0.125;
dp = 0.125;


%{
x1p = ydata(:,1);
y1p = ydata(:,2);
phi1p = ydata(:,3);

x2p = ydata(:,4);
y2p = ydata(:,5);
phi2p = ydata(:,6);

fast = 20;

figure;

for i = 1:(length(x1p)/fast)  %% Divided by 10 for faster animation
    clf
    
    %{
    hold on
    %plot([xstart xend],[ystart yend])
    plot([x1p(i), x1p(i) + bp*cos(phi1p(i))],[y1p(i), y1p(i) + bp*sin(phi1p(i))],'LineWidth',2,'Color','r');   
    plot([x1p(i), x1p(i) - bp*cos(phi1p(i))],[y1p(i), y1p(i) - bp*sin(phi1p(i))],'LineWidth',2,'Color','r'); 
        
    plot([x2p(i), x2p(i) + cp*cos(phi2p(i))],[y2p(i), y2p(i) + cp*sin(phi2p(i))],'LineWidth',2,'Color','b');        
    plot([x2p(i), x2p(i) - dp*cos(phi2p(i))],[y2p(i), y2p(i) - dp*sin(phi2p(i))],'LineWidth',2,'Color','b');
         
    axis([-2 2 -2 2])
    daspect([1 1 1])
    drawnow
    %pause(0.5)
    %}

    hold on
    %plot([xstart xend],[ystart yend])
    %plot([x1p(10*i), x1p(10*i) + bp*cos(phi1p(10*i))],[y1p(10*i), y1p(10*i) + bp*sin(phi1p(10*i))],'LineWidth',2,'Color','r');   
    %plot([x1p(10*i), x1p(10*i) - bp*cos(phi1p(10*i))],[y1p(10*i), y1p(10*i) - bp*sin(phi1p(10*i))],'LineWidth',2,'Color','r'); 
        
    %plot([x2p(10*i), x2p(10*i) + cp*cos(phi2p(10*i))],[y2p(10*i), y2p(10*i) + cp*sin(phi2p(10*i))],'LineWidth',2,'Color','b');        
    %plot([x2p(10*i), x2p(10*i) - dp*cos(phi2p(10*i))],[y2p(10*i), y2p(10*i) - dp*sin(phi2p(10*i))],'LineWidth',2,'Color','b');
         
    plot([x1p(fast*i), x1p(fast*i) + bp*cos(phi1p(fast*i))],[y1p(fast*i), y1p(fast*i) + bp*sin(phi1p(fast*i))],'LineWidth',2,'Color','r');   
    plot([x1p(fast*i), x1p(fast*i) - bp*cos(phi1p(fast*i))],[y1p(fast*i), y1p(fast*i) - bp*sin(phi1p(fast*i))],'LineWidth',2,'Color','r'); 
        
    plot([x2p(fast*i), x2p(fast*i) + cp*cos(phi2p(fast*i))],[y2p(fast*i), y2p(fast*i) + cp*sin(phi2p(fast*i))],'LineWidth',2,'Color','b');        
    plot([x2p(fast*i), x2p(fast*i) - dp*cos(phi2p(fast*i))],[y2p(fast*i), y2p(fast*i) - dp*sin(phi2p(fast*i))],'LineWidth',2,'Color','b');
         
    
    axis([-1 12 -1 12])
    daspect([1 1 1])
    drawnow
      
end
%}


%% Assignments and plots
tijd = linspace(0,time,time/h);
tijd = tijd';


% Question f.
figure
grid on
hold on

% Path of point A
xA_plot = ydata(:,1) - ap*cos(ydata(:,3));   %xA = x1 - a*cos(phi1);
yA_plot = ydata(:,2) - ap*cos(ydata(:,3));   %yA = y1 - a*sin(phi1);
plot(xA_plot,yA_plot,'LineWidth',1.1,'Color','r')

% Path of point C
xC_plot = ydata(:,4) + cp*cos(ydata(:,6));   %xC = x2 + c*cos(phi2);
yC_plot = ydata(:,5) + cp*sin(ydata(:,6));   %yC = y2 + c*sin(phi2);
plot(xC_plot,yC_plot,'LineWidth',1.1,'Color','b')
title("Path of point A and point B");
xlabel('x position (m)')
ylabel('y position (m)')
legend('point A','point C')


% Question g.
phi1d_plot = ydata(:,9);
phi2d_plot = ydata(:,12);

v_1 = sqrt(ydata(:,7).^2 + ydata(:,8).^2);
v_2 = sqrt(ydata(:,10).^2 + ydata(:,11).^2);

figure
subplot(2,1,1);
plot(tijd,v_1,'LineWidth',1.5,'Color','r');
hold on
plot(tijd,v_2,'LineWidth',1.5,'Color','b');
title("Linear speed of CM's");
xlabel('time (s)')
ylabel('speed (m/s)')
grid on
legend('v1','v2')

subplot(2,1,2);
plot(tijd,phi2d_plot,'LineWidth',1.5,'Color','r');
hold on
plot(tijd,phi1d_plot,'LineWidth',1.5,'Color','b');
title("Angular velocities");
xlabel('time (s)')
ylabel('angular velocity (rad/s)')
grid on
legend('phi2d','phi1d')


figure
subplot(2,1,1);
hold on
plot(tijd,ydata(:,7),'LineWidth',1.5,'Color','r');
plot(tijd,ydata(:,8),'LineWidth',1.5,'Color','b');
title("x and y velocities of body 1 CM");
xlabel('time (s)')
ylabel('speed (m/s)')
grid on
legend('x1d','y1d')

subplot(2,1,2);
hold on
plot(tijd,ydata(:,10),'LineWidth',1.5,'Color','g');
plot(tijd,ydata(:,11),'LineWidth',1.5,'Color','m');
title("x and y velocities of body 2 CM");
xlabel('time (s)')
ylabel('speed (m/s)')
grid on
legend('x2d','y2d')


% Question h.

% constants = [m1;  I1; m2; I2;   a;   b;     c;     d];
% vals =      [ 1; 0.1;  0;  0; 0.5; 0.5; 0.125; 0.125];
% E_kin = 0.5*1*(ydata(:,7).^2 + ydata(:,8).^2)  +  0.5*0.1*ydata(:,9).^2;  
E_kin = 0.5*1*v_1.^2  +  0.5*0.1*ydata(:,9).^2;  

figure
plot(tijd,E_kin,'LineWidth',2.4,'Color','r')
title("Kinetic energy and work done");
xlabel('time (s)')
ylabel('Energy and Work (J)')
grid on
hold on


% Question i.
torque1 = torque_data(1:end-1);                                 % torque applied to body 1  (2:end)?
torque2 = -torque_data(1:end-1);                                % torque applied to body 2
dphi1 = diff(ydata(:,3));                                       % increments of phi1
dphi2 = diff(ydata(:,6));                                       % increments of phi2
dwork = torque1.*dphi1 + torque2.*dphi2;                        % increments of work
workdone = zeros(size(dwork,1),1);                              % total work done at every point in time

for i = 1:(size(dwork))
    workdone(i) = sum(dwork(1:i));                              % Running sum of work increments (integral)
end

plot(tijd(2:end),workdone,'LineWidth',0.8,'Color','y') 
legend('Kinetic energy','Work done')


%% Functions

function yd = eomeval(y,torque)
m1 = 1;
I1 = 0.1;
m2 = 0;
I2 = 0;
a = 0.5;
b = 0.5;
c = 0.125;
d = 0.125;

% X = [x1;y1;phi1;x2;y2;phi2];
% Xd = [x1d;y1d;phi1d;x2d;y2d;phi2d];
% y = [X; Xd];

x1 = y(1);
y1 = y(2);
phi1 = y(3);
x2 = y(4);
y2 = y(5);
phi2 = y(6);
x1d = y(7);
y1d = y(8);
phi1d = y(9);
x2d = y(10);
y2d = y(11);
phi2d = y(12);


A = [    m1,         0,                               0,          0,         0,                             0,            1,           0,                      -sin(phi1),                             0;
          0,        m1,                               0,          0,         0,                             0,            0,           1,                       cos(phi1),                             0;
          0,         0,                              I1,          0,         0,                             0, -b*sin(phi1), b*cos(phi1), - a*sin(phi1)^2 - a*cos(phi1)^2,                             0;
          0,         0,                               0,         m2,         0,                             0,           -1,           0,                               0,                    -sin(phi2);
          0,         0,                               0,          0,        m2,                             0,            0,          -1,                               0,                     cos(phi2);
          0,         0,                               0,          0,         0,                            I2, -d*sin(phi2), d*cos(phi2),                               0, c*cos(phi2)^2 + c*sin(phi2)^2;
          1,         0,                    -b*sin(phi1),         -1,         0,                  -d*sin(phi2),            0,           0,                               0,                             0;
          0,         1,                     b*cos(phi1),          0,        -1,                   d*cos(phi2),            0,           0,                               0,                             0;
 -sin(phi1), cos(phi1), - a*sin(phi1)^2 - a*cos(phi1)^2,          0,         0,                             0,            0,           0,                               0,                             0;
          0,         0,                               0, -sin(phi2), cos(phi2), c*cos(phi2)^2 + c*sin(phi2)^2,            0,           0,                               0,                             0];
      
      
      
Z =    [                                                                          0;
                                                                                  0;
                                                                                  torque;
                                                                                  0;
                                                                                  0;
                                                                                  -torque;
                                                b*cos(phi1)*phi1d^2 + d*cos(phi2)*phi2d^2;
                                                b*sin(phi1)*phi1d^2 + d*sin(phi2)*phi2d^2;
        phi1d*(cos(phi1)*(x1d + a*phi1d*sin(phi1)) + sin(phi1)*(y1d - a*phi1d*cos(phi1)));
        phi2d*(cos(phi2)*(x2d - c*phi2d*sin(phi2)) + sin(phi2)*(y2d + c*phi2d*cos(phi2)))];

B = A\Z;
yd = [y(7:12); B(1:6)];
end


% Coordinate projection
function y = projection(X)

%constants = [m1; I1; m2; I2; a; b; c; d];
vals = [1; 0.1; 0; 0; 0.5; 0.5; 0.125; 0.125];
m1 = vals(1);
I1 = vals(2);
m2 = vals(3);
I2 = vals(4);
a = vals(5);
b = vals(6);
c = vals(7);
d = vals(8);

x1 = X(1);
y1 = X(2);
phi1 = X(3);
x2 = X(4);
y2 = X(5);
phi2 = X(6);

iter = 0;
tol = 1e-12;
qb_n1 = X;
q_n1 = qb_n1;
 
Cc = [x1 - x2 + b*cos(phi1) + d*cos(phi2)
      y1 - y2 + b*sin(phi1) + d*sin(phi2)];

% C1 = [ 1, 0, -b*sin(phi1), -1,  0, -d*sin(phi2);
%        0, 1,  b*cos(phi1),  0, -1,  d*cos(phi2)];

while (max(abs(Cc))>tol && iter<11)
    Cc = [qb_n1(1) - qb_n1(4) + b*cos(qb_n1(3)) + d*cos(qb_n1(6));
          qb_n1(2) - qb_n1(5) + b*sin(qb_n1(3)) + d*sin(qb_n1(6))];
    
    C1 = [ 1, 0, -b*sin(qb_n1(3)), -1,  0, -d*sin(qb_n1(6));
           0, 1,  b*cos(qb_n1(3)),  0, -1,  d*cos(qb_n1(6))];
       
    e = -Cc;
    pseudo = pinv(C1);
    delta = pseudo*e;
    q_n1 = qb_n1 + delta;
    qb_n1 = q_n1;
    iter = iter + 1;
end
y = q_n1;
end


% Velocity correction
function yd = velocity(X,Xd) 

% constants = [m1; I1; m2; I2; a; b; c; d];
vals = [1; 0.1; 0; 0; 0.5; 0.5; 0.125; 0.125];
m1 = vals(1);
I1 = vals(2);
m2 = vals(3);
I2 = vals(4);
a = vals(5);
b = vals(6);
c = vals(7);
d = vals(8);


x1 = X(1);
y1 = X(2);
phi1 = X(3);
x2 = X(4);
y2 = X(5);
phi2 = X(6);

x1d = Xd(1);
y1d = Xd(2);
phi1d = Xd(3);
x2d = Xd(4);
y2d = Xd(5);
phi2d = Xd(6);

qd_n1 = [Xd];


% Cd = 0
% Cd_bar + jacobian(Cd_bar,Xd) * delta_Xd = 0

% S = 0
% S_bar + jacobian(S,Xd) * delta_Xd = 0

% [Cd; S] = 0;
% [Cd_bar; S_bar] + jacobian([Cdr; S],Xd) * delta_Xd = 0

% New system    ( jacobian([Cd; S], Xd) = K )
% delta_q + mu*K  =  0
% K*delta_q       = -[Cd_bar; S_bar]

% pinv( jacobian(K,Xd) )

%C1 = jacobian(C,X);
C1 = [ 1, 0, -b*sin(phi1), -1,  0, -d*sin(phi2);
       0, 1,  b*cos(phi1),  0, -1,  d*cos(phi2)];  
     
%Cd = C1 * qd_n1;
Cd = [x1d - x2d - b*phi1d*sin(phi1) - d*phi2d*sin(phi2);
      y1d - y2d + b*phi1d*cos(phi1) + d*phi2d*cos(phi2)];
      
S = [cos(phi1)*(y1d - a*phi1d*cos(phi1)) - sin(phi1)*(x1d + a*phi1d*sin(phi1));
     cos(phi2)*(y2d + c*phi2d*cos(phi2)) - sin(phi2)*(x2d - c*phi2d*sin(phi2))];

% K = jacobian([Cd; S], Xd);
K = [     1,         0,                    -b*sin(phi1),         -1,         0,                  -d*sin(phi2);
          0,         1,                     b*cos(phi1),          0,        -1,                   d*cos(phi2);
 -sin(phi1), cos(phi1), - a*sin(phi1)^2 - a*cos(phi1)^2,          0,         0,                             0;
          0,         0,                               0, -sin(phi2), cos(phi2), c*cos(phi2)^2 + c*sin(phi2)^2];
 
e = -[Cd; S];

deltaqd = pinv(K)*e;
      
      
yd = qd_n1 + deltaqd; % Corrected velocities!


% Test to see if the correction works as it's supposed to

%{
x1d = yd(1);
y1d = yd(2);
phi1d = yd(3);
x2d = yd(4);
y2d = yd(5);
phi2d = yd(6);

Cd_TEST = [x1d - x2d - b*phi1d*sin(phi1) - d*phi2d*sin(phi2);
           y1d - y2d + b*phi1d*cos(phi1) + d*phi2d*cos(phi2)]
S_TEST = [cos(phi1)*(y1d - a*phi1d*cos(phi1)) - sin(phi1)*(x1d + a*phi1d*sin(phi1));
          cos(phi2)*(y2d + c*phi2d*cos(phi2)) - sin(phi2)*(x2d - c*phi2d*sin(phi2))]
%}
end

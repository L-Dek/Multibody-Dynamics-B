%% Levi Dekker 4224175
% Homework set 6 for Multibody Dynamics B
% 01-05-2018

%% Constants etc.
clear all
format short e

% constants
g_val = 9.81;
l1_val = 0.55;
l2_val = 0.55;
rho_val = 1180;
wi_val = 0.05;
th_val = 0.004;
m1_val = rho_val*l1_val*wi_val*th_val;
I1_val = (1/12)*m1_val*(l1_val^2 + wi_val^2);
m2_val = rho_val*l2_val*wi_val*th_val;
I2_val = (1/12)*m2_val*(l2_val^2 + wi_val^2);

syms g l1 l2 rho wi th m1 I1 m2 I2 real
para_sym = [g l1 l2 rho wi th m1 I1 m2 I2];
para_val = [g_val l1_val l2_val rho_val wi_val th_val m1_val I1_val m2_val I2_val];

%Time = 3;
Time = 11;

%q0 = [pi/2; pi/2];
q0 = [0.2*pi; 0.2*pi];
u0 = [0; 0];


%% This section is copied from hwsym.m by A.L. Schwab, TUDelft 2018
% define symbolic variables
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


%% Convert equations to state space form and substitute parameters
u = qd;

sol1 = solve(Equations==0,qdd);
ud = [sol1.phidd1; sol1.phidd2];
ud = simplify(ud);

yd = [u; ud];
yd = vpa(subs(yd,para_sym,para_val)); %evaluate for parameters like g, m1, m2


%num = subs(evalpara,q,[0;0]);
%num = subs(num,qd,[0;0]);

y = [q; u];
y0 = subs(y,q,q0);
y0 = subs(y0,u,u0);
y0 = double(y0);


%% Euler (Global error calculation)
%yold = y0;

%{
nmin = 6;
nmax = 17;

resvec = zeros(nmax-nmin+1,1);
hvec = zeros(nmax-nmin+1,1);
E1 = zeros(1,1);
ydata = zeros(1,4);

% No symbolic toolbox
for n = nmin:nmax
h = Time/(2^n);

    yold = y0;
    
    for i = 1:(2^n)
        %ydata(i,:) = yold(1:4);
        ydold = eomeval(yold(1),yold(2),yold(3),yold(4));
        ynew = yold + h*ydold;
        yold = ynew; 
    end
    
    resvec(n-nmin+1) = yold(2)
    hvec(n-nmin+1) = Time/2^n;
      
end

hvec = hvec(2:(nmax-nmin+1));

for i = 1:(nmax-nmin)
    E1(i) = abs(resvec(i+1) - resvec(i))
end

loglog(hvec,E1,'-o')
grid on
hold on

%}

% n = 12;
% h = Time/(2^n);
% 
% for i = 1:(2^n)
%     ydata(i,:) = yold(1:4);
%     ydold = eomeval(yold(1),yold(2),yold(3),yold(4));
%     ynew = yold + h*ydold;
%     yold = ynew; 
% end


%% Heun

nmin = 9;
nmax = 9;

resvec2 = zeros(nmax-nmin+1,1);
hvec2 = zeros(nmax-nmin+1,1);
E2 = zeros(1,1);
ydata2 = zeros(1,4);

% Heun scheme
for n = nmin:nmax
h = Time/(2^n);

    yold = y0;
    
    for i = 1:(2^n)
        ydata(i,:) = yold(1:4);
        ydold = eomeval(yold(1),yold(2),yold(3),yold(4));
        
        ystar = yold + h*ydold;
        ydstar = eomeval(ystar(1),ystar(2),ystar(3),ystar(4));
        ynew = yold + 0.5*h*(ydold + ydstar);
        
        yold = ynew; 
    end
    
    resvec2(n-nmin+1) = yold(2)
    hvec2(n-nmin+1) = Time/2^n;
      
end

hvec2 = hvec2(2:(nmax-nmin+1));

for i = 1:(nmax-nmin)
    E2(i) = abs(resvec2(i+1) - resvec2(i))
end

%{
loglog(hvec,E2,'-o')
grid on
hold on
%}


%% ANIMATION
%{
x1 = l1_val * cos(ydata(:,1));
y1 = l1_val * sin(ydata(:,1));
x2 = x1 + l2_val * cos(ydata(:,2));
y2 = y1 + l2_val * sin(ydata(:,2));
%}

x1 = l1_val * cos(ydata(:,1) - pi/2);
y1 = l1_val * sin(ydata(:,1) - pi/2);
x2 = x1 + l2_val * cos(ydata(:,2) - pi/2);
y2 = y1 + l2_val * sin(ydata(:,2) - pi/2);

for i = 1:length(x1)
    %plot([xstart xend],[ystart yend])
    plot([0 x1(i)],[0, y1(i)],[x1(i) x2(i)],[y1(i) y2(i)],'LineWidth',2);
    axis([-1.1 1.1 -1.1 1.1])
    daspect([1 1 1])
    drawnow
    pause(0.01)
end


%% RK4
%{
resvec3 = zeros(nmax-nmin+1,1);
hvec3 = zeros(nmax-nmin+1,1);
E3 = zeros(1,1);
ydata3 = zeros(1,4);


% RK4 scheme
for n = nmin:nmax
h = Time/(2^n);

    yold = y0;
    
    for i = 1:(2^n)
        k1 = eomeval(yold(1),yold(2),yold(3),yold(4));     
        yk1 = yold + 0.5*h*k1;
        
        k2 = eomeval(yk1(1),yk1(2),yk1(3),yk1(4));
        yk2 = yold + 0.5*h*k2;
        
        k3 = eomeval(yk2(1),yk2(2),yk2(3),yk2(4));
        yk3 = yold + h*k3;
        
        k4 = eomeval(yk3(1),yk3(2),yk3(3),yk3(4));
        ynew = yold + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
        
        yold = ynew; 
    end
    
    resvec3(n-nmin+1) = yold(2)
    hvec3(n-nmin+1) = Time/2^n;
      
end

hvec3 = hvec3(2:(nmax-nmin+1));

for i = 1:(nmax-nmin)
    E3(i) = abs(resvec3(i+1) - resvec3(i))
end

loglog(hvec,E3,'-o')
grid on
hold on

%}


%% Built in ODE functions
function x = eomeval(phi1,phi2,phid1,phid2)
x = zeros(4,1);
x(1) = phid1;
x(2) = phid2;
x(3) = (1.1*(0.002803092655*sin(2.0*phi1 - 2.0*phi2)*phid1^2 + 0.003745178891666666977657147441505*sin(phi1 - 1.0*phi2)*phid2^2 + 0.049996979901*sin(phi1 - 2.0*phi2) + 0.15040413788400001664081790764271*sin(phi1)))/(0.0030834019205*cos(2.0*phi1 - 2.0*phi2) - 0.0079081306751944455856852693484724);
x(4) = -(0.14278*(0.11523535833333333572925383236907*sin(phi1 - 1.0*phi2)*phid1^2 + 0.021595475*sin(2.0*phi1 - 2.0*phi2)*phid2^2 - 0.89982552000000004273450926461919*sin(phi2) + 1.155554235*sin(2.0*phi1 - 1.0*phi2)))/(0.0030834019205*cos(2.0*phi1 - 2.0*phi2) - 0.0079081306751944455856852693484724);
end

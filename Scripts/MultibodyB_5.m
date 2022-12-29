%% Levi Dekker 4224175
% Homework set 5 for Multibody Dynamics B

clear all

%% Constants and symbolic variables
syms l m M I g gamma real

syms alpha beta real
syms alphad betad real
syms alphadd betadd real

syms t real

q = [alpha; beta];
qd = [alphad; betad];
qdd = [alphadd; betadd];

Massvec = [M, M, m, m, I, m, m, I];
Mass = diag(Massvec)

F = [0; -M*g; 0; -m*g; 0; 0; -m*g; 0];

%% Coordinates and velocites
% Generalised coordinates are alpha and beta

% Rotation matrix because of incline with angle gamma:
C = [cos(gamma) sin(gamma); -sin(gamma) cos(gamma)];

% Positions of CoM's in rotated frame
Mpos_prime = [l*sin(alpha); l*cos(alpha)];
Apos_prime = (1/2) * Mpos_prime;
Bpos_prime = Mpos_prime + (1/2)*l*[sin(-beta); -cos(-beta)];

% Positions of CoM's in corrected (gravity pointing downwards) frame
Mpos = C*Mpos_prime;
Apos = C*Apos_prime;
Bpos = C*Bpos_prime;

x = [Mpos(1); Mpos(2); Apos(1); Apos(2); -alpha; Bpos(1); Bpos(2); -beta]; % -alpha, -beta ?
xd = simplify(jacobian(x,q)*qd);
xdd = simplify(jacobian(xd,q)*qd + jacobian(xd,qd)*qdd);

% Velocites
Mposd = simplify(jacobian(Mpos,alpha)*alphad);
Aposd = simplify(jacobian(Apos,alpha)*alphad);
Bposd = simplify(jacobian(Bpos,alpha)*alphad + jacobian(Bpos,beta)*betad);


%% Kinetic energies
%Leg A
TA = (1/2)*m*Aposd'*Aposd + (1/2)*I*alphad^2;

%Leg B
TB = (1/2)*m*Bposd'*Bposd + (1/2)*I*betad^2;

%Mass M
TM = (1/2)*M*Mposd'*Mposd;

T = TA + TB + TM;

%% Potential energies
%Leg A
VA = m*g*Apos(2);

%Leg B
VB = m*g*Bpos(2);

%Mass M
VM = M*g*Mpos(2);

V = VA + VB + VM;

%% Derivation of EoM using Lagrange (method from Farbod Alijani from Engineering Dynamics course)

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
Equations = simplify(DdT_Dtdqd - dT_dq + dV_dq);

T
V
Equations


%% TMT-method
Tmatrix = jacobian(x,q)

TMTeqn = simplify(Tmatrix' * Mass * Tmatrix * qdd == Tmatrix' * (F - Mass * jacobian(xd,q)*qd))

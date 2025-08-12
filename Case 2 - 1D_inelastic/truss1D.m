function [fe,me,hsv] = truss1D(xl,ul,A,E,rho,fy,H,hsv)

% function to conduct the explicit state determination of a 2-node truss

% INPUT:
%
% xl:  2-component vector, containing x1 and x2
% ul:  2-component vector, containing x1 and x2
% A:    cross-sectional area
% E:    elastic modulus
% rho: mass density
% fy: initial yield stress
% H: hardening modulus
% hsv: history variable vector: hsv(1) = stress of previous step, hsv(2) =
% strain of previous step, hsv(3) = effective plastic strain

% OUTPUT:
%
% fe:   internal force vector
% me:   mass matrix

le = xl(2)-xl(1);

eps = (ul(2)-ul(1))/le;

sig0 = hsv(1);
eps0 = hsv(2);
epeff = hsv(3);

sig = E*(eps-eps0)+sig0;    % trial elastic stress
sy = fy + H*epeff;

if sy < 0
   sy = 0; 
end

if abs(sig) > sy
    r1 = abs(sig)/sig;
    dlam = (abs(sig) - sy)/(E+H);
    sy = sy + H*dlam;
    sig = r1*sy;
    epeff = epeff + dlam;
end

hsv(1) = sig;
hsv(2) = eps;


hsv(3) = epeff;


%dt = 2.e-5;
%er = (eps-eps0)/dt;

%Cv = 0.0000001*E;

%fe(1) = -(sig+Cv*er)*A;
%fe(2) = (sig+Cv*er)*A;

fe(1) = -sig*A;
fe(2) = sig*A;

me(1) = rho*A*le/2;
me(2) = me(1);



end
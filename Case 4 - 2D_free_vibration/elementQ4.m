function [fe,me]=elementQ4(xy,E,v,rho,uv)

% Routine to calculate the lumped mass matrix {me}, 
% and the right-hand-side vector, {fe}, of a 4Q element, for explicit
% analysis

% INPUT:
%
% xy:   (4x2) array containing the nodal coordinates in PHYSICAL SPACE
% E: elastic modulus
% rho: density
% v: Poisson ratio
% uv:   nodal solution (8x1)
%
%
%
%
% OUTPUT:
%
% me:   element lumped mass matrix (as a column vector)
% fe:   element RHS vector

ksi=zeros(4,1);
eta=zeros(4,1);

me=zeros(8,1);
fe=zeros(8,1); 

ksi(1)=-((1./3.)^0.5); eta(1)=-((1./3.)^0.5);
ksi(2)=((1./3.)^0.5); eta(2)=-((1./3.)^0.5);
ksi(3)=((1./3.)^0.5); eta(3)=((1./3.)^0.5); 
ksi(4)=-((1./3.)^0.5); eta(4)=((1./3.)^0.5); 
mv = zeros(4,1);

Dmat=E/(1-v*v)*[1 v 0;v 1 0;0 0 0.5*(1.-v)];

for i1=1:4   % loop over the 4 Gauss points
    %fprintf('\n\n***  Quadrature point %d  ***\n',i1);
    N4Q=[0.25*(1-ksi(i1))*(1-eta(i1)) 0.25*(1+ksi(i1))*(1-eta(i1)) 0.25*(1+ksi(i1))*(1+eta(i1)) 0.25*(1-ksi(i1))*(1+eta(i1))];
    %fprintf('Derivatives of shape functions with parametric coordinates:\n');
    DN4Q=[-0.25*(1-eta(i1)) 0.25*(1-eta(i1)) 0.25*(1+eta(i1)) -0.25*(1+eta(i1))
          -0.25*(1-ksi(i1)) -0.25*(1+ksi(i1)) 0.25*(1+ksi(i1)) 0.25*(1-ksi(i1))];

    % calculate coordinates of G.P. in physical space:
      %fprintf('Physical coordinates for quadrature point\n');
      %xi=N4Q*xy(:,1)
      %yi=N4Q*xy(:,2)
      
      %    Calculate Jacobian Array: 
    Jmat=zeros(2,2);
    for k1=1:2
        for k2=1:2
    Jmat(k1,k2)=DN4Q(k1,:)*xy(:,k2);
        end
    end
    %fprintf('Jacobian Array:\n');
    %Jmat
    %fprintf('Jacobian Determinant:\n');
    J=det(Jmat);   % Jacobian Determinant
    Jinv=inv(Jmat); % Inverse of the Jacobian Array
     
    DNxmat=Jinv*DN4Q;  % 2x4 array containing in the 1st row the derivatives of the shape functions with x, and in the 2nd row the derivatives of the shape functions with y 
                
            %fprintf('B-matrix:\n');
    Bmat=[DNxmat(1,1) 0 DNxmat(1,2) 0 DNxmat(1,3) 0 DNxmat(1,4) 0; 
          0 DNxmat(2,1) 0 DNxmat(2,2) 0 DNxmat(2,3) 0 DNxmat(2,4);
          DNxmat(2,1) DNxmat(1,1) DNxmat(2,2) DNxmat(1,2) DNxmat(2,3) DNxmat(1,3) DNxmat(2,4) DNxmat(1,4)];            % [B] array  

  
    eps = Bmat*uv;
    sig = Dmat*eps;
      
    mv = mv + N4Q'*rho*J;
    fe = fe + Bmat'*sig*J;
      
end
for i1 = 1:4
    me(2*i1-1) = mv(i1);
    me(2*i1) = mv(i1);
end 

L0 = 0;     % set this to the reference size used for the computation of dt_cr!

if L0 > 0
% mass scaling (used for irregular meshes)

Jm0 = zeros(2,2);

Jm0(1,1) = -0.25*xy(1,1)+0.25*xy(2,1)+0.25*xy(3,1)-0.25*xy(4,1);
Jm0(1,2) = -0.25d0*xy(1,2)+0.25d0*xy(2,2)+0.25d0*xy(3,2)-0.25d0*xy(4,2);
Jm0(2,1) = -0.25d0*xy(1,1)-0.25d0*xy(2,1)+0.25d0*xy(3,1)+0.25d0*xy(4,1);
Jm0(2,2) = -0.25d0*xy(1,2)-0.25d0*xy(2,2)+0.25d0*xy(3,2)+0.25d0*xy(4,2);
      
      
Jdet = Jm0(1,1)*Jm0(2,2)-Jm0(1,2)*Jm0(2,1);
elsiz = sqrt(4.0*Jdet);
     
               
xksi = zeros(2,1);
xeta = xksi;
          
for jj = 1:2
    xksi(jj) = Jm0(1,jj);
    xeta(jj) = Jm0(2,jj);
end 
          
          
mag1 = 0.0;
for jj = 1:2
    mag1 = mag1 + xksi(jj)*xksi(jj);    
end 
          
mag1 = 2.0*sqrt(mag1);
          
if mag1 > 0. && mag1 < elsiz
    elsiz = mag1;
end
          
mag1 = 0.0;
for jj = 1:2
    mag1 = mag1 + xeta(jj)*xeta(jj);    
end 
          
mag1 = 2.0*sqrt(mag1);
          
if mag1 > 0. && mag1 < elsiz
    elsiz = mag1;
end


SF = L0/elsiz*L0/elsiz;

me = me*SF;

end

end





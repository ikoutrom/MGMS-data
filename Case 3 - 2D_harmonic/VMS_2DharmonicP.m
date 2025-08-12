% input file for analysis with ACM and selective Mass Scaling
% 2D model without mass block on top, with harmonic loading
% (see Fig. 6 of manuscript) 
% and with CONFORMING Auxiliary coarse mesh
% the script runs the analysis and produces the base reaction and displacement
% history for the top of the model, stored in vectors Reac and Utip2, respectively. 
% the applied load values are stored in the vector Load
% the time values are stored in the vector tval
% the time history of Utip2 and of Reac is plotted in the end.

% print the CURRENT time (we do the same after the end of the script, so
% that we can compute total computation time)
tv1 = now;
d = datetime(tv1,'ConvertFrom','datenum') 

% *************** INPUT parameters ***************************************

xmax = 10; ymax = 80;  % dimensions of domain

necx = 2; % coarse discretization along the x- axis
necy = 16; % coarse discretization along the y- axis


% material properties
rho = 7.5e-7;
v = 0.25;
E = 29000;
% end material properties


Nsub = 20;      % subdivision factor for finer mesh

Pval = 10; % total horizontal force at the top
Tval = 0.2; % period of external harmonic force


dt = 2.e-5;  % timestep for analysis
tend = 0.4;

% **************************************** END INPUT PARAMETERS



%===========================================================  coarse mesh

nelc = necx*necy;
netip = necx*necx;


dxc = xmax/necx;
dyc = ymax/necy;

nncx = necx + 1;
nncy = necy + 1;

nc = nncx*nncy; % total number of nodes for the COARSE mesh

XY=zeros(2,nc);


 for iny = 1:nncy
   for inx = 1:nncx
        in1 = (iny-1)*nncx+inx;
        XY(1,in1) = (inx-1)*dxc;
        XY(2,in1) = (iny-1)*dyc;
   end
end



Uc = zeros(nc*2,1);
Vc = Uc;
Ac = Uc;

IDc = zeros(nc*2,1);

for i1 = 1:2*nncx   % all the bottom node DOFs are restrained
   IDc(i1) = -1;       % the value of -1 corresponds to zero prescribed displacement
end


LMc = zeros(4,necx*necy);
xl = zeros(4,2);
ul = zeros(8,1);


for iely = 1:necy
for ielx = 1:necx
iel = ielx + (iely-1)*necx;
LMc(1,iel) = ielx + (iely-1)*nncx;
LMc(2,iel) = ielx + 1 + (iely-1)*nncx;
LMc(3,iel) = ielx + 1 + iely*nncx;
LMc(4,iel) = ielx + iely*nncx;
end
end

% =============================  end setup of coarse mesh


%  **************************************  FINE mesh


Mscalef = Nsub*Nsub; % this is the factor by which we must scale the mass
%definition
nefx = necx*Nsub; 
nefy = necy*Nsub;
nelf = nefx*nefy;

dxf = xmax/nefx;
dyf = ymax/nefy;

nnfx = nefx + 1;
nnfy = nefy + 1;

nf = nnfx*nnfy; % total number of nodes for the FINE mesh

XYf=zeros(2,nf);



 for iny = 1:nnfy 
   for inx = 1:nnfx
        in1 = (iny-1)*nnfx+inx;
        XYf(1,in1) = (inx-1)*dxf;
        XYf(2,in1) = (iny-1)*dyf;
   end
end

Nv = zeros(4,1);

Tf = zeros(2*nf,4);  % transformation matrix values
Ifc = zeros(2*nf,4);  % mapping from fine to coarse DOFs


% now, populate Tf and Ifc!

   for in1 = 1:nf 
        
        for iel = 1:nelc
            x1 = XY(1,LMc(1,iel)); 
            y1 = XY(2,LMc(1,iel));
            x2 = XY(1,LMc(3,iel)); 
            y2 = XY(2,LMc(3,iel));
            
            
            if (XYf(1,in1) >= x1) && (XYf(1,in1) <= x2) && (XYf(2,in1) >= y1) && (XYf(2,in1) <= y2)
                Nv(1) = (XYf(1,in1)-x2)/(x1-x2)*(XYf(2,in1)-y2)/(y1-y2);
                Nv(2) = (XYf(1,in1)-x1)/(x2-x1)*(XYf(2,in1)-y2)/(y1-y2);
                Nv(3) = (XYf(1,in1)-x1)/(x2-x1)*(XYf(2,in1)-y1)/(y2-y1);
                Nv(4) = (XYf(1,in1)-x2)/(x1-x2)*(XYf(2,in1)-y1)/(y2-y1);
                
                for i1 = 1:4
                    Tf(2*in1-1,i1) = Nv(i1);
                    Tf(2*in1,i1) = Nv(i1);
                    Ifc(2*in1-1,i1) = 2*LMc(i1,iel)-1;
                    Ifc(2*in1,i1) = 2*LMc(i1,iel);
                end  
            end  
        end 
   end


Uf = zeros(nf*2,1);
Vf = Uf;
Af = Uf;

Ff = zeros(nf*2,1);

IDf = zeros(nf*2,1);

for i1 = 1:2*nnfx   % all the bottom node DOFs are restrained
   IDf(i1) = -1;       % the value of -1 corresponds to zero prescribed displacement
end

LMf = zeros(4,nefx*nefy);


for iely = 1:nefy
for ielx = 1:nefx
iel = ielx + (iely-1)*nefx;
LMf(1,iel) = ielx + (iely-1)*nnfx;
LMf(2,iel) = ielx + 1 + (iely-1)*nnfx;
LMf(3,iel) = ielx + 1 + iely*nnfx;
LMf(4,iel) = ielx + iely*nnfx;
end
end


% ============================================  assemble Mlf
Mf = zeros(nf*2,1);
    for iel = 1:nelf
        for in1 =1:4
            xl(in1,1) = XYf(1,LMf(in1,iel));
            xl(in1,2) = XYf(2,LMf(in1,iel));
            ul(2*in1-1) = Uf(2*LMf(in1,iel)-1);
            ul(2*in1) = Uf(2*LMf(in1,iel));
        end 
    %[fe,me] = call the Q4 function here!
    [fe,me] = elementQ4(xl,E,v,rho,ul);
        for in1 =1:4
           Mf(2*LMf(in1,iel)-1) = Mf(2*LMf(in1,iel)-1) + me(2*in1-1);
           Mf(2*LMf(in1,iel)) = Mf(2*LMf(in1,iel)) + me(2*in1);
        end 
    end
  
   %==============================================================

%  ***************************************** END FINE MESH




% ============================================  assemble Mlc from the
% contrained finer mesh Mlf:
Mc = zeros(2*nc,1);
for i1 = 1:2*nf
    for i2 = 1:4
        Mc(Ifc(i1,i2)) = Mc(Ifc(i1,i2)) + Mf(i1)*Tf(i1,i2);
    end
end
  
   %==============================================================
    
  Mf = Mf*Mscalef;  
    %Mc = [rho*A*Ltot/4; rho*A*Ltot/2; rho*A*Ltot/4];
%Mf = zeros(nnc,1);
%for i1 = 1:nnodf
%     Mc(Ifc(i1,1)) = Mc(Ifc(i1,1)) + Mf(i1)*Tf(i1,1);
%     Mc(Ifc(i1,2)) = Mc(Ifc(i1,2)) + Mf(i1)*Tf(i1,2);
%end


absfac = 1;
SFfac = 1;

dt = dt*absfac*SFfac;

ttim = 0.;
Nprint = 1000;

Nstep = int32((tend/dt)/Nprint);
tval = zeros(Nprint,1);
Utip = zeros(Nprint,1);
Utip2 = Utip;
Reac = zeros(Nprint,1);
Load = Reac;

for iprint = 1:Nprint
for istep=1:Nstep
    ttim = ttim + dt;
    
    
    Ff = zeros(nf*2,1);
    
    P0 = Pval*sin(2*pi/Tval*ttim);
    Pext = P0/nefx;
    
    % apply concentrated external force at the row of nodes at the top:
    for i1 = ((nnfy-1)*nnfx+2):(nnfy*nnfx-1)
        Ff(2*i1-1) = -Pext;
    end
    
    % we also put half the value of force for the first and last nodes of
    % the top row!
    i1 = (nnfy-1)*nnfx+1;
    Ff(2*i1-1) = -Pext/2;
    Ff(2*nf-1) = -Pext/2;
    
    for iel = 1:nelf
        for in1 =1:4
            xl(in1,1) = XYf(1,LMf(in1,iel));
            xl(in1,2) = XYf(2,LMf(in1,iel));
            ul(2*in1-1) = Uf(2*LMf(in1,iel)-1);
            ul(2*in1) = Uf(2*LMf(in1,iel));
            
            % add contribution of coarse-scale displacements:
            
            for in2 = 1:4
                ul(2*in1-1) = ul(2*in1-1) + Uc(Ifc(2*LMf(in1,iel)-1,in2))*Tf(2*LMf(in1,iel)-1,in2);
                ul(2*in1) = ul(2*in1) + Uc(Ifc(2*LMf(in1,iel),in2))*Tf(2*LMf(in1,iel),in2);
            end  
        end 
    %[fe,me] = call the Q4 function here!
    [fe,me] = elementQ4(xl,E,v,rho,ul);
        for in1 =1:4
           Ff(2*LMf(in1,iel)-1) = Ff(2*LMf(in1,iel)-1) + fe(2*in1-1);
           Ff(2*LMf(in1,iel)) = Ff(2*LMf(in1,iel)) + fe(2*in1);
        end 
    end
    
    % now, we assemble the global force vector for the coarse-scale mesh:
    
    Fc = zeros(nc*2,1);
    for i1 = 1:2*nf
        for i2 = 1:4
            Fc(Ifc(i1,i2)) = Fc(Ifc(i1,i2)) + Ff(i1)*Tf(i1,i2);
        end
    end
    
    % finally, we conduct update: first for the fine-scale dofs...
    
   for i1 = 1:2*nf  
   if IDf(i1) == 0
       Af(i1) = -Ff(i1)/Mf(i1);
       Vf(i1) = Vf(i1) + Af(i1)*dt;
       Uf(i1) = Uf(i1) + Vf(i1)*dt;
   end
   end
    
    % then, we update the coarse-scale dofs... NOTE: the essential BCs are
    % applied to the coarse-scale dofs!!
    
    
   for i1 = 1:2*nc  
   if IDc(i1) == 0
       Ac(i1) = -Fc(i1)/Mc(i1);
       Vc(i1) = Vc(i1) + Ac(i1)*dt;
       Uc(i1) = Uc(i1) + Vc(i1)*dt;
   elseif IDc(i1) == -2
      Vc(i1) = vtip;
      Uc(i1) = Uc(i1) + Vc(i1)*dt;
   end
   end
   
   Ptip1 = 0;
   for i1 = 1:nncx
       Ptip1 = Ptip1 + Fc(2*i1-1);
   end
   
end
   tval(iprint) = ttim;
   Utip(iprint) = Uc(2*nc-1);
   Utip2(iprint) = Uc(2*nc-1)+Uf(2*nf-1);
   Load(iprint) = P0;
   Reac(iprint) = Ptip1;
end

figure(3); plot(tval,Reac);
figure(4); plot(tval,Utip2);


% REprint the CURRENT time now that the script is complete
tv1 = now;
d = datetime(tv1,'ConvertFrom','datenum')
% input file for analysis with SINGLE MESH
% 2D model without mass block on top, with harmonic loading
% (see Fig. 6 of manuscript) 
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

necx = 20; % discretization along the x- axis
necy = 160; % discretization along the y- axis


% material properties
rho = 7.5e-5;
v = 0.25;
E = 29000;
% end material properties


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

% =============================  end setup of mesh



% ============================================  assemble Mlc
Mc = zeros(nc*2,1);
    for iel = 1:nelc
        for in1 =1:4
            xl(in1,1) = XY(1,LMc(in1,iel));
            xl(in1,2) = XY(2,LMc(in1,iel));
            ul(2*in1-1) = Uc(2*LMc(in1,iel)-1);
            ul(2*in1) = Uc(2*LMc(in1,iel));
        end 
    %[fe,me] = call the Q4 function here!
    [fe,me] = elementQ4(xl,E,v,rho,ul);
        for in1 =1:4
           Mc(2*LMc(in1,iel)-1) = Mc(2*LMc(in1,iel)-1) + me(2*in1-1);
           Mc(2*LMc(in1,iel)) = Mc(2*LMc(in1,iel)) + me(2*in1);
        end 
    end
  
   %==============================================================

%  ***************************************** END FINE MESH


ttim = 0.;
Nprint = 1000;

Nstep = int32((tend/dt)/Nprint);
tval = zeros(Nprint,1);
Utip = zeros(Nprint,1);
Reac = zeros(Nprint,1);
Load = zeros(Nprint,1);

for iprint = 1:Nprint
for istep=1:Nstep
    ttim = ttim + dt;
    
    
    Fc = zeros(nc*2,1);
    
    P0 = Pval*sin(2*pi/Tval*ttim);
    Pext = P0/necx;
    
    % apply concentrated external force at the row of nodes at the top:
    for i1 = ((nncy-1)*nncx+2):(nncy*nncx-1)
        Fc(2*i1-1) = -Pext;
    end
    
    % we also put half the value of force for the first and last nodes of
    % the top row!
    i1 = (nncy-1)*nncx+1;
    Fc(2*i1-1) = -Pext/2;
    Fc(2*nc-1) = -Pext/2;
    
    for iel = 1:nelc
        for in1 =1:4
            xl(in1,1) = XY(1,LMc(in1,iel));
            xl(in1,2) = XY(2,LMc(in1,iel));
            ul(2*in1-1) = Uc(2*LMc(in1,iel)-1);
            ul(2*in1) = Uc(2*LMc(in1,iel));
            
        end 
    %[fe,me] = call the Q4 function here!
    [fe,me] = elementQ4(xl,E,v,rho,ul);
        for in1 =1:4
           Fc(2*LMc(in1,iel)-1) = Fc(2*LMc(in1,iel)-1) + fe(2*in1-1);
           Fc(2*LMc(in1,iel)) = Fc(2*LMc(in1,iel)) + fe(2*in1);
        end 
    end
    
    % finally, we conduct update: first for the fine-scale dofs...
    Ptip1 = 0;
   for i1 = 1:2*nc  
   if IDc(i1) == 0
       Ac(i1) = -Fc(i1)/Mc(i1);
       Vc(i1) = Vc(i1) + Ac(i1)*dt;
       Uc(i1) = Uc(i1) + Vc(i1)*dt;
   end
   end
   
   for i1 = 1:nncx
       Ptip1 = Ptip1 + Fc(2*i1-1);
   end
   
end
   tval(iprint) = ttim;
   Utip(iprint) = Uc(2*nc-1);
   Reac(iprint) = Ptip1;
   Load(iprint) = P0;
end

figure(1); plot(tval,Utip);
figure(2); plot(tval,Load);
figure(3); plot(tval,Reac);

% REprint the CURRENT time now that the script is complete
tv1 = now;
d = datetime(tv1,'ConvertFrom','datenum')

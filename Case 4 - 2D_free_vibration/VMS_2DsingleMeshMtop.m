% input file for analysis with SINGLE MESH,
% 2D models with mass block on top (see Fig. 6 of manuscript) 
% the script runs the analysis and produces the load and displacement
% history for the top of the model, stored in vectors Ptip and Utip, respectively. 
% the time values are stored in the vector tval
% the time history of Utip is plotted in the end.


% IMPORTANT NOTE: for fine meshes with nonrectangular element shapes, 
% we may need to select the timestep by trial and error to ensure stability!


% *********************************** INPUT PARAMETERS
xmax = 10; ymax = 90;



% print the CURRENT time (we do the same after the end of the script, so
% that we can compute total computation time)
tv1 = now;
d = datetime(tv1,'ConvertFrom','datenum')


necx = 1; % number of elements in x direction
necy = 9; % number of elements in y direction


rho = 7.5e-7;  % mass density (change as desired to obtain uniform mas scaling)
v = 0.25;       % Poisson ratio
E = 29000;      % Elastic modulus

dt = 4.e-5; % analysis timestep
ev = 2; % imperfection (if we want non-rectangular elements!)
tend = 1.2;  % time for analysis termination






% *********************************** END INPUT PARAMETERS



%===========================================================  coarse mesh
%definition

nelc = necx*necy;

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

% add imperfection to mesh shape (use val1 = 0 for zero imperfection):
 val1 = -ev;
 for iny = 2:nncy-2
     inx = 1;
     in1 = (iny-1)*nncx+inx;
     XY(2,in1) = XY(2,in1)+val1;
     val1 = -val1;
 end 

 
 val1 = ev;
 for iny = 2:nncy-2
     inx = 2;
     in1 = (iny-1)*nncx+inx;
     XY(2,in1) = XY(2,in1)+val1;
     val1 = -val1;
 end 


Uc = zeros(nc*2,1);
Vc = Uc;
Ac = Uc;

IDc = zeros(nc*2,1);

for i1 = 1:2*nncx   % all the bottom node DOFs are restrained
   IDc(i1) = -1;       % the value of -1 corresponds to zero prescribed displacement
end

for i1 = nc-nncx+1:nc   % all the top node horizontal displ. are restrained
   IDc(2*i1-1) = -2;  % the value of -2 corresponds to prescribed displacement 
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

Nsub = 1;      % subdivision factor for finer mesh
Mscalef = Nsub*Nsub; % this is the factor by which we must scale the mass
%definition
nefx = necx*Nsub; 
nefy = necy*Nsub;
nelf = nefx*nefy;
netop = nefx*nefx;

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


nm1 = zeros(4,1); 
for iecy = 1:necy
  for iecx = 1:necx
      
      iec = (iecy-1)*necx + iecx;
      
    for i1 = 1:4
        nm1(i1) = LMc(i1,iec);
        xl(i1,1) = XY(1,LMc(i1,iec));
        xl(i1,2) = XY(2,LMc(i1,iec));
    end
    
    for in1y = 1:Nsub+1
    for in1x = 1:Nsub+1
      ifx =  (iecx-1)*Nsub + in1x; 
      ify =  (iecy-1)*Nsub + in1y;
        inodf = (ify-1)*nnfx + ifx;
        
        ksi = -1 + (in1x-1)*2/Nsub;
        eta = -1 + (in1y-1)*2/Nsub;
        
        Nv(1) = 1/4*(1-ksi)*(1-eta);
        Nv(2) = 1/4*(1+ksi)*(1-eta);
        Nv(3) = 1/4*(1+ksi)*(1+eta);
        Nv(4) = 1/4*(1-ksi)*(1+eta);
        
        XYf(1,inodf)=0;
        XYf(2,inodf)=0;
        for i2 = 1:4
            Ifc(2*inodf-1,i2) = 2*nm1(i2)-1;
            Ifc(2*inodf,i2) = 2*nm1(i2);
            Tf(2*inodf-1,i2) = Nv(i2);
            Tf(2*inodf,i2) = Nv(i2);
            XYf(1,inodf)= XYf(1,inodf) + Nv(i2)*xl(i2,1);
            XYf(2,inodf)= XYf(2,inodf) + Nv(i2)*xl(i2,2);
        end
        
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

for i1 = nf-nnfx+1:nf   % all the top node horizontal displ. are restrained
   IDf(2*i1-1) = -1;  % the value of -1 corresponds to zero prescribed FINE-SCALE displacement
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
    if iel > (nelf-netop)
        fe = fe*200; % the top-block elements have a thickness of 200, so that we get a lumped mass at top!
        me = me*200;
    end
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
Nprint = 500;

Nstep = int32((tend/dt)/Nprint);
tval = zeros(Nprint*1.2,1);
Utip = zeros(Nprint*1.2,1);
Ptip = zeros(Nprint*1.2,1);
vtip = 0.;

for iprint = 1:Nprint
for istep=1:Nstep
    ttim = ttim + dt;
    
    Fc = zeros(nc*2,1);
    
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
    
    if ttim<=1
        vtip = ttim;
    else
        vtip = 1.;
    end
    
    Ptip1 = 0.;
   for i1 = 1:2*nc  
   if IDc(i1) == 0
       Ac(i1) = -Fc(i1)/Mc(i1);
       Vc(i1) = Vc(i1) + Ac(i1)*dt;
       Uc(i1) = Uc(i1) + Vc(i1)*dt;
   elseif IDc(i1) == -2
      Vc(i1) = vtip;
      Uc(i1) = Uc(i1) + Vc(i1)*dt;
      Ptip1 = Ptip1 + Fc(i1);
   end
   end
   
   
end
   tval(iprint) = ttim;
   Utip(iprint) = Uc(2*nc-1);
   Ptip(iprint) = Ptip1;
end

% now, we remove the restraint at the top and allow free vibration

for i1 = nc-nncx+1:nc   % all the top node horizontal displ. are restrained
   IDc(2*i1-1) = 0;  % we free the originally restrained nodes at the top 
end


for iprint = Nprint+1:1.2*Nprint
for istep=1:Nstep
    ttim = ttim + dt;
    
    Fc = zeros(nc*2,1);
    
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
    
    if ttim<=1
        vtip = ttim;
    else
        vtip = 1.;
    end
    
    Ptip1 = 0.;
   for i1 = 1:2*nc  
   if IDc(i1) == 0
       Ac(i1) = -Fc(i1)/Mc(i1);
       Vc(i1) = Vc(i1) + Ac(i1)*dt;
       Uc(i1) = Uc(i1) + Vc(i1)*dt;
   elseif IDc(i1) == -2
      Vc(i1) = vtip;
      Uc(i1) = Uc(i1) + Vc(i1)*dt;
      Ptip1 = Ptip1 + Fc(i1);
   end
   end
   
   
end
   tval(iprint) = ttim;
   Utip(iprint) = Uc(2*nc-1);
   Ptip(iprint) = Ptip1;
end
    
% REprint the CURRENT time now that the script is complete
tv1 = now;
d = datetime(tv1,'ConvertFrom','datenum')


% generate plots

%% Figure properties
Fig1=figure('Name','Tip displacement','Position',[50 50 700 950],'color','w');
xlabel('time (sec)','FontSize',16,'FontAngle','normal','FontWeight','bold','FontName','Times New Roman');
ylabel('displacement (in.)','FontSize',16,'FontAngle','normal','FontWeight','bold','FontName','Times New Roman');
xtickformat('%,.2f')
ytickformat('%,.1f')
axis([0 1.1*ttim 1.1*min(Utip) 1.1*max(Utip)])
    axis on
    hold all
    box on

set(gca,'fontsize', 16)  

set(0,'CurrentFigure',Fig1)

%% Figure Generation
plot(tval,Utip,'LineStyle','-','LineWidth',2,'Color',[0 0 0]);

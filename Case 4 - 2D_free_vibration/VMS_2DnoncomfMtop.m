% input file for analysis with ACM and selective Mass Scaling
% 2D models with mass block on top (see Fig. 6 of manuscript) 
% with NON-CONFORMING Auxiliary coarse mesh
% the script runs the analysis and produces the load and displacement
% history for the top of the model, stored in vectors Ptip and Utip2, respectively. 
% the time values are stored in the vector tval 
% the time history of Utip2 is plotted in the end.


% Note: the mesh is created as follows:
% Step 1: create a temporary, "conforming coarse mesh"
% Step 2: subdivide the "conforming coarse mesh" so that we obtain the
% fine, ACTUAL mesh
% Step 3: Revise the coarse mesh, so that it obtains the final,
% non-conforming shape of the ACM in the manuscript.

% IMPORTANT NOTE: in the specific analysis, we apply mass scaling 
% INDIVIDUALLY for the elements! The reason is that the elements may have
% irregular shapes, so some elements may have smaller size than simply
% Lc/Nsub, where Lc is the "coarse mesh size"!
% this is why the element computations for the unscaled mass use a Q4 routine "elementQ4",
% while the scaled mass matrix involves computations calling another Q4 routine "elementQ4scaleMass"

xmax = 10; ymax = 90;   % domain dimensions in the x- and y- direction

% print the CURRENT time (we do the same after the end of the script, so
% that we can compute total computation time)
tv1 = now;
d = datetime(tv1,'ConvertFrom','datenum') 



%************************************** input parameters


Nsub = 20;      % subdivision factor for finer mesh

ev = 2; % imperfection value for mesh (use zero for perfect mesh)

dt = 4.e-5; % analysis step

tend = 1;   % termination time for prescribed-displacement stage


Nprint = 500;   % number of steps to be printed for prescribed-displacement stage

%************************************** end input parameters


%===========================================================  coarse mesh
%definition
necx = 1; %2 % number of ACM elements in x-direction
necy = 9; %16   % number of ACM elements in y-direction
nelc = necx*necy; % total number of ACM elements
netip = necx*necx; % number of ACM elements in the mass-block region at the top 


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
 
 % end of adding imperfection
 

rho = 7.5e-7;
v = 0.25;
E = 29000;

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
            XYf(1,inodf)= XYf(1,inodf) + Nv(i2)*xl(i2,1);
            XYf(2,inodf)= XYf(2,inodf) + Nv(i2)*xl(i2,2);
        end
        
    end
    end
    
end   
end


% revise coordinates of coarse mesh (so that coarse mesh is rectangular!)

 for iny = 1:nncy
   for inx = 1:nncx
        in1 = (iny-1)*nncx+inx;
        XY(1,in1) = (inx-1)*dxc;
        XY(2,in1) = (iny-1)*dyc;
   end
 end

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
    
  %Mf = Mf*Mscalef;
  
  
  % ============================================  reassemble Mlf, applying
  % mass scaling on the ELEMENT LEVEL
Mf = zeros(nf*2,1);
    for iel = 1:nelf
        for in1 =1:4
            xl(in1,1) = XYf(1,LMf(in1,iel));
            xl(in1,2) = XYf(2,LMf(in1,iel));
            ul(2*in1-1) = Uf(2*LMf(in1,iel)-1);
            ul(2*in1) = Uf(2*LMf(in1,iel));
        end 
    %[fe,me] = call the Q4 function here!
    [fe,me] = elementQ4scaleMass(xl,E,v,rho,dxc,ul);
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
  
  
  
  
  
    %Mc = [rho*A*Ltot/4; rho*A*Ltot/2; rho*A*Ltot/4];
%Mf = zeros(nnc,1);
%for i1 = 1:nnodf
%     Mc(Ifc(i1,1)) = Mc(Ifc(i1,1)) + Mf(i1)*Tf(i1,1);
%     Mc(Ifc(i1,2)) = Mc(Ifc(i1,2)) + Mf(i1)*Tf(i1,2);
%end

ttim = 0.;


Nstep = int32((tend/dt)/Nprint);
tval = zeros(Nprint*1.2,1);
Utip = zeros(Nprint*1.2,1);
Utip2 = Utip;
Ptip = zeros(Nprint*1.2,1);
vtip = 0.;

for iprint = 1:Nprint
for istep=1:Nstep
    ttim = ttim + dt;
    
    
    Ff = zeros(nf*2,1);
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
    if iel > (nelf-netop)
        fe = fe*200; % the top-block elements have a thickness of 200!
    end
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
   Utip2(iprint) = Uc(2*nc-1)+Uf(2*nf-1);
   Ptip(iprint) = Ptip1;
end

% now, we remove the restraint at the top and allow free vibration

for i1 = nc-nncx+1:nc   % all the top node horizontal displ. are restrained
   IDc(2*i1-1) = 0;  % we free the originally restrained nodes at the top 
end

for i1 = nf-nnfx+1:nf   % all the top node horizontal displ. are restrained
   IDf(2*i1-1) = 0;  % we free the originally restrained nodes at the top 
end

    
for iprint = Nprint+1:1.2*Nprint
for istep=1:Nstep
    ttim = ttim + dt;
    Ff = zeros(nf*2,1);
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
   Utip2(iprint) = Uc(2*nc-1)+Uf(2*nf-1);
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
axis([0 1.1*ttim 1.1*min(Utip2) 1.1*max(Utip2)])
    axis on
    hold all
    box on

set(gca,'fontsize', 16)  

set(0,'CurrentFigure',Fig1)

%% Figure Generation
plot(tval,Utip2,'LineStyle','-','LineWidth',2,'Color',[0 0 0]);

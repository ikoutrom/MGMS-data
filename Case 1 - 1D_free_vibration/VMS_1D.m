nec = nec0;
dxc = xmax/nec;

nc = nec + 1;

XY=zeros(nc,1);


 for in1 = 1:nc
    XY(in1) = (in1-1)*dxc;
 end




Uc = zeros(nc,1);
Vc = Uc;
Ac = Uc;

IDc = zeros(nc,1);

IDc(1) = -1;

% uncomment the following line for cases involving prescribed displacement!
%IDc(nc) = -2; % -2 means prescribed motion

LMc = zeros(2,nec);
xl = zeros(2,1);
ul = zeros(2,1);


for iel = 1:nec
LMc(1,iel) = iel;
LMc(2,iel) = iel + 1;
end

% initial displacement
for i1 = 1:nc
   Uc(i1) = (1-idispsh)*ampl*xmax*sin(pi/2*XY(i1)/xmax) + idispsh*ampl*XY(i1);
end

% =============================  end setup of coarse mesh


%  **************************************  FINE mesh

Mscalef = Nsub*Nsub; % this is the factor by which we must scale the mass
%definition
nef = nec*Nsub; 

hsg = zeros(nef,nhv);
hl = zeros(nhv,1);

dxf = xmax/nef;

nf = nef + 1;

XYf=zeros(nf,1);

   for in1 = 1:nf
        XYf(in1) = (in1-1)*dxf;
   end

Nv = zeros(2,1);

Tf = zeros(nf,2);  % transformation matrix values
Ifc = zeros(nf,2);  % mapping from fine to coarse DOFs


% now, populate Tf and Ifc!

   for in1 = 1:nf 
        
        for iel = 1:nec
            x1 = XY(LMc(1,iel)); 
            x2 = XY(LMc(2,iel)); 
            
            
            if (XYf(in1) >= x1) && (XYf(in1) <= x2)
                Nv(1) = (XYf(in1)-x2)/(x1-x2);
                Nv(2) = (XYf(in1)-x1)/(x2-x1);
                
                for i1 = 1:2
                    Tf(in1,i1) = Nv(i1);
                    Ifc(in1,i1) = LMc(i1,iel);
                end  
            end  
        end 
   end


Uf = zeros(nf,1);
Vf = Uf;
Af = Uf;

Ff = zeros(nf,1);

IDf = zeros(nf,1);

IDf(1) = -1;

% uncomment the following line, for the case that we have prescribed
% motion:
%IDf(nf) = -1;

LMf = zeros(2,nef);


for iel = 1:nef
    LMf(1,iel) = iel;
    LMf(2,iel) = iel + 1;
end

Area = zeros(nef,1);
for i1 = 1:nef
   Area(i1) = A; 
end


for i1 = 1:Nsub
   Area((nec/2-1)*Nsub+i1) = area2; 
   Area((nec-1)*Nsub+i1) = area2;
end

% ============================================  assemble Mlf
Mf = zeros(nf,1);
    for iel = 1:nef
        for in1 =1:2
            xl(in1) = XYf(LMf(in1,iel));
            xl(in1) = XYf(LMf(in1,iel));
        end 
        [fe,me] = truss1D(xl,ul,Area(iel),E,rho,fy,Hmod,hsg(iel,:));
        for in1 =1:2
           Mf(LMf(in1,iel)) = Mf(LMf(in1,iel)) + me(in1);
        end 
    end
   %==============================================================

%  ***************************************** END FINE MESH




% ============================================  assemble Mlc from the
% contrained finer mesh Mlf:
Mc = zeros(nc,1);
for i1 = 1:nf
    for i2 = 1:2
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
Nprint = int32((tend/dt)/Nstep);
tval = zeros(Nprint,1);
Utip = zeros(Nprint,1);
Utip2 = Utip;
Ptip = zeros(Nprint,1);
vtip = 0.;


% conduct a first half-step

    Ff = zeros(nf,1);
    for iel = 1:nef
        for in1 =1:2
            xl(in1) = XYf(LMf(in1,iel));
            ul(in1) = Uf(LMf(in1,iel));
            
            % add contribution of coarse-scale displacements:
            
            for in2 = 1:2
                ul(in1) = ul(in1) + Uc(Ifc(LMf(in1,iel),in2))*Tf(LMf(in1,iel),in2);
            end  
        end 
        [fe,me] = truss1D(xl,ul,Area(iel),E,rho,fy,Hmod,hsg(iel,:));
        for in1 =1:2
           Ff(LMf(in1,iel)) = Ff(LMf(in1,iel)) + fe(in1);
        end 
    end
    
    % now, we assemble the global force vector for the coarse-scale mesh:
    
    Fc = zeros(nc,1);
    for i1 = 1:nf
        for i2 = 1:2
            Fc(Ifc(i1,i2)) = Fc(Ifc(i1,i2)) + Ff(i1)*Tf(i1,i2);
        end
    end
    
    % finally, we conduct update: first for the fine-scale dofs...
    
   for i1 = 1:nf  
   if IDf(i1) == 0
       Af(i1) = -Ff(i1)/Mf(i1);
       Vf(i1) = Vf(i1) + 0.5*Af(i1)*dt;
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
    
   for i1 = 1:nc  
   if IDc(i1) == 0
       Ac(i1) = -Fc(i1)/Mc(i1);
       Vc(i1) = Vc(i1) + 0.5*Ac(i1)*dt;
       Uc(i1) = Uc(i1) + Vc(i1)*dt;
   elseif IDc(i1) == -2
      Vc(i1) = vtip;
      Uc(i1) = Uc(i1) + Vc(i1)*dt;
   end
   end
   
   



for iprint = 1:Nprint
for istep=1:Nstep
    ttim = ttim + dt;
    Ff = zeros(nf,1);
    for iel = 1:nef
        for in1 =1:2
            xl(in1) = XYf(LMf(in1,iel));
            ul(in1) = Uf(LMf(in1,iel));
            
            % add contribution of coarse-scale displacements:
            
            for in2 = 1:2
                ul(in1) = ul(in1) + Uc(Ifc(LMf(in1,iel),in2))*Tf(LMf(in1,iel),in2);
            end  
        end 
        [fe,me] = truss1D(xl,ul,Area(iel),E,rho,fy,Hmod,hsg(iel,:));
        for in1 =1:2
           Ff(LMf(in1,iel)) = Ff(LMf(in1,iel)) + fe(in1);
        end 
    end
    
    % now, we assemble the global force vector for the coarse-scale mesh:
    
    Fc = zeros(nc,1);
    for i1 = 1:nf
        for i2 = 1:2
            Fc(Ifc(i1,i2)) = Fc(Ifc(i1,i2)) + Ff(i1)*Tf(i1,i2);
        end
    end
    
    % finally, we conduct update: first for the fine-scale dofs...
    
   for i1 = 1:nf  
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
    
   for i1 = 1:nc  
   if IDc(i1) == 0
       Ac(i1) = -Fc(i1)/Mc(i1);
       Vc(i1) = Vc(i1) + Ac(i1)*dt;
       Uc(i1) = Uc(i1) + Vc(i1)*dt;
   elseif IDc(i1) == -2
      Vc(i1) = vtip;
      Uc(i1) = Uc(i1) + Vc(i1)*dt;
   end
   end
   
   
end
   tval(iprint) = ttim;
   Utip(iprint) = Uc(nc);
   Utip2(iprint) = Uc(nc)+Uf(nf);
   Ptip(iprint) = Fc(nc);
end
  

nec = nec0*Nsub; %total number of coarse-scale elements
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


Area = zeros(nec,1);
for i1 = 1:nec
   Area(i1) = A; 
end


for i1 = 1:Nsub
   Area((nec0/2-1)*Nsub+i1) = area2; 
   Area((nec0-1)*Nsub+i1) = area2;
end

hsg = zeros(nec,3);

% =============================  end setup of coarse mesh

for i1 = 1:nc
   Uc(i1) = (1-idispsh)*ampl*xmax*sin(pi/2*XY(i1)/xmax) + idispsh*ampl*XY(i1);
end


% ============================================  assemble Mlf
Mc = zeros(nc,1);
    for iel = 1:nec
        for in1 =1:2
            xl(in1) = XY(LMc(in1,iel));
            xl(in1) = XY(LMc(in1,iel));
        end 
        [fe,me] = truss1D(xl,ul,Area(iel),E,rho,fy,Hmod,hsg(iel,:));
        for in1 =1:2
           Mc(LMc(in1,iel)) = Mc(LMc(in1,iel)) + me(in1);
        end 
    end
  
   %==============================================================



absfac = 1;
SFfac = 1;

dt = dt*absfac*SFfac;



ttim = 0.;
Nprint = int32((tend/dt))/Nstep; % we set up in a way that we print all steps!
tval = zeros(Nprint,1);
Utip = zeros(Nprint,1);
Vtip = Utip;
Utip2 = Utip;
Ptip = zeros(Nprint,1);
vtip = 0.;

Fc = zeros(nc,1);
    for iel = 1:nec
        for in1 =1:2
            xl(in1) = XY(LMc(in1,iel));
            ul(in1) = Uc(LMc(in1,iel)); 
        end 
        [fe,me] = truss1D(xl,ul,Area(iel),E,rho,fy,Hmod,hsg(iel,:));
        for in1 =1:2
           Fc(LMc(in1,iel)) = Fc(LMc(in1,iel)) + fe(in1);
        end 
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
    Fc = zeros(nc,1);
    for iel = 1:nec
        for in1 =1:2
            xl(in1) = XY(LMc(in1,iel));
            ul(in1) = Uc(LMc(in1,iel)); 
        end 
        [fe,me] = truss1D(xl,ul,Area(iel),E,rho,fy,Hmod,hsg(iel,:));
        for in1 =1:2
           Fc(LMc(in1,iel)) = Fc(LMc(in1,iel)) + fe(in1);
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
   Vtip(iprint) = Vc(nc);
   Ptip(iprint) = Fc(nc);
end
 

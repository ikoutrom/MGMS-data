% model for 1D mesh, monotonically increasing displacement at tip
% viscoplastic material behavior with softening
% (see Figs. 4, 10 and 11 of manuscript) 
% the script runs the analysis and produces the plots of the quantities
% presented in Figures 10 and 11 of the manuscript. 
% important output vectors: Utip = end displacement history, Ptip = end
% force history, tval = time value labels


% print the CURRENT time (we do the same after the end of the script, so
% that we can compute total computation time)
tv1 = now;
d = datetime(tv1,'ConvertFrom','datenum') 

% INPUT DEFINITIONS
xmax = 50;  % total length
nec = 10; %2 total number of coarse-scale elements
Nsub = 20;      % subdivision factor for finer mesh
%============================================= material properties

rho = 7.5e-7;
E = 29000;
A = 1;
eta = 0.04*E;  % viscoplasticity parameter
fy = 50; % initial yield stress
Hmod = -50; % hardening modulus (negative means SOFTENING)
nhv = 3;    % number of history variables per element

%============================================= end material properties


dt = 2.e-5;  % solution timestep
tend = 0.04;  % end of analysis time


%*****************************************  END OF INPUT parameters






%===========================================================  coarse mesh


dxc = xmax/nec;

nc = nec + 1;

XY=zeros(nc,1);


 for in1 = 1:nc
    XY(in1) = (in1-1)*dxc;
 end


 %dxc = 0;
 


Uc = zeros(nc,1);
Vc = Uc;
Ac = Uc;

IDc = zeros(nc,1);

IDc(1) = -1;

% uncomment the following line for cases involving prescribed displacement!
IDc(nc) = -2; % -2 means prescribed motion

LMc = zeros(2,nec);
xl = zeros(2,1);
ul = zeros(2,1);


for iel = 1:nec
LMc(1,iel) = iel;
LMc(2,iel) = iel + 1;
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
IDf(nf) = -1; % the final node of the fine mesh corresponds to a location with prescribed displacement, so we set its value to 0!

LMf = zeros(2,nef);


for iel = 1:nef
    LMf(1,iel) = iel;
    LMf(2,iel) = iel + 1;
end


% ============================================  assemble Mlf
Mf = zeros(nf,1);
    for iel = 1:nef
        for in1 =1:2
            xl(in1) = XYf(LMf(in1,iel));
            xl(in1) = XYf(LMf(in1,iel));
        end 
        [fe,me,hsg(iel,:)] = truss1Dvp(xl,ul,A,E,rho,fy,Hmod,eta,hsg(iel,:));
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
Nprint = int32((tend/dt)); % we set up in a way that we print all steps!
Nstep = 1;
tval = zeros(Nprint,1);
Utip = zeros(Nprint,1);
Utip2 = Utip;
Ptip = zeros(Nprint,1);
vtip = 0.;




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
        [fe,me,hsg(iel,:)] = truss1Dvp(xl,ul,A,E,rho,fy,Hmod,eta,hsg(iel,:));  

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
    
    if ttim<=0.05
        vtip = ttim/0.05;
    else
        vtip = 1.;
    end
    
    vtip = vtip*xmax;
    
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
  

% generate plots

%% Figure properties
Fig1=figure('Name','Tip displacement','Position',[50 50 700 950],'color','w');
xlabel('time sec','FontSize',16,'FontAngle','normal','FontWeight','bold','FontName','Times New Roman');
ylabel('displacement (in.)','FontSize',16,'FontAngle','normal','FontWeight','bold','FontName','Times New Roman');
xtickformat('%,.2f')
ytickformat('%,.1f')
axis([0 1.1*ttim 0 1.1*max(Utip)])
    axis on
    hold all
    box on

set(gca,'fontsize', 16)  

set(0,'CurrentFigure',Fig1)
%==========================================================================
%% Figure Generation
plot(tval,Utip,'LineStyle','-','LineWidth',2,'Color',[0 0 0]);

Fig2=figure('Name','Force-Displacement Response','Position',[700 50 700 950],'color','w');
xlabel('displacement (in.)','FontSize',16,'FontAngle','normal','FontWeight','bold','FontName','Times New Roman');
ylabel('force (kip)','FontSize',16,'FontAngle','normal','FontWeight','bold','FontName','Times New Roman');
xtickformat('%,.1f')
ytickformat('%,.0f')
axis([0 1.1*max(Utip) 0 1.1*max(Ptip)])
    axis on
    hold all
    box on

set(gca,'fontsize', 16)  
set(0,'CurrentFigure',Fig2)
plot(Utip,Ptip,'LineStyle','-','LineWidth',2,'Color',[0 0 0]);
%figure(1); plot(tval,Utip)
%figure(2); plot(Utip,Ptip)

epeff = zeros(nef,1);
xax1 = zeros(nef,1);
for i1=1:nef
   epeff(i1) = hsg(i1,3); 
end

xax1(1) = 0.5/nef;
for i1=2:nef
   xax1(i1) = xax1(i1-1)+1/nef; 
end


Fig3=figure('Name','Effective Plastic Strain','Position',[1350 50 700 950],'color','w');
xlabel('element number','FontSize',16,'FontAngle','normal','FontWeight','bold','FontName','Times New Roman');
ylabel('effective plastic strain','FontSize',16,'FontAngle','normal','FontWeight','bold','FontName','Times New Roman');
%xtickformat('%,.1f')
ytickformat('%,.2f')
axis([0 200 0 1.1*max(epeff)])
    axis on
    hold all
    box on

set(gca,'fontsize', 16)  
set(0,'CurrentFigure',Fig3)
plot(epeff,'LineStyle','-','LineWidth',2,'Color',[0 0 0]);

% REprint the CURRENT time now that the script is complete
tv1 = now;
d = datetime(tv1,'ConvertFrom','datenum')

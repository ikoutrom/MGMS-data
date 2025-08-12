% model for SINGLE 1D mesh (no ACM), monotonically increasing displacement at tip
% elastoplastic material behavior
% (see Figs. 4, 10 and 11 of manuscript) 
% the script runs the analysis and produces the plots of the quantities
% presented in Figures 10 and 11 of the manuscript. 
% important output vectors: Utip = end displacement history, Ptip = end
% force history, tval = time value labels

% print the CURRENT time (we do the same after the end of the script, so
% that we can compute total computation time)
tv1 = now;
d = datetime(tv1,'ConvertFrom','datenum') 

% *********************************** INPUT DEFINITIONS
xmax = 50;  % total length
nec = 10; % total number of coarse-scale elements
Nsub = 20;      % subdivision factor for finer mesh
%============================================= material properties

rho = 7.5e-7;
E = 29000;
A = 1;
fy = 50; % initial yield stress
Hmod = 0.015*E; % hardening modulus (negative means SOFTENING)
nhv = 3;    % number of history variables per element

%============================================= end material properties


dt = 2.e-5;  % solution timestep
tend = 0.04;  % end of analysis time

Nstep = 1; % how many substeps (between prints) we want to have - if =1, then we print ALL steps!

%******************************************************* END INPUT PARAM.

%===========================================================  coarse mesh

nec = nec*Nsub; 

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
IDc(nc) = -2; % -2 means prescribed motion

LMc = zeros(2,nec);
xl = zeros(2,1);
ul = zeros(2,1);


for iel = 1:nec
LMc(1,iel) = iel;
LMc(2,iel) = iel + 1;
end


hsg = zeros(nec,3);

% =============================  end setup of coarse mesh


% ============================================  assemble Mlf
Mc = zeros(nc,1);
    for iel = 1:nec
        for in1 =1:2
            xl(in1) = XY(LMc(in1,iel));
            xl(in1) = XY(LMc(in1,iel));
        end 
        [fe,me,hsg(iel,:)] = truss1D(xl,ul,A,E,rho,fy,Hmod,hsg(iel,:));
        for in1 =1:2
           Mc(LMc(in1,iel)) = Mc(LMc(in1,iel)) + me(in1);
        end 
    end
  
   %==============================================================


absfac = 1;
SFfac = 1;

dt = dt*absfac*SFfac;


ttim = 0.;
Nprint = int32((tend/dt)/Nstep); 
tval = zeros(Nprint,1);
Utip = zeros(Nprint,1);
Vtip = Utip;
Ptip = zeros(Nprint,1);
vtip = 0.;


for iprint = 1:Nprint
for istep=1:Nstep
    ttim = ttim + dt;
    Fc = zeros(nc,1);
    for iel = 1:nec
        for in1 =1:2
            xl(in1) = XY(LMc(in1,iel));
            ul(in1) = Uc(LMc(in1,iel)); 
        end 
        [fe,me,hsg(iel,:)] = truss1D(xl,ul,A,E,rho,fy,Hmod,hsg(iel,:));
        for in1 =1:2
           Fc(LMc(in1,iel)) = Fc(LMc(in1,iel)) + fe(in1);
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
   Vtip(iprint) = Vc(nc);
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

epeff = zeros(nec,1);
for i1=1:nec
   epeff(i1) = hsg(i1,3); 
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
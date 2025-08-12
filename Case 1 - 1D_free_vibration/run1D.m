% model for free vibration solution of 1D mesh - 
% analyses in Figs. 4 and 5 of the manuscript
% see Fig. 4 for geometry and model setup
% linear or sinunsoidal initial displacement pattern

clear; 
clc;

%***************************************** INPUT DEFINITIONS

xmax = 50; % total length of domain

nec0 = 10; % total number of coarse-scale elements
Nsub = 20;      % subdivision factor for finer mesh
% - NOTE: the total number of elements in the ACTUAL mesh will be
% nec0*Nsub!


%------------- material/sectional properties

rho = 7.5e-7; % density
E = 29000; % elastic modulus
A = 1; % cross-sectional area
fy = 1.e32; % initial yield stress (use huge value if we want linear elasticity)
Hmod = 290; % hardening modulus (zero = perfect plasticity, negative = softening)
nhv = 3;    % number of history variables Gauss point
area2 = 1; % % increased cross-sectional area for tip and middle regions (put "1" for having the same value as others)
% end material properties

ampl = 0.001; % initial displacement amplitude parameter
idispsh = 1; % initial displacement shape (give "0" for sinusoidal, "1" for linear!)

tend = 0.01;  % analysis termination time

Nstep = 1; % number of cycles run between prints: =1 means we print ALL steps!

%********************************************* END INPUT DEFINITIONS



% first, run single mesh without mass scaling:
dt = 1.e-6;  % solution timestep
VMS_1DsingleMesh;
tvala = tval; Ua = Utip/(ampl*xmax);


% next, run single mesh with mass scaling:
rho = 20^2*rho;
dt = 2.e-5;  % solution timestep (we increase it by a factor of 20!)

VMS_1DsingleMesh;
tvalb = tval; Ub = Utip/(ampl*xmax);

% finally, run with Selective MS scheme in manuscript, for the increased
% solution timestep:

rho = 7.5e-7; % density
VMS_1D;
tvalc = tval; Uc = Utip2/(ampl*xmax);

% generate plots

%% Figure properties
Fig1=figure('Name','Tip displacement','Position',[150 50 500 400],'color','w');
xlabel('time (sec)','FontSize',16,'FontAngle','normal','FontWeight','bold','FontName','Times New Roman');
ylabel('Normalized Tip Displacement','FontSize',16,'FontAngle','normal','FontWeight','bold','FontName','Times New Roman');
xtickformat('%,.3f')
ytickformat('%,.1f')
axis([0 0.011 -1.1 1.1])
    axis on
    hold all
    box on
%legend('Coarse Mesh','Fine Mesh, no MS','Fine Mesh, with MS','MBMS, fine mesh','Location','northwest');

%legend('Coarse Mesh','Fine Mesh, no MS');
%legend('Location','northwest');

%set(ch(3),'FontSize',12,'FontAngle','normal','FontWeight','bold','FontName','Times New Roman');

set(gca,'fontsize', 16)  

set(0,'CurrentFigure',Fig1)
%==========================================================================
%% Figure Generation
plot(tvala,Ua,'LineStyle','-','LineWidth',4,'Color',[0.75 0.75 0.75]);
plot(tvalb,Ub,'LineStyle',':','LineWidth',2,'Color',[0 0 0]);
plot(tvalc,Uc,'LineStyle','-','LineWidth',1,'Color',[0 0 1]);

lgd = legend({'no MS','with MS','MBMS'},...
    'FontSize',16, 'FontAngle','normal','FontWeight','normal','FontName','Times New Roman','Location','best');

%h=legend('Experiment','Analysis','Location','northwest');
%ch=get(h,'Children');


clear; close all; clc

%% ---- Propeller Geometry ----
% Select propeller (update if needed for your case)
PropellerName = 'APC20x8';

% Define radial positions from 0.1R to tip (10 elements as required)
r1 = linspace(0.1, 1, 10);

% Get chord and pitch distributions from provided function
[chord, pitch] = PropellerDimension(r1, PropellerName);

%% ---- Colour setup (just for nicer plots) ----
myred       = [216 30 49]/255;
myblue      = [27 99 157]/255;
myblack     = [0 0 0]/255;
mygreen     = [0 128 0]/255;
mycyan      = [2 169 226]/255;
myyellow    = [251 194 13]/255;
mygray      = [89 89 89]/255;

% Set default plotting colour order
set(groot,'defaultAxesColorOrder',[myblack;myblue;myred;mygreen;myyellow;mycyan;mygray]);

%% ---- General plot formatting ----
alw = 1;        % axis line width
fsz = 22;       % font size (kept large for report figures)
lw  = 2.5;      % line thickness
msz = 8;        % marker size

set(0,'defaultLineLineWidth',lw);
set(0,'defaultLineMarkerSize',msz);

%% ---- Figure 1: Pitch Distribution ----
hFig = figure(1); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

% Convert pitch from radians to degrees for plotting
plot(r1, pitch * 180/pi)

xlabel('Normalised radius r/R')
ylabel('Pitch angle (deg)')
xlim([0.1 1])
box on

% Make the plot look clean and consistent
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

% Export high-quality image for report
print -dpng -r300 Pitch_Distribution

%% ---- Figure 2: Chord Distribution ----
hFig = figure(2); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(r1, chord)

xlabel('Normalised radius r/R')
ylabel('Chord (m)')
xlim([0.1 1])
box on

set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

print -dpng -r300 Chord_Distribution

%% ---- Read XFOIL Data ----
% Load airfoil data for different Reynolds numbers
Re500  = readmatrix('Naca4412Re500');
Re650  = readmatrix('Naca4412Re650');
Re800  = readmatrix('Naca4412Re800');
Re1000 = readmatrix('Naca4412Re1000');

% Extract columns (alpha, CL, CD)
alpha500  = Re500(:,1);  cl500  = Re500(:,2);  cd500  = Re500(:,3);
alpha650  = Re650(:,1);  cl650  = Re650(:,2);  cd650  = Re650(:,3);
alpha800  = Re800(:,1);  cl800  = Re800(:,2);  cd800  = Re800(:,3);
alpha1000 = Re1000(:,1); cl1000 = Re1000(:,2); cd1000 = Re1000(:,3);

%% ---- Figure 3: CL vs Alpha ----
hFig = figure(3); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

% Plot lift curves (higher Re should generally give better performance)
plot(alpha650,  cl650)
plot(alpha800,  cl800)
plot(alpha1000, cl1000)

xlabel('Angle of attack \alpha (deg)')
ylabel('Lift coefficient C_L')
xlim([-15 15])
box on

set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend('Re = 650k','Re = 800k','Re = 1,000k');
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','northwest')

print -dpng -r300 Cl_vs_alpha

%% ---- Figure 4: CD vs Alpha ----
hFig = figure(4); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

% Plot drag curves (expect lower CD at higher Re, especially pre-stall)
plot(alpha650,  cd650)
plot(alpha800,  cd800)
plot(alpha1000, cd1000)

xlabel('Angle of attack \alpha (deg)')
ylabel('Drag coefficient C_D')
ylim([0 0.06])     % limit to relevant range for clarity
xlim([-15 15])
box on

set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend('Re = 650k','Re = 800k','Re = 1,000k');
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','north')

print -dpng -r300 Cd_vs_alpha

%% ---- Reset defaults ----
% Clean up so future plots aren’t affected
set(groot,'defaultAxesColorOrder','remove')
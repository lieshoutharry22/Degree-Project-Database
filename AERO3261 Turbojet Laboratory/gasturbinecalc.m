clear
clc
close all

% Defining constants
Pa = 103300;          % ambient pressure [Pa]
Ta = 297;             % ambient temperature [K]

gamma_c = 1.4;        % cold
gamma_h = 1.33;       % hot
cp_c = 1005;          % J/kgK
cp_h = 1156;          % J/kgK
R = 287;              % J/kgK

rho_fuel = 0.797;     % kg/L
LHV = 43.008e6;       % J/kg

A9 = 2500e-6;         % nozzle area [m^2] = 2500 mm^2
A3 = 3194e-6;         % compressor exit area [m^2]

file65k = '65k RPM';
file70k = '70k RPM';
file75k = '75k RPM';

data65k = readtable(file65k,'FileType','text','Delimiter','\t','VariableNamingRule','preserve');
data70k = readtable(file70k,'FileType','text','Delimiter','\t','VariableNamingRule','preserve');
data75k = readtable(file75k,'FileType','text','Delimiter','\t','VariableNamingRule','preserve');

$ Setting up empty data arrays

RPM = zeros(3,1);
Tm = zeros(3,1);

Tt2 = zeros(3,1);
Tt3 = zeros(3,1);
Tt4 = zeros(3,1);
Tt5 = zeros(3,1);

Pt2 = zeros(3,1);
Pt3 = zeros(3,1);
Pt4 = zeros(3,1);
Pt5 = zeros(3,1);

mdot_f = zeros(3,1);

% Extracting data

% 65K dataset 
RPM(1) = mean(data65k.("RPM"));
Tm(1) = mean(data65k.("Thrust (N)"));

Tt2(1) = mean(data65k.("T1 (C)")) + 273.15;
Tt3(1) = mean(data65k.("T2 (C)")) + 273.15;
Tt4(1) = mean(data65k.("T3 (C)")) + 273.15;
Tt5(1) = mean(data65k.("T4 (C)")) + 273.15;

Pt2(1) = mean(data65k.("P1 (kPa)"))*1000 + Pa;
Pt3(1) = mean(data65k.("P2 (kPa)"))*1000 + Pa;
Pt4(1) = mean(data65k.("P3 (kPa)"))*1000 + Pa;
Pt5(1) = mean(data65k.("P4 (kPa)"))*1000 + Pa;

mdot_f(1) = mean(data65k.("Fuel Flow  (L/hr)"))*rho_fuel/3600;

% 70K dataset 
RPM(2) = mean(data70k.("RPM"));
Tm(2) = mean(data70k.("Thrust (N)"));

Tt2(2) = mean(data70k.("T1 (C)")) + 273.15;
Tt3(2) = mean(data70k.("T2 (C)")) + 273.15;
Tt4(2) = mean(data70k.("T3 (C)")) + 273.15;
Tt5(2) = mean(data70k.("T4 (C)")) + 273.15;

Pt2(2) = mean(data70k.("P1 (kPa)"))*1000 + Pa;
Pt3(2) = mean(data70k.("P2 (kPa)"))*1000 + Pa;
Pt4(2) = mean(data70k.("P3 (kPa)"))*1000 + Pa;
Pt5(2) = mean(data70k.("P4 (kPa)"))*1000 + Pa;

mdot_f(2) = mean(data70k.("Fuel Flow  (L/hr)"))*rho_fuel/3600;

% 75K dataset 
RPM(3) = mean(data75k.("RPM"));
Tm(3) = mean(data75k.("Thrust (N)"));

Tt2(3) = mean(data75k.("T1 (C)")) + 273.15;
Tt3(3) = mean(data75k.("T2 (C)")) + 273.15;
Tt4(3) = mean(data75k.("T3 (C)")) + 273.15;
Tt5(3) = mean(data75k.("T4 (C)")) + 273.15;

Pt2(3) = mean(data75k.("P1 (kPa)"))*1000 + Pa;
Pt3(3) = mean(data75k.("P2 (kPa)"))*1000 + Pa;
Pt4(3) = mean(data75k.("P3 (kPa)"))*1000 + Pa;
Pt5(3) = mean(data75k.("P4 (kPa)"))*1000 + Pa;

mdot_f(3) = mean(data75k.("Fuel Flow  (L/hr)"))*rho_fuel/3600;

% Calculating additional values
pi_c = zeros(3,1);
eta_c = zeros(3,1);
eta_t = zeros(3,1);
eta_cc = zeros(3,1);

M9 = zeros(3,1);
T9 = zeros(3,1);
V9 = zeros(3,1);
rho9 = zeros(3,1);
mdot9 = zeros(3,1);
mdot_air = zeros(3,1);

FAR = zeros(3,1);
Tt = zeros(3,1);
deltaT = zeros(3,1);
TSFC = zeros(3,1);
Fs = zeros(3,1);

for i = 1:3

    % Compressor
    pi_c(i) = Pt3(i)/Pt2(i);
    Tt3s = Tt2(i)*(pi_c(i))^((gamma_c-1)/gamma_c);
    eta_c(i) = (Tt3s - Tt2(i)) / (Tt3(i) - Tt2(i));

    % Turbine
    Tt5s = Tt4(i)*(Pt5(i)/Pt4(i))^((gamma_h-1)/gamma_h);
    eta_t(i) = (Tt4(i) - Tt5(i)) / (Tt4(i) - Tt5s);

    % Stage 9 values
    M9(i) = sqrt((2/(gamma_h-1)) * ((Pt5(i)/Pa)^((gamma_h-1)/gamma_h) - 1));
    T9(i) = Tt5(i) * (Pa/Pt5(i))^((gamma_h-1)/gamma_h);
    V9(i) = M9(i)*sqrt(gamma_h*R*T9(i));
    rho9(i) = Pa/(R*T9(i));
    mdot9(i) = rho9(i)*A9*V9(i);

    % Mass flow
    mdot_air(i) = mdot9(i) - mdot_f(i);

    % Combustion chamber 
    eta_cc(i) = ((mdot_air(i)+mdot_f(i))*cp_h*(Tt4(i)-288) - mdot_air(i)*cp_c*(Tt3(i)-288)) ...
                / (mdot_f(i)*LHV);

    % FAR
    FAR(i) = mdot_f(i)/mdot_air(i);

    % Theoretical thrust
    Tt(i) = mdot9(i)*V9(i);

    % Difference in thrust
    deltaT(i) = ((Tt(i) - Tm(i))/Tm(i))*100;

    % TSFC
    TSFC(i) = (mdot_f(i)/Tm(i))*3600*10;

    % Specific thrust
    Fs(i) = Tm(i)/mdot_air(i);

end

$ Plots

%% define colours
myred       = [216 30 49]/255;
myblue      = [27 99 157]/255;
myblack     = [0 0 0]/255;
mygreen     = [0 128 0]/255;
mycyan      = [2 169 226]/255;
myyellow    = [251 194 13]/255;
mygray      = [89 89 89]/255;

set(groot,'defaultAxesColorOrder',[myblack;myblue;myred;mygreen;myyellow;mycyan;mygray]);

%% general plot parameters
alw = 1;        % axes line width
fsz = 20;       % font size
lw  = 2.5;      % line width
msz = 8;       % marker size

set(0,'defaultLineLineWidth',lw);
set(0,'defaultLineMarkerSize',msz);

%% 1) Low temperature group
hFig = figure(1); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(RPM,Tt2,'o')
plot(RPM,Tt3,'o')

xlabel('RPM')
ylabel('Total Temperature [K]')
box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend('T_{t1}','T_{t2}');
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','best')

print -dpng -r300 fig_temp_12

%% 2) High temperature group
hFig = figure(2); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(RPM,Tt4,'o')
plot(RPM,Tt5,'o')

xlabel('RPM')
ylabel('Total Temperature [K]')
box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend('T_{t3}','T_{t4}');
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','best')

print -dpng -r300 fig_temp_34

%% 3) Low pressure group
hFig = figure(3); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(RPM,Pt2/1000,'o')
plot(RPM,Pt5/1000,'o')

xlabel('RPM')
ylabel('Total Pressure [kPa]')
box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend('P_{t1}','P_{t4}');
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','best')

print -dpng -r300 fig_press_14

%% 4) High pressure group
hFig = figure(4); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(RPM,Pt3/1000,'o')
plot(RPM,Pt4/1000,'o')

xlabel('RPM')
ylabel('Total Pressure [kPa]')
box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend('P_{t2}','P_{t3}');
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','best')

print -dpng -r300 fig_press_23

%% 5) Efficiencies
hFig = figure(5); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(RPM,eta_c,'o')
plot(RPM,eta_t,'o')
plot(RPM,eta_cc,'o')

xlabel('RPM')
ylabel('Efficiency')
box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend('\eta_c','\eta_t','\eta_{cc}');
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','best')

print -dpng -r300 fig_efficiency

%% 6) Jet velocity
hFig = figure(6); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(RPM,V9,'o')

xlabel('RPM')
ylabel('Jet Velocity [m/s]')
box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

print -dpng -r300 fig_V9

%% 7) Thrust comparison
hFig = figure(7); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(RPM,Tm,'-o')
plot(RPM,Tt,'-o')

xlabel('RPM')
ylabel('Thrust [N]')
box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend('Measured','Theoretical');
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','best')

print -dpng -r300 fig_thrust

%% 8) Relative thrust difference
hFig = figure(8); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(RPM,deltaT,'o')

xlabel('RPM')
ylabel('\Delta T [%]')
box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

print -dpng -r300 fig_deltaT

%% 9) TSFC
hFig = figure(9); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(RPM,TSFC,'o')

xlabel('RPM')
ylabel('TSFC [kg/hr/daN]')
box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

print -dpng -r300 fig_TSFC

%% 10) Specific thrust
hFig = figure(10); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(RPM,Fs,'o')

xlabel('RPM')
ylabel('Specific Thrust [N/(kg/s)]')
box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

print -dpng -r300 fig_Fs

%% 11) FAR
hFig = figure(11); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(RPM,FAR,'o')

xlabel('RPM')
ylabel('FAR')
box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

print -dpng -r300 fig_FAR

%% 12) Air mass flow rate
hFig = figure(12); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

plot(RPM,mdot_air,'o')

xlabel('RPM')
ylabel('Air Mass Flow Rate [kg/s]')
box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

print -dpng -r300 fig_mdot_air

%% reset colours back to matlab defaults
set(groot,'defaultAxesColorOrder','remove')

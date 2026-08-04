%% Comments
% This m-file is based on the file created by Doug Auld for the
% Aerodynamics for students.
% Additional comments have been added by Dries Verstraete to assist students in
% AERO3261/AERO9261 with their propeller assignment

clear;
close all
clc

filenameprop        = 'testprop.mat';

%% Figure formatting
myred           = [216 30 49]/255;
myblue          = [27 99 157]/255;
myblack         = [0 0 0]/255;
mygreen         = [0 128 0]/255;
mycyan          = [2 169 226]/255;
myyellow        = [251 194 13]/255;
mygray          = [89 89 89]/255;
set(groot,'defaultAxesColorOrder',[myblack;myblue;myred;mygreen;myyellow;mycyan;mygray]);
alw             = 1;
fsz             = 26;
lw              = 2.5;
msz             = 40;
set(0,'defaultLineLineWidth',lw);
set(0,'defaultLineMarkerSize',msz);

%% Propeller Geometry
NrOfElements        = 10;
chordvec            = [0.027977708561170,0.034299337193030,0.040541124695084,0.042004958933961,0.040990306964408,0.037439210939482,0.031863802521854,0.025691220389119,0.019210138835871,0.012400000000000];
pitchvec            = [0.592297851384980,0.586232564288562,0.453445862572999,0.328392057826073,0.250049604543341,0.222368315231064,0.204746151583549,0.204091102397062,0.202581316907807,0.201000000000000];

dia                 = 0.5;
R                   = dia/2.0;

%% Blade segments
xs                  = 0.1*R;
xt                  = R;
rstep               = (xt-xs)/(NrOfElements-1);
r1                  = (xs:rstep:xt);

%% Atmospheric conditions
p                   = 82750;
T                   = 277.23;
R                   = 287.05;
rho                 = p/R/T;
mu                  = 1.458e-6*T^1.5/(T+110.4);

%% Flight speed settings
Vmin                = 0.1;
NrOfVs              = 100;

%% Aerodynamic data
filenames.Re50      = 'Naca4412Re50';
filenames.Re100     = 'Naca4412Re100';
filenames.Re200     = 'Naca4412Re200';
filenames.Re350     = 'Naca4412Re350';
filenames.Re500     = 'Naca4412Re500';
filenames.Re650     = 'Naca4412Re650';
filenames.Re800     = 'Naca4412Re800';
filenames.Re1000    = 'Naca4412Re1000';

astall = -8.9; % in DEGREES

[alfa500,cl500,cd500,Remat,Cd0mat] = dragdata(filenames,astall);

%% Convergence criterion
accuracy            = 1e-3;

%% RPM Loop
RPM_list = [2000, 5000, 8000];

Results(length(RPM_list)) = struct(...
    'RPM',[],'Vvec',[],'J',[],'t',[],'q',[],'eff',[],'T',[],'Q',[],'P',[]);

for rpmIdx = 1:length(RPM_list)

    RPM   = RPM_list(rpmIdx);
    n     = RPM / 60;
    omega = n * 2 * pi;
    Vmax  = n * dia * 1.5;
    Vvec  = linspace(Vmin, Vmax, NrOfVs);

    t      = NaN(1,length(Vvec));
    q      = NaN(1,length(Vvec));
    J      = NaN(1,length(Vvec));
    eff    = NaN(1,length(Vvec));
    icheck = NaN(1,length(Vvec));
    T      = NaN(1,length(Vvec));
    Q      = NaN(1,length(Vvec));
    P      = NaN(1,length(Vvec));

    for cntrV = 1:length(Vvec)
        V      = Vvec(cntrV);
        thrust = 0.0;
        torque = 0.0;

        for j = 1:size(r1,2)
            rad       = r1(j);
            pitch     = pitchvec(j);
            chord     = chordvec(j);
            theta     = pitch;
            a         = 0.1;
            b         = 0.01;
            finished  = false;
            sum       = 1;
            itercheck = 0;

            while (~finished)
                V0     = V*(1+a);
                V2     = omega*rad*(1-b);
                phi    = atan2(V0,V2);
                alpha  = theta-phi;
                Vlocal = sqrt(V0*V0+V2*V2);
                Re     = rho*chord*Vlocal/mu;

                cl = interp1(alfa500,cl500,alpha,'spline','extrap');
                cd = interp1(alfa500,cd500,alpha,'spline','extrap');
                if alpha < alfa500(1)
                    cl = cl500(1);
                    cd = cd500(1) + 0.1;
                elseif alpha > alfa500(end)
                    cl = cl500(end);
                    cd = cd500(end) + 0.1;
                end
                factor = interp1(Remat,Cd0mat,Re/500000,'spline','extrap');
                cd     = cd*factor;

                DtDr = 0.5*rho*Vlocal*Vlocal*2.0*chord*(cl*cos(phi)-cd*sin(phi));
                DqDr = 0.5*rho*Vlocal*Vlocal*2.0*chord*rad*(cd*cos(phi)+cl*sin(phi));

                tem1 = DtDr/(4.0*pi*rad*rho*V*V*(1+a));
                tem2 = DqDr/(4.0*pi*rad*rad*rad*rho*V*(1+a)*omega);
                anew = 0.5*(a+tem1);
                bnew = 0.5*(b+tem2);

                if (abs((anew-a)/a)<accuracy)
                    if (abs((bnew-b)/b)<accuracy)
                        finished = true;
                    end
                end
                a   = anew;
                b   = bnew;
                sum = sum+1;
                if (sum>500)
                    finished  = true;
                    itercheck = 1;
                end
            end

            if (isinf(a)||isinf(b)||isnan(a)||isnan(b))
                thrust = thrust;
                torque = torque;
            else
                thrust = thrust + DtDr*rstep;
                torque = torque + DqDr*rstep;
            end
        end

        t(cntrV)      = thrust/(rho*n*n*dia*dia*dia*dia);
        q(cntrV)      = torque/(rho*n*n*dia*dia*dia*dia*dia);
        J(cntrV)      = V/(n*dia);
        eff(cntrV)    = J(cntrV)/2.0/pi*t(cntrV)/q(cntrV);
        icheck(cntrV) = itercheck;
        T(cntrV)      = thrust;
        Q(cntrV)      = torque;
        P(cntrV)      = 2*pi*n*torque;
    end

    ind1       = find(t<0);
    ind2       = find(q<0);
    ind        = union(ind1,ind2);
    t(ind)     = []; q(ind)     = []; J(ind)    = [];
    eff(ind)   = []; icheck(ind)= []; T(ind)    = [];
    Q(ind)     = []; P(ind)     = []; Vvec(ind) = [];

    Results(rpmIdx).RPM  = RPM;
    Results(rpmIdx).Vvec = Vvec;
    Results(rpmIdx).J    = J;
    Results(rpmIdx).t    = t;
    Results(rpmIdx).q    = q;
    Results(rpmIdx).eff  = eff;
    Results(rpmIdx).T    = T;
    Results(rpmIdx).Q    = Q;
    Results(rpmIdx).P    = P;

end

%% Plotting
colors       = [myblack; myblue; myred];
legendLabels = {'2000 RPM','5000 RPM','8000 RPM'};

%% Figure 1 - Thrust & Power vs Flight Speed
hFig = figure(1); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

yyaxis left
for i = 1:3
    plot(Results(i).Vvec, Results(i).T, '--', 'Color', colors(i,:));
end
ylabel('Thrust (N)')
yyaxis right
for i = 1:3
    plot(Results(i).Vvec, Results(i).P, '-', 'Color', colors(i,:));
end
ylabel('Power (W)')
xlabel('Flight Speed V (m/s)')

box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend([strcat(legendLabels,' - T'), strcat(legendLabels,' - P')]);
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','best')

print -dpng -r300 Thrust_Power_vs_Speed

%% Figure 2 - Thrust & Power vs Advance Ratio
hFig = figure(2); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

yyaxis left
for i = 1:3
    plot(Results(i).J, Results(i).T,'--' , 'Color', colors(i,:));
end
ylabel('Thrust (N)')
yyaxis right
for i = 1:3
    plot(Results(i).J, Results(i).P, '-', 'Color', colors(i,:));
end
ylabel('Power (W)')
xlabel('Advance Ratio J')

box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend([strcat(legendLabels,' - T'), strcat(legendLabels,' - P')]);
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','best')

print -dpng -r300 Thrust_Power_vs_J

%% Figure 3 - Thrust Coefficient vs Advance Ratio
hFig = figure(3); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

for i = 1:3
    plot(Results(i).J, Results(i).t, 'Color', colors(i,:));
end
xlabel('Advance Ratio J')
ylabel('Thrust Coefficient C_t')

box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend(legendLabels);
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','best')

print -dpng -r300 Ct_vs_J

%% Figure 4 - Torque Coefficient vs Advance Ratio
hFig = figure(4); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

for i = 1:3
    plot(Results(i).J, Results(i).q, 'Color', colors(i,:));
end
xlabel('Advance Ratio J')
ylabel('Torque Coefficient C_q')

box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend(legendLabels);
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','best')

print -dpng -r300 Cq_vs_J

%% Figure 5 - Efficiency vs Advance Ratio
hFig = figure(5); clf
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
hold on
grid on
grid minor

for i = 1:3
    plot(Results(i).J, Results(i).eff, 'Color', colors(i,:));
end
xlabel('Advance Ratio J')
ylabel('Propulsive Efficiency \eta')

box on
set(gca,'LineWidth',alw)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle','-')
set(gca,'GridColor','k')
set(gca,'MinorGridColor','k')
set(findall(hFig,'-property','FontSize'),'FontSize',fsz)

hleg = legend(legendLabels);
set(hleg,'EdgeColor',hleg.Color)
set(hleg,'Location','best')

print -dpng -r300 Efficiency_vs_J

%% Reset default colour order
set(groot,'defaultAxesColorOrder','remove')

%% Save propeller data for engine matching
for i = 1:3
    Propeller.AdvanceRatio = Results(i).J;
    Propeller.Efficiency   = Results(i).eff;
    Propeller.ThrustCoeff  = Results(i).t;
    Propeller.TorqueCoeff  = Results(i).q;
    Propeller.Diameter     = dia;
    Propeller.RotSpeed     = Results(i).RPM / 60;
    fname = sprintf('testprop_%dRPM.mat', Results(i).RPM);
    save(fname, 'Propeller');
end
%% Assignment 1
clear
close all
clc

%% ---- Flight conditions ----
% Set atmospheric conditions (given / assumed)
FlightCondition.Pressure        = 82750;      % Pa
FlightCondition.Temperature     = 277.23;     % K
FlightCondition.R               = 287.05;     % gas constant
FlightCondition.Density         = FlightCondition.Pressure / ...
                                  FlightCondition.R / FlightCondition.Temperature;
FlightCondition.g               = 9.80665;    % gravity
FlightCondition.RoC             = 0;          % zero for cruise case

%% ---- Speed range ----
% Define flight speed sweep for performance analysis
FlightSpeed.Vmin                = 5;
FlightSpeed.Vmax                = 50;
FlightSpeed.NrOfPoints          = 100;
FlightCondition.Speed           = linspace(FlightSpeed.Vmin,FlightSpeed.Vmax,FlightSpeed.NrOfPoints);

% Basic aero quantities
FlightCondition.ClimbAngle      = atan(FlightCondition.RoC ./ FlightCondition.Speed);
FlightCondition.DynPressure     = 0.5 * FlightCondition.Density .* FlightCondition.Speed.^2;

%% ---- Aircraft properties ----
Aircraft.m                      = 13.5;      % mass
Aircraft.Sw                     = 0.55;      % wing area
Aircraft.b                      = 2.9;       % wingspan
Aircraft.CD0                    = 0.03;      % parasitic drag
Aircraft.e                      = 0.90;      % Oswald efficiency
Aircraft.CL0                    = 0.28;      % lift offset

Aircraft.W                      = Aircraft.m * FlightCondition.g;
Aircraft.AR                     = Aircraft.b^2 / Aircraft.Sw;

% Lift coefficient required for steady flight
FlightCondition.CL              = Aircraft.W ./ ...
    (FlightCondition.DynPressure .* Aircraft.Sw .* cos(FlightCondition.ClimbAngle));

% Drag model (parabolic)
FlightCondition.CD              = Aircraft.CD0 + ...
    (FlightCondition.CL - Aircraft.CL0).^2 ./ (pi * Aircraft.e * Aircraft.AR);

% Drag force
FlightCondition.Drag            = FlightCondition.DynPressure .* Aircraft.Sw .* FlightCondition.CD;

% Thrust required (includes climb term, currently zero here)
FlightCondition.ThrustRequired  = FlightCondition.Drag + ...
    Aircraft.W * sin(FlightCondition.ClimbAngle);

%% ---- Climb case ----
% Introduce a realistic climb rate to compare against cruise
RoC_climb = 2.5;   

ClimbAngle_climb = atan(RoC_climb ./ FlightCondition.Speed);

% Extra thrust needed due to climb
ThrustRequired_climb = FlightCondition.Drag + ...
    Aircraft.W .* sin(ClimbAngle_climb);

%% ---- Quick sanity check ----
figure
plot(FlightCondition.Speed,FlightCondition.ThrustRequired,'LineWidth',2)
xlabel('Flight Speed (m/s)')
ylabel('Required Thrust (N)')
title('Required Thrust vs Flight Speed')
grid on
box on

%% ---- Propeller data ----
% Load precomputed propeller performance map
filenameprop = 'testprop_8000RPM.mat';
load(filenameprop)

n       = Propeller.RotSpeed;
dia     = Propeller.Diameter;
J       = Propeller.AdvanceRatio;
eff     = Propeller.Efficiency;
t       = Propeller.ThrustCoeff;
q       = Propeller.TorqueCoeff;
rho     = FlightCondition.Density;

% Remove invalid entries (important for interpolation)
valid   = ~(isnan(J) | isnan(t) | isnan(q) | isnan(eff));
J       = J(valid);
t       = t(valid);
q       = q(valid);
eff     = eff(valid);

%% ---- Engine data ----
% Get engine performance map
[EngineMap] = EngineSpecs(FlightCondition);

ThrottleVec = [40 60 80 100];
MatchAll    = struct();

%% ---- Engine–prop matching ----
% Loop through throttle settings
for iT = 1:length(ThrottleVec)

    Throttle = ThrottleVec(iT);

    % Select appropriate engine row (scaled by throttle)
    RPMrow          = round(length(EngineMap.RPM(:,1))/100 * Throttle);
    EngineRPM       = EngineMap.RPM(RPMrow,:);
    EnginePower     = EngineMap.Power(RPMrow,:);
    EngineTorque    = EngineMap.Torque(RPMrow,:);
    EngineFuelFlow  = EngineMap.FuelFlow(RPMrow,:);
    nvec            = EngineRPM / 60;

    % Preallocate
    Match.J                 = NaN(size(FlightCondition.Speed));
    Match.RPM               = NaN(size(FlightCondition.Speed));
    Match.EngineTorque      = NaN(size(FlightCondition.Speed));
    Match.PropTorque        = NaN(size(FlightCondition.Speed));
    Match.PropThrust        = NaN(size(FlightCondition.Speed));
    Match.PropEff           = NaN(size(FlightCondition.Speed));
    Match.EngineFuelFlow    = NaN(size(FlightCondition.Speed));
    Match.EnginePower       = NaN(size(FlightCondition.Speed));
    Match.PropPower         = NaN(size(FlightCondition.Speed));
    indQmatch               = NaN(size(FlightCondition.Speed));

    for cntrV = 1:length(FlightCondition.Speed)

        Speed = FlightCondition.Speed(cntrV);

        % Convert speed to advance ratio
        Jvec = Speed ./ (nvec * dia);

        % Interpolate prop performance at this J
        qvec   = interp1(J,q,Jvec,'spline',NaN);
        tvec   = interp1(J,t,Jvec,'spline',NaN);
        effvec = interp1(J,eff,Jvec,'spline',NaN);

        % Convert coefficients to physical values
        PropTorqueVec = qvec .* rho .* nvec.^2 .* dia.^5;
        PropThrustVec = tvec .* rho .* nvec.^2 .* dia.^4;
        PropPowerVec  = 2*pi*nvec .* PropTorqueVec;

        % Find operating point where prop torque matches engine torque
        ind = find(PropTorqueVec >= EngineTorque,1,'first');

        if ~isempty(ind)
            indQmatch(cntrV)            = ind;
            Match.J(cntrV)              = Jvec(ind);
            Match.RPM(cntrV)            = EngineRPM(ind);
            Match.EngineTorque(cntrV)   = EngineTorque(ind);
            Match.PropTorque(cntrV)     = PropTorqueVec(ind);
            Match.PropThrust(cntrV)     = PropThrustVec(ind);
            Match.PropEff(cntrV)        = effvec(ind);
            Match.EngineFuelFlow(cntrV) = EngineFuelFlow(ind);
            Match.EnginePower(cntrV)    = EnginePower(ind);
            Match.PropPower(cntrV)      = PropPowerVec(ind);
        end
    end

    % Remove invalid matches (edge artefacts)
    indFalse = find(indQmatch == 1);
    Match.J(indFalse)              = NaN;
    Match.RPM(indFalse)            = NaN;
    Match.EngineTorque(indFalse)   = NaN;
    Match.PropTorque(indFalse)     = NaN;
    Match.PropThrust(indFalse)     = NaN;
    Match.PropEff(indFalse)        = NaN;
    Match.EngineFuelFlow(indFalse) = NaN;
    Match.EnginePower(indFalse)    = NaN;
    Match.PropPower(indFalse)      = NaN;

    MatchAll(iT).Throttle = Throttle;
    MatchAll(iT).Match    = Match;
end

%% ---- Performance extraction ----
% Extract key performance metrics (cruise + climb)
nT = length(ThrottleVec);

Vmax_all       = NaN(nT,1);
RoCmax_all     = NaN(nT,1);
Vbest_all      = NaN(nT,1);

for iT = 1:nT

    Match = MatchAll(iT).Match;

    % Cruise condition (where thrust balances drag)
    diffThrust = Match.PropThrust - FlightCondition.ThrustRequired;
    indCruise  = find(diffThrust >= 0,1,'last');

    if ~isempty(indCruise)
        Vmax_all(iT) = FlightCondition.Speed(indCruise);
    end

    % Rate of climb calculation
    RoCvec = (Match.PropThrust - FlightCondition.Drag) .* ...
              FlightCondition.Speed / Aircraft.W;

    [RoCmax_all(iT), indRoC] = max(RoCvec);

    Vbest_all(iT) = FlightCondition.Speed(indRoC);
end

%% ---- Plot: thrust vs speed ----
% Shows where T = D and climb margin
for iT = 1:length(ThrottleVec)
    plot(FlightCondition.Speed, MatchAll(iT).Match.PropThrust)
end
plot(FlightCondition.Speed, FlightCondition.Drag, '--')
plot(FlightCondition.Speed, ThrustRequired_climb, '--')

xlabel('Flight Speed (m/s)')
ylabel('Thrust (N)')
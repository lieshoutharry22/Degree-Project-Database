function [C,y] = NumericalFTCS_time(time)
% Input :
%    time - time of soln (in hours)
% Output :
%    C - concentration of phosphate (mg/L) at end time
% Set constant paramters
L = 2000;                    % Length of domain [m]
D = 5E-1;                   % Diffusivity constant [m^2/s]
ny = 129;
tf = time*60*60;
% Generate range
y = linspace(0,L,ny);      % y-domain [m]
dy = y(2) - y(1);          % y-domain spacing [m]

MStepSize = ((L/256)^2)/(2*D);
% Calculate number of timesteps
nt = tf/MStepSize;

% Coefficient in FTCS scheme
sigma = D*MStepSize/dy^2;

% Initialise solution vecto1
Cn = zeros(1,ny);              % Solution at timestep n
Cnp1 = zeros(1,ny);            % Solution at timestep n+1

% Set initial condition and boundary conditions
Cn(2:end) = 0;
Cn(1) = 1;
Cnp1(1) = 1;
Cnp1(end) = 0;

% Loop over time
for n = 2:nt
    % Loop over internal points
    for i = 2:ny-1
        Cnp1(i) = Cn(i) + sigma*(Cn(i-1) - 2*Cn(i) + Cn(i+1));
        Cnp1(ny) = Cn(ny) + sigma*(Cn(ny-1) - Cn(ny));
    end
    % Update solution for next timestep
    Cn = Cnp1;
end

% Set output
C = Cnp1;
return

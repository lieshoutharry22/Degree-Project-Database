function [C, y] = Analytical(ny)
% Input:
%   ny - number of points in domain
% Output:
%   C - concentration of phosphate (mg/L) at end time
%   y - spatial domain

% Set constant parameters
L = 2000;        % Length of domain [m]
D = 5E-1;        % Diffusivity constant [m^2/s]
tf = 30*60;      % End time [s]

% Generate spatial domain
y = linspace(0, L, ny);  % y-domain [m]
dy = y(2) - y(1);        % y-domain spacing [m]

% Calculate analytical solution
C = zeros(1, ny);
for i = 1:ny
    C(i) = 1 - erf(y(i) / (2 * sqrt(D * tf)));
end

end
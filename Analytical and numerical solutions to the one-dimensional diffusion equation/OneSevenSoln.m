% Parameters
L = 2000; % Length of the mine shaft in meters
D = 5e-1; % Diffusivity of CO in m^2/s
t = 0; % t=0 as specified in Q1 .7
nfArray = [1 , 10 , 100]; % Three different nf values
% Spatial domain (129 Spatial grid points )
x = linspace (0 , L , 129) ;
% Modified initial condition
Y_initial = @( x ) x / L - 1;
% Plot the initial condition
figure ;
plot (x , Y_initial ( x ) , 'k', 'LineWidth', 2 , 'DisplayName', 'Mine shaft') ;
hold on ;
for nf = nfArray
    % Fourier series approximation
    Y_approx = zeros ( size ( x ) ) ;
    for n = 1: nf
        % Calculating Bn
        Bn = -2 * mod(n,2) / ( n * pi ) ;
        % Fourier series term contribution at tEnd
        Y_approx = Y_approx + Bn * sin ( n * pi * x / L ) * exp ( - D * ( n * pi / L ) ^2 * t ) ;
    end
    plot (x , Y_approx , 'DisplayName', sprintf('Number of Fourier terms = %d', nf ) ) ;
end
% Customize the plot
title ('Fourier Series Approximation of the Modified Initial Condition ') ;
xlabel ('Position along the mine shaft , x (m)') ;
ylabel ('CO Mass Fraction , C_H (x ,0) ') ;
legend ('show', 'Location', 'north') ;
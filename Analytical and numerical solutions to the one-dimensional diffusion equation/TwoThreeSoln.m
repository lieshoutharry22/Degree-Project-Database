clear all
clc
close all


nx = [9,17,33,65,129.];

for i = 1:length(nx)
    x = linspace(0,2000,nx(i));      % x-domain [m]
    Y = NumericalFTCS(nx(i));
    
    Plot = figure(1);
    plot(x,Y)
    title ('CO Mass Fraction in the Mine Shaft at t = 30 min ') ;
    xlabel ('Distance along the shaft (m)') ;
    ylabel ('CO Mass Fraction ') ;
    legend ('Nx = 9', 'Nx = 17 ', 'Nx = 33 ', 'Nx = 65 ', 'Nx = 129 ', ' Analytical Solution ', 'Location ', 'best ') ;

    hold on
end


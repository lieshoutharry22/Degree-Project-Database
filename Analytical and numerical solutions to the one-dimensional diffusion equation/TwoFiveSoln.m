clear all
clc
close all


time = [1,6,12,18,24];

for i = 1:length(time)
    x = linspace(0,2000,129);      % x-domain [m]
    Y = NumericalFTCS_time(time(i));
    
    Plot = figure(1);
    plot(x,Y)
    title ('CO Mass Fraction in the Mine Shaft at Different Times ') ;
    xlabel ('Distance along the shaft (m)') ;
    ylabel ('CO Mass Fraction , C') ;
    legend ('1 hours', '6 hours', '12 hours', '18 hours', '24 hours', 'Location ', 'best ') ;

    hold on
end



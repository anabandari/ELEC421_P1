%%%%%%%%%%%%%%%%%%%%%
% Project 1 Main
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
% Part 1 - Section 2 Plots
%%%%%%%%%%%%%%%%%%%%%
%{ Running Naive DFT For different values of n
% exec_t = zeros(100,1);
% samples = linspace(1000,50000,100);
% for n=1:100
%   exec_t(n) = timeit(@() NaiveDFT(n)); % Measure execution time for Naive DFT
% end
% 
% figure;
% hold on 
% plot(samples, exec_t); 
% extimecomp = 25e-14*samples.^2;
% plot(samples, extimecomp);
% hold off
%}
clc;

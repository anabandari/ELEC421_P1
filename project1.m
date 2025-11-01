%%%%%%%%%%%%%%%%%%%%%
% Project 1 Main
%%%%%%%%%%%%%%%%%%%%%

% Prompt user for custom signal
prompt = "Would you like to enter your own analog sinusoidal signal? Y/N [Y]: ";
txt = input(prompt,"s");
if isempty(txt); txt = 'Y'; end

if upper(txt) == 'Y'
    % User-defined sinusoidal parameters
    A = input('Enter amplitude (e.g., 5): ');
    if isempty(A); A = 5; end
    
    f = input('Enter frequency in Hz (e.g., 100): ');
    if isempty(f); f = 100; end
    
    Fs = input('Enter sampling frequency in Hz (e.g., 1000): ');
    if isempty(Fs); Fs = 1000; end
    
    T = input('Enter duration in seconds (e.g., 1): ');
    if isempty(T); T = 1; end
    
    t = 0:1/Fs:T-1/Fs; % time vector
    Signal = A * sin(2*pi*f*t);
else
    % Default signal
    Fs = 1000; T = 1; t = 0:1/Fs:T-1/Fs;
    Signal = 5 * sin(200*pi*t);
end

%%%%%%%%%%%%%%%%%%%%%
% Part 1 - Section 2 Plots
%%%%%%%%%%%%%%%%%%%%%
% Running Naive DFT for different values of N

N_low = 1;
N_high = 100;
N_samples = 10;

N_values = linspace(N_low, N_high, N_samples); % example sample sizes
exec_t = zeros(length(N_values),1);

for idx = 1:length(N_values)
    N = round(N_values(idx));
    
    % Ensure Signal length matches N (truncate or pad)
    if length(Signal) < N
        x = [Signal zeros(1, N - length(Signal))];
    else
        x = Signal(1:N);
    end
    
    exec_t(idx) = timeit(@() naive_dft(x, N)); % Measure execution time
end

figure;
hold on
plot(N_values, exec_t, 'LineWidth', 1.5);
xlabel('Number of samples N');
ylabel('Execution Time (s)');
title('Execution Time of Naive DFT vs Sample Size');
grid on
hold off

clc;

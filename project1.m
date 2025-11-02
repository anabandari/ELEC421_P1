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
exec_t_ndft = zeros(length(N_values),1);

for idx = 1:length(N_values)
    N = round(N_values(idx));
    
    % Ensure Signal length matches N (truncate or pad)
    if length(Signal) < N
        x = [Signal zeros(1, N - length(Signal))];
    else
        x = Signal(1:N);
    end
    
    exec_t_ndft(idx) = timeit(@() naive_dft(x, N)); % Measure execution time
end



figure;
hold on
plot(N_values, exec_t_ndft, 'LineWidth', 1.5);
xlabel('Number of samples N', 'FontSize', 20);
ylabel('Execution Time (s)', 'FontSize', 20);
title('Execution Time of Naive DFT vs Sample Size', 'FontSize', 20);
grid on
hold off

clc;

%%%%%%%%%%%%%%%%%%%%%
% Part 3.2 - Compare DIT FFT and DIF FFT
%%%%%%%%%%%%%%%%%%%%%

% 3.2.1 Compare DIT FFT and DIF FFT for the same input signal

% We will compare the two FFTs against the first ten discrete values of
% our input signal

x_n_test = Signal(1:10).'; 

X_dit_fft = dit_fft(x_n_test);
X_dif_fft = bitrevorder(dif_fft(x_n_test)); % order of X_dif_fft is bit reversed
                                            % to match order of X_dit_fft

% using a table, we can see that the values of the 
% DIT FFT and DIF FFT are identical
compareFFTTable = table(X_dit_fft, X_dif_fft);
disp(compareFFTTable) % table is printed to console

% 3.2.1 Compare computational cost of DIT FFT and DIF FFT 

% create set of 2^n values to iterate through while testing ex_time
n2_values = 2:10;
N2_values = 2 .^ n2_values;

% vars to hold ex_time values
exec_t_dif = zeros(length(N2_values), 1);
exec_t_dit = zeros(length(N2_values), 1);

for idx = 1:length(N2_values)
    N2 = N2_values(idx);

    % Ensure Signal length matches N (truncate or pad)
    if length(Signal) < N2
        x = [Signal(:); zeros(N2 - length(Signal), 1)];
    else
        x = Signal(1:N2).';
    end

     % Measure DIF FFT and DIT FFT execution time
    exec_t_dif(idx) = timeit(@() dif_fft(x));
    exec_t_dit(idx) = timeit(@() dit_fft(x));
end

figure;
loglog(N2_values, exec_t_dif, '-s', 'LineWidth', 1.5, 'DisplayName', 'DIF FFT');
hold on;
loglog(N2_values, exec_t_dit, '-s', 'LineWidth', 1.5, 'DisplayName', 'DIT FFT');
xlabel('Number of Samples (N)', 'FontSize', 20);
ylabel('Execution Time (seconds)', 'FontSize', 20);
title('Execution Time of DIF FFT vs DIT FFT', 'FontSize', 20);
legend('Location','northwest');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%
% Part 4 - Frequency Analysis of Audio Signal
%%%%%%%%%%%%%%%%%%%%%

% 4.1: Load audio signals

clean_signal_filename = "clean_signal.wav";
noisy_signal_filename = "noisy_signal.wav";

% use audioread function to convert .WAV file to data values
% source: https://www.mathworks.com/help/matlab/ref/audioread.html
[cs_data, cs_sample_rate] = audioread(clean_signal_filename);
[ns_data, ns_sample_rate] = audioread(noisy_signal_filename);

% extra feature: apply bandpass, Butterworth filter to noisy audio file
% used filterDesigner tool from L15 to create filter with:
ns_filtered_data = filter(Hd2, ns_data);
audiowrite('noisy_signal_filtered.wav', ns_filtered_data, ns_sample_rate);
% sound(ns_filtered_data, ns_sample_rate); % to test sound in MATLAB

% 4.2: Apply FFT and Plot Spectrum

% take DIF FFT of clean signal and noisy signal
cs_X = bitrevorder(dif_fft(cs_data));
ns_X = bitrevorder(dif_fft(ns_data));
nsf_X = bitrevorder(dif_fft(ns_filtered_data));

% find lengths of DIF FFTs for future use
cs_N = length(cs_X);
ns_N = length(ns_X);
nsf_N = length(nsf_X);

% find values in the magnitude spectrum
cs_mag = abs(cs_X);
ns_mag = abs(ns_X);
nsf_mag = abs(nsf_X);

% find corresponding frequencies
cs_f = (0:cs_N-1)*(cs_sample_rate/cs_N);
ns_f = (0:ns_N-1)*(ns_sample_rate/ns_N);
nsf_f = (0:nsf_N-1)*(ns_sample_rate/nsf_N);

% plot magnitude spectrum vs its corresponding frequencies
figure;
plot(cs_f, cs_mag, 'b', 'LineWidth', 1);
title('Clean Signal Magnitude Spectrum', 'FontSize', 20);
xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
grid on;

figure;
plot(cs_f, cs_mag, 'b', 'LineWidth', 1);
title('Zoomed In Clean Signal Magnitude Spectrum', 'FontSize', 20);
xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
ylim([-0.1, 10]);
xlim([0, 4000]);
grid on;

figure;
plot(ns_f, ns_mag, 'b', 'LineWidth', 1);
title('Noisy Signal Magnitude Spectrum', 'FontSize', 20);
xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
grid on;

figure;
plot(ns_f, ns_mag, 'b', 'LineWidth', 1);
title('Zoomed In Noisy Signal Magnitude Spectrum', 'FontSize', 20);
xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
ylim([-0.1, 70]);
xlim([0, 4000]);
grid on;

figure;
plot(nsf_f, nsf_mag, 'b', 'LineWidth', 1);
title('Filtered Noisy Signal Magnitude Spectrum', 'FontSize', 20);
xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
ylim([-0.1, 70]);
xlim([0, 4000]);
grid on;

figure;
plot(nsf_f, nsf_mag, 'b', 'LineWidth', 1);
title('Zoomed-In Filtered Noisy Signal Magnitude Spectrum', 'FontSize', 20);
xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
ylim([-0.1, 500]);
xlim([0, 650]);
grid on;

% 4.3: Identify Dominant Frequencies

% the "dominant frequency" is taken as the frequency in the signal with
% the greatest amplitude

% we will find the max amplitude in each signal using the max fn
[cs_max, cs_max_index] = max(cs_mag);
[ns_max, ns_max_index] = max(ns_mag);

% use max index to find the dominant frequency 
cs_max_freq = cs_f(cs_max_index); % 439.9414 Hz
ns_max_freq = ns_f(ns_max_index); % 439.9414 Hz





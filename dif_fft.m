%%%%%%%%%%%%%%%%%%%%%
% Part 3 - Decimation-in-Frequency (DIF) FFT Implementation
%%%%%%%%%%%%%%%%%%%%%

% x_sample = [ones(5,  1)]; % sample signal for now
x_sample = [5,3,2,1,0,0].';

% declare global counter value to use in function
global counter;
counter = 0;

% call function
[X_k, X_k_bit_reversed, X_k_FFT, DIF_diff, t] = dif_fft_func(x_sample);


% FUNCTION: Takes DIT FFT of discrete input signal
% param: x_n, discrete input signal
% return: 
%   X_k, DIF FFT discrete output signal
%   X_k_bit_reversed, bit reversed DIF FFT discrete output signal
%   X_k_fft, MATLAB FFT discrete output signal
%   DIF_diff, difference between DIF FFT and MATLAB FFT outputs
function [X_k, X_k_bit_reversed, X_k_FFT, DIF_diff, t] = dif_fft_func(x_n)
tic;

global counter;
if (counter == 0) % only check on first recursive loop
    x_n_len = length(x_n);

    % check that signal is a power of 2
    % if not, zero pad x[n] to power of 2
    if (mod(x_n_len, log2(x_n_len)) ~= 0)
        next_exp_2 = ceil(log2(x_n_len));
        zeros_needed = [zeros(2^next_exp_2 - x_n_len, 1)];
        x_n = vertcat(x_n, zeros_needed);
    end

    % check that input signal is not empty
    if (x_n_len == 0)
        error('Input signal is empty!')
    end

    % increment counter only once so if block doesn't run again
    counter = counter + 1;
end


% check for base case N = 1
if (length(x_n) == 1)
    X_k = x_n;
    return 
end

N = length(x_n);

% "reverse" butterfly
% step 1: add the two elements
X_even = x_n(1:N/2) + x_n(N/2 + 1: N);

% step 2: subtract the two elements 
X_odd = x_n(1:N/2) - x_n(N/2 + 1: N);

% step 3: multiply by the twiddle factor
% 3.a: find twiddle factor
n = (0:N/2-1);
w_N = exp(-1j * n * 2 * pi / N).'; % convert to column vector

% 3.b multiply
X_odd = X_odd .* w_N; % multiplication by element in matlab

% Recursion
X_even = dif_fft_func(X_even);
X_odd = dif_fft_func(X_odd);

% Compile results
X_k = [X_even; X_odd];

t = toc;

X_k_bit_reversed = bitrevorder(X_k);
X_k_FFT = fft(x_n);
DIF_diff = abs(X_k_bit_reversed) - abs(X_k_FFT);

end


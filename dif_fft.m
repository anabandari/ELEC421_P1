%%%%%%%%%%%%%%%%%%%%%
% Part 3 - Decimation-in-Frequency (DIF) FFT Implementation
%%%%%%%%%%%%%%%%%%%%%


% FUNCTION: Takes DIF FFT of discrete input signal
% param: x_n, discrete input signal
% return: X_k, DIF FFT discrete output signal
function X_k = dif_fft(x_n)

x_n_len = length(x_n);

if (mod(x_n_len, log2(x_n_len)) ~= 0)
    next_exp_2 = ceil(log2(x_n_len));
    zeros_needed = [zeros(2^next_exp_2 - x_n_len, 1)];
    x_n = vertcat(x_n, zeros_needed);
end

% check for recursion base case N = 1
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
X_odd = X_odd .* w_N; % multiplication by element in MATLAB

% Recursion
X_even = dif_fft(X_even);
X_odd = dif_fft(X_odd);

% Compile results
X_k = [X_even; X_odd];

% Note: bit reversal is not done within this function but is applied when
% the function is invoked in the main project1.m file

end


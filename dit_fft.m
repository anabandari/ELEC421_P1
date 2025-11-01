%%%%%%%%%%%%%%%%%%%%%
% Part 2 - Decimation-in-Time (DIT) FFT Implementation
% 
% With reference from https://youtu.be/Ty0JcR6Dvis?si=RcP4nY_l-D2J5zA-
%%%%%%%%%%%%%%%%%%%%%

function X = dit_fft_test(x)
    x = x(:);  % ensure column vector
    N = length(x);

    % If N is not a power of 2, pad with zeros
    if mod(log2(N), 1) ~= 0
        N_next = 2^nextpow2(N);
        x = [x; zeros(N_next - N, 1)];
        N = N_next;
    end

    % Base case
    if N == 1
        X = x;
        return;
    end

    % Split into even and odd indices
    x_even = x(1:2:end);
    x_odd  = x(2:2:end);

    % Recursive calls
    X_even = dit_fft_test(x_even);
    X_odd  = dit_fft_test(x_odd);

    % Compute N/2 twiddle factors
    W = exp(-1i * 2 * pi * (0:(N/2 - 1))' / N);

    % Butterfly combination
    X = zeros(N,1);
    X(1:N/2)     = X_even + W .* X_odd;
    X(N/2+1:N)   = X_even - W .* X_odd;
end

test1 = fft([5,3,2,1,0,0])

test2 = dit_fft_test([5,3,2,1,0,0])
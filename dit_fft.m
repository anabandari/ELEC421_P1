%%%%%%%%%%%%%%%%%%%%%
% Part 2 - Decimation-in-Time (DIT) FFT Implementation
%%%%%%%%%%%%%%%%%%%%%


function X = dit_fft(x)
    x = x(:);  % ensure column vector
    N = length(x);

    % If N is not a power of 2, pad with zeros
    if mod(log2(N),1) ~= 0
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
    X_even = ditfft_butterfly_pad_twiddle(x_even);
    X_odd  = ditfft_butterfly_pad_twiddle(x_odd);

    % Compute twiddle factors using separate function
    W = twiddle(N);

    % Butterfly combination
    X = zeros(N,1);
    X(1:N/2)     = X_even + W .* X_odd;
    X(N/2+1:N)   = X_even - W .* X_odd;
end

function W = twiddle(N)
    % Compute the DIT FFT twiddle factors for N-point FFT
    k = (0:N/2-1).';
    W = exp(-1j*2*pi*k/N);
end
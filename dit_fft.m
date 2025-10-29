%%%%%%%%%%%%%%%%%%%%%
% Part 2 - Decimation-in-Time (DIT) FFT Implementation
%%%%%%%%%%%%%%%%%%%%%
function decintime = dit_fft(Signal, Samples)
    decintime = zeros(1, Samples);
    ditodd = Signal(2:2:Samples);
    diteven = Signal(1:2:Samples);
    twiddlehalf = Twiddle(Samples/2);
    twiddlefull = Twiddle(Samples);
    
    for k=1:Samples
        for r=1:(Samples/2)
            decintime(k) = decintime(k) + diteven(r)*twiddlehalf^(r*k) + twiddlefull^(k)*ditodd(r)*twiddlehalf^(r*k);
        end 
    end
end 

function twid = Twiddle(N)
    twid = zeros(1,N);
    for k = 0:N-1
        twid(k+1) = exp(-2*pi*1i*k/N); % Compute the twiddle factors
    end
end 
%%%%%%%%%%%%%%%%%%%%%
% Part 2 - Decimation-in-Time (DIT) FFT Implementation
%%%%%%%%%%%%%%%%%%%%%
function decintime = dit_fft(ditsignal)
    ditodd = ditsignal(2:2:end);
    diteven = ditsignal(1:2:end);
    twiddlehalf = Twiddle(size(ditsignal/2));
    twiddlefull = Twiddle(size(ditsignal));
    
    for k=0:size(ditsignal)
        for r=0:size(ditsignal/2)
            diteven(r)*Twiddlehalf(r)
        end 
    end
end 

function twid = Twiddle(N)
    twid = zeros(1,N);
    for k = 0:N-1
        twid(k+1) = exp(-2*pi*1i*k/N); % Compute the twiddle factors
    end
end 
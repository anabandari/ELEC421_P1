%%%%%%%%%%%%%%%%%%%%%
% Part 1 - Naive DFT Implementation and Cost Analysis
%%%%%%%%%%%%%%%%%%%%%
function ndft = naive_dft(Samples, Signal)
ndft = zeros(1, Samples); % Initialize the output array
    for k = 1:Samples   
        for n = 0:Samples
            ndft(k) = ndft(k) + Signal*exp(-2*pi*1i*k*n/Samples);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%
% Part 1 - Naive DFT Implementation and Cost Analysis
%
% Inputs: A discrete signal in terms of n, Number of samples 'N' 
% Outputs: The Naive DFT
% 
%%%%%%%%%%%%%%%%%%%%%

function ndft = naive_dft(Signal, Samples)
    ndft = zeros(1, Samples); % Initialize the output array
    for k = 1:Samples
        for n = 1:Samples
            ndft(k) = ndft(k) + Signal(n) * exp(-1i * 2 * pi * (k-1) * (n-1) / Samples);
        end
    end
end



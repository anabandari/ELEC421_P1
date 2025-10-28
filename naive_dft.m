%%%%%%%%%%%%%%%%%%%%%
% Part 1 - Naive DFT Implementation and Cost Analysis
%
% Inputs: A discrete signal in terms of n, Number of samples 'N' 
% Outputs: The Naive DFT
% 
%%%%%%%%%%%%%%%%%%%%%

function ndft = naive_dft(Samples, Signal)
ndft = zeros(1, Samples); % Initialize the output array
    for k = 1:Samples   
        for n = 1:Sampless
            ndft(k) = ndft(k) + Signal(n)*exp(-1i*k*(2*pi/Samples)*n);
        end
    end
end


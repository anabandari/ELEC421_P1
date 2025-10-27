%%%%%%%%%%%%%%%%%%%%%
% Part 1 - Naive DFT Implementation and Cost Analysis
%%%%%%%%%%%%%%%%%%%%%
function ndft = NaiveDFT(n)
ndft = zeros(1, n); % Initialize the output array
for k = 0:n-1
    for m = 0:n-1
        ndft(k+1) = ndft(k+1) + exp(-2*pi*1i*k*m/n); % Compute the DFT
    end
end
end


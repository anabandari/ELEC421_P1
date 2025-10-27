%%%%%%%%%%%%%%%%%%%%%
% Project 1
%%%%%%%%%%%%%%%%%%%%%

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

% Running Naive DFT For different values of n
exec_t = zeros(10,1);
samples = 1000;
for n=0:10
    exec_t(n+1) = timeit(@() NaiveDFT(samples)); % Measure execution time for Naive DFT
    samples = samples+1000;
end

figure;
subplot()

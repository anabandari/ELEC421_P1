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
exec_t = zeros(100,1);
samples = linspace(1000,50000,100);
for n=1:100
    exec_t(n) = timeit(@() NaiveDFT(n)); % Measure execution time for Naive DFT
end

figure;
hold on 
plot(samples, exec_t); 
extimecomp = 25e-14*samples.^2;
plot(samples, extimecomp);
hold off

clc;

%%%%%%%%%%%%%%%%%%%%%
% Part 2 - Decimation-in-Time (DIT) FFT Implementation
%%%%%%%%%%%%%%%%%%%%%
function twid = Twiddle(N)
    twid = zeros(1,N);
    for k = 0:N-1
        twid(k+1) = exp(-2*pi*1i*k/N); % Compute the twiddle factors
    end
end 

function decintime = dit(ditsignal)
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
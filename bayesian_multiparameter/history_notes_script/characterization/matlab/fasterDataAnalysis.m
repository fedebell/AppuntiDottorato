Num1 = 30000;
numPhase = 17;

sizeSpace = 5;

sparseMSEMean = zeros(Num1, 1);
sparseNean = zeros(Num1, 1);
sparseMSEError = zeros(Num1, 1);
sparseNError = zeros(Num1, 1);

step = 2;

for N=1:step:100
    counter = 0;
    sommaError = 0;
    sommaN = 0;
    sommaError2 = 0;
    sommaN2 = 0;
    for p=0:1:step-1
        for k=1:1:numPhase
            if(results(k, N+p, sizeSpace+1)~=0)
                sommaError = sommaError + results(k, N+p, sizeSpace+1);
                sommaError2 = sommaError2 + results(k, N+p, sizeSpace+1)^2;
                sommaN = sommaN + N+p;
                sommaN2 = sommaN2 + (N+p)^2;
                counter = counter + 1;
            end
        end
    end
    sparseMSEMean(round(sommaN/counter)) = sommaError/counter;
    sparseMSEError(round(sommaN/counter)) = sqrt((sommaError2/counter)-(sommaError/counter)^2);
    sparseNError(round(sommaN/counter)) = sqrt((sommaN2/counter)-(sommaN/counter)^2);
end


step = 20;

for N=100:step:1000
    counter = 0;
    sommaError = 0;
    sommaN = 0;
    sommaError2 = 0;
    sommaN2 = 0;
    for p=0:1:step-1
        for k=1:1:numPhase
            if(results(k, N+p, sizeSpace+1)~=0)
                sommaError = sommaError + results(k, N+p, sizeSpace+1);
                sommaError2 = sommaError2 + results(k, N+p, sizeSpace+1)^2;
                sommaN = sommaN + N+p;
                sommaN2 = sommaN2 + (N+p)^2;
                counter = counter + 1;
            end
        end
    end
    sparseMSEMean(round(sommaN/counter)) = sommaError/counter;
    sparseMSEError(round(sommaN/counter)) = sqrt((sommaError2/counter)-(sommaError/counter)^2);
    sparseNError(round(sommaN/counter)) = sqrt((sommaN2/counter)-(sommaN/counter)^2);
end

step = 100;

for N=1000:step:Num1-step
    counter = 0;
    sommaError = 0;
    sommaN = 0;
    sommaError2 = 0;
    sommaN2 = 0;
    for p=0:1:step-1
        for k=1:1:numPhase
            if(results(k, N+p, sizeSpace+1)~=0)
                sommaError = sommaError + results(k, N+p, sizeSpace+1);
                sommaError2 = sommaError2 + results(k, N+p, sizeSpace+1)^2;
                sommaN = sommaN + N+p;
                sommaN2 = sommaN2 + (N+p)^2;
                counter = counter + 1;
            end
        end
    end
    sparseMSEMean(round(sommaN/counter)) = sommaError/counter;
    sparseMSEError(round(sommaN/counter)) = sqrt((sommaError2/counter)-(sommaError/counter)^2);
    sparseNError(round(sommaN/counter)) = sqrt((sommaN2/counter)-(sommaN/counter)^2);
end


%Some of the entries of RMSE and xaxis are useless, because not all the entries of the vector results were codified (we jump of 100 resources at a time for example)
%Count the the number of non-null positions
counterNonZero = 0;
for N=1:1:Num1
    if(~(sparseMSEMean(N) == 0))
        counterNonZero = counterNonZero + 1;
    end
end

%Array with the RMSE and teh total number of resources with only non-trivial entries
MSEMean = zeros(counterNonZero, 1);
resources = zeros(counterNonZero, 1);
MSEMeanError = zeros(counterNonZero, 1);
resourcesError = zeros(counterNonZero, 1);

%Loading of the significat entries of RMSE and xaxis in arrayRMSE and arrayResources

index = 1;

for i=1:1:counterNonZero
    
    while((sparseMSEMean(index) == 0))
        index = index + 1;
    end
    
    MSEMean(i) = sparseMSEMean(index);
    resources(i) = index;
    
    MSEMeanError(i) = sparseMSEError(index);
    resourcesError(i) = sparseNError(index);
    
    index = index + 1;
    
end

%min(1/(4*floor(resources(i)/1))*(4-3*0.900^2), 1/12+(0.5-0.900)^2);
%min(1/(4*floor(resources(i)/2))*(4-3*0.875^2), 1/12+(0.5-0.875)^2);
%min(1/(4*floor(resources(i)/11))*(4-3*0.850^2), 1/12+(0.5-0.850)^2);
%min(1/(4*floor(resources(i)/51))*(4-3*0.765^2), 1/12+(0.5-0.765)^2);


%Lower bound on the precision
hs = zeros(1, counterNonZero);
for i=1:1:counterNonZero
    hs(i) = (pi/resources(i))^2 + min(1/(4*floor(resources(i)/51))*(4-3*0.765^2), 1/12+(0.5-0.765)^2)+min(1/(4*floor(resources(i)/11))*(4-3*0.850^2), 1/12+(0.5-0.850)^2)+...
        min(1/(4*floor(resources(i)/2))*(4-3*0.875^2), 1/12+(0.5-0.875)^2)+min(1/(4*floor(resources(i)/1))*(4-3*0.900^2), 1/12+(0.5-0.900)^2);
end

%SQL and HS as references
loglog(resources, hs)
hold on;

%{
%Plot of the precision of the Bayesian algorithm
loglog(arrayResources, arrayRMSE)
hold on
%}

e1 = plot(resources, MSEMean);

%e2 = errorbar(resources, MSEMean, MSEMeanError, MSEMeanError, resourcesError, resourcesError, 'CapSize', 0);
e1.LineWidth = 2.0;
%set(gca, 'XScale','log', 'YScale','log')

set(gcf,'position',[200, 200, 600, 400]);

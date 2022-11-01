Num1 = 30000;
numPhase = 17;

visibilities = [0.900, 0.875, 0.850, 0.765];

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
    
    if(counter ~= 0)
    
    sparseMSEMean(round(sommaN/counter)) = sommaError/counter;
    sparseMSEError(round(sommaN/counter)) = sqrt((sommaError2/counter)-(sommaError/counter)^2);
    sparseNError(round(sommaN/counter)) = sqrt((sommaN2/counter)-(sommaN/counter)^2);
    
    end
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
    
    if(counter ~= 0)
    
    sparseMSEMean(round(sommaN/counter)) = sommaError/counter;
    sparseMSEError(round(sommaN/counter)) = sqrt((sommaError2/counter)-(sommaError/counter)^2);
    sparseNError(round(sommaN/counter)) = sqrt((sommaN2/counter)-(sommaN/counter)^2);
    
    end
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
    
    if(counter ~= 0)
    
    sparseMSEMean(round(sommaN/counter)) = sommaError/counter;
    sparseMSEError(round(sommaN/counter)) = sqrt((sommaError2/counter)-(sommaError/counter)^2);
    sparseNError(round(sommaN/counter)) = sqrt((sommaN2/counter)-(sommaN/counter)^2);
    
    end
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

s = 1;

%Lower bound on the precision
hs = zeros(1, counterNonZero);
for i=1:1:counterNonZero
    hs(i) = (pi/resources(i))^2 + s/(4*resources(i));
end


%SQL and HS as references
loglog(resources, hs)
hold on;

%{
%Plot of the precision of the Bayesian algorithm
loglog(arrayResources, arrayRMSE)
hold on
%}

%e = plot(arrayResources, arrayRMSEerror);

e = errorbar(resources, MSEMean, MSEMeanError, MSEMeanError, resourcesError, resourcesError, 'CapSize', 0);
e.LineWidth = 2.0;
%set(gca, 'XScale','log', 'YScale','log')

set(gcf,'position',[200, 200, 600, 400]);


hold on;
e.LineWidth = 2.0;

set(gcf,'position',[200, 200, 600, 400]);
set(gca, 'XScale','log', 'YScale','log')

%Q-plates charges
s = [1 2 11 51];

%Computation of the lower bound
Itheta = zeros(1, 4);
Ivis = zeros(1, 4);

for i=1:1:4
    Itheta(i) = 2*s(i)^2*(1-sqrt(1-(visibilities(i))^2));
    Ivis(i) = 2*(1-sqrt(1-(visibilities(i))^2))/((visibilities(i))^2*sqrt(1-(visibilities(i))^2));
end

G = Q;

cvx_begin
variables nu(4)
minimize G(1, 1)*inv_pos(Itheta(1)*nu(1)+Itheta(2)*nu(2)+Itheta(3)*nu(3)+Itheta(4)*nu(4))+...
    G(2, 2)*inv_pos(Ivis(1)*nu(1))+G(3, 3)*inv_pos(Ivis(2)*nu(2))+G(4, 4)*inv_pos(Ivis(3)*nu(3))+G(5, 5)*inv_pos(Ivis(4)*nu(4));
subject to
nu >= 0;
2*(s(1)*nu(1)+s(2)*nu(2)+s(3)*nu(3)+s(4)*nu(4)) == 1;
cvx_end

CRbound = G(1, 1)*inv_pos(Itheta(1)*nu(1)+Itheta(2)*nu(2)+Itheta(3)*nu(3)+Itheta(4)*nu(4))+...
    G(2, 2)*inv_pos(Ivis(1)*nu(1))+G(3, 3)*inv_pos(Ivis(2)*nu(2))+G(4, 4)*inv_pos(Ivis(3)*nu(3))+G(5, 5)*inv_pos(Ivis(4)*nu(4));

%SQL plot
sql = zeros(1, counterNonZero);
for i=1:1:counterNonZero
    sql(i) = CRbound/resources(i);
end

loglog(resources, sql)


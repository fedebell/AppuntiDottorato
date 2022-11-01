%Simulated RMSE with the iterative algorithm
fileSource = fopen("data/iterativeTrueViscorrectederror", "r");

Num1 = 30000;
Num2 = 100;

phaseNum = 1;

thetaArray = [0.00235
0.06145
0.29345
0.21635
0.38000
0.49620
0.46280
1.13975
1.26430
1.44890
1.66450
1.75520
1.87500
2.12710
2.58995
2.74000
2.96000];

visibilities = [0.900, 0.875, 0.850, 0.765];

%Vector containing the RMSE of Num2 (=30) and the resource of experiments performed with the same number of resources
RMSE1 = zeros(Num1, phaseNum);
RMSE2 = zeros(Num1, phaseNum);
xaxis = zeros(Num1, phaseNum);

RMSEerror1 = zeros(Num1, phaseNum);
RMSEerror2 = zeros(Num1, phaseNum);
xaxiserror = zeros(Num1, phaseNum);

for l=1:1:phaseNum
      
    step = 2;
    
    for N=1:step:100
        
        somma = 0;
        somma4 = 0;
        counter = 0;
        
        Var = 0;
        errorRMSE = 0;
        errorN = 0;
        
        meanRes = 0;
        meanRes2 = 0;
        
        listN = zeros(1000, 1);
        
        for p=0:1:step-1
            for k=1:1:Num2
                if(~(results(l, N + p, k, 1) == 0))
                    counter = counter + 1;
                    somma = somma + (angularDist(results(l, N+p, k, 1), thetaArray(l)))^2;
                    somma4 = somma4 + (angularDist(results(l, N+p, k, 1), thetaArray(l)))^4;
                    listN(counter) = N+p;
                end
            end
        end
                
        if(counter ~= 0)
            
            RMSE1(N, l) = sqrt(somma/counter);
            
            sigma = sqrt((somma4/counter) - (somma/counter)^2);
            
            errorRMSE = sigma/(2*sqrt(somma));

            Nmedio = mean(mean(listN(1:counter)));
            
            errorN = sqrt(1/counter*sum((listN(1:counter)-Nmedio*ones(counter, 1)).^2));

            RMSEerror1(N, l) = errorRMSE;
            
            xaxiserror(N, l) = errorN;
            
            %We take as number of resources used the average value of the number of resources used in the Num2 (=30) repetitions of the experiment
            xaxis(N, l) = Nmedio;
        end
        
    end
    
    step = 2;
    
    for N=1:step:100
        
        somma = 0;
        somma4 = 0;
        counter = 0;
        
        Var = 0;
        errorRMSE = 0;
        errorN = 0;
        
        meanRes = 0;
        meanRes2 = 0;
        
        listN = zeros(1000, 1);
        
        for p=0:1:step-1
            for k=1:1:Num2
                if(~(results(l, N + p, k, 2) == 0))
                    counter = counter + 1;
                    somma = somma + (results(l, N+p, k, 2) - visibilities(2))^2;
                    somma4 = somma4 + (angularDist(results(l, N+p, k, 2), thetaArray(l)))^4;
                    listN(counter) = N+p;
                end
            end
        end
                
        if(counter ~= 0)
            
            RMSE2(N, l) = sqrt(somma/counter);
            
            sigma = sqrt((somma4/counter) - (somma/counter)^2);
            
            errorRMSE = sigma/(2*sqrt(somma));

            Nmedio = mean(mean(listN(1:counter)));
            
            errorN = sqrt(1/counter*sum((listN(1:counter)-Nmedio*ones(counter, 1)).^2));

            RMSEerror2(N, l) = errorRMSE;
            
            xaxiserror(N, l) = errorN;
            
            %We take as number of resources used the average value of the number of resources used in the Num2 (=30) repetitions of the experiment
            xaxis(N, l) = Nmedio;
        end
        
    end
    
    
    step = 20;
    
    for N=100:step:1000
        
        somma = 0;
        somma4 = 0;
        counter = 0;
        
        Var = 0;
        errorRMSE = 0;
        errorN = 0;
        
        meanRes = 0;
        meanRes2 = 0;
        
        listN = zeros(1000, 1);
        
        for p=0:1:step-1
            for k=1:1:Num2
                if(~(results(l, N + p, k, 1) == 0))
                    counter = counter + 1;
                    somma = somma + (angularDist(results(l, N+p, k, 1), thetaArray(l)))^2;
                    somma4 = somma4 + (angularDist(results(l, N+p, k, 1), thetaArray(l)))^4;
                    listN(counter) = N+p;
                end
            end
        end
                
        if(counter ~= 0)
            
            RMSE1(N, l) = sqrt(somma/counter);
            
            sigma = sqrt((somma4/counter) - (somma/counter)^2);
            
            errorRMSE = sigma/(2*sqrt(somma));

            Nmedio = mean(mean(listN(1:counter)));
            
            errorN = sqrt(1/counter*sum((listN(1:counter)-Nmedio*ones(counter, 1)).^2));

            RMSEerror1(N, l) = errorRMSE;
            
            xaxiserror(N, l) = errorN;
            
            %We take as number of resources used the average value of the number of resources used in the Num2 (=30) repetitions of the experiment
            xaxis(N, l) = Nmedio;
            
        end
        
    end
    
    step = 20;
    
    for N=100:step:1000
        
        somma = 0;
        somma4 = 0;
        counter = 0;
        
        Var = 0;
        errorRMSE = 0;
        errorN = 0;
        
        meanRes = 0;
        meanRes2 = 0;
        
        listN = zeros(1000, 1);
        
        for p=0:1:step-1
            for k=1:1:Num2
                if(~(results(l, N + p, k, 1) == 0))
                    counter = counter + 1;
                    somma = somma + (results(l, N+p, k, 2) - visibilities(2))^2;
                    somma4 = somma4 + (angularDist(results(l, N+p, k, 1), thetaArray(l)))^4;
                    listN(counter) = N+p;
                end
            end
        end
                
        if(counter ~= 0)
            
            RMSE2(N, l) = sqrt(somma/counter);
            
            sigma = sqrt((somma4/counter) - (somma/counter)^2);
            
            errorRMSE = sigma/(2*sqrt(somma));

            Nmedio = mean(mean(listN(1:counter)));
            
            errorN = sqrt(1/counter*sum((listN(1:counter)-Nmedio*ones(counter, 1)).^2));

            RMSEerror2(N, l) = errorRMSE;
            
            xaxiserror(N, l) = errorN;
            
            %We take as number of resources used the average value of the number of resources used in the Num2 (=30) repetitions of the experiment
            xaxis(N, l) = Nmedio;
            
        end
        
    end
    
    step = 50;
    
    for N=1000:step:Num1-step
        
        somma = 0;
        somma4 = 0;
        counter = 0;
        
        Var = 0;
        errorRMSE = 0;
        errorN = 0;
        
        meanRes = 0;
        meanRes2 = 0;
        
        listN = zeros(1000, 1);
        
        for p=0:1:step-1
            for k=1:1:Num2
                if(~(results(l, N + p, k, 1) == 0))
                    counter = counter + 1;
                    somma = somma + (angularDist(results(l, N+p, k, 1), thetaArray(l)))^2;
                    somma4 = somma4 + (angularDist(results(l, N+p, k, 1), thetaArray(l)))^4;
                    listN(counter) = N+p;
                end
            end
        end
                
        if(counter ~= 0)
            
            RMSE1(N, l) = sqrt(somma/counter);
            
            sigma = sqrt((somma4/counter) - (somma/counter)^2);
            
            errorRMSE = sigma/(2*sqrt(somma));

            Nmedio = mean(mean(listN(1:counter)));
            
            errorN = sqrt(1/counter*sum((listN(1:counter)-Nmedio*ones(counter, 1)).^2));

            RMSEerror1(N, l) = errorRMSE;
            
            xaxiserror(N, l) = errorN;
            
            %We take as number of resources used the average value of the number of resources used in the Num2 (=30) repetitions of the experiment
            xaxis(N, l) = Nmedio;
        end
        
    end
    
    step = 50;
    
    for N=1000:step:Num1-step
        
        somma = 0;
        somma4 = 0;
        counter = 0;
        
        Var = 0;
        errorRMSE = 0;
        errorN = 0;
        
        meanRes = 0;
        meanRes2 = 0;
        
        listN = zeros(1000, 1);
        
        for p=0:1:step-1
            for k=1:1:Num2
                if(~(results(l, N + p, k, 1) == 0))
                    counter = counter + 1;
                    somma = somma + (results(l, N+p, k, 2) - visibilities(2))^2;
                    somma4 = somma4 + (angularDist(results(l, N+p, k, 1), thetaArray(l)))^4;
                    listN(counter) = N+p;
                end
            end
        end
                
        if(counter ~= 0)
            
            RMSE2(N, l) = sqrt(somma/counter);
            
            sigma = sqrt((somma4/counter) - (somma/counter)^2);
            
            errorRMSE = sigma/(2*sqrt(somma));

            Nmedio = mean(mean(listN(1:counter)));
            
            errorN = sqrt(1/counter*sum((listN(1:counter)-Nmedio*ones(counter, 1)).^2));

            RMSEerror2(N, l) = errorRMSE;
            
            xaxiserror(N, l) = errorN;
            
            %We take as number of resources used the average value of the number of resources used in the Num2 (=30) repetitions of the experiment
            xaxis(N, l) = Nmedio;
        end
        
    end
    
end

xaxisMean = zeros(Num1, 1);
RMSEMeans = zeros(Num1, 1);

xaxisErrorMean = zeros(Num1, 1);
RMSEErrorMean = zeros(Num1, 1);

for i=1:1:N
    
    RMSEMeans(i, 1) = mean(RMSE1(i, :)+RMSE2(i, :));
    xaxisMean(i, 1) = mean(xaxis(i, :));
    
    xaxisErrorMean(i, 1) = sqrt(sum(xaxiserror(i, :).^2))/phaseNum;
    RMSEErrorMean(i, 1) = sqrt(sum(RMSEerror1(i, :).^2))/phaseNum + sqrt(sum(RMSEerror2(i, :).^2))/phaseNum;
    
end

%Some of the entries of RMSE and xaxis are useless, because not all the entries of the vector results were codified (we jump of 100 resources at a time for example)
%Count the the number of non-null positions
counterNonZero = 0;
for N=1:1:Num1
    if(~(RMSEMeans(N) == 0))
        counterNonZero = counterNonZero + 1;
    end
end

step = 2;

%Array with the RMSE and teh total number of resources with only non-trivial entries
arrayRMSE = zeros(counterNonZero, 1);
arrayResources = zeros(counterNonZero, 1);
arrayRMSEerror = zeros(counterNonZero, 1);
arrayResourcesError = zeros(counterNonZero, 1);

%Loading of the significat entries of RMSE and xaxis in arrayRMSE and arrayResources

index = 1;

for i=1:1:counterNonZero
    
    while((RMSEMeans(index) == 0))
        index = index + 1;
    end
    
    arrayRMSEerror(i) = RMSEErrorMean(index);
    arrayResourcesError(i) = xaxisErrorMean(index);
    
    arrayRMSE(i) = RMSEMeans(index);
    arrayResources(i) = xaxisMean(index);
    
    index = index + 1;
    
end

%HS plot
hs = zeros(1, counterNonZero);
for i=1:1:counterNonZero
    hs(i) = pi/arrayResources(i);
end

%SQL plot
sql = zeros(1, counterNonZero);
for i=1:1:counterNonZero
    sql(i) = 1/sqrt(arrayResources(i));
end

%Loading of the simulated data with the iterative estimation algorithm (refering to the mean value visibilities) from 'midVis'
resNum = [];
RMSE = [];
Max = [];
Min = [];

%Loading of the data from the file 'MidVis'
while(~feof(fileSource))
    str = convertCharsToStrings(fgetl(fileSource));
    tokenResources = extractBetween(str, "", ": [");
    str = str + "@";
    tokenRMSE = extractBetween(str, "RMSE: ", ", ");
    tokenMax = extractBetween(str, "+Sigma: ", ", ");
    tokenMin = extractBetween(str, "-Sigma: ", "@");
    resNum = [resNum, str2num(tokenResources)];
    RMSE = [RMSE, str2num(tokenRMSE)];
    Max = [Max, str2num(tokenMax)];
    Min = [Min, str2num(tokenMin)];
end

%SQL and HS as references
loglog(arrayResources, hs)
hold on;
loglog(arrayResources, sql)
hold on;

%Plot of the iterative algorithm with error bars
p = scatter(resNum, RMSE, '.');
hold on;
e = errorbar(resNum, RMSE, RMSE - Min, Max - RMSE, 'CapSize', 0);
e.Color = '#7F00FF';
e.LineWidth = 1.0;
hold on;

%{
%Plot of the precision of the Bayesian algorithm
loglog(arrayResources, arrayRMSE)
hold on
%}

e2 = errorbar(arrayResources, arrayRMSE, arrayRMSEerror, arrayRMSEerror, arrayResourcesError, arrayResourcesError, 'CapSize', 0);
e2.LineWidth = 2.0;
%set(gca, 'XScale','log', 'YScale','log')

set(gcf,'position',[200, 200, 600, 400]);

%Function that computes the angular distance 
function dist = angularDist(ang1, ang2)
    dist = pi - abs(mod((ang1-ang2), 2*pi)-pi);
end

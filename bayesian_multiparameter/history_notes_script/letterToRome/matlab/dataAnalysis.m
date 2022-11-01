%Simulated RMSE with the iterative algorithm
fileSource = fopen("midVis", "r");

Num1 = 30000;
Num2 = 30;
theta = 2;

%Vector containing the RMSE of Num2 (=30) and the resource of experiments performed with the same number of resources
RMSE = zeros(Num1, 1);
xaxis = zeros(Num1, 1);

for N=1:1:Num1
    RMSE(N) = sqrt(sum((results(N, :, 1)-theta).^2)/Num2);
    %We take as number of resources used the average value of the number of resources used in the Num2 (=30) repetitions of the experiment
    xaxis(N) = mean(results(N, :, 6));
end

%Some of the entries of RMSE and xaxis are useless, because not all the entries of the vector results were codified (we jump of 100 resources at a time for example)
%Count the the number of non-null positions
counterNonZero = 0;
for N=1:1:Num1
    if(~(results(N, :, 1) == 0))
        counterNonZero = counterNonZero + 1;
    end
end

%Array with the RMSE and teh total number of resources with only non-trivial entries
arrayRMSE = zeros(counterNonZero, 1);
arrayResources = zeros(counterNonZero, 1);

%Loading of the significat entries of RMSE and xaxis in arrayRMSE and arrayResources
index = 1;
for i=1:1:counterNonZero
    while(results(index, :, 1) == 0)
        index = index + 1;
    end
    arrayRMSE(i) = RMSE(index);
    arrayResources(i) = xaxis(index);
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

%Plot of the precision of the Bayesian algorithm
loglog(arrayResources, arrayRMSE)
hold on;

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


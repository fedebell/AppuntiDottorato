%0 -> only theta
%1, 2, 3, 4 -> theta + a visibility
%5 -> all the parameters

numVis = 3;

%Maximum number of used resources
Num1 = 30000;

%Number of repetition of each experiment for a given phase
Num2 = 50;

%Q-plates charges
s = [1 2 11 51];

%Number of phases
phaseNum = 17;

%Loading of the true values of the parameters
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

trueValues = zeros(17, 5);

trueValues(:, 1) = thetaArray;

for i=1:1:4
    trueValues(:, i+1) = visibilities(i)*ones(17, 1);
end

%Vector containing the RMSE of Num2 and the resource of experiments performed with the same number of resources
RMSE = zeros(Num1, phaseNum);
xaxis = zeros(Num1, phaseNum);

RMSEerror = zeros(Num1, phaseNum);
xaxiserror = zeros(Num1, phaseNum);

for param=1:5
    
    for l=1:1:phaseNum
        
        step = 2;
        
        N = 1;
        
        while(N <= Num1-step)
            
            if((1 <= N) && (N <= 100))
                step = 2;
            elseif((100 <= N) && (N <= 1000))
                step = 20;
            else
                step = 50;
            end
            
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
                    if(~(results(l, N + p, k, param) == 0))
                        counter = counter + 1;
                        if(param == 1)
                            somma = somma + (angularDist(results(l, N+p, k, param), trueValues(l, param)))^2;
                            somma4 = somma4 + (angularDist(results(l, N+p, k, param), trueValues(l, param)))^4;
                        else
                            somma = somma + (results(l, N+p, k, param) - trueValues(l, param))^2;
                            somma4 = somma4 + (results(l, N+p, k, param) - trueValues(l, param))^2;
                        end
                        listN(counter) = N+p;
                    end
                end
            end
            
            if(counter ~= 0)
                
                RMSE(N, l) = sqrt(somma/counter);
                
                sigma = sqrt((somma4/counter) - (somma/counter)^2);
                
                errorRMSE = sigma/(2*sqrt(somma));
                
                Nmedio = mean(mean(listN(1:counter)));
                
                errorN = sqrt(1/counter*sum((listN(1:counter)-Nmedio*ones(counter, 1)).^2));
                
                RMSEerror(N, l) = errorRMSE;
                
                xaxiserror(N, l) = errorN;
                
                %We take as number of resources used the average value of the number of resources used in the Num2 (=30) repetitions of the experiment
                xaxis(N, l) = Nmedio;
            end
            
            N = N + step;
        end
        
    end
    
    xaxisMean = zeros(Num1, 1);
    RMSEMeans = zeros(Num1, 1);
    
    xaxisErrorMean = zeros(Num1, 1);
    RMSEErrorMean = zeros(Num1, 1);
    
    for i=1:1:N
        
        RMSEMeans(i, 1) = mean(RMSE(i, :));
        xaxisMean(i, 1) = mean(xaxis(i, :));
        
        xaxisErrorMean(i, 1) = sqrt(sum(xaxiserror(i, :).^2))/phaseNum;
        RMSEErrorMean(i, 1) = sqrt(sum(RMSEerror(i, :).^2))/phaseNum;
        
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
    
    if(param == 1)
        %Array with the RMSE and the total number of resources with only non-trivial entries
        arrayRMSE = zeros(counterNonZero, 5);
        arrayResources = zeros(counterNonZero, 5);
        arrayRMSEerror = zeros(counterNonZero, 5);
        arrayResourcesError = zeros(counterNonZero, 5);
    end
    
    %Loading of the significat entries of RMSE and xaxis in arrayRMSE and arrayResources
    
    index = 1;
    
    for i=1:1:counterNonZero
        
        while((RMSEMeans(index) == 0))
            index = index + 1;
        end
        
        arrayRMSEerror(i, param) = RMSEErrorMean(index);
        arrayResourcesError(i, param) = xaxisErrorMean(index);
        
        arrayRMSE(i, param) = RMSEMeans(index);
        arrayResources(i, param) = xaxisMean(index);
        
        index = index + 1;
        
    end
    
end


%ACHTUNG: Now the vector arrayRMSE containe the RMSE of the phase and
%visibilities, but we need the MSE (and its error), therefore we need to
%transform them

for i=1:1:counterNonZero
    for j=1:1:5
        arrayRMSEerror(i, j) = 2*arrayRMSE(i, j)*arrayRMSEerror(i, j);
        arrayRMSE(i, j) = (arrayRMSE(i, j))^2;
    end
end

G = zeros(5, 5);

if(numVis == 0)
    e = errorbar(arrayResources(:, 1), arrayRMSE(:, 1), arrayRMSEerror(:, 1), arrayRMSEerror(:, 1), arrayResourcesError(:, 1), arrayResourcesError(:, 1), 'CapSize', 0);
    G(1, 1) = 1;
elseif((numVis) > 0 && (numVis < 5))
    e = errorbar(arrayResources(:, 1)+arrayResources(:, numVis), arrayRMSE(:, 1)+arrayRMSE(:, numVis), arrayRMSEerror(:, 1)+...
        arrayRMSEerror(:, numVis), arrayRMSEerror(:, 1)+arrayRMSEerror(:, numVis), arrayResourcesError(:, 1)+arrayResourcesError(:, numVis), ...
        arrayResourcesError(:, 1)+arrayResourcesError(:, numVis), 'CapSize', 0);
    G(1, 1) = 1;
    G(numVis+1, numVis+1) = 1;
else
    e = errorbar(arrayResources(:, 1)+arrayResources(:, 2)+arrayResources(:, 3)+arrayResources(:, 4)+arrayResources(:, 5), ...
        arrayRMSE(:, 1)+arrayRMSE(:, 2)+arrayRMSE(:, 3)+arrayRMSE(:, 4)+arrayRMSE(:, 5), ...
        arrayRMSEerror(:, 1)+arrayRMSEerror(:, 2)+arrayRMSEerror(:, 3)+arrayRMSEerror(:, 4)+arrayRMSEerror(:, 5),  ...
        arrayRMSEerror(:, 1)+arrayRMSEerror(:, 2)+arrayRMSEerror(:, 3)+arrayRMSEerror(:, 4)+arrayRMSEerror(:, 5), ...
        arrayResourcesError(:, 1)+arrayResourcesError(:, 2)+arrayResourcesError(:, 3)+arrayResourcesError(:, 4)+arrayResourcesError(:, 5), ...
        arrayResourcesError(:, 1)+arrayResourcesError(:, 2)+arrayResourcesError(:, 3)+arrayResourcesError(:, 4)+arrayResourcesError(:, 5), ...
        'CapSize', 0);
    G = eye(5);
end

hold on;
e.LineWidth = 2.0;

set(gcf,'position',[200, 200, 600, 400]);
set(gca, 'XScale','log', 'YScale','log')

%Computation of the lower bound
Itheta = zeros(1, 4);
Ivis = zeros(1, 4);

for i=1:1:4
    Itheta(i) = 2*s(i)^2*(1-sqrt(1-(visibilities(i))^2));
    Ivis(i) = 2*(1-sqrt(1-(visibilities(i))^2))/((visibilities(i))^2*sqrt(1-(visibilities(i))^2));
end

CRbound = zeros(counterNonZero, 1);

for i=1:1:counterNonZero
    cvx_begin
    variables nu(4)
    minimize G(1, 1)*inv_pos(Itheta(1)*nu(1)+Itheta(2)*nu(2)+Itheta(3)*nu(3)+Itheta(4)*nu(4))+...
        G(2, 2)*inv_pos(Ivis(1)*nu(1))+G(3, 3)*inv_pos(Ivis(2)*nu(2))+G(4, 4)*inv_pos(Ivis(3)*nu(3))+G(5, 5)*inv_pos(Ivis(4)*nu(4));
    subject to
    nu >= 0;
    2*(s(1)*nu(1)+s(2)*nu(2)+s(3)*nu(3)+s(4)*nu(4)) == arrayResources(i, 1);
    cvx_end

    CRbound(i) = G(1, 1)*inv_pos(Itheta(1)*nu(1)+Itheta(2)*nu(2)+Itheta(3)*nu(3)+Itheta(4)*nu(4))+...
    G(2, 2)*inv_pos(Ivis(1)*nu(1))+G(3, 3)*inv_pos(Ivis(2)*nu(2))+G(4, 4)*inv_pos(Ivis(3)*nu(3))+G(5, 5)*inv_pos(Ivis(4)*nu(4));
end

loglog(arrayResources, CRbound);

%Function that computes the angular distance 
function dist = angularDist(ang1, ang2)
    dist = pi - abs(mod((ang1-ang2), 2*pi)-pi);
end
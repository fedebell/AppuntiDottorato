%This is the final algorithm for the simulation of the experiment in Rome.

phaseNum = 17;

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

stat = 30;

upp = 2000;

red = 1;

error = ones(1, 200);

error(1) = 0.900;
error(2) = 0.875;
error(11) = 0.850;
error(51) = 0.765;

alt = 0;

fileSource = fopen("sfile", "r");
if(fileSource == -1)
    disp("Achtung: file sfile not found!");
else
    while(~feof(fileSource))
        
        lines = convertCharsToStrings(fgetl(fileSource));
        token = extractBetween(lines,"(", ")");
        name = "(" + token + ")";
        token = extractBetween(lines,"[", "]");
        totalM = str2num(token); %#ok<ST2NM>

        resourceDistribution1 = phaseEstimation(name + "_1", 100, totalM, phaseNum, thetaArray, stat, 2, 2, upp, red, estrazioni, error, 0, alt);
        resourceDistribution2 = phaseEstimation(name + "_2", 1000, totalM, phaseNum, thetaArray, stat, 4, 100, upp, red, estrazioni, error, 0, alt);
        resourceDistribution3 = phaseEstimation(name + "_3", 10000, totalM, phaseNum, thetaArray, stat, 20, 1000, upp, red, estrazioni, error, 0, alt);
        resourceDistribution4 = phaseEstimation(name + "_4", 30000, totalM, phaseNum, thetaArray, stat, 1000, 10000, upp, red, estrazioni, error, 0, alt);

        %{
        %Simplified executions
        resourceDistribution1 = phaseEstimation(name + "_1", 100, totalM, phaseNum, thetaArray, stat, 20, 2, upp, red, estrazioni, error, 0, alt);
        resourceDistribution2 = phaseEstimation(name + "_2", 1000, totalM, phaseNum, thetaArray, stat, 100, 100, upp, red, estrazioni, error, 0, alt);
        resourceDistribution3 = phaseEstimation(name + "_3", 10000, totalM, phaseNum, thetaArray, stat, 500, 1000, upp, red, estrazioni, error, 0, alt);
        resourceDistribution4 = phaseEstimation(name + "_4", 30000, totalM, phaseNum, thetaArray, stat, 2000, 10000, upp, red, estrazioni, error, 0, alt);
        %}
        
    end
    
end

function resourceDistribution = phaseEstimation(name, RES, totalM, phaseNum, thetaArray, stat, step, startPoint, upp, red, estrazioni, error, strategies, alt)

    %Nummber of point for which the algorithm is performed
    numPoints = floor((RES-startPoint)/step)+1;

    Ktot = size(totalM, 2);
    
    %Il primo valore è l'errore medio su tutti gli angoli, il secondo è
    %l'errore massimo e il terzo l'errore miimi al variare dell'angolo che
    %si misura.
    resourceDistribution = zeros(numPoints, Ktot, 50);
    
    var = zeros(Ktot, numPoints, 4);
    
    if(strategies == 0)
        strategies = optimalStrategies(totalM, error, alt);
    end
    
    %For every maximal entanglement size
    for c=1:1:Ktot
        %For every number of resources.
        for l=1:1:numPoints
            
            fin = -1;
            for i=1:1:Ktot
                if(strategies(c, i, 1) > 0)
                    fin = i;
                end
            end
            
            %We load the strategies.

            [nu, f] = divisionLim((l-1)*step+startPoint, strategies(c, 1:fin, 1), strategies(c, 1:fin, 2), red, error, alt);
 
            %Find the maximum value of M, notice that -1 als starting is ok.
            MAX = -1;
            for i=1:1:fin
                if (nu(i) > MAX)
                    MAX = nu(i);
                end
            end
            
            %The division functions return -1 in the first position when the procedure fails.
            %In this line we check for the success of the division procedure. We also check that 
            %M(i) is in the range for which we have simulated measurements results.
            if((nu(1) > 0) && (MAX <= upp))
                [valore, max, min] = algorithm(nu, estrazioni, phaseNum, thetaArray, stat, strategies(c, 1:fin, 1), strategies(c, 1:fin, 2));
                lennu = size(nu, 2);
                for p=1:1:lennu
                    resourceDistribution(l, c, p) = nu(p);
                end
                fsum = sqrt(sum(f));
            else
                valore = 0;
                max = 0;
                min = 0;
                fsum = 0;
            end
            var(c, l, 1) = valore;
            var(c, l, 2) = max;
            var(c, l, 3) = min;
            var(c, l, 4) = fsum;
        end
    end
    
    mkdir simulations
    mkdir simulations/plot;
    mkdir simulations/strategies;
    mkdir simulations/distributions;
    
    %Building of the strategies string
    
    fileID = fopen('./simulations/strategies/'+name, 'wt');
    
    for i=1:1:Ktot
        fprintf(fileID, "[");
        for j=1:1:Ktot-1
            fprintf(fileID, "%i ", strategies(i, j, 1));
        end
        fprintf(fileID, "%i", strategies(i, Ktot, 1));
        fprintf(fileID, "]; ");
        fprintf(fileID, "[");
        for j=1:1:Ktot-1
            fprintf(fileID, "%f ", strategies(i, j, 2));
        end
        fprintf(fileID, "%f", strategies(i, Ktot, 2));
        fprintf(fileID, "]\n");
    end

    %close(fileID);
    
    %Plot of the variance without redistribution.
    analysisStage(name, var, RES, Ktot, step, startPoint, strategies, red, error, alt);
end

%Performs a logarithmic plot of the precision of the algorithm obtained by
%chosing the best results for every possible startegy in the vector startegies
%We compare our results against the SQL and the reachable HS (pi/N).
function toPlot = analysisStage(name, var, RES, Ktot, step, startPoint, strategies, red, error, alt)
    
    %Nummber of point for which the algorithm is performed
    numPoints = floor((RES-startPoint)/step)+1;
    
    %SQL array
    sql = zeros(1, numPoints);
    for i=1:1:numPoints
        sql(i) = 1/sqrt((i-1)*step+startPoint);
    end
    
    %HS array
    hs = zeros(1, numPoints);
    for i=1:1:numPoints
        hs(i) = pi/((i-1)*step+startPoint);
    end
    
    %x axis used in the plot
    xaxis = zeros(1, numPoints);
    for i=1:1:numPoints
        xaxis(i) = step*(i-1)+startPoint;
    end
    
    fileID = fopen('./simulations/distributions/'+name, 'wt');
    
    toPlot = zeros(1, numPoints);

    for i=1:1:numPoints
        
        varMin = 9999;
        flag = 0;
        for k=1:1:Ktot
            if((var(k, i, 1) > 0) && (var(k, i, 1) < varMin))
                flag = k;
                varMin = var(k, i, 1);
            end
        end
        
        if(flag ~= 0)
            toPlot(i) = varMin;
            
            fin = -1;
            for j=1:1:Ktot
                if(strategies(flag, j, 1) > 0)
                    fin = j;
                end
            end
            
            [M, ~] = divisionLim((i-1)*step+startPoint, strategies(flag, 1:fin, 1), strategies(flag, 1:fin, 2), red, error, alt);
            
            fprintf(fileID, "%i: [", (i-1)*step+startPoint);
            for p=1:1:fin-1
                fprintf(fileID, "%i ", strategies(flag, p, 1));
            end
            fprintf(fileID, "%i], [", strategies(flag, fin, 1));
             for p=1:1:fin-1
                fprintf(fileID, "%i ", M(p));
            end
            fprintf(fileID, "%i], RMSE: ", M(fin));
            fprintf(fileID, "%f", toPlot(i));
            fprintf(fileID, ", +Sigma: ");
            fprintf(fileID, "%f", var(flag, i, 2));
            fprintf(fileID, ", -Sigma: ");
            fprintf(fileID, "%f", var(flag, i, 3));
            fprintf(fileID, "\n");
        else
            toPlot(i) = 0;
        end
    end
    
    %close(fileID);
    
    hFig = figure('Visible', 'off'); 
    
    scatter(xaxis, toPlot, '.');
    hold on;
    %We add the SQL and the HS
    plot(xaxis, sql);
    hold on;
    plot(xaxis, hs);

    set(gca,'xscale','log');
    set(gca,'yscale','log');
    set(gca,'TickLabelInterpreter','latex')

    xlabel("$N$", 'interpreter','latex');
    ylabel("$\overline{\Delta \hat{\theta}}$", 'interpreter','latex');
    set(hFig, 'Position',  [200, 200, 1000, 600]);
    
    %Set CreateFcn callback
    set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
   
    set(hFig,'Units','inches');
    screenposition = get(gcf,'Position');
    set(hFig, 'PaperPosition',[0 0 screenposition(3:4)], 'PaperSize',screenposition(3:4));

    % Save Fig file
    saveas(hFig, "./simulations/plot/" + name + '.pdf');
    
    close(hFig);
end

function [varFinaleM, maxError, minError] = algorithm(nu, estrazioni, phaseNum, thetaArray, stat, M, n) 

    %Table of partial RMSE
    varianceTable = zeros(2, phaseNum);
    estimators = zeros(1, 16);
    
    lenNu = size(nu, 2);
    
    maxI = -1;
    for i=1:1:lenNu
        if (nu(i) > 0)
            maxI = i;
        end
    end

    %Phase estimation algorithm
    for counter=1:1:phaseNum
    
        theta = thetaArray(counter);
    
        sum = 0;
        
        sum4 = 0;
    
        for k=1:1:stat
         
            currentTheta = 0;
    
            for i=1:1:maxI
                
                if(i==1)
                    %k0 = 0;
                    k0 = floor((currentTheta - pi/(1))/(2*pi/(M(i))));
                else
                    k0 = floor((currentTheta - pi/(M(i-1)*n(i-1)))/(2*pi/(M(i))));
                end
                
                vett =  [1, 2, 11, 51];

                pos = find(vett == M(i));

                measure1 = 0;
                for t=1:1:nu(i)
                    measure1 = measure1 + estrazioni(counter, pos, t, 1, k);
                end

                measure2 = 0;
                for t=1:1:nu(i)
                    measure2 = measure2 + estrazioni(counter, pos, t, 2, k);
                end
                
                x = 2*measure1/nu(i) - 1;
                y = 2*measure2/nu(i) - 1;
                
                estimators(i) = atan2(y, x);
        	 
                %\xi estimators, they are casted into the [0, 2*pi] interval
                if(estimators(i) < 0)
                    estimators(i) = estimators(i) + 2*pi;
                end
        
                estimators(i) = estimators(i)/M(i);
                
                estimator0 = (k0 - 1)*(2*pi)/(M(i))+estimators(i);
                estimator1 = k0*(2*pi)/(M(i))+estimators(i);
                estimator2 = (k0 + 1)*(2*pi)/(M(i))+estimators(i);
                          
                %The algorithm can be expressed almost without reference to the
                %number n by saying that we chose always the estimator
                %which is closer to the old one. We still use n in the
                %definition of k_0. But the role of n is just to compute
                %the error, indeed the algorithm can be completely
                %reformulated withput mentioning n, just by saying that we
                %take the nearest value. This is done in algorithmWrong2.
                
                if(i==1)
                    %currentTheta = estimator1;
                    if((estimator0 + pi/(M(i)*n(i)) > currentTheta - pi/(1)) && ...
                            ((estimator0 - pi/(M(i)*n(i)) < currentTheta + pi/(1))) )
                        currentTheta = estimator0;
                    elseif((estimator2 + pi/(M(i)*n(i)) > currentTheta - pi/(1)) && ...
                            ((estimator2 - pi/(M(i)*n(i)) < currentTheta + pi/(1))) )
                        currentTheta = estimator2;
                    else
                        currentTheta = estimator1;
                    end
                else
                    if((estimator0 + pi/(M(i)*n(i)) > currentTheta - pi/(M(i-1)*n(i-1))) && ...
                            ((estimator0 - pi/(M(i)*n(i)) < currentTheta + pi/(M(i-1)*n(i-1)))) )
                        currentTheta = estimator0;
                    elseif((estimator2 + pi/(M(i)*n(i)) > currentTheta - pi/(M(i-1)*n(i-1))) && ...
                            ((estimator2 - pi/(M(i)*n(i)) < currentTheta + pi/(M(i-1)*n(i-1)))) )
                        currentTheta = estimator2;
                    else
                        currentTheta = estimator1;
                    end
                
                end
                    
                while (currentTheta  >= 2*pi)
                    currentTheta = currentTheta - 2*pi;
                end
            
                while( currentTheta < 0.0)
                    currentTheta = currentTheta + 2*pi;
                end  
                
            end
      
            %Angular distance
            error2 = (pi - abs(mod(currentTheta - theta, 2*pi)-pi))^2;
            error4 = (pi - abs(mod(currentTheta - theta, 2*pi)-pi))^4;

            sum = sum + error2;
            sum4 = sum4 + error4;
            
        end

        error = sqrt(sum/stat);
        
        sigma = sqrt(sum4/stat-(sum/stat)^2);
        errorOnError = (1/2)*sigma/sqrt(sum);
        
        varianceTable(1, counter) = error;
        varianceTable(2, counter) = errorOnError;
        
    end

    %Calcolo varianza media
    
    media = 0;
    for counter = 1:1:phaseNum
        media = media + varianceTable(1, counter);
    end
    media = media/phaseNum;

    %Inserisco la varianza dell'errore.
    
    sum = 0;
    for counter = 1:1:phaseNum
        sum = sum + varianceTable(2, counter)^2;
    end
    
    sum = sqrt(sum)/phaseNum;
    
    maxError = media + sum;
    minError = media - sum;
    
    %Calculation of the mean variance across all the angles.
    
    varFinaleM = media;
end

%Risommazione delle risorse.
function N = resummation(nu, M, alt)
    N=0;
    for j=1:1:size(nu, 2)
        if(alt == 1)
            N = N + 2*nu(j);
        else
            N = N + 2*nu(j)*M(j);
        end
    end
end

%Resource distribution according to the regular startegy.
function [array, f] = divisionLimAlt(N, M, n, red, error)
    
    b = 0.727;

    K = size(M, 2);
    
    D = zeros(1, K);

    D(1) = 1/2;
    
    for i=2:1:K
        for j=i:1:K-2
            D(i) = D(i) + 1/(M(j)*n(j));
        end
        D(i) = D(i) + 1/(M(K-1)*n(K-1)*2);
        D(i) = D(i) + 1/(M(K)*2);
        D(i) = D(i)*M(i-1)*n(i-1);
        D(i) = 1 + D(i);
    end
    
    %Aggiungo correzzione a questi D in modo che l'errore non possa mai eccedere pi.
    for j=2:1:K
        if(((2*pi*D(j))/(n(j-1)*M(j-1))) > pi)
            D(j) = n(j-1)*M(j-1)/(2);
        end
    end
    
    %Strategy arrays (raw and rounded), to be filled.
    array = zeros(1, K);
    arrayn = zeros(1, K);
    
    %Calcolo \Omega e \alpha
    
    if(K>1)
        omega = 2*(log(D(1)^2)+log(log(C(n(1), M(1), error))))/(log(C(n(1), M(1), error)));
        for j=2:1:K-1
            omega = omega + 2*(log(D(j)^2)+log(log(C(n(j), M(j), error))) - ...
            2*log((n(j-1))) - 2*log(M(j-1)))/(log(C(n(j), M(j), error)));
        end
    else
        omega = 0;
    end
    
    den = 0;
    for j=1:1:K-1
        den = den + 1/(log(C(n(j), M(j), error)));
    end
    
    syms a;
    syms x;
    [alpha, nuK] = vpasolve([2*x + omega - 2*a*den == N, -pi^2/(4*b*error(M(K))^2*x^2*(M(K))^2)-b*(error(M(K)))^2*((3*pi^2/4)/(M(K))^2)*exp(-b*(error(M(K)))^2*x)+((2*pi)^2)*exp(a) == 0], [a x]);
    
    %Loading of the strategy
    arrayn(1) = (log(D(1)^2)+log(log(C(n(1), M(1), error))) - alpha)/(log(C(n(1), M(1), error)));
    array(1) = floor(arrayn(1));
    for j=2:1:K-1
        arrayn(j) = (log(D(j)^2)+log(log(C(n(j), M(j), error))) - ...
            2*log((n(j-1))) - 2*log(M(j-1)) - alpha)/(log(C(n(j), M(j), error)));
        array(j) = floor(arrayn(j));
    end
    
    arrayn(K) = nuK;
    array(K) = floor(arrayn(K));
    
    if(red == 1)
        
        %Total resources effectively used.
        T = resummation(array, M, 1);
        
        %Remaining resources to be redistributed
        Delta = N - T;
        
        %A single probe cannot be used.
        if(mod(Delta, 2) == 1)
            if(Delta > 0)
                Delta = Delta - 1;
            end
        end
        
        f = zeros(1, K);
        
        array(K) = array(K) + Delta/2;
        
        f(1) = (2*pi)^2*(D(1)^2)*C(n(1), M(1), error)^(-array(1));
        for i=2:1:K-1
            f(i) = (2*pi*D(i))^2*C(n(i), M(i), error)^(-array(i))/(n(i-1)*M(i-1))^2;
        end
        f(K) = pi^2/(4*b*(error(M(K))^2)*array(K)*(M(K))^2)+(3*pi^2/4)/(M(K))^2*exp(-b*(error(M(K)))^2*array(K));

    end
    
    %If at the end there is not even a single probe in the final state (K-th) then it means that such
    flag = 0;
    for ind=1:1:K
        if(array(ind) < 1)
            flag = 1;
        end
    end
    if(flag == 1)
        array(1) = -1;
    end
end

%Resource distribution according to the regular startegy.
function [array, f] = divisionLimTrad(N, M, n, red, error)
    
    b = 0.727;

    K = size(M, 2);
    
    D = zeros(1, K);
    %D(1) = 1;
    D(1) = 1/2;
    for i=2:1:K
        for j=i:1:K-2
            D(i) = D(i) + 1/(M(j)*n(j));
        end
        D(i) = D(i) + 1/(M(K-1)*n(K-1)*2);
        D(i) = D(i) + 1/(M(K)*2);
        D(i) = D(i)*M(i-1)*n(i-1);
        D(i) = 1 + D(i);
    end
    
    %Aggiungo correzzione a questi D in modo che l'errore non possa mai eccedere pi.
    
    for j=2:1:K
        if(((2*pi*D(j))/(n(j-1)*M(j-1))) > pi)
            D(j) = n(j-1)*M(j-1)/(2);
        end
    end
    
    %Strategy arrays (raw and rounded), to be filled.
    
    array = zeros(1, K);
    arrayn = zeros(1, K);
    
    %Calcolo \Omega e \alpha
    
    if(K>1)
        omega = 2*(log(D(1)^2)+log(log(C(n(1), M(1), error))))/(log(C(n(1), M(1), error)));
        for j=2:1:K-1
            omega = omega + 2*(log(D(j)^2)+log(log(C(n(j), M(j), error))) - ...
            2*log((n(j-1))) - 2*log(M(j-1)) - log(M(j)))/(log(C(n(j), M(j), error)));
        end
    else
        omega = 0;
    end
    
    den = 0;
    for j=1:1:K-1
        den = den + M(j)/(log(C(n(j), M(j), error)));
    end
    
    syms a;
    syms x;
    [alpha, nuK] = vpasolve([2*x*M(K) + omega - 2*a*den == N, -pi^2/(4*b*error(M(K))^2*x^2*(M(K))^2)-b*(error(M(K)))^2*((3*pi^2/4)/(M(K))^2)*exp(-b*(error(M(K)))^2*x)+((2*pi)^2)*exp(a)*M(K) == 0], [a x]);
    
    %Loading of the strategy
    arrayn(1) = (log(D(1)^2)+log(log(C(n(1), M(1), error))) - alpha)/(log(C(n(1), M(1), error)));
    array(1) = floor(arrayn(1));
    for j=2:1:K-1
        arrayn(j) = (log(D(j)^2)+log(log(C(n(j), M(j), error))) - ...
            2*log((n(j-1))) - 2*log(M(j-1)) - log(M(j)) - alpha)/(log(C(n(j), M(j), error)));
        array(j) = floor(arrayn(j));
    end
    
    arrayn(K) = nuK;
    array(K) = floor(arrayn(K));
    
    if(red == 1)
        
        %Total resources effectively used.
        T = resummation(array, M, 0);
        
        %Remaining resources to be redistributed
        Delta = N - T;
        
        %A single probe cannot be used.
        if(mod(Delta, 2) == 1)
            if(Delta > 0)
                Delta = Delta -1;
            end
        end
        
        f = zeros(1, K);
        
        f(1) = (2*pi)^2*(D(1)^2)*C(n(1), M(1), error)^(-array(1));
        for i=2:1:K-1
            f(i) = (2*pi*D(i))^2*C(n(i), M(i), error)^(-array(i))/(n(i-1)*M(i-1))^2;
        end
        f(K) = pi^2/(4*b*(error(M(K))^2)*array(K)*(M(K))^2)+(3*pi^2/4)/(M(K))^2*exp(-b*(error(M(K)))^2*array(K));

        [~, I] = sort(f, 'descend');
        
        while(Delta > 0)
            c = 1;
            while(Delta < 2*M(I(c)))
                c = c + 1;
            end
            Delta = Delta - 2*M(I(c));
            array(I(c)) = array(I(c)) + 1;
            if(I(c) == K)
                f(K) = pi^2/(4*b*(error(M(K))^2)*array(K)*(M(K))^2)+(3*pi^2/4)/(M(K))^2*exp(-b*(error(M(K)))^2*array(K));
            else
                f(I(c)) = f(I(c))/(C(n(I(c)), M(I(c)), error));
            end
            [~, I] = sort(f, 'descend');
        end
    end
    
    %If at the end there is not even a single probe in the final state (K-th) then it means that such
    flag = 0;
    for ind=1:1:K
        if(array(ind) < 1)
            flag = 1;
        end
    end
    if(flag == 1)
        array(1) = -1;
    end
end

function [array, f] = divisionLim(N, M, n, red, error, alt)
    if(alt == 1)
        [array, f] = divisionLimAlt(N, M, n, red, error);
    else
        [array, f] = divisionLimTrad(N, M, n, red, error);
    end
end

function strategies = optimalStrategies(totalM, error, alt)

    Ktot = size(totalM, 2);
    
    strategies = zeros(Ktot, Ktot, 2);
    
    for K=1:1:Ktot
    
        MIN = 9999;
        Mminindex = 0;
        n1ott = 1;
        
        for i=0:1:2^(K-2)-1
            
            M = decode(i, totalM(1, 1:K));
            
            [a, b] = compn1(M);
            extr1 = 1/b;
            extr2 = 1/a;
            
            for counter=1:1:39
                
                n1= extr1 + (extr2 - extr1)/40*counter;
                n = listn(n1, M);
                
                f = figureMerit(n, M, error, alt);
                
                if(f < MIN)
                    MIN = f;
                    Mminindex = i;
                    n1ott = n1;
                end
            end
        end
        
        Mmin = decode(Mminindex, totalM(1, 1:K));
        n = listn(n1ott, Mmin);
        
        for i=1:1:size(Mmin, 2)
            strategies(K, i, 1) = Mmin(i);
            strategies(K, i, 2) = n(i);
        end
    end

end

function [a, b] = compn1(M)

    K = size(M, 2);

    if(K == 1)
        a = 1;
        b = 1;
    elseif(K == 2)
        a = 1/(M(2)+1);
        b = 1/(M(2)+1);
    else
        %epsilon(1) = 0 is not used
        epsilon = zeros(1, K);
        for j=2:1:K
            epsilon(j) = M(j)/M(j-1);
        end
        
        a = 1/epsilon(K-1)*(1-1/epsilon(K));
        b = 1/epsilon(K-1);
        
        for j=K-2:-1:2
            a = 1/epsilon(j)*(1-a);
            b = 1/epsilon(j)*(1-b);
        end
        
        if(b < a)
            tmp = b;
            b = a;
            a = tmp;
        end
    end
end

function n = listn(n1, M)
    K = size(M, 2);
    
    %epsilon(1) = 0 is not used
    epsilon = zeros(1, K);
    for j=2:1:K
    	epsilon(j) = M(j)/M(j-1);
    end
    
    n = zeros(1, K);
    n(1) = n1;
    for j=2:1:K
        n(j) = n(j-1)/(n(j-1)-epsilon(j));
    end
end

function f = figureMerit(n, M, error, alt)

    K = size(M, 2);
    
    if(alt == 1)
        f = 1/(log(C(n(1), M(1), error)));
        for i=2:1:K-1
            f = f + 1/((n(i-1))^2*(log(C(n(i), M(i), error))));
        end
    else
        f = M(1)/(log(C(n(1), M(1), error)));
        for i=2:1:K-1
            f = f + M(i)/((n(i-1))^2*(log(C(n(i), M(i), error))));
        end
    end
    
end


function C = C(n, s, error)
    if(n>=2)
        C = exp(0.727*(sin(pi/n))^2*(error(s)^2));
    else
        C = exp(0.727*(sin(pi/2))^2*(error(s)^2));
    end
end

%{
General model:
     f(x) = exp(b*(sin(pi/x)^2))
Coefficients (with 95% confidence bounds):
       b =       0.727  (0.7243, 0.7297)

Goodness of fit:
  SSE: 0.04569
  R-square: 0.9972
  Adjusted R-square: 0.9972
  RMSE: 0.01326
%}

function M = decode(c, totalM)

    Ktot = size(totalM, 2);
  
    flagArray = zeros(1, Ktot);
    flagArray(1) = 1;
    flagArray(Ktot) = Ktot;

    if(Ktot > 1)
        counter = 2;
    else
        counter = 1;
    end
    for j=2:1:(Ktot-1)
        if(mod(c, 2)==1)
            flagArray(j) = j;
            counter = counter + 1;
            c = c - 1;
        end
        c = c/2;
    end
    
    M = zeros(1, counter);
    
    j=1;
    counter = 1;
    while(j<=Ktot)
        if(flagArray(j))
            M(counter) = totalM(j);
            counter = counter + 1;
        end
        j = j+1;
    end
            
end
%M = [1 3 7 11 24 51 141];
M = [1 2 11 51];
alt = 0;

K = size(M, 2);

%Random visibilities
error = rand(2000);

error(1) = 0.900;
error(2) = 0.875;
error(11) = 0.850;
error(51) = 0.765;

n = optimaln(M, error, alt);

theta = 1;

%Algorithm

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

num = 50;
start = 20;
step = 20;

resultsPhaseEstimation = zeros(30000/step, 50);

for N=start:step:30000
    
    for counter = 1:1:num
        
        
        arrayP = zeros(K, 1);
        b = 0.25;
        
        %inizializzazione
        arrayP(1) = D(1)^2/(2);
        for i=2:1:K
            arrayP(i) = D(i)^2/(n(i-1)^2*M(i-1)^2*M(i));
        end
        
        estVis2 = ones(K, 1);
        dist = zeros(K, 1);
        results = zeros(K, 2);
        
        ris = 0;
        while (ris < N)
            
            [~, pos] = max(arrayP);
            dist(pos) = dist(pos)+1;
            
            arrayP(1) = D(1)^2/(2)*exp(-dist(1)*b*estVis2(1)*sin(pi/n(1))^2)*estVis2(1)*sin(pi/n(1))^2;
            for i=2:1:K
                arrayP(i) = D(i)^2/(n(i-1)^2*M(i-1)^2*M(i))*exp(-dist(i)*b*estVis2(i)*sin(pi/n(i))^2)*estVis2(i)*sin(pi/n(i))^2;
            end
            
            ris = ris + M(pos);
            
            %Simulazione
            p0 = (1+error(M(pos))*cos(M(pos)*theta))/2;
            
            if(rand(1)<p0)
                results(pos, 1) = results(pos, 1) + 1;
            end
            
            p1 = (1+error(M(pos))*sin(M(pos)*theta))/2;
            
            if(rand(1)<p1)
                results(pos, 2) = results(pos, 2) + 1;
            end
            
            %Stima delle visibilità
            
            if(mod(ris, 20) == 0)
                
                for j=1:1:K
                    if(dist(j) > 1)
                        estVis2(j) = min(max(0, (dist(j)*((2*(results(j, 1)/dist(j))-1)^2+(2*(results(j, 2)/dist(j))-1)^2)-1)/(dist(j)-1)), 1);
                        %estVis2(j) = (2*(results(j, 1)/dist(j))-1)^2+(2*(results(j, 2)/dist(j))-1)^2;
                    end
                end
                
            end
            
        end
        
        resultsPhaseEstimation(N/step, counter) = algorithm(dist, results, 1, 1, M, n);
        
    end
    
end


function nlopt = optimaln(M, error, alt)

[a, b] = compn1(M);
extr1 = 1/b;
extr2 = 1/a;

div = 100;

MIN = 9999;
for counter=1:1:(div-1)
    
    n1 = extr1 + (extr2 - extr1)/div*counter;
    n = listn(n1, M);
    
    f = figureMerit(n, M, error, alt);
    
    if(f < MIN)
        MIN = f;
        n1ott = n1;
    end
end

nlopt = listn(n1ott, M);

end

    
function f = figureMerit(n, M, error, alt)

    K = size(M, 2);
    
    if(alt == 1)
        f = 1/(log(C(n(1), M(1), error)));
        for i=2:1:K
            f = f + 1/((n(i-1))^2*(log(C(n(i), M(i), error))));
        end
    else
        f = M(1)/(log(C(n(1), M(1), error)));
        for i=2:1:K
            f = f + M(i)/((n(i-1))^2*(log(C(n(i), M(i), error))));
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



function C = C(n, s, error)
    b = 0.25;
    if(n>=2)
        C = exp(b*(sin(pi/n))^2*(error(s)^2));
    else
        C = exp(b*(sin(pi/2))^2*(error(s)^2));
    end
end


function currentTheta = algorithm(nu, results, iterations, stat, M, n) 


    estimators = zeros(1, 16);
    lenNu = size(nu, 2);
    
    maxI = -1;
    for i=1:1:lenNu
        if (nu(i) > 0)
            maxI = i;
        end
    end

    %Phase estimation algorithm
    for counter=1:1:iterations
    
        for k=1:1:stat
         
            currentTheta = 0;
    
            for i=1:1:maxI
                
                if(i==1)
                    %k0 = 0;
                    k0 = floor((currentTheta - pi/(1))/(2*pi/(M(i))));
                else
                    k0 = floor((currentTheta - pi/(M(i-1)*n(i-1)))/(2*pi/(M(i))));
                end
                
                vett  = [1, 2, 11, 51];

                pos = find(vett == M(i));

                x = 2*results(pos, 1)/nu(i) - 1;
                y = 2*results(pos, 2)/nu(i) - 1;
                
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
        
        end
      
    end
    
end
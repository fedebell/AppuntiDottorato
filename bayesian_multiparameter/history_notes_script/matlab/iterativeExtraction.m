%This file generates the raw noisy data for the simulations, and it allows
%an error dependent on the stage.

sizeE  = [1, 2, 11, 51];

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

error = zeros(1, 200);

error(1) = 0.900;
error(2) = 0.875;
error(11) = 0.850;
error(51) = 0.765;

stat = 30;

startPoint = 2;

upp = 2000;

K = size(sizeE, 2);

rng('shuffle');

estrazioni = zeros(phaseNum, K, upp, 2, stat);

for counter=1:1:phaseNum
    
    theta = thetaArray(counter);
    
    for o=1:1:K
        
        errorTerm = error(sizeE(o));
        p1 = (1+errorTerm*cos(sizeE(o)*theta))/2;
        p2 = (1+errorTerm*sin(sizeE(o)*theta))/2;
        
        estrazioni(counter, o, :, 1, :) = binornd(1, p1*ones(upp, stat));
        estrazioni(counter, o, :, 2, :) = binornd(1, p2*ones(upp, stat));
        
    end
end

stringa = "extrazioniIterative.mat";

save(stringa, 'estrazioni', '-v7.3');
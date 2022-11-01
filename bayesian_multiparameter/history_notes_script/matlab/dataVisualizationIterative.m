hFig = figure('Visible', 'on');

fileSource1 = fopen("simulations/distributions/(1, 1, 1)_1", "r");
fileSource2 = fopen("simulations/distributions/(1, 1, 1)_2", "r");
fileSource3 = fopen("simulations/distributions/(1, 1, 1)_3", "r");
fileSource4 = fopen("simulations/distributions/(1, 1, 1)_4", "r");

fileID = fopen('iterative', 'wt');

if((fileSource1 == -1) || (fileSource2 == -1) || (fileSource3 == -1) || (fileSource4 == -1))
    disp("Achtung: files not found!");
else
    resNum1 = [];
    RMSE1 = [];
    Max1 = [];
    Min1 = [];
    while(~feof(fileSource1))
        str = convertCharsToStrings(fgetl(fileSource1));
        fprintf(fileID, str + '\n');
        tokenResources = extractBetween(str, "", ": [");
        str = str + "@";
        tokenRMSE = extractBetween(str, "RMSE: ", ", ");
        %tokenMax = extractBetween(str, "Max: ", ", ");
        %tokenMin = extractBetween(str, "Min: ", "@");
        tokenMax = extractBetween(str, "+Sigma: ", ", ");
        tokenMin = extractBetween(str, "-Sigma: ", "@");
        resNum1 = [resNum1, str2num(tokenResources)];
        RMSE1 = [RMSE1, str2num(tokenRMSE)];
        Max1 = [Max1, str2num(tokenMax)];
        Min1 = [Min1, str2num(tokenMin)];
    end
    
    resNum2 = [];
    RMSE2 = [];
    Max2 = [];
    Min2 = [];
    while(~feof(fileSource2))
        str = convertCharsToStrings(fgetl(fileSource2));
        fprintf(fileID, str + '\n');
        tokenResources = extractBetween(str, "", ": [");
        str = str + "@";
        tokenRMSE = extractBetween(str, "RMSE: ", ", ");
        %tokenMax = extractBetween(str, "Max: ", ", ");
        %tokenMin = extractBetween(str, "Min: ", "@");
        tokenMax = extractBetween(str, "+Sigma: ", ", ");
        tokenMin = extractBetween(str, "-Sigma: ", "@");
        resNum2 = [resNum2, str2num(tokenResources)];
        RMSE2 = [RMSE2, str2num(tokenRMSE)];
        Max2 = [Max2, str2num(tokenMax)];
        Min2 = [Min2, str2num(tokenMin)];
    end
    
    resNum3 = [];
    RMSE3 = [];
    Max3 = [];
    Min3 = [];
    while(~feof(fileSource3))
        str = convertCharsToStrings(fgetl(fileSource3));
        fprintf(fileID, str + '\n');
        tokenResources = extractBetween(str, "", ": [");
        str = str + "@";
        tokenRMSE = extractBetween(str, "RMSE: ", ", ");
        %tokenMax = extractBetween(str, "Max: ", ", ");
        %tokenMin = extractBetween(str, "Min: ", "@");
        tokenMax = extractBetween(str, "+Sigma: ", ", ");
        tokenMin = extractBetween(str, "-Sigma: ", "@");
        resNum3 = [resNum3, str2num(tokenResources)];
        RMSE3 = [RMSE3, str2num(tokenRMSE)];
        Max3 = [Max3, str2num(tokenMax)];
        Min3 = [Min3, str2num(tokenMin)];
    end
    
    resNum4 = [];
    RMSE4 = [];
    Max4 = [];
    Min4 = [];
    while(~feof(fileSource4))
        str = convertCharsToStrings(fgetl(fileSource4));
        fprintf(fileID, str + '\n');
        tokenResources = extractBetween(str, "", ": [");
        str = str + "@";
        tokenRMSE = extractBetween(str, "RMSE: ", ", ");
        %tokenMax = extractBetween(str, "Max: ", ", ");
        %tokenMin = extractBetween(str, "Min: ", "@");
        tokenMax = extractBetween(str, "+Sigma: ", ", ");
        tokenMin = extractBetween(str, "-Sigma: ", "@");
        resNum4 = [resNum4, str2num(tokenResources)];
        RMSE4 = [RMSE4, str2num(tokenRMSE)];
        Max4 = [Max4, str2num(tokenMax)];
        Min4 = [Min4, str2num(tokenMin)];
    end
    
    resNum = [resNum1, resNum2, resNum3, resNum4];
    %resNum = [resNum1, resNum2, resNum3];
    
    RMSE = [RMSE1, RMSE2, RMSE3, RMSE4];
    %RMSE = [RMSE1, RMSE2, RMSE3];
    
    Max = [Max1, Max2, Max3, Max4];
    Min = [Min1, Min2, Min3, Min4];
    
    resNumlog = log(resNum);
    RMSElog = log(RMSE);
    
    %Max = [Max1, Max2, Max3];
    %Min = [Min1, Min2, Min3];
    
    %From Max and Min we determine the size of the error bar, which is
    
    N = size(RMSE, 2);
    
    %SQL array
    sql = zeros(1, N);
    for i=1:1:N
        sql(i) = 1/sqrt(resNum(i));
    end
    
    %HS array
    hs = zeros(1, N);
    for i=1:1:N
        hs(i) = pi/(resNum(i));
    end
    
    scatter(resNum, RMSE, '.');
    hold on;
    
    %We add the SQL and the HS
    plot(resNum, sql);
    hold on;
    plot(resNum, hs);
    hold on;
    errorbar(resNum, RMSE, RMSE-Min, Max - RMSE, 'CapSize', 0);
    hold on;
    
    
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    set(gca,'TickLabelInterpreter','latex')
    title("Plot");
    xlabel("Resources");
    ylabel("$\Delta \hat{\theta}$", 'interpreter','latex');
    set(hFig, 'Position',  [200, 200, 1000, 600]);
    
    %Set CreateFcn callback
    set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    
    set(hFig,'Units','inches');
    screenposition = get(gcf,'Position');
    set(hFig, 'PaperPosition',[0 0 screenposition(3:4)], 'PaperSize',screenposition(3:4));
    
    %{
        resNumLog = zeros(N, 1);
        RMSELog = zeros(N, 1);
        
        for j=1:1:N
            resNumLog(j) = log(resNum(j));
        end
        
        for j=1:1:N
            RMSELog(j) = log(RMSE(j));
        end
        
        first = 0;
        for i=1:1:N
           if(RMSE(i) < sql(i))
              first = i;
              break;
           end
        end
        
        arrayC = zeros(N+1, 1);
        resNormArray = zeros(N+1, 1);
        
        for a=0:1:(N-first)
            
            F = @(k, x) k - x;
            k0 = 5;
            
            [k,resnorm,~,exitflag,output] = lsqcurvefit(F,k0,resNumLog(first:N-a),RMSELog(first:N-a));
            
            
            arrayC(a+1) = k;
            
            SStot = sum((RMSELog-mean(RMSELog)).^2);
            SSres = sum((RMSELog-F(k,resNumLog)).^2);
            Rsq = 1-SSres/SStot;
            
            resNormArray(a+1) = Rsq;
            
        end
        
        minr = 9999;
        mink = 9999;
        indexa = 9999;
        
        for a=0:1:N-first
           if(resNormArray(a+1) <= minr)
               mink = arrayC(a+1);
               minr = resNormArray(a+1);
               indexa = a;
           end
        end
       
        hold on
           plot(resNum, exp(F(mink,resNumLog)));
        hold off
    %}
    
end

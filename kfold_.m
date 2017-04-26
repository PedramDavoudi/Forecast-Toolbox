function [out_fis,outputf,Err,RMSE] = kfold_(Y,X,horizon)
% 1-fold cross validation, ANFIS-based cost function estimation

nf = size(X,2);% number of variables


% X=[X;nan(horizon,nf)];
% Y=[Y(:);nan(horizon,1)];
% X=lagmatrix(X,horizon);

EV=[X,Y];

% EV(1:horizon,:)=[];

nf=size(EV,2);
ns = size(EV,1); % Sample size

nepo= 150;  % number of learning epochs
x = 0.5*ones(1,nf);
%%%%%%%%%%%%%%%%%20 fold cross validation loop.%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:ns
    % ith sample is selected for test and others are choose for learning
    k=0;
    for j = i+1:ns
        k=k+1;
        EVt(k,:) = EV(j,:);
    end
    while k<ns-1
        k=k+1;
        EVt(k,:) = EV(ns - k,:);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ANFIS learning
    [ans,Err,out_fis] = main_ANFIS(nepo,x,EVt,EV);
    % Test ith sample
    output(i,1)=evalfis(EV(i,1:nf-1),out_fis);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Err = (output(1:20,:)-EV(1:20,nf));
RMSE = mse(Err)^.5;
outputf = output;
outputf=outputf(end-horizon+1:end);
end
function [output,Err,out_fis] = main_ANFIS(nepo,x,P_T,P_Tall)
% main function for ANFIS training 
trnData = P_T;
nf = size(P_T,2)-1;
EV = P_T(:,1:nf);

in_fis = genfis2(double(trnData(:,1:nf)),double(trnData(:,nf+1)),x,...
    [min(P_Tall,[],'omitnan'); max(P_Tall,[],'omitnan')]); 
    
out_fis = anfis(trnData,in_fis,nepo);
output=evalfis(trnData(:,1:nf),out_fis);

Err = mean(abs(output-trnData(:,nf).^2));
end



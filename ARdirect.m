function [yforecast,p]=ARdirect(AR_Var,dum,one_step,horizon)
% AR(p) Forecasting inflation bye direct approach 
%--------------------------------------------------------------------------
MaxLag=5; % Maximum Lag in Model Selection
IC=2;    %Information ceriterion 1=AIC, 2=SIC, 3=HQC
% i dont know why this is used. i must check later
% one_step=nan;

if sum(sum(dum(1:end-horizon,:)))==0
    d=0;% this variable check whether there is any dummy variables.
else
    d=1;
    s=sum(dum(1:end-horizon,:),1);
    m=1;
    for i=1:size(s,2)
        if s(1,i)~=0
            cn(1,m)=i;% this variable store the available dummies
            m=m+1;
        end
        
    end
end

%--------------------Preallocating Variables-------------------------------
%beta=cell(MaxLag,1);
% LagIn=zeros(MaxLag,1);% the index store the best lag by AIC SIC HQC, which is choosen below
yforecast=zeros(1,horizon);
%--------------------------------------------------------------------------



    beta=cell(MaxLag,1);
    LagIn=zeros(MaxLag,1);% the index store the best lag by AIC SIC HQC, which is choosen below

    for i=1:MaxLag;
        laginf=lagmatrix(AR_Var, 1:i);
        y=AR_Var(i+1:end);
        switch d
            case 0
                xx=[ones(length(y),1), laginf(i+1:end,:)];
            case 1
                xx=[ones(length(y),1), laginf(i+1:end,:), dum(i+1:i+length(y),cn)];
        end
        beta{i,1}=xx\y;
        yhat=xx*beta{i,1};
        e=y-yhat;
        sigma2=(e'*e)/(length(y)-i-1);
        switch IC
            case 1 % AIC
                LagIn(i,1)=-2*log(sigma2)+2*(i);
            case 2 %SIC
                LagIn(i,1)=log(sigma2)+(i)*log(length(y))/length(y);
            case 3 %HQC
                LagIn(i,1)=log(sigma2)+2*(i)*log(log(length(y)))/length(y);
        end
        
        clear y xx yhat e
    end
    
    [p, ra]=min(LagIn);
    %     [min_sic rs]=min(sic);
    %     [min_hqc rh]=min(hqc);
   
    %     laginff=lagmatrix(AR_Var, 0:ra-1);            % shif data for forecasting one-step ahead %
    Indx=length(AR_Var)-ra+1:length(AR_Var); % which rows would be used in calculation
    Indx=sort(Indx,'descend');
%     laginff=AR_Var(Indx);
for h=1:horizon
    yf=fliplr(yforecast(max(h-ra,1):h-1));
    if d==0
        yforecast(1,h)=[1 yf, AR_Var(Indx(1:end-h+1)).' ]*beta{ra,1};
    else
        yforecast(1,h)=[1, yf,AR_Var(Indx(1:end-h+1)).', dum(Indx(end)+h,cn)]*beta{ra,1};
    end
    
end 

clear laginff
if ~isnan(one_step)
    yforecast(1,1)=one_step;
end
end


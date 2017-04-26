function output = ANFIS_(Y,X,horizon)
%  adaptive neuro-fuzzy inference system
nf = size(X,2);% number of variables
ns = size(Y,1); % Sample size

X=[X;nan(horizon,nf)];
Y=[Y(:);nan(horizon,1)];
X=lagmatrix(X,horizon);

P_Tall=[X,Y];

P_Tall(1:horizon,:)=[];

nepo=150;% number of learning epochs
TraR=80; % Size of training sample in percent
x=0.5*ones(1,nf+1);%relative importance:radii = [0.6020, 0.6200, 0.5890, 0.4990];
nsTra = round(TraR*ns/100);

for i = 1:nsTra
    trnData(i,1:nf) = P_Tall(i,1:nf);
    trnData(i,nf+1) = P_Tall(i,nf+1);
end




mfType = 'gbellmf'; %'gauss2mf';%;'gbellmf''gaussmf'

in_fis = genfis2(double(trnData(:,1:nf)),double(trnData(:,nf+1)),x,...
[min(P_Tall,[],'omitnan'); max(P_Tall,[],'omitnan')]);
out_fis = anfis(trnData,in_fis,nepo);
output=evalfis(P_Tall(:,1:nf),out_fis);
output=output(end-horizon+1:end);
%{
close all;
figure(1);
plot(P_Tall(:,nf+1),'--r');
hold on
plot(output);
hold off;

Err = mean(abs(output-P_Tall(:,nf).^2));

figure; plotregression(P_Tall(:,nf+1),evalfis(P_Tall(:,1:nf),in_fis))
 
%}
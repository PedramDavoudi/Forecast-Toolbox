function DoAggregate(Weight,Names,Input_Adr)
T=length(Names);
Best=nan(T,5+4);
Min=Best;
Max=Best;

for i=1:T
    Data=dataset('xls',[Input_Adr, Names{i}, '.xlsx']); % ,'Range','C1:C2'
    Data = sortrows(Data,'RMSEF','ascend');
    Best(i,:)=double(Data(1,1:9))*Weight.(Names{i})/100;
    Min(i,:)=min(double(Data(1:20,1:9)))*Weight.(Names{i})/100;
    Max(i,:)=max(double(Data(1:20,1:9)))*Weight.(Names{i})/100;
end

Best=sum(Best,1);
Min=sum(Min,1);
Max=sum(Max,1);
Result_Database=mat2dataset([Best;Min;Max],'varnames',Data.Properties.VarNames(1:9),'obsnames',{'Best';'Min';'Max'});
export(Result_Database,'xlsfile',['Output\MainGroupAggr.xlsx']);
end
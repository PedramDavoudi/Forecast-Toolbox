function GetData2(includedate)
if ~exist('includedate','var')
   includedate=1; 
end
% data = Quandl.get('NSE/OIL');
%
% %mydata = Quandl.get('NSE/OIL', 'start_date','yyyy-mm-dd','end_date','yyyy-mm-dd');
% mydata = Quandl.get('NSE/OIL', 'start_date','2010-01-1','end_date','2015-01-04');
% mydata = Quandl.get('NSE/OIL', 'collapse' ,'annual');% ("weekly"|"monthly"|"quarterly"|"annual")
%
% % Transformations: `
% mydata = Quandl.get('NSE/OIL', 'transformation' ,'rdiff'); % ("diff"|"rdiff"|"normalize"|"cumulative")
%  %* Return only n number of rows: `
%  mydata = Quandl.get('NSE/OIL','rows',5);
%
% %% ## Available Data Types ##
% % There are four options for which datatype you would like your data returned as, you choose your type as follows:
% %
% % 	Quandl.get('NSE/OIL','type','ts')
% %
% % * **Timeseries (default)**: returns a timeseries if only 1 column in data, tscollection if more. `('type','ts')`
% % * **Financial timeseries** :`('type','fints')`
% % * **CSV string**: `('type','ASCII')`
% % * **DataMatrix**: `('type','data')`
% % * **Cell Strings**: `('type','cellstr')`
%
%  mydata = Quandl.get('NSE/OIL','rows',5,'type','fints');
%  %
%     data = Quandl.datatable('ZACKS/FE');
%
%     Quandl.search('crude oil');
%%
formatOut = 'yy.dd';
SD='2009-01-1';
ED=date;%'2017-12-31';
IndxBank=cell(100,1);
%comodities
% Agriculture
Ind=1;
Com(Ind)=Quandl.get('COM/WLD_IAGRICULTURE', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'Agriculture'};
Ind=Ind+1;
% Food
Com(Ind)=Quandl.get('COM/WLD_IFOOD', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'Food'};
Ind=Ind+1;
% Energy
Com(Ind)=Quandl.get('COM/WLD_IENERGY', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'Energy'};
Ind=Ind+1;
% Metals & Minerals Index
Com(Ind)=Quandl.get('COM/WLD_IMETMIN', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'MetalsMinerals'};
Ind=Ind+1;
% Raw Materials Index
Com(Ind)=Quandl.get('COM/WLD_IRAW_MATERIAL', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'RawMaterials'};
Ind=Ind+1;
% % Steel Index,(2005=100)
% Com(Ind)=Quandl.get('COM/WLD_ISTL_JP_INDX', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
% IndxBank(Ind)={'Steel'};
% Ind=Ind+1;
% Non-energy Index
Com(Ind)=Quandl.get('COM/WLD_INONFUEL', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'NonEnergy'};
Ind=Ind+1;
% Natural gas index,(2010=100)
Com(Ind)=Quandl.get('COM/WLD_INATGAS', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'NaturalGas'};
Ind=Ind+1;
% Base Metals (ex. iron ore) Index
Com(Ind)=Quandl.get('COM/WLD_IBASEMET', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'BaseMetals'};
Ind=Ind+1;
% All Commodity Price Index, 2005 = 100, includes both Fuel and Non-Fuel Price Indices
Com(Ind)=Quandl.get('COM/PALLFNF_INDEX', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'Commodity'};
Ind=Ind+1;
% Metals Price Index, 2005 = 100, includes Copper, Aluminum, Iron Ore, Tin, Nickel, Zinc, Lead, and Uranium Price Indices
Com(Ind)=Quandl.get('COM/PMETA_INDEX', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'Metals'};
Ind=Ind+1;

% Metals Price Index, 2005 = 100, includes Copper, Aluminum, Iron Ore, Tin, Nickel, Zinc, Lead, and Uranium Price Indices
Com(Ind)=Quandl.get('COM/PMETA_INDEX', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'Metals'};
Ind=Ind+1;

% Coal, Australian,($/mt)
Com(Ind)=Quandl.get('COM/WLD_COAL_AUS', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'Coal'};
Ind=Ind+1;
% Crude oil, WTI,($/bbl)
Com(Ind)=Quandl.get('COM/WLD_CRUDE_WTI', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'WTI'};
Ind=Ind+1;

% Crude oil, Brent,($/bbl)
Com(Ind)=Quandl.get('COM/WLD_CRUDE_BRENT', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'Brent'};
Ind=Ind+1;

% Copper,($/mt)
Com(Ind)=Quandl.get('COM/WLD_COPPER', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'Copper'};
Ind=Ind+1;

% Gold Price - London PM Fixing
Com(Ind)=Quandl.get('COM/AU_LPM', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'Gold'};
Ind=Ind+1;


% % Steel rebar,($/mt)
% Com(Ind)=Quandl.get('COM/WLD_STL_JP_REBAR', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
% IndxBank(Ind)={'rebar'};
% Ind=Ind+1;

% % Steel wire rod,($/mt)
% Com(Ind)=Quandl.get('COM/WLD_STL_JP_WIROD', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
% IndxBank(Ind)={'wirerod'};
% Ind=Ind+1;

% % Steel, hot rolled coilsheet,($/mt)
% Com(Ind)=Quandl.get('COM/WLD_STL_JP_HROLL', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
% IndxBank(Ind)={'hotRolled'};
% Ind=Ind+1;
% 
% % Steel, cold rolled coilsheet,($/mt)
% Com(Ind)=Quandl.get('COM/WLD_STL_JP_CROLL', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
% IndxBank(Ind)={' coldRolled'};
% Ind=Ind+1;


% Aluminium
Com(Ind)=Quandl.get('COM/WLD_ALUMINUM', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
IndxBank(Ind)={'Aluminium'};
Ind=Ind+1;

% % Iron Ore 62% Fe CFR China,CME
% Com(Ind)=Quandl.get('COM/FE_TJN', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
% IndxBank(Ind)={'IronOre'};
% Ind=Ind+1;

% Aluminium
% Com(Ind)=Quandl.get('COM/WLD_ALUMINUM', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
% IndxBank(Ind)={'Aluminium'};
% Ind=Ind+1;

% US CPI Market indexes


for i=1:Ind-1
    res=dataset(datenum(getabstime(Com(i))),Com(i).Data,'varnames',{'Date',IndxBank{i}});
    if i==1
        All=res;
    else
        All=join(All,res,'Type','outer', 'MergeKeys',true);
    end
end
if includedate==1
All.Date2=datestr(All.Date);
end
%{
% Opec Oil Price
oil=Quandl.get('OPEC/ORB', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
%oil.TimeInfo.StartDate;
oil=dataset(datenum(getabstime(oil)),oil.Data,'varnames',{'Date','oil'});
% Copper
CU=Quandl.get('LME/PR_CU', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
% CU.TimeInfo.getTimeStr;
CU=dataset(CU.CashBuyer.Data,datenum(getabstime(CU)),'varnames',{'CU','Date'});
% Aluminium
Al= Quandl.get('LME/PR_AA', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);

Al=dataset(Al.CashBuyer.Data,datenum(getabstime(Al)),'varnames',{'Al','Date'});
%  Al.TimeInfo.getTimeStr;
% Coal
Co=Quandl.get('ODA/PCOALAU_USD', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
Co=dataset(Co.Data,datenum(getabstime(Co)),'varnames',{'Co','Date'});
% Co=Quandl.get('EIA/COAL', 'collapse' ,'monthly', 'start_date','2012-01-1','end_date','2016-12-31');
%Co.TimeInfo.getTimeStr ;
% getabstime
% ORE Iron
Or=Quandl.get('ODA/PIORECR_USD', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);

Or=dataset(Or.Data,datenum(getabstime(Or)),'varnames',{'Or','Date'});

% US Energy price index
%Erg=Quandl.get('FRED/USACPIENGMINMEI', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
%Erg2=dataset(Erg.Data,datenum(getabstime(Erg),'varnames',{'Energy','Date'});

% Nasdaaq Index
NsQ=Quandl.get('NASDAQOMX/NQGI', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
NsQ= dataset(datenum(getabstime(NsQ)),NsQ.IndexValue.Data,'varnames',{'Date','NsQ'});


CpiUs= Quandl.get('BLSI/CUSR0000SA0', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
CpiUs=dataset(CpiUs.Data,datenum(getabstime(CpiUs)),'varnames',{'CPI','Date'});

% Commodity index
CMdi=Quandl.get('ODA/PALLFNF_INDEX', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
CMdi=dataset(CMdi.Data,datenum(getabstime(CMdi)),'varnames',{'Commodity','Date'});

clear d1

d1=join(NsQ,CpiUs,'Type','outer', 'MergeKeys',true);
d1=join(d1,CMdi,'Type','outer', 'MergeKeys',true);
d1=join(d1,oil,'Type','outer', 'MergeKeys',true);
d1 = join(d1,CU,'Type','outer', 'MergeKeys',true);
d1 = join(d1,Al,'Type','outer', 'MergeKeys',true);
d1 = join(d1,Co,'Type','outer', 'MergeKeys',true);
Data = join(d1,Or,'Type','outer', 'MergeKeys',true);
%}
%Data.DateIndex=datenum(Data.Date);
Data = sortrows(All,'Date','ascend');
delete('Input\Mar.xlsx')
export(Data,'xlsfile','Input\Mar.xlsx');
% Steel Bilet
warning('Input\Mar.xlsx created');
%Quandl.get('LME/PR_FM')

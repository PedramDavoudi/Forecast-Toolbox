function [Target, Exp_Var]=deseasonal(Target, Exp_Var,Date,Periodicity)
% [T1,N1]=size(Target);
% [T2,N2]=size(Exp_Var);

Data=[Target, Exp_Var];


% x12  Access to X12 seasonal adjustment program.
%
% Syntax with a single type of output requested
% ==============================================

%
% Syntax with mutliple types of output requested
% ===============================================
%
%     [Y1,Y2,...,OutpFile,ErrFile,Model,X] = x12(X,Range,...)
%
% See the option `'output='` for the types of output data available from
% X12.
%
% Input arguments
% ================
%
% * `X` [ tseries ] - Input data that will seasonally adjusted or filtered
% by the Census X12 Arima.
%
% * `Range` [ numeric ] - Date range on which the X12 will be run; if not
% specified or Inf the entire available range will be used.
%
% Output arguments
% =================
%
% * `Y`, `Y1`, `Y2`, ... [ tseries ] - Requested output data, by default
% only one type of output is returned, the seasonlly adjusted data; see the
% option `'output='`.
%
% * `OutpFile` [ cellstr ] - Contents of the output log files produced by
% X12; each cell contains the log file for one type of output requested.
%
% * `ErrFile` [ cellstr ] - Contents of the error files produced by X12;
% each cell contains the error file for one type of output requested.
%
% * `Model` [ struct ] - Struct array with model specifications and parameter
% estimates for each of the ARIMA models fitted; `Model` matches the size
% of `X` is 2nd and higher dimensions.
%
% * `X` [ tseries ] - Original input data with forecasts and/or backcasts
% appended if the options `'forecast='` and/or `'backcast='` are used.
%
% Options
% ========
%
% * `'backcast='` [ numeric | *`0`* ] - Run a backcast based on the fitted
% ARIMA model for this number of periods back to improve on the seasonal
% adjustment; see help on the `x11` specs in the X12-ARIMA manual. The
% backcast is included in the output argument `X`.
%
% * `'cleanup='` [ *`true`* | `false` ] - Delete temporary X12 files
% when done; the temporary files are named `iris_x12a.*`.
%
% * `'log='` [ `true` | *`false`* ] - Logarithmise the input data before,
% and de-logarithmise the output data back after, running `x12`.
%
% * `'forecast='` [ numeric | *`0`* ] - Run a forecast based on the fitted
% ARIMA model for this number of periods ahead to improve on the seasonal
% adjustment; see help on the `x11` specs in the X12-ARIMA manual. The
% forecast is included in the output argument `X`.
%
% * `'display='` [ `true` | *`false`* ] - Display X12 output messages in
% command window; if false the messages will be saved in a TXT file.
%
% * `'dummy='` [ tseries | *empty* ] - User dummy variable or variables (in
% case of a multivariate tseries object) used in X12-ARIMA regression; the
% dummy variables can also include values for forecasts and backcasts if
% you request them; the type of the dummy can be specified in the option
% `'dummyType='`.
%
% * `'dummyType='` [ `'ao'` | *`'holiday'`* | `'td'` ] - Type of the user
% dummy (which is specified through the option `'dummy='`); the three basic
% types of dummies are additive outlier (`'ao'`), holiday flows
% (`'holiday'`), and trading days (`'td'`); see the X12-ARIMA or X13-ARIMA
% documentation for more details (available from the U.S.Census Bureau
% website), look for the section on the REGRESSION spec, options 'user' and
% 'usertype'.
%
% * `'mode='` [ *`'auto'`* | `'add'` | `'logadd'` | `'mult'` | 
% `'pseudoadd'` | `'sign'` ] - Seasonal adjustment mode (see help on the
% `x11` specs in the X12-ARIMA manual); `'auto'` means that series with
% only positive or only negative numbers will be adjusted in the `'mult'`
% (multiplicative) mode, while series with combined positive and negative
% numbers in the `'add'` (additive) mode.
%
% * `'maxIter='` [ numeric | *`1500`* ] - Maximum number of iterations for
% the X12 estimation procedure. See help on the `estimation` specs in the
% X12-ARIMA manual.
%
% * `'maxOrder='` [ numeric | *`[2,1]`* ] - A 1-by-2 vector with maximum
% order for the regular ARMA model (can be `1`, `2`, `3`, or `4`) and
% maximum order for the seasonal ARMA model (can be `1` or `2`). See help
% on the `automdl` specs in the X12-ARIMA manual.
%
% * 'missing=' [ `true` | *`false`* ] - Allow for in-sample missing
% observations, and fill in values predicted by an estimated ARIMA process;
% if `false`, the seasonal adjustment will not run and a warning will be
% thrown.
%
% * `'output='` [ char | cellstr | *`'SA'`* ] - List of requested output
% data; the cellstr or comma-separated list can combine `'IR'` for the
% irregular component, `'SA'` for the final seasonally adjusted series,
% `'SF'` for seasonal factors, `'TC'` for the trend-cycle, and `'MV'` for
% the original series with missing observations replaced with ARIMA
% estimates. See also help on the `x11` specs in the X12-ARIMA manual.
%
% * `'saveAs='` [ char | *empty* ] - Name (or a whole path) under which
% X12-ARIMA output files will be saved.
%
% * `'specFile='` [ char | *`'default'`* ] - Name of the X12-ARIMA spec
% file; if `'default'` the IRIS default spec file will be used, see
% description.
%
% * `'tdays='` [ `true` | *`false`* ] - Correct for the number of trading
% days. See help on the `x11regression` specs in the X12-ARIMA manual.
%
% * `'tolerance='` [ numeric | *`1e-5`* ] - Convergence tolerance for the
% X12 estimation procedure. See help on the `estimation` specs in the
% X12-ARIMA manual.
%
% Description
% ============
%
% Missing observations
% ---------------------
%
% If you keep `'missing=' false` (this is the default for backward
% compatibility), `x12` will not run on series with in-sample missing
% observations, and a warning will be thrown.
%
% If you set `'missing=' true`, you allow for in-sample missing
% observations. The X12-ARIMA program handles missing observations by
% filling in values predicted by the estimated ARIMA process. You can
% request the series with missing values filled in by including `MV` in the
% option `'output='`.
%
% Spec file
% ----------
%
% The default X12-ARIMA spec file is `+thirdparty/x12/default.spc`. You can
% create your own spec file to include options that are not available
% through the IRIS interface. You can use the following pre-defined
% placeholders letting IRIS fill in some of the information needed (check
% out the default file):
%
% * `$series_data$` is replaced with a column vector of input observations;
% * `$series_freq$` is replaced with a number representing the date
% frequency: either 4 for quarterly, or 12 for monthly (other frequencies
% are currenlty not supported by X12-ARIMA);
% * `$series_startyear$` is replaced with the start year of the input
% series;
% * `$series_startper$` is replaced with the start quarter or month of the
% input series;
% * `$transform_function$` is replaced with `log` or `none` depending on
% the mode selected by the user;
% * `$forecast_maxlead$` is replaced with the requested number of ARIMA
% forecast periods used to extend the series before seasonal adjustment.
% * `$forecast_maxlead$` is replaced with the requested number of ARIMA
% forecast periods used to extend the series before seasonal adjustment.
% * `$tolerance$` is replaced with the requested convergence tolerance in
% the `estimation` spec.
% * `$maxiter$` is replaced with the requested maximum number of iterations
% in the `estimation` spec.
% * `$maxorder$` is replaced with two numbers separated by a blank space:
% maximum order of regular ARIMA, and maximum order of seasonal ARIMA.
% * `$x11_mode$` is replaced with the requested mode: `'add'` for additive,
% `'mult'` for multiplicative, `'pseudoadd'` for pseudo-additive, or
% `'logadd'` for log-additive;
% * `$x12_save$` is replaced with the list of the requested output
% series: `'d10'` for seasonals, `'d11'` for final seasonally adjusted
% series, `'d12'` for trend-cycle, `'d13'` for irregular component.
%
% Two of the placeholders, `'$series_data$` and `$x12_output$`, are
% required; if they are not found in the spec file, IRIS throws an error.
%
% Estimates of ARIMA model parameters
% ------------------------------------
%
% The ARIMA model specification, `Model`, is a struct with three fields:
%
% * `.spec` - a cell array with the first cell giving the structure of the
% non-seasonal ARIMA, and the second cell giving the
% structure of the seasonal ARIMA; both specifications follow the usual
% Box-Jenkins notation, e.g. `[0 1 1]`.
%
% * `.ar` - a numeric array with the point estimates of the AR coefficients
% (non-seasonal and seasonal).
%
% * `.ma` - a numeric array with the point estimates of the MA coefficients
% (non-seasonal and seasonal).
%
% Example 1
% ==========
%
% If you wish to run `x12` on the entire range on which the input time
% series is defined, and do not use any options, you can omit the second
% input argument (the date range). These following three calls to `x12` do
% exactly the same:
%
%     xsa = x12(x);
%     xsa = x12(x,Inf);
%     xsa = x12(x,get(x,'range'));
%
% Example 2
% ==========
%
% If you wish to specify some of the options, you have to enter a date
% range or use `Inf`:
%
%     xsa = x12(x,Inf,'mode=','add');
%

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

opt = passvalopt('tseries.x12',varargin{:});

if strcmp(opt.mode,'sign')
    opt.mode = 'auto';
end

if nargin > 1 && ~isnumeric(range)
    error('Incorrect type of input argument(s).');
end

if nargin < 2
    range = Inf;
end
range = specrange(x,range);

doOutput();

% Forecasts and backcasts.
if islogical(opt.arima)
    % Obsolete option `'arima='`.
    if opt.arima
        opt.backcast = 8;
        opt.forecast = 8;
    else
        opt.backcast = 0;
        opt.forecast = 0;
    end
elseif ~isempty(opt.arima)
    opt.backcast = opt.arima(1);
    opt.forecast = opt.arima(end);
end

%--------------------------------------------------------------------------

co = comment(x);
tmpsize = size(x.data);
x.data = x.data(:,:);
[data,range] = rangedata(x,range);

% Extended range with backcasts and forecasts.
xRange = range(1)-opt.backcast : range(end)+opt.forecast;

% Fill in zeros for NaNs in dummy variables on the extended range.
dummy = [];
if ~isempty(opt.dummy) && isa(opt.dummy,'tseries')
    dummy = rangedata(opt.dummy,xRange);
    dummy = dummy(:,:);
    dummy(isnan(dummy)) = 0;
end

if opt.log
    data = log(data);
end

outp = regexp(opt.output,'[a-zA-Z]\d\d','match');
nOutp = length(outp);

[data,varargout{1:nOutp},varargout{nOutp+(1:3)}] = ...
    thirdparty.x12.x12(data,range(1),dummy,opt);

for i = 1 : nOutp
    if opt.log
        varargout{i} = exp(varargout{i});
    end
    varargout{i} = ...
        reshape(varargout{i},[size(varargout{i},1),tmpsize(2:end)]);
    varargout{i} = replace(x,varargout{i},range(1),co);
end

% Reshape the model spec struct to match the dimensions and size of input
% and output tseries.
if length(tmpsize) > 2
    varargout{nOutp+3} = ...
        reshape(varargout{nOutp+3},[1,tmpsize(2:end)]);
end

% Return original series with forecasts and backcasts.
if nargout >= nOutp+3
    x.start = x.start - opt.backcast;
    x.data = data;
    if length(tmpsize) > 2
        x.data = reshape(x.data,[size(x.data,1),tmpsize(2:end)]);
    end
    x = mytrim(x);
    varargout{nOutp+4} = x;
end


% Nested functions.


%**************************************************************************
    function doOutput()
        % Map aliases for output arguments to X12 codes.
        opt.output = regexp(opt.output,'\w+','match');
        map = {...
            'SF','d10', ...
            'seasonals','d10', ...
            'seasonal','d10', ...
            'seasfactors','d10', ...
            'seasfact','d10', ...
            ...
            'SA','d11', ...
            'seasadj','d11', ...
            ...
            'TC','d12', ...
            'trendcycle','d12', ...
            'trend','d12', ...
            ...
            'IR','d13', ...
            'irregular','d13', ...
            'irreg','d13', ...
            ...
            'MV','mv', ...
            'missingvaladj','mv', ...
            'missing','mv', ...
            };
        for ii = 1 : 2 : length(map)
            opt.output = strrep(opt.output,map{ii},map{ii+1});
        end
    end % doOutput()


end
return
%% Two method are available now: 13 Movibng Average and Quadratic Parametric Trend Estimation  
method=1;

% Copyright 2015 The MathWorks, Inc.
if method==1
%% Seasonal Adjustment Using S(n,m) Seasonal Filters  
% This example shows how to apply $S_{n\times m}$ seasonal filters to deseasonalize
% a time series (using a multiplicative decomposition). The time series
% is monthly international airline passenger counts from 1949 to 1960.   
[T N]= size(Data);
% Copyright 2015 The MathWorks, Inc.
for n=1:N
y = Data(:,n);

%%
% The data shows an upward linear trend and a seasonal component with periodicity
% 12.  

%% Detrend the data using a 13-term moving average. 
% Before estimating the seasonal component, estimate and remove the linear
% trend. Apply a 13-term symmetric moving average, repeating the first and
% last observations six times to prevent data loss. Use weight 1/24 for
% the first and last terms in the moving average, and weight 1/12 for all
% interior terms. 
%
% Divide the original series by the smoothed series to detrend the data.
% Add the moving average trend estimate to the observed time series plot. 
%sW13 = [1/24;repmat(1/12,11,1);1/24];
sW13 = [1/(2*Periodicity);repmat(1/Periodicity,Periodicity-1,1);1/(2*Periodicity)];

yS = conv(y,sW13,'same');
%yS(1:6) = yS(7); yS(T-5:T) = yS(T-6);
yS(1:Periodicity/2) = yS(Periodicity/2+1); yS(T-(Periodicity/2-1):T) = yS(T-Periodicity/2);

xt = y./yS;
% 
% h = plot(yS,'r','LineWidth',2);
% legend(h,'13-Term Moving Average')
% hold off     

%% Create seasonal indices. 
% Create a cell array, |sidx|, to store the indices corresponding to each
% period. The data is monthly, with periodicity 12, so the first element
% of |sidx| is a vector with elements 1, 13, 25,...,133 (corresponding to
% January observations). The second element of |sidx| is a vector with elements
% 2, 14, 16,...,134 (corresponding to February observations). This is repeated
% for all 12 months. 
s = Periodicity;
sidx = cell(s,1); % Preallocation

for i = 1:s
 sidx{i,1} = i:s:T;
end

% sidx{1:2} 

%%
% Using a cell array to store the indices allows for the possibility that
% each period does not occur the same number of times within the span of
% the observed series.   

%% Apply an S(3,3) filter. 
% Apply a 5-term $S_{3\times 3}$ seasonal moving average to the
% detrended series |xt|. That is, apply a moving average to the January
% values (at indices 1, 13, 25,...,133), and then apply a moving average
% to the February series (at indices 2, 14, 26,...,134), and so on for the
% remaining months.  
%
% Use asymmetric weights at the ends of the moving average (using |conv2|).
% Put the smoothed values back into a single vector. 
%
% To center the seasonal component around one, estimate, and then divide
% by, a 13-term moving average of the estimated seasonal component.

% S3x3 seasonal filter
% Symmetric weights
%sW3 = [1/9;2/9;1/3;2/9;1/9];
sW3 = [1/9;2/9;1/3;2/9;1/9];
% Asymmetric weights for end of series
aW3 = [.259 .407;.37 .407;.259 .185;.111 0];

% Apply filter to each month
shat = NaN*y;
for i = 1:s
    ns = length(sidx{i});
    first = 1:4;
    last = ns - 3:ns;
    dat = xt(sidx{i});
    
    sd = conv(dat,sW3,'same');
    sd(1:2) = conv2(dat(first),1,rot90(aW3,2),'valid');
    sd(ns  -1:ns) = conv2(dat(last),1,aW3,'valid');
    shat(sidx{i}) = sd;
end

% 13-term moving average of filtered series
sW13 = [1/24;repmat(1/12,11,1);1/24];
sb = conv(shat,sW13,'same');
sb(1:6) = sb(s+1:s+6); 
sb(T-5:T) = sb(T-s-5:T-s);

% Center to get final estimate
s33 = shat./sb;

% figure
% plot(s33)
% h2 = gca;
% h2.XLim = [0,T];
% h2.XTick = 1:12:T;
% h2.XTickLabel = datestr(dates(1:12:T),10);
% title 'Estimated Seasonal Component'; 

%%
% Notice that the seasonal level changes over the range of the data. This
% illustrates the difference between an $S_{n\times m}$ seasonal filter and
% a stable seasonal filter. A stable seasonal filter assumes that the
% seasonal level is constant over the range of the data.

%% Apply a 13-term Henderson filter. 
% To get an improved estimate of the trend component, apply a 13-term Henderson
% filter to the seasonally adjusted series. The necessary symmetric and
% asymmetric weights are provided in the following code. 

% Deseasonalize series
dt = y./s33;
%{
% Henderson filter weights
sWH = [-0.019;-0.028;0;.066;.147;.214;
      .24;.214;.147;.066;0;-0.028;-0.019];
% Asymmetric weights for end of series
aWH = [-.034  -.017   .045   .148   .279   .421;
       -.005   .051   .130   .215   .292   .353;
        .061   .135   .201   .241   .254   .244;
        .144   .205   .230   .216   .174   .120;
        .211   .233   .208   .149   .080   .012;
        .238   .210   .144   .068   .002  -.058;
        .213   .146   .066   .003  -.039  -.092;
        .147   .066   .004  -.025  -.042  0    ;
        .066   .003  -.020  -.016  0      0    ;
        .001  -.022  -.008  0      0      0    ;
       -.026  -.011   0     0      0      0    ;
       -.016   0      0     0      0      0    ];

% Apply 13-term Henderson filter
first = 1:12;
last = T-11:T;
h13 = conv(dt,sWH,'same');
h13(T-5:end) = conv2(dt(last),1,aWH,'valid');
h13(1:6) = conv2(dt(first),1,rot90(aWH,2),'valid');

% New detrended series
xt = y./h13;

% figure
% plot(y)
% h3 = gca;
% h3.XLim = [0,T];
% h3.XTick = 1:12:T;
% h3.XTickLabel = datestr(dates(1:12:T),10);
% title 'Airline Passenger Counts';
% hold on
% plot(h13,'r','LineWidth',2);
% legend('13-Term Henderson Filter')
% hold off     

%% Apply an S(3,5) seasonal filter. 
% To get  6. an improved estimate of the seasonal component, apply a 7-term
% $S_{3\times 5}$ seasonal moving average to the newly detrended
% series. The symmetric and asymmetric weights are provided in the following
% code. Center the seasonal estimate to fluctuate around 1. 
%
% Deseasonalize the original series by dividing it by the centered seasonal
% estimate. 

% S3x5 seasonal filter 
% Symmetric weights
sW5 = [1/15;2/15;repmat(1/5,3,1);2/15;1/15];
% Asymmetric weights for end of series
aW5 = [.150 .250 .293;
       .217 .250 .283;
       .217 .250 .283;
       .217 .183 .150;
       .133 .067    0;
       .067   0     0];

% Apply filter to each month
shat = NaN*y;
for i = 1:s
    ns = length(sidx{i});
    first = 1:6;
    last = ns-5:ns;
    dat = xt(sidx{i});
    
    sd = conv(dat,sW5,'same');
    sd(1:3) = conv2(dat(first),1,rot90(aW5,2),'valid');
    sd(ns-2:ns) = conv2(dat(last),1,aW5,'valid');
    shat(sidx{i}) = sd;
end

% 13-term moving average of filtered series
sW13 = [1/24;repmat(1/12,11,1);1/24];
sb = conv(shat,sW13,'same');
sb(1:6) = sb(s+1:s+6); 
sb(T-5:T) = sb(T-s-5:T-s);


% Center to get final estimate
s35 = shat./sb;

% Deseasonalized series
dt = y./s35;

figure
plot(dt)
h4 = gca;
h4.XLim = [0,T];
h4.XTick = 1:12:T;
h4.XTickLabel = datestr(dates(1:12:T),10);
title 'Deseasonalized Airline Passenger Counts';    

%%
% The deseasonalized series consists of the long-term trend and irregular
% components. With the seasonal component removed, it is easier to see turning
% points in the trend.  

%% Plot the components and the original series. 
% Compare the original series to a series reconstructed using the component
% estimates. 
figure
plot(y,'Color',[.85,.85,.85],'LineWidth',4)
h5 = gca;
h5.XLim = [0,T];
h5.XTick = 1:12:T;
h5.XTickLabel = datestr(dates(1:12:T),10);
title 'Airline Passenger Counts';
hold on
plot(h13,'r','LineWidth',2)
plot(h13.*s35,'k--','LineWidth',1.5)
legend('Original Series','13-Term Henderson Filter',...
       'Trend and Seasonal Components')
hold off     

%% Estimate the irregular component. 
% Detrend and deseasonalize the original series. Plot the remaining estimate
% of the irregular component. 
Irr = dt./h13;

figure
plot(Irr)
h6 = gca;
h6.XLim = [0,T];
h6.XTick = 1:12:T;
h6.XTickLabel = datestr(dates(1:12:T),10);
title 'Airline Passenger Counts Irregular Component';

%%
% You can optionally model the detrended and deseasonalized series using
% a stationary stochastic process model.   
%}
Data(:,n)=dt;
end

elseif method==2
    %{
%% Step 1: Load the Data 
% Load the accidental deaths data set.
y = Data;
T = length(y);
 
%%
% The data shows a potential quadratic trend and a strong seasonal component
% with periodicity 12.  

%% Step 2: Fit Quadratic Trend
% Fit the polynomial
%
% $$T_t = \beta_0 + \beta_1t + \beta_2t^2$$
%
% to the observed series.  
t = (1:T)';
X = [ones(T,1) t t.^2];

b = X\y;
tH = X*b;
%  
% h2 = plot(tH/1000,'r','LineWidth',2);
% legend(h2,'Quadratic Trend Estimate')
% hold off
 
%% Step 3. Detrend Original Series. 
% Subtract the fitted quadratic line from the original data. 
xt = y - tH;  

%% Step 4. Estimate Seasonal Indicator Variables
% Create indicator (dummy) variables for each month. The first indicator
% is equal to one for January observations, and zero otherwise. The second
% indicator is equal to one for February observations, and zero otherwise.
% A total of 12 indicator variables are created for the 12 months. Regress
% the detrended series against the seasonal indicators. 
mo = repmat((1:12)',6,1);
sX = dummyvar(mo);
  
bS = sX\xt;
st = sX*bS;

figure
plot(st/1000)
title 'Parametric Estimate of Seasonal Component (Indicators)';
h3 = gca;
h3.XLim = [0,T];
ylabel 'Number of Deaths (in thousands)';
h3.XTick = 1:12:T;
h3.XTickLabel = datestr(dates(1:12:T),10);
%%
% In this regression, all 12 seasonal indicators are included in the design
% matrix. To prevent collinearity, an intercept term is not included (alternatively,
% you can include 11 indicators and an intercept term).  

%% Step 5. Deseasonalize Original Series
% Subtract the estimated seasonal component from the original series. 
dt = y - st;

% figure
% plot(dt/1000)
% title 'Monthly Accidental Deaths (Deseasonalized)';
% h4 = gca;
% h4.XLim = [0,T];
% ylabel 'Number of Deaths (in thousands)';
% h4.XTick = 1:12:T;
% h4.XTickLabel = datestr(dates(1:12:T),10);   
% %%
% The quadratic trend is much clearer with the seasonal component removed.  

%% Step 6. Estimate Irregular Component
% Subtract the trend and seasonal estimates from the original series. The
% remainder is an estimate of the irregular component. 
% bt = y - tH - st;

% figure
% plot(bt/1000)
% title('Irregular Component')
% h5 = gca;
% h5.XLim = [0,T];
% ylabel 'Number of Deaths (in thousands)';
% h5.XTick = 1:12:T;
% h5.XTickLabel = datestr(dates(1:12:T),10);    
    %}
end
%%
% You can optionally model the irregular component using a stochastic process
% model.
%%
% References: 
%
% Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. _Time Series Analysis: Forecasting and Control_. 3rd ed. Englewood Cliffs, NJ: Prentice Hall, 1994.  
Target=Data(:,1);
Exp_Var=Data(:,2:end);  
end
%{
q=Date(:);%qq(1900,1):qq(2000,1);
%q=q';


[T N]=size(data);
for i=1:N
    x1=tseries(q(1:T,1),data(1:T,i));
    x2=x12(x1);
    x3=x2.data;
    data_Adj(:,i)=x3;
end
Target=data_Adj(:,1:N1);
Exp_Var=data_Adj(:,N1+1:N1+N2);

end


%}
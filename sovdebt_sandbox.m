
% Jack Watson 
% Sovereign Debt project with Aayushi Mishra & Evan Kodra
% Initial sandbox for exploring data, building regression models, PCA, etc
dbstop if error

% root dir
rdir = "C:\Users\watson.jac\OneDrive - Northeastern University\";
% working dir
wdir = "C:\Users\watson.jac\Desktop\MATLAB scripts\Sovereign Debt";
% path to data directory
dpath = rdir + "\Research Projects\Sovereign Debt\Data\";
cd(wdir) % set working directory
addpath(rdir + "\Data\World map") % subfolder. for global mapping

%% Load data

% Load water stress data
% One scalar water stress value (WS.water_stress) per country (WS.country_name), for 1997 only (WS.year)
% Corresponds to x1 in Jack's notes
WS = readtable(dpath + "historic_water.csv");

% Load historic climate data
% One scalar per var per country, 1997 only
% vars disaggregated by extreme cold and heat and 100-year heatwaves
% Corresponds to x2, x3, x4, x5, x6, x7 in Jack's notes
HC = readtable(dpath + "historic_climate.csv");

% Load future climate data
% 81 values per variable per country corresponding to years 2020-2100 for each country
FC = readtable(dpath + "future_climate.csv");

% Load social impact score data
% 21 values of single var (SIS) for years 2000-2020 for each country
SI = readtable(dpath + "social_impact_score.csv");

% Load lagged financial data + climate scores
% 26 values per var for years 1995-2020, 7 financial vars, 1 PCRS, 1 is_default = 9 vars
% MATLAB puts an "x" in front of variable names that start with a number
LF = readtable(dpath + "lagged_full_data.csv");

%% Pre-process data
% Add 1997 water stress data to 1997 historical climate table
wsvec = zeros(height(HC),1);
for i = 1:length(wsvec)
    wsvec(i) = WS.water_stress(find(string(WS.iso3) == HC.iso3{i}));
end
HC.water_stress = wsvec;

% Create mapping of ISO3 country codes to full country names
ISO3NmT = table(WS.country_name, WS.iso3, VariableNames=["name", "iso3"]);


%% Scale/normalize data
% rescale defaults to scaling to [0,1] when given one argument
% NOTE: should probably scale these relative to min/max of indv variables,
% not to min/max of all vars in a given table. TODO later...

% use syntax: 
isscaleduqul = @(x) all(round(max(x,[],1)) == 1) && all(logical(sum(x==1, 1))) ...
    && all(round(min(x,[],1)) == 0); 
% (where x can be any of the below arrays) to check for correct normalization. 
% The "uqul" is "UniQue Upper Limit", for continuous variables like precise
% temperature and pressure. 
% The statement should return true if each variable (col) is scaled s.t.
% the max is 1, there is only one occurence of 1, and the min is zero. 
% Or just:
isscaled = @(x) all(round(max(x,[],1)) == 1) && all(round(min(x,[],1)) == 0);
% for discrete and discrete-derived variables.

ws = rescale(WS.water_stress); % just use rescale for 1D arrays/single vars
hc = normalize(table2array([HC(:,5:11)]), 1, "range", [0 1]);
fc = normalize(table2array([FC(:,5:11)]), 1, "range", [0 1]);
si = rescale(SI.social_impact_score);
lf = normalize(table2array([LF(:,4:10)]), 1, "range", [0 1]);

lfdvarnames = LF.Properties.VariableNames(4:11);
LFn = normalize(LF, "range", [0,1], "DataVariables", lfdvarnames);

hcdvarnames = HC.Properties.VariableNames(5:12);
HCn = normalize(HC, "range", [0,1], "DataVariables", hcdvarnames);


%% Logistic regression
% iso3ctry = "ARM";
% isdft = LFn.is_default(LFn.iso3 == iso3ctry);
% ctrl_fts = table2array(LFn(LFn.iso3 == iso3ctry, 4:11));

xtrain = LFn(LFn.year <= 2015,:);
xtest = LFn(LFn.year > 2015,:);

% logitmodelspec = "is_default ~ x5yr_avg_inflation*x5yr_avg_financial_dev" ...
%     + "*x5yr_avg_debt_to_gdp*x5yr_avg_gdp_deflator*x5yr_avg_gdp_growth" ...
%     + "*x5yr_avg_REER*x5yr_avg_terms_of_trade_index*x5yr_avg_climate_score" ...
%     + " - x5yr_avg_inflation:x5yr_avg_financial_dev:x5yr_avg_debt_to_gdp" ...
%     + ":x5yr_avg_gdp_deflator:x5yr_avg_gdp_growth:x5yr_avg_REER" ...
%     + ":x5yr_avg_terms_of_trade_index:x5yr_avg_climate_score";
% logitmodelspec = "is_default ~ x5yr_avg_inflation*x5yr_avg_financial_dev" ...
%     + "*x5yr_avg_debt_to_gdp*x5yr_avg_gdp_deflator*x5yr_avg_gdp_growth" ...
%     + "*x5yr_avg_REER*x5yr_avg_terms_of_trade_index*x5yr_avg_climate_score";

logitmodelspec = "is_default ~ x5yr_avg_inflation+x5yr_avg_financial_dev" ...
    + "+x5yr_avg_debt_to_gdp+x5yr_avg_gdp_deflator+x5yr_avg_gdp_growth" ...
    + "+x5yr_avg_REER+x5yr_avg_terms_of_trade_index+x5yr_avg_climate_score";

logitmdl = fitglm(xtrain, logitmodelspec, 'link','logit', 'Distribution', 'binomial');

% Evaluate logistic regression model on test set
ypredlogit = predict(logitmdl, xtest);
successtf = xtest.is_default == round(ypredlogit);
accy_logit = sum(successtf)/length(successtf);
% [xtest.is_default, round(ypredlogit)]


%% Principal component analysis (PCA)

%pca_X_HC = [table2array(HCn(:,6:8)), table2array(HCn(:,10:12))]; % without excess mortality
pca_X_HC = [table2array(HCn(:,6:8)), table2array(HCn(:,10:11))]; % without em or water stress
%pca_X_HC = table2array(HCn(:,5:12)); % with excess mortality

[pca_coeff_hc, pca_score_hc, pca_latent_hc] = pca(pca_X_HC, "NumComponents", 1);
pca_PCRS_HC = normalize(pca_score_hc, "range", [0,5]);


%% Visualize data

% Heatmaps ----------------------------------------------------------------
cmapHot = flipud(hot(256));
cmapRed = sky(256);
cmapRed = [cmapRed(:,3), cmapRed(:,1).*[1,1]];
% Water stress
%countryHeatmap(WS.water_stress, sNames, flipud(hot(256)))\
countryHeatmap(WS.water_stress, WS.country_name, cmapHot)

% PCA PCRS, HC ------------------------------------------------------------
HCsNames = iso32fullname(ISO3NmT, HCn.iso3);
LF97 = LF(LF.year == 1997,:);
LF97 = assign_by_iso3(HCn.iso3, pca_PCRS_HC, LF97, "PCRS_HC");
LF97.country_name = iso32fullname(ISO3NmT, LF97.iso3);
PCA_ICE_ratio = LF97.PCRS_HC./LF97.x5yr_avg_climate_score;
PCA_ICE_ratio(PCA_ICE_ratio > 38) = 38;
LF97_RMSE = rmse(LF97.PCRS_HC, LF97.x5yr_avg_climate_score);

countryHeatmap(pca_PCRS_HC, HCsNames, parula(256))
countryHeatmap(PCA_ICE_ratio, LF97.country_name, turbo(256))

% Histograms --------------------------------------------------------------

figure; 
histogram(PCA_ICE_ratio, 50); hold on

figure;
histogram(LF97.x5yr_avg_climate_score, 40); hold on
histogram(LF97.PCRS_HC, 40); hold on
legend({"ICE, Lagged full, 1997", "PCA, Historical Climate, 1997"})

figure;
h = plot(PCA_ICE_ratio, "o"); hold on; grid on
h.MarkerFaceColor = h.Color; % Change color from default to purple

% Time series -------------------------------------------------------------
sNamesTsFC = unique(FC.country_name);
tsFCbyState = cellfun(@(x) fc(string(FC.country_name) == string(x), :), ...
    sNamesTsFC, "UniformOutput",false);

hold on; arrayfun(@(x) plot(table2array(tsFCbyState{1}(x,:))), 5:11)
% or
hold on; plot(tsFCbyState{71})

% Purpose: calculate market risk economic capital (MEC) for the banking book using moving block bootstrapping 
% Input: historical yield curves and Greeks (key rate durations/Gammas)
% Output: stressed loss given time horizon and condidence interval 
% References: Rebonato et al., "Evolving yield curves in the real-world measures: a semi-parametric approach"  
% Written by: Enterprise Model Resources Group
% Version: Version 4
% Date: June 6, 2016 
% Update Date: October 10, 2017

clc; 
clear all;
close all;

% Control random number generation
rng('default'); 

% Load data from input file
filename = 'MEC_Input_2017Q2.xlsx';

[controls, val_date] = xlsread(filename, 'Control Page', 'B:B');
B = controls(1);                    % Number of Yield Curve Paths
jump_freq = controls(2);            % Jump Frequency
T = controls(3);                    % Time Horizon (months)for MEC
conf_interval = controls(4);        % Confidence Interval for MEC

% Historical data for bootstrapping of Interest Rates (%): 3m, 6m, 1y, 2Y, 3Y, 5Y, 7Y, 10Y, 15Y, 20Y and 30Y
[hist_data] = xlsread(filename, 'Historical Data', 'B2:L279'); % up to last term point as of valuation date 

% Yield Curve as of Valuation Date (e.g., Dec/1/2015)
[yc] = 0.01 * xlsread(filename, 'Initial Yield Curve', 'B:B');
yc = transpose(yc);

% Load Base Market Value (EVE) as of the valuation date 
mv0 = xlsread(filename, 'Base Market Value', 'B1');

% Load Key Rate Durations and Gamma from QRM (US Dollar per 1 BP)
[QRM] = xlsread(filename, 'Greeks','B3:B14');
[other_portfolios] = xlsread(filename, 'Greeks','E3:F4');

Dur_MSR_BP = other_portfolios(1,1);
Gamma_MSR_BP = other_portfolios(1,2);
Dur_BOLI_BP = other_portfolios(2,1);

KRD = (-10000)*QRM(1:11)';
Gamma = 100000000*QRM(12);
Dur_MSR = (-10000)*Dur_MSR_BP;
Gamma_MSR = 100000000*Gamma_MSR_BP;
Dur_BOLI = (-10000)*Dur_BOLI_BP;

[n, kz] = size(hist_data);

% Relative Change in Rates
ZZ = zeros(n-1,kz);

for j=1:n-1
    ZZ(j,:) = (hist_data(j+1,:) - hist_data(j,:))./hist_data(j,:); 
end

% Initialize Rate Shock Scenarios
rates = zeros(B,T,11);

% Loop over Bootstraps
for b = 1:B
    u = zeros(T,1);                 % Initialize Index

    % Number of Observations for Relative Change (ZZ): n-1 
    u(1) = ceil((n-1)*rand);
    
    for j = 1:11
        rates(b,1,j) = ZZ(u(1),j);
    end
    
    for t = 2:T  
        if rand <= jump_freq         % Skip to new series w/ prob jump_freq
           u(t) = ceil((n-1)*rand);
        elseif u(t-1) == n-1
           u(t) = 1;
        else
           u(t) = u(t-1) + 1;  
        end
        
        for j = 1:11
            rates(b,t,j) = ZZ(u(t),j);
        end
    end
end

% Cumulative change in rates [1+r(i)]
cumprods  = zeros(B,T,11);
for j = 1:11
    cumprods(:,:,j) = cumprod(1+rates(:,:,j),2);
end

cum_yc = zeros(B,11);
for j = 1:11
    cum_yc(:,j) = cumprods(:,end,j)-1;
end

yc_full = zeros(B,11,T);
for t = 1:T
    for j = 1:11
        yc_full(:,j,t) = cumprods(:,t,j)*yc(j);
    end
end
     
% EVE calculation for months 1 to 12
R_values = zeros(B,T);
mv_values = zeros(B,T);

yc_current = yc_full(:,:,1);

% Market Value Change at Day 1 End
for j = 1:11
    R_values(:,1) = R_values(:,1) - KRD(j)*(yc_current(:,j)-yc(:,j)); 
end
R_values(:,1) = R_values(:,1) -(Dur_MSR+Dur_BOLI)*(yc_current(:,6)-yc(:,6)) +...
    0.5*Gamma_MSR*(yc_current(:,6)-yc(:,6)).^2 + 0.5*Gamma*(yc_current(:,9)-yc(:,9)).^2;

mv_values(:,1) = mv0*(1+yc(1)/T)+R_values(:,1);

for t = 2:T %Market Value Change per Month
    yc_past = yc_current;
    yc_current = yc_full(:,:,t);
    
    for j = 1:11
        R_values(:,t) = R_values(:,t) - KRD(j)*(yc_current(:,j)-yc_past(:,j)); 
    end
    R_values(:,t) = R_values(:,t) -(Dur_MSR+Dur_BOLI)*(yc_current(:,6)-yc_past(:,6)) +...
        0.5*Gamma_MSR*(yc_current(:,6)-yc_past(:,6)).^2 + 0.5*Gamma*(yc_current(:,9)-yc_past(:,9)).^2;
    
    mv_values(:,t) = mv_values(:,t-1).*(1+yc_past(:,1)/T) + R_values(:,t);
end

% Change in EVE within 1-year horizon
delta_mv = mv_values-mv0;
PL = sort(delta_mv);

% Stressed loss at Selected Confidence Interval
loss_6m = PL(round(B*(1-conf_interval)),6);
mec_6m = max(0,-loss_6m);
mec_output = [loss_6m mec_6m];

% Stressed Yield Curve
k = find(mv_values(:,6)-mv0 == loss_6m); 
syc = yc_current(k,:);

filename1 = 'MEC Output with 6-Month Horizon_check';
xlswrite(filename1, mec_output);

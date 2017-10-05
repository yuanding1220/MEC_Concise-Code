% Purpose: Determine jump frequency in yield curve moving block bootstrapping
% References: Rebonato, Mahal, Joshi, Bucholz, and Nyholm (2005)
% Written by: Enterprise Model Resources Group

clc; 
clear all;
close all;

% Set random seed
rng('default')

% Load historical data for bootstrapping
filename = 'Historical Monthly Data_Aug 2017.xlsx';
Z = xlsread(filename,'C2:M281');
[n kz] = size(Z);

% Relative Change in Rates
ZZ = zeros(n-1,kz);
for j=1:n-1
    ZZ(j,:) = (Z(j+1,:) - Z(j,:))./Z(j,:); 
end

% Jump Frequency
JF = 1;

% Number of Simulated Interest Rate Paths/Yield Curve Scenarios
B = 1000000;

% Time Horizon on Monthly Basis (e.g., 12-month period)
T = 12;

% Initialize Rate Shock Scenarios
r3m = zeros(B,T);
r6m = zeros(B,T);
r12m = zeros(B,T);
r2y = zeros(B,T);
r3y = zeros(B,T);
r5y = zeros(B,T);
r7y = zeros(B,T);
r10y = zeros(B,T);
r15y = zeros(B,T);
r20y = zeros(B,T);
r30y = zeros(B,T);

% Loop over Bootstraps
for b = 1:B
    u = zeros(T,1);
    u(1) = ceil((n-1)*rand);
    r3m(b,1) = ZZ(u(1),1);
    r6m(b,1) = ZZ(u(1),2);
    r12m(b,1)= ZZ(u(1),3);
    r2y(b,1) = ZZ(u(1),4);
    r3y(b,1) = ZZ(u(1),5);
    r5y(b,1) = ZZ(u(1),6);
    r7y(b,1) = ZZ(u(1),7);
    r10y(b,1) = ZZ(u(1),8);
    r15y(b,1) = ZZ(u(1),9);
    r20y(b,1) = ZZ(u(1),10);
    r30y(b,1) = ZZ(u(1),11);
    
    for t = 2:T 
        if rand <= JF    
           u(t) = ceil((n-1)*rand);
        elseif u(t-1) == n-1
           u(t) = 1;
        else
           u(t) = u(t-1) + 1;  
        end
        r3m(b,t) = ZZ(u(t),1);
        r6m(b,t) = ZZ(u(t),2);
        r12m(b,t) =ZZ(u(t),3);
        r2y(b,t) = ZZ(u(t),4);
        r3y(b,t) = ZZ(u(t),5);
        r5y(b,t) = ZZ(u(t),6);
        r7y(b,t) = ZZ(u(t),7);
        r10y(b,t) = ZZ(u(t),8);
        r15y(b,t) = ZZ(u(t),9);
        r20y(b,t) = ZZ(u(t),10);
        r30y(b,t) = ZZ(u(t),11);
    end
end

% Autocorrelation of Simulated Sampling Data 
for b=1:B
    ac_3m_s = corrcoef(r3m(b,2:end),r3m(b,1:end-1));
    ac_3m_ss(b) = ac_3m_s(1,2);

    ac_6m_s = corrcoef(r6m(b,2:end),r6m(b,1:end-1));
    ac_6m_ss(b) = ac_6m_s(1,2);
    
    ac_12m_s = corrcoef(r12m(b,2:end),r12m(b,1:end-1));
    ac_12m_ss(b) = ac_12m_s(1,2);
    
    ac_2y_s = corrcoef(r2y(b,2:end),r2y(b,1:end-1));
    ac_2y_ss(b) = ac_2y_s(1,2);

    ac_3y_s = corrcoef(r3y(b,2:end),r3y(b,1:end-1));
    ac_3y_ss(b) = ac_3y_s(1,2);
    
    ac_5y_s = corrcoef(r5y(b,2:end),r5y(b,1:end-1));
    ac_5y_ss(b) = ac_5y_s(1,2);
    
    ac_7y_s = corrcoef(r7y(b,2:end),r7y(b,1:end-1));
    ac_7y_ss(b) = ac_7y_s(1,2);

    ac_10y_s = corrcoef(r10y(b,2:end),r10y(b,1:end-1));
    ac_10y_ss(b) = ac_10y_s(1,2);

    ac_15y_s = corrcoef(r15y(b,2:end),r15y(b,1:end-1));
    ac_15y_ss(b) = ac_15y_s(1,2);
    
    ac_20y_s = corrcoef(r20y(b,2:end),r20y(b,1:end-1));
    ac_20y_ss(b) = ac_20y_s(1,2);
    
    ac_30y_s = corrcoef(r30y(b,2:end),r30y(b,1:end-1));
    ac_30y_ss(b) = ac_30y_s(1,2);
end

ac_s = [ac_3m_ss' ac_6m_ss' ac_12m_ss' ac_2y_ss' ac_3y_ss' ac_5y_ss' ac_7y_ss' ac_10y_ss' ac_15y_ss' ...
        ac_20y_ss' ac_30y_ss'];

ac_ss = mean(ac_s);
xlswrite('Simulated Autocorrelation_JF=1000.xlsx',ac_ss');

% a = zeros(2,11);
% for i=1:11
%     a(1,i) = sum(ZZ(:,i)>0);
%     a(2,i) = sum(ZZ(:,i)<0);    
% end

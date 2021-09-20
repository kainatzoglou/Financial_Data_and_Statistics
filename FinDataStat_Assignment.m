%% LECTURE 1 EXERCISE OF THE ASSIGNMENT

%% OK READING THE DATA OF AMZN FROM THE CSV FILE
clear all
close all

AMZ = readtable('AMAZON.csv');
% Convert the table format to datetime type for the first column and double
% for the rest four columns
date = table2array(AMZ(:,1));
price = table2array(AMZ(:,2));
trd_vol = table2array(AMZ(:,3));
ofi_buy = table2array(AMZ(:,4));
ofi_sell = table2array(AMZ(:,5));
date_serial = datenum(date);

%% OK Reading the data from APPLE

APL = readtable('APPLE.csv');
% Convert the table format to datetime type for the first column and double
% for the rest four columns
date1 = table2array(APL(:,1));
price1 = table2array(APL(:,2));
trd_vol1 = table2array(APL(:,3));
ofi_buy1 = table2array(APL(:,4));
ofi_sell1 = table2array(APL(:,5));
date_serial1 = datenum(date1);


%% OK Reading the data from Google

GOOGL = readtable('GOOGLE.csv');
% Convert the table format to datetime type for the first column and double
% for the rest four columns
date2 = table2array(GOOGL(:,1));
price2 = table2array(GOOGL(:,2));
trd_vol2 = table2array(GOOGL(:,3));
ofi_buy2 = table2array(GOOGL(:,4));
ofi_sell2 = table2array(GOOGL(:,5));
date_serial2 = datenum(date2);

%% OK Reading the data from Microsoft

MICR = readtable('MICROSOFT.csv');
% Convert the table format to datetime type for the first column and double
% for the rest four columns
date3 = table2array(MICR(:,1));
price3 = table2array(MICR(:,2));
trd_vol3 = table2array(MICR(:,3));
ofi_buy3 = table2array(MICR(:,4));
ofi_sell3 = table2array(MICR(:,5));
date_serial3 = datenum(date3);

%% OK Reading the data from Tesla

TSL = readtable('TESLA.csv');
% Convert the table format to datetime type for the first column and double
% for the rest four columns
date4 = table2array(TSL(:,1));
price4 = table2array(TSL(:,2));
trd_vol4 = table2array(TSL(:,3));
ofi_buy4 = table2array(TSL(:,4));
ofi_sell4 = table2array(TSL(:,5));
date_serial4 = datenum(date4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   SECTION 1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OK Figure 1: Plotting the datasets for the whole Amazon dataset

figure (1)
subplot(2,2,1)
plot(date(1:length(date)),price(1:length(price)),'Linewidth',1.2);
xlabel("Dates")
ylabel("Amazon Prices")

%figure (2)
subplot(2,2,2)
plot(date(1:length(date)),trd_vol(1:length(trd_vol)));
xlabel("Dates")
ylabel("Amazon Traded Volume")

%figure (3)
subplot(2,2,3)
plot(date(1:length(date)),ofi_buy(1:length(trd_vol)));
xlabel("Dates")
ylabel("Amazon Buy Order Flow")

%figure (4)
subplot(2,2,4)
plot(date(1:length(date)),ofi_sell(1:length(trd_vol)));
xlabel("Dates")
ylabel("Amazon Sell Order Flow")

%% OK Figure 2: Plotting the datasets for the whole Apple dataset

figure(2)
subplot(2,2,1)
plot(date1(1:length(date1)),price1(1:length(price1)));
xlabel("Dates")
ylabel("Apple Prices")

subplot(2,2,2)
plot(date1(1:length(date1)),trd_vol1(1:length(trd_vol1)));
xlabel("Dates")
ylabel("Apple Traded Volume")

subplot(2,2,3)
plot(date1(1:length(date1)),ofi_buy1(1:length(trd_vol1)));
xlabel("Dates")
ylabel("Aple Buy Order Flow")
 
subplot(2,2,4)
plot(date1(1:length(date1)),ofi_sell1(1:length(trd_vol1)));
xlabel("Dates")
ylabel("Apple Sell Order Flow")
 
%% OK Figure 3: Plotting the datasets for the whole Google dataset

figure(3)
subplot(2,2,1)
plot(date2(1:length(date2)),price2(1:length(price2)));
xlabel("Dates")
ylabel("Google Prices")

subplot(2,2,2)
plot(date2(1:length(date2)),trd_vol2(1:length(trd_vol2)));
xlabel("Dates")
ylabel("Google Traded Volume")

subplot(2,2,3)
plot(date2(1:length(date2)),ofi_buy2(1:length(trd_vol2)));
xlabel("Dates")
ylabel("Google Buy Order Flow")
 
subplot(2,2,4)
plot(date2(1:length(date2)),ofi_sell2(1:length(trd_vol2)));
xlabel("Dates")
ylabel("Google Sell Order Flow")


%% OK Figure 4: Plotting the datasets for the whole Microsoft dataset

figure(4)
subplot(2,2,1)
plot(date3(1:length(date3)),price3(1:length(price3)));
xlabel("Dates")
ylabel("Microsoft Prices")

subplot(2,2,2)
plot(date3(1:length(date3)),trd_vol3(1:length(trd_vol3)));
xlabel("Dates")
ylabel("Microsoft Traded Volume")

subplot(2,2,3)
plot(date3(1:length(date3)),ofi_buy3(1:length(trd_vol3)));
xlabel("Dates")
ylabel("Microsoft Buy Order Flow")
 
subplot(2,2,4)
plot(date3(1:length(date3)),ofi_sell3(1:length(trd_vol3)));
xlabel("Dates")
ylabel("Microsoft Sell Order Flow")


%% OK Figure 5: Plotting the datasets for the whole Tesla dataset

figure(5)
subplot(2,2,1)
plot(date4(1:length(date4)),price4(1:length(price4)));
xlabel("Dates")
ylabel("Tesla Prices")

subplot(2,2,2)
plot(date4(1:length(date4)),trd_vol4(1:length(trd_vol4)));
xlabel("Dates")
ylabel("Tesla Traded Volume")

subplot(2,2,3)
plot(date4(1:length(date4)),ofi_buy4(1:length(trd_vol4)));
xlabel("Dates")
ylabel("Tesla Buy Order Flow")
 
subplot(2,2,4)
plot(date4(1:length(date4)),ofi_sell4(1:length(trd_vol4)));
xlabel("Dates")
ylabel("Tesla Sell Order Flow")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   SECTION 2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OK Compute Log - returns of the 5 stocks
l=568913;
r = log(price(52:982)./price(51:981)); % Log-returns
 
r0 = log(price(2:l)./price(1:l-1)); 
r1 = log(price1(2:l)./price1(1:l-1)); % Log-returns
r2 = log(price2(2:l)./price2(1:l-1)); % Log-returns
r3 = log(price3(2:l)./price3(1:l-1)); % Log-returns
r4 = log(price4(2:l)./price4(1:l-1)); % Log-returns

% OK Compute Trade Volume of the 5 stocks
volume= trd_vol(51:981);
v0 = trd_vol(1:l) ;
v1 = trd_vol1(1:l);
v2 = trd_vol2(1:l);
v3 = trd_vol3(1:l);
v4 = trd_vol4(1:l);


% OK Figure 6: Compute Volatility of the 5 stocks
x = 360;
y = 1;
vol = zeros(1608,1);
vol0 = zeros((length(r0)-mod(length(r0),360))/360,1);
vol1 = zeros((length(r1)-mod(length(r1),360))/360,1);
vol2 = zeros((length(r2)-mod(length(r2),360))/360,1);
vol3 = zeros((length(r3)-mod(length(r3),360))/360,1);
vol4 = zeros((length(r4)-mod(length(r4),360))/360,1);
meand = zeros(1608,1);

x = 360;
y = 1;
for i=1:(length(r0)-mod(length(r0),360))/360
    vol0(i) = std(r0(y:x,1));
    y = y + 360;
    x = x + 360;
end
x = 360;
y = 1;
for i=1:(length(r1)-mod(length(r1),360))/360
    
    vol1(i) = std(r1(y:x,1));
   
    y = y + 360;
    x = x + 360;
end
x = 360;
y = 1;
for i=1:(length(r2)-mod(length(r2),360))/360
    
    vol2(i) = std(r2(y:x,1));
    y = y + 360;
    x = x + 360;
end
x = 360;
y = 1;
for i=1:(length(r3)-mod(length(r3),360))/360
    
    vol3(i) = std(r3(y:x,1));
    y = y + 360;
    x = x + 360;
end
x = 360;
y = 1;
for i=1:(length(r4)-mod(length(r4),360))/360
   
    vol4(i) = std(r4(y:x,1));
    y = y + 360;
    x = x + 360;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% SECTION 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure (6)
subplot(2,3,1)
plot(date(1:l-1),r0);
xlabel("Dates")
ylabel("Amazon Log-Returns")
title("Amazon Log-Returns")

subplot(2,3,2)
plot(date(1:l-1),r2);
xlabel("Dates")
ylabel("Google Log-Returns")
title("Google Log-Returns")

subplot(2,3,3)
plot(date(1:l-1),r4);
xlabel("Dates")
ylabel("Tesla Log-Returns")
title("Tesla Log-Returns")

subplot(2,3,4)
plot((1:1580),vol0);
xlabel("Hours")
ylabel("Amazon Volatility")
title("Amazon Rolling Volatility")
 
subplot(2,3,5)
plot((1:1580),vol2);
xlabel("Hours")
ylabel("Google Volatility")
title("Google Rolling Volatility")

subplot(2,3,6)
plot((1:1580),vol4);
xlabel("Hours")
ylabel("Tesla Volatility")
title("Tesla Rolling Volatility")



figure ()
subplot(2,2,1)
plot(date(1:l-1),r1);
xlabel("Dates")
ylabel("Apple Log-Returns")
title("Apple Log-Returns")

subplot(2,2,2)
plot(date(1:l-1),r3);
xlabel("Dates")
ylabel("Microsoft Log-Returns")
title("Microsoft Log-Returns")

subplot(2,2,3)
plot((1:1580),vol1);
xlabel("Hours")
ylabel("Apple Volatility")
title("Apple Rolling Volatility")

subplot(2,2,4)
plot((1:1580),vol3);
xlabel("Hours")
ylabel("Microsoft Volatility")
title("Microsoft Rolling Volatility")




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OK Figure 7: Empirical distribution of log-returns Amazon vs Gaussian

x0 = linspace(min(r0),max(r0),1000);
g0 = exp(-(x0-mean(r0)).^2/(2*std(r0)^2))/sqrt(2*pi*std(r0)^2);
figure (7)
subplot(1,3,1)
histogram(r0,200,'Normalization','pdf')
set(gca,'YScale','log');
xlim([-0.04,0.02]);
ylim([1E-003,1E+003]);
hold on
plot(x0,g0,'linewidth',2)
xlim([-0.04,0.02]);
ylim([1E-003,1E+003]);
xlabel('$r0$','Interpreter','LaTex')
ylabel('$p(r0)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Amazon LogReturns')

%% OK Figure 7b: Empirical distribution of log-returns Google vs Gaussian

x2 = linspace(min(r2),max(r2),1000);
g2 = exp(-(x2-mean(r2)).^2/(2*std(r2)^2))/sqrt(2*pi*std(r2)^2);
subplot(1,3,2)
histogram(r2,200,'Normalization','pdf')
set(gca,'YScale','log');
xlim([-0.08,0.08]);
ylim([1E-003,1E+004]);
hold on
plot(x2,g2,'linewidth',2)
xlim([-0.08,0.08]);
ylim([1E-003,1E+004]);
xlabel('$r2$','Interpreter','LaTex')
ylabel('$p(r2)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Google LogReturns')

%% OK Figure 7c: Empirical distribution of log-returns Tesla vs Gaussian
subplot(1,3,3)
x4 = linspace(min(r4),max(r4),1000);
g4 = exp(-(x4-mean(r4)).^2/(2*std(r4)^2))/sqrt(2*pi*std(r4)^2);
histogram(r4,200,'Normalization','pdf')
set(gca,'YScale','log');
xlim([-0.08,0.08]);
ylim([1E-003,1E+004]);
hold on
plot(x4,g4,'linewidth',2)
xlim([-0.08,0.08]);
ylim([1E-003,1E+004]);
xlabel('$r4$','Interpreter','LaTex')
ylabel('$p(r4)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Tesla LogReturns')


%% OK Figure 7d: Empirical distribution of log-returns Google vs Gaussian

x1 = linspace(min(r1),max(r1),1000);
g1 = exp(-(x1-mean(r1)).^2/(2*std(r1)^2))/sqrt(2*pi*std(r1)^2);
figure()
subplot(1,2,1)
histogram(r1,200,'Normalization','pdf')
set(gca,'YScale','log');
xlim([-0.08,0.08]);
ylim([1E-003,1E+004]);
hold on
plot(x1,g1,'linewidth',2)
xlim([-0.08,0.08]);
ylim([1E-003,1E+004]);
xlabel('$r1$','Interpreter','LaTex')
ylabel('$p(r1)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Apple LogReturns')

%% OK Figure 7e: Empirical distribution of log-returns Tesla vs Gaussian
subplot(1,2,2)
x3 = linspace(min(r3),max(r3),1000);
g3 = exp(-(x3-mean(r3)).^2/(2*std(r3)^2))/sqrt(2*pi*std(r3)^2);
histogram(r3,200,'Normalization','pdf')
set(gca,'YScale','log');
xlim([-0.08,0.08]);
ylim([1E-003,1E+004]);
hold on
plot(x3,g3,'linewidth',2)
xlim([-0.08,0.08]);
ylim([1E-003,1E+004]);
xlabel('$r3$','Interpreter','LaTex')
ylabel('$p(r3)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Microsoft LogReturns')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OK Figure 8:Analyze the distribution of Amazon trade volume
y0 = linspace(min(v0),max(v0),1000);
f0 = exp(-(y0-mean(v0)).^2/(2*std(v0)^2))/sqrt(2*pi*std(v0)^2);

%%%%%%% this is for gpd %%%%%%%%%%%%%%%%%
posv0=v0(v0 ~= 0);
paramEsts0 = gpfit(posv0);
kHat0      = paramEsts0(1)   % Tail index parameter
sigmaHat0  = paramEsts0(2)

parmHat0 = wblfit(posv0);
pHat0 = lognfit(posv0);
phat0 = gamfit(posv0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8)
subplot(1,3,1)
histogram(v0,200,'Normalization','pdf')
xlim([0,8000]);
ylim([0,2E-03]);
hold on
plot(y0,f0,'linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
line(y0,gppdf(y0,kHat0,sigmaHat0),'linewidth',2,'color','blue');
% hold on
% line(y0,wblpdf(y0,parmHat0(1),parmHat0(2)+3),'linewidth',2,'color','yellow');
hold on
line(y0,lognpdf(y0,pHat0(1),pHat0(2)),'linewidth',2,'color','green');
hold on
line(y0,gampdf(y0,phat0(1),phat0(2)),'linewidth',2,'color','yellow');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim([0,8000]);
ylim([0,2E-03]);
xlabel('$v0$','Interpreter','LaTex')
ylabel('$p(v0)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Amazon Trading Volume')


%% OK Figure 8b: Analyze the distribution of Google trade volume
y2 = linspace(min(v2),max(v2),1000);
f2 = exp(-(y2-mean(v2)).^2/(2*std(v2)^2))/sqrt(2*pi*std(v2)^2);

%%%%%%% this is for gpd %%%%%%%%%%%%%%%%%
posv2=v2(v2 ~= 0);
paramEsts2 = gpfit(posv2);
kHat2      = paramEsts2(1)   % Tail index parameter
sigmaHat2  = paramEsts2(2)
pHat2 = lognfit(posv2);
phat2 = gamfit(posv2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,3,2)
histogram(v2,200,'Normalization','pdf')
xlim([0,4000]);
ylim([0,4E-03]);
hold on
plot(y2,f2,'linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
line(y2,gppdf(y2,kHat2,sigmaHat2),'linewidth',2,'color','blue');
hold on
line(y2,lognpdf(y2,pHat2(1),pHat2(2)),'linewidth',2,'color','green');
hold on
line(y2,gampdf(y2,phat2(1),phat2(2)),'linewidth',2,'color','yellow');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim([0,4000]);
ylim([0,4E-03]);
xlabel('$v2$','Interpreter','LaTex')
ylabel('$p(v2)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Google Trading Volume')


%% OK Figure 8c:Analyze the distribution of Tesla trade volume
y4 = linspace(min(v4),max(v4),1000);
f4 = exp(-(y4-mean(v4)).^2/(2*std(v4)^2))/sqrt(2*pi*std(v4)^2);

%%%%%%% this is for gpd %%%%%%%%%%%%%%%%%
posv4=v4(v4 ~= 0);
paramEsts4 = gpfit(posv4);
kHat4      = paramEsts4(1)   % Tail index parameter
sigmaHat4  = paramEsts4(2)

pHat4 = lognfit(posv4);
phat4 = gamfit(posv4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8)
subplot(1,3,3)
histogram(v4,200,'Normalization','pdf')
xlim([0,8000]);
ylim([0,2E-03]);
hold on
plot(y4,f4,'linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
line(y4,gppdf(y4,kHat4,sigmaHat4),'linewidth',2,'color','blue');

hold on
line(y4,lognpdf(y4,pHat4(1),pHat4(2)),'linewidth',2,'color','green');
hold on
line(y4,gampdf(y4,phat4(1),phat4(2)),'linewidth',2,'color','yellow');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim([0,8000]);
ylim([0,2E-03]);
xlabel('$v4$','Interpreter','LaTex')
ylabel('$p(v4)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Tesla Trading Volume')


%% OK Figure 8d: Analyze the distribution of Apple trade volume
y1 = linspace(min(v1),max(v1),1000);
f1 = exp(-(y1-mean(v1)).^2/(2*std(v1)^2))/sqrt(2*pi*std(v1)^2);

%%%%%%% this is for gpd %%%%%%%%%%%%%%%%%
posv1=v1(v1 ~= 0);
paramEsts1 = gpfit(posv1);
kHat1      = paramEsts1(1)   % Tail index parameter
sigmaHat1  = paramEsts1(2)
pHat1 = lognfit(posv1);
phat1 = gamfit(posv1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8)
subplot(1,2,1)
histogram(v1,200,'Normalization','pdf')
% xlim([0,4000]);
% ylim([0,4E-03]);
hold on
plot(y1,f1,'linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
line(y1,gppdf(y1,kHat1,sigmaHat1),'linewidth',2,'color','blue');

hold on
line(y1,lognpdf(y1,pHat1(1),pHat1(2)),'linewidth',2,'color','green');
hold on
line(y1,gampdf(y1,phat1(1),phat1(2)),'linewidth',2,'color','yellow');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim([0,4000]);
ylim([0,3.25E-04]);
xlabel('$v1$','Interpreter','LaTex')
ylabel('$p(v1)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Apple Trading Volume')


%% OK Figure 8e: Analyze the distribution of Microsoft trade volume
y3 = linspace(min(v3),max(v3),1000);
f3 = exp(-(y3-mean(v3)).^2/(2*std(v3)^2))/sqrt(2*pi*std(v3)^2);

%%%%%%% this is for gpd %%%%%%%%%%%%%%%%%
posv3=v3(v3 ~= 0);
paramEsts3 = gpfit(posv3);
kHat3     = paramEsts3(1)   % Tail index parameter
sigmaHat3  = paramEsts3(2)

pHat3 = lognfit(posv3);
phat3 = gamfit(posv3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)
histogram(v3,200,'Normalization','pdf')
xlim([0,4000]);
ylim([0,3.25E-04]);
hold on
plot(y3,f3,'linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
line(y3,gppdf(y3,kHat3,sigmaHat3),'linewidth',2,'color','blue');

hold on
line(y3,lognpdf(y3,pHat3(1),pHat3(2)),'linewidth',2,'color','green');
hold on
line(y3,gampdf(y3,phat3(1),phat3(2)),'linewidth',2,'color','yellow');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xlim([0,4000]);
% ylim([0,4E-03]);
xlabel('$v3$','Interpreter','LaTex')
ylabel('$p(v3)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Microsoft Trading Volume')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% SECTION 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% OK Figure 9: Analyze the distribution of Amazon volatility 
vol0=vol0';
z0 = linspace(min(vol0),max(vol0),1000);
h0 = exp(-(z0-mean(vol0)).^2/(2*std(vol0)^2))/sqrt(2*pi*std(vol0)^2);

%%%%%%% this is for lognormal %%%%%%%%%%%%%%%%%
posvol0=vol0(vol0 ~= 0);
pHat = lognfit(posvol0)
phat0 = gamfit(posvol0);
paramEsts0 = gpfit(posvol0);
kHat0     = paramEsts0(1)   % Tail index parameter
sigmaHat0  = paramEsts0(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure (9)
subplot(1,3,1)
histogram(vol0,200,'Normalization','pdf')
% % xlim([0,1E-03]);
% % ylim([1E+02,1E+04]);
hold on
plot(z0,h0,'linewidth',2)
% % xlim([0,1E-03]);
% % ylim([1E+02,1E+04]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
line(z0,lognpdf(z0,pHat(1),pHat(2)),'linewidth',2,'color','green');

hold on
line(z0,gppdf(z0,kHat0,sigmaHat0),'linewidth',2,'color','blue');
hold on
line(z0,gampdf(z0,phat0(1),phat0(2)),'linewidth',2,'color','yellow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('$vol0$','Interpreter','LaTex')
ylabel('$p(vol0)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Amazon Volatility')

%% OK Figure 9b: Analyze the distribution of Google Volatility

z2 = linspace(min(vol2),max(vol2),1000);
h2 = exp(-(z2-mean(vol2)).^2/(2*std(vol2)^2))/sqrt(2*pi*std(vol2)^2);

%%%%%%% this is for lognormal %%%%%%%%%%%%%%%%%
posvol2=vol2(vol2 ~= 0);
pHat2 = lognfit(posvol2)

phat2 = gamfit(posvol2);
paramEsts2 = gpfit(posvol2);
kHat2     = paramEsts2(1)   % Tail index parameter
sigmaHat2  = paramEsts2(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,3,2)
histogram(vol2,200,'Normalization','pdf')
% % xlim([0,1E-03]);
% % ylim([1E+02,1E+04]);
hold on
plot(z2,h2,'linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
line(z2,lognpdf(z2,pHat2(1),pHat2(2)),'linewidth',2,'color','green');

hold on
line(z2,gppdf(z2,kHat2,sigmaHat2),'linewidth',2,'color','blue');
hold on
line(z2,gampdf(z2,phat2(1),phat2(2)),'linewidth',2,'color','yellow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % xlim([0,1E-03]);
% % ylim([1E+02,1E+04]);
xlabel('$vol2$','Interpreter','LaTex')
ylabel('$p(vol2)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Google Volatility')


%% OK Figure 9c: Analyze the distribution of Tesla Volatility

z4 = linspace(min(vol4),max(vol4),1000);
h4 = exp(-(z4-mean(vol4)).^2/(2*std(vol4)^2))/sqrt(2*pi*std(vol4)^2);

%%%%%%% this is for lognormal %%%%%%%%%%%%%%%%%
posvol4=vol4(vol4 ~= 0);
pHat4 = lognfit(posvol4)

phat4 = gamfit(posvol4);
paramEsts4 = gpfit(posvol4);
kHat4     = paramEsts4(1)   % Tail index parameter
sigmaHat4  = paramEsts4(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,3,3)
histogram(vol4,200,'Normalization','pdf')
% % xlim([0,1E-03]);
% % ylim([1E+02,1E+04]);
hold on
plot(z4,h4,'linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
line(z4,lognpdf(z4,pHat4(1),pHat4(2)),'linewidth',2,'color','green');

hold on
line(z4,gppdf(z4,kHat4,sigmaHat4),'linewidth',2,'color','blue');
hold on
line(z4,gampdf(z4,phat4(1),phat4(2)),'linewidth',2,'color','yellow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % xlim([0,1E-03]);
% % ylim([1E+02,1E+04]);
xlabel('$vol4$','Interpreter','LaTex')
ylabel('$p(vol4)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Tesla Volatility')



%% OK Figure 9b: Analyze the distribution of Apple Volatility

z1 = linspace(min(vol1),max(vol1),1000);
h1 = exp(-(z1-mean(vol1)).^2/(2*std(vol1)^2))/sqrt(2*pi*std(vol1)^2);

%%%%%%% this is for lognormal %%%%%%%%%%%%%%%%%
posvol1=vol1(vol1 ~= 0);
pHat1 = lognfit(posvol1)

phat1 = gamfit(posvol1);
paramEsts1 = gpfit(posvol1);
kHat1     = paramEsts1(1)   % Tail index parameter
sigmaHat1  = paramEsts1(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1)
histogram(vol1,200,'Normalization','pdf')
% % xlim([0,1E-03]);
% % ylim([1E+02,1E+04]);
hold on
plot(z1,h1,'linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
line(z1,lognpdf(z1,pHat1(1),pHat1(2)),'linewidth',2,'color','green');

hold on
line(z1,gppdf(z1,kHat1,sigmaHat1),'linewidth',2,'color','blue');
hold on
line(z1,gampdf(z1,phat1(1),phat1(2)),'linewidth',2,'color','yellow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % xlim([0,1E-03]);
% % ylim([1E+02,1E+04]);
xlabel('$vol2$','Interpreter','LaTex')
ylabel('$p(vol2)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Apple Volatility')


%% OK Figure 9c: Analyze the distribution of Microsoft Volatility

z3 = linspace(min(vol3),max(vol3),1000);
h3 = exp(-(z3-mean(vol3)).^2/(2*std(vol3)^2))/sqrt(2*pi*std(vol3)^2);

%%%%%%% this is for lognormal %%%%%%%%%%%%%%%%%
posvol3=vol3(vol3 ~= 0);
pHat3 = lognfit(posvol3)

phat3 = gamfit(posvol3);
paramEsts3 = gpfit(posvol3);
kHat3     = paramEsts3(1)   % Tail index parameter
sigmaHat3  = paramEsts3(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)
histogram(vol3,200,'Normalization','pdf')
% % xlim([0,1E-03]);
% % ylim([1E+02,1E+04]);
hold on
plot(z3,h3,'linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
line(z3,lognpdf(z3,pHat3(1),pHat3(2)),'linewidth',2,'color','green');

hold on
line(z3,gppdf(z3,kHat3,sigmaHat3),'linewidth',2,'color','blue');
hold on
line(z3,gampdf(z3,phat3(1),phat3(2)),'linewidth',2,'color','yellow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % xlim([0,1E-03]);
% % ylim([1E+02,1E+04]);
xlabel('$vol4$','Interpreter','LaTex')
ylabel('$p(vol4)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
%set(gca,'FontSize',20)
title('PDF of Microsoft Volatility')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% SECTION 7%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% OK Descriptive statistics for log - returns Amazon %%%%%%%%%%
fprintf('Descriptive statistics for Amazon: \n','FontWeight','bold');
N = length(r0);
m = sum(r0)/N; % Compute mean and store value in variable
fprintf('\n')
fprintf('Mean = %4.3f\n',m)
s = sqrt(sum((r0-m).^2)/N); % Compute std. deviation and store value in variable
fprintf('Std. deviation = %4.3f\n',s)
fprintf('Skewness = %4.3f\n',sum((r0-m).^3)/(N*s^3))
fprintf('Excess kurtosis = %4.3f\n',sum((r0-m).^4)/(N*s^4)-3) 
fprintf('\n')

%% OK Descriptive statistics for trade volume Amazon %%%%%%%%%%
N_volume = length(v0); % Number of elements of trade volume
m_volume = sum(v0)/N_volume; % Compute mean and store value in variable
fprintf('\n')
fprintf('Volume mean = %4.3f\n',m_volume)
s_volume = sqrt(sum((v0-m_volume).^2)/N_volume); % Compute std. deviation and store value in variable
fprintf('Std. deviation of volume = %4.3f\n',s_volume)
fprintf('Skewness = %4.3f\n',sum((v0-m_volume).^3)/(N_volume*s_volume^3))
fprintf('Excess kurtosis = %4.3f\n',sum((v0-m_volume).^4)/(N_volume*s_volume^4)-3)
fprintf('\n')

%% OK Descriptive statistics for volatility Amazon%%%%%%%%%%

N_volatility = length(vol0); % Number of elements of trade volume
m_volatility = sum(vol0)/N_volatility; % Compute mean and store value in variable
fprintf('\n')
fprintf('Volatility mean = %4.3f\n',m_volatility)
s_volatility = sqrt(sum((vol0-m_volatility).^2)/N_volatility); % Compute std. deviation and store value in variable
fprintf('Std. deviation of volatility = %4.3f\n',s_volatility)
fprintf('Skewness = %4.3f\n',sum((vol0-m_volatility).^3)/(N_volatility*s_volatility^3))
fprintf('Excess kurtosis = %4.3f\n',sum((vol0-m_volatility).^4)/(N_volatility*s_volatility^4)-3)
fprintf('\n')

%% OK Descriptive statistics for log - returns Google %%%%%%%%%%

fprintf('Descriptive statistics for Google: \n','FontWeight','bold');
N2= length(r2);
m2= sum(r2)/N2; % Compute mean and store value in variable
fprintf('\n')
fprintf('Mean = %4.3f\n',m2)
s2= sqrt(sum((r2-m2).^2)/N2); % Compute std. deviation and store value in variable
fprintf('Std. deviation = %4.3f\n',s2)
fprintf('Skewness = %4.3f\n',sum((r2-m2).^3)/(N2*s2^3))
fprintf('Excess kurtosis = %4.3f\n',sum((r2-m2).^4)/(N2*s2^4)-3)
fprintf('\n')

%% OK Descriptive statistics for trade volume Google %%%%%%%%%%
N_volume2 = length(v2); % Number of elements of trade volume
m_volume2 = sum(v2)/N_volume2; % Compute mean and store value in variable
fprintf('\n')
fprintf('Volume mean = %4.3f\n',m_volume2)
s_volume2 = sqrt(sum((v2-m_volume2).^2)/N_volume2); % Compute std. deviation and store value in variable
fprintf('Std. deviation of volume = %4.3f\n',s_volume2)
fprintf('Skewness = %4.3f\n',sum((v2-m_volume2).^3)/(N_volume2*s_volume2^3))
fprintf('Excess kurtosis = %4.3f\n',sum((v2-m_volume2).^4)/(N_volume2*s_volume2^4)-3)
fprintf('\n')

%% OK Descriptive statistics for volatility Google%%%%%%%%%%
N_volatility2 = length(vol2); % Number of elements of trade volume
m_volatility2 = sum(vol2)/N_volatility2; % Compute mean and store value in variable
fprintf('\n')
fprintf('Volatility mean = %4.3f\n',m_volatility2)
s_volatility2 = sqrt(sum((vol2-m_volatility2).^2)/N_volatility2); % Compute std. deviation and store value in variable
fprintf('Std. deviation of volatility = %4.3f\n',s_volatility2)
fprintf('Skewness = %4.3f\n',sum((vol2-m_volatility2).^3)/(N_volatility2*s_volatility2^3)) 
fprintf('Excess kurtosis = %4.3f\n',sum((vol2-m_volatility2).^4)/(N_volatility2*s_volatility2^4)-3) 
fprintf('\n')


%% OK Descriptive statistics for log - returns Tesla %%%%%%%%%%

fprintf('Descriptive statistics for Tesla: \n','FontWeight','bold');
N4= length(r4);
m4= sum(r4)/N4; % Compute mean and store value in variable
fprintf('\n')
fprintf('Mean = %4.3f\n',m4)
s4= sqrt(sum((r4-m4).^2)/N4); % Compute std. deviation and store value in variable
fprintf('Std. deviation = %4.3f\n',s4)
fprintf('Skewness = %4.3f\n',sum((r4-m4).^3)/(N4*s4^3))
fprintf('Excess kurtosis = %4.3f\n',sum((r4-m4).^4)/(N4*s4^4)-3)
fprintf('\n')

%% OK Descriptive statistics for trade volume Tesla %%%%%%%%%%
N_volume4 = length(v4); % Number of elements of trade volume
m_volume4 = sum(v4)/N_volume4; % Compute mean and store value in variable
fprintf('\n')
fprintf('Volume mean = %4.3f\n',m_volume4)
s_volume4 = sqrt(sum((v4-m_volume4).^2)/N_volume4); % Compute std. deviation and store value in variable
fprintf('Std. deviation of volume = %4.3f\n',s_volume4)
fprintf('Skewness = %4.3f\n',sum((v4-m_volume4).^3)/(N_volume4*s_volume4^3))
fprintf('Excess kurtosis = %4.3f\n',sum((v4-m_volume4).^4)/(N_volume4*s_volume4^4)-3)
fprintf('\n')

%% OK Descriptive statistics for volatility Tesla%%%%%%%%%%
N_volatility4 = length(vol4); % Number of elements of trade volume
m_volatility4 = sum(vol4)/N_volatility4; % Compute mean and store value in variable
fprintf('\n')
fprintf('Volatility mean = %4.3f\n',m_volatility4)
s_volatility4 = sqrt(sum((vol4-m_volatility4).^2)/N_volatility4); % Compute std. deviation and store value in variable
fprintf('Std. deviation of volatility = %4.3f\n',s_volatility4)
fprintf('Skewness = %4.3f\n',sum((vol4-m_volatility4).^3)/(N_volatility4*s_volatility4^3)) 
fprintf('Excess kurtosis = %4.3f\n',sum((vol4-m_volatility4).^4)/(N_volatility4*s_volatility4^4)-3) 
fprintf('\n')


%% OK Descriptive statistics for log - returns Apple %%%%%%%%%%

fprintf('Descriptive statistics for Apple: \n','FontWeight','bold');
N1= length(r1);
m1= sum(r1)/N1; % Compute mean and store value in variable
fprintf('\n')
fprintf('Mean = %4.3f\n',m1)
s1= sqrt(sum((r1-m1).^2)/N1); % Compute std. deviation and store value in variable
fprintf('Std. deviation = %4.3f\n',s1)
fprintf('Skewness = %4.3f\n',sum((r1-m1).^3)/(N1*s1^3))
fprintf('Excess kurtosis = %4.3f\n',sum((r1-m1).^4)/(N1*s1^4)-3)
fprintf('\n')

%% OK Descriptive statistics for trade volume Apple %%%%%%%%%%
N_volume1 = length(v1); % Number of elements of trade volume
m_volume1 = sum(v1)/N_volume1; % Compute mean and store value in variable
fprintf('\n')
fprintf('Volume mean = %4.3f\n',m_volume1)
s_volume1 = sqrt(sum((v1-m_volume1).^2)/N_volume1); % Compute std. deviation and store value in variable
fprintf('Std. deviation of volume = %4.3f\n',s_volume1)
fprintf('Skewness = %4.3f\n',sum((v1-m_volume1).^3)/(N_volume1*s_volume1^3))
fprintf('Excess kurtosis = %4.3f\n',sum((v1-m_volume1).^4)/(N_volume1*s_volume1^4)-3)
fprintf('\n')

%% OK Descriptive statistics for volatility Apple%%%%%%%%%%
N_volatility1 = length(vol1); % Number of elements of trade volume
m_volatility1 = sum(vol1)/N_volatility1; % Compute mean and store value in variable
fprintf('\n')
fprintf('Volatility mean = %4.3f\n',m_volatility1)
s_volatility1 = sqrt(sum((vol1-m_volatility1).^2)/N_volatility1); % Compute std. deviation and store value in variable
fprintf('Std. deviation of volatility = %4.3f\n',s_volatility1)
fprintf('Skewness = %4.3f\n',sum((vol1-m_volatility1).^3)/(N_volatility1*s_volatility1^3)) 
fprintf('Excess kurtosis = %4.3f\n',sum((vol1-m_volatility1).^4)/(N_volatility1*s_volatility1^4)-3) 
fprintf('\n')


%% OK Descriptive statistics for log - returns Microsoft %%%%%%%%%%

fprintf('Descriptive statistics for Microsoft: \n','FontWeight','bold');
N3= length(r3);
m3= sum(r3)/N3; % Compute mean and store value in variable
fprintf('\n')
fprintf('Mean = %4.3f\n',m3)
s3= sqrt(sum((r3-m3).^2)/N3); % Compute std. deviation and store value in variable
fprintf('Std. deviation = %4.3f\n',s3)
fprintf('Skewness = %4.3f\n',sum((r3-m3).^3)/(N3*s3^3))
fprintf('Excess kurtosis = %4.3f\n',sum((r3-m3).^4)/(N3*s3^4)-3)
fprintf('\n')

%% OK Descriptive statistics for trade volume Microsoft %%%%%%%%%%
N_volume3 = length(v3); % Number of elements of trade volume
m_volume3 = sum(v3)/N_volume3; % Compute mean and store value in variable
fprintf('\n')
fprintf('Volume mean = %4.3f\n',m_volume3)
s_volume3 = sqrt(sum((v3-m_volume3).^2)/N_volume3); % Compute std. deviation and store value in variable
fprintf('Std. deviation of volume = %4.3f\n',s_volume3)
fprintf('Skewness = %4.3f\n',sum((v3-m_volume3).^3)/(N_volume3*s_volume3^3))
fprintf('Excess kurtosis = %4.3f\n',sum((v3-m_volume3).^4)/(N_volume3*s_volume3^4)-3)
fprintf('\n')

%% OK Descriptive statistics for volatility Microsoft%%%%%%%%%%
N_volatility3 = length(vol3); % Number of elements of trade volume
m_volatility3 = sum(vol3)/N_volatility3; % Compute mean and store value in variable
fprintf('\n')
fprintf('Volatility mean = %4.3f\n',m_volatility3)
s_volatility3 = sqrt(sum((vol3-m_volatility3).^2)/N_volatility3); % Compute std. deviation and store value in variable
fprintf('Std. deviation of volatility = %4.3f\n',s_volatility3)
fprintf('Skewness = %4.3f\n',sum((vol3-m_volatility3).^3)/(N_volatility3*s_volatility3^3)) 
fprintf('Excess kurtosis = %4.3f\n',sum((vol3-m_volatility3).^4)/(N_volatility3*s_volatility3^4)-3) 
fprintf('\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% SECTION 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% OK Figure 10: CCDF of log-returns for Amazon
x_ccdf0 = sort(r0); % Returns sorted in ascending order
y_ccdf0 = 1:1:length(r0); 
y_ccdf0 = 1 - y_ccdf0/(length(r0)+1); % Calculating CCDF as rank-frequency plot
N_ccdf0 = length(r0);
m_ccdf0 = sum(r0)/N_ccdf0;
s_ccdf0 = std(r0);

c0 = 0.5*(1 - erf((x_ccdf0-m_ccdf0)/(s_ccdf0*sqrt(2)))); % Gaussian CCDF
figure(10)
subplot(1,3,1)
semilogy(x_ccdf0,y_ccdf0,'b','LineWidth',2)
%semilogy(x_ccdf0,y_ccdf0,'v')
hold on
semilogy(x_ccdf0,c0,'r','LineWidth',2)
ylim([1e-4,1])
xlabel('$r$','Interpreter','LaTex')
ylabel('$C(r)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
set(gca,'FontSize',20)
title('CCDF for LogReturns- Amazon')

%% OK Figure 10b: CCDF for trade volume Amazon
x_tccdf0 = sort(v0); % Returns sorted in ascending order
y_tccdf0 = 1:1:length(v0); 
y_tccdf0 = 1 - y_tccdf0/(length(v0)+1); % Calculating CCDF as rank-frequency plot
N_tccdf0 = length(v0);
m_tccdf0 = sum(v0)/N_tccdf0;
s_tccdf0 = std(v0);
c1 = 0.5*(1 - erf((x_tccdf0-m_tccdf0)/(s_tccdf0*sqrt(2)))); % Gaussian CCDF 

subplot(1,3,2)
semilogy(x_tccdf0,y_tccdf0,'b','LineWidth',2)
hold on
semilogy(x_tccdf0,c1,'r','LineWidth',2)
ylim([1e-4,1])
xlabel('$v$','Interpreter','LaTex')
ylabel('$C(v)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
set(gca,'FontSize',20)
title('CCDF for Trading Volume-Amazon')

%% OK Figure 10c: CCDF for volatility Amazon
x_vccdf0 = sort(vol0); % Returns sorted in ascending order
y_vccdf0 = 1:1:length(vol0); 
y_vccdf0 = 1 - y_vccdf0/(length(vol0)+1); % Calculating CCDF as rank-frequency plot
N_vccdf0 = length(vol0);
m_vccdf0 = sum(vol0)/N_vccdf0;
s_vccdf0 = std(vol0);

c2 = 0.5*(1 - erf((x_vccdf0-m_vccdf0)/(s_vccdf0*sqrt(2)))); % Gaussian CCDF
 
subplot(1,3,3)
semilogy(x_vccdf0,y_vccdf0,'b','LineWidth',2)
hold on
semilogy(x_vccdf0,c2,'r','LineWidth',2)
ylim([1e-4,1])
xlabel('$v$','Interpreter','LaTex')
ylabel('$C(v)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
set(gca,'FontSize',20)
title('CCDF for Volatility-Amazon')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% SECTION 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OK Figure 11: Q-Q PLOTS for Amazon

% Q-Q PLOT FOR AMAZON RETURNS
figure (11)
subplot(1,3,1)
qqplot(r0)
title('Q-Q Plot for Amazon Returns')

% Q-Q PLOT FOR AMAZON TRADING VOLUME
subplot(1,3,2)
qqplot(v0)
title('Q-Q Plot for Amazon Trading Volume')

% Q-Q PLOT FOR VOLATILITY
subplot(1,3,3)
qqplot(vol0)
title('Q-Q Plot for Amazon Volatility')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% SECTION 10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OK Figure12: Fitting right & left tail via Maximum Likelihood  for log returns%%%%%%

p = 0.1; % Defining tails as top p% of returns (both positive and negative)

%%% Right tail %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s_r = sort(r0); % Sorting returns
r_right = s_r(round((1-p)*length(s_r)):end); % Selecting top p% of returns
N = length(r_right); % Number of returns selected as right tail
alpha_right = N/sum(log(r_right/min(r_right))); % Maximum-likelihood estimate for right tail exponent
fprintf('Right tail exponent: %4.3f\n',alpha_right)
x_right = linspace(min(r_right),max(r_right),100); % Grid of points between min and max values in right tail
y_right = alpha_right*(x_right/min(r_right)).^(-alpha_right-1)/min(r_right); % Values of power law distribution on grid of points
[b_right,a_right] = histnorm(r_right,20);


figure(12)
subplot(1,3,1)
loglog(a_right,b_right,'ob','MarkerSize',8,'MarkerFaceColor','b')
hold on
loglog(x_right,y_right,'r','LineWidth',2)
xlabel('$r$','Interpreter','LaTex')
ylabel('$p(r)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
set(gca,'FontSize',20)
title('Right tail')

%%% Left tail %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_left = s_r(1:round(p*length(s_r))); % Selecting bottom p% of returns
r_left = abs(r_left); % Converting negative returns to positive numbers 
N = length(r_left); % Number of returns selected as left tail
alpha_left = N/sum(log(r_left/min(r_left))); % Maximum-likelihood estimate for left tail exponent
fprintf('Left tail exponent: %4.3f\n',alpha_left) 
x_left = linspace(min(r_left),max(r_left),100);
y_left = alpha_left*(x_left/min(r_left)).^(-alpha_left-1)/min(r_left); % Power law distribution 
[b_left,a_left] = histnorm(r_left,20);
subplot(1,3,2)
loglog(a_left,b_left,'ob','MarkerSize',8,'MarkerFaceColor','b')
hold on
loglog(x_left,y_left,'r','LineWidth',2)
xlabel('$r$','Interpreter','LaTex')
ylabel('$p(r)$','Interpreter','LaTex')
set(gca,'TickLabelInterpreter','LaTex')
set(gca,'FontSize',20)
title('Left tail')

%% OK Bootstrap analysis 
bts = 0.8; % Fraction of data to be retained in each bootstrap sample
Nbts = 500; % Number of bootstrap samples
alpha = 0.9; % Significance level

%%% Right tail with bootstrap 
alpha_right_bts = []; % Vector to collect bootstrap estimates for right tail exponent

for i = 1:Nbts
   
    r_bts = s_r(randperm(length(s_r))); % Random permutation of returns
    r_bts = r_bts(1:round(bts*length(r_bts))); % Bootstrapping bts% of returns 
    r_bts = sort(r_bts); % Sorting bootstrapped returns
    r_right_bts = r_bts(round((1-p)*length(r_bts)):end); % Selecting top p% of returns
    
    N_bts = length(r_right_bts); % Number of bootstrapped returns
    
    alpha_right_bts = [alpha_right_bts; N_bts/sum(log(r_right_bts/min(r_right_bts)))];

end

alpha_right_bts = sort(alpha_right_bts); % Sorting bootstrap estimates for right tail exponent

fprintf('Right tail interval at %3.2f CL: [%4.3f; %4.3f] \n',alpha,alpha_right_bts(round(0.5*(1-alpha)*Nbts)),alpha_right_bts(round(0.5*(1+alpha)*Nbts)))
fprintf('\n')

%%% Left tail with bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_left_bts = []; % Vector to collect bootstrap estimates for left tail exponent

for i = 1:Nbts
   
    r_bts = s_r(randperm(length(s_r))); % Random permutation of returns
    r_bts = r_bts(1:round(bts*length(r_bts))); % Bootstrapping bts% of returns 
    r_bts = sort(r_bts); % Sorting bootstrapped returns
    r_left_bts = r_bts(1:round(p*length(r_bts))); % Selecting bottom p% of returns
    r_left_bts = abs(r_left_bts); % Converting returns to positive
    
    N_bts = length(r_left_bts); % Number of bootstrapped returns
    
    alpha_left_bts = [alpha_left_bts; N_bts/sum(log(r_left_bts/min(r_left_bts)))];

end

alpha_left_bts = sort(alpha_left_bts); % Sorting bootstrap estimates for right tail exponent

fprintf('Left tail interval at %3.2f CL: [%4.3f; %4.3f] \n',alpha,alpha_left_bts(round(0.5*(1-alpha)*Nbts)),alpha_left_bts(round(0.5*(1+alpha)*Nbts)))
fprintf('\n')

subplot(1,3,3)
boxplot([alpha_right_bts alpha_left_bts])
xticks([1 2])
xticklabels({'$\alpha_R$','$\alpha_L$'})
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','LaTex')
title('Boostrap tail exponent values')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% SECTION 11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Body-Tail fitting


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% SECTION 12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OK Covariance Matrices

covar_matr1 =zeros(5,5);
covar_matr2 =zeros(5,5);
covar_matr3 =zeros(5,5);
r_total = [r0,r1,r2,r3,r4];
volume_total = [v0,v1,v2,v3,v4];
volat_total = [vol0',vol1,vol2,vol3,vol4];
k=2;
for i=1:4
     for j=k:5
         x = cov(r_total(:,i),r_total(:,j));
         y = cov(volume_total(:,i),volume_total(:,j));
         z = cov(volat_total(:,i),volat_total(:,j));
         covar_matr1(i,j) = x(1,2);
         covar_matr2(i,j) = y(1,2);
         covar_matr3(i,j) = z(1,2);
     end
     k=k+1;
 end
     
 
 for i=1:5
     for j=1:5
         if i==j
             covar_matr1(i,j) = var(r_total(:,i));
             covar_matr2(i,j) = var(volume_total(:,i));
             covar_matr3(i,j) = var(volat_total(:,i));
         elseif i>j
             covar_matr1(i,j) = covar_matr1(j,i);
             covar_matr2(i,j) = covar_matr2(j,i);
             covar_matr3(i,j) = covar_matr3(j,i);
         end
     end
 end
 
%% OK Correlation Matrices
 
corr_matr1 =zeros(5,5);
corr_matr2 =zeros(5,5);
corr_matr3 =zeros(5,5);
r_total = [r0,r1,r2,r3,r4];
volume_total = [v0,v1,v2,v3,v4];
volat_total = [vol0',vol1,vol2,vol3,vol4];
k=2;
for i=1:4
     for j=k:5
         x = corr(r_total(:,i),r_total(:,j));
         y = corr(volume_total(:,i),volume_total(:,j));
         z = corr(volat_total(:,i),volat_total(:,j));
         corr_matr1(i,j) = x;
         corr_matr2(i,j) = y;
         corr_matr3(i,j) = z;
     end
     k=k+1;
 end
     
 
 for i=1:5
     for j=1:5
         if i==j
             corr_matr1(i,j) = 1;
             corr_matr2(i,j) = 1;
             corr_matr3(i,j) = 1;
         elseif i>j
            corr_matr1(i,j) = corr_matr1(j,i);
             corr_matr2(i,j) = corr_matr2(j,i);
             corr_matr3(i,j) = corr_matr3(j,i);
         end
     end
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION 13 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test Hypothesis on all previous results
%% MLE parameters for 
%phat=mle(r0);
%phat=mle(r0,'distribution','bino');

%% Test Hypothesis on Normalilty (Kolmogorov-Smirnov test)
% Test on log-returns

[hk_log0,pk_log0] = kstest(r0);
[hk_log2,pk_log2] = kstest(r2);
 
% Test on trade volume
[hk_tr0,pk_tr0] = kstest(v0);
[hk_tr2, pk_tr2]= kstest(v2);

% Test on volatility
[hk_volat0,pk_volat0] = kstest(vol0);
[hk_volat2,pk_volat2] = kstest(vol2);
 
%% Test Stationarity (KPSS test)
% Test on log-returns
[hs_log0,ps_log0] = kpsstest(r0);
[hs_log2,ps_log2] = kpsstest(r2);
[hs_log4,ps_log4] = kpsstest(r4);
[hs_log1,ps_log1] = kpsstest(r1);
[hs_log3,ps_log3] = kpsstest(r3);

% Test on trade volume
[hs_tr0,ps_tr0] = kpsstest(v0);
[hs_tr2,ps_tr2]= kpsstest(v2);
[hs_tr4,ps_tr4] = kpsstest(v4);
[hs_tr1,ps_tr1]= kpsstest(v1);
[hs_tr3,ps_tr3] = kpsstest(v3);

% Test on volatility
[hs_volat0,ps_volat0]= kpsstest(vol0);
[hs_volat2,ps_volat2] = kpsstest(vol2);
[hs_volat4,ps_volat4]= kpsstest(vol4);
[hs_volat1,ps_volat1] = kpsstest(vol1);
[hs_volat3,ps_volat3]= kpsstest(vol3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION 14 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Scatter plot for log-returns of Amazon and Google
l=568913;
ret0=(price(2:l)-price(1:l-1));
ret2=(price2(2:l)-price2(1:l-1));
ret0_perc=(price(2:l)-price(1:l-1))./price(1:l-1);
ret2_perc=(price2(2:l)-price2(1:l-1))./price2(1:l-1);

l1 = 142228;
figure(10)
scatter(price(1:l),price2(1:l));
xlabel('time')
ylabel('Log-returns')
title("Amazon and Google Log-Returns")

figure(13)
subplot(2,2,1)
scatter(price(1:l1),price2(1:l1));
xlabel('Amazon log-returns')
ylabel('Google log-returns')
title("First Part")
subplot(2,2,2)
scatter(price(l1+1:2*l1),price2(l1+1:2*l1))
xlabel('Amazon log-returns')
ylabel('Google log-returns')
title("Second Part")
subplot(2,2,3)
scatter(price(2*l1+1:3*l1),price2(2*l1+1:3*l1))
xlabel('Amazon log-returns')
ylabel('Google log-returns')
title("Third Part")
subplot(2,2,4)
scatter(price(3*l1+1:l),price2(3*l1+1:l))
xlabel('Amazon log-returns')
ylabel('Google log-returns')
title("Fourth Part")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION 15 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First Part - price/price dependency
part1_x = price(1:l1);
part1_y = price2(1:l1);
p=polyfit(part1_x,part1_y,1);
yfit = polyval(p,part1_x);
yresid = part1_y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(part1_y)-1) * var(part1_y);
rsq = 1 - SSresid/SStotal;
corr_rsq = sqrt(rsq);

%Pearson's correlation (linear)
 
[RHO1,PVAL1] = corr(part1_x,part1_y,'Type','Pearson');
[RHO2,PVAL2] = corr(part1_x,part1_y,'Type','Kendall');
[RHO3,PVAL3] = corr(part1_x,part1_y,'Type','Spearman');

% Second Part - price/price dependency
scatter(price(l1+1:2*l1),price2(l1+1:2*l1))
part2_x = price(l1+1:2*l1);
part2_y = price2(l1+1:2*l1);
p=polyfit(part2_x,part2_y,1);
yfit = polyval(p,part2_x);
yresid = part2_y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(part2_y)-1) * var(part2_y);
rsq = 1 - SSresid/SStotal;
corr_rsq = sqrt(rsq);

% Pearson's correlation (linear)

[RHO4,PVAL4] = corr(part2_x,part2_y,'Type','Pearson');
[RHO5,PVAL5] = corr(part2_x,part2_y,'Type','Kendall');
[RHO6,PVAL6] = corr(part2_x,part2_y,'Type','Spearman');
 
% Third Part - price/price dependency
scatter(price(2*l1+1:3*l1),price2(2*l1+1:3*l1))
part3_x = price(2*l1+1:3*l1);
part3_y = price2(2*l1+1:3*l1);
p=polyfit(part3_x,part3_y,1);
yfit = polyval(p,part3_x);
yresid = part3_y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(part3_y)-1) * var(part3_y);
rsq = 1 - SSresid/SStotal;
corr_rsq = sqrt(rsq);

%Pearson's correlation (linear)

[RHO7,PVAL7] = corr(part3_x,part3_y,'Type','Pearson');
[RHO8,PVAL8] = corr(part3_x,part3_y,'Type','Kendall');
[RHO9,PVAL9] = corr(part3_x,part3_y,'Type','Spearman');

% Fourth Part - price/price dependency
scatter(price(3*l1+1:l),price2(3*l1+1:l))
part4_x = price(3*l1+1:l);
part4_y = price2(3*l1+1:l);
p=polyfit(part4_x,part4_y,1);
yfit = polyval(p,part4_x);
yresid = part4_y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(part4_y)-1) * var(part4_y);
rsq = 1 - SSresid/SStotal;
corr_rsq = sqrt(rsq);

% Pearson's correlation (linear)

[RHO10,PVAL10] = corr(part4_x,part4_y,'Type','Pearson');
[RHO11,PVAL11] = corr(part4_x,part4_y,'Type','Kendall');
[RHO12,PVAL12] = corr(part4_x,part4_y,'Type','Spearman');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION 16 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CAUSALITY

part_xtmin1 = price(end-2000:end-1000);
part_ytmin1 = price2(end-2000:end-1000);
part_yt = price2(end-1000:end);

matrix_fit = [part_xtmin1 part_ytmin1];
mdl = fitlm(matrix_fit,part_yt);

 
matrix_fit2 = [part_ytmin1];
mdl2 = fitlm(matrix_fit2,part_yt);

b=mdl.Residuals(:,1);
c=table2array(b);
a=sum(c.^2);

var_e1 = mdl2.SSE;
var_e2 = mdl.SSE;
 
F_test = log(var_e1/var_e2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION 17 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Autocorrelation
% Amazon returns
l=568913;
returns = price(2:end) - price(1:end-1);
abs_returns = abs(returns);
sqr_returns = returns.^2;

aut_returns = autocorr(returns);
aut_absreturns = autocorr(abs_returns);
aut_sqrreturns = autocorr(sqr_returns);

figure(14)
subplot(3,2,1)
plot(returns)
ylabel('returns')
xlabel('time')
title("Returns")
subplot(3,2,2)
autocorr(returns)
subplot(3,2,3)
plot(abs_returns)
ylabel('abs-returns')
xlabel('time')
title("Absolute Returns")
subplot(3,2,4)
autocorr(abs_returns)
subplot(3,2,5)
plot(sqr_returns)
ylabel('sqr-returns')
xlabel('time')
title("Square Returns")
subplot(3,2,6)
autocorr(sqr_returns)

trading=trd_vol(1:l);
aut_trading=autocorr(trading);
aut_vol=autocorr(vol0);

figure(15)
subplot(2,2,1)
plot(trading)
ylabel('trade volume')
xlabel('time')
title("Trade Volume")
subplot(2,2,2)
autocorr(trading,'NumLags',50)
subplot(2,2,3)
plot(vol0)
ylabel('volatility')
xlabel('time')
title("Volatility")
subplot(2,2,4)
autocorr(vol0,'NumLags',50)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION 18 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Scatter plot for price and Trading Volume of Amazon
figure(16)
subplot(1,3,1)
scatter(price(1:l),trd_vol(1:l));
ylabel('Trading volume')
xlabel('Price')
title("Price-trd_vol Amazon")


%% Scatter plot for Volume between Amazon and Google
subplot(1,3,2)
scatter(trd_vol(1:l),trd_vol2(1:l)); 
ylabel('Trading volume Google')
xlabel('Trading volume Amazon')
title("Trd_vol Amazon-Google")

%% Scatter plot for Volatility between Amazon and Google
subplot(1,3,3)
scatter(vol0,vol2);
ylabel('Volatility Google')
xlabel('Volatility Amazon')
title("Volatility Amazon-Google")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION 19 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Linear Dependencies for log-returns of Amazon and Google
p=polyfit(r0,r2,1);
yfit = polyval(p,r0);
yresid = r2 - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(r2)-1) * var(r2);
rsq = 1 - SSresid/SStotal;

%% Non-linear Dependencies for log-returns of Amazon and Google
r0_sort=sort(r0);
r2_sort=sort(r2);
[RHO,PVAL] = corr(trd_vol(1:1000),trd_vol2(1:1000),'Type','Kendall');

%% Distribution of log-returns in different time scales
x0 = linspace(min(r0),max(r0),1000);
g0 = exp(-(x0-mean(r0)).^2/(2*std(r0)^2))/sqrt(2*pi*std(r0)^2);

r01=r0(1:360);
figure (12)
subplot(2,2,1)
histogram(r01,200,'Normalization','pdf')
set(gca,'YScale','log');
% xlim([-0.04,0.02]);
% ylim([1E-003,1E+003]);
hold on
xlabel('1 hour')
ylabel('Log-Returns')
% set(gca,'TickLabelInterpreter','LaTex')
% set(gca,'FontSize',20)
title('1 hour PDF-Amazon LogRet')

r02=r0(1:24*360);
subplot(2,2,2)
histogram(r02,200,'Normalization','pdf')
set(gca,'YScale','log');
% xlim([-0.04,0.02]);
% ylim([1E-003,1E+003]);
hold on
xlabel('1 day')
ylabel('Log-Returns')
% set(gca,'TickLabelInterpreter','LaTex')
% set(gca,'FontSize',20)
title('1 day PDF-Amazon LogRet')

r03=r0(1:7*24*360);
subplot(2,2,3)
histogram(r03,200,'Normalization','pdf')
set(gca,'YScale','log');
% xlim([-0.04,0.02]);
% ylim([1E-003,1E+003]);
hold on
xlabel('1 week')
ylabel('Log-Returns')
% set(gca,'TickLabelInterpreter','LaTex')
% set(gca,'FontSize',20)
title('1 week PDF-Amazon LogRet')

r04=r0(1:24*30*360);
subplot(2,2,4)
histogram(r04,200,'Normalization','pdf')
set(gca,'YScale','log');
% xlim([-0.04,0.02]);
% ylim([1E-003,1E+003]);
hold on
xlabel('1 month')
ylabel('Log-Returns')
% set(gca,'TickLabelInterpreter','LaTex')
% set(gca,'FontSize',20)
title('1 month PDF-Amazon LogRet')




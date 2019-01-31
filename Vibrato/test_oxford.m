% Script to test the convergence of the vibrato method, with browninan
% bridges over barriers. Using the function VMC_barr, written by Giles and
% Burgos.

UNDERLYING = [];
VALUE = [];
DELTA = [];
VEGA = [];
THETA = [];
RHO = [];

for s = 5:0.1:15
    [val,valD,valV,valR,valT,var,varD,varV,varR,varT] = VMC_barr(0.05,0.2,1,s,12,10,100000,10,10,0.7);
    UNDERLYING = [UNDERLYING s];
    VALUE = [VALUE val];
    DELTA = [DELTA valD];
    VEGA = [VEGA valV];
    THETA = [THETA valT];
    RHO = [RHO valR];
end
%%
figure
%plot(linspace(1,20,600),payoff(linspace(1,20,600),K),'k')
%hold on
subplot(3,2,[1,2])
plot(UNDERLYING,VALUE,'k')
grid
xlabel('Underlying')
ylabel('Value')
title('Price')
subplot(3,2,3)
plot(UNDERLYING,DELTA,'k')
grid
xlabel('Underlying')
ylabel('Value')
title('Delta')
subplot(3,2,4)
plot(UNDERLYING,VEGA,'k')
grid
xlabel('Underlying')
ylabel('Value')
title('Vega')
subplot(3,2,5)
plot(UNDERLYING,RHO,'k')
grid
xlabel('Underlying')
ylabel('Value')
title('Rho')
subplot(3,2,6)
plot(UNDERLYING,THETA,'k')
%legend('PAYOFF','PRICE','DELTA','VEGA','RHO','THETA')
grid
xlabel('Underlying')
ylabel('Value')
title('Theta')

% %% 
% figure
% plot(UNDERLYING,VALUE,'o--g')
% hold on
% plot(UNDERLYING,DELTA,'o--r')
% plot(UNDERLYING,VEGA,'o--m')
% plot(UNDERLYING,THETA,'o--c')
% plot(UNDERLYING,RHO,'o--b')
% legend('PRICE','DELTA','VEGA','THETA','RHO')
% xlabel('Underlying')
% ylabel('value')
% grid on

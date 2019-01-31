% Script for plotting of smoothed barriers

K = 10;
H = 14;
delta = 0.25;

x_values = [];
f_values = [];
derivatives = [];

for x = 1:0.01:20
    X = ADRev(x);
    f = adr_linbarrier_payoff(X,K,H,delta);
    rev = chainRule(f);
    x_values = [x_values x];
    f_values = [f_values f.value];
    derivatives = [derivatives X.derivative];
end

%figure
subplot(1,2,1)
plot(x_values,f_values,'b')
hold on
plot(x_values,derivatives,'--r')
grid
xlabel('x')
ylabel('Value')
title('Linear smoothing of Barrier')
legend('f(x)','df/dx')
axis([9 15 -7 4.2])


%%%%%
x_values = [];
f_values = [];
derivatives = [];
delta = 7;

for x = 1:0.01:20
    X = ADRev(x);
    f = adr_sigmoidbarrier_payoff(X,K,H,delta);
    rev = chainRule(f);
    x_values = [x_values x];
    f_values = [f_values f.value];
    derivatives = [derivatives X.derivative];
end

%figure
subplot(1,2,2)
plot(x_values,f_values,'b')
hold on
plot(x_values,derivatives,'--r')
grid
xlabel('x')
ylabel('Value')
title('Sigmoid smoothing of Barrier')
legend('f(x)','df/dx')
axis([9 15 -7 4.2])
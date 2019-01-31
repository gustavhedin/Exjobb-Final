% Test for testing of convergence of Reverse AD
clear all

MClimit = 100000;

smoothing_parameter = 10;
DELTA_mat = [];
RHO_mat = [];
VEGA_mat = [];
THETA_mat = [];
OPTIONPRICE_mat = [];

Q = 1000:1000:MClimit;

for q = Q

%UNDERLYING = [];
%OPTIONPRICE = [];
%DELTA = [];
   
    % inputs:
    % volatility sigma
    % interest rate r
    % maturity time tau
    start_der = 0;
    sigma_der = 0;
    r_der = 0;
    tau_der = 0;

    K = 10; % Strike price
    H = 14; % Barrrier

    % constants:
    N = 12; % # steps
    nbrMC = q; % # MC simulations
    X = zeros(N+1,1);

    for i = 1:nbrMC
        startvalue = ADRev(9); 
        sigma      = ADRev(0.2);
        r       =   ADRev(0.05);
        tau      =     ADRev(1);

        eRev_forward = geometric_brownian_adr(N,r,sigma,tau,startvalue,K,H,smoothing_parameter);
        reverse_sweep = chainRule(eRev_forward); 
        X = X + reverse_sweep.value(end);
        start_der = start_der + startvalue.derivative(end);
    end

    X = X/nbrMC;
    start_der = start_der/nbrMC;
    optionprice = X(end)*exp(-r.value*tau.value);

    DELTA = start_der*exp(-r.value*tau.value);    
    DELTA_mat = [DELTA_mat; DELTA];
    OPTIONPRICE_mat = [OPTIONPRICE_mat; optionprice(1)];

end
figure
plot(Q,DELTA_mat,'o--r')
%hold on
%plot(DELTA_mat(2,:),'o--g')
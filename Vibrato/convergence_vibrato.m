% Script to test the convergence of the vibrato method

MC_vec = [];
STRIKEPRICE = [];
DELTA = [];

PRICE = [];
x_mean = [];
MC = 1000:1000:100000;

for j = MC

    % User defined parameters:
    r = 0.05;           % risk free interest rate
    sigma = 0.2;        % volatility
    Startvalue = 9;     % Starting value for the underlying asset at time 0.
    T = 1;              % Time horizon
    N = 12;             % # simulation points on [0,T];
    K = 10;             % Strike price
    nbr_MC = j;    % # of Monte Carlo simulations
    nbrMC_z = 10;       % # of samples over the barrier


    % Declaring nedded variables:
    h = T/N;                            % Stepsize
    X = zeros(N,nbr_MC);                % Path matrix
    X(1,:) = Startvalue*ones(1,nbr_MC); % Initiate startvalue
    Y_delta = zeros(N,nbr_MC);          % Tangent process, Delta
    Y_delta(1,:) = ones(1,nbr_MC);      % Initiate Delta to be 1.
    Y_vega = zeros(N,nbr_MC);           % Tangent process, Vega
    Y_rho = zeros(N,nbr_MC);            % Tangent process, Rho
    Y_theta = zeros(N,nbr_MC);            % Tangent process, Theta

     for n=2:N
             Z = randn(1,nbr_MC);
            % Asset
            X(n,:) = X(n-1,:) + r*h*X(n-1,:) + sigma*sqrt(h)*X(n-1,:).*Z;

            % Tangent Processes: 
            % Delta: theta = X0:
            Y_delta(n,:) = Y_delta(n-1,:) + r*h*Y_delta(n-1,:) + 0*X(n-1,:) ...
                + (Y_delta(n-1,:)*sigma*sqrt(h)+ 0*X(n-1,:)).*Z; 

            % Vega: theta = sigma 
            Y_vega(n,:) = Y_vega(n-1,:) + r*h*Y_vega(n-1,:) + 0*X(n-1,:) ...
                + (Y_vega(n-1,:)*sigma*sqrt(h)+ sqrt(h)*X(n-1,:)).*Z; 

            % Rho: theta = sigma 
            Y_rho(n,:) = Y_rho(n-1,:) + r*h*Y_rho(n-1,:) + h*X(n-1,:) ...
                + (Y_rho(n-1,:)*sigma*sqrt(h)+ 0*X(n-1,:)).*Z; 
            
            % Theta: theta = sigma 
            Y_theta(n,:) = Y_theta(n-1,:) + r*h*Y_theta(n-1,:) + (r/N)*X(n-1,:) ...
                + (Y_theta(n-1,:)*sigma*sqrt(h)+ (sigma/sqrt(N*T))*X(n-1,:)).*Z; 
     end


     % Find price 
     Z = randn(1,nbr_MC);
     Xend = X(end,:)+r*h*X(end,:) + sigma*sqrt(h)*X(end,:).*Z;

     % Check asset at barrier checkpoints ?     
     % payoff = V(end)*prodsum(I1, I2 ...) , Ii = 1 if X > K else 0, at a specific barrier checkingpoint.  
     
     
     %Price = mean(max(Xend-K,0))*exp(-r*T);   

     % Generate random nubers for the step over the barrier.
     Z = randn(nbrMC_z,nbr_MC);

     % Obtain valued for the process and derivatives on the other side of
     % the barrier:

     firstpart = repmat(X(end,:)+r*h*X(end,:),nbrMC_z,1);
     lastpart = repmat(sigma*X(end,:)*sqrt(h),nbrMC_z,1);
     % First derivatives:
        X_Tplus = firstpart + Z.*lastpart;
        X_Tminus = firstpart - Z.*lastpart;
        X_Tdot = firstpart;

        V_Tplus = payoff(X_Tplus,K);
        V_Tminus = payoff(X_Tminus,K);
        V_Tdot = payoff(X_Tdot,K);
        V_now = (1/2)*(V_Tplus+V_Tminus);

        % Delta
         dmu_dtheta = repmat(Y_delta(end,:)*(1+r*h) + X(end,:)*0,nbrMC_z,1);
         dsig_dtheta = repmat(Y_delta(end,:)*sigma*sqrt(h) + X(end,:)*0,nbrMC_z,1);
         divfactor = repmat(1./(X(end,:)*sigma*sqrt(h)),nbrMC_z,1);
         Delta = mean(mean(dmu_dtheta.*(1/2).*(V_Tplus-V_Tminus).*(Z.*divfactor)...
                + dsig_dtheta.*(V_Tplus-2*V_Tdot+V_Tminus).*((Z.^2-1).*divfactor)))*exp(-r*T);

%         % Vega
%          dmu_dtheta =  repmat(Y_vega(end,:)*(1+r*h) + X(end,:)*0,nbrMC_z,1);
%          dsig_dtheta =  repmat(Y_vega(end,:)*sigma*sqrt(h) + X(end,:)*sqrt(h),nbrMC_z,1);
%          Vega = mean(mean(dmu_dtheta.*(1/2).*(V_Tplus-V_Tminus).*(Z.*divfactor)...
%                 + dsig_dtheta.*(V_Tplus-2*V_Tdot+V_Tminus).*((Z.^2-1).*divfactor)))*exp(-r*T);
% 
%         % Rho
%          dmu_dtheta = repmat(Y_rho(end,:)*(1+r*h) + X(end,:)*h,nbrMC_z,1);
%          dsig_dtheta =  repmat(Y_rho(end,:)*sigma*sqrt(h) + X(end,:)*0,nbrMC_z,1);
%          Rho = mean(mean(dmu_dtheta.*(1/2).*(V_Tplus-V_Tminus).*(Z.*divfactor)...
%                 + dsig_dtheta.*(V_Tplus-2*V_Tdot+V_Tminus).*((Z.^2-1).*divfactor)-T*V_now))*exp(-r*T);
% 
%         % Theta
%          dmu_dtheta =  repmat(Y_theta(end,:)*(1+r*h) + X(end,:)*0,nbrMC_z,1);
%          dsig_dtheta =  repmat(Y_theta(end,:)*sigma*sqrt(h) + X(end,:)*0 ,nbrMC_z,1);
%          Theta = mean(mean(dmu_dtheta.*(1/2).*(V_Tplus-V_Tminus).*(Z.*divfactor)...
%                 + dsig_dtheta.*(V_Tplus-2*V_Tdot+V_Tminus).*((Z.^2-1).*divfactor)-r*V_now))*exp(-r*T); 
        
        Price = mean(payoff(Xend,K))*exp(-r*T); 
    
    X = [X;Xend];
    
    MC_vec = [MC_vec j];
    DELTA = [DELTA Delta];
    PRICE = [PRICE Price];
    
    MEAN = mean(X');
    x_mean = [x_mean MEAN(end)];

end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    %plot(linspace(1,20,600),payoff(linspace(1,20,600),K),'k')
    %hold on
    plot(MC_vec,PRICE,'k')
    grid
    xlabel('Underlying')
    ylabel('Value')
    title('Price')
    figure
    plot(MC_vec,DELTA,'k')
    grid
    xlabel('Underlying')
    ylabel('Value')
    title('Delta')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
function P = payoff(x,K)
   %P=max(x-K,0);
   P = (x < 14).* max(x-K,0); % Strike K, Barrier at 12
end
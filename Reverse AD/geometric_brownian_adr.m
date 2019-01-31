function X_path = geometric_brownian_adr(N,r,sigma,T,startvalue,K,H,smoothing_parameter)

% Conceptual arlgorithm:
t = (0:1:N)'/N;               % t is the column vector [0 1/N 2/N ... 1]
W = [0; cumsum(randn(N,1))]/sqrt(N);%S is running sum of N(0,1/N) variables
t = t*T;
W = W*sqrt(T);
Y = (r-(sigma^2)/2)*t + sigma * W;
X = startvalue*exp(Y);

% Chose smoothing style of barrier:
%X_path = adr_linbarrier_payoff(X,K,H,smoothing_parameter);
X_path = adr_sigmoidbarrier_payoff(X,K,H,(1/smoothing_parameter));

end

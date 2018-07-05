function [ eta, u ] = matEvaluateExactSolution( obj, xb, zb, time )
%MATEVALUATEEXACTSOLUTION Summary of this function goes here
%   Detailed explanation goes here

I = sqrt( -1 );
tau0 = 3 * obj.miu0;
lambda = sqrt( I * obj.omega / obj.miu0 );
tau = tau0 / 3 * ( lambda ^ 2 * tanh( lambda ) / ...
    ( lambda + ( lambda^2 / obj.cD - 1 ) * tanh(lambda) ) );
beta2 = ( obj.omega ^ 2 - I * obj.omega * tau ) / obj.gra / obj.H0;
s1 = -0.5 + sqrt( 0.25 - beta2 );
s2 = -0.5 - sqrt( 0.25 - beta2 );

A = + obj.a*s2* obj.x1^s2 / (s2*obj.x2^s1 * obj.x1^s2 - s1*obj.x1^s1 * obj.x2^s2);
B = - obj.a*s1* obj.x1^s1 / (s2*obj.x2^s1 * obj.x1^s2 - s1*obj.x1^s1 * obj.x2^s2);

syms x z
eta = ( A * x .^ s1 + B * x .^ s2 );
u0 = - obj.gra / I / obj.omega * diff( eta, x );
u = u0 .* ( 1 - cosh(lambda .* z ) / cosh(lambda) / ( 1 + tanh(lambda)*lambda/obj.cD ) );

x = xb;
z = zb;

eta = eval( real( eta * exp( I * obj.omega * time ) ) );
u = eval( real( u * exp( I * obj.omega * time ) ) );
end


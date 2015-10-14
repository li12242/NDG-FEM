function [rhoP,rhouP,rhovP,EnerP] = ChannelBC2D(rhoP, rhouP, rhovP, EnerP, time)
  
% function [rho,rhou,rhov,Ener] = ChannelBC2D(rho, rhou, rhov, Ener, time)
% Purpose: Impose channel boundary conditions on 2D Euler equations on weak form

Globals2D;
gamma = 1.5; mu = 1e-2; pbar = 10;

xB = x(vmapP(mapB)); yB = y(vmapP(mapB));

% Quadratic shear flow, relies on gamma=1.5
rhoP(mapB)  = 1;
rhouP(mapB) = yB.^2;
rhovP(mapB) = 0;
EnerP(mapB) = (2*mu*xB + pbar)/(gamma-1) + .5*(yB.^4);
return

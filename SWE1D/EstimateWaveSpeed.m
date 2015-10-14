function [SM, SP] = EstimateWaveSpeed(mesh, hM, hP, uM, uP)
    gra = 9.8; hDelta = 10^-1;
    
    cs = 0.5*( sqrt(gra*hM) + sqrt(gra*hP) ) + 0.25*(uM - uP).*mesh.nx;
    us = ( 0.5*(uM + uP) ).*mesh.nx + (sqrt(gra*hM) - sqrt(gra*hP));

    us(cs<0) = 0; % middle is dry 
    cs(cs<0) = 0; % middle is dry
    hs = cs.^2/gra;
    
    SM = zeros(size(us)); SP = zeros(size(us));
    isWetM = hM>hDelta; isWetP = hP>hDelta;

    flag = ~isWetP & ~isWetM; % both is dry
    SM(flag) = 0;
    SP(flag) = 0;

    flag =  isWetP & ~isWetM; % left is dry
    SM(flag) = mesh.nx(flag) .* uP(flag) - 2*sqrt(gra.*hP(flag) ); % hM = 0
    flag = ~isWetP &  isWetM; % right is dry
    SP(flag) = mesh.nx(flag) .* uM(flag) + 2*sqrt(gra.*hM(flag) ); % hP = 0

    flag = isWetM & ( hs>hM ); % hs > hM >0
    SM(flag) = mesh.nx(flag).*uM(flag) - ...
        sqrt(gra*( hs(flag) + hM(flag) ).*hs(flag)./2./hM(flag) );

    flag = isWetP & ( hs>hP ); % hs > hP > 0
    SP(flag) = mesh.nx(flag).*uP(flag) + ...
        sqrt(gra*( hs(flag) + hP(flag) ).*hs(flag)./2./hP(flag) );

    flag = isWetM & ( hs<=hM );
    SM(flag) = mesh.nx(flag).*uM(flag) - sqrt(gra.*hM(flag));
    SM(flag) = min(SM(flag),...
        mesh.nx(flag).*us(flag) - sqrt(gra.*hs(flag)) );

    flag = isWetP & ( hs<=hP );
    SP(flag) = mesh.nx(flag).*uP(flag) + sqrt(gra.*hP(flag));
    SP(flag) = max(SP(flag),...
        mesh.nx(flag).*us(flag) + sqrt(gra.*hs(flag)) );
end% func
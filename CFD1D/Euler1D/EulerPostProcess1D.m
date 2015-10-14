figure('color', 'w')

subplot(2,2,1)
plot(x(:), rho(:))
xlabel('rho')

subplot(2,2,2)
plot(x(:), rhou(:)./rho(:))
xlabel('rho*u')
set(gca, 'Ylim', [-0.2, 1])

subplot(2,2,3)
gamma = 1.4;
press = (gamma - 1).*(Ener - 0.5*(rhou).^2./rho);
plot(x(:), press(:))
xlabel('pressure')

subplot(2,2,4)
cvel = sqrt(gamma*press./rho);
mach = (rhou./rho)./cvel;
plot(x(:), mach(:))
xlabel('Mach number')
set(gca, 'Ylim', [-0.2, 1])
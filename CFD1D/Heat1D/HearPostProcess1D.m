y = exp(-FinalTime)*sin(x);

figure('color', 'w')
subplot(2,1,1)
hold on
plot(x(:), u(:))
plot(x(:), y(:), 'o')

error = abs(y-u);
subplot(2,1,2)
plot(x(:), error(:))
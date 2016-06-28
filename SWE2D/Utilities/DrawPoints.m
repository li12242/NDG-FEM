function DrawPoints(mesh, h, qx, qy)

subplot(1, 3, 1)
plot3(mesh.x, mesh.y, h, 'b.'); hold on;
% plot(mesh.x(mesh.vmapM), mesh.y(mesh.vmapM));

subplot(1, 3, 2)
plot3(mesh.x, mesh.y, qx, 'r.');

subplot(1, 3, 3)
plot3(mesh.x, mesh.y, qy, 'r.');
drawnow;
end% func
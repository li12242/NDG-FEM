function DrawPoints(mesh, h, qx, qy)

subplot(2,2, [1,3])
plot3(mesh.x, mesh.y, h, 'b.'); hold on;
% plot(mesh.x(mesh.vmapM), mesh.y(mesh.vmapM));

subplot(2,2,2)
plot3(mesh.x, mesh.y, qx, 'r.');
subplot(2,2,4)
plot3(mesh.x, mesh.y, qy, 'r.');

drawnow;
end% func
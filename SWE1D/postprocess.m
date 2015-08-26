pix_left = 126; pix_right = 706;
pix_top = 72; pix_bottom = 601;
ratio.x = 1000/(pix_right - pix_left); ratio.y = 10./(pix_bottom-pix_top);
P(1,:) = [126-pix_left, pix_bottom-72].*[ratio.x, ratio.y];
P(2,:) = [240-pix_left, pix_bottom-72].*[ratio.x, ratio.y];
P(3,:) = [387-pix_left, pix_bottom-332].*[ratio.x, ratio.y];
P(4,:) = [578-pix_left, pix_bottom-332].*[ratio.x, ratio.y];
P(5,:) = [578-pix_left, pix_bottom-493].*[ratio.x, ratio.y];
P(6,:) = [706-pix_left, pix_bottom-493].*[ratio.x, ratio.y];

figure('Color', 'w')
hold on
plot(P(:,1), P(:,2), 'ko-');
% D = Q(:,:,1);
plot(mesh.x, Q(:,:,1), 'r'); 
set(gca, 'YLim', [0,11]);
legend('Theoretial','NDG')
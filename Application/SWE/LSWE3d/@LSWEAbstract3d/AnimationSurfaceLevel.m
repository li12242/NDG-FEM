function AnimationSurfaceLevel( obj )
%ANIMATIONSURFACELEVEL Summary of this function goes here
%   Detailed explanation goes here

visual = Visual2d( obj.mesh2d );
OutstepNum = obj.outputFile.outputStep;

video = VideoWriter( [obj.outputFile.casename, '/', ...
    obj.outputFile.casename, '.avi'] );
video.FrameRate = 30;
open( video );

% initialize axis with the first step output
figure('Color', 'w');
[ fphys2d, ~ ] = obj.outputFile.readOutputResult(1);
visual.drawResult( fphys2d(:, :, 1) );
zlim([-1.1, 1.1]);
zlabel('$\xi$ (m)', 'Interpreter', 'Latex', 'FontSize', 14);
xlabel('$x$ (m)', 'Interpreter', 'Latex', 'FontSize', 14);
ylabel('$y$ (m)', 'Interpreter', 'Latex', 'FontSize', 14);

for n = 1 : 2 : OutstepNum
    [ fphys2d, ~ ] = obj.outputFile.readOutputResult(n);
    eta = fphys2d(:, :, 1);
    
    visual.drawResult( eta );
    frame = getframe( gcf );
    writeVideo( video, frame );
end

close( video );

end


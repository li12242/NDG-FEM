function showParameter( obj )

w = 2 * pi / obj.T;
c = sqrt( obj.gra * obj.H );
k = w/c;
waveLen = obj.T * sqrt(obj.gra * obj.H);
fprintf('\nWave parameter:\n');
fprintf('\tWave amplitude (m): %f\n', obj.delta * obj.H);
fprintf('\tWave period (s): %f\n', obj.T);
fprintf('\tWave speed (m/s): %f\n', sqrt(obj.gra * obj.H ));
fprintf('\tWave length (m): %f\n', waveLen);
fprintf('\tWave k/H: %e\n', k/obj.H);

fprintf('\nChannel parameter\n');
fprintf('\tChannel depth (m): %f\n', obj.H);
fprintf('\tChannle length (m): %f\n', obj.ChLength);
fprintf('\tNumber of wave: %f\n', obj.ChLength / waveLen);
fprintf('\tTime to propage through the channel (h): %f\n', ...
    obj.ChLength/sqrt(obj.gra * obj.H ) / 3600 );

Ne = obj.meshUnion.K/2;
fprintf('\nModel parameter\n');
fprintf('\tElement number: %d\n', Ne );
fprintf('\tNumber of elements in each wave: %f\n',  Ne/(obj.ChLength / waveLen));
fprintf('\tSimulate time (h) = %f\n', obj.getOption('finalTime')/3600 );
fprintf('\tNumber of waves pass through the channel: %f\n', ...
    obj.getOption('finalTime')/( obj.ChLength/sqrt(obj.gra * obj.H ) ) );
fprintf('\tNumber of wave periods: %f\n', obj.getOption('finalTime')/obj.T);
fprintf('\tOutput interval: %d\n', obj.getOption('outputTimeInterval'));

end


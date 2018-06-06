%> \brief check scalar fields to evaluate the OBCs performance
%> \details The scalar fields include
%> * mechanical energy
%> * excess mass
%> * energy flux
%>
function checkResult( obj )

chpos = makeNdgPostProcessFromNdgPhys( obj );
bot = obj.fphys{1}(:, :, 4);

% mechnical energy
mechEng = zeros( chpos.Nt, 1 );
mechEngExt = zeros( chpos.Nt, 1 );

% excess mass
mass = zeros( chpos.Nt, 1 );
massExt = zeros( chpos.Nt, 1 );

% energy flux
engFlux = zeros( chpos.Nt, 1 );
engFluxExt = zeros( chpos.Nt, 1 );

line = StdLine( obj.meshUnion(1).cell.N );
for t = 1:chpos.Nt
    fphys = chpos.accessOutputResultAtStepNum(t);
    time = chpos.time{1}(t);
    
    for m = 1:obj.Nmesh
        mesh = obj.meshUnion(m);
        [h, hu] = obj.setExactSolution( mesh.x, time );
        
        % mechanical energy
        mechEng(t) = sum ( mesh.GetMeshAverageValue( ...
            obj.gra * (fphys{m}(:, :, 1) - bot).^2 ...
            + 0.5*( fphys{m}(:, :, 2).^2 + fphys{m}(:, :, 3).^2 )./fphys{m}(:, :, 1).^2 ) );
        mechEngExt(t) = sum ( ...
            mesh.GetMeshAverageValue( obj.gra * ( h - bot ).^2 ...
            + 0.5*( hu.^2 ) ./ h.^2 ) );
        
        % excess mass
        mass(t) = mean( mesh.GetMeshAverageValue( fphys{1}(:,:,1) ) );
        massExt(t) = mean( mesh.GetMeshAverageValue( h ) );
        
        % energy flux
        boundaryNodeId = mesh.cell.Fmask(:, 4);
        H = fphys{m}( boundaryNodeId, 1, 1 );
        eta = H - obj.H;
        U = fphys{m}( boundaryNodeId, 1, 2 ) ./ H;
        V = fphys{m}( boundaryNodeId, 1, 3 ) ./ H;
        temp = obj.H .* U .* ( obj.gra .* eta + 0.5 * V.^2 );
        engFlux(t) = sum( line.wq .* temp .* 0.5 );
        He = h( boundaryNodeId, 1 );
        Ue = hu( boundaryNodeId, 1 ) ./ He;
        eta = He - obj.H;
        engFluxExt(t) = sum( line.wq .* obj.H .* Ue .* ( obj.gra .* eta ) );
    end
end
figure('color', 'w', 'Position', [100, 170, 625, 635]);
subplot(3, 1, 1); hold on; grid on; box on;
plot( chpos.time{1}, mechEngExt, 'r', 'LineWidth', 2 );
plot( chpos.time{1}, mechEng, 'b', 'LineWidth', 2 );
xlabel('Time (s)', 'Interpreter', 'latex','FontSize', 16);
ylabel('Energy (J/$\mathrm{m}^2$)', 'Interpreter', 'latex','FontSize', 16);

subplot( 3, 1, 2 ); hold on; grid on; box on;
plot( chpos.time{1}, massExt, 'r', 'LineWidth', 2 );
plot( chpos.time{1}, mass, 'b', 'LineWidth', 2 );
xlabel('Time (s)', 'Interpreter', 'latex','FontSize', 16);
ylabel('Mean sea level (m)', 'Interpreter', 'latex','FontSize', 16);

subplot( 3, 1, 3 ); hold on; grid on; box on;
plot( chpos.time{1}, engFluxExt, 'r', 'LineWidth', 2 );
plot( chpos.time{1}, engFlux, 'b', 'LineWidth', 2 );
xlabel('Time (s)', 'Interpreter', 'latex','FontSize', 16);
ylabel('Energy flux (m)', 'Interpreter', 'latex','FontSize', 16);

title('Flather boundary for closed channed', 'Interpreter', 'latex','FontSize', 16);
end


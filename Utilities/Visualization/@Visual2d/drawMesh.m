function drawMesh(obj)
    mesh = obj.mesh;
    patch('Vertices', [mesh.vx(:), mesh.vy(:)], ...
        'Faces', mesh.EToV', 'FaceColor', [0.8, 0.9, 1]);
    box on;
    grid on;
end
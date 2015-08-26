% Err = [ 50  0.008242     0.044100   0.360795  
%  100  0.001851     0.011108   0.104890  
%  200  0.001121     0.009720   0.139260  ];

% rarefaction
Err = [50  0.005601 	 0.012645 	 0.046397 	
 100 0.002152 	 0.005301 	 0.021544 	
 200 0.001238 	 0.003462 	 0.019960 	];


L1 = Err(:,2);
L2 = Err(:,3);
Linf = Err(:,4);

mesh = Err(:,1); deltaX = 1./mesh;

ErrStr = {'L1', 'L2', 'Linf'};
for i = 1:3
    subplot(3,1,i)
    plot(log10(mesh), log10(Err(:, i+1)), 'ro-')
    xlabel(ErrStr{i})
end


ErrMius = Err(1:end-1, 2:4);
ErrPlus = Err(2:end, 2:4);

rate = log10(ErrMius./ErrPlus)./log10(2);

fprintf('|mesh |L1 \t| rate \t| L2 \t| rate \t| Linf \t| rate \t|\n')
fprintf('| | %f \t| ~ \t| %f \t| ~ \t| %f \t| ~ \t|\n', L1(1), ...
     L2(1), Linf(1))
for i = 1:(size(Err, 1) - 1)
    fprintf('| | %f \t| %f \t| %f \t| %f \t| %f \t| %f \t|\n', L1(i+1), ...
        rate(i,1), L2(i+1), rate(i,2), Linf(i+1), rate(i,3))
end

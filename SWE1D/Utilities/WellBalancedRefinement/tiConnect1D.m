function [EToE,EToF]= tiConnect1D(mesh, EToV)
% re-connect the element and faces
% 1D line element
Nfaces=2;
K = mesh.nElement;

fnodes = [EToV(:, 1); EToV(:, 2)];
EToE= (1:K)'*ones(1,Nfaces); EToF= ones(K,1)*(1:Nfaces);
id = fnodes(:);

spNodeToNode = [id, (1:Nfaces*K)', EToE(:), EToF(:)];
% Now we sort by global face number.
sorted=sortrows(spNodeToNode,1);

[indices,~]=find( sorted(1:(end-1),1)==sorted(2:end,1) );

% make links reflexive 
matchL = [sorted(indices,:)   ;sorted(indices+1,:)];
matchR = [sorted(indices+1,:) ;sorted(indices,:)];

% insert matches
EToE(matchL(:,2)) = matchR(:,3); EToF(matchL(:,2)) = matchR(:,4);
end% func
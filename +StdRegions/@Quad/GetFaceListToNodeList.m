function [ nodelist ] = GetFaceListToNodeList( obj )
%GETFACELISTTONODELIST Summary of this function goes here
%   Detailed explanation goes here

Np = obj.nOrder+1; % No of points at one face
nodelist = zeros(obj.nFaceNode, 1);

% node arrangement matrix
indMatrix = flip(reshape(1:Np^2, Np, Np)');

% face I, bottom
nodelist(1:Np) = indMatrix(Np, :);
% face II, right
nodelist(Np+1:2*Np) = flip(indMatrix(:, Np));
% face III, up
nodelist(2*Np+1:3*Np) = flip(indMatrix(1, :));
% face IV, left
nodelist(3*Np+1:4*Np) = indMatrix(:, 1);
end


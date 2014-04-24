function gphs = newGphUs(Pts, parGph)
% Connect nodes to generate edges.
%
% Remark
%   The edge is undirected and the edge feature is symmetric.
%
% Input
%   Pts     -  graph node, 1 x mG (cell), d x ni
%   parGph  -  parameter for computing the adjacency matrix
%              see gphEg.m for more details
%
% Output
%   gphs    -  graph, 1 x mG (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-24-2014

% dimension
mG = length(Pts);

% per graph
gphs = cell(1, mG);
for iG = 1 : mG
    gphs{iG} = newGphU(Pts{iG}, parGph);
end

function cellcomplex_annotate_nodes(S, col)
% CELLCOMPLEX_ANNOTATE_NODES
% Overlay node indices on current figure.
%
% INPUT
%   S.nodes : [K x 2] node coordinates
%   col:    color of node indexing
%
% The function writes node numbers at their spatial locations
% on the current axes. No new figure is created.
%
% (C) 2026 Moo K. Chung

if ~isstruct(S) || ~isfield(S,'nodes') || isempty(S.nodes)
    return
end

Nodes = S.nodes;
K     = size(Nodes,1);

hold on;

% small vertical offset for readability
offset = 0.028 * max(range(Nodes));  

for k = 1:K
    text(Nodes(k,1), Nodes(k,2) + offset, ...
        num2str(k), ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment','center', ...
        'Color',col);
end

hold off;

end
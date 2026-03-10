function S = cellcomplex_add_faces_random(S, p, seed)
% CELLCOMPLEX_ADD_FACES_RANDOM
% Randomly select faces on the existing rectangular grid (each with prob p),
% MERGE connected selected faces into larger regions, DELETE all interior
% shared edges, and rebuild only coherent boundary cycles as directed edges.
%
% Output is NOT a checkerboard. Orientation is assigned per connected
% component: each component boundary is randomly CW or CCW.
%
% NOTE: Here S.faces is NOT a true list of 2-cells of a cubical complex.
% It is metadata used to group and color boundary cycles associated with each
% merged region component. The topological object returned is the boundary
% 1-skeleton in S.edges and S.cycles.
%
% INPUT:
%   S      : struct with field
%            .nodes : [K x 2] grid coordinates on a rectangular lattice
%            Optional existing fields .edges, .flows, .cycles, .faces are overwritten.
%   p      : Bernoulli probability for selecting each unit square face
%   seed   : optional random seed
%
% OUTPUT:
%   S.nodes  : unchanged input node coordinates
%   S.edges  : [E x 2] directed boundary edges only
%   S.flows  : [E x 1] positive unit flow along stored edge direction
%   S.cycles : cell array; each entry is one directed boundary cycle
%   S.faces  : struct describing merged regions:
%              .cycles : same as S.cycles
%              .sign   : [Nc x 1] +1 CW, -1 CCW
%              .mask   : [n x m] face selection mask
%
% (C) 2026 Moo K. Chung

if nargin < 3 || isempty(seed)
    seed = 1;
end
rng(seed);

nodes = S.nodes;

% ---------------------------------------------------------------------
% infer rectangular grid size from node coordinates                     % <— FIX
% ---------------------------------------------------------------------
xv = unique(nodes(:,1));
yv = unique(nodes(:,2));

nColV = numel(xv);
nRowV = numel(yv);

m = nColV - 1;
n = nRowV - 1;

idxV = @(r,c) (r-1)*nColV + c;

% ---------------------------------------------------------------------
% face mask: select faces independently                                 % <— FIX
% ---------------------------------------------------------------------
mask = (rand(n,m) < p);

% guarantee at least one face
if ~any(mask(:))
    mask(randi(n), randi(m)) = true;
end

% ---------------------------------------------------------------------
% connected components (4-neighborhood)                                 % <— same as build_fractal
% ---------------------------------------------------------------------
visited = false(n,m);

cycles = {};
signs  = [];
allEdges = zeros(0,2);

for r0 = 1:n
    for c0 = 1:m

        if ~mask(r0,c0) || visited(r0,c0)
            continue
        end

        % BFS to collect one component
        Q = [r0 c0];
        visited(r0,c0) = true;
        comp = zeros(0,2);

        while ~isempty(Q)
            r = Q(1,1);
            c = Q(1,2);
            Q(1,:) = [];

            comp(end+1,:) = [r c]; %#ok<AGROW>

            % 4-neighbors
            nb = [r-1 c; r+1 c; r c-1; r c+1];
            for t = 1:4
                rr = nb(t,1);
                cc = nb(t,2);

                if rr>=1 && rr<=n && cc>=1 && cc<=m
                    if mask(rr,cc) && ~visited(rr,cc)
                        visited(rr,cc) = true;
                        Q(end+1,:) = [rr cc]; %#ok<AGROW>
                    end
                end
            end
        end

        % -----------------------------------------------------------------
        % boundary segments in CW orientation around region                 % <— same as build_fractal
        % -----------------------------------------------------------------
        seg = zeros(0,2);

        for t = 1:size(comp,1)
            r = comp(t,1);
            c = comp(t,2);

            TL = idxV(r,   c);
            TR = idxV(r,   c+1);
            BL = idxV(r+1, c);
            BR = idxV(r+1, c+1);

            % If neighbor NOT selected, that side is on the boundary.
            % We orient boundary consistently CW around the component:
            % top:    TL->TR
            % right:  TR->BR
            % bottom: BR->BL
            % left:   BL->TL
            if r==1 || ~mask(r-1,c)
                seg(end+1,:) = [TL TR]; %#ok<AGROW>
            end
            if c==m || ~mask(r,c+1)
                seg(end+1,:) = [TR BR]; %#ok<AGROW>
            end
            if r==n || ~mask(r+1,c)
                seg(end+1,:) = [BR BL]; %#ok<AGROW>
            end
            if c==1 || ~mask(r,c-1)
                seg(end+1,:) = [BL TL]; %#ok<AGROW>
            end
        end

        if isempty(seg)
            continue
        end

        % -----------------------------------------------------------------
        % choose random orientation per component                           % <— same as build_fractal
        % -----------------------------------------------------------------
        if rand < 0.5
            sgn = 1;                      % keep CW as built
        else
            sgn = -1;                     % flip to CCW
            seg = seg(:,[2 1]);           % reverse every directed segment
        end

        % -----------------------------------------------------------------
        % order segments into cycles (touching / multiple cycles allowed)   % <— same as build_fractal
        % -----------------------------------------------------------------
        used = false(size(seg,1),1);
        cycList = {};

        while any(~used)
            i0 = find(~used, 1, 'first');
            e0 = seg(i0,:);
            used(i0) = true;

            vstart = e0(1);
            vcur   = e0(2);

            cyc = e0;

            while true
                cand = find(~used & seg(:,1)==vcur, 1, 'first');

                if isempty(cand)
                    break
                end

                enext = seg(cand,:);
                used(cand) = true;

                cyc(end+1,:) = enext; %#ok<AGROW>
                vcur = enext(2);

                if vcur == vstart
                    break
                end
            end

            cycList{end+1} = cyc; %#ok<AGROW>
        end

        % store ALL cycles found for this component
        for q = 1:numel(cycList)
            cycles{end+1} = cycList{q}; %#ok<AGROW>
            signs(end+1,1) = sgn; %#ok<AGROW>
            allEdges = [allEdges; cycList{q}]; %#ok<AGROW>
        end
    end
end

% ---------------------------------------------------------------------
% build S.edges as union of boundary directed edges only                 % <— same as build_fractal
% ---------------------------------------------------------------------
if isempty(allEdges)
    S.edges  = zeros(0,2);
    S.flows  = zeros(0,1);
    S.cycles = {};
    S.faces  = struct('cycles', {{}}, 'sign', [], 'mask', mask);
    return
end

[edgesU, ~] = unique(allEdges, 'rows', 'stable');

S.edges  = edgesU;
S.flows  = ones(size(edgesU,1),1);
S.cycles = cycles;

S.faces  = struct('cycles', {cycles}, 'sign', signs, 'mask', mask);

end
function V = validate_suite(folder, varargin)
%VALIDATE_SUITE  Independent check of a suite folder and per-instance INDEX.csv.
%
%   V = validate_suite('core/D')
%   V = validate_suite(folder, 'n_random', 200, 'seed', 1, 'index_file', 'INDEX.csv')
%
% For every wrapper of the folder (version 2 fields required):
%   1. loads, has the version-2 fields (cardinality, func_idx, idx2cell,
%      cell2idx, f_star, domain_size); the optimal solution comes from the
%      folder's reference_optima.csv (attach_reference_optimum) when the
%      wrapper does not carry it, as in the core;
%   2. f = f_star at EVERY point of x_star_all (relative tolerance 1e-9) and
%      x_star_all_idx agrees with cell2idx(x_star_all);
%   3. x* is a local minimizer of the finite landscape: no neighbour has a
%      lower value, where the neighbours of x* are, for every ordered
%      discrete coordinate the levels j+-1 and for every categorical
%      coordinate ALL the other categories (categories have no order), the
%      other coordinates fixed — a failure here means the documented
%      minimizer is rounded and a grid point is better;
%   3b. certification by enumeration (fully finite domains only, i.e. no
%      continuous variable, and domain_size <= 'enumerate_max', default 2e6):
%      f_cert = min over the WHOLE domain of func_idx, computed blockwise;
%      x0 is then compared with f_cert, so a trivial instance (column
%      x0_is_optimum) is detected even without a documented minimizer;
%      the number of global minimizers of the finite problem is counted;
%      if f_cert < f_star the wrapper's f_star is NOT the optimum (flagged);
%      larger or continuous domains stay explicitly 'not certified';
%   4. on n_random Halton-distributed index points: func_idx == func_F(idx2cell),
%      cell2idx(idx2cell(z)) == z, two evaluations give the same value;
%   5. fraction of non-finite values among those points and the 2n+1 points
%      of the reference initial list (centre and axial vertices);
%   6. the protocol's starting point x0 (PB.x0_idx; the centre when absent) is not an optimum
%      (problem), unless the wrapper declares it (PB.x0_is_optimum, possible
%      in the controlled-cardinality family: column x0_is_optimum, the
%      instance is trivial and data_profiles excludes it by default); the
%      2n axial vertices are only counted (column
%      axial_initial_optima: how many of them are optima — relevant to a
%      solver that starts from the centre + axial vertices, as some
%      direct-search methods do; unavoidable with 2-4 levels, where every index is
%      extreme or central);
%   7. mean evaluation time.
% Writes <folder>/INDEX.csv (one row per instance) and prints a summary.
opts = struct('n_random', 200, 'seed', 1, 'index_file', 'INDEX.csv', 'enumerate_max', 2e6, 'enumerate', true);
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
addpath(folder);
L = dir(fullfile(folder, '*.m'));
L = L(~ismember({L.name}, {'add_problems.m', 'add_strategy.m', 'add_all_strategies.m'}));   % every wrapper of the folder
names = sort(cellfun(@(s) s(1:end-2), {L.name}, 'UniformOutput', false));
fields = {'name','nvar','var_types','ProblemData','func_F','f_star','x_star','x_star_all', ...
          'cardinality','func_idx','idx2cell','cell2idx','x_star_idx','x_star_all_idx','f_star_documented','domain_size'};
% random index points: Halton points (deterministic, portable; opts.seed is the index offset)
V = struct([]); nfail = 0;
fid = fopen(fullfile(folder, opts.index_file), 'w');
fprintf(fid, 'instance,base_function,n,pattern,cardinalities,domain_size,f_star,f_star_documented,n_minimizers,n_minimizers_documented,x_star_local_min,certified,f_star_certified,n_certified_minimizers,f_x0,x0_is_optimum,axial_initial_optima,nonfinite_fraction,eval_us,problems\n');
for k = 1:numel(names)
    nm = names{k}; prob = {};
    try
        PB = load_quiet(nm);
    catch ME
        fprintf('FAIL %s: does not load (%s)\n', nm, ME.message); nfail = nfail + 1; continue;
    end
    PB = attach_reference_optimum(PB, folder);           % core wrappers keep x* in reference_optima.csv
    miss = fields(~isfield(PB, fields));
    if ~isempty(miss), prob{end+1} = ['missing fields: ' strjoin(miss, ' ')]; end %#ok<AGROW>
    n = PB.nvar; card = PB.cardinality; known = isfinite(PB.f_star);
    % 2. optima
    nmin = 0; localmin = NaN;
    if known && isfield(PB, 'x_star_all')
        nmin = numel(PB.x_star_all);
        for q = 1:nmin
            fq = PB.func_F(PB.x_star_all{q});
            if ~(abs(fq - PB.f_star) <= 1e-9 * max(1, abs(PB.f_star))), prob{end+1} = sprintf('f(x*_%d)=%.10g ~= f*=%.10g', q, fq, PB.f_star); end %#ok<AGROW>
            if isfield(PB, 'x_star_all_idx') && ~isequal(PB.cell2idx(PB.x_star_all{q}), PB.x_star_all_idx{q}), prob{end+1} = sprintf('x_star_all_idx{%d} inconsistent', q); end %#ok<AGROW>
        end
        % 3. local minimality in index space
        localmin = true;
        if isfield(PB, 'x_star_idx')
            z = PB.x_star_idx;
            for i = 1:n
                if isinf(card(i)), continue; end
                if strcmp(PB.var_types{i}, 'K'), alt = setdiff(1:card(i), z(i));   % every other category
                else, alt = z(i) + [-1 1]; alt = alt(alt >= 1 & alt <= card(i)); end  % ordered: the two adjacent levels
                for a = alt
                    zz = z; zz(i) = a;
                    if PB.func_idx(zz) < PB.f_star - 1e-12 * max(1, abs(PB.f_star)), localmin = false; end
                end
            end
            if ~localmin, prob{end+1} = 'x* is not a local minimizer in index space (rounded minimizer?)'; end %#ok<AGROW>
        end
    end
    % 3b. certification by enumeration of fully finite domains
    certified = false; fcert = NaN; ncert = NaN; x0opt_cert = false;
    if opts.enumerate && isfield(PB, 'func_idx') && all(isfinite(card)) && prod(card) <= opts.enumerate_max
        [fcert, ncert] = enumerate_domain(PB, card);
        certified = true;
        if known && fcert < PB.f_star - 1e-9 * max(1, abs(PB.f_star))
            prob{end+1} = sprintf('enumeration finds f = %.12g < f_star = %.12g: f_star is not the optimum', fcert, PB.f_star); %#ok<AGROW>
        end
        % trivial instance detected from the certified optimum (also when no minimizer is documented)
        fx0c = PB.func_idx(P_x0_idx(PB, card));
        if isfinite(fcert) && abs(fx0c - fcert) <= 1e-9 * max(1, abs(fcert)), x0opt_cert = true; end
    end
    if x0opt_cert
        x0opt = true;
        if known && ~(isfield(PB, 'x0_is_optimum') && PB.x0_is_optimum)
            prob{end+1} = 'x0 (centre) attains the certified optimum but PB.x0_is_optimum is false'; %#ok<AGROW>
        end
    end
    % 4-5. random points
    nnf = 0; ntot = 0; tsum = 0;
    if isfield(PB, 'func_idx')
        U = halton_point(opts.seed * 1000 + (1:opts.n_random), n);
        for t = 1:opts.n_random
            u = U(t, :); z = zeros(1, n);
            for i = 1:n
                if isinf(card(i)), z(i) = PB.lb_C_work(i) + u(i) * (PB.ub_C_work(i) - PB.lb_C_work(i));
                else, z(i) = 1 + floor(u(i) * card(i)); end
            end
            t0 = tic; f1 = PB.func_idx(z); tsum = tsum + toc(t0); ntot = ntot + 1;
            f2 = PB.func_F(PB.idx2cell(z)); f3 = PB.func_idx(z);
            if ~(isequaln(f1, f2) || abs(f1 - f2) <= 1e-12 * max(1, abs(f1))), prob{end+1} = 'func_idx ~= func_F'; break; end %#ok<AGROW>
            if ~isequaln(f1, f3), prob{end+1} = 'not deterministic'; break; end %#ok<AGROW>
            z2 = PB.cell2idx(PB.idx2cell(z));
            if any(abs(z2 - z) > 1e-12 * max(1, abs(z))), prob{end+1} = 'cell2idx(idx2cell(z)) ~= z'; break; end %#ok<AGROW>
            if ~isfinite(f1), nnf = nnf + 1; end
        end
    end
    % 5-6. reference initial list: centre and axial vertices
    fx0 = NaN; n_axial = 0; x0opt = false;
    if isfield(PB, 'x0_is_optimum') && PB.x0_is_optimum && known
        % declared trivial: check that the declaration is true
        if abs(PB.func_idx(P_x0_idx(PB, card)) - PB.f_star) > 1e-9 * max(1, abs(PB.f_star)), prob{end+1} = 'PB.x0_is_optimum is true but f(x0) ~= f_star'; end %#ok<AGROW>
    end
    if isfield(PB, 'func_idx')
        Z = [0.5*ones(n,1), eye(n), 1-eye(n)];
        for c = 1:size(Z, 2)
            z = zeros(1, n);
            for i = 1:n
                if isinf(card(i)), z(i) = PB.lb_C_work(i) + Z(i,c) * (PB.ub_C_work(i) - PB.lb_C_work(i));
                else, z(i) = round(Z(i,c) * (card(i) - 1)) + 1; end
            end
            if c == 1, z = P_x0_idx(PB, card); end            % the protocol's x0 (Halton in the core, centre otherwise)
            f = PB.func_idx(z); ntot = ntot + 1; if ~isfinite(f), nnf = nnf + 1; end
            if c == 1, fx0 = f; end
            if known && abs(f - PB.f_star) <= 1e-9 * max(1, abs(PB.f_star))
                if c == 1
                    x0opt = true;                                 % trivial for the protocol: recorded, never a failure
                    if isfield(PB, 'x0_is_optimum') && ~PB.x0_is_optimum && ~isfield(PB, 'core')
                        prob{end+1} = 'x0 is an optimum but PB.x0_is_optimum is false'; %#ok<AGROW>
                    end
                else
                    n_axial = n_axial + 1;
                end
            end
        end
    end
    pattern = ''; if isfield(PB, 'pattern'), pattern = PB.pattern; end
    fdoc = NaN; if isfield(PB, 'f_star_documented'), fdoc = PB.f_star_documented; end
    ds = NaN; if isfield(PB, 'domain_size'), ds = PB.domain_size; end
    cstr = strjoin(arrayfun(@(c) num2str(c), card, 'UniformOutput', false), ';');
    ndoc = nmin; if isfield(PB, 'n_minimizers_documented'), ndoc = PB.n_minimizers_documented; end
    basef = ''; if isfield(PB, 'base_function'), basef = PB.base_function; end
    fprintf(fid, '%s,%s,%d,%s,%s,%.16g,%.16g,%.16g,%d,%d,%d,%d,%.16g,%g,%.16g,%d,%d,%.4f,%.1f,%s\n', nm, basef, n, pattern, cstr, ds, PB.f_star, fdoc, nmin, ndoc, localmin, certified, fcert, ncert, fx0, x0opt, n_axial, nnf/max(1,ntot), 1e6*tsum/max(1,ntot), strjoin(prob, ' | '));
    V(end+1).instance = nm; V(end).n = n; V(end).f_star = PB.f_star; V(end).n_minimizers = nmin; %#ok<AGROW>
    V(end).certified = certified; V(end).f_star_certified = fcert; V(end).n_certified_minimizers = ncert; V(end).axial_initial_optima = n_axial; V(end).x0_is_optimum = x0opt;
    V(end).x_star_local_min = localmin; V(end).nonfinite_fraction = nnf/max(1,ntot); V(end).problems = prob;
    if ~isempty(prob), nfail = nfail + 1; fprintf('PROBLEM %s: %s\n', nm, strjoin(prob, ' | ')); end
end
fclose(fid);
rmpath(folder);
known = sum(isfinite([V.f_star]));
fprintf('\n%s: %d instances, %d with known optimum (%d with several minimizers), %d certified by enumeration, %d with problems; instances with non-finite values: %d; trivial (x0 is an optimum): %d; with an optimum at an axial initial vertex: %d; INDEX written to %s\n', ...
    folder, numel(V), known, sum([V.n_minimizers] > 1), sum([V.certified]), nfail, sum([V.nonfinite_fraction] > 0), sum([V.x0_is_optimum]), sum([V.axial_initial_optima] > 0), fullfile(folder, opts.index_file));
end

function z = P_x0_idx(PB, card)
if isfield(PB, 'x0_idx'), z = PB.x0_idx(:)'; return; end
n = PB.nvar; z = zeros(1, n);
for i = 1:n
    if isinf(card(i)), z(i) = 0.5 * (PB.lb_C_work(i) + PB.ub_C_work(i)); else, z(i) = round(0.5*(card(i)-1)) + 1; end
end
end

function [fmin, nmin] = enumerate_domain(PB, card)
% minimum of func_idx over the whole finite domain, blockwise (mixed radix
% counter), and the number of points attaining it (relative tolerance 1e-12)
n = numel(card); N = prod(card); fmin = Inf; nmin = 0;
z = ones(1, n);
for q = 1:N
    f = PB.func_idx(z);
    if isfinite(f)
        if ~isfinite(fmin) || f < fmin - 1e-12 * max(1, abs(fmin)), fmin = f; nmin = 1;
        elseif abs(f - fmin) <= 1e-12 * max(1, abs(fmin)), nmin = nmin + 1;
        end
    end
    i = 1;                                   % advance the mixed-radix counter
    while i <= n
        z(i) = z(i) + 1;
        if z(i) <= card(i), break; end
        z(i) = 1; i = i + 1;
    end
end
end

function PB = load_quiet(name)
try
    [~, PB] = evalc('feval(name)');
catch
    evalc('PB = feval(name);');
end
end

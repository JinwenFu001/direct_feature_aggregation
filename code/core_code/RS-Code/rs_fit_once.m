function rs_fit_once(input_dir, output_dir)
% Fit selected RS methods from CSV files written by R.
%
% Required input files in input_dir:
%   X_train.csv, y_train.csv, X_valid.csv, y_valid.csv, X_test.csv, y_test.csv
%   taxonomy.csv
%
% Optional input files:
%   gamma_range.csv, seed.csv, selection.txt, model_mode.txt, mu.csv
%   method_list.txt with requested methods, one per line or comma/space separated
%   For model_mode = article_aligned:
%     A_full.csv
%     C_RS_DL2.csv, gidx_RS_DL2.csv, cnorm_RS_DL2.csv
%     C_RS_CL2.csv, gidx_RS_CL2.csv, cnorm_RS_CL2.csv
%     C_RS_L1.csv, gidx_RS_L1.csv, cnorm_RS_L1.csv
%
% Output files in output_dir:
%   rs_summary.csv, rs_betas.csv

if nargin < 2 || isempty(output_dir)
    output_dir = input_dir;
end

set(0, 'DefaultFigureVisible', 'off');
rs_code_dir = fileparts(mfilename('fullpath'));
addpath(genpath(rs_code_dir));

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

X_train = readmatrix(fullfile(input_dir, 'X_train.csv'));
y_train = readmatrix(fullfile(input_dir, 'y_train.csv'));
X_valid = readmatrix(fullfile(input_dir, 'X_valid.csv'));
y_valid = readmatrix(fullfile(input_dir, 'y_valid.csv'));
X_test = readmatrix(fullfile(input_dir, 'X_test.csv'));
y_test = readmatrix(fullfile(input_dir, 'y_test.csv'));
taxonomy = round(readmatrix(fullfile(input_dir, 'taxonomy.csv')));

X_train = sparse(X_train);
X_valid = sparse(X_valid);
X_test = sparse(X_test);

y_train = y_train(:);
y_valid = y_valid(:);
y_test = y_test(:);

gamma_file = fullfile(input_dir, 'gamma_range.csv');
if exist(gamma_file, 'file')
    gamma_range = readmatrix(gamma_file);
    gamma_range = gamma_range(:)';
else
    gamma_range = logspace(log10(1e-4), log10(5e-2), 50);
end

seed_file = fullfile(input_dir, 'seed.csv');
if exist(seed_file, 'file')
    seed = readmatrix(seed_file);
    seed = seed(1);
else
    seed = 1;
end

selection_file = fullfile(input_dir, 'selection.txt');
if exist(selection_file, 'file')
    selection = strtrim(fileread(selection_file));
else
    selection = 'validation';
end

model_file = fullfile(input_dir, 'model_mode.txt');
if exist(model_file, 'file')
    model_mode = strtrim(fileread(model_file));
else
    model_mode = 'paper_rootless';
end

option.mu = 1e-2;
mu_file = fullfile(input_dir, 'mu.csv');
if exist(mu_file, 'file')
    mu_value = readmatrix(mu_file);
    option.mu = mu_value(1);
end
option.gammarange = gamma_range;
option.fig = 0;
option.verbose = false;
option.display_iter = false;
option.maxiter = 10000;
option.tol = 1e-7;
option.nfold = 5;

all_methods = {'RS-DL2', 1; 'RS-CL2', 2; 'RS-L1', 3};
all_beta_names = {'RS_DL2', 'RS_CL2', 'RS_L1'};
all_penalty_keys = {'RS_DL2', 'RS_CL2', 'RS_L1'};

method_file = fullfile(input_dir, 'method_list.txt');
if exist(method_file, 'file')
    method_text = strtrim(fileread(method_file));
    requested_methods = regexp(method_text, '[,\s]+', 'split');
    requested_methods = requested_methods(~cellfun('isempty', requested_methods));
    [is_known, method_idx] = ismember(requested_methods, all_methods(:, 1));
    if any(~is_known)
        error('Unknown requested RS method(s): %s', strjoin(requested_methods(~is_known), ', '))
    end
    method_idx = method_idx(is_known);
else
    method_idx = 1:size(all_methods, 1);
end

methods = all_methods(method_idx, :);
beta_names = all_beta_names(method_idx);
penalty_keys = all_penalty_keys(method_idx);
n_methods = size(methods, 1);
p = size(X_train, 2);

summary_method = cell(n_methods, 1);
summary_penalty = zeros(n_methods, 1);
summary_gamma = zeros(n_methods, 1);
summary_valid = zeros(n_methods, 1);
summary_test = zeros(n_methods, 1);
summary_intercept = zeros(n_methods, 1);
summary_time = zeros(n_methods, 1);
summary_status = cell(n_methods, 1);
betas = zeros(p, n_methods);

y_mean = mean(y_train);
y_train_c = y_train - y_mean;

article_aligned = strcmpi(model_mode, 'article_aligned');
if article_aligned
    A_full = sparse(readmatrix(fullfile(input_dir, 'A_full.csv')));
    node_order = round(readmatrix(fullfile(input_dir, 'node_order.csv')));
    root_node = round(readmatrix(fullfile(input_dir, 'root_node.csv')));
    root_col = find(node_order(:) == root_node(1), 1);
    if isempty(root_col)
        error('Could not find root node in node_order.')
    end
end

for imethod = 1:n_methods
    method_name = methods{imethod, 1};
    penalty = methods{imethod, 2};
    penalty_key = penalty_keys{imethod};
    current_intercept = NaN;
    method_timer = tic;
    try
        if article_aligned
            A_tree = A_full;
            C_tree = sparse(readmatrix(fullfile(input_dir, ['C_', penalty_key, '.csv'])));
            g_idx_tree = round(readmatrix(fullfile(input_dir, ['gidx_', penalty_key, '.csv'])));
            CNorm_tree = readmatrix(fullfile(input_dir, ['cnorm_', penalty_key, '.csv']));
            CNorm_tree = CNorm_tree(1);
            y_model = y_train;
            current_intercept = 0;
            use_exact_l1 = penalty == 3;
        else
            [A_tree, C_tree, CNorm_tree, g_idx_tree] = mat2SPGtree_rootless(taxonomy, penalty);
            A_tree = sparse(A_tree);
            C_tree = sparse(C_tree);
            y_model = y_train_c;
            current_intercept = y_mean;
            use_exact_l1 = false;
        end

        X_train_expanded = X_train * A_tree;

        if strcmpi(selection, 'cv')
            if use_exact_l1
                error('article_aligned exact RS-L1 currently supports validation selection only.')
            end
            rng(seed);
            [coef_rs, ~, gamma_opt, CV] = cv_SPG_cvrt( ...
                'group', y_model, X_train_expanded, [], ...
                C_tree, CNorm_tree, option, g_idx_tree);
            beta_rs = A_tree * coef_rs;
            valid_mse = NaN;
        else
            valid_mse_by_gamma = zeros(length(gamma_range), 1);
            beta_by_gamma = zeros(p, length(gamma_range));
            coef_init = zeros(size(A_tree, 2), 1);
            penalty_factor = ones(size(A_tree, 2), 1);
            if use_exact_l1
                penalty_factor(root_col) = 0;
            end
            for igamma = 1:length(gamma_range)
                gamma = gamma_range(igamma);
                if use_exact_l1
                    lasso_option = option;
                    lasso_option.b_init = coef_init;
                    coef_rs = weighted_lasso_fista( ...
                        y_model, X_train_expanded, gamma, penalty_factor, lasso_option);
                    coef_init = coef_rs;
                else
                    spg_option = option;
                    spg_option.b_init = coef_init;
                    coef_rs = SPG( ...
                        'group', y_model, X_train_expanded, gamma, 0, ...
                        C_tree, CNorm_tree, spg_option, g_idx_tree);
                    coef_init = coef_rs;
                end
                beta_by_gamma(:, igamma) = A_tree * coef_rs;
                pred_valid = current_intercept + X_valid * beta_by_gamma(:, igamma);
                valid_mse_by_gamma(igamma) = mean((y_valid - pred_valid).^2);
            end
            [valid_mse, best_idx] = min(valid_mse_by_gamma);
            gamma_opt = gamma_range(best_idx);
            beta_rs = beta_by_gamma(:, best_idx);
            CV = valid_mse_by_gamma(:)';
        end

        pred_test = current_intercept + X_test * beta_rs;
        test_mse = mean((y_test - pred_test).^2);

        summary_method{imethod} = method_name;
        summary_penalty(imethod) = penalty;
        summary_gamma(imethod) = gamma_opt;
        summary_valid(imethod) = valid_mse;
        summary_test(imethod) = test_mse;
        summary_intercept(imethod) = current_intercept;
        summary_time(imethod) = toc(method_timer);
        summary_status{imethod} = 'ok';
        betas(:, imethod) = beta_rs;

        cv_file = fullfile(output_dir, [method_name, '_cv.csv']);
        writematrix([gamma_range(:), CV(:)], cv_file);
    catch ME
        summary_method{imethod} = method_name;
        summary_penalty(imethod) = penalty;
        summary_gamma(imethod) = NaN;
        summary_valid(imethod) = NaN;
        summary_test(imethod) = NaN;
        summary_intercept(imethod) = current_intercept;
        summary_time(imethod) = toc(method_timer);
        summary_status{imethod} = ME.message;
        betas(:, imethod) = NaN;
    end
end

summary_table = table( ...
    summary_method, summary_penalty, summary_gamma, summary_valid, ...
    summary_test, summary_intercept, summary_time, summary_status, ...
    'VariableNames', {'method', 'penalty', 'gamma', 'valid_mse', ...
    'test_mse', 'intercept', 'elapsed_seconds', 'status'});

writetable(summary_table, fullfile(output_dir, 'rs_summary.csv'));
beta_table = array2table(betas, 'VariableNames', beta_names);
writetable(beta_table, fullfile(output_dir, 'rs_betas.csv'));

end

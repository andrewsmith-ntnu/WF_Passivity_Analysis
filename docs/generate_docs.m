% generate_docs.m
% Scans .m files in the repository, extracts the leading comment blocks
% (function/script descriptions) and writes Markdown files under docs/.

root = fileparts(mfilename('fullpath'));
if isempty(root)
    root = pwd; % when run from MATLAB path
end
repo_root = fileparts(root);
addpath(repo_root);

m_files = dir(fullfile(repo_root, '**', '*.m'));
functions_md = fullfile(repo_root, 'docs', 'functions.md');
scripts_md = fullfile(repo_root, 'docs', 'scripts.md');

% Collect entries first for index ordering
func_entries = struct('name', {}, 'path', {}, 'summary', {}, 'header', {});
script_entries = struct('name', {}, 'path', {}, 'summary', {}, 'header', {});

for k = 1:numel(m_files)
    full = fullfile(m_files(k).folder, m_files(k).name);
    % skip this generator file
    if strcmp(full, mfilename('fullpath'))
        continue
    end
    % read file
    txt = fileread(full);
    % split lines
    lines = regexp(txt, '\r?\n', 'split');
    % extract leading comment block (lines starting with % before first non-comment non-empty)
    header = {};
    in_header = true;
    for i = 1:min(numel(lines), 200)
        line = strtrim(lines{i});
        if startsWith(line, '%') && in_header
            header{end+1} = regexprep(line, '^%\s?', '');
        elseif isempty(line) && in_header
            header{end+1} = '';
        else
            break;
        end
    end

    % find one-line summary: first non-empty header line
    summary = '';
    for h = 1:numel(header)
        if ~isempty(strtrim(header{h}))
            summary = header{h};
            break;
        end
    end

    % detect if file contains a function definition within first 20 lines
    is_function = false;
    for i = 1:min(numel(lines), 20)
        if contains(lines{i}, 'function')
            is_function = true;
            break;
        end
    end

    entry.name = m_files(k).name;
    entry.path = full;
    entry.summary = summary;
    entry.header = header;

    if is_function
        func_entries(end+1) = entry; %#ok<SAGROW>
    else
        script_entries(end+1) = entry; %#ok<SAGROW>
    end
end

% write functions.md with index
fid_f = fopen(functions_md, 'w');
fid_s = fopen(scripts_md, 'w');

fprintf(fid_f, '# Functions\n\n');
if ~isempty(func_entries)
    fprintf(fid_f, '## Index\n\n');
    for i = 1:numel(func_entries)
        name = func_entries(i).name;
        summary = func_entries(i).summary;
        if isempty(summary)
            summary = '_No summary found._';
        end
        fprintf(fid_f, '- [%s](#%s) — %s\n', name, lower(strrep(name, '.', '')), summary);
    end
    fprintf(fid_f, '\n---\n\n');
end

for i = 1:numel(func_entries)
    e = func_entries(i);
    fprintf(fid_f, '## %s\n\n', e.name);
    if ~isempty(e.summary)
        fprintf(fid_f, '%s\n\n', e.summary);
    end
    if ~isempty(e.header)
        for h = 1:numel(e.header)
            fprintf(fid_f, '%s\n', e.header{h});
        end
    else
        fprintf(fid_f, '_No header found._\n');
    end
    fprintf(fid_f, '\n---\n\n');
end

% write scripts.md with index
fprintf(fid_s, '# Scripts\n\n');
if ~isempty(script_entries)
    fprintf(fid_s, '## Index\n\n');
    for i = 1:numel(script_entries)
        name = script_entries(i).name;
        summary = script_entries(i).summary;
        if isempty(summary)
            summary = '_No summary found._';
        end
        fprintf(fid_s, '- [%s](#%s) — %s\n', name, lower(strrep(name, '.', '')), summary);
    end
    fprintf(fid_s, '\n---\n\n');
end

for i = 1:numel(script_entries)
    e = script_entries(i);
    fprintf(fid_s, '## %s\n\n', e.name);
    if ~isempty(e.summary)
        fprintf(fid_s, '%s\n\n', e.summary);
    end
    if ~isempty(e.header)
        for h = 1:numel(e.header)
            fprintf(fid_s, '%s\n', e.header{h});
        end
    else
        fprintf(fid_s, '_No header found._\n');
    end
    fprintf(fid_s, '\n---\n\n');
end

fclose(fid_f);
fclose(fid_s);

fprintf('Documentation generated: %s, %s\n', functions_md, scripts_md);

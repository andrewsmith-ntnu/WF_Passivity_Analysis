function generate_docs_public(targetDir)
% generate_docs_public  Generate documentation for the Public folder
%
% Usage:
%   generate_docs_public()         % auto-detect the Public folder (parent of this docs folder)
%   generate_docs_public(targetDir) % provide the folder to document
%
% This function scans .m files in the target directory, extracts one-line
% summaries and short help blocks (including comment blocks immediately
% following a function signature), and writes three Markdown files into a
% `docs` subfolder:
%   - functions.md   : entries for files that contain a function
%   - scripts.md     : entries for script files
%   - Public_overview.md : a compact list of all files with one-line summaries

if nargin < 1 || isempty(targetDir)
    % Determine targetDir as the parent of this file's folder
    thisfull = mfilename('fullpath');
    if isempty(thisfull)
        % fallback: use current working directory
        targetDir = pwd;
    else
        docsFolder = fileparts(thisfull);
        targetDir = fileparts(docsFolder);
    end
end

% Normalize path
targetDir = char(targetDir);
if ~isfolder(targetDir)
    error('Target directory does not exist: %s', targetDir);
end

% List .m files in target dir (non-recursive)
fl = dir(fullfile(targetDir, '*.m'));
if isempty(fl)
    warning('No .m files found in %s', targetDir);
end

outDir = fullfile(targetDir, 'docs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

functionsList = {};
scriptsList = {};

for k = 1:numel(fl)
    fname = fl(k).name;
    fpath = fullfile(targetDir, fname);
    try
        txt = fileread(fpath);
    catch
        warning('Could not read %s', fpath);
        continue
    end
    lines = regexp(txt, '\r\n|\n|\r', 'split');
    % find the function signature and any comment block immediately after it
    func_idx = [];
    rawtext = txt;
    % Try to capture: (signature line) followed by an optional contiguous block of comment lines
    token = regexp(rawtext, '(^\s*function[^\r\n]*)(?:\r\n|\n|\r)((?:[ \t]*%[^\r\n]*(?:\r\n|\n|\r))*)', 'once', 'tokens');
    sig_line = '';
    header_block = '';
    if ~isempty(token)
        sig_line = strtrim(token{1});
        header_block = token{2};
        % compute the approximate line index of the signature
        prefix = rawtext(1:regexp(rawtext, ['(^\s*function[^\r\n]*)'], 'once')-1);
        func_idx = numel(regexp(prefix, '\r\n|\n|\r')) + 1;
    else
        % fallback: check if file contains any function definition
        if ~isempty(regexp(rawtext, '^\s*function\b', 'once')) || ~isempty(regexp(rawtext, '\r\n\s*function\b', 'once'))
            % find line index by scanning lines for the function keyword
            for i = 1:numel(lines)
                if ~isempty(regexp(strtrim(lines{i}), '^function\b', 'once'))
                    func_idx = i; sig_line = strtrim(lines{i}); break
                end
            end
        end
    end
    % Extract top-of-file leading comment block (for scripts or fallback)
    top_comments = {};
    for i = 1:min(numel(lines),50)
        ln = strtrim(lines{i});
        if startsWith(ln, '%')
            top_comments{end+1} = regexprep(ln, '^%\s?', ''); %#ok<AGROW>
        elseif isempty(ln)
            % allow initial blank lines before comments
            continue
        else
            break
        end
    end
    % Initialize summary/help
    summary = '';
    signature = '';
    help_text = '';
    if ~isempty(func_idx)
        % It's a function file
        % capture the function signature line (first function line)
        signature = strtrim(lines{func_idx});
        % 1) Try to find a comment block immediately after signature
        header_after = {};
        jstart = func_idx + 1;
        % skip initial blanks
        while jstart <= min(numel(lines), func_idx+200) && isempty(strtrim(lines{jstart}))
            jstart = jstart + 1;
        end
        if jstart <= min(numel(lines), func_idx+200) && startsWith(strtrim(lines{jstart}), '%')
            seen_comment = false;
            for j = jstart:min(numel(lines), func_idx+200)
                ln = strtrim(lines{j});
                if startsWith(ln, '%')
                    header_after{end+1} = regexprep(ln, '^%\s?', ''); %#ok<AGROW>
                    seen_comment = true;
                elseif isempty(ln) && seen_comment
                    header_after{end+1} = ''; %#ok<AGROW>
                else
                    break;
                end
            end
        end
        if ~isempty(header_after)
            % take the first non-empty line as one-line summary
            for m=1:numel(header_after)
                if ~isempty(strtrim(header_after{m}))
                    summary = strtrim(header_after{m});
                    break
                end
            end
        end
        % build full help_text from the header_after block if present
        if ~isempty(header_after)
            help_text = strjoin(header_after, '\n');
        end
        % 2) fallback: trailing inline comment on signature line
        if isempty(summary)
            sigln = lines{func_idx};
            pct = strfind(sigln, '%');
            if ~isempty(pct)
                % take text after first % on the signature line
                inline_text = strtrim(regexprep(sigln(pct(1):end), '^%\s?', ''));
                summary = inline_text;
                if isempty(help_text)
                    help_text = inline_text;
                end
            end
        end
        % 3) fallback: top-of-file leading comments
        if isempty(summary) && ~isempty(top_comments)
            % find first non-empty top comment
            for m=1:numel(top_comments)
                if ~isempty(strtrim(top_comments{m}))
                    summary = strtrim(top_comments{m});
                    break
                end
            end
        end
        % if no help_text yet, but top_comments exist, use them
        if isempty(help_text) && ~isempty(top_comments)
            help_text = strjoin(top_comments, '\n');
        end
        if isempty(summary)
            summary = '_No summary found._';
        end
        functionsList{end+1} = struct('file', fname, 'signature', signature, 'summary', summary, 'help', help_text); %#ok<AGROW>
    else
        % Script file: use top-of-file leading comments as summary
        if ~isempty(top_comments)
            % first non-empty top comment
            for m=1:numel(top_comments)
                if ~isempty(strtrim(top_comments{m}))
                    summary = strtrim(top_comments{m});
                    break
                end
            end
            if isempty(summary)
                summary = '_No summary found._';
            end
        else
            summary = '_No summary found._';
        end
        help_text = '';
        if ~isempty(top_comments)
            help_text = strjoin(top_comments, '\n');
        end
        scriptsList{end+1} = struct('file', fname, 'summary', summary, 'help', help_text); %#ok<AGROW>
    end
end

% Write functions.md
fid = fopen(fullfile(outDir, 'functions.md'), 'w', 'n', 'UTF-8');
if fid == -1, error('Could not open functions.md for writing'); end
fprintf(fid, '# Functions in %s\n\n', targetDir);
if isempty(functionsList)
    fprintf(fid, '_No function files found._\n');
else
    for i = 1:numel(functionsList)
        entry = functionsList{i};
        fprintf(fid, '### %s\n\n', entry.file);
        if ~isempty(entry.signature)
            fprintf(fid, '**Signature:** `%s`\n\n', strtrim(entry.signature));
        end
        fprintf(fid, '**Summary:** %s\n\n', entry.summary);
        if isfield(entry, 'help') && ~isempty(entry.help)
            fprintf(fid, '#### Help\n\n');
            fprintf(fid, '```matlab\n');
            % write help text preserving newlines
            help_lines = regexp(entry.help, '\n', 'split');
            for hh = 1:numel(help_lines)
                fprintf(fid, '%s\n', help_lines{hh});
            end
            fprintf(fid, '```\n\n');
        end
    end
end
fclose(fid);

% Write scripts.md
fid = fopen(fullfile(outDir, 'scripts.md'), 'w', 'n', 'UTF-8');
if fid == -1, error('Could not open scripts.md for writing'); end
fprintf(fid, '# Scripts in %s\n\n', targetDir);
if isempty(scriptsList)
    fprintf(fid, '_No script files found._\n');
else
    for i = 1:numel(scriptsList)
        entry = scriptsList{i};
        fprintf(fid, '### %s\n\n', entry.file);
        fprintf(fid, '**Summary:** %s\n\n', entry.summary);
        if isfield(entry, 'help') && ~isempty(entry.help)
            fprintf(fid, '#### Help\n\n');
            fprintf(fid, '```matlab\n');
            help_lines = regexp(entry.help, '\n', 'split');
            for hh = 1:numel(help_lines)
                fprintf(fid, '%s\n', help_lines{hh});
            end
            fprintf(fid, '```\n\n');
        end
    end
end
fclose(fid);

% Write Public_overview.md
fid = fopen(fullfile(outDir, 'Public_overview.md'), 'w', 'n', 'UTF-8');
if fid == -1, error('Could not open Public_overview.md for writing'); end
fprintf(fid, '# Public folder overview\n\n');
fprintf(fid, 'This page lists the functions and scripts in the Public folder with a one-line summary. For full per-file details see `functions.md` and `scripts.md`.\n\n');
if ~isempty(functionsList)
    fprintf(fid, '## Functions\n\n');
    for i = 1:numel(functionsList)
        entry = functionsList{i};
        fprintf(fid, '- **%s** — %s\n', entry.file, entry.summary);
    end
    fprintf(fid, '\n');
end
if ~isempty(scriptsList)
    fprintf(fid, '## Scripts\n\n');
    for i = 1:numel(scriptsList)
        entry = scriptsList{i};
        fprintf(fid, '- **%s** — %s\n', entry.file, entry.summary);
    end
    fprintf(fid, '\n');
end
fclose(fid);

fprintf('Generated docs in %s\n', outDir);
end

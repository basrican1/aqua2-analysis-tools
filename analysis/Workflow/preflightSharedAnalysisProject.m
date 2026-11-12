function preflight = preflightSharedAnalysisProject(manifestDir, varargin)
% Validate a complete SharedAnalysisMetadata v1 declaration directory.
% This is a structural guard only: it does not load data or alter computation.
%
% Example:
%   p = preflightSharedAnalysisProject('/path/to/manifests/post2');

p = inputParser;
p.addParameter('Required', true, @islogical);
p.parse(varargin{:});
opts = p.Results;

preflight = struct('requested', strlength(string(manifestDir)) > 0, ...
    'validated', false, 'manifest_dir', char(string(manifestDir)), ...
    'validator', 'validateSharedAnalysisProject', 'schema_version', '', ...
    'shared_metadata_root', '', 'shared_metadata_revision', '', ...
    'shared_metadata_dirty', [], ...
    'report', struct(), 'message', 'No project manifest directory supplied.');

if ~preflight.requested
    if opts.Required
        error('post2:MissingProjectManifest', ...
            'A SharedAnalysisMetadata declaration directory is required for this run.');
    end
    return;
end

manifestDir = char(string(manifestDir));
if ~isfolder(manifestDir)
    error('post2:MissingProjectManifest', 'Project manifest directory does not exist: %s', manifestDir);
end
if exist('validateSharedAnalysisProject', 'file') ~= 2
    error('post2:MissingSharedValidator', ...
        ['validateSharedAnalysisProject is unavailable. Run BWA_setupPaths() or add ', ...
         'SharedAnalysisMetadata/adapters/matlab to the MATLAB path.']);
end

report = validateSharedAnalysisProject(manifestDir, 'ThrowOnError', false);
preflight.report = report;
preflight.schema_version = char(string(report.schema_version));
preflight.shared_metadata_root = char(string(report.shared_root));
[preflight.shared_metadata_revision, preflight.shared_metadata_dirty] = ...
    sharedGitProvenance(preflight.shared_metadata_root);
preflight.validated = logical(report.passed);
if ~report.passed
    error('post2:InvalidProjectManifest', ...
        'Shared project manifest validation failed with %d error(s): %s', ...
        report.n_errors, char(strjoin(unique(string(report.issues.message(report.issues.severity == "error"))), ' | ')));
end
preflight.message = sprintf('Validated SharedAnalysisMetadata schema v%s.', preflight.schema_version);
end

function [revision, dirty] = sharedGitProvenance(sharedRoot)
% Best-effort provenance: an exported result can retain this even if the
% shared resource is later updated. Non-Git copies remain valid but unknown.
revision = '';
dirty = [];
if strlength(string(sharedRoot)) == 0 || ~isfolder(sharedRoot)
    return;
end
quotedRoot = strrep(char(sharedRoot), '"', '\"');
[status, output] = system(sprintf('git -C "%s" rev-parse HEAD', quotedRoot));
if status == 0
    revision = strtrim(output);
end
[status, output] = system(sprintf('git -C "%s" status --porcelain', quotedRoot));
if status == 0
    dirty = ~isempty(strtrim(output));
end
end

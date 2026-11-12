function root = setup_aqua2_tools()
% Add only this curated export to the MATLAB path. External tools are opt-in.
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'ingest')));
addpath(genpath(fullfile(root, 'analysis')));
addpath(genpath(fullfile(root, 'anatomy')));
addpath(genpath(fullfile(root, 'experimental')));
addpath(fullfile(root, 'examples'));
addpath(fullfile(root, 'hpc'));
end

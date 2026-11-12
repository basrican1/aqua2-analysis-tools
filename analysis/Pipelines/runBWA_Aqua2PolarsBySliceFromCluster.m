function loaded = runBWA_Aqua2PolarsBySliceFromCluster(pathIn, varargin)
    % Convenience wrapper:
    %   - loads BWA_AQ2 + canonical tables
    %   - runs BWA_Aqua2PolarsBySlice with filesTbl fallback wired
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(repoRoot));

    loaded = loadClusterOutputs(pathIn, 'attachLegacyMeta', false);
    BWA_Aqua2PolarsBySlice(loaded.BWA_AQ2, 'filesTbl', loaded.files, 'eventsTbl', loaded.events, varargin{:});
end

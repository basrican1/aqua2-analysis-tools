function BWA_AQ2 = attachLegacyMetaFromCanonical(BWA_AQ2, filesTbl)
    % Attach legacy-style opts/DB fields to BWA_AQ2 from canonical files table.
    % Safe no-op if opts/DB already exist.

    if isfield(BWA_AQ2, 'opts') && isfield(BWA_AQ2, 'DB')
        return
    end
    [opts, DB] = buildLegacyMetaFromCanonical(filesTbl, BWA_AQ2.fts1);
    BWA_AQ2.opts = opts;
    BWA_AQ2.DB = DB;
end

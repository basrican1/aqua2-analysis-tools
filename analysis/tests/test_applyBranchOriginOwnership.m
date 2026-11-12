function tests = test_applyBranchOriginOwnership
tests = functiontests(localfunctions);
end

function testBranchOriginOwnershipReassignsPrimaryCell(testCase)
imageSz = [8 8 1];
f1 = sub2ind([8 8], 4, 4);
f2 = [f1; sub2ind([8 8], 4, 5)];
f3 = [f2; sub2ind([8 8], 2, 2); sub2ind([8 8], 2, 3); sub2ind([8 8], 3, 2)];
f4 = [f3; sub2ind([8 8], 2, 4); sub2ind([8 8], 3, 3)];
xSpaTemp = {double([f1; 64 + f2; 128 + f3; 192 + f4])};

mask1 = false(8, 8); mask1(4:5, 4:5) = true;
mask2 = false(8, 8); mask2(1:4, 1:4) = true;

BWA_AQ2 = struct();
BWA_AQ2.fts1.loc.xSpaTemp = xSpaTemp;
BWA_AQ2.fts1.loc.t0 = 1;
BWA_AQ2.fts1.loc.t1 = 4;
BWA_AQ2.fts1.fileID = {'f1'};
BWA_AQ2.fts1.region.cell.mask = {mask1, mask2};
BWA_AQ2.fts1.region.cell.memberIdx = [1 2 NaN];
BWA_AQ2.opts.sz = imageSz;
BWA_AQ2.opts.spatialRes = 1;

canonicalOut = struct();
canonicalOut.subset.events = table(string("E1"), string("f1"), string("C1"), ...
    'VariableNames', {'event_uid','file_id','primary_cell_uid'});
canonicalOut.subset.cells = table(string(["C1";"C2"]), string(["f1";"f1"]), ...
    'VariableNames', {'cell_uid','file_id'});
canonicalOut.subset.event_cell_map = table(string(["E1";"E1"]), string(["C1";"C2"]), logical([true;false]), ...
    'VariableNames', {'event_uid','cell_uid','is_primary_cell'});

[out, info] = apply_branch_origin_ownership(canonicalOut, BWA_AQ2, ...
    'minBranchPixels', 3, 'minBranchPersistence', 2);

verifyTrue(testCase, info.applied);
verifyEqual(testCase, string(out.subset.events.onset_cell_uid), "C1");
verifyEqual(testCase, string(out.subset.events.branch_origin_cell_uid), "C2");
verifyEqual(testCase, string(out.subset.events.primary_cell_uid), "C2");
verifyTrue(testCase, logical(out.subset.event_cell_map.is_primary_cell(2)));
verifyFalse(testCase, logical(out.subset.event_cell_map.is_primary_cell(1)));
end

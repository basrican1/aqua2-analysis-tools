function tests = test_wavefrontMetricsFromGeometry
tests = functiontests(localfunctions);
end

function testSharedGeometryAndWavefrontMetrics(testCase)
imageSz = [5 5 1];
t0 = 1;
t1 = 4;

% Frame 1: (3,3)
% Frame 2: add (4,3)
% Frame 3: add (5,3)
% Frame 4: add (5,4)
f1 = sub2ind([5 5], 3, 3);
f2 = [f1; sub2ind([5 5], 3, 4)];
f3 = [f2; sub2ind([5 5], 3, 5)];
f4 = [f3; sub2ind([5 5], 4, 5)];
xSpaTemp_ID = [f1; 25 + f2; 50 + f3; 75 + f4];

xSpaTemp = {double(xSpaTemp_ID)};
allt0 = t0;
allt1 = t1;
spatialRes = 1;

out = CalcCentroidAndWavefrontMetrics(xSpaTemp, allt0, allt1, imageSz, spatialRes);

verifyEqual(testCase, numel(out.geometry), 1);
verifyEqual(testCase, out.geometry{1}.firstFrontierIdx, 1);
verifyEqual(testCase, out.centroid.distTot(1), 1.35355339059327, 'AbsTol', 1e-6);
verifyEqual(testCase, out.centroid.netDisplacementUm(1), 1.2747548783982, 'AbsTol', 1e-6);
verifyGreaterThan(testCase, out.centroid.tortuosity(1), 1);

W = out.wavefront{1};
verifyEqual(testCase, W.wavefront_max_reach_px, sqrt(5), 'AbsTol', 1e-6);
verifyEqual(testCase, W.wavefront_net_displacement_px, sqrt(5), 'AbsTol', 1e-6);
verifyGreaterThan(testCase, W.wavefront_frontier_length_px, W.wavefront_max_reach_px);
verifyGreaterThan(testCase, W.wavefront_tortuosity, 1);
verifyEqual(testCase, W.peakFrame, 4);
verifyTrue(testCase, isfinite(W.branchOriginFrame));
verifyEqual(testCase, W.branchPeakFrame, 4);
verifyLessThanOrEqual(testCase, W.branchOriginFrame, W.branchPeakFrame);
end

function testWavefrontIgnoresDisconnectedRoguePixels(testCase)
imageSz = [7 7 1];
t0 = 1;
t1 = 4;

% Main event expands rightward from the center.
f1 = sub2ind([7 7], 4, 4);
f2 = [f1; sub2ind([7 7], 4, 5)];
f3 = [f2; sub2ind([7 7], 4, 6); sub2ind([7 7], 1, 1)];
f4 = [f3; sub2ind([7 7], 4, 7)];
xSpaTemp_ID = [f1; 49 + f2; 98 + f3; 147 + f4];

xSpaTemp = {double(xSpaTemp_ID)};
out = CalcCentroidAndWavefrontMetrics(xSpaTemp, t0, t1, imageSz, 1);
W = out.wavefront{1};

verifyEqual(testCase, W.onsetOriginXY, [4 4], 'AbsTol', 1e-6);
verifyEqual(testCase, W.peakPointXY, [7 4], 'AbsTol', 1e-6);
verifyLessThan(testCase, W.wavefront_max_reach_px, 4);
verifyEqual(testCase, W.peakFrame, 4);
end

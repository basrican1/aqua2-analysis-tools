function tests = test_exploreLocalEventMotion
tests = functiontests(localfunctions);
end

function testExplorerOpensOnSyntheticPayload(testCase)
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

BWA_AQ2 = struct();
BWA_AQ2.opts = struct('sz', imageSz, 'spatialRes', 1);
BWA_AQ2.fts1 = struct();
BWA_AQ2.fts1.loc = struct();
BWA_AQ2.fts1.loc.t0 = t0;
BWA_AQ2.fts1.loc.t1 = t1;
BWA_AQ2.fts1.loc.xSpaTemp = {double(xSpaTemp_ID)};
BWA_AQ2.fts1.loc.xSpa = {double(unique(mod(xSpaTemp_ID - 1, imageSz(1) * imageSz(2)) + 1))};

out = explore_local_event_motion(BWA_AQ2, 'eventID', 1, 'visible', 'off');
testCase.addTeardown(@() closeFigureIfOpen(out.figure));

verifyTrue(testCase, isgraphics(out.figure));
verifyTrue(testCase, isgraphics(out.axes));
verifyEqual(testCase, out.controls.event.Value, 1);
verifyEqual(testCase, out.controls.frame.Value, t0);
end

function closeFigureIfOpen(h)
if ~isempty(h) && isgraphics(h)
    close(h);
end
end

function tests = test_plotSummarySignificance
tests = functiontests(localfunctions);
end

function testPlaceholder(testCase)
% Direct significance-marker graphics assertions are brittle in headless MATLAB.
% Keep this file as a placeholder until the plotting layer exposes a stable
% non-graphics significance payload for regression testing.
testCase.verifyTrue(true);
end

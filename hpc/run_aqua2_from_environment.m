function result = run_aqua2_from_environment()
% Invoke a separately installed, reviewed runner. No engine code is bundled.
runner = strtrim(getenv('AQUA2_RUNNER'));
assert(isvarname(runner), 'AQUA2_RUNNER must be a MATLAB function name.');
inputTif = getenv('AQUA2_INPUT'); presetCsv = getenv('AQUA2_PRESET');
assert(isfile(inputTif) && isfile(presetCsv), 'Explicit input and preset files required.');
assert(exist(runner,'file') == 2, 'Install the reviewed external runner first.');
result = feval(runner, inputTif, presetCsv);
end

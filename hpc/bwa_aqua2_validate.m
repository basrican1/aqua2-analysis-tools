function info = bwa_aqua2_validate(inputTif, cellMat, somaMat)
% Validate OME metadata and optional AQuA2 region MAT files.

arguments
    inputTif (1,:) char
    cellMat (1,:) char = ''
    somaMat (1,:) char = ''
end

bioformatsRoot = getenv('BIOFORMATS_ROOT');
helpersRoot = getenv('AQUA2_HELPERS_ROOT');
assert(isfolder(bioformatsRoot), 'Set BIOFORMATS_ROOT to your Bio-Formats MATLAB installation.');
assert(isfolder(helpersRoot), 'Set AQUA2_HELPERS_ROOT to your OME helper installation.');
addpath(genpath(bioformatsRoot), '-end');
addpath(genpath(helpersRoot), '-end');

assert(isfile(inputTif), 'Input TIFF does not exist: %s', inputTif);
[frameInterval, spatialRes] = BWAOmeHeaderRead(inputTif);

reader = bfGetReader(inputTif);
cleanup = onCleanup(@() reader.close());
meta = reader.getMetadataStore();
width = double(meta.getPixelsSizeX(0).getValue());
height = double(meta.getPixelsSizeY(0).getValue());
depth = double(meta.getPixelsSizeZ(0).getValue());
frames = double(meta.getPixelsSizeT(0).getValue());

info = struct( ...
    'inputTif', inputTif, ...
    'size', [height, width, depth, frames], ...
    'frameIntervalSeconds', frameInterval, ...
    'spatialResolutionMicrons', spatialRes, ...
    'cellCount', validateRegionFile(cellMat, 'cell'), ...
    'somaCount', validateRegionFile(somaMat, 'soma'));

fprintf('Input: %s\n', inputTif);
fprintf('OME size [Y X Z T]: %d %d %d %d\n', info.size);
fprintf('OME frame interval: %.9g seconds\n', frameInterval);
fprintf('OME spatial resolution: %.9g microns/pixel\n', spatialRes);
fprintf('Regions: %d cells, %d somas\n', info.cellCount, info.somaCount);
end

function count = validateRegionFile(path, label)
if isempty(path)
    count = 0;
    return
end
assert(isfile(path), '%s region MAT does not exist: %s', label, path);
data = load(path, 'bd0');
assert(isfield(data, 'bd0') && iscell(data.bd0), ...
    '%s region MAT must contain cell array bd0: %s', label, path);
count = numel(data.bd0);
end

function dataType = bwaDatTreeToDataType(datTree)
%BWADATTREETODATATYPE Map legacy datTree names to rule-set data types.
switch char(string(datTree))
    case 'basic.area'
        dataType = 'area';
    case 'curve.duration'
        dataType = 'duration';
    case 'curve.rise19'
        dataType = 'rise';
    case 'curve.fall91'
        dataType = 'fall';
    case 'curve.decayTau'
        dataType = 'decayTau';
    case 'curve.dffAUC'
        dataType = 'AUC';
    case 'curve.dffMax2'
        dataType = 'amplitude';
    case 'propagation.maxPropSpeed'
        dataType = 'maxspeed';
    otherwise
        dataType = '';
end
end

function scenarioRecord(scenarioFile, fileName)
% Generates the scenario input file

% Copyright 2018 The MathWorks, Inc.

helperSaveScenarioToFile(scenarioFile, fileName);
load(fileName);
RoadNetwork.RoadCenters = scenarioFile.RoadCenters;
RoadNetwork.NumLanes = scenarioFile.RoadSegments.NumLanes;
save(fileName, '-regexp', '^(?!(scenarioFile|fileName)$).')
end
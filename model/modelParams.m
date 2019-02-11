
% Copyright 2018 The MathWorks, Inc.

scenarioFileName = 'scenarioInput';
% Workaround for passing string to MATLAB function
scenarioFileNameU = uint8(scenarioFileName);
set_param([bdroot, '/Scene Management/Actor Pose Reader'], 'ScenarioFileName', scenarioFileName)
%% Initialization of workspace variables to be used by various model blocks
smartActors = 0; % Actor IDs of actors nominated as smart
numSmartActors = 0; % Number of nominated smart actors
numSmartActorsDM = 0; % Number of actors with driver model
numAddedRogueActors = 0; % Number of added rogue actors
smartActorsDM = 0; % Actor IDs of actors with driver model

%% Configurable parameters
maxVehicles = 25;
standStillGap = 10;
safeTimeGap = 1.4; % Vehicles must maintain atleast this time-gap w.r.t vehicle ahead for safety.
sensorRange = 600; % Vicinity range that smart vehicle monitors around.
numPlans = 2; % Number of Plans. Must match with number of implemented plans in plan-library
maxLanes = 10; % Maximum number of lanes
% Vehicle must be atleast this much close to the vehicle in front of it
% before even considering the feasibility of overtake maneuver.
minOvertakeProximity = 50;
% All externally added rogue actors have actor IDs > rogueActorIdBase
rogueActorIdBase = 4000;
%% Read from input simulation logs
scene = load(char(scenarioFileNameU));% Load the input mat file
% Set simulation time as timestamp of last log in the input .mat file
simTime = scene.vehiclePoses(end).SimulationTime;
numActors = length(scene.vehiclePoses(1).ActorPoses);
roadWidth = scene.RoadNetwork.RoadWidth;% Road width
laneWidth = roadWidth / (scene.RoadNetwork.NumMarkings - 1); % Lane width
numLanes = scene.RoadNetwork.NumLanes; % Num of lanes in the scenario.
laneLocationsY = scene.RoadNetwork.LaneMarkingLocation ;
% Read Init Poses info of actors from input simulation logs.
for i = 1 : numActors
    actorPoseInit(scene.vehiclePoses(1).ActorPoses(i).ActorID) = scene.vehiclePoses(1).ActorPoses(i);
end
clear i;
% Placing smart vehicles in middle of the right lane.
initYPos = (scene.RoadNetwork.RoadCenters(1, 2) - laneWidth / 2);
% Placing rogue vehicles in middle of the right lane.
initYPosRogue = (scene.RoadNetwork.RoadCenters(1, 2) - laneWidth / 2);
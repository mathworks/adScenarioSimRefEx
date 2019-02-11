classdef OvertakeQueries < Simulink.IntEnumType
    
    % Copyright 2018 The MathWorks, Inc.
    
    enumeration
        VehiclesAheadSameLane(0)
        VehiclesAheadLeftLane(1)
    end
    methods (Static)
        function retVal = getDefaultValue()
            retVal = OvertakeQueries.VehiclesAheadSameLane;
        end
    end
end
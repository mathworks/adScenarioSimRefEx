classdef visualizer < matlab.System ...
        & matlab.system.mixin.CustomIcon
    % Visualizer Creates and updates a bird's-eye plot.
    % You can opt for 'World-coordinate' view or 'Self-centered' view.
    % The Visualization block is for display purposes and you will have
    % to disable it to set the model in 'Rapid Accelerator' mode.
    %
    % See also: birdsEyePlot
    
    % Copyright 2017-2018 The MathWorks, Inc.
    
    %#codegen
    
    properties(Nontunable)
        %XLim Minimum and maximum distance in the longitudinal axis
        XLim = [-1 150]
        
        %YLim Minimum and maximum distance in the lateral axis
        YLim = [-30 30]
        
        %Number of smart actors
        NumSmartActors = 0
        
        %Visualization type
        View = 'Global view'
    end
    
    properties
        %Smart actor index
        smartActorIndex = 1
    end
    
    properties(Hidden, Constant)
        %To create the combo box list items
        ViewSet = ...
            matlab.system.StringSet({'Global view', 'Self-Centered view'})
    end
    
    properties(Nontunable, Logical)
        %HasActors Display actors data
        HasActors = true
    end
    
    properties(Nontunable, Logical)
        %HasRoads Display road boundary data
        HasRoads = false
        
        %HasLaneMarkings Display lane markings
        HasLaneMarkings = false
    end
    
    properties(Access=private)
        pFig
        pBEP
        pActorPlotter
        pActorProfile
        pLaneBoundaryPlotter
        pSimulinkUIToolbar
        pBlockName
        pIsLegendOn
        pLaneMarkingPlotter
    end
    
    methods
        function obj = visualizer(varargin)
            % Constructor
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Access=protected)
        
        function setupImpl(obj,varargin)
            wasFigureClosed = (isempty(obj.pFig) || ~ishghandle(obj.pFig));
            obj.pBlockName = gcb;
            if wasFigureClosed
                % Find the hidden figure handle
                root = groot;
                shh = get(root,'ShowHiddenHandles');
                set(root,'ShowHiddenHandles','on');
                hfig = findobj('Tag',obj.pBlockName);
                set(root,'ShowHiddenHandles',shh);
                if isempty(hfig)
                    hfig = figure('Name','Visualization','Tag',obj.pBlockName, ...
                        'Visible','off','NumberTitle','off');
                    obj.pIsLegendOn = true;
                else % Hide the figure while we build the bird's-eye plot
                    set(hfig,'Visible','off')
                    obj.pIsLegendOn = ~isempty(get(hfig.CurrentAxes,'Legend'));
                end
                
                obj.pFig = hfig;
            end
            
            % Create BEP before the toolbar because clf clears toolbars
            isBEPNeeded = (isempty(obj.pBEP) || wasFigureClosed);
            if isBEPNeeded
                clf(obj.pFig);
                hax = axes(obj.pFig);
                obj.pBEP = birdsEyePlot('Parent',hax,'XLim',obj.XLim,'YLim',obj.YLim);
                
                if obj.HasActors
                    obj.pActorPlotter = outlinePlotter(obj.pBEP);
                    obj.pActorProfile = struct( ...
                        'Length', 4.7, ...
                        'Width', 1.8, ...
                        'OriginOffset', [-1.35 0]);
                end
                
                if obj.HasRoads
                    grey = 0.3*ones(1,3);
                    obj.pLaneBoundaryPlotter = laneBoundaryPlotter(obj.pBEP,'DisplayName','Roads', 'Color', grey);
                end
                
                if obj.HasLaneMarkings
                    obj.pLaneMarkingPlotter = laneMarkingPlotter(obj.pBEP,'DisplayName','Lane Markings', 'FaceColor', [0.6 0.6 0.6]);
                end
                
            end
        end
        
        function resetImpl(obj)
            modelName = bdroot;
            
            % Create scope toolbar
            if isempty(obj.pFig.UserData) % Toolbar got cleared
                if isempty(obj.pSimulinkUIToolbar) % Toolbar was never created
                    t = findall(obj.pFig,'Type','uitoolbar');
                    if isempty(t)
                        t = uitoolbar(obj.pFig);
                    end
                    obj.pSimulinkUIToolbar = driving.internal.SimulinkUIToolbar(...
                        'Toolbar', t,...
                        'ModelName', modelName, ...
                        'BlockName', obj.pBlockName);
                else % Make sure that the toolbar is registered to model events
                    registerToModelButtonEvents(obj.pSimulinkUIToolbar);
                end
                userData.SimulinkUIToolbar = obj.pSimulinkUIToolbar;
                obj.pFig.UserData  = userData;
            else
                obj.pSimulinkUIToolbar = obj.pFig.UserData.SimulinkUIToolbar;
                registerToModelButtonEvents(obj.pSimulinkUIToolbar);
            end
            
            % Turn off the legend if it was off earlier
            if ~obj.pIsLegendOn
                legend(obj.pFig.CurrentAxes,'off')
            end
            
            % Bring the figure to front, set it to visible
            isDirty = get_param(bdroot,'Dirty');
            set_param(obj.pBlockName,'UserData',obj.pFig)
            set_param(obj.pBlockName,'OpenFcn','helperOpenFcn');
            set_param(bdroot,'Dirty',isDirty);
            set(obj.pFig,'Visible','on');
        end
        
        function releaseImpl(obj)
            % Release resources, such as file handles
            if ~isempty(obj.pFig.UserData)
                modelName = bdroot;
                isFastRestart = strcmp(get_param(modelName,'FastRestart'),'on');
                if isFastRestart %In fast restart mode, SimulationStatus is still 'running' at the end
                    setStoppedIcon(obj.pFig.UserData.SimulinkUIToolbar);
                else
                    update(obj.pFig.UserData.SimulinkUIToolbar);
                end
                release(obj.pFig.UserData.SimulinkUIToolbar);
            end
            dirtyFlag = get_param(bdroot,'Dirty');
            set_param(obj.pBlockName,'OpenFcn','');
            set_param(bdroot,'Dirty',dirtyFlag);
        end
        
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object
            % configuration, for the command line and System block dialog
            flag = false;
            
            hasActor = strcmpi(obj.View, 'Self-Centered view');
            if strcmp(prop, 'smartActorIndex') && hasActor
                flag = true;
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            
            if isLocked(obj)
                s.pFig = obj.pFig;
                s.pBEP = obj.pBEP;
                s.pActorPlotter         = obj.pActorPlotter;
                s.pActorProfile         = obj.pActorProfile;
                s.pLaneBoundaryPlotter  = obj.pLaneBoundaryPlotter;
                s.pIsLegendOn           = obj.pIsLegendOn;
                s.pSimulinkUIToolbar    = saveobj(obj.pSimulinkUIToolbar);
                s.pLaneMarkingPlotter   = obj.pLaneMarkingPlotter;
                s.pLaneMarkingPlotter   = obj.pLaneMarkingPlotter;
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked)
            if wasLocked
                obj.pFig = s.pFig;
                obj.pBEP = s.pBEP;
                obj.pActorPlotter           = s.pActorPlotter;
                obj.pActorProfile           = s.pActorProfile;
                obj.pLaneBoundaryPlotter    = s.pLaneBoundaryPlotter;
                obj.pIsLegendOn             = s.pIsLegendOn;
                obj.pSimulinkUIToolbar      = loadobj(s.pSimulinkUIToolbar);
                obj.pLaneMarkingPlotter     = s.pLaneMarkingPlotter;
                
                s = rmfield(s,{'pFig','pBEP','pActorPlotter','pActorProfile',...
                    'pLaneBoundaryPlotter','pSimulinkUIToolbar', ...
                    'pLaneMarkingPlotter'});
            end
            loadObjectImpl@matlab.System(obj,s,wasLocked);
        end
        
        function stepImpl(obj,varargin)
            %Update the Simulink control toolbar
            update(obj.pFig.UserData.SimulinkUIToolbar);
            
            % Update the bird's-eye plot if it is visible
            if strcmp(get(obj.pFig, 'Visible'), 'on')
                idx = 1;
                
                if obj.HasActors
                    actors = varargin{idx};
                    idx = idx+1;
                    plotActors(obj, actors);
                end
                
                if obj.HasRoads
                    % Road boundaries are in global view by default
                    rbSmart = varargin{idx};
                    if strcmp(obj.View, 'Self-Centered view')
                        % Self-Centered Road boundaries
                        rbSmart = driving.scenario.roadBoundariesToEgo({rbSmart}, actors(obj.smartActorIndex));
                    end
                    plotLanes(obj, rbSmart);
                    idx = idx+1;
                end
                
                if obj.HasLaneMarkings
                    % Lane markings are in global view by default
                    lm = varargin{idx};
                    if strcmp(obj.View, 'Self-Centered view')
                        % Self-Centered Lane markings
                        lm.LaneMarkingVertices = driving.scenario.laneMarkingsToEgo(lm.LaneMarkingVertices, actors(obj.smartActorIndex));
                    end
                    plotLaneMarkings(obj, lm);
                end
            end
        end
    end
    
    methods(Access=private)
        function plotActors(obj, actors)
            if (isnumeric(actors) && actors == 0)
                % Input is disconnected, so return
                return
            end
            
            numActors = numel(actors);
            actorPoses = actors(1 : numActors);
            actorProfile = obj.pActorProfile;
            pos = NaN(numActors, 2);
            yaw = NaN(numActors, 1);
            
            if strcmp(obj.View, 'Self-Centered view')
                % Self-Centered view of actors
                for m = 1:numActors
                    pos(m, :) = actorPoses(m).Position(1 : 2) - actorPoses(obj.smartActorIndex).Position(1 : 2);
                    yaw(m) = actorPoses(m).Yaw - actorPoses(obj.smartActorIndex).Yaw;
                end
            else
                % Global View of actors
                for m = 1:numActors
                    pos(m, :) = actorPoses(m).Position(1 : 2);
                    yaw(m) = actorPoses(m).Yaw;
                end
            end
            
            % Dimensions of all the actors in visualization
            actorProfileLength = actorProfile.Length * ones(numActors, 1);
            actorProfileWidth = actorProfile.Width * ones(numActors, 1);
            
            % increasing the dimensions of the smart actors to make it
            % distinguishable
            for i=1:obj.NumSmartActors
                actorProfileLength(i) = 4.7;
                actorProfileWidth(i) = 2.5;
            end
            
            plotOutline(obj.pActorPlotter, pos, yaw, ...
                actorProfileLength, ...
                actorProfileWidth, ...
                'OriginOffset', ones(numActors, 1) * actorProfile.OriginOffset);
        end
        
        function plotLanes(obj, rbSmart)
            if (isnumeric(rbSmart) && isscalar(rbSmart) && rbSmart == 0)
                % Input is disconnected, so return
                return
            end
            plotLaneBoundary(obj.pLaneBoundaryPlotter, {rbSmart});
        end
        
        function plotLaneMarkings(obj,lm)
            if (isnumeric(lm) && isscalar(lm) && lm == 0)
                % Input is disconnected, so return
                return
            end
            plotLaneMarking(obj.pLaneMarkingPlotter,lm.LaneMarkingVertices,lm.LaneMarkingFaces);
        end
        
    end
    
    % Simulink interface
    methods(Access=protected)
        function str = getIconImpl(~)
            str = sprintf('Visualizer');
        end
        
        function num = getNumInputsImpl(obj)
            num = 0;
            
            if obj.HasActors
                num = num+1;
            end
            if obj.HasRoads
                num = num+1;
            end
            if obj.HasLaneMarkings
                num = num+1;
            end
        end
        function varargout = getInputNamesImpl(obj)
            varargout = {};
            
            if obj.HasActors
                varargout = {varargout{:} 'Actors'};
            end
            if obj.HasRoads
                varargout = {varargout{:} 'Roads'};
            end
            if obj.HasLaneMarkings
                varargout = {varargout{:} 'Lane Markings'};
            end
        end
    end
    
    methods(Access = protected, Static)
        function header = getHeaderImpl
            % Define header panel for System block dialog
            header = matlab.system.display.Header(...
                'Title', 'Visualization',...
                'Text', getHeaderText());
        end
        
        function simMode = getSimulateUsingImpl
            % Return only allowed simulation mode in System block dialog
            simMode = 'Interpreted execution';
        end
        
        function flag = showSimulateUsingImpl
            % Return false if simulation mode hidden in System block dialog
            flag = false;
        end
    end
end

function str = getHeaderText
str = sprintf([...
    'The Visualization block creates and maintains a bird''s-eye plot.',...
    'You can opt for ''World-coordinate'' view or ''Self-centered'' view.']);
end
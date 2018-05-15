function outputs = GARD_ExecuteFlightPlan(Arguments)
% function outputs = GARD_ExecuteFlightPlan(Arguments)
%
%
%
% Flight Phases
% =============
% 1 - Taxi
% 2 - Takeoff / Initial Climb
% 3 - Climb
% 4 - En-Route
% 5 - Decent / Terminal
% 6 - Approach
% 7 - Missed Approach
%
% Arguments
% =========
% CurrentWaypoint  = Arguments(1);
% CurrentLatitude = Arguments(2);
% CurrentLongitude = Arguments(3);
% CurrentHeight = Arguments(4);
% CurrentFlightPhase = Arguments(5);

% Setup initial conditions
InitialConditions.Position = [-0.4726279 2.67082098 4000/3.28]';
InitialConditions.Velocity = [-51.6 28.8 0.0]';
InitialConditions.Attitude = EulerToQuat([0,0,152*pi/180])';
InitialConditions.AngularRates = [0 0 0]';
InitialConditions.FuelMass = 247;
InitialConditions.EngineSpeed = (2400/60) * 2*pi;
InitialConditions.SimulationDate = [15 01 2005];

% start in the climb
InitialConditions.FlightPhase = 3;

CurrentWaypoint  = Arguments(1);
CurrentLatitude = Arguments(2);
CurrentLongitude = Arguments(3);
CurrentHeight = Arguments(4);
CurrentFlightPhase = Arguments(5);














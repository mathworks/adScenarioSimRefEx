# Automated Driving Scenario Simulation Reference Example (Overtake Maneuver)
Overtake maneuver strategy allows a smart vehicle (an independent agent), to safely overtake the vehicles ahead of it by taking decisions based on the belief/knowledge it has about its environment.

This demo showcases a Simulink model architecture for creating and simulating synthetic scenarios. It reads as input the scenario file saved using the Driving Scenario Designer (DSD) application. This architecture allows creation of **synthetic scenarios**, by:
  * Marking an actor in the scenario as an autonomous smart actor.
  * Installing car-following (driver) model on some of the actors.
  * Introducing rogue actors (actors devoid of any intelligence) in the scenario.

**Driving scenario designer (DSD) application** is part of Automated Driving System Toolbox (ADST). Refer to the documentation [here](https://www.mathworks.com/help/driving/ref/drivingscenariodesigner-app.html) for more information.

Configuration parameters can be set for individual actors to observe the variations in the behavior. The plan algorithm in the smart actor supports overtake maneuver on straight road.

## MATLAB Toolbox dependencies
This model has been tested with MATLAB R2018b. 

To run this model, you need: MATLAB, Automated Driving System Toolbox (ADST), Model Predictive Control Toolbox, Simulink, Simulink Coder, and Stateflow.

## Running the model
The model needs a scenario file saved using the DSD application. A sample scenario file *scenarioInput.mat* is already included in the repository. 

Run the Simulink model **mOvertakeManeuver.slx** and see the visualization.

## Further documentation
See the file [Overtake-Maneuver.pdf](https://github.com/mathworks/Overtake-Maneuver/blob/master/Overtake-Maneuver.pdf) for detailed documentation covering - Scenario description, new scenario generation, Simulink model description, configuration and scenario visualization.

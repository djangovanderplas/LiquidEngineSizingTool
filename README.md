# LiquidEngineSizingTool
Welcome to the Liquid Engine Sizing Tool, a project aimed at creating a comprehensive engine design tool for liquid rocket engines. This tool provides a complete solution for sizing the combustion chamber, injector and nozzle geometry of your engine.

My goal with this project is not only to create a tool that can be useful for others, but also to learn more about the design of liquid rocket engines myself. By implementing this tool, I hope to gain deeper insights into the complex process of engine design.

## Roadmap
This is a small roadmap of stuff I want to implement:

**Phase 0: Design Selection**

This is basically the phase before any calculations where you look into the requirements of the engine. This is stuff like selecting propellant, defining thrust level etc.
- [x] Implement a way to find the optimal O/F ratio for Isp for a propellant combination
- [x] Implement a way to find an appropriate value for the L*
- [ ] Implement a way to find propellant combinations
- [ ] Improve the L* availablility
- [ ] Make a simple vehicle dynamics analyzer so you can find what Exit Pressure is most optimal for your engine
- [ ] Maybe add some kind of tradeoff tool so you can make choices

Known values: T, Pe, Pa_design, O/F, L*, Fuel, Oxidizer

**Phase 1: Preliminary Design**

This phase involves the initial layout of the engine in such a way that all engine dimensions are known. 
- [x] Implement initial sizing of engine characteristics

Known values: pressure ratio, mw, gamma, expansion ratio, combustion temperature, throat temperature, exit velocity, mass flow rate, Isp, throat area, exit area
- [x] Implement engine geometry
- [x] Implement various options for nozzles like bell and conical
- [x] Calculate engine volume
- [x] Calculate combustion chamber dimensions
- [ ] Plot the engine geometry (WIP)

**Phase 2: Initial analysis**

This phase involves analysing the requirements of the subsystems of the engine like the cooling and injectors.
- [ ] Simulate the heat flux through the engine
- [ ] Find out cooling channel requirements
- [ ] Simulate the injector requirements
- [ ] Find pressure and thermal stresses on material

**Phase 3: Initial Design**

This phase involves coming up with an initial design from the initial analysis
- [ ] Implement some kind of material selection method
- [ ] Implement the ability for multiple materials (multi material printing)
- [ ] Mass estimation for various materials
- [ ] Allow for material trade-off
- [ ] Allow for injector trade-off
- [ ] Analyse combustion stability in the injector

**Phase 4: Detail Analysis**

- [ ] Implement FEM methods to analyse stresses on materials
- [ ] Simulate fluid flow through the injector
- [ ] Simulate combustion instability
- [ ] Simulate combustion chamber efficiencies

**Phase 5: Detail Design**

- [ ] Allow for easy export to CAD software
- [ ] Allow for export of heat fluxes for FEM

**Phase 6: Validation**

- [ ] Allow for checking simulation with tested data
- [ ] Allow for determining pressures/temperatures at certain positions so it can be compared to sensor data

**Phase 7: Feed system simulation**
- [ ] Allow for more advanced feed system design simulations
- [ ] Allow for basic tank requirements

**Phase 8: Turbomachinery**
- [ ] Add turbomachinery design capabilities

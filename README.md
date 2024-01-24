# Amphiphiles
by R M L Evans. Initial coding 20/07/2010. Repository initiated 24/1/24.
A simple 2D lattice model of amphiphilic phase behaviour, useful for teaching about thermodynamic equilibrium and soft-matter systems.



An amphiphile is a molecule with a hydrophobic part and a hydrophilic part. When dissolved in water, amphiphiles spontaneously form non-trivial structures that balance entropy maximization against the energetic advantage of hiding the hyrdophbic parts from the water. 

This is a Metropolis Monte Carlo simulation of thermodynamic equilibrium (canonical NVT ensemble) in a simple 2D system that has the idealized features of an amphiphilic solution, such as and aqueous solution of surfactant (like washing-up liquid) or phospholipids (like in biological cell mambranes). 

The lattice is populated with a number of square particles with two hydrophobic (A) and two hydrophilic (B) sides. Unoccupied sites are assumed to contain solvent. Energy cost of each AA, AB and BB contact is defined. By default, AA contacts are advantageous (carrying one negative energy unit), mimicking the attractive hydrophobic interaction.

A configuration is shown in which two neighbouring sites are occupied with particles in an energetically expensive configuration (with one AB contact).

![image](https://github.com/RMLEvans/Amphiphiles/assets/30809235/5edc7890-ecc8-4725-b334-61a3eb7dc78b)


A Metropolis Monte Carlo trial move involves either swapping the occupancy of two randomly chosen sites or randomly rotating a particle at a randomly chosen occupied site.

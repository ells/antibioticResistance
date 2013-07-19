antibioticResistance
====================

A population model evaluating the consequences of the evolution of double-resistance and tradeoffs on the benefits of two-drug antibiotic treatments


	Model Classes:

	1.) SimulationManager  

		- user must input the desired treatment type delineated by string associated integers
                       
		- increments multi-resistant genotype's tradeoff and performs the selected simulation
                      
	2.) SimulationSettings 
		- handles declaration, setting, and getting of simulation constants


	3.) StaticSimulation   
		- responsible for all simulations where treatment regimens remain static
                       		----> no treatment control
                       		----> no mutation control
                       		----> single drug A
                       		----> single drug B
			----> separate cocktail
                       		----> combined cocktail
                       
	4.) CycleSimulation    
		- responsible for cycling simulation, cycling periodicity must be set from Manager or Settings

	4.) MixSimulation      
		- responsible for mixing simulation, switching periodicity must be set from Manager or Settings
     

                      
                      

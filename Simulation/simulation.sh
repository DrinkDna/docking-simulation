#pdb2gmx: protein preperation

grep -v HOH 1aki.pdb > 1AKI_clean.pdb

#protein topology
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce

#solvation (simulating a simple aqueous system)
#define the box and fill with solvent 
#editconf and genbox module

#define the box using editconf
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic

#fill the box with solvent using solvate 
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top

#adding ions 
#Since life does not exist at a net charge, we must add ions to our system = to make stable system 
#using genion
#genion read topology file and replace water mol with the ions as we specified 
#input file = .tpr = grompp module (GROMACS pre-processor). grompp use during 1st sim
#grompp module process the coordinate file and topology 
#topology file = describes the molcule = how atoms connected to each other (bonds, angles, dihedrals) 

#.tpr = topology parameter file = contains all atoms and their interactions with each other = combines topology file + coordinates + all ffp + simulation setting (in mdp file)
#.mdp = simulation setting = used to run energy minimization or run md simulation 

#grompp make .tpr from .mdp 

#assemble .tpr file with  .mdp .gro .top
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr

#now we have .tpr file. now replace water mol with ions 
#we have atomic level discription of our system (solute and solvent) in binary file .tpr 
#pass .tpr file to genion 

#-pname =positive ion 
#-nname =negative ion 
#-neutral the system 
#why we add topol.top = bcz every time it will be modify after each steps or process 
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

  #after this you will see the change in the topol.top file...no of SOL mol will be reduced than previous (10644 - 10636)s
  
#Energy minimization 

#we have solvated and electroneutral system assembled 
#make system is free from steric clashesh or inappropriate geometry 
#we make system relaxed through EM = make system free from steric clashesh or inappropriate geometry 

#Run the energy minimization through the GROMACS MD engine, mdrun
#you need minim.mdp file

#assemble minim.mdp + 1AKI_solv_ions.gro + topol.top into em.tpr
#and run Energy minimization using mdrun 
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr

gmx mdrun -v -deffnm em

#-v = prints the progress to the screen at every step of em

#use mdrun -s flag = if name is not rm.tpr 




#to dtermine EM succesfull or not / access whether EM is successfull or not = 2 key indicators Epot and MaxForce 
#2 factors = Epot = should be negative (10*5-10*6 magnitude)
#and MaxFore = The largest force acting on any single atom in the system after EM.

#our emtol (energy minimization tolerance) = 1000 (means the goal is to reduce the maximum force to below 1000 kJ/mol/nm.)

#high Fmax indicates 
#1 There are unresolved steric clashes
#2 system is not relaxed enough, making it unstable for MD simulation
#Fmax < emtol = good = go for md run 

#we have em.edr (binary enrgy file)  em.trr (trajectory file)  em.gro (energy minimized structure)

#extract enegy and then plot using Xmgrace

gmx energy -f em.edr -o potential.xvg

#enter 10 to select potential to get energy in potential.xvg file 

#Now that our system is at an energy minimum, we can begin real dynamics.


#___________________________________________________________________________________________________________________________________________________________________________________________

#EM ensured that we have a reasonable starting structure, in terms of geometry and solvent orientation

#Equilibrate the solvent and ions around the protein = statbilize solvent and ions without letting the protein change too much 

#make relax the Solvent allow to adjust the protein surface into it
#system will colapse without equilibration , as we heat the system solvent molecules reacts violently , explosive behaviour, huge repulsive forces

#apply position restraints to restrain heavy atoms of the protein except H = prevents distorting while the solvent equilibrates
#Restraints mean atoms can move only if the force is strong enough to overcome the restraint

#The protein is restrained so that only the solvent moves and "settles in" around it.

##system can be stable in two ways 
#1 NVT = controlled temperature = NVT allows us to gradually thermalize the system under fixed volume, preventing large fluctuations in energy or atomic positions
#prevent temperature fluctuation and system instabilities due to sudden high temp....thermostat gently increase temp...

#NVT stabilize tempt
#NPT stabilize pressure and density 

#make nvt.mdp file -c = coordinates -r = position restraint coordinates 
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -v 



#NVT equilibration, 

	#The system is heated to the desired temperature.

	#Atom velocities are generated.

	#The temperature is stabilized using V-rescale thermostat.

	#Volume and pressure remain fixed (pressure coupling off).

	#Position restraints are active to allow only solvent to equilibrate around the protein.

#analyze temp progression 
gmx energy -f nvt.edr -o temperature.xvg

#the temperature of the system quickly reaches the target value (300 K), and remains stable over the remainder of the equilibration



#NPT stabilize the pressure and density = isothermal-isobaric ensemble 

gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx_mpi mdrun -deffnm npt -v 

#then analyze the pressure and density progression 

gmx_mpi energy -f npt.edr -o pressure.xvg
gmx_mpi energy -f npt.edr -o density.xvg


#Production MD
#release the position restraints and run production MD for data collection
#make md.mdp file 

gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr

#-.cpt = continue from previous step 

gmx_mpi mdrun -deffnm md_0_1 -v 



#analysis 
#we have simulated our protein

#tools :
	#trjconv = post-processing tool to strip out coordinates
	#trajectory = records position of all atoms in the system at every time step during simulation 
	
	
gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center

	#-s md_0_1.tpr	Gives the reference structure and topology.
	#-f md_0_1.xtc	The raw trajectory (with broken protein).
	
#after above 
	#The protein is now whole in every frame.
	#centered in the box.
	 


# look at structural stability first
gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns

#calculate rmsd 

gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns


#rg plot 
#The radius of gyration of a protein is a measure of its compactness. If a protein is stably folded, 
#it will likely maintain a relatively steady value of Rg. If a protein unfolds, its Rg will change over time.



gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg








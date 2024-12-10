#!/bin/bash
Help()
{
  ans=$(zenity --text-info --title 'gmx_qk 1.0.1' \
      --html \
      --ok-label="Quit" \
      --extra-button="Tutorial" \
      --extra-button="Workflow" \
      --extra-button Simulation \
      --extra-button "Trajectory Analysis" \
      --extra-button "FEL (Free Energy Landscape)" \
      --width=600 \
      --height=600 \
      --filename=<(echo "
<html>
    <body>
        <div style='text-align:center; font-family: Arial, sans-serif;'>
            <h1 style='color: blue;'>Welcome to gmx_qk</h1>
            <img src="/usr/share/icons/gmx_qk" alt='gmx_qk Logo' style='display: block; margin: 20px auto; width: 200px; height: auto;'>
            <p style='font-size: 18px; text-align: justify;'>
                The gmx_qk is a <b>Zenity, Python, Gromacs, and g_mmpbsa</b> dependent bash program. 
                It's designed for beginners to GROMACS who wish to simulate <b>protein or protein-ligand complexes, 
                including MM/PBSA calculations</b>. 
                <br><br>
                gmx_qk is a fully automated program, efficiently compatible with GROMACS 5.0 and newer versions such as 2021.4. 
                Informative widgets are supported by <b>Zenity (GUI)</b>. 
                Additional functionality for <b>post-MD simulation trajectory analysis and Free Energy Landscape (FEL)</b> analysis has also been embedded.
                
            </p>
            <strong>Please cite this: J. Chem. Inf. Model. 2023, 63, 9, 2603–2608</strong>
            <p>
                <strong>Author:</strong><br>
                Dr. Harvinder Singh<br>
                Research Associate<br>
                Department of Translational and Regenerative Medicine, PGIMER, Chandigarh (160012)<br>
                <a href='mailto:harvindermaan4@gmail.com'>harvindermaan4@gmail.com</a>
            </p>
        </div>
    </body>
</html>
"))
  rc=$?
#  echo "${rc}-${ans}"
  echo $ans
  if [[ $ans = "Tutorial" ]]
  then
    theurl="https://github.com/harry-maan/gmx_qk/blob/main/gmx_qk.pdf"
    zenity --text-info --title="gmx_qk 1.0.0" --html --url=$theurl \
       --checkbox="I read it...and I'm good to go" --width 800 --height 600;
       Help;
  elif [[ $ans = "Workflow" ]]
  then
    theurl="https://github.com/harry-maan/gmx_qk/blob/main/workflow.png"
    zenity --text-info --title="gmx_qk 1.0.0" --html --url=$theurl \
       --checkbox="I read it...and I'm good to go" --width 800 --height 600;
       Help;
  elif [[ $ans = "License" ]]
  then
    theurl="https://github.com/harry-maan/gmx_qk/blob/main/license.png"
    zenity --text-info --title="gmx_qk 1.0.0" --html --url=$theurl \
       --checkbox="I read it...and I'm good to go" --width 800 --height 600;
       Help;
  elif [[ $ans = "Simulation" ]]
  then
    version_gmx;
  elif [[ $ans = "Trajectory Analysis" ]]
  then
    gmx_trj_qk;  
  elif [[ $ans = "FEL (Free Energy Landscape)" ]]
  then
    gmx_fel;
  fi
   
}

function simulation_type(){
  go_no_go=$(zenity --forms --title "Topology/Simulation/Simulation-MMPBSA" --text "Combo name" \
  --add-combo "Select a simulation type" --combo-values "Simulate the Protein-Ligand (P-L complex)|Simulate Protein-Protein (P-P complex)|P-L complex MD & Thermal MM/PBSA calculation|P-P complex MD & Thermal MM/PBSA calculation|Protein in water|Done")
  
  case "${go_no_go}" in
  "Simulate the Protein-Ligand (P-L complex)")
  choose1;;
  "Simulate Protein-Protein (P-P complex)")
  ppcomplex;;
  "P-L complex MD & Thermal MM/PBSA calculation")
  choose1;;
  "P-P complex MD & Thermal MM/PBSA calculation")
  ppcomplex;;
  "Protein in water")
  ppcomplex;;
  "Done") 
  break;;
esac
}
function version_gmx() {
  version=$(zenity --forms --title "Gromacs version" --text "For Smooth working Select Version /opt" --add-combo "Select a version do you like to select" --combo-values "/usr/local|/opt")
  case $? in
           0)
                  dir=$(zenity  --file-selection --title="Choose a working directory" --directory);cd $dir;simulation_type;;
           1)
                  zenity --question \
                         --title="Please select a gromacs version!" \
                         --text="No gromacs version selected. Do you want to select one?!" \
                         && version_gmx || exit;;
          -1)
                  echo "An unexpected error has occurred."; exit;;
  esac
}

function ppcomplex() {
  protein="$(zenity --file-selection --title='Select a Protein-Protein complex file(.pdb)')"
  case $? in
           0)
                  ppcomplex_simu;;
           1)
                  zenity --question \
                         --title="Please select a Protein-Protein complex file(.pdb)" \
                         --text="No Protein-Protein complex file(.pdb) selected. Do you want to select one?" \
                         && choose1 || exit;;
          -1)
                  echo "An unexpected error has occurred."; exit;;
  esac
}
function ppcomplex_simu() {
source $version/gromacs/bin/GMXRC || source $version/gromacs/bin/GMXRC
	gmx pdb2gmx -f "$protein" -ignh  -ff charmm27 -water tip3p |is_ok=$(zenity --text-info --title "Is there everything OK?" --width 600 --height 300)
	ok=$?
	if ((ok==0)); then
	source $version/gromacs/bin/GMXRC || source $version/gromacs/bin/GMXRC
	gmx=gmx
	        case "${go_no_go}" in
          "Simulate Protein-Protein (P-P complex)")
            parameters;;
            "P-P complex MD & Thermal MM/PBSA calculation")
            parameters;;
            "Protein in water")
            parameters;;          
  
  esac
  else
    zenity --text-info --title "Please check the Protein/Protein-Protein complex file (.pdb) for its correction" --width 600 --height 300
  fi
}

function choose1() {
  protein="$(zenity --file-selection --title='Select a receptor File(.pdb)')"
  case $? in
           0)
                  choose2;;
           1)
                  zenity --question \
                         --title="Please select a receptor File(.pdb)" \
                         --text="No receptor File(pdb) selected. Do you want to select one?" \
                         && choose1 || exit;;
          -1)
                  echo "An unexpected error has occurred."; exit;;
  esac
}
function choose2() {
  ligand="$(zenity --file-selection --title='Please select a ligand File(.pdb)')"
  case $? in
           0)
                  choose3;;
           1)
                  zenity --question \
                         --title="Select a ligand File(pdb)" \
                         --text="No ligand File(pdb). Do you want to select one?" \
                         && choose2 || exit;;
          -1)
                  echo "An unexpected error has occurred."; exit;;
  esac
}
function choose3() {
  ligitpf="$(zenity --file-selection --title='Select a ligand ITP File(itp)')"
  case $? in
           0)
                  topology;;
           1)
                  zenity --question \
                         --title="Select a ligand ITP File(itp)" \
                         --text="No ligand ITP File(itp). Do you want to select one?" \
                         && choose3 || exit;;
          -1)
                  echo "An unexpected error has occurred."; exit;;
  esac
}
function topo_simu(){
  cp $ligitpf $protein $ligand .
  base_name=$(basename ${ligitpf})
  go_no_go=$(zenity --forms --title "Topology/Simulation/Simulation-MMPBSA" --text "Combo name" \
  --add-combo "Select what you want to do" --combo-values "Generate Topology|Simulate the P-L complex|MD & Thermal MM/PBSA analysis|Done")
  
  case "${go_no_go}" in
  "Generate Topology")
  topology;;
  "Simulate the P-L complex")
  topology;;
  "MD & Thermal MM/PBSA analysis")
  topology;;
  "Done") 
  break;;
esac
}


#
function topology() {
  cp $ligitpf $protein $ligand .
  base_name=$(basename ${ligitpf})
	source $version/gromacs/bin/GMXRC || source $version/gromacs/bin/GMXRC
	gmx pdb2gmx -f "$protein" -ignh  -ff charmm27 -water tip3p |is_ok=$(zenity --text-info --title "Is there everything OK?" --width 600 --height 300)
	ok=$?
	if ((ok==0)); then
		gmx editconf -f "$ligand" -o lig.gro
cat << EOF > gmx_input_edit.py
#!/usr/bin/python3
import fnmatch
import os
def replace_line(file_name, line_num, text): # to replace the count of atoms
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()
with open('conf.gro', 'r') as file:
    with open('lig.gro', 'r') as foo:
    # read a list of lines into data
        dat2 = foo.readlines()
    data = file.readlines()
    
with open('conf.gro', 'w') as file:
    dat3 = data[0:-1] + dat2[2:-1]# concatenate the conf.gro and lig.gro
    file.writelines(dat3)
with open('conf.gro', 'r') as fp:
    for count, line in enumerate(fp): #count the lines in conf.gro file
        pass
    acount = str(count - 1) # conversion of integer into text
    acount1 = ' ' + acount + '\n' # text supplied with space and insert a new line
def replace_line(file_name, line_num, text): # to replace the count of atoms
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()
replace_line('conf.gro', 1, acount1)
def search_string_in_file(file_name, string_to_search):
    """Search for the given string in file and return lines containing that string,
    along with line numbers"""
    line_number = 0
    list_of_results = []
    # Open the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            line_number += 1
            if string_to_search in line:
                # If yes, then add the line number & line as a tuple in the list
                list_of_results.append((line_number, line.rstrip()))
    # Return list of tuples containing line numbers and lines where string is found
    return list_of_results
for file1 in os.listdir('.'):
    if fnmatch.fnmatch(file1, '$base_name'): #file1 for lig.itp 
        ligitp = '#include "' + file1 + '"' + "\n"
        matched_lines = search_string_in_file(file1, '[ moleculetype ]')
        print('Total Matched lines : ', len(matched_lines))
        for elem2 in matched_lines:
            elm3 = elem2[0] + 1
            replace_line(file1, elm3, 'ligand     3\n')
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*ol.top'):
        matched_lines = search_string_in_file(file, '.ff/forcefield.itp')
        print('Total Matched lines : ', len(matched_lines))
        for elem in matched_lines:
            elm1 = elem[0] #+ 1
        with open(file, 'r') as tip:
            data = tip.readlines()
            #print(file)
            #ligitp = '#include "' + file + '"' + data[19] + "\n"
            #ligitp =  '#include "lig_sam.itp"' + data[19] + '\n'
            LIG =  data[-1] + 'ligand              1\n'#replace_line(file, -1, LIG)
            replace_line(file, elm1, ligitp)
            replace_line(file, -1, LIG )
EOF
chmod 755 gmx_input_edit.py
./gmx_input_edit.py
rm -rf gmx_input_edit.py
	source $version/gromacs/bin/GMXRC || source $version/gromacs/bin/GMXRC
	gmx=gmx
	        case "${go_no_go}" in
          "Generate Topology")
            output;;
            "Simulate the Protein-Ligand (P-L complex)")
            parameters;;
            "P-L complex MD & Thermal MM/PBSA calculation")
            parameters;;
          
  
  esac
  else
    zenity --text-info --title "Please check the Protein and Ligand file for its correction" --width 600 --height 300
  fi
}
function output() {
  mkdir output
  cp *.itp *.top *.pdb *.gro output/
  tar -czf topology.tar.gz output/
  rm -rf *.itp *.top *.pdb *.gro *.mdp *.tpr *.cpt *.log *.ndx *.edr *.trr *.xtc *.py output/
}
function parameters(){
  data=$(zenity --forms --separator="," \
--title="MD Parameters" \
--text="1. Please enter digits only in picoseconds. \n 2. Please eneter digits only in picoseconds \n 3. String for output file name \n" \
--add-entry="Simulation Course (picoseconds)" \
--add-entry="Energy minimization" \
--add-entry="Name of output" --width 600 --height 300) ;

picoseconds=$( echo $data | awk -F ',' '{print $1}' )
Energy_minimization=$( echo $data | awk -F ',' '{print $2}' )
fname=$( echo $data | awk -F ',' '{print $3}' )
nsteps=$((picoseconds * 500))
frames=$((nsteps / 50000))
	        case "${go_no_go}" in
          "Simulate Protein-Protein (P-P complex)")
            simulate_protein;;
            "P-P complex MD & Thermal MM/PBSA calculation")
            parameters_mmpba;;
            "Simulate the Protein-Ligand (P-L complex)")
            simulate;;
            "P-L complex MD & Thermal MM/PBSA calculation")
            parameters_mmpba;;
            "Protein in water")
            simulate_protein;output_md;;
          
  
  esac
  
}
function parameters_mmpba(){
  if zenity --question  --title="Expected MD trajectory frames" --text="Total $frames will be produced. Do you want to reduce no. of frames to save time and computational power?" --width=300 --height=100
then 
    zenity --info --title="The total number of frames will be produced" --text="$frames frames, please sample the frames for mmpbsa calculations.\nDo you like to reduce no. of frames" --width=300 --height=100
  data=$(zenity --forms --separator="," \
--title="MM/PBSA Parameters" \
--text="1. Please enter digits only for sampling the fames.\n 2. Please enter digits ranges between 2-4.\n 3. Please enter digits only (mb)" \
--add-entry="Factor for reducing frames" \
--add-entry="Dielectric constant" \
--add-entry="Memory available for run");

factor=$( echo $data | awk -F ',' '{print $1}' )
dieC=$( echo $data | awk -F ',' '{print $2}' )
memory=$( echo $data | awk -F ',' '{print $3}' )
echo $factor 
echo $dieC
echo $memory
  case "${go_no_go}" in
  "Simulate the Protein-Ligand (P-L complex)")
  simulate;reduce_frames;output_mmpbsa;;
  "Simulate Protein-Protein (P-P complex)")
  simulate_protein;reduce_frames_protein;output_md;;
  "P-L complex MD & Thermal MM/PBSA calculation")
  simulate;reduce_frames;output_mmpbsa;;
  "P-P complex MD & Thermal MM/PBSA calculation")
  simulate_protein;reduce_frames_protein;output_mmpbsa;;
  "Protein in water")
  simulate_protein
  esac
  else
      data=$(zenity --forms --separator="," \
      --title="MM/PBSA Parameters" \
      --text="1. Please enter digits ranges between 2-4.\n2. Please enter digits only (MB)" \
      --add-entry="Dielectric constant" \
      --add-entry="Memory available for run" --width 300 --height 100);

dieC=$( echo $data | awk -F ',' '{print $1}' )
memory=$( echo $data | awk -F ',' '{print $2}' )
echo $dieC
echo $memory
  case "${go_no_go}" in
  "Simulate the Protein-Ligand (P-L complex)")
  simulate;no_reduction;output_md;;
  "Simulate Protein-Protein (P-P complex)")
  simulate_protein;no_reduction_protein;output_md;;
  "P-L complex MD & Thermal MM/PBSA calculation")
  simulate;no_reduction;output_mmpbsa;;
  "P-P complex MD & Thermal MM/PBSA calculation")
  simulate_protein;no_reduction_protein;output_mmpbsa;;
  "Protein in water")
  simulate_protein;output_md;;
  esac
fi
}

function simulate(){
em=$(cat << EOF > em.mdp
;
;	GROMACS
;	Energy Minimization Script
;
;
define		= -DFLEXIBLE	; pass to preprocessor
cpp		= usr/bin/cpp	; location of preprocessor
constraints	= none
integrator	= steep		; steepest decents minimum (else cg)
nsteps		= $Energy_minimization
;
;	Energy Minimizing Stuff
;
emtol		= 10  		; convergence total force(kJ/mol/nm) is smaller than
emstep		= 0.01		; initial step size (nm)
nstcomm		= 100		; frequency or COM motion removal
ns_type		= grid
rlist		= 1.4		; cut-off distance for short range neighbors
rcoulomb	= 1.1		; distance for coulomb cut-off
coulombtype	= PME		; electrostatics (Particle Mesh Ewald method)
fourierspacing	= 0.12		; max grid spacing when using PPPM or PME
vdw-type	= Shift
rvdw		= 1.1		; Verlet cutoff
Tcoupl		= no		; temperature coupling
Pcoupl		= no		; pressure coupling
gen_vel		= no
EOF
)
nvt=$(cat << EOF > nvt.mdp
title		= OPLS Lysozyme NVT equilibration 
define		= -DPOSRES	; position restrain the protein
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= 50000		; 2 * 50000 = 100 ps
dt		    = 0.002		; 2 fs
; Output control
nstxout		= 500		; save coordinates every 1.0 ps
nstvout		= 500		; save velocities every 1.0 ps
nstenergy	= 500		; save energies every 1.0 ps
nstlog		= 500		; update log file every 1.0 ps
; Bond parameters
continuation	        = no		; first dynamics run
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10		; 20 fs, largely irrelevant with Verlet
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME	; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		; cubic interpolation
fourierspacing	= 0.16	; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1           ; time constant, in ps
ref_t		= 300 	  300           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl		= no 		; no pressure coupling in NVT
; Periodic boundary conditions
pbc		= xyz		    ; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= yes		; assign velocities from Maxwell distribution
gen_temp	= 300		; temperature for Maxwell distribution
gen_seed	= -1		; generate a random seed
EOF
)
npt=$(cat << EOF > npt.mdp
title		= OPLS Lysozyme NPT equilibration 
define		= -DPOSRES	; position restrain the protein
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= 50000		; 2 * 50000 = 100 ps
dt		    = 0.002		; 2 fs
; Output control
nstxout		= 500		; save coordinates every 1.0 ps
nstvout		= 500		; save velocities every 1.0 ps
nstenergy	= 500		; save energies every 1.0 ps
nstlog		= 500		; update log file every 1.0 ps
; Bond parameters
continuation	        = yes		; Restarting after NVT 
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		    ; cubic interpolation
fourierspacing	= 0.16		; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1	        ; time constant, in ps
ref_t		= 300 	  300	        ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT
pcoupltype	        = isotropic	            ; uniform scaling of box vectors
tau_p		        = 2.0		            ; time constant, in ps
ref_p		        = 1.0		            ; reference pressure, in bar
compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com
; Periodic boundary conditions
pbc		= xyz		; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= no		; Velocity generation is off 
EOF
)

mdp=$(cat << EOF > md.mdp
title		= OPLS Lysozyme MD simulation 
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= $nsteps     ; 2 * 50000000 = 100000 ps (200 ns)
dt		= 0.002	; 2 fs
; Output control
nstxout		        = 0		; save coordinates every 10.0 ps
nstvout		        = 0		; save velocities every 10.0 ps
nstenergy	        = 500000		; save energies every 10.0 ps
nstlog		        = 500000		; update log file every 10.0 ps
nstxout-compressed  = 500000      ; save compressed coordinates every 10.0 ps
                                ; nstxout-compressed replaces nstxtcout
compressed-x-grps   = System    ; replaces xtc-grps
; Bond parameters
continuation	        = yes		; Restarting after NPT 
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		    ; cubic interpolation
fourierspacing	= 0.16		; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1	        ; time constant, in ps
ref_t		= 300 	  300	        ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT
pcoupltype	        = isotropic	            ; uniform scaling of box vectors
tau_p		        = 2.0		            ; time constant, in ps
ref_p		        = 1.0		            ; reference pressure, in bar
compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc		= xyz		; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= no		; Velocity generation is off
EOF
)





  gmx editconf -f conf.gro -d 1.0 -bt triclinic -o box.gro
  gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o box_sol.gro 
  gmx grompp -f em.mdp -c box_sol.gro -p topol.top -o ION.tpr -maxwarn 4
  gmx genion -s ION.tpr -p topol.top -conc 0.1 -neutral -o box_sol_ion.gro << EOF
15
EOF
  gmx grompp -f em.mdp -c box_sol_ion.gro -p topol.top -o EM.tpr -maxwarn 4
  gmx mdrun -v -deffnm EM
  gmx make_ndx -f lig.gro -o index_lig.ndx << EOF
0 & ! a H*
q
EOF
  gmx genrestr -f lig.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000 << EOF
3
EOF
cat << EOF > posre_inclusion.py
#!/usr/bin/python3
import fnmatch
import os
def search_string_in_file(file_name, string_to_search):
    """Search for the given string in file and return lines containing that string,
    along with line numbers"""
    line_number = 0
    list_of_results = []
    # Open the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            line_number += 1
            if string_to_search in line:
                # If yes, then add the line number & line as a tuple in the list
                list_of_results.append((line_number, line.rstrip()))
    # Return list of tuples containing line numbers and lines where string is found
    return list_of_results
matched_lines = search_string_in_file('topol.top', '; Include water topology')
print('Total Matched lines : ', len(matched_lines))
for elem in matched_lines:
    #print('Line Number = ', elem[0], ' :: Line = ', elem[1])
    #print(elem[0])
    elm1 = elem[0] -2
    print(elm1)
    #txt=str(elem[1])
    #txtf=txt[0:21] + '\n' + txt[21:42]
    
def replace_line(file_name, line_num, text): # to replace the count of atoms
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()
replace_line('topol.top', elm1, '#include "posre_lig.itp"\n' )
EOF
chmod 755 posre_inclusion.py
  python posre_inclusion.py
  rm -rf posre_inclusion.py
  gmx make_ndx -f EM.gro -o index.ndx <<EOF
1 | 13
q
EOF
  gmx grompp -f nvt.mdp -c E*.gro -r E*.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 4
  gmx mdrun -deffnm nvt
  gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 4
  gmx mdrun -deffnm npt
  gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_10.tpr -maxwarn 4
  source /usr/local/gromacs/bin/GMXRC
  gmx mdrun -deffnm md_10
}
no_reduction(){
  gmx trjconv -s md_*.tpr -f md_*.xtc -o md_10_center.xtc -center -pbc mol -ur compact << EOF
1 | 0
EOF
mm_pbsa_calculation
}
reduce_frames(){
  gmx trjconv -s md_*.tpr -f md_*.xtc -o md_10_center.xtc -center -pbc mol -ur compact -skip $factor << EOF
1 | 0
EOF
mm_pbsa_calculation
}
mm_pbsa_calculation(){
source /usr/local/gromacs/bin/GMXRC
gmx trjconv -s md_*.tpr -f *_center.xtc -o start_frame.pdb -dump 0 << EOF
0
EOF
cat << EOF > pbsa.mdp
;Polar calculation: "yes" or "no"
polar		= yes

;=============
;PSIZE options
;=============
;Factor by which to expand molecular dimensions to get coarsegrid dimensions.
cfac 		= 1.5

;The desired fine mesh spacing (in A)
gridspace 	= 0.5

:Amount (in A) to add to molecular dimensions to get fine grid dimensions.
fadd 		= 5

;Maximum memory (in MB) available per-processor for a calculation.
gmemceil 	= $memory

;=============================================
;APBS kwywords for polar solvation calculation
;=============================================
;Charge of positive ions
pcharge 	= 1

;Radius of positive charged ions
prad		= 0.95

;Concentration of positive charged ions
pconc           = 0.150 

;Charge of negative ions
ncharge 	= -1

;Radius of negative charged ions
nrad		= 1.81

;Concentration of negative charged ions
nconc 		= 0.150

;Solute dielectric constant
pdie 		= 2

;Solvent dielectric constant
sdie 		= 80

;Reference or vacuum dielectric constant
vdie 		= 1

;Solvent probe radius
srad 		= 1.4

;Method used to map biomolecular charges on grid. chgm = spl0 or spl2 or spl4
chgm            = spl4

;Model used to construct dielectric and ionic boundary. srfm = smol or spl2 or spl4
srfm            = smol

;Value for cubic spline window. Only used in case of srfm = spl2 or spl4.
swin 		= 0.30

;Numebr of grid point per A^2. Not used when (srad = 0.0) or (srfm = spl2 or spl4)
sdens 		= 10

;Temperature in K
temp 		= 300

;Type of boundary condition to solve PB equation. bcfl = zero or sdh or mdh or focus or map
bcfl 		= mdh

;Non-linear (npbe) or linear (lpbe) PB equation to solve
PBsolver 	= lpbe


;========================================================
;APBS kwywords for Apolar/Non-polar solvation calculation
;========================================================
;Non-polar solvation calculation: "yes" or "no"
apolar		= yes

;Repulsive contribution to Non-polar 
;===SASA model ====

;Gamma (Surface Tension) kJ/(mol A^2)
gamma           = 0.0226778

;Probe radius for SASA (A)
sasrad          = 1.4

;Offset (c) kJ/mol
sasaconst       = 3.84982

;===SAV model===
;Pressure kJ/(mol A^3)
press           = 0

;Probe radius for SAV (A)
savrad          = 0

;Offset (c) kJ/mol
savconst        = 0

;Attractive contribution to Non-polar
;===WCA model ====
;using WCA method: "yes" or "no"
WCA             = no

;Probe radius for WCA
wcarad          = 1.20

;bulk solvent density in A^3
bconc		= 0.033428

;displacment in A for surface area derivative calculation
dpos		= 0.05

;Quadrature grid points per A for molecular surface or solvent accessible surface
APsdens		= 20

;Quadrature grid spacing in A for volume integral calculations
grid            = 0.45 0.45 0.45

;Parameter to construct solvent related surface or volume
APsrfm          = sacc

;Cubic spline window in A for spline based surface definitions
APswin          = 0.3

;Temperature in K
APtemp          = 300
EOF
chmod 755 pbsa.mdp
g_mmpbsa -s md_*.tpr -f *_center.xtc -n index.ndx -i pbsa.mdp -pdie $dieC -pbsa -decomp << EOF
1
13
EOF
cat << EOF > stats.py
#!/usr/bin/python3
import re
import numpy as np
from scipy import stats
import argparse
import os
import math
import scipy.stats as spstat

def main():
    args = ParseOptions()
    CheckInput(args)
    #File => Frame wise component energy
    frame_wise = open(args.outfr, 'w')
    frame_wise.write('#Time E_VdW_mm(Protein)\tE_Elec_mm(Protein)\tE_Pol(Protein)\tE_Apol(Protein)\tE_VdW_mm(Ligand)\tE_Elec_mm(Ligand)\tE_Pol(Ligand)\tE_Apol(Ligand)\tE_VdW_mm(Complex)\tE_Elec_mm(Complex)\tE_Pol(Complex)\tE_Apol(Complex)\tDelta_E_mm\tDelta_E_Pol\tDelta_E_Apol\tDelta_E_binding\n')
    #Complex Energy
    c = []
    if args.multiple:
        MmFile, PolFile, APolFile = ReadMetafile(args.metafile)
        for i in range(len(MmFile)):
            cTmp = Complex(MmFile[i],PolFile[i],APolFile[i])
            cTmp.CalcEnergy(args,frame_wise,i)
            c.append(cTmp)
    else:
        cTmp = Complex(args.molmech,args.polar,args.apolar)
        cTmp.CalcEnergy(args,frame_wise,0)
        c.append(cTmp)
    #Summary in output files => "--outsum" and "--outmeta" file options
    Summary_Output_File(c, args)

class Complex():
    def __init__(self,MmFile,PolFile,APolFile):
        self.TotalEn = []
        self.Vdw, self.Elec, self.Pol, self.Sas, self.Sav, self.Wca =[], [], [], [], [], []
        self.MmFile = MmFile
        self.PolFile = PolFile
        self.APolFile = APolFile
        self.AvgEnBS = []
        self.CI = []
        self.FinalAvgEnergy = 0
        self.StdErr = 0
    
    def CalcEnergy(self,args,frame_wise,idx):
        mmEn = ReadData(self.MmFile,n=7)
        polEn = ReadData(self.PolFile,n=4)
        apolEn = ReadData(self.APolFile,n=10)
        CheckEnData(mmEn,polEn,apolEn)
    
        time, MM, Vdw, Elec, Pol, Apol, Sas, Sav, Wca = [], [], [], [], [], [], [], [], []    
        for i in range(len(mmEn[0])):
            #Vacuum MM
            Energy = mmEn[5][i] + mmEn[6][i] - (mmEn[1][i] + mmEn[2][i] + mmEn[3][i] + mmEn[4][i])
            MM.append(Energy)
            Energy = mmEn[5][i] - (mmEn[1][i] + mmEn[3][i])
            Vdw.append(Energy)
            Energy = mmEn[6][i] - (mmEn[2][i] + mmEn[4][i])
            Elec.append(Energy)
            # Polar
            Energy = polEn[3][i] - (polEn[1][i] + polEn[2][i])
            Pol.append(Energy)
            #Non-polar
            Energy = apolEn[3][i] + apolEn[6][i] + apolEn[9][i] - (apolEn[1][i] + apolEn[2][i] + apolEn[4][i] + apolEn[5][i] + apolEn[7][i] + apolEn[8][i])
            Apol.append(Energy)
            Energy = apolEn[3][i] - (apolEn[1][i] + apolEn[2][i])
            Sas.append(Energy)
            Energy = apolEn[6][i] - (apolEn[4][i] + apolEn[5][i])
            Sav.append(Energy)
            Energy = apolEn[9][i] - (apolEn[7][i] + apolEn[8][i])
            Wca.append(Energy)
            #Final Energy
            time.append(mmEn[0][i])
            Energy = MM[i] + Pol[i] + Apol[i]
            self.TotalEn.append(Energy)

        # Writing frame wise component energy to file
        frame_wise.write('\n#Complex %d\n' % ( (idx+1)))
        for i in range(len(time)):
            frame_wise.write('%15.3lf %15.3lf %15.3lf %15.3lf %15.3lf' % (time[i], mmEn[1][i], mmEn[2][i], polEn[1][i], (apolEn[1][i] + apolEn[4][i] + apolEn[7][i])))
            frame_wise.write('%15.3lf %15.3lf %15.3lf %15.3lf'         %          (mmEn[3][i], mmEn[4][i], polEn[2][i], (apolEn[2][i] + apolEn[5][i] + apolEn[8][i])))
            frame_wise.write('%15.3lf %15.3lf %15.3lf %15.3lf'         %          (mmEn[5][i], mmEn[6][i], polEn[3][i], (apolEn[3][i] + apolEn[6][i] + apolEn[9][i])))
            frame_wise.write('%15.3lf %15.3lf %15.3lf %15.3lf\n'         % (MM[i], Pol[i], Apol[i], self.TotalEn[i]))

        #Bootstrap analysis energy components
        if(args.bootstrap):
            bsteps = args.nbstep
            avg_energy, error = BootStrap(Vdw,bsteps)
            self.Vdw.append(avg_energy)
            self.Vdw.append(error)
            avg_energy, error = BootStrap(Elec,bsteps)
            self.Elec.append(avg_energy)
            self.Elec.append(error)
            avg_energy, error = BootStrap(Pol,bsteps)
            self.Pol.append(avg_energy)
            self.Pol.append(error)    
            avg_energy, error = BootStrap(Sas,bsteps)
            self.Sas.append(avg_energy)
            self.Sas.append(error)
            avg_energy, error = BootStrap(Sav,bsteps)
            self.Sav.append(avg_energy)
            self.Sav.append(error)
            avg_energy, error = BootStrap(Wca,bsteps)
            self.Wca.append(avg_energy)
            self.Wca.append(error)
            #Bootstrap => Final Average Energy
            self.AvgEnBS, AvgEn, EnErr, CI = ComplexBootStrap(self.TotalEn,bsteps)
            self.FinalAvgEnergy = AvgEn
            self.StdErr = EnErr
            self.CI = CI
        #If not bootstrap then average and standard deviation
        else:
            self.Vdw.append(np.mean(Vdw))
            self.Vdw.append(np.std(Vdw))
            self.Elec.append(np.mean(Elec))
            self.Elec.append(np.std(Elec))
            self.Pol.append(np.mean(Pol))
            self.Pol.append(np.std(Pol))
            self.Sas.append(np.mean(Sas))
            self.Sas.append(np.std(Sas))
            self.Sav.append(np.mean(Sav))
            self.Sav.append(np.std(Sav))
            self.Wca.append(np.mean(Wca))
            self.Wca.append(np.std(Wca))
            self.FinalAvgEnergy = np.mean(self.TotalEn)
            self.StdErr = np.std(self.TotalEn)
                

def Summary_Output_File(AllComplex,args):
        fs = open(args.outsum,'w')

        if args.multiple:
            fm = open(args.outmeta,'w')
            fm.write('# Complex_Number\t\tTotal_Binding_Energy\t\tError\n')
    
        for n in range(len(AllComplex)):
            fs.write('\n\n#Complex Number: %4d\n' % (n+1))    
            fs.write('===============\n   SUMMARY   \n===============\n\n')
            fs.write('\n van der Waal energy      = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].Vdw[0], AllComplex[n].Vdw[1]))
            fs.write('\n Electrostattic energy    = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].Elec[0],AllComplex[n].Elec[1]))
            fs.write('\n Polar solvation energy   = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].Pol[0], AllComplex[n].Pol[1]))
            fs.write('\n SASA energy              = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].Sas[0], AllComplex[n].Sas[1]))
            fs.write('\n SAV energy               = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].Sav[0], AllComplex[n].Sav[1]))
            fs.write('\n WCA energy               = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].Wca[0], AllComplex[n].Wca[1]))
            fs.write('\n Binding energy           = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].FinalAvgEnergy, AllComplex[n].StdErr))
            fs.write('\n===============\n    END     \n===============\n\n')

        if args.multiple:
                fm.write('%5d %15.3lf %7.3lf\n' % (n+1 , AllComplex[n].FinalAvgEnergy, AllComplex[n].StdErr))

def CheckEnData(mmEn,polEn,apolEn):
    frame = len(mmEn[0])
    for i in range(len(mmEn)):
        if(len(mmEn[i]) != frame):
            print ("Number of row is not same for all column")
            exit(1)
    for i in range(len(polEn)):
        if(len(polEn[i]) != frame):
            print ("Number of row is not same for all column")
            exit(1)
    for i in range(len(polEn)):
        if(len(polEn[i]) != frame):
            print ("Number of row is not same for all column")
            exit(1)

def CheckInput(args):
    if args.multiple:
        if not os.path.exists(args.metafile):
            print ('\n{0} not found....\n') .format(args.metafile)
            exit(1)
    else:
        if not os.path.exists(args.molmech):
            print ('\n{0} not found....\n') .format(args.molmech)
            exit(1)
        if not os.path.exists(args.polar):
            print ('\n{0} not found....\n') .format(args.polar)
            exit(1)
        if not os.path.exists(args.apolar):
            print ('\n{0} not found....\n') .format(args.polar)
            exit(1)
    
def ParseOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument("-mt", "--multiple", help='If given, calculate for multiple complexes. Need Metafile containing path of energy files', action="store_true")
    parser.add_argument("-mf", "--metafile", help='Metafile containing path to energy files of each complex in a row obtained from g_mmpbsa in following order: \
                                                       [MM file] [Polar file] [ Non-polar file] ',action="store", default='metafile.dat', metavar='metafile.dat')
    parser.add_argument("-m", "--molmech", help='Vacuum Molecular Mechanics energy file obtained from g_mmpbsa',action="store", default='energy_MM.xvg', metavar='energy_MM.xvg')
    parser.add_argument("-p", "--polar", help='Polar solvation energy file obtained from g_mmpbsa',action="store",default='polar.xvg', metavar='polar.xvg')
    parser.add_argument("-a", "--apolar", help='Non-Polar solvation energy file obtained from g_mmpbsa',action="store",default='apolar.xvg',metavar='polar.xvg')
    parser.add_argument("-bs", "--bootstrap", help='If given, Enable Boot Strap analysis',action="store_true")
    parser.add_argument("-nbs", "--nbstep", help='Number of boot strap steps for average energy calculation',action="store", type=int, default=500, metavar=500)
    parser.add_argument("-of", "--outfr", help='Energy File: All energy components frame wise',action="store",default='full_energy.dat', metavar='full_energy.dat')
    parser.add_argument("-os", "--outsum", help='Final Energy File: Full Summary of energy components',action="store",default='summary_energy.dat', metavar='summary_energy.dat')
    parser.add_argument("-om", "--outmeta", help='Final Energy File for Multiple Complexes: Complex wise final binding nergy',action="store",default='meta_energy.dat',metavar='meta_energy.dat')
    args = parser.parse_args()
    return args

def ReadData(FileName,n=2):
    infile = open(FileName,'r')
    x, data = [],[]
    for line in infile:
        line = line.rstrip('\n')
        if not line.strip():
            continue
        if(re.match('#|@',line)==None):
            temp = line.split()
            data.append(np.array(temp))
    for j in range(0,n):
        x_temp =[]
        for i in range(len(data)):
            x_temp.append(float(data[i][j]))
        x.append(x_temp)
    return x

def ComplexBootStrap(x,step=1000):
    avg =[]
    x = np.array(x)
    n = len(x)
    idx = np.random.randint(0,n,(step,n))
    sample_x = x[idx]
    avg = np.sort(np.mean(sample_x,1))
    CI_min = avg[int(0.005*step)]
    CI_max = avg[int(0.995*step)]
    #print('Energy = %13.3f; Confidance Interval = (-%-5.3f / +%-5.3f)\n' % (np.mean(avg), (np.mean(avg)-CI_min), (CI_max-np.mean(avg))))
    return avg, np.mean(avg), np.std(avg), [(np.mean(avg)-CI_min), (CI_max-np.mean(avg))]

def BootStrap (x,step=1000):
    if(np.mean(x)) == 0:
        return 0.000, 0.000
    else:
        avg =[]
        x = np.array(x)
        n = len(x)
        idx = np.random.randint(0,n,(step,n))
        sample_x = x[idx]
        avg = np.sort(np.mean(sample_x,1))
        return np.mean(avg),np.std(avg)

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def ReadMetafile(metafile):
    MmFile,PolFile, APolFile = [], [], []
    FileList = open(metafile,'r')
    for line in FileList:
        line = line.rstrip('\n')
        if not line.strip():
            continue
        temp = line.split()    
        MmFile.append(temp[0])
        PolFile.append(temp[1])
        APolFile.append(temp[2])
        
        
        if not os.path.exists(temp[0]):
            print ('\n{0} not found....\n') .format(temp[0])
            exit(1)
        if not os.path.exists(temp[1]):
            print ('\n{0} not found....\n') .format(temp[1])
            exit(1)
        if not os.path.exists(temp[2]):
            print ('\n{0} not found....\n') .format(temp[2])
            exit(1)
                        
    return MmFile, PolFile, APolFile

if __name__=="__main__":
    main()
EOF
chmod 755 stats.py
python3 stats.py -bs -nbs 2000 -m energy_MM.xvg -p polar.xvg -a apolar.xvg
sed '1d;2d;3d' full_energy.dat | sed '1i\,Time(Corresponding_frame),E_VdW_mm(Protein),E_Elec_mm(Protein),E_Pol(Protein),E_Apol(Protein),E_VdW_mm(Ligand),E_Elec_mm(Ligand),E_Pol(Ligand),E_Apol(Ligand),E_VdW_mm(Complex),E_Elec_mm(Complex),E_Pol(Complex),E_Apol(Complex),Delta_E_mm,Delta_E_Pol,Delta_E_Apol,Delta_E_binding' |  sed "s/ \{1,\}/,/g" > full_energy_corrected.csv 
sed -i 's/^,//' full_energy_corrected.csv
rm -rf full_energy.dat stats.py
cat << EOF > csv_plot_fe.py
#!/usr/bin/python3
import os, sys, matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import io
import pandas as pd
import csv
df = pd.read_csv("full_energy_corrected.csv",index_col=False)
columns_name = list(df.columns)
mean = df[columns_name[16]].mean()
sd = df[columns_name[16]].std()
sd_string=str(round(sd, 2))
mean_string=str(round(mean, 2))
text = '(nm)' + ' "Mean +/- Sd (' + mean_string + ' +/- ' + sd_string + ')"'
df.plot(kind = 'line', x = columns_name[0], y = columns_name[13:17])
plt.legend()
#x = round(df[columns_name[1]].max(),2)
#plt.text(3.5, x,text, fontsize = 10, color = 'b')
plt.grid(visible=True, which='major', color='grey', linestyle='--')
plt.ylabel('∆G binding free energy (kJ/mol)')
file_name = 'full_energy_corrected' + '.tiff' 
plt.xlabel(columns_name[0])
plt.suptitle('MM/PBSA-based ∆G Calculation')
plt.title('MM/PBSA-base ∆G')
plt.savefig(file_name, dpi=300)
EOF
chmod 755 csv_plot_fe.py
python3 csv_plot_fe.py



cat << EOF > decomp.py
#!/usr/bin/python3
import re
import numpy as np
import argparse
import sys
import os
import math

def main():
	args = ParseOptions()
	CheckInput(args)
	MMEnData,resnameA = ReadData(args.molmech)
	polEnData,resnameB = ReadData(args.polar)
	apolEnData,resnameC = ReadData(args.apolar)
	resname = CheckResname(resnameA,resnameB,resnameC)
	print("Total number of Residue: {0}\n" . format(len('resname')+1))
	Residues = []
	for i in range(len(resname)):
		CheckEnData(MMEnData[i],polEnData[i],apolEnData[i])
		r = Residue()
		r.CalcEnergy(MMEnData[i],polEnData[i],apolEnData[i],args)
		Residues.append(r)
		print(' %8s %8.4f %8.4f' % (resname[i], r.TotalEn[0], r.TotalEn[1]))
	fout = open(args.output,'w')
	fmap = open(args.outmap,'w')
	fout.write('#Residues  MM Energy(+/-)dev/error  Polar Energy(+/-)dev/error APolar Energy(+/-)dev/error Total Energy(+/-)dev/error\n')
	for i in range(len(resname)):
		if (args.cutoff == 999):
			fout.write("%-8s  %4.4f  %4.4f    %4.4f  %4.4f    %4.4f  %4.4f    %4.4f  %4.4f \n" %(resname[i],Residues[i].FinalMM[0],Residues[i].FinalMM[1], 
                                      Residues[i].FinalPol[0], Residues[i].FinalPol[1], Residues[i].FinalAPol[0], Residues[i].FinalAPol[1], Residues[i].TotalEn[0],Residues[i].TotalEn[1] ))
		elif (args.cutoff <= Residues[i].TotalEn[0]) or ( (-1 *args.cutoff) >= Residues[i].TotalEn[0]):
			fout.write("%-8s  %4.4f  %4.4f    %4.4f  %4.4f    %4.4f  %4.4f    %4.4f  %4.4f \n" % (resname[i],Residues[i].FinalMM[0],Residues[i].FinalMM[1], 
                                  Residues[i].FinalPol[0], Residues[i].FinalPol[1], Residues[i].FinalAPol[0], Residues[i].FinalAPol[1], Residues[i].TotalEn[0],Residues[i].TotalEn[1] ))
			
		fmap.write("%-8d     %4.4f \n" %((i+1), Residues[i].TotalEn[0]))


class Residue():
	def __init__(self):
		self.FinalMM, self.FinalPol, self.FinalAPol, self.TotalEn = [], [], [], []
		
	def BootStrap(self,x,step):
		avg =[]
		x = np.array(x)
		n = len(x)
		idx = np.random.randint(0,n,(step,n))
		sample_x = x[idx]
		avg = np.sort(np.mean(sample_x,1))
		return np.mean(avg), np.std(avg)

	def CalcEnergy(self,MM,Pol,APol,args):
		TotalEn = np.sum([MM,Pol,APol],axis=0)
		if(args.bootstrap):
			self.FinalMM = self.BootStrap(MM,args.nbstep)
			self.FinalPol = self.BootStrap(Pol,args.nbstep)
			if (np.mean(APol) == 0):
				self.FinalAPol = [0.0,0.0]
			else:
				self.FinalAPol = self.BootStrap(APol,args.nbstep)
			self.TotalEn = self.BootStrap(TotalEn,args.nbstep)
		else:
			self.FinalMM = [np.mean(MM),np.std(MM)]
			self.FinalPol = [np.mean(Pol),np.std(Pol)]
			if (np.mean(APol) == 0):
				self.FinalAPol = [0.0,0.0]
			else:
				self.FinalAPol = [np.mean(APol),np.std(APol)]
			self.TotalEn = [np.mean(TotalEn),np.std(TotalEn)]
		self.FinalMM = np.round(self.FinalMM,4)
		self.FinalPol = np.round(self.FinalPol,4)
		self.FinalAPol = np.round(self.FinalAPol,4)
		self.TotalEn = np.round(self.TotalEn,4)

def CheckEnData(MM,Pol,APol):
	if(len(Pol) != len(MM)):
		print("Times or Frames Mismatch between files")
		exit(1)
	if(len(APol) != len(Pol)):
		print("Times or Frames Mismatch between files")
		exit(1)
	if(len(APol) != len(MM)):
		print("Times or Frames Mismatch between files")
		exit(1)


def ParseOptions():
        parser = argparse.ArgumentParser()
        parser.add_argument("-m", "--molmech", help='Molecular Mechanics energy file',action="store", default='contrib_MM.dat', metavar='contrib_MM.dat')
        parser.add_argument("-p", "--polar", help='Polar solvation energy file',action="store",default='contrib_pol.dat', metavar='contrib_pol.dat')
        parser.add_argument("-a", "--apolar", help='Non-Polar solvation energy file',action="store",default='contrib_apol.dat',metavar='contrib_apol.dat')
        parser.add_argument("-bs", "--bootstrap", help='Switch for Error by Boot Strap analysis',action="store_true")
        parser.add_argument("-nbs", "--nbstep", help='Number of boot strap steps',action="store", type=int,default=500, metavar=500)
        parser.add_argument("-ct", "--cutoff", help='Absolute Cutoff: energy output above and below this value',action="store",type=float,default=999, metavar=999)
        parser.add_argument("-o", "--output", help='Final Decomposed Energy File',action="store",default='final_contrib_energy.dat', metavar='final_contrib_energy.dat')
        parser.add_argument("-om", "--outmap", help='energy2bfac input file: to map energy on structure for visualization',action="store",default='energyMapIn.dat', metavar='energyMapIn.dat')
        return parser.parse_args()

def CheckResname(resA,resB,resC):
	if(len(resA) != len(resB)):
		print("ERROR: Total number of residues mismatch between files")
		exit(1)
	if(len(resB) != len(resC)):
		print("ERROR: Total number of residues mismatch between files")
		exit(1)
	if(len(resC) != len(resA)):
		print("ERROR: Total number of residues mismatch between files")
		exit(1)
	for i in range(len(resA)):
		if (resA[i] != resB[i]):
			print("ERROR: Residue mismatch between files")
			exit(1)
	for i in range(len(resB)):
		if (resB[i] != resC[i]):
			print("ERROR: Residue mismatch between files")
			exit(1)
	for i in range(len(resA)):
		if (resA[i] != resC[i]):
			print("ERROR: Residue mismatch between files")
			exit(1)
	return resA

def CheckInput(args):
	if not os.path.exists(args.molmech):
		print('\n{0} not found....\n' .format(args.molmech))
		exit(1)
	if not os.path.exists(args.polar):
		print('\n{0} not found....\n' .format(args.polar))
		exit(1)
	if not os.path.exists(args.apolar):
		print('\n{0} not found....\n' .format(args.polar))
		exit(1)


def ReadData(FileName):
        infile = open(FileName,'r')
        x, data,resname = [],[],[]
        for line in infile:
                line = line.rstrip('\n')
                if not line.strip():
                        continue
                if(re.match('#|@',line)==None):
                        temp = line.split()
                        data.append(np.array(temp))
                if(re.match('#',line)):
                        resname = line.split()
                        n = len(resname[1:])
        for j in range(1,n):
                x_temp =[]
                for i in range(len(data)):
                        x_temp.append(float(data[i][j]))
                x.append(x_temp)
        return x, resname[2:]

if __name__=="__main__":
        main()
EOF
chmod 755 decomp.py
python3 decomp.py -bs -nbs 2000 -m contrib_MM.dat -p contrib_pol.dat -a contrib_apol.dat
sed '1d' final_contrib_energy.dat | sed '1i\Residues,MM_Energy,dev-error,Polar_Energy,dev-error,APolar_Energy,dev-error,Total_Energy,dev-error' | sed "s/ \{1,\}/,/g" > final_contrib_energy_corrected.csv
sed -i 's/^,//' final_contrib_energy_corrected.csv
rm -rf final_contrib_energy.dat decomp.py
cat << EOF > csv_plot.py
#!/usr/bin/python3
import os, sys, matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import io
import pandas as pd
import csv
df = pd.read_csv("final_contrib_energy_corrected.csv",index_col=False)
columns_name = list(df.columns)
#mean = df[columns_name[7]].mean()
#sd = df[columns_name[7]].std()
#sd_string=str(round(sd, 2))
#mean_string=str(round(mean, 2))
#text = '(nm)' + ' "Mean +/- Sd (' + mean_string + ' +/- ' + sd_string + ')"'
df.plot(kind = 'line', x = columns_name[0], y = columns_name[7])
plt.legend()
#x = round(df[columns_name[1]].max(),2)
#plt.text(3.5, x,text, fontsize = 10, color = 'b')
plt.grid(visible=True, which='major', color='grey', linestyle='--')
plt.ylabel('∆G binding free energy (kJ/mol)')
file_name = 'final_contrib_energy_corrected' + '.tiff' 
plt.xlabel(columns_name[0])
plt.suptitle('MM/PBSA-base ∆G calculation(P-L complex)')
plt.title('Energy Decomposition Analysis')
plt.savefig(file_name, dpi=300)
EOF
chmod 755 csv_plot.py
python3 csv_plot.py
}
function output_md(){
  mkdir output_md
  mv *.itp *.top *.pdb *.gro *.mdp *.tpr *.xtc *.cpt *.log *.ndx *.edr *.trr *.csv output_md/
  tar -czf $fname.tar.gz output_md/
  rm -rf *.itp *.top *.pdb *.gro *.mdp *.tpr *.cpt *.log *.ndx *.edr *.trr *.xtc output_md/
  zenity --info --title="Successfull!!!" --text="Hey buddy you have done it. All the files are in $fname.tar.gz" --width 300 --height 100
}
function output_mmpbsa(){
  mkdir output_mmpbsa
  mv *.itp *.top *.pdb *.gro *.mdp *.tpr *.xtc *.cpt *.log *.ndx *.edr *.trr *.xvg *.dat *.csv *.tiff *.py output_mmpbsa/
  tar -czf $fname.tar.gz output_mmpbsa/
  rm -rf *.itp *.top *.pdb *.gro *.mdp *.tpr *.xtc *.cpt *.log *.ndx *.edr *.trr *.dat *.tiff output_mmpbsa/
  zenity --info --title="Successfull!!!" --text="Hey buddy you have done it. All the files are in $fname.tar.gz" --width 300 --height 100
}
##################################################################################################################################
##################################################################################################################################
#######Protein-Protein complex simulation###################################Protein in water######################################
#########################################P-P complex MD and MMPBSA calculation####################################################
##################################################################################################################################
##################################################################################################################################
function simulate_protein(){
em=$(cat << EOF > em.mdp
;
;	GROMACS
;	Energy Minimization Script
;
;
define		= -DFLEXIBLE	; pass to preprocessor
cpp		= usr/bin/cpp	; location of preprocessor
constraints	= none
integrator	= steep		; steepest decents minimum (else cg)
nsteps		= $Energy_minimization
;
;	Energy Minimizing Stuff
;
emtol		= 10  		; convergence total force(kJ/mol/nm) is smaller than
emstep		= 0.01		; initial step size (nm)
nstcomm		= 100		; frequency or COM motion removal
ns_type		= grid
rlist		= 1.4		; cut-off distance for short range neighbors
rcoulomb	= 1.1		; distance for coulomb cut-off
coulombtype	= PME		; electrostatics (Particle Mesh Ewald method)
fourierspacing	= 0.12		; max grid spacing when using PPPM or PME
vdw-type	= Shift
rvdw		= 1.1		; Verlet cutoff
Tcoupl		= no		; temperature coupling
Pcoupl		= no		; pressure coupling
gen_vel		= no
EOF
)
nvt=$(cat << EOF > nvt.mdp
title		= OPLS Lysozyme NVT equilibration 
define		= -DPOSRES	; position restrain the protein
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= 50000		; 2 * 50000 = 100 ps
dt		    = 0.002		; 2 fs
; Output control
nstxout		= 500		; save coordinates every 1.0 ps
nstvout		= 500		; save velocities every 1.0 ps
nstenergy	= 500		; save energies every 1.0 ps
nstlog		= 500		; update log file every 1.0 ps
; Bond parameters
continuation	        = no		; first dynamics run
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10		; 20 fs, largely irrelevant with Verlet
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME	; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		; cubic interpolation
fourierspacing	= 0.16	; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1           ; time constant, in ps
ref_t		= 300 	  300           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl		= no 		; no pressure coupling in NVT
; Periodic boundary conditions
pbc		= xyz		    ; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= yes		; assign velocities from Maxwell distribution
gen_temp	= 300		; temperature for Maxwell distribution
gen_seed	= -1		; generate a random seed
EOF
)
npt=$(cat << EOF > npt.mdp
title		= OPLS Lysozyme NPT equilibration 
define		= -DPOSRES	; position restrain the protein
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= 50000		; 2 * 50000 = 100 ps
dt		    = 0.002		; 2 fs
; Output control
nstxout		= 500		; save coordinates every 1.0 ps
nstvout		= 500		; save velocities every 1.0 ps
nstenergy	= 500		; save energies every 1.0 ps
nstlog		= 500		; update log file every 1.0 ps
; Bond parameters
continuation	        = yes		; Restarting after NVT 
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		    ; cubic interpolation
fourierspacing	= 0.16		; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1	        ; time constant, in ps
ref_t		= 300 	  300	        ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT
pcoupltype	        = isotropic	            ; uniform scaling of box vectors
tau_p		        = 2.0		            ; time constant, in ps
ref_p		        = 1.0		            ; reference pressure, in bar
compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com
; Periodic boundary conditions
pbc		= xyz		; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= no		; Velocity generation is off 
EOF
)

mdp=$(cat << EOF > md.mdp
title		= OPLS Lysozyme MD simulation 
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= $nsteps     ; 2 * 50000000 = 100000 ps (200 ns)
dt		= 0.002	; 2 fs
; Output control
nstxout		        = 0		; save coordinates every 10.0 ps
nstvout		        = 0		; save velocities every 10.0 ps
nstenergy	        = 500000		; save energies every 10.0 ps
nstlog		        = 500000		; update log file every 10.0 ps
nstxout-compressed  = 500000      ; save compressed coordinates every 10.0 ps
                                ; nstxout-compressed replaces nstxtcout
compressed-x-grps   = System    ; replaces xtc-grps
; Bond parameters
continuation	        = yes		; Restarting after NPT 
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		    ; cubic interpolation
fourierspacing	= 0.16		; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1	        ; time constant, in ps
ref_t		= 300 	  300	        ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT
pcoupltype	        = isotropic	            ; uniform scaling of box vectors
tau_p		        = 2.0		            ; time constant, in ps
ref_p		        = 1.0		            ; reference pressure, in bar
compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc		= xyz		; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= no		; Velocity generation is off
EOF
)





  gmx editconf -f conf.gro -d 1.0 -bt triclinic -o box.gro
  gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o box_sol.gro 
  gmx grompp -f em.mdp -c box_sol.gro -p topol.top -o ION.tpr -maxwarn 4
  gmx genion -s ION.tpr -p topol.top -conc 0.1 -neutral -o box_sol_ion.gro << EOF
13
EOF
  gmx grompp -f em.mdp -c box_sol_ion.gro -p topol.top -o EM.tpr -maxwarn 4
  gmx mdrun -v -deffnm EM
    case "${go_no_go}" in
  "Simulate Protein-Protein (P-P complex)")
  printf "splitch 1 & 2\nq\n"|gmx make_ndx -f EM.gro -o index.ndx;;
  "Protein in water")
  printf "1\nq\n"|gmx make_ndx -f EM.gro -o index.ndx;;
  "P-P complex MD & Thermal MM/PBSA calculation")
  printf "splitch 1 & 2\nq\n"|gmx make_ndx -f EM.gro -o index.ndx;;
esac
  gmx grompp -f nvt.mdp -c E*.gro -r E*.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 4
  gmx mdrun -deffnm nvt
  gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 4
  gmx mdrun -deffnm npt
  gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_10.tpr -maxwarn 4
  source /usr/local/gromacs/bin/GMXRC
  gmx mdrun -deffnm md_10
}
no_reduction_protein(){
  gmx trjconv -s md_*.tpr -f md_*.xtc -o md_10_center.xtc -center -pbc mol -ur compact << EOF
1 | 0
EOF
mm_pbsa_calculation_protein
}
reduce_frames_protein(){
  gmx trjconv -s md_*.tpr -f md_*.xtc -o md_10_center.xtc -center -pbc mol -ur compact -skip $factor << EOF
1 | 0
EOF
mm_pbsa_calculation_protein
}
mm_pbsa_calculation_protein(){
gmx trjconv -s md_*.tpr -f *_center.xtc -o start_frame.pdb -dump 0 << EOF
0
EOF
cat << EOF > pbsa.mdp
;Polar calculation: "yes" or "no"
polar		= yes

;=============
;PSIZE options
;=============
;Factor by which to expand molecular dimensions to get coarsegrid dimensions.
cfac 		= 1.5

;The desired fine mesh spacing (in A)
gridspace 	= 0.5

:Amount (in A) to add to molecular dimensions to get fine grid dimensions.
fadd 		= 5

;Maximum memory (in MB) available per-processor for a calculation.
gmemceil 	= $memory

;=============================================
;APBS kwywords for polar solvation calculation
;=============================================
;Charge of positive ions
pcharge 	= 1

;Radius of positive charged ions
prad		= 0.95

;Concentration of positive charged ions
pconc           = 0.150 

;Charge of negative ions
ncharge 	= -1

;Radius of negative charged ions
nrad		= 1.81

;Concentration of negative charged ions
nconc 		= 0.150

;Solute dielectric constant
pdie 		= 2

;Solvent dielectric constant
sdie 		= 80

;Reference or vacuum dielectric constant
vdie 		= 1

;Solvent probe radius
srad 		= 1.4

;Method used to map biomolecular charges on grid. chgm = spl0 or spl2 or spl4
chgm            = spl4

;Model used to construct dielectric and ionic boundary. srfm = smol or spl2 or spl4
srfm            = smol

;Value for cubic spline window. Only used in case of srfm = spl2 or spl4.
swin 		= 0.30

;Numebr of grid point per A^2. Not used when (srad = 0.0) or (srfm = spl2 or spl4)
sdens 		= 10

;Temperature in K
temp 		= 300

;Type of boundary condition to solve PB equation. bcfl = zero or sdh or mdh or focus or map
bcfl 		= mdh

;Non-linear (npbe) or linear (lpbe) PB equation to solve
PBsolver 	= lpbe


;========================================================
;APBS kwywords for Apolar/Non-polar solvation calculation
;========================================================
;Non-polar solvation calculation: "yes" or "no"
apolar		= yes

;Repulsive contribution to Non-polar 
;===SASA model ====

;Gamma (Surface Tension) kJ/(mol A^2)
gamma           = 0.0226778

;Probe radius for SASA (A)
sasrad          = 1.4

;Offset (c) kJ/mol
sasaconst       = 3.84982

;===SAV model===
;Pressure kJ/(mol A^3)
press           = 0

;Probe radius for SAV (A)
savrad          = 0

;Offset (c) kJ/mol
savconst        = 0

;Attractive contribution to Non-polar
;===WCA model ====
;using WCA method: "yes" or "no"
WCA             = no

;Probe radius for WCA
wcarad          = 1.20

;bulk solvent density in A^3
bconc		= 0.033428

;displacment in A for surface area derivative calculation
dpos		= 0.05

;Quadrature grid points per A for molecular surface or solvent accessible surface
APsdens		= 20

;Quadrature grid spacing in A for volume integral calculations
grid            = 0.45 0.45 0.45

;Parameter to construct solvent related surface or volume
APsrfm          = sacc

;Cubic spline window in A for spline based surface definitions
APswin          = 0.3

;Temperature in K
APtemp          = 300
EOF
printf "Protein_chain1\nProtein_chain2\n" |g_mmpbsa -s md_*.tpr -f *_center.xtc -n index.ndx -i pbsa.mdp -pdie $dieC -pbsa -decomp
cat << EOF > stats.py
#!/usr/bin/python3
import re
import numpy as np
from scipy import stats
import argparse
import os
import math
import scipy.stats as spstat

def main():
    args = ParseOptions()
    CheckInput(args)
    #File => Frame wise component energy
    frame_wise = open(args.outfr, 'w')
    frame_wise.write('#Time E_VdW_mm(Protein)\tE_Elec_mm(Protein)\tE_Pol(Protein)\tE_Apol(Protein)\tE_VdW_mm(Ligand)\tE_Elec_mm(Ligand)\tE_Pol(Ligand)\tE_Apol(Ligand)\tE_VdW_mm(Complex)\tE_Elec_mm(Complex)\tE_Pol(Complex)\tE_Apol(Complex)\tDelta_E_mm\tDelta_E_Pol\tDelta_E_Apol\tDelta_E_binding\n')
    #Complex Energy
    c = []
    if args.multiple:
        MmFile, PolFile, APolFile = ReadMetafile(args.metafile)
        for i in range(len(MmFile)):
            cTmp = Complex(MmFile[i],PolFile[i],APolFile[i])
            cTmp.CalcEnergy(args,frame_wise,i)
            c.append(cTmp)
    else:
        cTmp = Complex(args.molmech,args.polar,args.apolar)
        cTmp.CalcEnergy(args,frame_wise,0)
        c.append(cTmp)
    #Summary in output files => "--outsum" and "--outmeta" file options
    Summary_Output_File(c, args)

class Complex():
    def __init__(self,MmFile,PolFile,APolFile):
        self.TotalEn = []
        self.Vdw, self.Elec, self.Pol, self.Sas, self.Sav, self.Wca =[], [], [], [], [], []
        self.MmFile = MmFile
        self.PolFile = PolFile
        self.APolFile = APolFile
        self.AvgEnBS = []
        self.CI = []
        self.FinalAvgEnergy = 0
        self.StdErr = 0
    
    def CalcEnergy(self,args,frame_wise,idx):
        mmEn = ReadData(self.MmFile,n=7)
        polEn = ReadData(self.PolFile,n=4)
        apolEn = ReadData(self.APolFile,n=10)
        CheckEnData(mmEn,polEn,apolEn)
    
        time, MM, Vdw, Elec, Pol, Apol, Sas, Sav, Wca = [], [], [], [], [], [], [], [], []    
        for i in range(len(mmEn[0])):
            #Vacuum MM
            Energy = mmEn[5][i] + mmEn[6][i] - (mmEn[1][i] + mmEn[2][i] + mmEn[3][i] + mmEn[4][i])
            MM.append(Energy)
            Energy = mmEn[5][i] - (mmEn[1][i] + mmEn[3][i])
            Vdw.append(Energy)
            Energy = mmEn[6][i] - (mmEn[2][i] + mmEn[4][i])
            Elec.append(Energy)
            # Polar
            Energy = polEn[3][i] - (polEn[1][i] + polEn[2][i])
            Pol.append(Energy)
            #Non-polar
            Energy = apolEn[3][i] + apolEn[6][i] + apolEn[9][i] - (apolEn[1][i] + apolEn[2][i] + apolEn[4][i] + apolEn[5][i] + apolEn[7][i] + apolEn[8][i])
            Apol.append(Energy)
            Energy = apolEn[3][i] - (apolEn[1][i] + apolEn[2][i])
            Sas.append(Energy)
            Energy = apolEn[6][i] - (apolEn[4][i] + apolEn[5][i])
            Sav.append(Energy)
            Energy = apolEn[9][i] - (apolEn[7][i] + apolEn[8][i])
            Wca.append(Energy)
            #Final Energy
            time.append(mmEn[0][i])
            Energy = MM[i] + Pol[i] + Apol[i]
            self.TotalEn.append(Energy)

        # Writing frame wise component energy to file
        frame_wise.write('\n#Complex %d\n' % ( (idx+1)))
        for i in range(len(time)):
            frame_wise.write('%15.3lf %15.3lf %15.3lf %15.3lf %15.3lf' % (time[i], mmEn[1][i], mmEn[2][i], polEn[1][i], (apolEn[1][i] + apolEn[4][i] + apolEn[7][i])))
            frame_wise.write('%15.3lf %15.3lf %15.3lf %15.3lf'         %          (mmEn[3][i], mmEn[4][i], polEn[2][i], (apolEn[2][i] + apolEn[5][i] + apolEn[8][i])))
            frame_wise.write('%15.3lf %15.3lf %15.3lf %15.3lf'         %          (mmEn[5][i], mmEn[6][i], polEn[3][i], (apolEn[3][i] + apolEn[6][i] + apolEn[9][i])))
            frame_wise.write('%15.3lf %15.3lf %15.3lf %15.3lf\n'         % (MM[i], Pol[i], Apol[i], self.TotalEn[i]))

        #Bootstrap analysis energy components
        if(args.bootstrap):
            bsteps = args.nbstep
            avg_energy, error = BootStrap(Vdw,bsteps)
            self.Vdw.append(avg_energy)
            self.Vdw.append(error)
            avg_energy, error = BootStrap(Elec,bsteps)
            self.Elec.append(avg_energy)
            self.Elec.append(error)
            avg_energy, error = BootStrap(Pol,bsteps)
            self.Pol.append(avg_energy)
            self.Pol.append(error)    
            avg_energy, error = BootStrap(Sas,bsteps)
            self.Sas.append(avg_energy)
            self.Sas.append(error)
            avg_energy, error = BootStrap(Sav,bsteps)
            self.Sav.append(avg_energy)
            self.Sav.append(error)
            avg_energy, error = BootStrap(Wca,bsteps)
            self.Wca.append(avg_energy)
            self.Wca.append(error)
            #Bootstrap => Final Average Energy
            self.AvgEnBS, AvgEn, EnErr, CI = ComplexBootStrap(self.TotalEn,bsteps)
            self.FinalAvgEnergy = AvgEn
            self.StdErr = EnErr
            self.CI = CI
        #If not bootstrap then average and standard deviation
        else:
            self.Vdw.append(np.mean(Vdw))
            self.Vdw.append(np.std(Vdw))
            self.Elec.append(np.mean(Elec))
            self.Elec.append(np.std(Elec))
            self.Pol.append(np.mean(Pol))
            self.Pol.append(np.std(Pol))
            self.Sas.append(np.mean(Sas))
            self.Sas.append(np.std(Sas))
            self.Sav.append(np.mean(Sav))
            self.Sav.append(np.std(Sav))
            self.Wca.append(np.mean(Wca))
            self.Wca.append(np.std(Wca))
            self.FinalAvgEnergy = np.mean(self.TotalEn)
            self.StdErr = np.std(self.TotalEn)
                

def Summary_Output_File(AllComplex,args):
        fs = open(args.outsum,'w')

        if args.multiple:
            fm = open(args.outmeta,'w')
            fm.write('# Complex_Number\t\tTotal_Binding_Energy\t\tError\n')
    
        for n in range(len(AllComplex)):
            fs.write('\n\n#Complex Number: %4d\n' % (n+1))    
            fs.write('===============\n   SUMMARY   \n===============\n\n')
            fs.write('\n van der Waal energy      = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].Vdw[0], AllComplex[n].Vdw[1]))
            fs.write('\n Electrostattic energy    = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].Elec[0],AllComplex[n].Elec[1]))
            fs.write('\n Polar solvation energy   = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].Pol[0], AllComplex[n].Pol[1]))
            fs.write('\n SASA energy              = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].Sas[0], AllComplex[n].Sas[1]))
            fs.write('\n SAV energy               = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].Sav[0], AllComplex[n].Sav[1]))
            fs.write('\n WCA energy               = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].Wca[0], AllComplex[n].Wca[1]))
            fs.write('\n Binding energy           = %15.3lf   +/-  %7.3lf kJ/mol\n' % (AllComplex[n].FinalAvgEnergy, AllComplex[n].StdErr))
            fs.write('\n===============\n    END     \n===============\n\n')

        if args.multiple:
                fm.write('%5d %15.3lf %7.3lf\n' % (n+1 , AllComplex[n].FinalAvgEnergy, AllComplex[n].StdErr))

def CheckEnData(mmEn,polEn,apolEn):
    frame = len(mmEn[0])
    for i in range(len(mmEn)):
        if(len(mmEn[i]) != frame):
            print ("Number of row is not same for all column")
            exit(1)
    for i in range(len(polEn)):
        if(len(polEn[i]) != frame):
            print ("Number of row is not same for all column")
            exit(1)
    for i in range(len(polEn)):
        if(len(polEn[i]) != frame):
            print ("Number of row is not same for all column")
            exit(1)

def CheckInput(args):
    if args.multiple:
        if not os.path.exists(args.metafile):
            print ('\n{0} not found....\n') .format(args.metafile)
            exit(1)
    else:
        if not os.path.exists(args.molmech):
            print ('\n{0} not found....\n') .format(args.molmech)
            exit(1)
        if not os.path.exists(args.polar):
            print ('\n{0} not found....\n') .format(args.polar)
            exit(1)
        if not os.path.exists(args.apolar):
            print ('\n{0} not found....\n') .format(args.polar)
            exit(1)
    
def ParseOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument("-mt", "--multiple", help='If given, calculate for multiple complexes. Need Metafile containing path of energy files', action="store_true")
    parser.add_argument("-mf", "--metafile", help='Metafile containing path to energy files of each complex in a row obtained from g_mmpbsa in following order: \
                                                       [MM file] [Polar file] [ Non-polar file] ',action="store", default='metafile.dat', metavar='metafile.dat')
    parser.add_argument("-m", "--molmech", help='Vacuum Molecular Mechanics energy file obtained from g_mmpbsa',action="store", default='energy_MM.xvg', metavar='energy_MM.xvg')
    parser.add_argument("-p", "--polar", help='Polar solvation energy file obtained from g_mmpbsa',action="store",default='polar.xvg', metavar='polar.xvg')
    parser.add_argument("-a", "--apolar", help='Non-Polar solvation energy file obtained from g_mmpbsa',action="store",default='apolar.xvg',metavar='polar.xvg')
    parser.add_argument("-bs", "--bootstrap", help='If given, Enable Boot Strap analysis',action="store_true")
    parser.add_argument("-nbs", "--nbstep", help='Number of boot strap steps for average energy calculation',action="store", type=int, default=500, metavar=500)
    parser.add_argument("-of", "--outfr", help='Energy File: All energy components frame wise',action="store",default='full_energy.dat', metavar='full_energy.dat')
    parser.add_argument("-os", "--outsum", help='Final Energy File: Full Summary of energy components',action="store",default='summary_energy.dat', metavar='summary_energy.dat')
    parser.add_argument("-om", "--outmeta", help='Final Energy File for Multiple Complexes: Complex wise final binding nergy',action="store",default='meta_energy.dat',metavar='meta_energy.dat')
    args = parser.parse_args()
    return args

def ReadData(FileName,n=2):
    infile = open(FileName,'r')
    x, data = [],[]
    for line in infile:
        line = line.rstrip('\n')
        if not line.strip():
            continue
        if(re.match('#|@',line)==None):
            temp = line.split()
            data.append(np.array(temp))
    for j in range(0,n):
        x_temp =[]
        for i in range(len(data)):
            x_temp.append(float(data[i][j]))
        x.append(x_temp)
    return x

def ComplexBootStrap(x,step=1000):
    avg =[]
    x = np.array(x)
    n = len(x)
    idx = np.random.randint(0,n,(step,n))
    sample_x = x[idx]
    avg = np.sort(np.mean(sample_x,1))
    CI_min = avg[int(0.005*step)]
    CI_max = avg[int(0.995*step)]
    #print('Energy = %13.3f; Confidance Interval = (-%-5.3f / +%-5.3f)\n' % (np.mean(avg), (np.mean(avg)-CI_min), (CI_max-np.mean(avg))))
    return avg, np.mean(avg), np.std(avg), [(np.mean(avg)-CI_min), (CI_max-np.mean(avg))]

def BootStrap (x,step=1000):
    if(np.mean(x)) == 0:
        return 0.000, 0.000
    else:
        avg =[]
        x = np.array(x)
        n = len(x)
        idx = np.random.randint(0,n,(step,n))
        sample_x = x[idx]
        avg = np.sort(np.mean(sample_x,1))
        return np.mean(avg),np.std(avg)

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def ReadMetafile(metafile):
    MmFile,PolFile, APolFile = [], [], []
    FileList = open(metafile,'r')
    for line in FileList:
        line = line.rstrip('\n')
        if not line.strip():
            continue
        temp = line.split()    
        MmFile.append(temp[0])
        PolFile.append(temp[1])
        APolFile.append(temp[2])
        
        
        if not os.path.exists(temp[0]):
            print ('\n{0} not found....\n') .format(temp[0])
            exit(1)
        if not os.path.exists(temp[1]):
            print ('\n{0} not found....\n') .format(temp[1])
            exit(1)
        if not os.path.exists(temp[2]):
            print ('\n{0} not found....\n') .format(temp[2])
            exit(1)
                        
    return MmFile, PolFile, APolFile

if __name__=="__main__":
    main()
EOF
chmod 755 stats.py
python3 stats.py -bs -nbs 2000 -m energy_MM.xvg -p polar.xvg -a apolar.xvg
sed '1d;2d;3d' full_energy.dat | sed '1i\,Time,E_VdW_mm(Protein),E_Elec_mm(Protein),E_Pol(Protein),E_Apol(Protein),E_VdW_mm(Ligand),E_Elec_mm(Ligand),E_Pol(Ligand),E_Apol(Ligand),E_VdW_mm(Complex),E_Elec_mm(Complex),E_Pol(Complex),E_Apol(Complex),Delta_E_mm,Delta_E_Pol,Delta_E_Apol,Delta_E_binding' |  sed "s/ \{1,\}/,/g" > full_energy_corrected.csv 
rm -rf full_energy.dat stats.py
sed -i 's/^,//' full_energy_corrected.csv
cat << EOF > csv_plot_fe.py
#!/usr/bin/python3
import os, sys, matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import io
import pandas as pd
import csv
df = pd.read_csv("full_energy_corrected.csv",index_col=False)
columns_name = list(df.columns)
df.plot(kind = 'line', x = columns_name[0], y = columns_name[13:17])
plt.legend()
x = round(df[columns_name[1]].max(),2)
#plt.text(3.5, x,text, fontsize = 10, color = 'b')
plt.grid(visible=True, which='major', color='grey', linestyle='--')
plt.ylabel('∆G binding free energy (kJ/mol)')
file_name = 'full_energy_corrected' + '.tiff' 
plt.xlabel(columns_name[0])
plt.suptitle('MM/PBSA-based ∆G Calculation')
plt.title('MM/PBSA-base ∆G')
plt.savefig(file_name, dpi=300)
EOF
chmod 755 csv_plot_fe.py
python3 csv_plot_fe.py
cat << EOF > decomp.py
#!/usr/bin/python3
import re
import numpy as np
import argparse
import sys
import os
import math

def main():
	args = ParseOptions()
	CheckInput(args)
	MMEnData,resnameA = ReadData(args.molmech)
	polEnData,resnameB = ReadData(args.polar)
	apolEnData,resnameC = ReadData(args.apolar)
	resname = CheckResname(resnameA,resnameB,resnameC)
	print("Total number of Residue: {0}\n" . format(len('resname')+1))
	Residues = []
	for i in range(len(resname)):
		CheckEnData(MMEnData[i],polEnData[i],apolEnData[i])
		r = Residue()
		r.CalcEnergy(MMEnData[i],polEnData[i],apolEnData[i],args)
		Residues.append(r)
		print(' %8s %8.4f %8.4f' % (resname[i], r.TotalEn[0], r.TotalEn[1]))
	fout = open(args.output,'w')
	fmap = open(args.outmap,'w')
	fout.write('#Residues  MM Energy(+/-)dev/error  Polar Energy(+/-)dev/error APolar Energy(+/-)dev/error Total Energy(+/-)dev/error\n')
	for i in range(len(resname)):
		if (args.cutoff == 999):
			fout.write("%-8s  %4.4f  %4.4f    %4.4f  %4.4f    %4.4f  %4.4f    %4.4f  %4.4f \n" %(resname[i],Residues[i].FinalMM[0],Residues[i].FinalMM[1], 
                                      Residues[i].FinalPol[0], Residues[i].FinalPol[1], Residues[i].FinalAPol[0], Residues[i].FinalAPol[1], Residues[i].TotalEn[0],Residues[i].TotalEn[1] ))
		elif (args.cutoff <= Residues[i].TotalEn[0]) or ( (-1 *args.cutoff) >= Residues[i].TotalEn[0]):
			fout.write("%-8s  %4.4f  %4.4f    %4.4f  %4.4f    %4.4f  %4.4f    %4.4f  %4.4f \n" % (resname[i],Residues[i].FinalMM[0],Residues[i].FinalMM[1], 
                                  Residues[i].FinalPol[0], Residues[i].FinalPol[1], Residues[i].FinalAPol[0], Residues[i].FinalAPol[1], Residues[i].TotalEn[0],Residues[i].TotalEn[1] ))
			
		fmap.write("%-8d     %4.4f \n" %((i+1), Residues[i].TotalEn[0]))


class Residue():
	def __init__(self):
		self.FinalMM, self.FinalPol, self.FinalAPol, self.TotalEn = [], [], [], []
		
	def BootStrap(self,x,step):
		avg =[]
		x = np.array(x)
		n = len(x)
		idx = np.random.randint(0,n,(step,n))
		sample_x = x[idx]
		avg = np.sort(np.mean(sample_x,1))
		return np.mean(avg), np.std(avg)

	def CalcEnergy(self,MM,Pol,APol,args):
		TotalEn = np.sum([MM,Pol,APol],axis=0)
		if(args.bootstrap):
			self.FinalMM = self.BootStrap(MM,args.nbstep)
			self.FinalPol = self.BootStrap(Pol,args.nbstep)
			if (np.mean(APol) == 0):
				self.FinalAPol = [0.0,0.0]
			else:
				self.FinalAPol = self.BootStrap(APol,args.nbstep)
			self.TotalEn = self.BootStrap(TotalEn,args.nbstep)
		else:
			self.FinalMM = [np.mean(MM),np.std(MM)]
			self.FinalPol = [np.mean(Pol),np.std(Pol)]
			if (np.mean(APol) == 0):
				self.FinalAPol = [0.0,0.0]
			else:
				self.FinalAPol = [np.mean(APol),np.std(APol)]
			self.TotalEn = [np.mean(TotalEn),np.std(TotalEn)]
		self.FinalMM = np.round(self.FinalMM,4)
		self.FinalPol = np.round(self.FinalPol,4)
		self.FinalAPol = np.round(self.FinalAPol,4)
		self.TotalEn = np.round(self.TotalEn,4)

def CheckEnData(MM,Pol,APol):
	if(len(Pol) != len(MM)):
		print("Times or Frames Mismatch between files")
		exit(1)
	if(len(APol) != len(Pol)):
		print("Times or Frames Mismatch between files")
		exit(1)
	if(len(APol) != len(MM)):
		print("Times or Frames Mismatch between files")
		exit(1)


def ParseOptions():
        parser = argparse.ArgumentParser()
        parser.add_argument("-m", "--molmech", help='Molecular Mechanics energy file',action="store", default='contrib_MM.dat', metavar='contrib_MM.dat')
        parser.add_argument("-p", "--polar", help='Polar solvation energy file',action="store",default='contrib_pol.dat', metavar='contrib_pol.dat')
        parser.add_argument("-a", "--apolar", help='Non-Polar solvation energy file',action="store",default='contrib_apol.dat',metavar='contrib_apol.dat')
        parser.add_argument("-bs", "--bootstrap", help='Switch for Error by Boot Strap analysis',action="store_true")
        parser.add_argument("-nbs", "--nbstep", help='Number of boot strap steps',action="store", type=int,default=500, metavar=500)
        parser.add_argument("-ct", "--cutoff", help='Absolute Cutoff: energy output above and below this value',action="store",type=float,default=999, metavar=999)
        parser.add_argument("-o", "--output", help='Final Decomposed Energy File',action="store",default='final_contrib_energy.dat', metavar='final_contrib_energy.dat')
        parser.add_argument("-om", "--outmap", help='energy2bfac input file: to map energy on structure for visualization',action="store",default='energyMapIn.dat', metavar='energyMapIn.dat')
        return parser.parse_args()

def CheckResname(resA,resB,resC):
	if(len(resA) != len(resB)):
		print("ERROR: Total number of residues mismatch between files")
		exit(1)
	if(len(resB) != len(resC)):
		print("ERROR: Total number of residues mismatch between files")
		exit(1)
	if(len(resC) != len(resA)):
		print("ERROR: Total number of residues mismatch between files")
		exit(1)
	for i in range(len(resA)):
		if (resA[i] != resB[i]):
			print("ERROR: Residue mismatch between files")
			exit(1)
	for i in range(len(resB)):
		if (resB[i] != resC[i]):
			print("ERROR: Residue mismatch between files")
			exit(1)
	for i in range(len(resA)):
		if (resA[i] != resC[i]):
			print("ERROR: Residue mismatch between files")
			exit(1)
	return resA

def CheckInput(args):
	if not os.path.exists(args.molmech):
		print('\n{0} not found....\n' .format(args.molmech))
		exit(1)
	if not os.path.exists(args.polar):
		print('\n{0} not found....\n' .format(args.polar))
		exit(1)
	if not os.path.exists(args.apolar):
		print('\n{0} not found....\n' .format(args.polar))
		exit(1)


def ReadData(FileName):
        infile = open(FileName,'r')
        x, data,resname = [],[],[]
        for line in infile:
                line = line.rstrip('\n')
                if not line.strip():
                        continue
                if(re.match('#|@',line)==None):
                        temp = line.split()
                        data.append(np.array(temp))
                if(re.match('#',line)):
                        resname = line.split()
                        n = len(resname[1:])
        for j in range(1,n):
                x_temp =[]
                for i in range(len(data)):
                        x_temp.append(float(data[i][j]))
                x.append(x_temp)
        return x, resname[2:]

if __name__=="__main__":
        main()
EOF
chmod 755 decomp.py
python3 decomp.py -bs -nbs 2000 -m contrib_MM.dat -p contrib_pol.dat -a contrib_apol.dat
sed '1d' final_contrib_energy.dat | sed '1i\Residues,MM_Energy,dev-error,Polar_Energy,dev-error,APolar_Energy,dev-error,Total_Energy,dev-error' | sed "s/ \{1,\}/,/g" > final_contrib_energy_corrected.csv
rm -rf final_contrib_energy.dat decomp.py
sed -i 's/^,//' final_contrib_energy_corrected.csv
rm -rf final_contrib_energy.dat decomp.py
cat << EOF > csv_plot.py
#!/usr/bin/python3
import os, sys, matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import io
import pandas as pd
import csv
df = pd.read_csv("final_contrib_energy_corrected.csv",index_col=False)
columns_name = list(df.columns)
df.plot(kind = 'line', x = columns_name[0], y = columns_name[7])
plt.legend()
#x = round(df[columns_name[1]].max(),2)
#plt.text(3.5, x,text, fontsize = 10, color = 'b')
plt.grid(visible=True, which='major', color='grey', linestyle='--')
plt.ylabel('∆G binding free energy (kJ/mol)')
file_name = 'final_contrib_energy' + '.tiff' 
plt.xlabel(columns_name[0])
plt.suptitle('MM/PBSA-based ∆G Calculation (P-P Complex)')
plt.title('Energy Decomposition Analysis')
plt.savefig(file_name, dpi=300)
EOF
chmod 755 csv_plot.py
python3 csv_plot.py
}



Help

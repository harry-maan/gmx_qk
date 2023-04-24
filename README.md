# gmx_qk
gmx_qk is an automated bash workflow that enables users to execute simulations of **protein in water/protein-ligand complexes** with only a basic understanding of **UNIX or command-line tools**.
After receiving the input files and parameters, it starts the MD simulations in a matter of seconds, whereas a command-line-based protocol would take **20 to 30 minutes**.
Through the use of **g_mmpbsa**, it acts as a bridge between the **MD engine** and **MM/PBSA-based binding free energy calculations** (gromacs version 2021.4)
Simplicity and lack of expertise in complicated tool installation or compilation processes, such as Gromacs command line-based utilities.
Users are assisted in producing reproducible research outcomes by a single workflow.
  
It is now available in public doamain with GPL 3.0 Licence. (It can be used free of cost, modified and distributed with no warranty).
Please cite the publication 
Gmx_qk: An Automated Protein/Protein–Ligand Complex Simulation Workflow Bridged to MM/PBSA, Based on Gromacs and Zenity-Dependent GUI for Beginners in MD Simulation Study

Harvinder Singh, Anupam Raja, Ajay Prakash, and Bikash Medhi
Journal of Chemical Information and Modeling Article ASAP
DOI: 10.1021/acs.jcim.3c00341

**NOTE: Make sure that internet connection is working during installation to download its all dependencies from different repositories i.e., Python modules, and build-essentials**
## Download
	git clone https://github.com/harry-maan/gmx_qk.git
## Installation
	cd gmx_qk/gmx_qk
	ls
	sudo bash configure.sh
## usage 
	bash gmx_qk
or or look application section

Thank you for taking the time to read this post. If you found the article useful, please consider citing it in your own work. Your support would be greatly appreciated and would help recognize the scientific contribution made in this research.
[Support with Paytm](9587677525@paytm)

Harvinder Singh,
PhD scholar,
Dept. of Pharmacology
PGIMER, Chandigarh
harvindermaan4@gmail.com

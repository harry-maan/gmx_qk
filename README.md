# gmx_qk
gmx_qk is an automated bash workflow that enables users to execute simulations of **protein in water/protein-ligand complexes** with only a basic understanding of **UNIX or command-line tools**.
After receiving the input files and parameters, it starts the MD simulations in a matter of seconds, whereas a command-line-based protocol would take **20 to 30 minutes**.
Through the use of **g_mmpbsa**, it acts as a bridge between the **MD engine** and **MM/PBSA-based binding free energy calculations** (gromacs version 2021.4)
Simplicity and lack of expertise in complicated tool installation or compilation processes, such as Gromacs command line-based utilities.
Users are assisted in producing reproducible research outcomes by a single workflow.
  
As its code under development phase, no publication has been made, so it is not available in public doamain. (comming soon after any publication record).


**NOTE: Make sure that internet connection is working during installation to download its all dependencies from different repositories i.e., Python modules, and build-essentials**
## Download
	wget --no-check-certificate 'https://drive.google.com/file/d/1iqEs7Qb-1YS6XEEc7qjj9XF5cevd1Z42/view?usp=share_link' -o gmx_qk
## Installation
	cd gmx_qk
	ls
	sudo bash configure.sh
## usage
	bash gmx_qk or look application section 

Thanks for being here.

Harvinder Singh,
		
PhD scholar,

Dept. of Pharmacology

PGIMER, Chandigarh

harvindermaan4@gmail.com

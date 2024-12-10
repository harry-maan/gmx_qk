#!/bin/bash
Help()
{
  ans=$(zenity --text-info --title 'gmx_trj_qk 1.0.0' \
      --html \
      --ok-label="Quit" \
      --extra-button="Tutorial" \
      --extra-button="Workflow" \
      --extra-button "Trajectory Analysis" \
      --extra-button "FEL (Free Energy Landscape)" \
      --extra-button "Back to gmx_qk" \
      --width=600 \
      --height=600 \
      --filename=<(echo "
<html>
    <body>
        <div style='text-align:center; font-family: Arial, sans-serif;'>
            <h1 style='color: blue;'>Welcome to gmx_trj_qk</h1>
            <img src="/usr/share/icons/gmx_qk" alt='gmx_qk Logo' style='display: block; margin: 20px auto; width: 200px; height: auto;'>
            <p style='font-size: 18px; text-align: justify;'>
                The gmx_trj_qk is a <b>Zenity, Python, and Gromacs</b> dependent bash program. 
                It's designed for beginners to GROMACS who wish to analyze the <b>protein or protein-ligand complexes Molecular Dynamics simulation trajectory, including RMSD, RMSF, Radius of Gyration, H-bond and contacts calculations</b>. 
                <br><br>
                gmx_trj_qk is a fully automated program, efficiently compatible with GROMACS 5.0 and newer versions such as 2021.4. 
                Informative widgets are supported by <b>Zenity (GUI)</b>. 
                Additional functionality for <b>post-MD simulation trajectory Free Energy Landscape (FEL)</b> analysis has also been embedded.
                
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
  elif [[ $ans = "Back to gmx_qk" ]]
  then
    gmx_qk;
  elif [[ $ans = "Trajectory Analysis" ]]
  then
    simulation_type;
  elif [[ $ans = "FEL (Free Energy Landscape)" ]]
  then
    gmx_fel;
  fi
   
}
function simulation_type(){
  dir=$(zenity  --file-selection --title="Choose a working directory" --directory)
  cd $dir
  go_no_go=$(zenity --forms --title "P-L complex Simulation/P-P complex Simulation/Protein in water" --text "Simulation type" \
  --add-combo "Select a simulation type" --combo-values "Protein-Ligand (P-L complex) simulation analysis|Protein-Protein (P-P complex) simulation analysis|Protein in water simulation analysis|Done")
  
  case "${go_no_go}" in
  "Protein-Ligand (P-L complex) simulation analysis")
  want_to_do;;
  "Protein-Protein (P-P complex) simulation analysis")
  want_to_do_pp1;;
  "Protein in water simulation analysis")
  want_to_do_p;;
  "Done") 
  break;;
esac
}


function want_to_do(){
      data=$(zenity --list --checklist --column "TRUE", --column "Analysis"\
      --title="Checklis for Analysis want to do?" \
      --text="1. Available analysis modules are listed please tick the box.\n2. It will be resulted publishable plots and concerned raw data (.csv)" \
      TRUE RMSD \
      TRUE RMSF \
      TRUE 'Radius_of_Gyration' \
      TRUE 'Contacts_and_Distanse' \
      TRUE 'H_bond_Calculations' --width 300 --height 400) ;

choose1;
}
function want_to_do_pp1(){
      data=$(zenity --list --checklist --column "TRUE", --column "Analysis"\
      --title="Checklis for Analysis want to do?" \
      --text="1. Available analysis modules are listed please tick the box.\n2. It will be resulted publishable plots and concerned raw data (.csv)" \
      TRUE RMSD \
      TRUE RMSF \
      TRUE 'Radius_of_Gyration' \
      TRUE 'Contacts_and_Distanse' \
      TRUE 'H_bond_Calculations' --width 300 --height 400) ;

choose1;
}
function want_to_do_p(){
      data=$(zenity --list --checklist --column "TRUE", --column "Analysis"\
      --title="Checklis for Analysis want to do?" \
      --text="1. Available analysis modules are listed please tick the box.\n2. It will be resulted publishable plots and concerned raw data (.csv)" \
      TRUE RMSD \
      TRUE RMSF \
      TRUE 'Radius_of_Gyration' --width 300 --height 400) ;

choose1;
}




function choose1() {
  tpr="$(zenity --file-selection --title='Select an input md.tpr(.tpr)')"
  case $? in
           0)
                  choose2;;
           1)
                  zenity --question \
                         --title="Select an input md.tpr(.tpr)" \
                         --text="No md.tpr(.tpr) selected. Do you want to select one?" \
                         && choose1 || exit;;
          -1)
                  echo "An unexpected error has occurred."; exit;;
  esac
}
function choose2() {
  trajectory="$(zenity --file-selection --title='Please select a trajectory file(.xtc)')"
  case $? in
           0)
                  choose3;;
           1)
                  zenity --question \
                         --title="Select a trajectory file (.xtc)" \
                         --text="No trajectory file (.xtc). Do you want to select one?" \
                         && choose2 || exit;;
          -1)
                  echo "An unexpected error has occurred."; exit;;
  esac
}
function choose3() {
  index="$(zenity --file-selection --title='Select an index file(.ndx)')"
  case $? in
           0)
                  analysis;RMSD;RMSF;Radius_of_Gyration;Contacts_and_Distanse;H_bond_Calculations;plotting;;
           1)
                  zenity --question \
                         --title="Select an index file(.ndx)" \
                         --text="No index ile(.ndx). Do you want to select one?" \
                         && choose3 || exit;;
          -1)
                  echo "An unexpected error has occurred."; exit;;
  esac
}




function analysis() {
	cp $tpr $index .
	source /usr/local/gromacs/bin/GMXRC
	printf "Protein\nSystem\n" | gmx trjconv -s $tpr -f $trajectory -center -ur compact -pbc mol -o md_center.xtc
	printf "Backbone\nSystem\n" | gmx trjconv -s $tpr -f md_center.xtc -fit rot+trans -o md_fit.xtc
	printf "keep 1\na CA\nname 1 Calpha\nq\n" |gmx make_ndx -f $tpr -o CA.ndx
	
}

# 	 	         	 	 		gmx covar -f md_fit.xtc -s $tpr -v eigenvec.trr #gmx anaeig -v eigenvec.trr -f md_fit.xtc -eig eigenval.xvg -s $tpr -first 1 -last 2 -2d
# 	 	         	 	 		fi


function RMSD(){
 case "${data}" in
          *"RMSD"*)
            printf "Calpha\nCalpha\n" | gmx rms -s $tpr -f md_fit.xtc -n CA.ndx -o rmsd.xvg -fit rot+trans;file1='rmsd.xvg';;
     esac
}
function RMSF(){
 case "${data}" in
            *"RMSF"*)
          printf "Calpha\n" | gmx rmsf -s $tpr -f md_fit.xtc -n CA.ndx -o rmsf.xvg -fit;file2='rmsf.xvg';;
     esac
}
function Radius_of_Gyration(){
 case "${data}" in
            *"Radius_of_Gyration"*)
          echo Protein | gmx gyrate -s $tpr -f md_fit.xtc -o gyrate.xvg;file3='gyrate.xvg';;
     esac
}

function Contacts_and_Distanse(){
 case "${data}" in
          *"Contacts_and_Distanse"*)
          if [[ $go_no_go == "Protein-Ligand (P-L complex) simulation analysis" ]]
          then
           printf "1\n13\n" |gmx mindist -f md_fit.xtc -s $tpr -n $index -od distance.xvg -on count.xvg -or residue.xvg
          elif [[ $go_no_go == "Protein-Protein (P-P complex) simulation analysis" ]]
          then
           printf "Protein_chain1\nProtein_chain2\n" |gmx mindist -f md_fit.xtc -s $tpr -n $index -od distance.xvg -on count.xvg -or residue.xvg
           fi;file4='distance.xvg';file5='count.xvg';file6='residue.xvg';;
     esac
}

function H_bond_Calculations(){
 	        case "${data}" in
          *"H_bond_Calculations"*)
          if [[ $go_no_go == "Protein-Ligand (P-L complex) simulation analysis" ]]
          then
          printf "1\n13\n" | gmx hbond -f md_fit.xtc -s $tpr -n $index -num hbnum.xvg
          elif [[ $go_no_go == "Protein-Protein (P-P complex) simulation analysis" ]]
          then
           printf "Protein_chain1\nProtein_chain2\n" | gmx hbond -f md_fit.xtc -s $tpr -n $index -num hbnum.xvg
          fi;file7='hbnum.xvg';;
     esac
}




function plotting(){
for f in $file1 $file2 $file3 $file4 $file5 $file6 $file7  #*.xvg #count.xvg residue.xvg distance.xvg rmsd.xvg rmsf.xvg gyrate.xvg 
do
	myvar="$(grep -w 'yaxis' "$f")";yaxis=$(echo $myvar |grep -o -P '(?<=").*(?=")');echo $yaxis # extraction of ylabel
	myvar="$(grep -w 'xaxis' "$f")";xaxis=$(echo $myvar |grep -o -P '(?<=").*(?=")');echo $xaxis # xlabel extraction
	myvar="$(grep -w 'title' "$f")";suptitle=$(echo $myvar |grep -o -P '(?<=").*(?=")');echo $suptitle # title extraction
	myvar="$(grep -w 's0 legend' "$f")";title=$(echo $myvar |grep -o -P '(?<=").*(?=")');echo $title
	myvar="$(grep -w 's1 legend' "$f")";s1=$(echo $myvar |grep -o -P '(?<=").*(?=")');echo $s1
	myvar="$(grep -w 's2 legend' "$f")";s2=$(echo $myvar |grep -o -P '(?<=").*(?=")');echo $s2
	myvar="$(grep -w 's3 legend' "$f")";s3=$(echo $myvar |grep -o -P '(?<=").*(?=")');echo $s3
	sed -i '/#/d' "$f";sed -i '/@/d' "$f"  # remove header of the xvg file
	sed "s/ \{1,\}/,/g" "$f" > "${f%%.*}"2.csv #converted into csv file
	sed -i 's/^,//' "${f%%.*}"2.csv #;rm -rf "${f%%.*}"2.csv "${f%%.*}".xvg
	echo ",$xaxis,$yaxis,$s1,$s2,$s3" > "${f%%.*}".csv && cat "${f%%.*}"2.csv >> "${f%%.*}".csv #insert headers in csv file
	sed -i 's/^,//' "${f%%.*}".csv ;rm -rf "${f%%.*}"2.csv "${f%%.*}".xvg
	cat << EOF > csv_plot.py
#!/usr/bin/python3
import os, sys, matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import io
import pandas as pd
import csv
df = pd.read_csv("${f%%.*}.csv",index_col=False)
columns_name = list(df.columns)
mean = df[columns_name[1]].mean()
sd = df[columns_name[1]].std()
sd_string=str(round(sd, 2))
mean_string=str(round(mean, 2))
text = '$yaxis' + ' "Mean +/- Sd (' + mean_string + ' +/- ' + sd_string + ')"'
df.plot(kind = 'line', x = columns_name[0], y = columns_name[1], label = text)
plt.legend(fontsize=7)
x = round(df[columns_name[1]].max(),2)
#plt.text(3.5, x,text, fontsize = 10, color = 'b')
plt.grid(visible=True, which='major', color='grey', linestyle='--')
plt.ylabel('$yaxis')
file_name = '${f%%.*}' + '.tiff' 
plt.xlabel('$xaxis')
plt.suptitle('$suptitle')
plt.title('$title')
plt.savefig(file_name, dpi=300)
EOF


	chmod 755 csv_plot.py
	python3 csv_plot.py
done
rm -rf csv_plot.py
}
Help

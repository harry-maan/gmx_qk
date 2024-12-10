#!/bin/bash
Help()
{
  ans=$(zenity --text-info --title 'gmx_fel 1.0.0' \
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
            <h1 style='color: blue;'>Welcome to gmx_fel</h1>
            <img src="/usr/share/icons/gmx_qk" alt='gmx_qk Logo' style='display: block; margin: 20px auto; width: 200px; height: auto;'>
            <p style='font-size: 18px; text-align: justify;'>
                The gmx_fel is a <b>Zenity, Python, Gromacs, and g_mmpbsa</b> dependent bash program. 
                It's designed for beginners to GROMACS who wish conduct <b> Free Energy Landscape analysis for protein in water or protein-ligand complexes Molecular Dynamics Simulation trajectory including Gibbs free energy, Entropy, Enthalpy and Probability Distribution calculations</b>. 
                <br><br>
                gmx_fel is a fully automated program, efficiently compatible with GROMACS 5.0 and newer versions such as 2021.4. 
                Informative widgets are supported by <b>Zenity (GUI)</b>. 
                
                
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
    zenity --text-info --title="gmx_fel 1.0.0" --html --url=$theurl \
       --checkbox="I read it...and I'm good to go" --width 800 --height 600;
       Help;
  elif [[ $ans = "Workflow" ]]
  then
    theurl="https://github.com/harry-maan/gmx_qk/blob/main/workflow.png"
    zenity --text-info --title="gmx_fel 1.0.0" --html --url=$theurl \
       --checkbox="I read it...and I'm good to go" --width 800 --height 600;
       Help;
  elif [[ $ans = "License" ]]
  then
    theurl="https://github.com/harry-maan/gmx_qk/blob/main/license.png"
    zenity --text-info --title="gmx_fel 1.0.0" --html --url=$theurl \
       --checkbox="I read it...and I'm good to go" --width 800 --height 600;
       Help;
  elif [[ $ans = "Let's go!!" ]]
  then
    simulation_type;
  fi
   
}
function simulation_type(){
  dir=$(zenity  --file-selection --title="Choose a working directory" --directory)
  cd $dir
  go_no_go=$(zenity --forms --title "P-L complex Simulation/Protein in water" --text "Simulation type" \
  --add-combo "Select a simulation type" --combo-values "Protein-Ligand (P-L complex) simulation anal
ysis|Protein in water simulation analysis|Done")
  
  case "${go_no_go}" in
  "Protein-Ligand (P-L complex) simulation analysis")
  want_to_do;;
  "Protein in water simulation analysis")
  want_to_do;;
  "Done") 
  break;;
esac
}
function want_to_do(){
      data=$(zenity --list --checklist --column "TRUE", --column "Analysis"\
      --title="Checklis for Analysis want to do?" \
      --text="1. Available analysis modules are listed please tick the box.\n2. It will be resulted p
ublishable plots and concerned raw data (.csv)" \
      TRUE Gibbs \
      TRUE Entropy \
      TRUE 'PDF'  --width 300 --height 400) ;

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
                  step_data;analysis;Gibbs;; #Entropy;PDF;plotting;;
           1)
                  zenity --question \
                         --title="Select a trajectory file (.xtc)" \
                         --text="No trajectory file (.xtc). Do you want to select one?" \
                         && choose2 || exit;;
          -1)
                  echo "An unexpected error has occurred."; exit;;
  esac
}

function step_data(){
data_step=$(zenity --forms --separator="," \
--title="Time Intervals you would like to slice trajectory (ns)" \
--text="1. Please enter digits only to slice the trajectory.\n 2. Total duration of simulation in nanoseconds. \n 3. Tiles format of images i.e., 5X2 for 10 slices of trajectory. \n 4. Figure Quality in dpi (numerics only)" \
--add-entry="Factor for slice trajectory" \
--add-entry="Total Duration of Trajectory (ns)"  \
--add-entry="Tiles of images" \
--add-entry="DPI");

step=$( echo $data_step | awk -F ',' '{print $1}' )
total_t=$( echo $data_step | awk -F ',' '{print $2}' )
tiles=$( echo $data_step | awk -F ',' '{print $3}' )
dpi=$( echo $data_step | awk -F ',' '{print $4}' )
echo ${step}
echo ${total_t}
echo ${tiles}
echo ${dpi} 
}


function analysis() {
        cp $tpr .
        source /usr/local/gromacs/bin/GMXRC
        gmx make_ndx -f $tpr -o index.ndx <<EOF
4 | 13
q
EOF
        printf "Protein\nSystem\n" | gmx trjconv -s $tpr -f $trajectory -center -ur compact -pbc mol -o md_center.xtc
        printf "Backbone\nSystem\n" | gmx trjconv -s $tpr -f md_center.xtc -fit rot+trans -o md_fit.xtc
        printf "keep 1\na CA\nname 1 Calpha\nq\n" |gmx make_ndx -f $tpr -o CA.ndx
        
}



function Gibbs(){
 case "${data}" in
          *"Gibbs"*)
          if [[ $go_no_go == "Protein-Ligand (P-L complex) simulation analysis" ]]
          then
           G_pl;
          elif [[ $go_no_go == "Protein in water simulation analysis" ]]
          then
           G_pw;
           fi;;
     esac
}
function Entropy(){
 case "${data}" in
            *"Entropy"*)
          printf "Calpha\n" | gmx Entropy -s $tpr -f md_fit.xtc -n CA.ndx -o Entropy.xvg -fit;file2='Entropy.xvg';;
     esac
}
function PDF(){
 case "${data}" in
            *"PDF"*)
          echo Protein | gmx gyrate -s $tpr -f md_fit.xtc -o gyrate.xvg;file3='gyrate.xvg';;
     esac
}

function Contacts_and_Distanse(){
 case "${data}" in
          *"Contacts_and_Distanse"*)
          if [[ $go_no_go == "Protein-Ligand (P-L complex) simulation analysis" ]]
          then
           G_pl;
          elif [[ $go_no_go == "Protein in water simulation analysis" ]]
          then
           G_pw;
           fi;;
     esac
}

function G_pl(){
STEP=$step

# Loop 
for ((i=0; i<${total_t}; i+=STEP)); do
  # Calculate the start and end of the range
  start=$i
  end=$((i + STEP))

  # Create the filenames based on the range
  eigenval="md_${start}_${end}_eigenval.xvg"
  eigenvec="md_${start}_${end}_eigenvec.trr"
  pca_2dproj="PCA_2dproj_md_${start}_${end}.xvg"
  fel_sham="FEL_PCA_sham_md_${start}_${end}"

  # covar
  echo 21 21 | gmx covar -f md_fit.xtc -s $tpr -n index.ndx -o $eigenval -v $eigenvec -b $start -e $end -tu ns

  # anaeig
  echo 21 21 |gmx anaeig -v $eigenvec -f md_fit.xtc -eig $eigenval -s $tpr -first 1 -last 2 -2d $pca_2dproj -n index.ndx

  # sham
  gmx sham -f $pca_2dproj -ls ${fel_sham}.xpm -notime
  
  #change the title
  sed -i "s|/\* title:.*|/* title:   \"${start}_${end}ns\" */|" ${fel_sham}.xpm
  sed -i "s|/\* title:.*|/* title:   \"${start}_${end}ns\" */|" ${fel_sham}.xpm
  sed -i "s|/\* title:.*|/* title:   \"${start}_${end}ns\" */|" ${fel_sham}.xpm
  sed -i "s|/\* title:.*|/* title:   \"${start}_${end}ns\" */|" entropy.xpm
  sed -i "s|/\* title:.*|/* title:   \"${start}_${end}ns\" */|" enthalpy.xpm
  sed -i "s|/\* title:.*|/* title:   \"${start}_${end}ns\" */|" prob.xpm

  # Run xpm2ps for Gibbs free energy (G), Entropy (TDS), Enthalpy (H), Prob 
  gmx xpm2ps -f ${fel_sham}.xpm -o ${fel_sham}.eps -rainbow red
  gmx xpm2ps -f entropy.xpm -o entropy_${fel_sham}.eps -rainbow red
  gmx xpm2ps -f enthalpy.xpm -o enthalpy_${fel_sham}.eps -rainbow red
  gmx xpm2ps -f prob.xpm -o prob_${fel_sham}.eps -rainbow red


  rm *#
done
montage -tile 5X2 -geometry +2+2 -density 300 FEL_PCA_sham_md_*.eps combinde_landscape_G.eps
convert -density 300 combinde_landscape_G.eps combinde_landscape_G.tif
rm FEL_PCA_sham_md_*.eps
rm *#

# combine the plots using  montage and converting to .tif file (Entropy)
montage -tile 5X2 -geometry +2+2 -density 300 entropy_FEL_PCA_sham_md_*.eps combinde_landscape_entropy.eps
convert -density 300 combinde_landscape_entropy.eps combinde_landscape_entropy.tif
rm entropy_FEL_PCA_sham_md_*.eps


# combine the plots using  montage and converting to .tif file (Enthalpy)
montage -tile 5X2 -geometry +2+2 -density 300 enthalpy_FEL_PCA_sham_md_*.eps combinde_landscape_enthalpy.eps
convert -density 300 combinde_landscape_enthalpy.eps combinde_landscape_enthalpy.tif
rm enthalpy_FEL_PCA_sham_md_*.eps


## combine the plots using  montage and converting to .tif file (prob)
montage -tile 5X2 -geometry +2+2 -density 300 prob_FEL_PCA_sham_md_*.eps combinde_landscape_prob.eps
convert -density 300 combinde_landscape_prob.eps combinde_landscape_prob.tif
rm prob_FEL_PCA_sham_md_*.eps
rm -rf *.xpm *.eps
}

function G_pw(){
STEP=${step}

# Loop 
for ((i=0; i<${total_t}; i+=STEP)); do
  # Calculate the start and end of the range
  start=$i
  end=$((i + STEP))

  # Create the filenames based on the range
  eigenval="md_${start}_${end}_eigenval.xvg"
  eigenvec="md_${start}_${end}_eigenvec.trr"
  pca_2dproj="PCA_2dproj_md_${start}_${end}.xvg"
  fel_sham="FEL_PCA_sham_md_${start}_${end}"

  # covar
  echo 4 4 | gmx covar -f md_fit.xtc -s $tpr -n index.ndx -o $eigenval -v $eigenvec -b $start -e $end -tu ns

  # anaeig
  echo 4 4 |gmx anaeig -v $eigenvec -f md_fit.xtc -eig $eigenval -s $tpr -first 1 -last 2 -2d $pca_2dproj -n index.ndx

  # sham
  gmx sham -f $pca_2dproj -ls ${fel_sham}.xpm -notime
  
  #change the title
  sed -i "s|/\* title:.*|/* title:   \"${start}_${end}ns\" */|" ${fel_sham}.xpm
  sed -i "s|/\* title:.*|/* title:   \"${start}_${end}ns\" */|" ${fel_sham}.xpm
  sed -i "s|/\* title:.*|/* title:   \"${start}_${end}ns\" */|" entropy.xpm
  sed -i "s|/\* title:.*|/* title:   \"${start}_${end}ns\" */|" enthalpy.xpm
  sed -i "s|/\* title:.*|/* title:   \"${start}_${end}ns\" */|" prob.xpm

  # Run xpm2ps for Gibbs free energy (G), Entropy (TDS), Enthalpy (H), Prob 
  gmx xpm2ps -f ${fel_sham}.xpm -o ${fel_sham}.eps -rainbow red
  gmx xpm2ps -f entropy.xpm -o entropy_${fel_sham}.eps -rainbow red
  gmx xpm2ps -f enthalpy.xpm -o enthalpy_${fel_sham}.eps -rainbow red
  gmx xpm2ps -f prob.xpm -o prob_${fel_sham}.eps -rainbow red


  rm *#
done
montage -tile ${tiles} -geometry +2+2 -density ${dpi} FEL_PCA_sham_md_*.eps combinde_landscape_G.eps
convert -density 300 combinde_landscape_G.eps combinde_landscape_G.tif
rm FEL_PCA_sham_md_*.eps
rm *#

# combine the plots using  montage and converting to .tif file (Entropy)
montage -tile ${tiles} -geometry +2+2 -density ${dpi} entropy_FEL_PCA_sham_md_*.eps combinde_landscape_entropy.eps
convert -density 300 combinde_landscape_entropy.eps combinde_landscape_entropy.tif
rm entropy_FEL_PCA_sham_md_*.eps


# combine the plots using  montage and converting to .tif file (Enthalpy)
montage -tile ${tiles} -geometry +2+2 -density ${dpi} enthalpy_FEL_PCA_sham_md_*.eps combinde_landscape_enthalpy.eps
convert -density 300 combinde_landscape_enthalpy.eps combinde_landscape_enthalpy.tif
rm enthalpy_FEL_PCA_sham_md_*.eps


## combine the plots using  montage and converting to .tif file (prob)
montage -tile ${tiles} -geometry +2+2 -density ${dpi} prob_FEL_PCA_sham_md_*.eps combinde_landscape_prob.eps
convert -density 300 combinde_landscape_prob.eps combinde_landscape_prob.tif
rm prob_FEL_PCA_sham_md_*.eps
rm -rf *.xpm *.eps
}

# combine the plots using  montage and converting to .tif file (Gibbs free energy)


Help




#!/bin/bash
function buil_essential(){
  sudo apt install -y build-essential
  sudo apt-get install -y zenity
  sudo apt install python2
  sudo apt-get install python3.8
  sudo apt-get install montage
  #sudo apt-get -y install python3-pip
  alias python=python3.8
  curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
  python3.8 get-pip.py
  pip3 install -r requirements.txt
  Help;
}

function install_gmx_qk(){
sudo cp gmx_qk.sh /usr/local/bin/gmx_qk
sudo cp gmx_trj_qk.sh /usr/local/bin/gmx_trj_qk
sudo cp gmx_fel.sh /usr/local/bin/gmx_fel
sudo chmod +x /usr/local/bin/gmx_qk
sudo chmod +x /usr/local/bin/gmx_fel
sudo chmod +x /usr/local/bin/gmx_trj_qk
sudo cp logo.png /usr/share/icons/gmx_qk
sudo chmod +x /usr/share/icons/gmx_qk
cat << EOF > gmx_qk.desktop
[Desktop Entry]
Version=1.0
Exec=/usr/local/bin/gmx_qk
Name=gmx_qk
Icon=/usr/share/icons/gmx_qk
Terminal=true
Type=Application
Categories=Application
EOF
cat << EOF > gmx_trj_qk.desktop
[Desktop Entry]
Version=1.0
Exec=/usr/local/bin/gmx_trj_qk
Name=gmx_trj_qk
Icon=/usr/share/icons/gmx_trj_qk
Terminal=true
Type=Application
Categories=Application
EOF
sudo cp gmx_qk.desktop /usr/share/applications
sudo chmod +x /usr/share/applications/gmx_qk.desktop
tar -zxvf g_mmpbsa.tar.gz
sudo cp g_mmpbsa/bin/g_mmpbsa /usr/local/bin/
sudo cp g_mmpbsa/bin/energy2bfac /usr/local/bin/
zenity --info --title="Successfull!!!" --text="Hey buddy you have done it!!!"
}
Help()
{
  ans=$(zenity --info --title 'gmx_qk 1.0.0' \
      --text 'Is there CUDA support available for MD simulation acceleration?' \
      --ok-label Quit \
      --extra-button Yes \
      --extra-button No --width 300 --height 100 )
  echo $ans
  if [[ $ans = "Yes" ]]
  then
    tar -zxvf gromacs-2021.4.tar.gz
    cd gromacs-2021.4
    mkdir build
    cd build
    cmake .. -DGMX_BUILD_OWN_FFTW=ON  -DGMX_GPU=CUDA -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda
    make -j8
    #make check -j8
    sudo make install
    cd ../../
    tar -zxvf gromacs-5.0.tar.gz
    cd gromacs-5.0
    mkdir build
    cd build
    cmake ..  -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=/opt/gromacs -DGMX_GPU=OFF
    make -j8
    #make check -j8
    sudo make install
    cd ../../
    install_gmx_qk;
  elif [[ $ans = "No" ]]
  then
    tar -zxvf gromacs-2021.4.tar.gz
    cd gromacs-2021.4
    mkdir build
    cd build
    cmake .. -DGMX_BUILD_OWN_FFTW=ON
    make -j8
    #make check -j8
    sudo make install
    cd ../../
    tar -zxvf gromacs-5.0.tar.gz
    cd gromacs-5.0
    mkdir build
    cd build
    cmake ..  -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=/opt/gromacs -DGMX_GPU=OFF
    make -j8
    #make check -j8
    sudo make install
    cd ../../
    install_gmx_qk;
  fi
  }
  
buil_essential

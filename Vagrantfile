# -*- mode: ruby -*-
# vi: set ft=ruby :

#Definir variable config
Vagrant.configure("2") do |config|
  #--- INICIO MAQUINA MPI ---
  #Configurar para iniciar ssh de servidor usando: vagrant ssh cl 
  config.vm.define "mpi" do |mpi|
    mpi.vm.box = "ubuntu/bionic64"
    mpi.vm.network "private_network", ip: "192.168.0.100"
    mpi.vm.provision "shell", inline: <<-SHELL
      echo "mpi" > /etc/hostname
      #Confiugrar hostname de la maquina
      hostname -b mpi
      #Instalar esenciales y actualizar
      sudo apt-get update -y
      sudo apt-get install build-essential -y
      #Compilador de c
      sudo apt-get install gcc -y
      #Librerias para MPI
      sudo apt-get install libopenmpi-dev -y
      sudo apt-get install openmpi-bin -y
      #Compilar
      #mpicc PmstMPI.c -o Pmst
      #Correr
      #mpirun -np 4 ./Pmst

    SHELL
  end
end

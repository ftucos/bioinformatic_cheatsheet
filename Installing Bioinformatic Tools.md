# Installing Bioinformatic Tools

## IGV

Download igv from https://software.broadinstitute.org/software/igv/download; if you are installing it on a remote server, install the command line version, and be sure to have java 11 installed 

```bash
# Donwnload IGV in the directory you want to keep all of your tools
cd ~/software
wget https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_2.8.10.zip
unzip IGV_2.8.10.zip
# run igv for the first time, it will create a new directory ~/.igv for the settings
./IGV_2.8.10/igv.sh
# edit ~/.igv/java_arguments to increase the startup ram
nano ~/.igv/java_arguments
# uncomment "-Xmx8G" and customize the startup ram required, usually the defauklt 4GB is fine for samll projects
```

Add igv to the path

Note: you need to specify absolute path, not relative like ```./igv.sh```

```bash
# make a softlink to /usr/local/bin/ or any other directory in the path
sudo ln -s ~/path/to/IGV_2.8.10/igv.sh /usr/local/bin/igv
```

### Using IGV on a server

To allow display forwarding and use the visual interface of igv, you should start your ssh connection enabling X11 forwarding with 

```bash
ssh -X user@host
# sometimes you may need to start a trusted forwarding to make it work
ssh -Y user@host
# you can verify that the forwarding is configured by checking the $DISPLAY variable was defined
```

#### Adjust screen resolution

you can check the screen resolution on the X11 server ```xdpyinfo | grep -e dimensions -e resolution```

and customize it with the *xrandr* package

```bash
# match your screen resolution
xrandr -q 
# if not working it will list you all tha availables options and you can pick one manualy 
```

## fast-dump and derivatives

Fast-dump, faster-dump and parallel-fasterq-dump to download reads from SRA by default will chache the downloaded files in ~/ncbi/public. This can be a problem if you have the home directory on a different and limited storage drive. Saturating that drive could eventually even lead you to not being able to access anymore trought ssh.

To change the cache location there is a GUI tool:

https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration

```bahs
vdb-config -i
```

## MISO

It requires a venv with python2.7 only

```bash
conda create env --name miso python=2.7
conda activate miso
pip install misopy
```

## RSEM

You need to compile RSEM from the github source code.

The install argument will move the installation in the `/usr/local/bin` directory

```bash
git clone https://github.com/deweylab/RSEM.git
cd RSEM
sudo make install
```

Yo may have problems with installation due to a custom compiler installed in conda (eg gxx_linux-64 required for compiling velocyto). In order to compiler RSEM it's recomended to `conda deactivate` in order to user system compiler

you may also need to install the developmental version of zlib ``sudo apt-get install -y zlib1g-dev``

## velocyto 

First of all install all the python dependencies

```bash
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
```

than install velocyto

`pip install --user velocyto`

If it fails, you may need a custom compiler

`conda install gxx_linux-64`

## MACS3

https://github.com/macs3-project/MACS/blob/master/docs/INSTALL.md

```bash
conda create -n macs3 python=3.8
conda activate macs3

conda install -c conda-forge numpy Cython cykhash
conda install -c bioconda fermi-lite
pip install macs3
```

## PECA

https://github.com/SUwonglab/PECA

1. Install whatever you need from conda

2. Install Matlab (full version, you need command line interpreted because the script runs code not compiled). You also need to install the optimization toolbox for the fsolve function.

3. If on MacOS, Add matlab to the path (in the `~/.bashrc` or  `~/.bash_profile`: 

   `export PATH="/Applications/MATLAB_R2022a.app/bin:$PATH"` and than comment in the bash script (PECA.sh,...) the `module load matlab` section (that makes sense only on linux)

   also added ` -p 8` to HOMER scripts to enable multicore processing.

4. chek to have the required genome (eg hg38) installed in homer: `less /Users/tucos/opt/anaconda3/envs/macs3/share/homer-4.10-0/.//config.txt` otherwise install it `perl /Users/tucos/opt/anaconda3/envs/macs3/share/homer-4.10-0/.//configureHomer.pl -install hg38`


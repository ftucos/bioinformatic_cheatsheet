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


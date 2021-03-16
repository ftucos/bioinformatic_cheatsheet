# Set Up new google cloud server

## Set yout password

configure your password with `sudo passwd`

## Mount DISK 1 in home folder

```bash
#find the disk
sudo lsblk
# format it as ext4
sudo mkfs.ext4 -m 0 -F -E lazy_itable_init=0,lazy_journal_init=0,discard /dev/sdb
# make the directory to mount it
sudo mkdir -p /mnt/disks/data
# mount it 
sudo mount -o discard,defaults /dev/sdb /mnt/disks/data
# give read-write permission
sudo chmod a+w /mnt/disks/data

# make a softlink to access from the home
ln -s /mnt/disks/data/ ~/userData
```

To automate the mount at the every boot of the machine, go on google cloud, edit the instance and add a custo metadata:

- Set `startup-script` as name and add the following script as value

  ```bash
  #!/usr/bin/bash
  sudo mount -o discard,defaults /dev/sdb /mnt/disks/data
  ```

## Enable ssh access from terminal or vscode

1. generate a key for the username you have on the server

   ```bash
   ssh-keygen -t rsa -f ~/.ssh/gcloud_bioinformatics -C [USERNAME ON THE SERVER]
   ```

2. copy the content of `~/.ssh/gcloud_bioinformatics.pub` in the public ssh keys on google cloud

3. now you can ssh specifing the custom key

   `ssh -i ~/.ssh/gcloud_bioinformatics tucci_bioinformatics@34.89.214.89`

4. to get access through VSCode you need to install the remote extension, click on it in the bottom-left corner and customize the configuration file

   ```
   Host gcloud_bioinformatics
     HostName 34.89.214.89
     User tucci_bioinformatics
     ForwardAgent yes
     ForwardX11Trusted yes
     IdentityFile ~/.ssh/gcloud_bioinformatics
   ```

## Change .bashrc

```bash
# Makes the prompt much more user friendly.
# But I do agree that the command to set it up looks a bit crazy.
export PS1='\[\e]0;\w\a\]\n\[\e[32m\]\u@\h \[\e[33m\]\w\[\e[0m\]\n\$ '

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions
alias ll='ls -lah'

## Safe versions of the default commands.
## Will ask permissions before overwriting files.
alias rm='rm -i'
alias mv='mv -i'
alias cp='cp -i'
```



## Install softwares

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# press space untill you have to input yes
rm -f Miniconda3-latest-Linux-x86_64.sh
source .bashrc

# add channels
conda config --append channels bioconda
conda config --append channels conda-forge

# install required packages
conda install r numpy scipy cython numba matplotlib gxx_linux-64 scikit-learn h5py click bamtools parallel samtools bowtie jupyter nextflow nf-core notebook picard sra-tools vcftools igv loompy trim-galore fastqc r-base seaborn statsmodels numba pytables -y

pip install multiqc scvelo velocyto fastqp scanpy MulticoreTSNE

# Install R-studio server
sudo apt-get update
sudo apt-get install gdebi-core
wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-1.4.1103-amd64.deb
# Confirm pressing y
sudo gdebi rstudio-server-1.4.1103-amd64.deb

# Install tldr, alternative to man
 sudo apt-get install -y tldr
```

#### Install Docker and enable it at startup

```bash
sudo apt install docker.io
# start it 
sudo systemctl start docker
# this would make docker autostart, but it would not work if you set the container image on an attached disk that gets mounted every time
## sudo systemctl enable docker
# add your user to the docker group otherwise you will also need to run it with sudo
sudo usermod -a -G docker $USER
# test it  
docker run hello-world
```

**Move docker image file location to your bif disk**

follow the tutorial

https://www.guguweb.com/2019/02/07/how-to-move-docker-data-directory-to-another-location-on-ubuntu/

## copy file from bucket

```bash
gsutil -m cp -r gs://bulk-rnaseq/* /mnt/SW480_variant_call/raw_data
```

## Organize the tool folder

```bash
mkdir -p resources/tools
mkdir -p resources/genomes
```

## Download reference genomes

- Don't include haplotype, this will simply increase the multimappers since it's not handled well by most aligner (except BWA?)
- Never use masked genomes, use softmasked instead
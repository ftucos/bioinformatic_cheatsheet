# GCLOUD computing RNAseq

1. Mount DISK 1 in home folder

   ```bash
   sudo mount /dev/sdb /mnt/
   # check with 
   df -lh
   # grand read/write permission
   sudo chmod -R a+w /mnt/
   ```

2. copy file from bucket

   ```bash
   mkdir -p /mnt/SW480_variant_call/raw_data
   gsutil -m cp -r gs://bulk-rnaseq/* /mnt/SW480_variant_call/raw_data
   
   mkdir -p /mnt/SW480_variant_call/sandbox /mnt/SW480_variant_call/logs /mnt/SW480_variant_call/scripts
   ```

3. Staff to install

   ```bash
   # GNU parallel
   sudo apt-get install parallel
   # PIP
   sudo apt install python3-pip
   # Cutadapt
   python3 -m pip install --user --upgrade cutadapt
   ```

4. Add local pip installation to the PATH

   ```bash
   # edit bashrc
   nano ~/.bashrc 
   # add the following string at the end of the file
   export PATH=$PATH:~/.local/bin
   ```

   
# Bash Notes

Useful resource: https://github.com/denysdovhan/bash-handbook

## Make a shell script

1. Never use sh for execution, allways bash. Sh was completely absorbed into bash, running a script with sh it's just hipster.

2. You don't need to specify the interpreted with file extension like "script.sh", you just need to specify the interpreter in the shebang line

   ```bash
   #!/usr/bin/bash
   echo "hello world"
   ```

   than make it executable  ```chmod +x my_script``` or ```chmod 755 my_script``` (if you don't want to give write permission to everybody)

   ![chmod Cheatsheet : linux](https://i.redd.it/vkxuqbatopk21.png)

#### Positional variables

```bash
#!bin/bash
echo "The script is: $0"
echo "Argument 1 is: $1"
# the brackets are required otherwise it will be interpreted as $1 followed by a 0
echo "Argument 10 is: ${10}"
echo "Passed a total number of $# arguments"
```

#### Named variables

Example from https://github.com/ludvigla/VisiumTrim/blob/main/TSO_polyA_trimming.sh

```bash
# Get script name
SCRIPT=$(basename "$0")

FASTQ=''
while (( "$#" )); do
  case "$1" in
    -l|--homopolymer-length)
      LENGTH=$2
      shift 2
      ;;
    -o|--output)
      OUTPUT=$2
      shift 2
      ;;
    -e|--error-tolerance)
      ETOL=$2
      shift 2
      ;;
    --overlap)
      OVERLAP=$2
      shift 2
      ;;
    -h|--help)
cat << EOF
$SCRIPT [-h] [-l -o -e --overlap n] 
...
```





#### Mandatory and optional variables

```bash
# mandatory argument
NAME1=${1?Error: missing name}
# optional argument (default value is "friend")
NAME2=${2:-friend}

echo "Hi $NAME1 and ${NAME2}!"
```

#### Make all the scripts in a folder executable

```bash
chmod +x *
```

#### Relative vs absolute shebang

```bash
#!/usr/bin/env python

#!/whatever/absolute/path/is/my/python
```

It makes no sense to use it with bash since bash is allways in `/bin/bash` but it's really usefull with python since everybody has multiple and diferent python installation and you cannot guess the path (especially since most bioinformatician (and me) prefer using conda python which is in the user home directory thus the path is user specific)

#### Writing safe shell scripts

https://sipb.mit.edu/doc/safe-shell/

Quits the whole script if any errors happens

```bash
#!/usr/bin/env bash
set -eux -o pipefail

command_allowed_to_fail || true
# or simply
command_allowed_to_fail ||
```

If you need to run 

```bahs
source ~/.bashrc
conda acrivate my_env
```

You will need to set -ue after that 



### ask confirmation before executing

From: https://stackoverflow.com/questions/1885525/how-do-i-prompt-a-user-for-confirmation-in-bash-script

```bash
read -p "Do you want to run variant filtration [y/N]" -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1 # handle exits from shell or function but don't exit interactive shell
fi
```

**Explanation:**

The `read` command outputs the prompt (`-p "prompt"`) then accepts one character (`-n 1`) and accepts backslashes literally (`-r`) (otherwise `read` would see the backslash as an escape and wait for a second character). The default variable for `read` to store the result in is `$REPLY` if you don't supply a name like this: `read -p "my prompt" -n 1 -r my_var`

The `if` statement uses a regular expression to check if the character in `$REPLY` matches (`=~`) an upper or lower case "Y". The regular expression used here says "a string starting (`^`) and consisting solely of one of a list of characters in a bracket expression (`[Yy]`) and ending (`$`)". The anchors (`^` and `$`) prevent matching longer strings. In this case they help reinforce the one-character limit set in the `read` command.

The negated form uses the logical "not" operator (`!`) to match (`=~`) any character that is not "Y" or "y". An alternative way to express this is less readable and doesn't as clearly express the intent in my opinion in this instance. However, this is what it would look like: `if [[ $REPLY =~ ^[^Yy]$ ]]`



## Globing

* ```*``` matches every character one or more times excepth ```\```. Useful for catching file names but noth paths
* ```**``` matches everything including ```\```. Useful for catching also the path
* ```?``` matches everything excepth ```\``` only one time, useful for catching file names but noth paths
* ```photo[1234].jpg```  or ```photo[1-4].jpg``` matches photo1.jpg, photo2.jpg... photo4.jpg

## Brace expansion

```bash
touch note1.txt note2.txt note3.txt
# or
touch note{1,2,3}.txt
# or
touch note{1..3}.txt

# Only pari numbers
touch note{00..4..2}.txt 
# note00.txt note02.txt note04.txt
```

## Command substitution

```bash
now=`date +%T`
# or
now=$(date +%T)

echo $now # 19:08:26
```

![img](https://wizardzines.com/comics/brackets-cheatsheet/brackets-cheatsheet.png)



## Variable expansion and quotes

quotes are usevul for avoid ambiguity with spaces

```bash
# This makes two separate folders
mkdir empy folder
# to escape the space you have several options
mkdir "empty folder"
# don't get confused with ` type of single quote
mkdir 'empty folder'
mkdir empty\ folder
```

the difference between ```'``` and ```"``` is that only double quotes allows variable and command expansion (but not globing) while single quotes protects everything

```bash
NAME="Francesco"
echo "My name is $NAME"
#=> My name is Francesco
echo 'My name is $NAME'
#=> My name is $NAME
echo "$NAME has \$5 in his pocket"
#=> Francesco has $5 in his pocket
echo "${NAME}one"
#=> Francescone
```

`${VARIABLE}other characters It is optional but serve to protect the variable to be expanded from characters immediately following it which could be interpreted as part of the name.



Remember all types of quotes protects brace expansion (both double and single)

```bash
touch file 1.txt
ls
#=> file      1.txt
touch "file 1.txt" "file 2.txt" "file 3.txt"
ls
#=> file 1.txt	file 2.txt	file 3.txt
touch "file {1..3}.txt"
ls
#=> file {1..3}.txt
touch file\ {1..3}.txt
ls
#=> file 1.txt	file 2.txt	file 3.txt
```

and globing

```bash
ls file*.txt
#=> file 1.txt	file 2.txt	file 3.txt
ls "file*.txt"
#=> ls: file*.txt: No such file or director"
$FILE_NUMBER=3
ls "file 1.txt" "file 2.txt" "file ${FILE_NUMBER}.txt"
#=> file 1.txt	file 2.txt	file 3.txt
```



## Create a softlink in the $PATH to allow and app/script to be lunched by the shell

```
ln -s /Users/tucos/Tools/igv/igvtools /usr/local/bin/igvtools
```

This operation is permitte only for the /usr/local/bin/ and not for /usr/bin/

## Loop over a list of files

```bash
for f in file1 file2 file3 file5
do
 echo "Processing $f"
done
```

```bash
FILES="file1
/path/to/file2
/etc/resolv.conf"
for f in $FILES
do
	echo "Processing $f"
done
```

**Alternative formatting**

```bash
for f in *.bam
do echo "Processing $f"
done

for f in *.bam; do
  echo "Processing $f"
done

# semicolon can be used in place of line break
for f in *.bam; do echo "Processing $f file.."; done
```

#### Parallelize on multicore

```bash
for f in *.bam; do
	echo "Processing $f" & done
```

This approach will take advantage of all the cores available and overload your machine, a better and more controlled approach is to use xargs or GNU parallel

## Parallelize for a list of file with GNU Parallel

GNU parallel is a bit slower (+3 ms/job vs +0.3 ms/job) but more complete alternative to `xargs`. Specifically, GNU parallel handles better special characters as excapes, pipes and newlines.

documentation: https://www.gnu.org/software/parallel/parallel_tutorial.html

``` bash
find ./processed/bam -name "*.sorted.bam" | parallel samtools index {}

# run long script
find ./processed/bam -name "*.sorted.bam" | parallel 'samtools index"{}" | head'
# Print to consol the output of the first job in progress
find ./processed/bam -name "*.sorted.bam" | parallel -k --lb 'samtools index"{}" | head'

{.} # Filename without extension
{/} # Filename without the path
{//} # Path without filename
{/.} # Remove path and extension
{#} # job number: 1 2 3 ...
{=s/pattern/replacement/=}	# sed replacement
```

#### string replacement using regex characters (and capturing groups)

```bash
VAR="file.edited_07.txt"
parallel echo '{=s/\.edited_([0-9]+)\.txt/final_\1/=}.tsv' ::: $VAR
```



### Quoting in GNU Parallel

GNU Parallel is very liberal with quoting, you need to quote only special shell characters as ```  ( ) $ ` ' " < > ; | \``` and depending on the context ```  ~ & # ! ? space * {```

```bash
parallel “zcat {} | bzip2 >{.}.bz2” ::: *
# Can also be writteng quoting only shell special characters
parallel zcat {} “|” bzip2 “>”{.}.bz2 ::: *
```

-q ??

This is really important when you want to execute mutliple commands

```bash
parallel “command1 {}; command2 {}” ::: *
# or
parallel command1 {}";" command2 {} ::: *
```

### Remove Path from find

```bash
find ./ -name "*Aligned.sortedByCoord.out.bam" | sed "s/.*\///" | parallel echo {}
```

## Tricks with variables

### Edit a variable with SED

```bash
# remove the extension
FILE_NAME=file.txt
#{object//pattern/replacement}
FILE_NAME=${FILE_NAME//txt/tsv}
echo $FILE_NAME
```

### Assign default value to a variable

Set a variable only if it's unset

```bash
# set ${VAR=value} is the same as VAR=value, but allows
# you to use also other assignment operators as :=
set ${NAME:=sconosciuto}
echo $NAME
# sconosciuto
NAME=Mario
echo ${NAME:=sconosciuto}
# Mario

# yes, you can assign a variable on the fly!
echo ${NAME:=Luca}
# Luca
```

**Tip: ${var:-defaultValue} vs ${var:=defaultValue}**

```bash
echo ${SURNAME:-Rossi}
# Rossi

# it returns a default value without assigning it to the variable
echo ${SURNAME}
#

# thus it makes no sense to set with :-
```

### Remove Pattern 

**Front of $VAR**

```bash
var=path/to/my/file.txt
# one # removes the shortest match
echo ${var#*/}
# to/my/file.txt

# Two ## removes the longest match
echo ${var##*/}
# file.txt

# this is an alternative approach to using POSIX basename
echo $(basename -- "$var")
```

**End of $VAR**

use `%` and `%%` the same way

### Summary: String manipulation

| Pattern                       | Function                                                     |
| ----------------------------- | ------------------------------------------------------------ |
| ${parameter:-defaultValue}    | Get default shell variables value                            |
| ${parameter:=defaultValue}    | Set default shell variables value                            |
| ${parameter:?"Error Message"} | Display an error message if parameter is not set             |
| ${#var}                       | Find the length of the string                                |
| ${var%pattern}                | Remove from shortest rear (end) pattern                      |
| ${var%%pattern}               | Remove from longest rear (end) pattern                       |
| ${var:num1:num2}              | Substring                                                    |
| ${var#pattern}                | Remove from shortest front pattern                           |
| ${var##pattern}               | Remove from longest front pattern                            |
| ${var/pattern/string}         | Find and replace (only replace first occurrence)             |
| ${var//pattern/string}        | Find and replace all occurrences                             |
| ${!prefix*}                   | Expands to the names of variables whose names begin with prefix. |
| ${var,} ${var,pattern}        | Convert first character to lowercase.                        |
| ${var,,} ${var,,pattern}      | Convert all characters to lowercase.                         |
| ${var^} ${var^pattern}        | Convert first character to uppercase.                        |
| ${var^^} ${var^^pattern}      | Convert all character to uppercase..                         |

## Rename multiple fiels

```bash
for i in *string_to_edit*.txt
do
  mv $i ${i/string_to_edit/string_edited}
done
```

## Add a location to the PATH

```bash
# add this function to the ~/.zshrc so you don't have to call it every time you open a new terminal
export PATH=$PATH:$HOME/.local/bin
```

## mkdir complete path

```bash
 mkdir -p ./variant_call/output/FastQC2
```

#### check if a directory exist

```bash
DIR="/etc/httpd/"
if [ -d "$DIR" ]; then
  ### Take action if $DIR exists ###
  echo "Installing config files in ${DIR}..."
else
  ###  Control will jump here if $DIR does NOT exists ###
  echo "Error: ${DIR} not found. Can not continue."
  exit 1
fi
```

## Write console output to a file

#### Overview:

> *Please note that the `n.e.` in the syntax column means "not existing".*
> There is a way, but it's too complicated to fit into the column. You can find a helpful link in the List section about it.

```
          || visible in terminal ||   visible in file   || existing
  Syntax  ||  StdOut  |  StdErr  ||  StdOut  |  StdErr  ||   file   
==========++==========+==========++==========+==========++===========
    >     ||    no    |   yes    ||   yes    |    no    || overwrite
    >>    ||    no    |   yes    ||   yes    |    no    ||  append
          ||          |          ||          |          ||
   2>     ||   yes    |    no    ||    no    |   yes    || overwrite
   2>>    ||   yes    |    no    ||    no    |   yes    ||  append
          ||          |          ||          |          ||
   &>     ||    no    |    no    ||   yes    |   yes    || overwrite
   &>>    ||    no    |    no    ||   yes    |   yes    ||  append
          ||          |          ||          |          ||
 | tee    ||   yes    |   yes    ||   yes    |    no    || overwrite
 | tee -a ||   yes    |   yes    ||   yes    |    no    ||  append
          ||          |          ||          |          ||
 n.e. (*) ||   yes    |   yes    ||    no    |   yes    || overwrite
 n.e. (*) ||   yes    |   yes    ||    no    |   yes    ||  append
          ||          |          ||          |          ||
|& tee    ||   yes    |   yes    ||   yes    |   yes    || overwrite
|& tee -a ||   yes    |   yes    ||   yes    |   yes    ||  append
```

#### List:

- `command > output.txt`

  The standard output stream will be redirected to the file only, it will not be visible in the terminal. If the file already exists, it gets overwritten.

- `command >> output.txt`

  The standard output stream will be redirected to the file only, it will not be visible in the terminal. If the file already exists, the new data will get appended to the end of the file.

- `command 2> output.txt`

  The standard error stream will be redirected to the file only, it will not be visible in the terminal. If the file already exists, it gets overwritten.

- `command 2>> output.txt`

  The standard error stream will be redirected to the file only, it will not be visible in the terminal. If the file already exists, the new data will get appended to the end of the file.

- `command &> output.txt`

  Both the standard output and standard error stream will be redirected to the file only, nothing will be visible in the terminal. If the file already exists, it gets overwritten.

- `command &>> output.txt`

  Both the standard output and standard error stream will be redirected to the file only, nothing will be visible in the terminal. If the file already exists, the new data will get appended to the end of the file..

- `command | tee output.txt`

  The standard output stream will be copied to the file, it will still be visible in the terminal. If the file already exists, it gets overwritten.

- `command | tee -a output.txt`

  The standard output stream will be copied to the file, it will still be visible in the terminal. If the file already exists, the new data will get appended to the end of the file.

- **(\*)**

  Bash has no shorthand syntax that allows piping only StdErr to a second command, which would be needed here in combination with `tee` again to complete the table. If you really need something like that, please look at ["How to pipe stderr, and not stdout?" on Stack Overflow](https://stackoverflow.com/q/2342826/4464570) for some ways how this can be done e.g. by swapping streams or using process substitution.

- `command |& tee output.txt`

  Both the standard output and standard error streams will be copied to the file while still being visible in the terminal. If the file already exists, it gets overwritten.

- `command |& tee -a output.txt`

  Both the standard output and standard error streams will be copied to the file while still being visible in the terminal. If the file already exists, the new data will get appended to the end of the file.

## Remove multiple files in fancy ways

E.g. we have to remove some files from this list:

```bash
abc.log.2012-03-14
abc.log.2012-03-27
abc.log.2012-03-28
abc.log.2012-03-29
abc.log.2012-03-30
abc.log.2012-04-02
abc.log.2012-04-04
abc.log.2012-04-05
abc.log.2012-04-09
abc.log.2012-04-10
```

#### Wildcard

```bash 
rm -f abc.log.*             # Remove them all
rm -f abc.log.2012*         # Remove all logs from 2012
rm -f abc.log.2012-0[123]*  # Remove all files from the first quarter of 2012
```

#### RegEX

Regular expressions are more powerful than wildcards; you can feed the output of `grep` to `rm -f`. For example, if some of the file names start with `"abc.log"` and some with `"ABC.log"`, `grep`lets you do a case-insensitive match:

```sh
rm -f $(ls | grep -i '^abc\.log\.')
# equivalent to
 rm -f `ls | grep -i '^abc\.log\.'`
```

This will cause problems if any of the file names contain funny characters, including spaces. Be careful.

When I do this, I run the `ls | grep ...` command first and check that it produces the output I want -- *especially* if I'm using `rm -f`:

```sh
$ ls | grep -i '^abc\.log\.'
(check that the list is correct)
$ rm -f $(!!)
```

where `!!` expands to the previous command. Or I can type up-arrow or Ctrl-P and edit the previous line to add the `rm -f` command.

#### File List

```bash
ls > list
rm -f $(<list)
# or
rm -f `cat list`
```

## RegEX string replacement with SED

```<<<``` is required to pass the path as a string otherwise sed will try to open the file with the given path

```bash
# extracts Index from Illumina TruSeq sequencing file name
INDEX=`sed -E 's/.*\([ACTG]\{6\}\).*/\1/' <<< "Sample_SW480_spheres/SW480_spheres_CGATGT_L006_R2_001.fastq.gz"`
```

## Extract filename/path using POSIX

```bash
echo $(basename -- "/path/to/file.txt")
# returns: file.txt

echo $(dirname -- "/path/to/file.txt")
# returns: /path/to
```

Note: the `--` are not mandatory but recomended to avoid problem in case  the variable starts with a dash (otherwise it would be interpreted as an option and not an argument)

## The use of parenthesis 

https://dev.to/rpalo/bash-brackets-quick-reference-4eh6

#### ( single parentheses )

runs script in a subshell, the variables are deleted at the end of the script

```bash
( a=banana; mkdir $a )
echo $a
# => 
```

#### <( process substitution )

Creates a temporary file `/dev/fd/<n>` deleted after the end of the process.

```bash
gunzip -c file.txt.gz | grep "my_string"
grep <(gunzip -c file.txt.gz) "my_string"
```

#### if [] if [[]] for evaluation of truth

```
# Single-line
if [[ 2 -ne 1 ]]; then echo "true"; else echo "false"; fi

# Multi-line
if [[ 2 -ne 1 ]]; then
  echo "true"
else
  echo "false"
fi
```

Expressions enclosed inside `[[ ]]` (or `[ ]` for `sh`) are called **test commands** or **primaries**. These expressions help us to indicate results of a conditional. In the tables below, we are using `[ ]`, because it works for `sh` too. Here is an answer about [the difference between double and single square brackets in bash](http://serverfault.com/a/52050).

**Working with the file system:**

| Primary               | Meaning                                                     |
| --------------------- | ----------------------------------------------------------- |
| `[ -e FILE ]`         | True if `FILE` **e**xists.                                  |
| `[ -f FILE ]`         | True if `FILE` exists and is a regular **f**ile.            |
| `[ -d FILE ]`         | True if `FILE` exists and is a **d**irectory.               |
| `[ -s FILE ]`         | True if `FILE` exists and not empty (**s**ize more than 0). |
| `[ -r FILE ]`         | True if `FILE` exists and is **r**eadable.                  |
| `[ -w FILE ]`         | True if `FILE` exists and is **w**ritable.                  |
| `[ -x FILE ]`         | True if `FILE` exists and is e**x**ecutable.                |
| `[ -L FILE ]`         | True if `FILE` exists and is symbolic **l**ink.             |
| `[ FILE1 -nt FILE2 ]` | FILE1 is **n**ewer **t**han FILE2.                          |
| `[ FILE1 -ot FILE2 ]` | FILE1 is **o**lder **t**han FILE2.                          |

**Working with strings:**

| Primary            | Meaning                                          |
| ------------------ | ------------------------------------------------ |
| `[ -z STR ]`       | `STR` is empty (the length is **z**ero).         |
| `[ -n STR ]`       | `STR` is not empty (the length is **n**on-zero). |
| `[ STR1 == STR2 ]` | `STR1` and `STR2` are equal.                     |
| `[ STR1 != STR2 ]` | `STR1` and `STR2` are not equal.                 |

**Arithmetic binary operators:**

| Primary             | Meaning                                            |
| ------------------- | -------------------------------------------------- |
| `[ ARG1 -eq ARG2 ]` | `ARG1` is **eq**ual to `ARG2`.                     |
| `[ ARG1 -ne ARG2 ]` | `ARG1` is **n**ot **e**qual to `ARG2`.             |
| `[ ARG1 -lt ARG2 ]` | `ARG1` is **l**ess **t**han `ARG2`.                |
| `[ ARG1 -le ARG2 ]` | `ARG1` is **l**ess than or **e**qual to `ARG2`.    |
| `[ ARG1 -gt ARG2 ]` | `ARG1` is **g**reater **t**han `ARG2`.             |
| `[ ARG1 -ge ARG2 ]` | `ARG1` is **g**reater than or **e**qual to `ARG2`. |

Conditions may be combined using these **combining expressions:**

| Operation            | Effect                                                   |
| -------------------- | -------------------------------------------------------- |
| `[ ! EXPR ]`         | True if `EXPR` is false.                                 |
| `[ (EXPR) ]`         | Returns the value of `EXPR`.                             |
| `[ EXPR1 -a EXPR2 ]` | Logical *AND*. True if `EXPR1` **a**nd `EXPR2` are true. |
| `[ EXPR1 -o EXPR2 ]` | Logical *OR*. True if `EXPR1` **o**r `EXPR2` are true.   |

#### More

* ```((...))```used of arithmetic operations (outputs nothing)
* ```$(...)``` let's integrate the result of a subprocess into a string

* ```$((...))``` outputs the result of an arithmetic operation

* ```$(...)$?``` outputs the exit code of a subprocess



### File test operation

```bash
file=~/Download/file.txt
if [-e "$file"]
then
	echo "file.txt was downloaded"
else
	echo "file.txt missing"
fi
```

Multiple evaluation options: 

* `-e` a file/folder/device exist; `-f` for files; `-d` for directories; 

* `-s` file size is not 0
* `-r`, `-w` and `-x` to evaluate the read, write and execute permission

## Keep track of changes in a directory in real time

show the current directory and file sizes updating every 1 second

```bash
watch -n1 tree . --du -h
```

## Convert pictures

Batch convert keeping the same file structure

```bash
# requires imagemagick
for f in **/*.tif
do
echo "Converting $f"
convert $f $(dirname "$f")/$(basename "$f" .tif).jpg
rm $f
done
```

## Proces management

```bash
# Start a long running process just for the purpose of testing
sleep 300
# STOP a process
^C
# SUSPEND a process
^X
# list jobs (suspendend or in background), -l show the ProcesID in an additional column
jobs -l
# or
ps -f
# Restart a process in foreground
fg
# or specify the job number if you have multiple
fg %2
# restart it in background
bg
# start a process directly in bacground
sleep 300 &
```

all the jobs live withing current terminal session, if you close it you will close all the running ans suspended job started in that terminal

### how to terminate a job/program

```bash
# SIGTERM: sends to the program a request of terminating is process, the program can thus handle the termination closing some files handling requests, removing some temporary files etc, or even ignoring the termination request.
kill -15 <PID> 
# or simply 
kill <PID>

# SIGKILL: the kernel kills the process in an harsh and inevitable way. If you kill the writing of a file, that file can end up to be corrupted
kill -9 <PID>
```

Actually kill is not only terminating process but can also suspend and restart them

```bash
# SIGTSTP (the signal sent when you press ^Z)
kill -20 <PID>
# SIGCONT (equivalent to bg on the last process)
kill -18 <PID>
# list all the kill options
kill -l
```

## Use screen sessions

Screen session are purposeless on your local machine with a GUI since you can already open multiple terminal and or TABs. It becames crucial when you are on an ssh connection because if you loose the connection the ongoing process will be arrested. If you start a screen session instead, if you loose the connection/close the terminal you get only detached and the prosess keeps going on

Some useful commands:

```bash
# start a session
screen 
# or specify a name for it
screen -S new_session

# list all the available commands once you are in the session: press CTRL(^)+A followed by the question mark 
^A ? 

# Detatch the session
^A d
# List all the sessions
screen -ls
# recover (attach) the session: you can use the session id or the session name (if you specified one)
screen -r new_session
# recover a screen session aleready attached somewhere else
screen -rd new_session

# Create a new window into the same session
^A c
# Rename the current window
^A A
# List available windows into the current session and switch to a different one
^A "
# cycle between available windows
ˆA space
# split horizontally the screen and than use ^A " to select the new screen
^A S
# swith between the different parts of the splitted screen
^A TAB
# close the current slice
^A X
# close all the slices aside the current one
^A Q
# exit and close the current session
ˆD

#Close a dedatched session
screen -XS test_screen quit
```

To scroll in screen mode:

```bash
^A Esc
```

## Functions

```bash
my_func () {
  # statements
}

my_func # call my_func
```

```bash
# function with params
greeting () {
  if [[ -n $1 ]]; then
    echo "Hello, $1!"
  else
    echo "Hello, unknown!"
  fi
  return 0
}

greeting Denys  # Hello, Denys!
greeting        # Hello, unknown!
```

### Keep track of occpied space in a folder

```bash
watch -n1 du -h --max-depth 1
```

## Google Cloud

### Download a bucket

1. Be sure to have guts installed

2. Be shure to be logged in the account owner of the bucket (it not public) otherwise

   ```bash
   gcloud init --console-only
   ```

3. Download the content of your bucket

   ```bash
   # -r to download all the files into the directory
   gutils cp -r gs://my-bucket ~/destination/directory/
   # or specify only which type of files to download
   gutils cp -r *.txt gs://my-bucket ~/destination/directory/
   ```

   

## Batch convert PPTX to PDF

Cannot parallelize `unoconv` because it will mess with temporary files generated in the destination folder*. To make it work with GNU parallel, you shuld export everything in subdirectories*

```bash
for i in PPTX/*.ppt*
do
	file=${i//**\//}
	unoconv -f pdf -o "PDF/${file//ppt*/pdf}" PPTX/"$file"
done
```

## Git

#### start tracking an existing directory

```bash
git init
git add ./scripts/
git commit -m 'message'
git remote add origin <url>
git push -u origin master
```

#### Remove .DS_store from your repo

https://stackoverflow.com/questions/107701/how-can-i-remove-ds-store-files-from-a-git-repository

Remove existing files from the repository:

```
find . -name .DS_Store -print0 | xargs -0 git rm -f --ignore-unmatch
```

Add the line `.DS_Store` to the file `.gitignore`, which can be found at the top level of your repository (or created if it isn't there already). You can do this easily with this command in the top directory

```
echo .DS_Store >> .gitignore
```

Then

```
git add .gitignore
git commit -m '.DS_Store banished!'
```

### Setup remote key

```bash
# generate an rsa key (dsa not supported by GitHub)
ssh-keygen -t rsa
# a prompt will ask you for the name: e.g. 
> github_ftucos
# copy the public key to the clipboard
pbcopy < github_ftucos.pub
# Add it to your public access key in github (in settings)
# edit your ~/.ssh/config by adding:
Host github_as_ftucos
  HostName github.com
  User git
  IdentityFile ~/.ssh/github_ftucos
  IdentitiesOnly yes
```

When you add a remote to your github repo, replace github.com with the host name in the `.ssh/config`

```bash
# replace
git remote add origin git@github.com:ftucos/repo-name.git
# with
git remote add origin git@github_as_ftucos:ftucos/repo-name.git
```


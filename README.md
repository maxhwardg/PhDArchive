# Max Ward's PhD Code Archive
This repository is an archive of code and data associated with Max Ward's PhD thesis. Some of the same code and data 
has appeared previously at https://github.com/maxhwardg/advanced_multiloops which is associated with some publications.

Several programs and a data set are included. The programs are folding algorithms, energy calculators, and parameter 
training algorithms for the models explored in my thesis. These include the linear, logarithmic, Aalberts & Nandagopal, 
average asymmetry, and stem length penalties models. These programs can be found in the programs/ directory. They use 
a library of useful code found in librnary/. The data set is contained in the data_set directory.

The code contains some parts of RNAstructure 5.8, which I have been given permission to use, and which is licensed under
 GPL. The data set source is described in my thesis. It is a modification of the ArchiveII data set from the Mathews 
lab (found at http://rna.urmc.rochester.edu/pub/archiveII.tar.gz). The tRNAs have been replaced by those from 
RNAstralign (https://rna.urmc.rochester.edu/pub/RNAStralign.tar.gz). Extensive cleaning has been done to the data 
set to remove duplicate sequences, partial structures, and suspicious structures.

## Compilation & Requirements
To compile the programs you will need a few things. The first is a C++ compiler that supports the C++11 standard, 
and the second is CMake (version >= 3.5.x), which is the build system used. The code has only been tested using the 
default GCC (G++), and CMake installation under Ubuntu 16.04, though anything meeting the previous requirements should 
work. With this in mind, compilation can be done the usual CMake way. For example, let's assume we are using GNU Linux 
system, and we're in the root directory for this project:

```
mkdir bin
cd bin
cmake ../ -DCMAKE_BUILD_TYPE=Release
make all
```

Now, let's run the tests to make sure everything is working. The following assumes you are still in the bin directory.

```
cd librnary_tests/
./librnary_tests
```

## A Note about data_tables
A common usage mistake is incorrectly providing the relative path to the data_tables directory to the executables. 
If you are getting strange looking structures with odd free energies, this is the most likely cause.

## Folding Programs
All the folding programs have the form fold_*. They all support usual command line flags, and usage information via 
"--help" can be seen. Let's run through an example usage.

Go to the root directory of this repository. Also, I will assume the project was built in the bin directory.

```
./bin/programs/fold_linear --help
```

This should bring up some useful information. Now let's try folding an RNA. Try:

```
./bin/programs/fold_linear
AUCGUAGCUAGUCUGAGCUAUCGUAGUCUGAGUCUGCA
```

This should yield the output:
```
...((((((......)))))).((((.......)))).
MFE: -9.8 (kcal/mol)
```

Now, we will try folding an RNA specified by a CT file. To do this, I use the read_cts program, which is useful for 
extracting CT information to be piped to other programs.

```
echo  "data_set/ct_files/tRNA_tdbD00012070.ct" | ./bin/programs/read_cts  | ./bin/programs/fold_linear
```

Should result in:

```
((((((((.((...))..((((....(((((((.......)))))))))))((((.......)))))))))))).
MFE: -36.9 (kcal/mol)
```

Other folding programs work in the same way. However, the various model specific parameters can be modified via 
command line.

## Energy Calculators
All the energy calculator programs have the form energy_*. They all support usual command line flags, and usage 
information via "--help" can be seen. Let's run through an example usage.

Go to the root directory of this repository. Also, I will assume the project was built in the bin directory.

```
./bin/programs/energy_logarithmic --help
```

This should bring up some useful information. Now let's try some free energy calculation an a demo RNA. I will use the 
folding output from the linear folding algorithm in the previous section:

```
./bin/programs/energy_logarithmic
AUCGUAGCUAGUCUGAGCUAUCGUAGUCUGAGUCUGCA
...((((((......)))))).((((.......)))).
```

Should result in:

```
Free Energy Change: -9.8 (kcal/mol)
```

Now, let's try out the verbose mode, which gives a complete breakdown of the energy calculation. This includes 
accoutrements.

```
./bin/programs/energy_logarithmic -v
AUCGUAGCUAGUCUGAGCUAUCGUAGUCUGAGUCUGCA
...((((((......)))))).((((.......)))).
```

Should result in:

```
Free Energy Change: -9.8 (kcal/mol)
External-loop: -32/-98
AU/GU penalties: 5
Stacking: -37
Stacking(22, 36, mismatch-mediated coax (3,20)) 
 Two-loop closed by (3, 20) and (4, 19): -14/-54
  Two-loop closed by (4, 19) and (5, 18): -13/-40
   Two-loop closed by (5, 18) and (6, 17): -21/-27
    Two-loop closed by (6, 17) and (7, 16): -34/-6
     Two-loop closed by (7, 16) and (8, 15): -21/28
      One-loop closed by (8, 15): 49/49
 Two-loop closed by (22, 36) and (23, 35): -25/-12
  Two-loop closed by (23, 35) and (24, 34): -10/13
   Two-loop closed by (24, 34) and (25, 33): -21/23
    One-loop closed by (25, 33): 44/44
```

Finally, we will use the energy calculator on a CT file. This works similarly as for the folding programs. However, 
the -s flag is used to get the secondary structure from a CT file.

Doing,

```
echo  "data_set/ct_files/tRNA_tdbD00012070.ct" | ./bin/programs/read_cts -s | ./bin/programs/energy_logarithmic -v
```

Should output,

```
Free Energy Change: -35 (kcal/mol)
External-loop: -17/-350
AU/GU penalties: 0
Stacking: -17
Stacking(0, 73, 3' dangle) 
 Two-loop closed by (0, 73) and (1, 72): -33/-333
  Two-loop closed by (1, 72) and (2, 71): -24/-300
   Two-loop closed by (2, 71) and (3, 70): -22/-276
    Two-loop closed by (3, 70) and (4, 69): -33/-254
     Two-loop closed by (4, 69) and (5, 68): -33/-221
      Two-loop closed by (5, 68) and (6, 67): -24/-188
       Multi-loop closed by (6, 67): 28/-164
       AU/GU penalties: 5
       Multi-loop closure: 71
       Multi-loop closure featues: branches=4 unpaired=6 
       Stacking: -48
       Stacking(9, 26, flush coax (27,45)) Stacking(6, 67, flush coax (50,66)) 
        Two-loop closed by (9, 26) and (10, 25): -34/-29
         Two-loop closed by (10, 25) and (11, 24): -21/5
          Two-loop closed by (11, 24) and (12, 23): -24/26
           One-loop closed by (12, 23): 50/50
        Two-loop closed by (27, 45) and (28, 44): -25/-84
         Two-loop closed by (28, 44) and (29, 43): -21/-59
          Two-loop closed by (29, 43) and (30, 42): -21/-38
           Two-loop closed by (30, 42) and (31, 41): -33/-17
            Two-loop closed by (31, 41) and (32, 40): -34/16
             One-loop closed by (32, 40): 50/50
        Two-loop closed by (50, 66) and (51, 65): -33/-79
         Two-loop closed by (51, 65) and (52, 64): -24/-46
          Two-loop closed by (52, 64) and (53, 63): -33/-22
           Two-loop closed by (53, 63) and (54, 62): -33/11
            One-loop closed by (54, 62): 44/44
```

As with the folding algorithm programs, the energy programs can be given flags to adjust the parameters of their 
associated models.

## Parameter Training Algorithms
The parameter training programs are train_linear, train_logarithmic, and train_an. They all have similar input requirements. Instructions for flags can be found by calling a program with the flag "-h".

### Basic Usage
All parameter training algorithms need three pieces of input to function. The path to the data tables, which are the parameters for the entire nearest neighbor model as per RNAstructure; the path to the data set, which is a folder containing all the .ct files for RNAs used in training; a list of .ct file names to use for training. If run from the root directory of this repository, the default paths for the data tables and data set will work out of the box. To change them, use the flag "-h" to see how. The list of .ct files to use should be given via standard input. The folder "data_sets/" contains several .ctset files. These are the lists of .ct files used for the paper. The folder "data_sets/ct_files" contains the archive of training data used in the paper.

The range of parameters to optimize over is hard coded in each program. They are those used in the paper. To change the parameter ranges, the code for a program must be changed.

### An Example
Let us consider an example run of parameter optimization. For this example, I will assume that we are in the root directory of the project. Also, I will assume that "path_to_programs/" is the path to the folder the programs are in. The following should bring up usage information for the linear parameter training program.

```
./bin/programs/train_linear -h
```

It should look a bit like this:

```
Trains the parameters of the linear model using IBF. Expects a .ctset file as input on standard in. For some .ctset files, see the data_set directory.
Usage:
  Train Linear Model [OPTION...]

  -d, --data_path arg  Path to data_tables (default: data_tables/)
  -c, --ct_path arg    Path to the folder of CTs (default:
                       data_set/ct_files/)
  -t, --threads arg    Number of threads to use (default: 16)
  -h, --help           Print help
```

Now, let us run the program using the small training set. By default, this will use all the cores on your machine, so be wary! The following invocation should starting training:

```
./bin/programs/train_linear < data_set/small.ctset
```

It can take a while. You should see something like the following as a result:

```
Seed #1 init = 32 branch = -58 unpaired = 3: 0.101219
Seed #2 init = 112 branch = -38 unpaired = -4: 0.472647
Seed #3 init = 183 branch = -9 unpaired = -16: 0.096862
Seed #4 init = 125 branch = -44 unpaired = 2: 0.463108
Seed #5 init = 193 branch = 5 unpaired = -5: 0.464895
Epoch #0: 
	Average F-Score = 0.0342677
	Best score = 0.768073
	init = 109 branch = 30 unpaired = 11
Epoch #1: 
	Average F-Score = 0.417572
	Best score = 0.630142
	init = 75 branch = -7 unpaired = 0
Epoch #2: 
	Average F-Score = 0.617731
	Best score = 0.636124
	init = 163 branch = -25 unpaired = -2
Epoch #3: 
	Average F-Score = 0.590204
	Best score = 0.629207
	init = 115 branch = -15 unpaired = 0
Epoch #4: 
	Average F-Score = 0.631208
	Best score = 0.633857
	init = 111 branch = -18 unpaired = 1
Epoch #5: 
	Average F-Score = 0.625122
	Best score = 0.631451
	init = 110 branch = -10 unpaired = -1
Epoch #6: 
	Average F-Score = 0.632576
	Best score = 0.635146
	init = 124 branch = -14 unpaired = -1
Epoch #7: 
	Average F-Score = 0.631371
	Best score = 0.633531
	init = 109 branch = -13 unpaired = 0
Epoch #8: 
	Average F-Score = 0.633684
	Best score = 0.634308
	init = 109 branch = -13 unpaired = 0
Best parameters: init = 109 branch = -13 unpaired = 0
```

Each epoch is a training set. We can see the F-Score steadily increases until we converge on a parameter set, which is given at the end. Note that the parameters found are different from those in the paper, since we training on the small data rather than the large.

The same training process will work for any of the parameter optimization programs. Please keep in mind that parameter training uses all cores and can therefore use a lot of memory!
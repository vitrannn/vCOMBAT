@Written by Vi Tran on June 2020
All the commands below are Unix commands and can be performed with Mac Terminal or Unix Shell or Cygwin in Windows.

Before running the code, the system needs to have a list of software properly installed. The list of software and how to install them are described in "installation_requirement.docx" document.

1.	How to compile (build the executable files):

•	Get to the code folder

$ ls
C_Code
bin
compileCode.sh
inputsample.txt
output_from_C
output_from_R

•	Build and Compile the code
$ ./compileCode.sh
Expected to see:
	changing directory..
	erasing stuff..
	cmaking..
	-- The C compiler identification is GNU 8.3.0
	-- The CXX compiler identification is GNU 8.3.0

Make sure the C and CXX compiler is GNU version 8 or higher. 

If you are using Mac and the working C compiler is the default clang gcc by MacOS, we need to link the working gcc with the newly installed GNU gcc. We do so by changing CC and CXX to the installed gcc location. For example, link gcc to the installed gcc version 8  at /usr/local/bin/gcc-8 and /usr/local/bin/g++-8 as in the following commands. 

$ export CC=/usr/local/bin/gcc-8
$ export CXX=/usr/local/bin/g++-8
$ export PATH=/usr/bin:/usr/local/bin:$PATH

2.	How to run
Run executable file tuberculosis_simulation from the C code folder 
$ cd ..
$./bin/tuberculosis_simulation -d 1000000 -p 100000  -t 100:1 -n 100 -R 0  -K 0  -A 1  -D 0.01  -V 0.000000000000001  -C 1000000 -k 50 -r 50 -m output3_100_t100_sp1.csv -i A_input3_100_t100_sp1.csv -v


Each parameter called by the C code is explained in ModelParameterinCode.docx document. The simulation takes the external concentrations from A_input3_100_t100_sp1.csv. The mathematical output results are stored in output3_100_t100_sp1.csv in the same folder.

3.	Expected results running with above parameters in command line console

Starting simulation with arguments
------------------------------------

Starting population     	100000
Antibiotic dose         	1e+06
Time of simulation      	100
Step size               	1

Target molecules        	100
Maximum kill rate       	0
Killing threshold       	50
Baseline replication    	0
Target association rate 	1
Target dissociation rate	0.01
Drug Molecular Weight    	555.5
Carrying capacity       	1e+06
Intracellular volume    	1e-15
Outputing compartmentBoundComplexState matrix to simCode/output3_100_t100_sp1.csv
Output File order is as follows:
L0 Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Li Ln tm BP An AT 

creating system with 104 free variables
Results readout
---------------

Final population 100000

It took me (9.858000 milliseconds).

=============================

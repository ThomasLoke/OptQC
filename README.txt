README for OptQC:

- This program is designed to use simulated annealing to find an optimal permutations P and Q such that
  the decomposition of the unitary matrix U = Q^T.P^T.U'.P.Q requires less gates.

- Requires a multiprocessor system running MPI and LAPACK.

- The makefile is provided in "Makefile" - to compile, simply run "make". The makefile assumes
  that the modules for MPI, Intel MKL libraries (for the LAPACK subroutines) and the Intel Fortran
  compiler have been preloaded to the system.

- The command to execute the program (inside a PBS script) follows the format:

	mpirun ./OptQC <Name-of-matrix-file> <Real/Complex-matrix> <Number-of-iterations-for-SA> <Number-of-iterations-for-QS> <Alpha-value>

  where five parameters are supplied to the program. The parameters are:
-> <Name-of-matrix-file>: The file in which the matrix is stored. For example, if the matrix was
   stored in a file "FileA.txt", then the argument would be supplied as "FileA" (without the
   quotation marks). It is assumed that all matrices will be stored in a file with a .txt extension.
   The file containing the matrix follows the following format:

   	<Dimension of matrix U>
   	U(1,1)
   	U(1,2)
   	...
   	U(1,n)
   	U(2,1)
   	...
   	U(n,n)

   where n is the dimension of the matrix U.
   For a real matrix, the element U(i,j) is just a real number in decimal form.
   For a complex matrix, the elements U(i,j) are written as:
	
	<Real-part>, <Imaginary-part>

   where both the real and imaginary parts are real numbers in decimal form.
-> <Real/Complex-matrix>: This argument indicates whether the matrix is to be interpreted as
   a real or complex matrix.
   For a real matrix, this argument is 0.
   For a complex matrix, this argument is 1.
-> <Number-of-iterations-for-SA>: This argument specifies the number of iterations to use when performing
   simulated annealing. Hence, a positive integer must be supplied.
-> <Number-of-iterations-for-QS>: This argument specifies the number of iterations to use when performing
   the search for a different qubit permutation. Hence, a positive integer must be supplied.
-> <Alpha-value>: This argument specifies the alpha value, which determines the acceptance t
   threshold in the simulated annealing process. Must be a floating point value in the range 0 <= alpha < 1
- Several output files are generated from this program (* indicates the <Name-of-matrix-file> argument):
-> "*_gates.txt": This file gives the optimal decomposition of the matrix in text.
-> "*_plot.tex": This file, together with the Qcircuit TeX package
   (see http://www.cquic.org/Qcircuit/Qcircuit.tex), can be used to draw the quantum circuit of the
   optimal decomposition using any TeX installation.
-> "*_history.dat": This file provides the time-series of the simulated annealing process for the process
   that reaches the optimal decomposition.
-> "*_perm.dat": This file gives the chosen qubit permutation q in list form, and the chosen permutation p in
   list form.

- On the Fornax supercomputer (access provided by iVEC@UWA), the PBS script (job.pbs) used to call mpirun 
  is as follows:

#!/bin/tcsh
#PBS -l nodes=8:ppn=12,mem=64gb
#PBS -l walltime=01:00:00
#PBS -W group_list=xxx
#PBS -q workq
#PBS -j oe

cd $PBS_O_WORKDIR
mpirun ./OptQC RandReal 0 40000 1000 0.010

  This PBS script executes the program on 8 nodes (96 threads), reading a complex unitary matrix from the file RandReal.txt, 
  using 40000 iterations for simulated annealing and 1000 iterations to find a different qubit permutation. The output files 
  "RandReal_plot.txt", "RandReal_plot.tex", "RandReal_history.dat" and "RandReal_perm.dat" are produced. See the file 
  "job.pbs.o1688661" for the console output produced by the program.
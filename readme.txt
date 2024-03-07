Matrix Transpose Program

Description:
Our program performs the transposition of a square matrix using parallel processing with MPI. It offers three different methods for the all-to-all communication primitive: MPI_Alltoall, hypercubic permutation-based all-to-all, and arbitrary permutation-based all-to-all.

Custom all-to-all communications:
The hypercubic permutation-based method involves performing a series of exchanges between pairs of processes whose ranks differ in only one bit position. Each process starts with a segment of the original matrix and exchanges data with its partners in each stage of the hypercube, gradually assembling the transposed matrix.
The arbitrary permutation-based method makes each process send and receive data from any other process, not just those that differ by one bit position. Each process sends a segment of its local matrix to every other process and receives the corresponding segments to assemble the transposed matrix.

Compilation:
A Makefile is provided for easy compilation. Simply run the command `make` in the directory containing the source code and Makefile. This will generate the executable named `transpose`.

Execution:
To run the program, use the following command:
mpirun -np <num_procs> ./transpose <input_file> <output_file> <all_to_all_choice> <matrix_size>
- <num_procs>: Number of processes to use (should be a power of 2)
- <input_file>: Path to the file containing the input matrix.
- <output_file>: Path to the file where the transposed matrix will be written.
- <all_to_all_choice>: Choice of all-to-all communication method ('a' for arbitrary, 'h' for hypercubic, 'm' for MPI_Alltoall)
- <matrix_size>: Size of the square matrix (n x n)

Input File Format:
The input file should contain the square matrix with rows separated by new lines and elements within each row separated by spaces.

Output File Format:
The output file will contain the transposed matrix in the same format as the input file.

Example:
mpirun -np 8 ./transpose matrix.txt transpose.txt a 24
This will transpose a 24x24 matrix using the arbitrary permutation-based all-to-all communication method with 8 processes.

Notes:
- The number of processes (<num_procs>) must be a power of 2 and should divide the matrix size (<matrix_size>)
- Ensure that the input matrix file exists and is properly formatted before running the program.

Machine Used:
- Model: MacBook Pro
- Size: 14-inch, 2021
- Chip: Apple M1 Pro
- Memory: 16 GB
- Operating System: macOS Sonoma 14.1.1
- Filesystem Capacity: 460 Gi with 43% used
- MPI Version: Open MPI 5.0.1
- Compiler: Apple clang version 15.0.0 (clang-1500.0.40.1), arm64-apple-darwin23.1.0

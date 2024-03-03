#include <iostream>
#include <mpi.h>
#include <fstream>

// mpirun -np 8 ./transpose matrix.txt transpose.txt a 24

int HPC_Alltoall_H(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
    // hypercubic permutation
    return MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
}

int HPC_Alltoall_A(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
    // arbitrary permutation
    return MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
}

int alltoall(char choice, const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
    switch (choice)
    {
    case 'm':
        return MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);    
    case 'a':
        return HPC_Alltoall_A(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);    
    case 'h':
        return HPC_Alltoall_H(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);    
    default:
        std::cerr << "Enter a valid All to All choice: 'a', 'h', 'm'" << std::endl;
        return -1;
    }
}

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (!world_rank)
        if (argc != 5) {
            std::cerr << "Usage: " << argv[0] << " <input.txt> <output.txt> <permutation_choice> <matrix_dimension>" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

    const char* input_file = argv[1];
    const char* output_file = argv[2];
    char perm_choice = argv[3][0];
    int matrix_dim = std::stoi(argv[4]);

    int rows_per_procs = matrix_dim / world_size;

    std::vector<int> matrix(matrix_dim * matrix_dim);
    std::vector<int> local_matrix(rows_per_procs * matrix_dim);
    std::vector<int> transp(matrix_dim * matrix_dim);
    std::vector<int> local_transp(rows_per_procs* matrix_dim);

    std::fstream infile(input_file);
    if (!world_rank)
        for (int i = 0; i < matrix_dim; i++)
            for (int j = 0; j < matrix_dim; j++) infile >> matrix[i * matrix_dim + j];

    MPI_Scatter(matrix.data(), rows_per_procs * matrix_dim, MPI_INT, local_matrix.data(), rows_per_procs * matrix_dim, MPI_INT, 0, MPI_COMM_WORLD);

    // Used this to check and it works
    // if (world_rank == 1) {
    //     for (int k = 0; k < matrix_dim * rows_per_procs; k++) {
    //         std::cout << local_matrix[k] << " ";
    //     }
    //     std::cout << world_rank << std::endl;
    // }

    double start_time = MPI_Wtime();

    // we need to modify this so that it sends by column not row
    alltoall(perm_choice, local_matrix.data(), rows_per_procs, MPI_INT, local_transp.data(), rows_per_procs, MPI_INT, MPI_COMM_WORLD);

    //
    if (world_rank == 0) {
        for (int k = 0; k < matrix_dim; k++) {
            std::cout << local_transp[k] << " ";
        }
    }
    
    double end_time = MPI_Wtime();
    double time_taken_ms = (end_time - start_time) * 1000.0; 

    MPI_Gather(local_transp.data(), rows_per_procs * matrix_dim, MPI_INT, transp.data(), rows_per_procs * matrix_dim, MPI_INT, 0, MPI_COMM_WORLD);

    std::ofstream outfile(output_file);
    if (!world_rank) {
        std::ofstream outfile(output_file);
        for (int i = 0; i < matrix_dim; i++) {
            for (int j = 0; j < matrix_dim; j++) {
                outfile << transp[j * matrix_dim + i] << " ";
            }
            outfile << std::endl;
        }
        printf("Time taken: %.6f milliseconds\n", time_taken_ms);
    }
    MPI_Finalize();
    
    return 0;
}

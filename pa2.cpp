#include <iostream>
#include <mpi.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <unordered_map>

int MPI_Type_size_wrapper(MPI_Datatype datatype) {
    int size;
    MPI_Type_size(datatype, &size);
    return size;
}

int log2_int(int value) {
    int logbase2 = 0;
    while (value >>= 1) ++logbase2;
    return logbase2;
}

int HPC_Alltoall_H(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
    int world_size, world_rank;
    MPI_Comm_rank(comm, &world_rank);
    MPI_Comm_size(comm, &world_size);

    // number of stages of the hypercube (stages)
    int stages = log2_int(world_size);

    std::vector<MPI_Request> send_requests(stages);
    std::vector<MPI_Request> recv_requests(stages);
    std::vector<MPI_Status> send_statuses(stages);
    std::vector<MPI_Status> recv_statuses(stages);

    int modified_sendcount;
    int modified_recvcount;
    int splits;
    int partner;
    int sendbuf_size;

    splits = (world_size / (1 << 1));
    modified_sendcount = sendcount * splits;
    modified_recvcount = recvcount * splits;
    sendbuf_size = sendcount * world_size * MPI_Type_size_wrapper(sendtype);

    char* temp_send_buffer = new char[sendbuf_size];
    memcpy(temp_send_buffer, static_cast<char*>(const_cast<void*>(sendbuf)), sendbuf_size);

    for (int j = 0; j < stages; j++) {

        partner = world_rank ^ (1 << (stages - j - 1));
        char* send_pos = (world_rank & (1 << (stages - j - 1))) ? 
                         temp_send_buffer :
                         temp_send_buffer + modified_sendcount * MPI_Type_size_wrapper(sendtype);

        char* recv_pos = static_cast<char*>(recvbuf);

        MPI_Isend(send_pos, modified_sendcount, sendtype, partner, 0, comm, &send_requests[j]);
        MPI_Irecv(recv_pos, modified_recvcount, recvtype, partner, 0, comm, &recv_requests[j]);

        MPI_Wait(&send_requests[j], &send_statuses[j]);
        MPI_Wait(&recv_requests[j], &recv_statuses[j]);

        char* storage_buffer = new char[sendbuf_size];
        for (int k = 0; k < splits; k++) {
            char* send_pos = (world_rank & (1 << (stages - j - 1))) ? 
                            temp_send_buffer + modified_sendcount * MPI_Type_size_wrapper(sendtype) :
                            temp_send_buffer;
    
            memcpy(storage_buffer + sendcount * MPI_Type_size_wrapper(sendtype) * k * 2, send_pos + sendcount * MPI_Type_size_wrapper(sendtype) * k, sendcount * MPI_Type_size_wrapper(sendtype));
            memcpy(storage_buffer + sendcount * MPI_Type_size_wrapper(sendtype) * (2 * k + 1), recv_pos + sendcount * MPI_Type_size_wrapper(sendtype) * k, sendcount * MPI_Type_size_wrapper(sendtype));
        }
        memcpy(temp_send_buffer, storage_buffer, sendbuf_size);
        delete[] storage_buffer;
    }

    int ref_world_rank = world_rank;

    char* reversed = new char[sendbuf_size];
    char* holder = new char[sendbuf_size];
    int packet_size = sendcount * MPI_Type_size_wrapper(sendtype);
    int tuple_packet_size = 2 * packet_size;

    if (world_rank % 2 != 0) {

        ref_world_rank = world_size - 1 - world_rank;

        for (int k = 0; k < sendbuf_size; k+=packet_size) {
            memcpy(reversed + k, temp_send_buffer + (sendbuf_size - k - packet_size), packet_size);
        }
        memcpy(temp_send_buffer, reversed, sendbuf_size);
    }

    int dist = tuple_packet_size * (ref_world_rank / 2); // go to the copy position

    std::unordered_map<int, int> visited;

    for (int k = 0, i = 0; k < sendbuf_size; k+=tuple_packet_size, i++) { // can be skipped for rank = 0

        if (visited.find(k) != visited.end() || visited.find((k + dist) % sendbuf_size) != visited.end()) {
            continue;
        }

        visited[k] = 1;
        visited[((k + dist) % sendbuf_size)] = 1;
        memcpy(holder + k, temp_send_buffer + ((k + dist) % sendbuf_size), tuple_packet_size);
        memcpy(holder + ((k + dist) % sendbuf_size), temp_send_buffer + k, tuple_packet_size);
    }
    memcpy(static_cast<char*>(recvbuf), holder, sendbuf_size);
    delete[] temp_send_buffer;
    delete[] reversed;
    delete[] holder;
    return MPI_SUCCESS;
}

int HPC_Alltoall_A(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {

    int world_size, world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    std::vector<MPI_Request> send_requests(world_size - 1);
    std::vector<MPI_Request> recv_requests(world_size - 1);

    std::vector<MPI_Status> send_statuses(world_size - 1);
    std::vector<MPI_Status> recv_statuses(world_size - 1);

    // copy the data that is supposed to be sent to itself directly into the recvbuf
    memcpy(static_cast<char*>(recvbuf) + world_rank * recvcount * MPI_Type_size_wrapper(recvtype), 
                static_cast<const char*>(sendbuf) + world_rank * sendcount * MPI_Type_size_wrapper(sendtype), 
                sendcount * MPI_Type_size_wrapper(sendtype));

    for (int j = 1; j < world_size; j++) {

        int send_partner = (world_rank + j) % world_size;
        int recv_partner = (world_rank - j + world_size) % world_size;

        MPI_Isend(static_cast<const char*>(sendbuf) + send_partner * sendcount * MPI_Type_size_wrapper(sendtype), sendcount, sendtype, send_partner, 0, comm, &send_requests[j - 1]);

        MPI_Irecv(static_cast<char*>(recvbuf) + recv_partner * recvcount * MPI_Type_size_wrapper(recvtype), recvcount, recvtype, recv_partner, 0, comm, &recv_requests[j - 1]);
    }

    // wait for all non-blocking operations to complete
    MPI_Waitall(world_size - 1, &send_requests[0], &send_statuses[0]);
    MPI_Waitall(world_size - 1, &recv_requests[0], &recv_statuses[0]);

    return MPI_SUCCESS;
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
    std::vector<int> temp_transp(rows_per_procs * matrix_dim);

    std::fstream infile(input_file);
    if (!world_rank)
        for (int i = 0; i < matrix_dim; i++)
            for (int j = 0; j < matrix_dim; j++) infile >> matrix[i * matrix_dim + j];

    MPI_Scatter(matrix.data(), rows_per_procs * matrix_dim, MPI_INT, local_matrix.data(), rows_per_procs * matrix_dim, MPI_INT, 0, MPI_COMM_WORLD);

    double start_time = MPI_Wtime();

    // transpose done in each processor to prepare for the all to all
    for (int i = 0; i < matrix_dim; ++i) {
        for (int j = 0; j < rows_per_procs; ++j) {
            temp_transp[i * rows_per_procs + j] = local_matrix[j * matrix_dim + i];
        }
    }

    int items_per_procs = (rows_per_procs * matrix_dim) / world_size;   

    alltoall(perm_choice, temp_transp.data(), items_per_procs, MPI_INT, local_transp.data(), items_per_procs, MPI_INT, MPI_COMM_WORLD);

    // reorganization needs to happen here to prepare for gathering
    std::vector<int> gatherready(rows_per_procs * matrix_dim);
    int count = 0;
    int shift = 0;
    for (int i = 0; i < matrix_dim; i++) {
        for (int j = 0; j < rows_per_procs; j++) {
            gatherready[i * rows_per_procs + j] = local_transp[(((i + count) * rows_per_procs + shift) % (matrix_dim * rows_per_procs)) + j];
        }

        count += rows_per_procs - 1;
        if ((i + 1) % (matrix_dim/rows_per_procs) == 0 && i) {shift += rows_per_procs;}
    }

    double end_time = MPI_Wtime();
    double time_taken_ms = (end_time - start_time) * 1000.0; 

    MPI_Gather(gatherready.data(), rows_per_procs * matrix_dim, MPI_INT, transp.data(), rows_per_procs * matrix_dim, MPI_INT, 0, MPI_COMM_WORLD);

    std::ofstream outfile(output_file);
    if (!world_rank) {
        std::ofstream outfile(output_file);
        for (int i = 0; i < matrix_dim; i++) {
            for (int j = 0; j < matrix_dim; j++) {
                outfile << transp[i * matrix_dim + j] << " ";
            }
            outfile << std::endl;
        }
        printf("Time taken: %.6f milliseconds,\n", time_taken_ms);
    }

    MPI_Finalize();
    return 0;
}

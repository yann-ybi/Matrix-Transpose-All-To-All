#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Check for proper command line argument
    if (argc < 2 && world_rank == 0) {
        std::cerr << "Usage: " << argv[0] << " <number_of_points>" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Total number of points to generate for estimating π
    long long int total_points = atoll(argv[1]);
    long long int points_per_process = total_points / world_size;
    long long int points_inside_circle = 0;

    // Seed the random number generator to vary results across different processes
    std::srand(static_cast<unsigned int>(time(nullptr)) + world_rank * 123);

    // Generate points and count how many fall inside the unit circle
    for (long long int i = 0; i < points_per_process; ++i) {
        double x = 2.0 * rand() / RAND_MAX - 1.0;
        double y = 2.0 * rand() / RAND_MAX - 1.0;
        if (x*x + y*y <= 1.0) {
            ++points_inside_circle;
        }
    }

    // Gather the counts of points inside the circle from all processes
    long long int total_points_inside_circle = 0;
    MPI_Reduce(&points_inside_circle, &total_points_inside_circle, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Master process computes and prints the estimated value of π
    if (world_rank == 0) {
        double pi_estimate = 4.0 * total_points_inside_circle / total_points;
        std::cout << "Estimated π = " << pi_estimate << std::endl;
    }

    MPI_Finalize();
    return 0;
}

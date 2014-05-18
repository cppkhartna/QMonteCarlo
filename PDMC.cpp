#include "PDMC.h"

using namespace std;
const double eh_to_ev = 27.211384523232323;

PDMC::PDMC()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
        start = MPI_Wtime();

    long utime, seed;
    utime = (long) time(NULL);

    seed = abs(((utime*181)*((rank-83)*359))%104729); 

    mt_seed32(seed);   MPI_Comm_size(MPI_COMM_WORLD, &size);
};

void PDMC::run(QModel* Q, int tau_max)
{
    double E_curr = 0, E_cum = 0, E_avg = 0;
    //if (E_array == NULL)
        //E_array = new double[tau_max+1];
    int N_0;
    for (int tau = 0; tau <= tau_max; tau++)
    {
        E_curr = (Q->getE_r() + Q->E_proton())*eh_to_ev;
        if (tau > 3000)
        {
            E_cum += E_curr;
            E_avg = E_cum/(tau-3000);
        }
        else
        {
            E_avg = E_curr;
        }
        //E_array[tau] = E_avg;

        if (tau % 100 == 0 && tau != 0)
        {
            cout.precision(17);
            cout << tau << " : " << N_0 << " : " << E_avg << endl;
        }

        Q->setV_avg(0.0);

        Q->walk();
        Q->branch();

        unsigned long long N_1 = Q->getN_1();
        double V_avg = Q->getV_avg();
        double E_r;

        double vals[size];
        unsigned long long nums[size];

        MPI_Gather(&V_avg, 1, MPI_DOUBLE, vals, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&N_1, 1, MPI_UNSIGNED_LONG_LONG, nums, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            for (int j = 0; j < size; j++)
            {
                vals[0] += vals[j];
                nums[0] += nums[j];
            }
            V_avg = vals[0];
            E_r = V_avg/double(nums[0]);
        }
        MPI_Bcast(&E_r, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        Q->setE_r(E_r);

        N_0 = N_1;
    }
}

//replica* QModel::getReplicas()
//{
    //return replicas;
//};

//double* QModel::getEnergies()
//{
    //return E_array;
//};

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    QMolecule *mol = new QMolecule();
    PDMC P;
    int N_0 = 4000;
    mol->add("H", 1.0, vec_3d(0, 0, 1.0));
    mol->add("H", 1.0, vec_3d(0, 0, -1.0));
    mol->setD(2);
    double x[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    mol->init_replicas(N_0, x);

    P.run(mol, 4000);

    MPI_Finalize();
}

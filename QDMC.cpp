#include "QDMC.h"

using namespace std;
const double eh_to_ev = 27.211384523232323;

QDMC::QDMC()
{
    long utime;
    utime = (long) time(NULL);
    mt_seed32(utime);
};

void QDMC::run(QModel* Q, int tau_max)
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

        int N_1 = Q->getN_1();
        double E_r = Q->getV_avg()/double(N_1);
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

int main()
{
    QMolecule *mol = new QMolecule();
    QDMC Q;
    int N_0 = 4000;
    mol->add("H", 1.0, vec_3d(0, 0, 1.0));
    mol->add("H", 1.0, vec_3d(0, 0, -1.0));
    //mol->add("H", 1.0, vec_3d(0, 0, -1.0));
    //mol->add("O", 1.0, vec_3d(0, 1, 0.0));
    mol->setD(2);
    //double x[3] = {0.0, 0.0, 0.0};
    double x[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    //std::cout << "Joba" << std::endl;
    //std::cout << mol->E_proton() << std::endl;
    mol->init_replicas(N_0, x);
    //std::cout << "Joba" << std::endl;
    Q.run(mol, 4000);
    //run(4000, 10000, 3, 0.0, NULL);
}

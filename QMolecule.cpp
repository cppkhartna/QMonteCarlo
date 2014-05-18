#include "QMolecule.h"

using namespace std;
const double eh_to_ev = 27.211384523232323;

QDMC::QDMC()
{
    long utime;
    utime = (long) time(NULL);
    mt_seed32(utime);
};

QModel::~QModel()
{
    replicas.clear();
};

void QModel::init_replicas(int N_0, double *x)
{
    this->N_0 = N_0;
    for (int i = 0; i < N_0; i++)
    {
        replica *current = new replica();

        for (int j = 0; j < d; j++)
        {
            current->coords.resize(d);
            for (int k = 0; k < MAXDIM; k++)
                current->coords[j].add_random(k, x[j], x_max, x_min);
        }

        replicas.push_front(*current);
    }
};

void QModel::walk()
{
    for (auto it = replicas.begin(); it != replicas.end(); it++)
    {
        for (int j = 0; j < d; j++)
            for (int k = 0; k < MAXDIM; k++)
            {
                double p = rd_normal(0.0, 1.0);
                double change = sqrt(dtau) * p;
                it->coords[j].add_random(k, change, x_max, x_min);
            }
    }
}

void QModel::branch()
{
    for (auto it = replicas.begin(); it != replicas.end(); ++it)
    {
        double u = rd_uniform(0.0, 1.0);
        int m_n = min(int(W(*it) + u), 3);
        if (m_n == 0)
        {
            replicas.erase(it++);
        }
        else if (m_n != 1)
        {
            if (N_0 < N_max)
            {
                for (int j = 1; j < m_n; j++)
                {
                    //copy
                    replica *current = new replica();

                    current->coords.resize(d);
                    for (int i = 0; i < d; i++)
                    {
                        current->coords[i] = it->coords[i];
                    }
                    replicas.push_front(*current);
                }
            }
        }
    }
    N_1 = replicas.size();
}

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

inline double QModel::W(replica &x)
{
    double V_x = V(x);
    V_avg += V_x;
    return exp(-(V_x - E_r)*dtau);
}

//replica* QModel::getReplicas()
//{
    //return replicas;
//};

//double* QModel::getEnergies()
//{
    //return E_array;
//};

inline double QMolecule::E_proton()
{
    double result = 0;
    for (auto it = nuclei.begin(); it < nuclei.end(); it++)
    {
        for (auto jt = it+1; jt < nuclei.end(); jt++)
            result += ((*jt).charge*(*it).charge/sub_len((*jt).pos, (*it).pos));
    }
    return result;
};

void QMolecule::add(string name, double charge, vec_3d pos)
{
   nucleus *nuc = new nucleus(name, charge, pos);
   nuclei.push_back(*nuc);
};

nucleus::nucleus(string name, double charge, vec_3d pos)
{
    this->name = name;
    this->charge = charge;
    this->pos = pos;
};

inline double QMolecule::V(replica &rep)
{
    double V_x = 0.0;
    for (int i = 0; i < d; i++)
    {
        if (d == 1 && nuclei.size() == 1) // atom
        {
            V_x -= (1/sqrt(rep.coords[i].length2()));
        }
        else 
        {
            // ion+
            for (auto it = nuclei.begin(); it < nuclei.end(); it++)
            {
                V_x -= ((*it).charge/sub_len(rep.coords[i], (*it).pos));
            }
            // molecule
            for (int j = i+1; j < d; j++)
            {
                //pos1 = rep.coords[j];
                V_x += (1/sub_len(rep.coords[i], rep.coords[j]));
            }
        }
    
    }
    return V_x;
};

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

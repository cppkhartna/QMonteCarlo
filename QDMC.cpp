#include "QDMC.h"

using namespace std;
QDMC Qh;

QDMC::QDMC()
{
    long utime;
    utime = (long) time(NULL);
    mt_seed32(utime);
};

QDMC::~QDMC()
{
    replica *p1 = replicas;
    replica *p0 = NULL;
    while (p1 != NULL)
    {
        delete[] p1->x;
        p0 = p1->next;
        delete p1;
        p1 = p0;
    }
};

void QDMC::init_replicas(int N_0, double x, double y, double z)
{
    this->N_0 = N_0;
    for (int i = 0; i < N_0; i++)
    {
        replica *current = new replica();

        current->x = new double[3];
        current->x[0] = x;
        current->x[1] = y;
        current->x[2] = z;

        current->next = replicas;
        replicas = current;
    }
};

void QDMC::walk()
{
    replica *p1 = replicas;
    while (p1 != NULL)
    {
        for (int j = 0; j < d; j++)
        {
            double p = rd_normal(0.0, 1.0);
            double change = sqrt(dtau) * p;
            p1->x[j] += change;
            if (p1->x[j] > x_max)
                p1->x[j] = x_max;
            if (p1->x[j] < x_min)
                p1->x[j] = x_min;
        }
        p1 = p1->next;
    }
}

void QDMC::branch()
{
    replica *p1 = replicas;
    replica *p0 = NULL;
    N_1 = 0;
    while (p1 != NULL)
    {
        N_1++;
        double u = rd_uniform(0.0, 1.0);
        int m_n = min(int(W(p1) + u), 3);
        if (m_n == 0)
        {
            if (p0 != NULL)
            {
                p0->next = p1->next;
            }
            else
            {
                p0 = p1->next;
            }

            delete[] p1->x;
            delete p1;
            p1 = p0;
        }
        else if (m_n != 1)
        {
            for (int j = 1; j < m_n; j++)
            {
                replica *current = new replica();

                current->x = new double[3];
                for (int i = 0; i < d; i++)
                    current->x[i] = p1->x[i];

                current->next = replicas;
                replicas = current;
            }
        }
        p0 = p1;
        p1 = p1->next;
    }
}

void QDMC::run(int N_0, int tau_max)
{
    //init_replicas(N_0);
    for (int tau = 0; tau <= tau_max; tau++)
    {
        if (tau % 100 == 0)
        {
            cout.precision(17);
            cout << tau << " : " << N_0 << " : " <<  E_r << endl;
        }
        V_avg = 0.0;

        walk();
        branch();

        E_r = V_avg/double(N_0);
        N_0 = N_1;
    }
}

inline double QDMC::W(replica *x)
{
    double V_x = model->V(x);
    V_avg += V_x;
    return exp(-(V_x - E_r)*dtau);
}

replica* QDMC::getReplicas()
{
    return replicas;
};

extern "C" replica* run(int N_0, int tau_max, bool is_atom)
{
    if (is_atom)
    {
        Qh.init_replicas(N_0, 0.0, 0.0, 1.0);
        QAtomH* atom = new QAtomH();
        Qh.setModel(atom);
        std::cout << "atom" << std::endl;
    }
    else
    {
        Qh.init_replicas(N_0, 0.0, 0.0, 0.0);
        QIonH* ion = new QIonH();
        ion->setR(2.0);
        Qh.setModel(ion);
        std::cout << "ion" << std::endl;
    }
    Qh.run(N_0, tau_max);
    return Qh.getReplicas();
};

inline double QAtomH::V(replica *rep)
{
    double V_x = 0.0;
    double x = rep->x[0];
    double y = rep->x[1];
    double z = rep->x[2];
    double r = sqrt(x*x + y*y + z*z);
    if (r > 0.0)
    {
        V_x = V_x - 1.0/r;
    }
    return V_x;
};

void QDMC::setModel(QModel* mod)
{
    model = mod;
};

void QIonH::setR(double R_proton)
{
    R = R_proton;
};

inline double QIonH::V(replica *rep)
{
    double V_x = 0.0;
    double x = rep->x[0];
    double y = rep->x[1];
    double z = rep->x[2];
    double r1 = sqrt(x*x + y*y + (z - 0.5*R)*(z - 0.5*R));
    double r2 = sqrt(x*x + y*y + (z + 0.5*R)*(z + 0.5*R));
    if (r1 > 0.0 && r2 > 0)
    {
       V_x  = V_x - 1.0/r1 - 1.0/r2;
    }
    return V_x;
};

int main()
{
    run(4000, 10000, false);
}

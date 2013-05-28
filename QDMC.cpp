#include "QDMC.h"

using namespace std;
const double eh_to_ev = 27.211384523232323;

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

void QDMC::init_replicas(int N_0, double *x)
{
    this->N_0 = N_0;
    for (int i = 0; i < N_0; i++)
    {
        replica *current = new replica();

        current->x = new double[d];
        for (int j = 0; j < d; j++)
            current->x[j] = x[j];

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
            //if (N_0 < N_max)
            {
                for (int j = 1; j < m_n; j++)
                {
                    replica *current = new replica();

                    current->x = new double[d];
                    for (int i = 0; i < d; i++)
                        current->x[i] = p1->x[i];

                    current->next = replicas;
                    replicas = current;
                }
            }
        }
        p0 = p1;
        p1 = p1->next;
    }
}

void QDMC::run(int tau_max)
{
    double E_curr = 0, E_cum = 0, E_avg = 0;
    if (E_array == NULL)
        E_array = new double[tau_max+1];
    for (int tau = 0; tau <= tau_max; tau++)
    {
        E_curr = (E_r + E_proton())*eh_to_ev;
        if (tau > 3000)
        {
            E_cum += E_curr;
            E_avg = E_cum/(tau-3000);
        }
        else
        {
            E_avg = E_curr;
        }
        E_array[tau] = E_avg;

        if (tau % 100 == 0 && tau != 0)
        {
            cout.precision(17);
            cout << tau << " : " << N_0 << " : " << E_avg << endl;
        }

        V_avg = 0.0;

        walk();
        branch();

        E_r = V_avg/double(N_1);
        N_0 = N_1;
    }
}

inline double QDMC::W(replica *x)
{
    double V_x = V(x);
    V_avg += V_x;
    return exp(-(V_x - E_r)*dtau);
}

replica* QDMC::getReplicas()
{
    return replicas;
};

double* QDMC::getEnergies()
{
    return E_array;
};

QDMC *Qh; 

extern "C" replica* run(int N_0, int tau_max, int is_atom, double R, double* es)
{
    if (Qh != NULL)
    {
        delete Qh;
    };

    if (is_atom == 1)
    {
        Qh = (QAtomH*) new QAtomH();
        double x[3] = {0.0, 0.0, 1.0};
        Qh->init_replicas(N_0, x);
    }
    else if (is_atom == 2)
    {
        Qh = (QIonH*) new QIonH();
        double x[3] = {0.0, 0.0, 1.0};
        Qh->init_replicas(N_0, x);
        if (R == 0)
            Qh->setR(2.0);
        else
            Qh->setR(R);
    }
    else
    {
        Qh = (QMoleculeH*) new QMoleculeH();
        double x[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        Qh->init_replicas(N_0, x);
        if (R == 0)
            Qh->setR(1.4);
        else
            Qh->setR(R);
    }
    Qh->setE_array(es);
    Qh->run(tau_max);
    return Qh->getReplicas();
};

void QIonH::setR(double R_proton)
{
    this->R_proton = R_proton;
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
        V_x = - 1.0/r;
    }
    return V_x;
};

inline double QIonH::V(replica *rep)
{
    double V_x = 0.0;
    double x = rep->x[0+ind];
    double y = rep->x[1+ind];
    double z = rep->x[2+ind];
    double r_minus = sqrt(x*x + y*y + (z - 0.5*R())*(z - 0.5*R()));
    double r_plus  = sqrt(x*x + y*y + (z + 0.5*R())*(z + 0.5*R()));
    if (r_minus > 0.0 && r_plus > 0)
    {
       V_x  = 0 - 1.0/r_plus - 1.0/r_minus;
    }
    return V_x;
};

inline double QMoleculeH::V(replica *rep)
{
    double V_x = 0.0;
    // first electron
    ind = 0;
    V_x += QIonH::V(rep);
    double x1 = rep->x[0];
    double y1 = rep->x[1];
    double z1 = rep->x[2];

    // second electron
    ind = 3;
    V_x += QIonH::V(rep);
    double x2 = rep->x[3];
    double y2 = rep->x[4];
    double z2 = rep->x[5];

    // r12
    double r12= sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    if (r12 > 0)
    {
       V_x  += 1.0/r12;
    }
    return V_x;
};

inline double QIonH::E_proton()
{
    return 1.0/(double)R();
};

inline double QAtomH::E_proton()
{
    return 0;
};

int main()
{
    run(4000, 10000, 3, 0.0, NULL);
}

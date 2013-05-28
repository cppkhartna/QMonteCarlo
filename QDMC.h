#include <time.h>
#include <iostream>
#include <cmath>
#include "mtwist/randistrs.h"

struct replica
{
public:
    double *x;
    replica *next;
};

class QDMC
{
    int N_0 = 4000;
    int N_1;
    int N_max = 10000;
    double dtau = 0.01;
    int x_min = -30;
    int x_max = 30;
    int n_b = 200;
    int d = 6; // 6
    double E_r = -1.0;
    double V_avg = 0.0;
    replica* replicas = NULL;
    double* E_array;
public:
    QDMC();
    virtual ~QDMC();
    void init_replicas(int N_0, double *x);
    void run(int tau_max);
    void walk();
    void branch();
    void count(){};
    replica* getReplicas();
    double* getEnergies();
    double W(replica *x);
    virtual double V(replica *rep) = 0;
    virtual double E_proton() = 0;
    virtual void setR(double R_proton){R_proton=0;};
};

class QAtomH: public QDMC
{
public:
    double V(replica *rep);
    double E_proton();
};

class QIonH: public QDMC
{
    double R_proton;
public:
    int ind = 0;
    double V(replica *rep);
    void setR(double R_proton);
    double R(){return R_proton;};
    double E_proton();
};

class QMoleculeH: public QIonH
{
public:
    double V(replica *rep);
};

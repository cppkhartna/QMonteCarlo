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

class QModel
{
public:
    virtual double V(replica *rep) = 0;
};

class QAtomH: public QModel
{
public:
    double V(replica *rep);
};

class QIonH: public QModel
{
    double R;
public:
    double V(replica *rep);
    void setR(double R_proton);
};

class QDMC
{
    int N_0 = 4000;
    int N_1;
    int N_max = 2000;
    double dtau = 0.01;
    int x_min = -30;
    int x_max = 30;
    int n_b = 200;
    int d = 3;
    double E_r = -1.0;
    double V_avg = 0.0;
    replica* replicas = NULL;
    QModel* model;
public:
    QDMC();
    ~QDMC();
    void init_replicas(int N_0, double x, double y, double z);
    void run(int N_0, int tau_max);
    void walk();
    void branch();
    void count(){};
    replica* getReplicas();
    double W(replica *x);
    void setModel(QModel* mod);
};

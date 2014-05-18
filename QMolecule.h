#include <time.h>
#include <iostream>
#include <cmath>
#include <list>
#include <vector>
#include <string>
#include "math.h"
#include "mtwist/randistrs.h"

using std::vector;
using std::string;
using std::list;

#define MAXDIM 3

struct replica
{
public:
    vector<vec_3d> coords;
};

class QModel
{
    int N_0 = 4000;
    int N_1;
    int N_max = 10000;
    double dtau = 0.01;
    double x_min = -30;
    double x_max = 30;
    int n_b = 200;
    double E_r = -1.0;
    double V_avg = 0.0;
    //double* E_array;
protected:
    int d = 2; // number of electrons
    list<replica> replicas;
public:
    QModel(){};
    virtual ~QModel();
    void init_replicas(int N_0, double *x);
    //void run(int tau_max);

    void walk();
    void branch();
    //void count(){};

    //replica* getReplicas();
    //double* getEnergies();

    double W(replica &x);
    virtual double V(replica &rep) = 0;
    virtual double E_proton() = 0;

    void setD(int dd){d = dd;};
    void setE_r(double _E_r){E_r = _E_r;};
    void setV_avg(double _V_avg){V_avg = _V_avg;};

    double getE_r(){return E_r;};
    double getV_avg(){return V_avg;};
    long getN_1(){return N_1;};
    //void setE_array(double* es){E_array = es;};
};

class QDMC
{
public:
    QDMC();
    void run(QModel* Q, int tau_max);
};

struct nucleus
{
public:
    string name;
    double charge;
    vec_3d pos;
    nucleus(string name, double charge, vec_3d pos);
    ~nucleus(){};
    friend std::ostream& operator<<(std::ostream& os, const nucleus& nl);
};

//: public QModel
class QMolecule : public QModel
{
    vector<nucleus> nuclei;
public:
    double V(replica &rep);
    void add(string name, double charge, vec_3d pos);
    double E_proton();
};

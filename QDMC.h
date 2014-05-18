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

class QDMC
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
    double* E_array;
protected:
    int d = 2; // number of electrons
    list<replica> replicas;
public:
    QDMC();
    virtual ~QDMC();
    void init_replicas(int N_0, double *x);
    void run(int tau_max);

    void walk();
    void branch();

    void count(){};
    //replica* getReplicas();
    double* getEnergies();
    double W(replica &x);
    virtual double V(replica &rep) = 0;
    virtual double E_proton() = 0;
    virtual void setR(double R_proton){R_proton=0;};
    void setD(int dd){d=dd;}
    void setE_array(double* es){E_array = es;};
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

//: public QDMC
class QMolecule : public QDMC
{
    vector<nucleus> nuclei;
public:
    double V(replica &rep);
    void add(string name, double charge, vec_3d pos);
    double E_proton();
};

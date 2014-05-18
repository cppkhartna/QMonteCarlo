#include <mpi.h>
#include "QMolecule.h"

class PDMC
{
    int rank, size;
    double end, start;
public:
    PDMC();
    void run(QModel* Q, int tau_max);
};


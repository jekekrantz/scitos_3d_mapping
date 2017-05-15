#ifndef DistanceWeightFunction2FILEtest_H
#define DistanceWeightFunction2FILEtest_H

#include <cmath>
#include <sys/time.h>
#include "DistanceWeightFunction2.h"

using namespace Eigen;
namespace reglib{


class DistanceWeightFunction2FILE : public DistanceWeightFunction2
{
public:
    std::string path;
    double maxd;
    double mind;
    unsigned int steps;
    double mul;

    double maxp;
    double minp;

    double * probs_vec;
    double * infront_vec;

    virtual DistanceWeightFunction2 * clone();

    DistanceWeightFunction2FILE(std::string path_);
    ~DistanceWeightFunction2FILE();

    virtual double getProbInfront(double d, bool debugg = false);
    virtual double getProb(double d, bool debugg = false);
    virtual std::string getString();
    virtual int getInd(double d);
};

}

#endif // DistanceWeightFunction2FILEtest_H

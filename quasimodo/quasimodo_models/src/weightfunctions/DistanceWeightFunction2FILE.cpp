#include "weightfunctions/DistanceWeightFunction2FILE.h"

namespace reglib{

DistanceWeightFunction2FILE::DistanceWeightFunction2FILE(std::string path_){
    path = path_;

    probs_vec   = 0;
    infront_vec = 0;

    maxp = 0.999;
    minp = 0.000;

    std::ifstream file (path, std::ios::in | std::ios::binary | std::ios::ate);
    if (file.is_open()){
      std::streampos size = file.tellg();
      char * buffer = new char [size];
      file.seekg (0, std::ios::beg);
      file.read (buffer, size);
      file.close();

      double* buffer_double         = (double *)buffer;
      unsigned long * buffer_long   = (unsigned long * )buffer;
      unsigned long current         = 0;
      steps = buffer_long   [current++];
      maxd  = buffer_double [current++];
      mind  = buffer_double [current++];

      printf("steps: %i maxd: %f mind: %f size: %i\n",steps,maxd,mind,size/sizeof(double));

      mul = 1.0/((maxd-mind)/double(steps-1));

      probs_vec   = new double[steps];
      infront_vec = new double[steps];

      for(unsigned int i = 0; i < steps; i++){
          probs_vec[i]   = buffer_double[current++];
      }
      for(unsigned int i = 0; i < steps; i++){infront_vec[i] = buffer_double[current++];}
      delete[] buffer;
    }else{
        std::cout << "Unable to open file";
    }
}
DistanceWeightFunction2FILE::~DistanceWeightFunction2FILE(){
    if(probs_vec != 0) {delete[] probs_vec;}
    if(infront_vec != 0) {delete[] infront_vec;}
}

inline int DistanceWeightFunction2FILE::getInd(double d){
    return (d-mind)*mul;
}

double  DistanceWeightFunction2FILE::getProb(double d,bool debugg){
    int ind = getInd(d);
    if(ind >= 0 && ind < steps){
        return std::min(maxp,std::max(minp,probs_vec[ind]));
    }else{
        return minp;
    }
}

double  DistanceWeightFunction2FILE::getProbInfront(double d,bool debugg){
    int ind = getInd(d);
    if(ind >= 0 && ind < steps){
        return std::min(maxp,std::max(minp,infront_vec[ind]));
    }else if(ind < 0){
        return minp;
    }else{
        return maxp;
    }
}

std::string DistanceWeightFunction2FILE::getString(){return "DistanceWeightFunction2FILE";}

DistanceWeightFunction2 * DistanceWeightFunction2FILE::clone(){
    DistanceWeightFunction2FILE * func  = new DistanceWeightFunction2FILE(path);
    return func;
}

}



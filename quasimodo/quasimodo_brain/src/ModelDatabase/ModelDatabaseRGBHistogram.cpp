#include "ModelDatabaseRGBHistogram.h"

ModelDatabaseRGBHistogram::ModelDatabaseRGBHistogram(int res_){
	res = res_;
	storage = new ModelStorageFile();
	printf("made a ModelDatabaseRGBHistogram(%i)\n",res);
}
ModelDatabaseRGBHistogram::~ModelDatabaseRGBHistogram(){}

std::vector< double > getDescriptor(int res, reglib::Model * model){
	std::vector< double > descriptor;
    int resi = res+1;
    descriptor.resize(resi*resi*resi);
    for(int i = 0; i < resi*resi*resi; i++){descriptor[i] = 0;}
	std::vector<reglib::superpoint> & points = model->points;
    for(unsigned int i = 0; i < points.size(); i++){
        reglib::superpoint & sp = points[i];
        int r = sp.r;
        int g = sp.g;
        int b = sp.b;
        int rind = int(double(res)*r/256.0);
        int gind = int(double(res)*g/256.0);
        int bind = int(double(res)*b/256.0);
        int ind = resi*resi*rind + resi*gind + bind;
        descriptor[ind]++;
    }
/*
    for(unsigned int i = 0; i < points.size(); i++){
        reglib::superpoint & sp = points[i];
        int r = sp.r;
        int g = sp.g;
        int b = sp.b;
        double rind = double(res)*r/256.0;
        double gind = double(res)*g/256.0;
        double bind = double(res)*b/256.0;

        double rind0 = int(rind);
        double gind0 = int(gind);
        double bind0 = int(bind);

        double rind1 = int(rind+1.0);
        double gind1 = int(gind+1.0);
        double bind1 = int(bind+1.0);

        double wr0   = 1.0- (rind - rind0);
        double wg0   = 1.0- (gind - gind0);
        double wb0   = 1.0- (bind - bind0);

        double wr1   = 1.0-wr0;
        double wg1   = 1.0-wg0;
        double wb1   = 1.0-wb0;

        descriptor[resi*resi*int(rind0) + resi*int(gind0) + int(bind0)]+=wr0*wg0*wb0;
        descriptor[resi*resi*int(rind0) + resi*int(gind0) + int(bind1)]+=wr0*wg0*wb1;
        descriptor[resi*resi*int(rind0) + resi*int(gind1) + int(bind0)]+=wr0*wg1*wb0;
        descriptor[resi*resi*int(rind0) + resi*int(gind1) + int(bind1)]+=wr0*wg1*wb1;
        descriptor[resi*resi*int(rind1) + resi*int(gind0) + int(bind0)]+=wr1*wg0*wb0;
        descriptor[resi*resi*int(rind1) + resi*int(gind0) + int(bind1)]+=wr1*wg0*wb1;
        descriptor[resi*resi*int(rind1) + resi*int(gind1) + int(bind0)]+=wr1*wg1*wb0;
        descriptor[resi*resi*int(rind1) + resi*int(gind1) + int(bind1)]+=wr1*wg1*wb1;
    }
    */
    for(int i = 0; i < resi*resi*resi; i++){descriptor[i] /= double(points.size());}
	return descriptor;
}

bool ModelDatabaseRGBHistogram::add(reglib::Model * model){  
    model->keyval = std::to_string(rand());
    storage->add(model);
    modelkeys.insert(model->keyval);

    std::vector< double > descriptor = getDescriptor(res,model);
    descriptors.push_back(std::make_pair(model->keyval,descriptor));
	return true;
}

bool ModelDatabaseRGBHistogram::remove(reglib::Model * model){
    storage->remove(model);
    modelkeys.erase(model->keyval);
    for(unsigned int i = 0; i < descriptors.size(); i++){
        if(descriptors[i].first.compare(model->keyval) == 0){
            descriptors[i] = descriptors.back();
            descriptors.pop_back();
            return true;
        }
    }
	return false;
}

double dist(std::vector< double > descriptorA, std::vector< double > descriptorB){
	double sum = 0;
	for(unsigned int i = 0; i < descriptorA.size(); i++){
		sum += std::min(descriptorA[i],descriptorB[i]);
	}
	return sum;
}

bool ModelDatabaseRGBHistogram_sort (std::pair<double,std::string> i, std::pair<double,std::string> j) { return i.first>j.first; }

std::vector<reglib::Model *> ModelDatabaseRGBHistogram::search(reglib::Model * model, int number_of_matches){
    std::vector<reglib::Model *> ret;
	std::vector< double > descriptor = getDescriptor(res,model);

    std::vector< std::pair<double,std::string> > distances;

    for(unsigned int i = 0; i < descriptors.size(); i++){
        std::string & key = descriptors[i].first;
        std::vector<double> & descriptor2 = descriptors[i].second;
        if(key.compare(model->keyval) != 0){
            double d = dist(descriptor,descriptor2);
            if(d > 0.5){
                distances.push_back(std::make_pair(d,key));
            }
        }
    }

    std::sort (distances.begin(), distances.end(), ModelDatabaseRGBHistogram_sort);

    for(unsigned int i = 0; i < distances.size() && i < number_of_matches; i++){
        printf("match to %s distance %f\n",distances[i].second.c_str(),distances[i].first);
        if(model == 0 || model->keyval.compare(distances[i].second) == 0){continue;}
        reglib::Model * mod = storage->fetch(distances[i].second);
        if(mod != 0 && mod->keyval.compare(model->keyval) != 0){ret.push_back(mod);}
        if(ret.size() == number_of_matches){return ret;}
    }
	return ret;
}

#ifndef ModelDatabase_H
#define ModelDatabase_H

#include <vector>

// PCL specific includes
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include "model/Model.h"
#include "../Util/Util.h"

#include "../ModelStorage/ModelStorage.h"
#include <set>

class ModelDatabase{
	public:
	std::vector<reglib::Model * > models;
	ModelStorage * storage;
    std::set<std::string> modelkeys;

	//Add pointcloud to database, return index number in database, weight is the bias of the system to perfer this object when searching
	virtual bool add(reglib::Model * model);
    virtual bool setStorage(ModelStorage * storage);
	// return true if successfull
    // return false if fail
    virtual bool remove(reglib::Model * model);
    virtual bool update(reglib::Model * model);
		
	//Find the number_of_matches closest matches in dabase to the pointcloud for index 
	virtual std::vector<reglib::Model * > search(reglib::Model * model, int number_of_matches);
	virtual std::vector<reglib::Model * > getBestModels(int nr = 20);
		
	ModelDatabase();
	~ModelDatabase();
};

#include "ModelDatabaseBasic.h"
#include "ModelDatabaseRGBHistogram.h"
#include "ModelDatabaseRetrieval.h"
#endif // ModelDatabase_H

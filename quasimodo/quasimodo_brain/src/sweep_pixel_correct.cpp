#include "ModelStorage/ModelStorage.h"
#include "Util/Util.h"
#include "CameraOptimizer/CameraOptimizer.h"

void update(std::string path, CameraOptimizer * cameraOptimizer){
	for(unsigned int i = 0; i < 8;i++){
		path.pop_back();
	}
	path += "/";
	printf("updating: %s\n",path.c_str());

	for(unsigned int i = 0; true; i++){
		char buf [1024];
		sprintf(buf,"intermediate_cloud%4.4i.pcd",i);
		std::string object = path+std::string(buf);
		sprintf(buf,"intermediate_cloud%4.4i_backup.pcd",i);
		std::string backupobject = path+std::string(buf);

		sprintf(buf,"%s/intermediate_cloud%4.4i_updated.txt",path.c_str(),i);
		//if(quasimodo_brain::fileExists(std::string(buf))){continue;}

		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGB> ());

		if(quasimodo_brain::fileExists(backupobject)){
			if (pcl::io::loadPCDFile (backupobject, *cloud) < 0)  {
				std::cout << "Error loading backup point cloud " << object << std::endl << std::endl;
				break;
			}
		}else if(quasimodo_brain::fileExists(object)){
			if (pcl::io::loadPCDFile (object, *cloud) < 0)  {
				std::cout << "Error loading point cloud " << object << std::endl << std::endl;
			}
			pcl::io::savePCDFileBinaryCompressed(backupobject,*cloud);
		}else{
			break;
		}

		unsigned long width = 640;
		unsigned long height = 480;
		for(unsigned long h = 0; h < height;h ++){
			for(unsigned long w = 0; w < width;w ++){
				unsigned long ind = h*width+w;
				pcl::PointXYZRGB & p = cloud->points[ind];
				if(p.z > 0 && !std::isnan(p.z)){
					double ratio = cameraOptimizer->getRange(float(w)/float(width),float(h)/float(height),p.z)/p.z;
					p.x *= ratio;
					p.y *= ratio;
					p.z *= ratio;
				}
			}
		}

		std::ofstream myfile;
		myfile.open (std::string(buf));
		myfile << "dummy";
		myfile.close();
		pcl::io::savePCDFileBinaryCompressed(object,*cloud);
	}
//	for (auto cloudFile : clouds){
//		std::string object = path+cloudFile.toStdString();
//		printf("cloud: %s\n",object.c_str());
//		std::string tmp = cloudFile.toStdString();
//		tmp.pop_back();
//		tmp.pop_back();
//		tmp.pop_back();
//		tmp.pop_back();
//		tmp += ".improved";

//		QStringList impoutput = QDir(path.c_str()).entryList(QStringList(tmp.c_str()));

//		if(impoutput.size() == 0){
//			pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGB> ());

//			  if (pcl::io::loadPCDFile (object, *cloud) < 0)  {
//				std::cout << "Error loading point cloud " << object << std::endl << std::endl;
//				continue;
//			  }

//			  unsigned int width = cloud->width;
//			  unsigned int height = cloud->height;
//			  printf("%i %i\n",width,height);

//			  //z_ratios[ind] = getRange(double(w)/double(width), double(h)/double(height),z,h % 20 == 0 && w % 20 == 0)/z;




////			processSweep(sweep_xml,"");

//		}
//	}
}

int main(int argc, char **argv){
	CameraOptimizer * cameraOptimizer = 0;

	bool visualize = false;
	int inputstate = 0;
	for(int i = 1; i < argc;i++){
		printf("input: %s\n",argv[i]);
			 if(std::string(argv[i]).compare("-v") == 0){						visualize = true; printf("visualize on\n");}
		else if(std::string(argv[i]).compare("-path") == 0){					inputstate = 0;}
		else if(std::string(argv[i]).compare("-optimizer") == 0){				inputstate = 1;}
		else if(inputstate == 0){
			std::string path = std::string(argv[i]);
			printf("path: %s\n",path.c_str());
			bool no_optimizer_set = cameraOptimizer == 0;
			if(no_optimizer_set){
				printf("loading camera optimizer: %s\n",(path+"/sweep_CameraOptimizerGridXYZ.bin").c_str());
				cameraOptimizer = CameraOptimizer::load(path+"/sweep_CameraOptimizerGridXYZ.bin");
			}

			if(visualize){cameraOptimizer->show(true);}

			std::vector<std::string> sweep_xmls = semantic_map_load_utilties::getSweepXmls<pcl::PointXYZRGB>(path);
			printf("number of sweep_xmls: %i\n",sweep_xmls.size());
			for (unsigned int i = 0; i < sweep_xmls.size(); i++) {
				printf("room: %i / %i\n",i+1,sweep_xmls.size());
				update(sweep_xmls[i],cameraOptimizer);
			}
			if(no_optimizer_set){delete cameraOptimizer;}
		}else if(inputstate == 1){
			if(cameraOptimizer != 0){delete  cameraOptimizer;}
			cameraOptimizer = CameraOptimizer::load(std::string(argv[i]));
		}

	}



	return 0;
}

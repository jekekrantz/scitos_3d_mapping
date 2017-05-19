#include "ModelDatabase/ModelDatabase.h"
#include "ModelStorage/ModelStorage.h"
#include "Util/Util.h"

//#include "/home/johane/koltuntest/FastGlobalRegistration/source/FastGlobalRegistration/app.cpp"
#include "reg2/RegistrationRandom2.cpp"
#include <pcl/features/fpfh_omp.h>


using namespace quasimodo_brain;

int state = -1;
void keyboardEventOccurred (const pcl::visualization::KeyboardEvent &event, void* viewer_void){
	pcl::visualization::PCLVisualizer *viewer = static_cast<pcl::visualization::PCLVisualizer *> (viewer_void);
	std::string s = event.getKeySym ();
	bool ed = event.keyDown ();
	if(ed){
		std::cout << s << std::endl;
	}
	if (s == "KP_7" && ed){state = 0;}
	if (s == "KP_8" && ed){state = 1;}
	if (s == "KP_4" && ed){state = 2;}
	if (s == "KP_5" && ed){state = 3;}
	if (s == "KP_0" && ed){state = 4;}
	if (s == "h" && ed){state = 5;}
	if (s == "j" && ed){state = 6;}
	if (s == "k" && ed){state = 7;}
	if (s == "l" && ed){state = 8;}
	if (s == "w" && ed){state = 9;}
//	if (event.getKeySym () == "a" && event.keyDown ()){state = 0;}
//	if (event.getKeySym () == "s" && event.keyDown ()){state = 1;}
//	if (event.getKeySym () == "d" && event.keyDown ()){state = 2;}
//	if (event.getKeySym () == "f" && event.keyDown ()){state = 3;}
//	if (event.getKeySym () == "g" && event.keyDown ()){state = 4;}
//	if (event.getKeySym () == "h" && event.keyDown ()){state = 5;}
//	if (event.getKeySym () == "j" && event.keyDown ()){state = 6;}
//	if (event.getKeySym () == "k" && event.keyDown ()){state = 7;}
//	if (event.getKeySym () == "l" && event.keyDown ()){state = 8;}
//	if (event.getKeySym () == " " && event.keyDown ()){state = 9;}
	//if (event.keyDown ()){printf("state: %i\n",state);}
}

Eigen::Affine3d randomrot = Eigen::Affine3d::Identity();

bool changed = false;

//void keyboardEventOccurred2 (const pcl::visualization::KeyboardEvent &event, void* viewer_void){
//	pcl::visualization::PCLVisualizer *viewer = static_cast<pcl::visualization::PCLVisualizer *> (viewer_void);
//	std::string s = event.getKeySym ();
//	bool ed = event.keyDown ();
//	if(ed){std::cout << s << std::endl;}

//	if (s == "KP_7" && ed){randomrot = Eigen::AngleAxisd(-0.01, Eigen::Vector3d::UnitX()) * randomrot;}
//	if (s == "KP_8" && ed){randomrot = Eigen::AngleAxisd( 0.01, Eigen::Vector3d::UnitX()) * randomrot;}

//	if (s == "KP_4" && ed){randomrot = Eigen::AngleAxisd(-0.01, Eigen::Vector3d::UnitY()) * randomrot;}
//	if (s == "KP_5" && ed){randomrot = Eigen::AngleAxisd( 0.01, Eigen::Vector3d::UnitY()) * randomrot;}

//	if (s == "KP_1" && ed){randomrot = Eigen::AngleAxisd(-0.01, Eigen::Vector3d::UnitZ()) * randomrot;}
//	if (s == "KP_2" && ed){randomrot = Eigen::AngleAxisd( 0.01, Eigen::Vector3d::UnitZ()) * randomrot;}
//	changed = true;
//}

pcl::PointCloud<pcl::PointXYZRGB>::Ptr color(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld, double r, double g, double b){
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr out (new pcl::PointCloud<pcl::PointXYZRGB>);
	out->points.resize(cld->points.size());

	for(unsigned int i = 0; i < cld->points.size(); i++){
		out->points[i].x = cld->points[i].x;
		out->points[i].y = cld->points[i].y;
		out->points[i].z = cld->points[i].z;
		out->points[i].r = r;
		out->points[i].g = g;
		out->points[i].b = b;
	}

	return out;
}

pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr color(pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cld, double r, double g, double b){
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr out (new pcl::PointCloud<pcl::PointXYZRGBNormal>);
	out->points.resize(cld->points.size());

	for(unsigned int i = 0; i < cld->points.size(); i++){
		out->points[i]	 = cld->points[i];
		out->points[i].r = r;
		out->points[i].g = g;
		out->points[i].b = b;
	}

	return out;
}

void center(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld,int offset){
	double x = 0;
	double y = 0;
	double z = 0;
	for(unsigned int i = 1; i < cld->points.size(); i++){
		x += cld->points[i].x;
		y += cld->points[i].y;
		z += cld->points[i].z;
	}
	x /= double(cld->points.size());
	y /= double(cld->points.size());
	z /= double(cld->points.size());
	for(unsigned int i = 0; i < cld->points.size(); i++){
		cld->points[i].x -= x-offset;
		cld->points[i].y -= y;
		cld->points[i].z -= z;
	}
}

void center(pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cld){
	double x = 0;
	double y = 0;
	double z = 0;
	for(unsigned int i = 1; i < cld->points.size(); i++){
		x += cld->points[i].x;
		y += cld->points[i].y;
		z += cld->points[i].z;
	}
	x /= double(cld->points.size());
	y /= double(cld->points.size());
	z /= double(cld->points.size());
	for(unsigned int i = 0; i < cld->points.size(); i++){
		cld->points[i].x -= x;
		cld->points[i].y -= y;
		cld->points[i].z -= z;
	}
}



void testRegister(reglib::Model * model1, reglib::Model * model2, boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer){
//	Eigen::Affine3d startrot = Eigen::Affine3d::Identity();
//	startrot = Eigen::AngleAxisd(-30*2*M_PI/360.0, Eigen::Vector3d::UnitX());

//	transformSuperPoints(model1->points,startrot.matrix());
//	transformSuperPoints(model2->points,startrot.matrix());

//	model1->recomputeModelPoints(startrot.matrix());
//	model2->recomputeModelPoints(startrot.matrix());


/*
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr model1Cloud = model1->getPCLnormalcloud(1,false);
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr model2Cloud = model2->getPCLnormalcloud(1,false);

	printf("%i %i -> ",model1Cloud->points.size(),model2Cloud->points.size());

	double feature_radius_ = 0.10;

	pcl::PointCloud<pcl::FPFHSignature33>::Ptr model1_features(new pcl::PointCloud<pcl::FPFHSignature33>());
	pcl::PointCloud<pcl::FPFHSignature33>::Ptr model2_features(new pcl::PointCloud<pcl::FPFHSignature33>());

	pcl::FPFHEstimationOMP<pcl::PointXYZRGBNormal, pcl::PointXYZRGBNormal, pcl::FPFHSignature33> fest;
	fest.setRadiusSearch(feature_radius_);

	fest.setInputCloud(model1Cloud);
	fest.setInputNormals(model1Cloud);
	fest.compute(*model1_features);

	fest.setInputCloud(model2Cloud);
	fest.setInputNormals(model2Cloud);
	fest.compute(*model2_features);


	printf("%i %i\n",model1Cloud->points.size(),model2Cloud->points.size());

	std::vector<Eigen::Vector3f> model1_xyx;
	std::vector<Eigen::VectorXf> model1_feature;
	for (int v = 0; v < model1Cloud->points.size(); v++) {
		const pcl::PointXYZRGBNormal & pt = model1Cloud->points[v];
		model1_xyx.push_back(Eigen::Vector3f(pt.x, pt.y, pt.z));
		Eigen::VectorXf f (33);
		const pcl::FPFHSignature33 &feature = model1_features->points[v];
		for (int j = 0; j < 33; j++) {f(j) = feature.histogram[j];}
		model1_feature.push_back(f);
	}

	std::vector<Eigen::Vector3f> model2_xyx;
	std::vector<Eigen::VectorXf> model2_feature;
	for (int v = 0; v < model2Cloud->points.size(); v++) {
		const pcl::PointXYZRGBNormal & pt = model2Cloud->points[v];
		model2_xyx.push_back(Eigen::Vector3f(pt.x, pt.y, pt.z));
		Eigen::VectorXf f (33);
		const pcl::FPFHSignature33 &feature = model2_features->points[v];
		for (int j = 0; j < 33; j++) {f(j) = feature.histogram[j];}
		model2_feature.push_back(f);
	}

	CApp app;
	app.LoadFeature(model1_xyx,model1_feature);
	app.LoadFeature(model2_xyx,model2_feature);
//	app.ReadFeature(argv[2]);
	app.NormalizePoints();
	app.AdvancedMatching();
	app.OptimizePairwise(true, ITERATION_NUMBER);
	Eigen::Matrix4f T = app.GetTrans();

	std::cout << T << std::endl;


*/



//	Eigen::Affine3d current_guess = Ymean*randomrot*Xmean.inverse();



//	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr model2Cloud = model2->getPCLnormalcloud(1,false);
//	center(model2Cloud);
//	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer2 = boost::shared_ptr<pcl::visualization::PCLVisualizer>(new pcl::visualization::PCLVisualizer ("Modelserver Viewer"));
//	viewer2->addCoordinateSystem(1.0);
//	viewer2->setBackgroundColor(1.0,1.0,1.0);
//	viewer2->registerKeyboardCallback (keyboardEventOccurred2, (void*)viewer.get ());

//	viewer2->removeAllPointClouds();
//	viewer2->addPointCloud<pcl::PointXYZRGBNormal>(model2Cloud,pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(model2Cloud), "cld2" );

//	while (!viewer->wasStopped () ){
//		viewer->spinOnce (10);
//		boost::this_thread::sleep (boost::posix_time::microseconds (1000));
//		if(changed){
//			Eigen::Affine3d randomrot = Eigen::Affine3d::Identity();
//			randomrot =	Eigen::AngleAxisd(rx, Eigen::Vector3d::UnitX()) *
//						Eigen::AngleAxisd(ry, Eigen::Vector3d::UnitY()) *
//						Eigen::AngleAxisd(rz, Eigen::Vector3d::UnitZ());
//			pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr tc (new pcl::PointCloud<pcl::PointXYZRGBNormal> ());
//			pcl::transformPointCloud (*model2Cloud, *tc, randomrot);
//			viewer2->removeAllPointClouds();
//			viewer2->addPointCloud<pcl::PointXYZRGBNormal>(tc,pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(tc), "cld2" );
//			changed = false;
//		}
//	}


//	pcl::PointCloud<pcl::PointXYZRGB>::Ptr model2Cloud = model2->submodels.front()->frames.front()->getSmallPCLcloud();//model2->getPCLnormalcloud(1,false);

//	//center(model2Cloud);
//	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer2 = boost::shared_ptr<pcl::visualization::PCLVisualizer>(new pcl::visualization::PCLVisualizer ("Modelserver Viewer"));
//	viewer2->addCoordinateSystem(1.0);
//	viewer2->setBackgroundColor(1.0,1.0,1.0);
////	viewer2->registerKeyboardCallback (keyboardEventOccurred2, (void*)viewer.get ());

//	viewer2->removeAllPointClouds();
//	viewer2->addPointCloud<pcl::PointXYZRGB>(model2Cloud,pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(model2Cloud), "cld2" );

//	while (!viewer->wasStopped () ){
//		viewer->spinOnce (10);
//		boost::this_thread::sleep (boost::posix_time::microseconds (1000));
//		if(changed){
//			pcl::PointCloud<pcl::PointXYZRGB>::Ptr tc (new pcl::PointCloud<pcl::PointXYZRGB> ());
//			pcl::transformPointCloud (*model2Cloud, *tc, randomrot);
//			viewer2->removeAllPointClouds();
//			viewer2->addPointCloud<pcl::PointXYZRGB>(tc,pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(tc), "cld2" );
//			changed = false;
//		}
//	}

//	return;

//	for(double rx = 0; rx < 360; rx += 1){
//		printf("rx: %f\n",rx);
//		Eigen::Affine3d randomrot = Eigen::Affine3d::Identity();
//		randomrot =	Eigen::AngleAxisd(-rx*2*M_PI/360.0, Eigen::Vector3d::UnitX()) *
//				Eigen::AngleAxisd(0, Eigen::Vector3d::UnitY()) *
//				Eigen::AngleAxisd(0, Eigen::Vector3d::UnitZ());
//		pcl::PointCloud<pcl::PointXYZRGB>::Ptr tc (new pcl::PointCloud<pcl::PointXYZRGB> ());
//		pcl::transformPointCloud (*model2Cloud, *tc, randomrot);
//		viewer2->removeAllPointClouds();
//		viewer2->addPointCloud<pcl::PointXYZRGB>(tc,pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(tc), "cld2" );
//		viewer2->spin();
//	}

//exit(0);
//	for(double ry = 0; ry < 2*M_PI; ry += 0.01){
//		printf("ry: %f\n",ry);
//		Eigen::Affine3d randomrot = Eigen::Affine3d::Identity();
//		randomrot =	Eigen::AngleAxisd(0, Eigen::Vector3d::UnitX()) *
//				Eigen::AngleAxisd(ry, Eigen::Vector3d::UnitY()) *
//				Eigen::AngleAxisd(0, Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(-30*2*M_PI/360.0, Eigen::Vector3d::UnitX());
//		pcl::PointCloud<pcl::PointXYZRGB>::Ptr tc (new pcl::PointCloud<pcl::PointXYZRGB> ());
//		pcl::transformPointCloud (*model2Cloud, *tc, randomrot);
//		viewer2->removeAllPointClouds();
//		viewer2->addPointCloud<pcl::PointXYZRGB>(tc,pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(tc), "cld2" );
//		viewer2->spin();
//	}

//	for(double rz = 0; rz < 2*M_PI; rz += 0.01){
//		printf("rz: %f\n",rz);
//		Eigen::Affine3d randomrot = Eigen::Affine3d::Identity();
//		randomrot =	Eigen::AngleAxisd(0, Eigen::Vector3d::UnitX()) *
//				Eigen::AngleAxisd(0, Eigen::Vector3d::UnitY()) *
//				Eigen::AngleAxisd(rz, Eigen::Vector3d::UnitZ());
//		pcl::PointCloud<pcl::PointXYZRGB>::Ptr tc (new pcl::PointCloud<pcl::PointXYZRGB> ());
//		pcl::transformPointCloud (*model2Cloud, *tc, randomrot);
//		viewer2->removeAllPointClouds();
//		viewer2->addPointCloud<pcl::PointXYZRGB>(tc,pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(tc), "cld2" );
//		viewer2->spin();
//	}


	reglib::RegistrationRandom2 *	reg	= new reglib::RegistrationRandom2(5);
	//reg->initTransform = (Eigen::AngleAxisd(-30*2*M_PI/360.0, Eigen::Vector3d::UnitX())).matrix();
	reg->steprx = 1;
	reg->stepry	= 36;
	reg->steprz	= 1;
//	reg->src_meantype	= 3;
//	reg->dst_meantype	= 3;
//	reg->start_rx	= 0;//2*M_PI*20.0/360;
//	reg->start_ry	= 0;//2*M_PI*20.0/360;
//	reg->start_rz	= 0;


	//reg->stop_rx	= reg->start_rx;
	//reg->stop_ry	= reg->start_ry;
	//reg->stop_rz	= reg->start_rz;


	reg->visualizationLvl				= 1;
	reg->viewer							= viewer;
	reglib::ModelUpdaterBasicFuse * mu	= new reglib::ModelUpdaterBasicFuse( model2, reg);
	mu->occlusion_penalty               = 15;
	mu->viewer							= viewer;
	mu->show_scoring					= true;//fuse scoring show


	reglib::FusionResults fr = mu->registerModel(model1);
	//		if(fr.score > 100){
	//			reglib::UpdatedModels ud = mu->fuseData(&fr, model2H, model1H);
	//		}
	delete mu;
	delete reg;
}





int main(int argc, char **argv){

	ros::init(argc, argv, "testmerge");
	ros::NodeHandle n;

	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = boost::shared_ptr<pcl::visualization::PCLVisualizer>(new pcl::visualization::PCLVisualizer ("Modelserver Viewer"));
	viewer->addCoordinateSystem(0.1);
	viewer->setBackgroundColor(1.0,0.0,1.0);
	viewer->registerKeyboardCallback (keyboardEventOccurred, (void*)viewer.get ());
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer2 = boost::shared_ptr<pcl::visualization::PCLVisualizer>(new pcl::visualization::PCLVisualizer ("Modelserver Viewer"));
	viewer2->addCoordinateSystem(1.0);
	viewer2->setBackgroundColor(1.0,1.0,1.0);

	int v1(0);
	viewer2->createViewPort(0.0, 0.0, 0.5, 1.0, v1);
	viewer2->setBackgroundColor (0.5, 0.5, 0.5, v1);

	int v2(0);
	viewer2->createViewPort(0.5, 0.0, 1.0, 1.0, v2);
	viewer2->setBackgroundColor (0.7, 0.7, 0.7, v2);

	ModelStorageFile *  storage = new ModelStorageFile();
	ModelDatabase * modeldatabase	= new ModelDatabaseRetrieval(n);
	modeldatabase->setStorage(storage);

	std::vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> models;
	std::vector<reglib::Model * > mods;

	int offset = 0;
	for (std::map<std::string,std::string>::iterator it=storage->keyPathMap.begin(); it!=storage->keyPathMap.end(); ++it){
		reglib::Model * model = storage->fetch(it->first);
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld	= model->getPCLcloud(1,false);
		center(cld,offset++);

		models.push_back(cld);
		mods.push_back(model);
		//storage->fullHandback();
	}

//	for(unsigned int i = 0; i < models.size(); i++){
//		pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr model2Cloud = mods[i]->getPCLnormalcloud(1,false);
//		center(model2Cloud);
//		boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer2 = boost::shared_ptr<pcl::visualization::PCLVisualizer>(new pcl::visualization::PCLVisualizer ("Modelserver Viewer"));
//		viewer2->addCoordinateSystem(1.0);
//		viewer2->setBackgroundColor(1.0,1.0,1.0);
//		viewer2->registerKeyboardCallback (keyboardEventOccurred2, (void*)viewer.get ());

//		viewer2->removeAllPointClouds();
//		viewer2->addPointCloud<pcl::PointXYZRGBNormal>(model2Cloud,pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(model2Cloud), "cld2" );

//		while (!viewer->wasStopped () ){
//			viewer->spinOnce (10);
//			boost::this_thread::sleep (boost::posix_time::microseconds (1000));
//			if(changed){
//				Eigen::Affine3d randomrot = Eigen::Affine3d::Identity();
//				randomrot =	Eigen::AngleAxisd(rx, Eigen::Vector3d::UnitX()) *
//							Eigen::AngleAxisd(ry, Eigen::Vector3d::UnitY()) *
//							Eigen::AngleAxisd(rz, Eigen::Vector3d::UnitZ());
//				pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr tc (new pcl::PointCloud<pcl::PointXYZRGBNormal> ());
//				pcl::transformPointCloud (*model2Cloud, *tc, randomrot);
//				viewer2->removeAllPointClouds();
//				viewer2->addPointCloud<pcl::PointXYZRGBNormal>(tc,pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(tc), "cld2" );
//				changed = false;
//			}
//		}
//	}


	char buf[1024];
	for(unsigned int i = 0; i < models.size(); i++){

		sprintf(buf,"model%i",i);
		viewer->addPointCloud<pcl::PointXYZRGB> (models[i], pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(models[i]), buf );
	}

	int current = 54;
	int other = 55;

    testRegister(mods[current],mods[other],viewer2);


	sprintf(buf,"model%i",current);
	viewer->removePointCloud(buf);
	viewer->addPointCloud<pcl::PointXYZRGB> (color(models[current],0,255,0), pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(color(models[current],0,255,0)), buf );

	sprintf(buf,"model%i",other);
	viewer->removePointCloud(buf);
	viewer->addPointCloud<pcl::PointXYZRGB> (color(models[other],0,0,255), pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(color(models[other],0,0,255)), buf );

	while (!viewer->wasStopped () ){
		viewer->spinOnce (10);
		boost::this_thread::sleep (boost::posix_time::microseconds (1000));
		if(state != -1){
			if(state == 0 && current > 0){
				sprintf(buf,"model%i",current);
				viewer->removePointCloud(buf);
				viewer->addPointCloud<pcl::PointXYZRGB> (models[current], pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(models[current]), buf );
				current--;
				sprintf(buf,"model%i",current);
				viewer->removePointCloud(buf);
				viewer->addPointCloud<pcl::PointXYZRGB> (color(models[current],0,255,0), pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(color(models[current],0,255,0)), buf );
			}


			if(state == 1 && current < (models.size()+1)){
				sprintf(buf,"model%i",current);
				viewer->removePointCloud(buf);
				viewer->addPointCloud<pcl::PointXYZRGB> (models[current], pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(models[current]), buf );
				current++;
				sprintf(buf,"model%i",current);
				viewer->removePointCloud(buf);
				viewer->addPointCloud<pcl::PointXYZRGB> (color(models[current],0,255,0), pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(color(models[current],0,255,0)), buf );
			}


			if(state == 2 && other > 0){
				sprintf(buf,"model%i",other);
				viewer->removePointCloud(buf);
				viewer->addPointCloud<pcl::PointXYZRGB> (models[other], pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(models[other]), buf );
				other--;
				sprintf(buf,"model%i",other);
				viewer->removePointCloud(buf);
				viewer->addPointCloud<pcl::PointXYZRGB> (color(models[other],0,255,0), pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(color(models[other],0,255,0)), buf );
			}


			if(state == 3 && other < (models.size()+1)){
				sprintf(buf,"model%i",other);
				viewer->removePointCloud(buf);
				viewer->addPointCloud<pcl::PointXYZRGB> (models[other], pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(models[other]), buf );
				other++;
				sprintf(buf,"model%i",other);
				viewer->removePointCloud(buf);
				viewer->addPointCloud<pcl::PointXYZRGB> (color(models[other],0,0,255), pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(color(models[other],0,0,255)), buf );
			}

			printf("current: %i other: %i",current,other);
			if(state == 4){
				testRegister(mods[current],mods[other],viewer2);
			}
			state = -1;
		}
	}


	//		std::vector<reglib::Model * > res = modeldatabase->search(model,5);


	//		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld1	= model->getPCLcloud(1,false);
	//		viewer->addPointCloud<pcl::PointXYZRGB> (cld1, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld1), "model");

	//		for(unsigned int i = 0; i < res.size(); i++){
	//			pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld2	= res[i]->getPCLcloud(1,false);
	//			for(unsigned int j = 0; j < cld2->points.size(); j++){
	//				cld2->points[j].x += i+2;
	//			}

	//			char buf[1024];
	//			sprintf(buf,"model%i",i);
	//			viewer->addPointCloud<pcl::PointXYZRGB> (cld2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld2), buf );
	//		}


	//		while (!viewer->wasStopped () ){
	//			viewer->spinOnce (100);
	//			boost::this_thread::sleep (boost::posix_time::microseconds (100000));
	//			if(state >= 0){
	//				printf("state: %i\n",state);
	//				if(state == 9){
	//					state = -1;
	//					break;
	//				}
	//				testRegister(model,res[state],viewer);
	//				state = -1;
	//			}
	//		}

	//		storage->fullHandback();



//	for (std::map<std::string,std::string>::iterator it=storage->keyPathMap.begin(); it!=storage->keyPathMap.end(); ++it){
//		reglib::Model * model = storage->fetch(it->first);
//		std::vector<reglib::Model * > res = modeldatabase->search(model,5);

//		viewer->removeAllPointClouds();
//		viewer->removeAllShapes();
//		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld1	= model->getPCLcloud(1,false);
//		viewer->addPointCloud<pcl::PointXYZRGB> (cld1, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld1), "model");

//		for(unsigned int i = 0; i < res.size(); i++){
//			pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld2	= res[i]->getPCLcloud(1,false);
//			for(unsigned int j = 0; j < cld2->points.size(); j++){
//				cld2->points[j].x += i+2;
//			}

//			char buf[1024];
//			sprintf(buf,"model%i",i);
//			viewer->addPointCloud<pcl::PointXYZRGB> (cld2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld2), buf );
//		}


//		while (!viewer->wasStopped () ){
//			viewer->spinOnce (100);
//			boost::this_thread::sleep (boost::posix_time::microseconds (100000));
//			if(state >= 0){
//				printf("state: %i\n",state);
//				if(state == 9){
//					state = -1;
//					break;
//				}
//				testRegister(model,res[state],viewer);
//				state = -1;
//			}
//		}

//		storage->fullHandback();

//	}
	/*
	reglib::Model * model1						= reglib::Model::loadFast(std::string(argv[1]));
	reglib::Model * model2						= reglib::Model::loadFast(std::string(argv[2]));

	reglib::Model * model1H = new reglib::Model();
	model1->parrent = model1H;
	model1H->submodels.push_back(model1);
	model1H->submodels_relativeposes.push_back(Eigen::Matrix4d::Identity());
	model1H->recomputeModelPoints();

	reglib::Model * model2H = new reglib::Model();
	model2->parrent = model2H;
	model2H->submodels.push_back(model2);
	model2H->submodels_relativeposes.push_back(Eigen::Matrix4d::Identity());
	model2H->recomputeModelPoints();

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld1	= model1H->getPCLcloud(1,false);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld2	= model2H->getPCLcloud(1,false);

	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = boost::shared_ptr<pcl::visualization::PCLVisualizer>(new pcl::visualization::PCLVisualizer ("Modelserver Viewer"));
	viewer->addCoordinateSystem(0.1);
	viewer->setBackgroundColor(1.0,1.0,1.0);
	viewer->addPointCloud<pcl::PointXYZRGB> (cld1, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld1), "model1");
	viewer->addPointCloud<pcl::PointXYZRGB> (cld2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld2), "model2");
	printf("models added\n");
	viewer->spin();
	{
		reglib::RegistrationRandom *	reg	= new reglib::RegistrationRandom(5);
		reglib::ModelUpdaterBasicFuse * mu	= new reglib::ModelUpdaterBasicFuse( model2H, reg);
		mu->occlusion_penalty               = 15;
		mu->viewer							= viewer;
		mu->show_scoring					= false;//fuse scoring show
		reg->visualizationLvl				= 0;

		reglib::FusionResults fr = mu->registerModel(model1H);
//		if(fr.score > 100){
//			reglib::UpdatedModels ud = mu->fuseData(&fr, model2H, model1H);

//		}
		delete mu;
		delete reg;
	}

	{
		reglib::RegistrationRandom *	reg	= new reglib::RegistrationRandom(5);
		reglib::ModelUpdaterBasicFuse * mu	= new reglib::ModelUpdaterBasicFuse( model1H, reg);
		mu->occlusion_penalty               = 15;
		mu->viewer							= viewer;
		mu->show_scoring					= false;//fuse scoring show
		reg->visualizationLvl				= 0;

		reglib::FusionResults fr = mu->registerModel(model2H);
//		if(fr.score > 100){
//			reglib::UpdatedModels ud = mu->fuseData(&fr, model1H, model2H);

//		}
		delete mu;
		delete reg;
	}
*/
	return 0;
}

#include "../../include/core/superpoint.h"

namespace reglib
{

//pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr superpoint::getCloud(){}

//pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr	cloud	(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
//cloud->width	= width;
//cloud->height	= height;
//cloud->points.resize(width*height);


double mysign(double v){
	if(v < 0){return -1;}
	return 1;
}

superpoint::superpoint(bool dummy){
	x	= 0;
	y	= 0;
	z	= 0;

	nx	= 0;
	ny	= 0;
	nz	= -1;

	r	= 0;
	g	= 0;
	b	= 0;

	point_information = 0;
	normal_information = 0;
	colour_information = 0;
	last_update_frame_id = 0;

	is_boundry = true;
}

superpoint::superpoint(Eigen::Vector3f p, Eigen::Vector3f n, Eigen::VectorXf f, double pi, double fi, int id, bool is_boundry_){
    x	= p(0);
    y	= p(1);
    z	= p(2);

    nx	= n(0);
    ny	= n(1);
    nz	= n(2);

    r	= f(0);
    g	= f(1);
    b	= f(2);

    point_information = pi;
    normal_information = pi;
    colour_information = fi;
    last_update_frame_id = id;

    is_boundry = is_boundry_;
}

superpoint::~superpoint(){}

void superpoint::merge(superpoint p, double weight){

    double newpweight = weight*p.point_information	+ point_information;
    x	= weight*p.point_information*p.x		+ point_information*x;
    y	= weight*p.point_information*p.y		+ point_information*y;
    z	= weight*p.point_information*p.z		+ point_information*z;
    x	/= newpweight;
    y	/= newpweight;
    z	/= newpweight;
    point_information = newpweight;

    double newnweight = weight*p.point_information	+ point_information;

    nx	= weight*p.normal_information*p.nx		+ normal_information*nx;
    ny	= weight*p.normal_information*p.ny		+ normal_information*ny;
    nz	= weight*p.normal_information*p.nz		+ normal_information*nz;
    double nnorm = sqrt(nx*nx+ny*ny+nz*nz);
    if(nnorm != 0){
        nx /= nnorm;
        ny /= nnorm;
        nz /= nnorm;
    }
    normal_information = newnweight;


    double newcweight = weight*p.colour_information	+ colour_information;
    r	= weight*p.colour_information*p.r		+ colour_information*r;
    g	= weight*p.colour_information*p.g		+ colour_information*g;
    b	= weight*p.colour_information*p.b		+ colour_information*b;
    r	/= newcweight;
    g	/= newcweight;
    b	/= newcweight;
    colour_information = newcweight;

    last_update_frame_id = std::max(p.last_update_frame_id,last_update_frame_id);

    if(!p.is_boundry){is_boundry = false;}
}

void superpoint::print(){
	printf("point: %5.5f %5.5f %5.5f normal: %5.5f %5.5f %5.5f Colour: %5.5f %5.5f %5.5f Infos: %5.5f %5.5f %5.5f\n",x,y,z,nx,ny,nz,r,g,b,point_information,normal_information,colour_information);
}

void superpoint::transform(Eigen::Matrix4d cp){
	double m00 = cp(0,0); double m01 = cp(0,1); double m02 = cp(0,2); double m03 = cp(0,3);
	double m10 = cp(1,0); double m11 = cp(1,1); double m12 = cp(1,2); double m13 = cp(1,3);
	double m20 = cp(2,0); double m21 = cp(2,1); double m22 = cp(2,2); double m23 = cp(2,3);

	double tx	= m00*x + m01*y + m02*z + m03;
	double ty	= m10*x + m11*y + m12*z + m13;
	double tz	= m20*x + m21*y + m22*z + m23;
    x = tx;
    y = ty;
    z = tz;

	double tnx	= m00*nx + m01*ny + m02*nz;
	double tny	= m10*nx + m11*ny + m12*nz;
	double tnz	= m20*nx + m21*ny + m22*nz;
    nx = tnx;
    ny = tny;
    nz = tnz;
}

double superpoint::angle(superpoint p){return nx*p.nx + ny*p.ny + nz*p.nz;}
double superpoint::distance(superpoint p){
	double signval = mysign(z-p.z);
	double dx = p.x-x;
	double dy = p.y-y;
	double dz = p.z-z;
	double pdx = p.nx*dx;
	double pdy = p.ny*dy;
	double pdz = p.nz*dz;

//	printf("superpoint::distance -> %1.1f\n",signval);
//	printf("dx %8.8f = (%8.8f - %8.8f)\n",dx,p.x,x);
//	printf("dy %8.8f = (%8.8f - %8.8f)\n",dy,p.y,y);
//	printf("dz %8.8f = (%8.8f - %8.8f)\n",dz,p.z,z);
//	printf("--\n");
//	printf("pdx %8.8f = (%8.8f * %8.8f)\n",pdx,p.nx,dx);
//	printf("pdy %8.8f = (%8.8f * %8.8f)\n",pdy,p.ny,dy);
//	printf("pdz %8.8f = (%8.8f * %8.8f)\n",pdz,p.nz,dz);
//	printf("sum: %8.8f\n",pdx+pdy+pdz);
//	printf("superpoint::distance -> %1.1f * (%8.8f + %8.8f + %8.8f) from %8.8f = (%8.8f - %8.8f), %8.8f = (%8.8f - %8.8f), %8.8f = (%8.8f - %8.8f)\n",signval,pdx,pdy,pdz,dx,x,p.x,dy,y,p.y,dz,z,p.z);
	return signval*fabs(pdx+pdy+pdz);
	//return mysign(z-p.z)*fabs(nx*(x-p.x) + ny*(y-p.y) + nz*(z-p.z));//dst_z-tz;//mysign(dst_z-tz)*fabs(tnx*(dst_x-tx) + tny*(dst_y-ty) + tnz*(dst_z-tz));
}

}


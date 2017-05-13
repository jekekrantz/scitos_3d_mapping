#!/usr/bin/python

import rospy
from soma_llsd_msgs.msg import Scene, Segment, Observation
from soma_llsd.srv import * #InsertScene, InsertSegment, AddObservationsToSegment
from quasimodo_conversions.srv import * #soma_insert_model, soma_insert_modelRequest, soma_insert_modelResponse 

class SomaInsertModelServer(object):

    def __init__(self):

        self.srv = rospy.Service("/quasimodo_conversions/insert_soma_model", soma_insert_model, self.callback)

    def callback(self, req):
        
         

        #string episode_id
        #string waypoint
        #string meta_data
        #uint32 timestamp

        #tf/tfMessage transform
        #sensor_msgs/PointCloud2 cloud
        #sensor_msgs/Image rgb_img
        #sensor_msgs/Image depth_img
        #sensor_msgs/CameraInfo camera_info
        #geometry_msgs/Pose robot_pose
        #---
        #bool result
        #soma_llsd_msgs/Scene response

        observations = []
        for cloud, mask, pose, frame in zip(req.model.clouds, req.model.masks, req.model.local_poses, req.model.frames):

            scene = InsertSceneRequest()
            scene.cloud = cloud
            scene.rgb_img = frame.rgb
            scene.depth_img = frame.depth
            scene.camera_info = frame.camera
            robot_pose = pose

            rospy.wait_for_service("/soma_llsd/insert_scene")
            try:
                insert_scene = rospy.ServiceProxy("/soma_llsd/insert_scene", InsertScene)
                resp = insert_scene(scene)
            except rospy.ServiceException, e:
                print "Service call failed: %s"%e
                return False

            #string id
            #int32 timestamp
            #string meta_data
            #string scene_id

            #geometry_msgs/Pose pose # centroid in map co-ordinates
            #sensor_msgs/PointCloud2 map_cloud   #cloud in map co-ordinates
            #sensor_msgs/PointCloud2 camera_cloud # cloud in camera co-ordinates
            #sensor_msgs/PointCloud2 room_cloud  # cloud aligned to meta-room

            #sensor_msgs/Image rgb_cropped     # optional, required only for web services
            #sensor_msgs/Image depth_cropped   # optional, required only for web services
            #sensor_msgs/Image image_mask

            observation = Observation()
            observation.pose = pose
            #observation.camera_cloud = cloud
            observation.image_mask = mask
            observation.scene_id = resp.response.id
            scene_id = resp.response.id
            observations.append(observation)
        
        #string meta_data
        #string scene_id
        #soma_llsd_msgs/Observation[] observations       # optional, can be added later
        #---
        #bool result
        #soma_llsd_msgs/Segment response
            
        segment = InsertSegmentRequest()
        segment.scene_id = scene_id
        segment.observations = observations
  
        rospy.wait_for_service("/soma_llsd/insert_segment")
        try:
            insert_segment = rospy.ServiceProxy("/soma_llsd/insert_segment", InsertSegment)
            segment_resp = insert_segment(segment)
        except rospy.ServiceException, e:
            print "Service call failed: %s"%e
            return False
 
        ##string segment_id
        ##soma_llsd_msgs/Observation[] observations
        ##string scene_id
        ##---
        ##bool result
        #
        #observations_req = AddObservationsToSegmentRequest()
        #observations_req.segment_id = segment_resp.response.id
        ##observations_req.scene_id = segment_resp.response.scene_id
        #observations_req.observations = observations
        #
        #rospy.wait_for_service("/soma_llsd/insert_segment")
        #try:
        #    insert_segment = rospy.ServiceProxy("/soma_llsd/insert_segment", InsertSegment)
        #    resp = insert_segment(segment)
        #except rospy.ServiceException, e:
        #    print "Service call failed: %s"%e
        #    return False

        soma_resp = soma_insert_modelResponse()
        soma_resp.segment = segment_resp.response

        return soma_resp


if __name__ == "__main__":

    rospy.init_node('soma_insert_model_server')
    sims = SomaInsertModelServer()
    rospy.spin()


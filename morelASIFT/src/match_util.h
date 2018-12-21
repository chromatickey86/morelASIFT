#include <opencv2/opencv.hpp>

#ifndef _MATCH_UNIT_H_
#define _MATCH_UNIT_H_

void determine_good_match_by_brute_force(const cv::Mat &descriptors_object, const cv::Mat &descriptors_scene,
                          std::vector<cv::DMatch> &good_matches);

void determine_good_match_by_brute_force(const cv::Mat &descriptors_object, const cv::Mat &descriptors_scene, const std::vector<cv::KeyPoint> &lkps, const std::vector<cv::KeyPoint> &rkps, std::vector<cv::DMatch> &good_matches);

void matched_keypoint_to_mat(const std::vector<cv::DMatch> &good_matches,
                             const std::vector<cv::KeyPoint> &keypoints_object, const std::vector<cv::KeyPoint> &keypoints_scene,
                             cv::Mat &obj_mat, cv::Mat &scene_mat);

void draw_transformed_box(const cv::Mat &corners, const cv::Point2f &offset, cv::Mat &img, const cv::Scalar &color=cv::Scalar(0, 255, 0), int thick=4);

void find_transform(const cv::Mat &lhs_image, const std::vector<cv::KeyPoint> &lhs_kps, const cv::Mat &lhs_desc,
		const cv::Mat &rhs_image, const std::vector<cv::KeyPoint> &rhs_kps, const cv::Mat &rhs_desc);

void match_keypoint(const cv::Mat &lhs_image, const std::vector<cv::KeyPoint> &lhs_kps,
		const cv::Mat &rhs_image, const std::vector<cv::KeyPoint> &rhs_kps);

#endif

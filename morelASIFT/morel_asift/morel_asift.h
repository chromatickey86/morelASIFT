#ifndef __MOREL_ASIFT_H__
#define __MOREL_ASIFT_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>
#include <opencv2/opencv.hpp>

#include "demo_lib_sift.h"

#include "library.h"
#include "frot.h"
#include "fproj.h"
#include "compute_asift_keypoints.h"
#include "MatchingModel.hpp"

void morel_asift(cv::Mat img, std::vector<cv::KeyPoint> &keypoints,
					  std::vector<KeyPoint_param> &prmKeypoints, cv::Mat &descriptors);

void morel_asift(cv::Mat img, std::vector<MatchingModel::View> &views,
					  std::vector<cv::KeyPoint> &keypoints,
					  std::vector<KeyPoint_param> &prmKeypoints, cv::Mat &descriptors);

#endif

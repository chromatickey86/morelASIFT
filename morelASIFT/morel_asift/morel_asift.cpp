#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>
#include <opencv2/opencv.hpp>
//#include "asift.h"
#include "morel_asift.h"
#include "demo_lib_sift.h"
#include "evaluation_kit.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void morel_asift(cv::Mat img, std::vector<cv::KeyPoint> &keypoints,
					  std::vector<KeyPoint_param> &prmKeypoints, cv::Mat &descriptors) {
	cv::Mat src;
	
	//////////////////////////////////////////////// Input
	// Read image1
	int w, h;
	w = img.cols;
	h = img.rows;
	cvtColor(img, src, CV_BGR2GRAY);
	std::vector<float> ipixels(w * h);
	
	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols; j++)
			ipixels[i*img.cols + j] = (float)src.at<uchar>(i, j);
	
	src.release();
	
	///// Resize the images to area wS*hW in remaining the apsect-ratio
	///// Resize if the resize flag is not set or if the flag is set unequal to 0
	
	float zoom = 0;
	int wS = 0, hS = 0;
	std::vector<float> ipixels_zoom;
	
	ipixels_zoom.resize(w*h);
	ipixels_zoom = ipixels;
	wS = w;
	hS = h;
	zoom = 1;
	
	///// Compute ASIFT keypoints
	// number N of tilts to simulate t = 1, \sqrt{2}, (\sqrt{2})^2, ..., {\sqrt{2}}^(N-1)
	// t = 1, sqrt(2), 2, 2*sqrt(2), 4, 4*sqrt(2), 8 (num_of_tilts = 7)
	// t = 1, sqrt(2), 2, 2*sqrt(2), 4, 4*sqrt(2) (num_of_tilts = 6, paper's param)
	int num_of_tilts = 6;
	int verb = 0;
	// Define the SIFT parameters
	siftPar siftparameters;
	default_sift_parameters(siftparameters);
	
	std::vector<std::vector<keypointslist>> tkeypoints;
	
	int num_keys = 0;
	
	
	std::cout << "Computing keypoints on a image..." << std::endl;
	std::time_t tstart, tend;
	tstart = time(0);
	
	num_keys = compute_asift_keypoints(ipixels_zoom, wS, hS, num_of_tilts, verb, tkeypoints, siftparameters);
	
	tend = time(0);
	std::cout << "Keypoints computation accomplished in " << std::difftime(tend, tstart) << " seconds." << std::endl;
	std::cout << "Detected keypoints: " << num_keys << std::endl;
	
	///// Convert the result of compute_asift_keypoints into a form of keypoints of openCV
	if(!keypoints.empty())
		keypoints.clear();
	if(!prmKeypoints.empty())
		prmKeypoints.clear();
	prmKeypoints = std::vector<KeyPoint_param>(num_keys);
	keypoints = std::vector<cv::KeyPoint>(num_keys);
	descriptors.create(num_keys, VecLength, cv::DataType<float>::type);
	int itKps = 0;
	for(int itTilt=0; itTilt<(int)tkeypoints.size(); itTilt++) {
		for(int itLgtd=0; itLgtd<(int)tkeypoints[itTilt].size(); itLgtd++) {
			for(int it=0; it<(int)tkeypoints[itTilt][itLgtd].size(); it++) {
				keypoints[itKps].pt.x = tkeypoints[itTilt][itLgtd][it].x;
				keypoints[itKps].pt.y = tkeypoints[itTilt][itLgtd][it].y;
				keypoints[itKps].octave = std::log2(tkeypoints[itTilt][itLgtd][it].octSize);
				keypoints[itKps].size = tkeypoints[itTilt][itLgtd][it].scale;
				keypoints[itKps].angle = tkeypoints[itTilt][itLgtd][it].angle;
				prmKeypoints[itKps].tilt = tkeypoints[itTilt][itLgtd][it].tilt;
				prmKeypoints[itKps].angle = tkeypoints[itTilt][itLgtd][it].longitude;
				for(int itv=0; itv<VecLength; itv++)
					descriptors.at<float>(itKps, itv) = tkeypoints[itTilt][itLgtd][it].vec[itv];
				itKps++;
			}
		}
	}
	std::cout << "Done" << std::endl;
}

void morel_asift(cv::Mat img, std::vector<MatchingModel::View> &views,
					  std::vector<cv::KeyPoint> &keypoints,
					  std::vector<KeyPoint_param> &prmKeypoints, cv::Mat &descriptors) {
	
	cv::Mat src;
	
	//////////////////////////////////////////////// Input
	// Read image1
	int w, h;
	w = img.cols;
	h = img.rows;
	cvtColor(img, src, CV_BGR2GRAY);
	std::vector<float> ipixels(w * h);
	
	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols; j++)
			ipixels[i*img.cols + j] = (float)src.at<uchar>(i, j);
	
	src.release();
	
	///// Resize the images to area wS*hW in remaining the apsect-ratio
	///// Resize if the resize flag is not set or if the flag is set unequal to 0
	
	float zoom = 0;
	int wS = 0, hS = 0;
	std::vector<float> ipixels_zoom;
	
	ipixels_zoom.resize(w*h);
	ipixels_zoom = ipixels;
	wS = w;
	hS = h;
	zoom = 1;
	
	///// Compute ASIFT keypoints
	// Define the SIFT parameters
	siftPar siftparameters;
	default_sift_parameters(siftparameters);
	
	std::vector<std::vector<keypointslist>> tkeypoints;
	
	int num_keys = 0;
	
	std::cout << "Computing keypoints on a image..." << std::endl;
	std::time_t tstart, tend;
	tstart = time(0);
	
	//num_keys = compute_asift_keypoints(ipixels_zoom, wS, hS, num_of_tilts, verb, tkeypoints, siftparameters);
	
	
	
	
	
	
	
	
	
	
	
	tend = time(0);
	std::cout << "Keypoints computation accomplished in " << std::difftime(tend, tstart) << " seconds." << std::endl;
	std::cout << "Detected keypoints: " << num_keys << std::endl;
	
	///// Convert the result of compute_asift_keypoints into a form of keypoints of openCV
	if(!keypoints.empty())
		keypoints.clear();
	if(!prmKeypoints.empty())
		prmKeypoints.clear();
	prmKeypoints = std::vector<KeyPoint_param>(num_keys);
	keypoints = std::vector<cv::KeyPoint>(num_keys);
	descriptors.create(num_keys, VecLength, cv::DataType<float>::type);
	int itKps = 0;
	for(int itTilt=0; itTilt<(int)tkeypoints.size(); itTilt++) {
		for(int itLgtd=0; itLgtd<(int)tkeypoints[itTilt].size(); itLgtd++) {
			for(int it=0; it<(int)tkeypoints[itTilt][itLgtd].size(); it++) {
				keypoints[itKps].pt.x = tkeypoints[itTilt][itLgtd][it].x;
				keypoints[itKps].pt.y = tkeypoints[itTilt][itLgtd][it].y;
				keypoints[itKps].octave = std::log2(tkeypoints[itTilt][itLgtd][it].octSize);
				keypoints[itKps].size = tkeypoints[itTilt][itLgtd][it].scale;
				keypoints[itKps].angle = tkeypoints[itTilt][itLgtd][it].angle;
				prmKeypoints[itKps].tilt = tkeypoints[itTilt][itLgtd][it].tilt;
				prmKeypoints[itKps].angle = tkeypoints[itTilt][itLgtd][it].longitude;
				for(int itv=0; itv<VecLength; itv++)
					descriptors.at<float>(itKps, itv) = tkeypoints[itTilt][itLgtd][it].vec[itv];
				itKps++;
			}
		}
	}
	std::cout << "Done" << std::endl;
}













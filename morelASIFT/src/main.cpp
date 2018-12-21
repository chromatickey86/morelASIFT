//
//  main.cpp
//  ASIFT_MatchingProbModel
//
//  Created by Yum Joohyuk on 2018. 5. 27..
//  Copyright © 2018년 Yum Joohyuk. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include "morel_asift.h"
#include "evaluation_kit.h"
#include "match_util.h"
#include "MatchingModel.hpp"

#define SW_VER "hessASIFT_v1.0"

int parse(int argc, const char * argv[]) {
	printf("%s <function> <Input_image> <output> <views.xml>\n", argv[0]);
	printf("<function>\n");
	printf("0 : compute origianl ASIFT features proposed by Morel et al.\n");
	printf("%s 0 <input> <output>\n\n", argv[0]);
	
	printf("1 : compute ASIFT features concerning input views.xml\n");
	printf("%s 0 <input> <output> <views.xml>\n\n", argv[0]);
	return 0;
}

int main_test(int argc, const char * argv[]) {
	
	// 1. Load source image
	std::string nmSrc1 = "data/adam1.png";
	std::string nmSrc2 = "data/adam2.png";
	cv::Mat image1 = cv::imread(nmSrc1);
	if(image1.empty()) {
		std::cout << "main() - [ERR] '" << nmSrc1 << "' does not exist\n" << std::endl;
		system("pause");
		exit(1);
	}
	cv::Mat image2 = cv::imread(nmSrc2);
	if(image2.empty()) {
		std::cout << "main() - [ERR] '" << nmSrc2 << "' does not exist\n" << std::endl;
		system("pause");
		exit(1);
	}
	
	// 2. Call ASIFT function
	std::vector<cv::KeyPoint> keypoints1, keypoints2;
	std::vector<KeyPoint_param> prmKeypoints1, prmKeypoints2;
	cv::Mat descriptors1, descriptors2;
	morel_asift(image1, keypoints1, prmKeypoints1, descriptors1);
	morel_asift(image2, keypoints2, prmKeypoints2, descriptors2);
	
	// 3. Matching
	find_transform(image1, keypoints1, descriptors1, image2, keypoints2, descriptors2);
	
	return 0;
}

int maintask_morelASIFT(int argc, const char * argv[]) {
	
	printf("function 0 [morelASIFT]\n");
	
	const int req_argc = 4; //program + function + input + output
	if(argc != req_argc) {
		printf("maintask_morelASIFT() - [Parsing error] # of argument(%d) is not satisfied\n", argc);
		parse(argc, argv);
	}
	
	// 1. Load source image
	cv::Mat src = cv::imread(argv[2]);
	if(src.empty()) {
		printf("maintask_morelASIFT() - [File i/o error] %s does not exist\n", argv[1]);
		std::system("pause");
		exit(1);
	}
	printf("src\t:%s\\n", argv[2]);
	
	// 2. Call ASIFT function
	std::vector<cv::KeyPoint> keypoints;
	std::vector<KeyPoint_param> prmKeypoints;
	cv::Mat descriptors;
	morel_asift(src, keypoints, prmKeypoints, descriptors);
	
	// 3. output
	char* ofn = new char[100];
	strcpy(ofn, argv[3]);
	evaluation_kit eval_kit(ofn, evaluation_kit::DOT_ASIFT, evaluation_kit::WRITE);
	eval_kit(keypoints, prmKeypoints, descriptors);
	delete[] ofn;
	
	return 0;
}

int maintask_customASIFT(int argc, const char * argv[]) {
	printf("function 0 [customASIFT]\n");
	
	const int req_argc = 5; // program + function + input_image + output + views.xml
	if(argc != req_argc) {
		printf("maintask_customASIFT() - [Parsing error] # of argument(%d) is not satisfied\n", argc);
		parse(argc, argv);
	}
	
	// 1. Load source image
	cv::Mat src = cv::imread(argv[2]);
	if(src.empty()) {
		printf("maintask_customASIFT() - [File i/o error] %s does not exist\n", argv[1]);
		std::system("pause");
		exit(1);
	}
	printf("src\t:%s\\n", argv[2]);
	
	// 2. Call ASIFT function
	std::vector<cv::KeyPoint> keypoints;
	std::vector<KeyPoint_param> prmKeypoints;
	cv::Mat descriptors;
	
	std::vector<MatchingModel::View> views;
	std::string fnm(argv[4]);
	MatchingModel::readViewsFromXml(fnm, views);
	if(views.empty()) {
		printf("main() - [Content error] %s has nothing\n", argv[4]);
		std::system("pause");
		exit(1);
	}
	morel_asift(src, views, keypoints, prmKeypoints, descriptors);
	
	// 3. output
	char* ofn = new char[100];
	strcpy(ofn, argv[3]);
	evaluation_kit eval_kit(ofn, evaluation_kit::DOT_ASIFT, evaluation_kit::WRITE);
	eval_kit(keypoints, prmKeypoints, descriptors);
	delete[] ofn;
	
	return 0;
}

int main(int argc, const char * argv[]) {
	
	printf("===========================================\n");
	printf("%s\n", SW_VER);
	
	const int req_argc = 2; //program + function
	if(argc < req_argc)
		return parse(argc, argv);
	int function = std::atoi(argv[1]);
	if(function == 0)
		return maintask_morelASIFT(argc, argv);
	else if(function == 1)
		return main_test(argc, argv);
	else {
		printf("main() - [Parsing error] # of argument(%d) is not satisfied\n", argc);
		parse(argc, argv);
	}
	
	return 0;
}








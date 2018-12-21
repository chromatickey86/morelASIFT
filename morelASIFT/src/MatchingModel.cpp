//
//  MatchingModel.cpp
//  MatchingModel
//
//  Created by JoohyukYum on 2016. 7. 18..
//  Copyright © 2016년 JoohyukYum. All rights reserved.
//

#include <vector>
#include "MatchingModel.hpp"

//Build//
//0. Display
void MatchingModel::display(const std::string &nmWnd, const std::vector<InfoMatch> &infoMatches) {

	const cv::Scalar white = cv::Scalar(255.0, 255.0, 255.0);
//	const cv::Scalar black = cv::Scalar(0.0, 0.0, 0.0);
//	const cv::Scalar gray = cv::Scalar(128.0, 128.0, 128.0);
//	const cv::Scalar red = cv::Scalar(0.0, 0.0, 255.0);
//	const cv::Scalar green = cv::Scalar(0.0, 255.0, 0.0);
	const int width = 600;
	const int hfWidth = width/2;
	const double unit = hfWidth/3.5;
	cv::Mat frameBuff = cv::Mat::zeros(width, width, CV_8UC3);
	
	//display coordinate system
	cv::arrowedLine(frameBuff, cv::Point(0,hfWidth), cv::Point(width-1, hfWidth), white, 1.0, CV_AA, 0, 0.02);
	cv::arrowedLine(frameBuff, cv::Point(hfWidth,width-1), cv::Point(hfWidth, 0), white, 1.0, CV_AA, 0, 0.02);
	for(double col=0.0; col<width; col+=(0.1*unit)) {
		cv::line(frameBuff, cv::Point(col, hfWidth-2), cv::Point(col, hfWidth+2), white, 2.0, CV_AA, 0);
		cv::line(frameBuff, cv::Point(hfWidth-2, col), cv::Point(hfWidth+2, col), white, 2.0, CV_AA, 0);
	}
	
	//display matched views
	for(int it=0; it<(int)infoMatches.size(); it++) {
		double ms = infoMatches[it].matchingScore;
		double sr = infoMatches[it].vwTrain.sr*unit;
		double tau = infoMatches[it].vwTrain.tau*unit;
		if(ms > 0)
			cv::circle(frameBuff, cv::Point(hfWidth+tau, hfWidth-sr), 2.0, colorCoding(ms, 0.0, 100.0), 1.5);
	}
	cv::imshow(nmWnd, frameBuff);
	cv::imwrite("capture.png", frameBuff);
	
	//DEBUG::parameter K
//	double idx=-2.0;
//	std::vector<double> histK;
//	calc_K(infoMatches, histK);
//	for(int it=0; it<(int)histK.size()/2; it++) {
//		printf("[%.1f](%.1f)%.3f\n", idx, histK[2*it], histK[2*it+1]);
//		idx += 0.1;
//	}
}

void MatchingModel::display(const std::string &nmWnd, const std::vector<View> &views) {

	const cv::Scalar white = cv::Scalar(255.0, 255.0, 255.0);
//	const cv::Scalar red = cv::Scalar(0.0, 0.0, 255.0);
	const cv::Scalar green = cv::Scalar(0.0, 255.0, 0.0);
	const int width = 800;
	const int hfWidth = width/2;
	const double unit = hfWidth/3.5;
	cv::Mat frameBuff = cv::Mat::zeros(width, width, CV_8UC3);
	
	//display coordinate system
	cv::arrowedLine(frameBuff, cv::Point(0,hfWidth), cv::Point(width-1, hfWidth), white, 1.0, CV_AA, 0, 0.02);
	cv::arrowedLine(frameBuff, cv::Point(hfWidth,width-1), cv::Point(hfWidth, 0), white, 1.0, CV_AA, 0, 0.02);
	for(double col=0.0; col<width; col+=(0.1*unit)) {
		cv::line(frameBuff, cv::Point(col, hfWidth-2), cv::Point(col, hfWidth+2), white, 2.0, CV_AA, 0);
		cv::line(frameBuff, cv::Point(hfWidth-2, col), cv::Point(hfWidth+2, col), white, 2.0, CV_AA, 0);
	}
	
	//display matched views
	for(int it=0; it<(int)views.size(); it++) {
		char cstr[4];
		sprintf(cstr, "%d", it);
		std::string str = cstr;
		double sr = views[it].sr*unit;
		double tau = views[it].tau*unit;
		cv::Point position = cv::Point(hfWidth+tau, hfWidth-sr);
		cv::putText(frameBuff, str, position, CV_FONT_HERSHEY_COMPLEX_SMALL, 0.8, white, 1, CV_AA);
		cv::circle(frameBuff, position, 2.0, green, 1.5);
	}
	cv::imshow(nmWnd, frameBuff);
}

void MatchingModel::displayHemisphere(const std::string &nmWnd, const int lat, const int lon) {
	double latitude = lat/180.0*CV_PI;
	double longitude = lon/180.0*CV_PI;
	const cv::Scalar white = cv::Scalar(255.0, 255.0, 255.0);
	const cv::Scalar green = cv::Scalar(0.0, 255.0, 0.0);
	const int rngLat = 15;
	const int width = 300;
	const int hfWidth = width/2;
	cv::Mat frameBuff = cv::Mat::zeros(width, width, CV_8UC3);
	
	//display coordinate system
	cv::arrowedLine(frameBuff, cv::Point(0,hfWidth), cv::Point(width-1, hfWidth), white, 1.0, CV_AA, 0, 0.02);
	cv::arrowedLine(frameBuff, cv::Point(hfWidth,width-1), cv::Point(hfWidth, 0), white, 1.0, CV_AA, 0, 0.02);
	for(int itLat=0; itLat<90; itLat+=rngLat)
		if(itLat!=0) cv::circle(frameBuff, cv::Point(hfWidth, hfWidth), hfWidth*std::sin((double)itLat/180.0*CV_PI), white, 1.0);
	char cstr1[5], cstr2[5];
	sprintf(cstr1, "%d", lat);
	sprintf(cstr2, "%d", lon);
	std::string str = (std::string)"lat:" + cstr1 + " lon:" + cstr2;
	cv::putText(frameBuff, str, cv::Point(hfWidth, 30), CV_FONT_HERSHEY_COMPLEX_SMALL, 0.8, white, 1, CV_AA);
	
	//display matched views
	double radius = std::sin(latitude) * hfWidth;
	cv::Point position = cv::Point(hfWidth, hfWidth) + cv::Point(std::cos(longitude)*radius, std::sin(longitude)*radius);
	cv::circle(frameBuff, position, 3.0, green, 2.0);
	
	//display
	cv::imshow(nmWnd, frameBuff);
}

cv::Scalar MatchingModel::colorCoding(const double &val, const double &min, const double &max) {
    const double MAX = 256.0;
    cv::Scalar ret(0.0, 0.0, 0.0);
    if(max <= min || val>max || val<min) {
        std::cerr << "::MatchingModel::clrCoding:: Scope of val is wrong" <<std::endl;
        return ret;
    }
    double V = val / (max - min) * 255.0 *  3.0;
    if(V<MAX) {
        ret[2] = V;
    }
    else if(V<MAX*2.0) {
        ret[2] = MAX*2.0 - V;
        ret[1] = V - MAX;
    }
    else {
        ret[1] = MAX*3.0 - V;
        ret[0] = V - MAX*2.0;
    }
    return ret;
}

//1. ParameterCalcuation
double MatchingModel::calc_K(const std::vector<InfoMatch> &infoMatches) {
    std::vector<double> history;
    return calc_K(infoMatches, history);
}

double MatchingModel::calc_K(const std::vector<InfoMatch> &infoMatches, std::vector<double> &history) {
	const double rngTau = 0.1;
	std::vector<std::vector<double>> histogramMS, histogramSR;
	for(double itTau=-2.0; itTau<=2.1; itTau+=rngTau) {
		histogramMS.push_back(std::vector<double>());
		histogramSR.push_back(std::vector<double>());
	}
	for(int it=0; it<(int)infoMatches.size(); it++) {
		int idx = (infoMatches[it].vwTrain.g+2.0)*10.0 + 0.5;
		histogramMS[idx].push_back(infoMatches[it].matchingScore);
		histogramSR[idx].push_back(infoMatches[it].vwTrain.sr);
		 
	}
	for(int it=0; it<(int)histogramMS.size(); it++) {
		double max=0.0;
		double idxMax=0;
		for(int iit=0; iit<(int)histogramMS[it].size(); iit++) {
			if(histogramMS[it][iit]>max) {
				max = histogramMS[it][iit];
//				idxMax = iit;
				idxMax = histogramSR[it][iit];
			}
		}
		history.push_back(max);
		history.push_back(idxMax);
	}
	return 0.0;
}

void MatchingModel::calcModifiedViews(std::vector<View> &views, const int lat, const int lon) {
	double latitude = lat/180.0*CV_PI;
	double longitude = lon/180.0*CV_PI;
	double cos, sin, t;
	double A, B, C, D, SX, SY, G;
	for(int it=0; it<(int)views.size(); it++) {
		t = 1.0/std::cos(latitude);
		cos = std::cos(longitude);
		sin = std::sin(longitude);
		A = views[it].sx*cos + views[it].g*views[it].sy*sin/t;
		B = -views[it].sx*sin + views[it].g*views[it].sy*cos/t;
		C = views[it].sy*sin/t;
		D = views[it].sy*cos/t;
		SX = std::sqrt(A*A + C*C);
		SY = (A*D - B*C) / SX;
		G = (A*B + C*D)/(A*D - B*C);
		views[it] = View(SX, SY, G);
	}
}


//2. ViewSelection
int MatchingModel::buildOriginalViews(std::vector<View> &views) {
	const int nTilts(6);
	const float minTilt(1.0);
	const float srTilt(std::sqrt(2.0));
	return buildOriginalViews(views, nTilts, minTilt, srTilt);
}

int MatchingModel::buildOriginalViews(std::vector<View> &views, int nTilts, float minTilt, double srTilt) {
	int nAngl;
	float tilt, angle, srAngl;
	views.clear();
	for (int it_tilt(0); it_tilt < nTilts; it_tilt++) {
		tilt = (float)minTilt * cv::pow(srTilt, it_tilt);
		if (it_tilt == 0) {
			views.push_back(MatchingModel::View(std::acos(1.0/tilt), 0.0));
		}
		else {
			nAngl = std::round(10.0 * tilt / 2.0);
			nAngl = (nAngl % 2 == 1 ? nAngl + 1 : nAngl) / 2;
			srAngl = 72.0 / tilt;
			if (it_tilt == 6 - 1) nAngl++;	//Is it right??
			for (int it_angl(0); it_angl < nAngl; it_angl++) {
				angle = srAngl * it_angl;
				views.push_back(MatchingModel::View(std::acos(1.0/tilt), angle/180.0*CV_PI));
			}
		}
	}
//	printf("===========================\n");
//	double area = 0.0;
//	for (int it(0); it < (int)views.size(); it++) {
//		printf("[%d] (%.3f, %.3f)\n", it, views[it][0]*180.0/CV_PI, views[it][1]*180.0/CV_PI);
//		area += std::cos(views[it][0]);
//	}
//	printf("area(%.3f)\n", area);
//	printf("===========================\n");
	return (int)views.size();
}

int MatchingModel::buildUniformDistrViews(std::vector<View> &views) {
	const int width = 3*2 + 1;
	const double rng = 0.375;
//	const double sigma = 0.562;
	const double sigma = 0.566;

	const double init_sr = (double)(-width/2) * rng;
	const double init_tau = (double)(-width/2) * rng;
	double sr, sxpsy, sx, sy, g, distance;
	double area=0.0;
	for(int it_tau=0; it_tau<width; it_tau++) {
		double tau = init_tau + rng*it_tau;
		for(int it_sr=0; it_sr<width; it_sr++) {
			sr = init_sr + rng*it_sr;
			distance = std::sqrt(sr*sr + tau*tau);
			sxpsy = std::pow(10.0, sr/W1);
			double unit_area = std::exp(-(sr*sr + tau*tau)/(2.0*sigma*sigma));
			sx = std::sqrt(unit_area * sxpsy);
			sy = sx / sxpsy;
			area += sx*sy;
			g = std::tan(tau)*W2;
			views.push_back(View(sx, sy, g));
			printf("(%.3f, %.3f) -> (%.3f, %.3f, %.3f)\n", sxpsy, tau, sx, sy, g);
		}
	}
	printf("total area : %.3f\n", area);
	
	return (int)views.size();
}

int MatchingModel::buildProposedViews(std::vector<View> &views) {
	
	//0. initial
	views.push_back(View(1.0, 1.0, 0.0));
	
	//1. Select viewpoints
	/*CONFIG*/ const double sigma = 0.566;
	/*CONFIG*/ const double rngRadius = 3.5;
	/*CONFIG*/ const int stRadius = 100;
	const double intvlRadius = (double)stRadius/rngRadius;
	/*CONFIG*/ const int stAngl = 360;
	/*CONFIG*/ const double intvlAngl = CV_PI / stAngl;
	double radius, angl, sr, tau, unit_area, sxpsy, sx, sy, g;
	for(int itRadius=0; itRadius<stRadius; itRadius++) {
		radius = itRadius * intvlRadius;
		for(int itAngl=0; itAngl<stAngl; itAngl++) {
			angl = itAngl * intvlAngl;
			sr = radius * std::cos(angl);
			tau = radius * std::sin(angl);
			unit_area = std::exp(-(sr*sr + tau*tau)/(2.0*sigma*sigma));
			sxpsy = std::pow(10.0, sr/W1);
			sx = std::sqrt(unit_area * sxpsy);
			sy = sx / sxpsy;
			g = std::tan(tau)*W2;
			View tview(sx, sy, g);
			for(int itViews=0; itViews<(int)views.size(); itViews++) {
				
			}
		}
	}
	
	
	return (int)views.size();
}


//Data//
//1. TestImageGen
void MatchingModel::imageSynthesis(const cv::Mat src, const int transType) {
	
	if(transType==AFFINE_TRANSFORM) {
		std::string nmFolder = "AffineImage/";
		std::string nmBatchFile = "script.bat";
		FILE *fscript = fopen(nmBatchFile.c_str(), "w");

//		for(int it_g=-10; it_g<=10; it_g++) {
		for(int it_g=-20; it_g<=20; it_g++) {
			for(int it_sy=1; it_sy<=10; it_sy++) {
				for(int it_sx=1; it_sx<=10; it_sx++) {
					double sx = (double)it_sx/10.0;
					double sy = (double)it_sy/10.0;
					double g = (double)it_g/10.0;
					
					cv::Mat tfImg = affineTransform(src, sx, sy, g);
					cv::imshow("display", tfImg);
					
					char numSx[7], numSy[7], numG[7];
					sprintf(numSx, "%.2f", sx);
					sprintf(numSy, "%.2f", sy);
					sprintf(numG, "%.2f", g);
					std::string nmInstance = (std::string)"sx"+numSx+(std::string)"sy"+numSy+(std::string)"g"+numG;
					std::string nmImgFile = nmFolder + nmInstance + ".png";
					
					//File writing
					cv::imwrite(nmImgFile.c_str(), tfImg);
					
					MatchingModel::View tgView(sx, sy, g);
					for(int it_base_g=0; it_base_g<=10; it_base_g+=5) {
						for(int it_base_sy=5; it_base_sy<=10; it_base_sy+=5) {
							for(int it_base_sx=5; it_base_sx<=10; it_base_sx+=5) {
								double bs_sx = (double)it_base_sx/10.0;
								double bs_sy = (double)it_base_sy/10.0;
								double bs_g = (double)it_base_g/10.0;
								MatchingModel::View bsView(bs_sx, bs_sy, bs_g);
								cv::Mat H = getHomography(src, bsView, tgView);
								sprintf(numSx, "%.2f", bs_sx);
								sprintf(numSy, "%.2f", bs_sy);
								sprintf(numG, "%.2f", bs_g);
								std::string nmInstance_detail = (std::string)"sx"+numSx+(std::string)"sy"+numSy+(std::string)"g"+numG;
								std::string nmHomoFile = nmFolder + nmInstance_detail + nmInstance + ".csv";
								wrtMat2csv(nmHomoFile, H);
							}
						}
					}
					
					//BatchFile writing
					std::string cmd = "ASIFT_REF_SW.exe " + nmImgFile + " result/" + nmInstance + ".sift\n";
					fprintf(fscript, "%s", cmd.c_str());
					cv::waitKey(10);
				}
			}
		}
		fclose(fscript);
	}
	else if(transType==PERSPECTIVE_TRANSFORM) {
		std::string nmFolder = "PerspectImage/";
		std::string nmBatchFile = "script.bat";
		FILE *fscript = fopen(nmBatchFile.c_str(), "w");
		const double PI_180 = CV_PI/180.0;
		const double distance = 10.0;
		const double zoom = 6.0;
		const int vertTilt = 1;
		const int stridLat = 5;
		const int stridLon = 10;
		
		cv::Mat blurImg, tfImg;
		for (int itLat=0; itLat<=80; itLat+=stridLat) {
			for(int itLon=0; itLon<=180; itLon+=stridLon) {
				double latitude = itLat * PI_180;
				double longitude = itLon * PI_180;
				double t = 1.0/std::cos(latitude);
				cv::GaussianBlur(src, blurImg, cv::Size(), t/3.0);
				tfImg = perspectiveTransform(distance, zoom, blurImg, latitude, longitude, vertTilt);
				cv::Mat H1 = getHomographyFromFront(distance, zoom, blurImg, latitude, longitude, vertTilt);
				cv::Mat H2 = getHomographyFromTilted(distance, zoom, blurImg, latitude, longitude, vertTilt);
				cv::imshow("display", tfImg);
				
				char numLat[5], numLon[5];
				sprintf(numLat, "%d", itLat);
				sprintf(numLon, "%d", itLon);
				std::string instanceName = (std::string)"lat" + numLat + "lon" + numLon;
				std::string fnImage = nmFolder + instanceName + ".png";
				std::string fnHomo1 =  nmFolder + instanceName + ".csv";
				std::string fnHomo2 =  nmFolder + instanceName + "_.csv";
				cv::imwrite(fnImage, tfImg);
				wrtMat2csv(fnHomo1, H1);
				wrtMat2csv(fnHomo2, H2);
				std::string cmd = "ASIFT_REF_SW.exe " + fnImage + " result/" + instanceName + ".asift\n";
				fprintf(fscript, "%s", cmd.c_str());
				cv::waitKey(10);
			}
		}
		fclose(fscript);
	}
}

//2. File I/O
int MatchingModel::wrtMat2csv(const std::string &fnm, const cv::Mat mat) {
	if(mat.empty())
		return -1;
	FILE *fp = fopen(fnm.c_str(), "w");
	for(int itRows=0; itRows<mat.rows; itRows++) {
		for(int itCols=0; itCols<mat.cols; itCols++) {
			fprintf(fp, "%lf,", mat.at<double>(itRows, itCols));
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fclose(fp);
	return 1;
}

cv::Mat MatchingModel::readXml(const std::string &fnm) {
	cv::FileStorage fs;
	cv::Mat ret;
	if(fs.open(fnm, cv::FileStorage::READ)) {
		fs["mat"] >> ret;
		fs.release();
	}
	return ret;
}

int MatchingModel::wrtXml(const std::string &fnm, cv::Mat mat) {
	cv::FileStorage fs;
	if(fs.open(fnm, cv::FileStorage::WRITE)) {
		fs << "mat" << mat;
		fs.release();
		return 0;
	}
	return -1;
}

void MatchingModel::readInfoMatches(std::vector<InfoMatch> &infoMatches) {
	std::string nmFolder = "infoMatches_civilwar/sx1.00sy1.00g0.50";
//	std::string nmFolder = "infoMatches_sierra/sx0.50sy0.50g1.00";
	cv::Mat t;
	infoMatches.clear();
	for(int it_g=-20; it_g<=20; it_g++) {
		double g = (double)it_g/10.0;
		char num_g[5];
		sprintf(num_g, "%.2f", g);
		std::string instanceName = nmFolder+(std::string)"g"+num_g+"_infoMatches.xml";
		t = readXml(instanceName);
		if(!t.empty()) {
			for(int r=0; r<t.rows; r++) {
				View vwQuery(t.at<float>(r,0), t.at<float>(r,1), t.at<float>(r,2));
				View vwTrain(t.at<float>(r,3), t.at<float>(r,4), t.at<float>(r,5));
				InfoMatch info(1, vwQuery, vwTrain, t.at<float>(r,6), t.at<float>(r,7), t.at<float>(r,8),
									t.at<float>(r,9), t.at<float>(r,10), t.at<float>(r,11));
				infoMatches.push_back(info);
			}
		}
		else
			std::cout << "ERROR::" << instanceName << " doesn't exist\n";
	}
}

void MatchingModel::wrtInfoMatches2csv(const std::string &fnm, const std::vector<InfoMatch> &infoMatches) {
	cv::Mat t((int)infoMatches.size(), /*#ofElements*/12, CV_64F);
//	double repeatability, matchingScore;
//	int correspondance, correctMatches, nQueryKps, nTrainKps;
	for(int i=0; i<(int)infoMatches.size(); i++) {
		t.at<double>(i, 0) = infoMatches[i].vwQuery.sx;
		t.at<double>(i, 1) = infoMatches[i].vwQuery.sy;
		t.at<double>(i, 2) = infoMatches[i].vwQuery.g;
		t.at<double>(i, 3) = infoMatches[i].vwTrain.sx;
		t.at<double>(i, 4) = infoMatches[i].vwTrain.sy;
		t.at<double>(i, 5) = infoMatches[i].vwTrain.g;
		t.at<double>(i, 6) = infoMatches[i].repeatability;
		t.at<double>(i, 7) = infoMatches[i].matchingScore;
		t.at<double>(i, 8) = infoMatches[i].correspondance;
		t.at<double>(i, 9) = infoMatches[i].correctMatches;
		t.at<double>(i, 10) = infoMatches[i].nQueryKps;
		t.at<double>(i, 11) = infoMatches[i].nTrainKps;
	}
	wrtMat2csv(fnm, t);
}

void MatchingModel::wrtViews2xml(const std::string & fnm, const std::vector<View> &views) {
	cv::Mat matViews((int)views.size(), 3, CV_32F);
	for(int it=0; it<(int)views.size(); it++) {
		matViews.at<float>(it, 0) = views[it].sx;
		matViews.at<float>(it, 1) = views[it].sy;
		matViews.at<float>(it, 2) = views[it].g;
	}
	wrtXml(fnm, matViews);
}

void MatchingModel::readViewsFromXml(const std::string & fnm, std::vector<View> &views) {
	cv::Mat matViews = readXml(fnm);
	cv::Mat txml = MatchingModel::readXml(fnm);
	for(int it=0; it<(int)txml.rows; it++) {
		MatchingModel::View tview(matViews.at<double>(it, 0), matViews.at<double>(it, 1), matViews.at<double>(it, 2));
		views.push_back(tview);
	}
}

//Transform//
void MatchingModel::bound(int x, int y, double cos, double sin, int &xmin, int &xmax, int &ymin, int &ymax) {
	int rx = (int)floor(cos * (double)x + sin * (double)y);
	int ry = (int)floor(-sin * (double)x + cos * (double)y);
	if(rx < xmin) xmin = rx;
	if(rx > xmax) xmax = rx;
	if(ry < ymin) ymin = ry;
	if(ry > ymax) ymax = ry;
}

void MatchingModel::rotation(cv::Vec3f &vec, const double rx, const double ry, const double rz) {
	float cx = std::cos(rx), cy = std::cos(ry), cz = std::cos(rz);
	float sx = std::sin(rx), sy = std::sin(ry), sz = std::sin(rz);
	
	cv::Vec3f newVec;
	newVec[0] = cy*(-sz*vec[1] + cz*vec[0]) + sy*vec[2];
	newVec[1] = -sx*(cy*vec[2] + sy*(sz*vec[1] - cz*vec[0])) + cx*(cz*vec[1] + sz*vec[0]);
	newVec[2] = cx*(cy*vec[2] + sy*(sz*vec[1] - cz*vec[0])) + sx*(cz*vec[1] + sz*vec[0]);
	vec = newVec;
}

cv::Mat MatchingModel::perspectiveTransform(double distance, double zoom, const cv::Mat img, const double lat, const double lon, const int vertTilt) {
	
	const double lenFocal = img.cols/2;
	const int hfImgWidth = img.cols/2, hfImgHeight = img.rows/2;
	
	//Find the surface equation
	// A * x + B * y + C * c + D = 0
	// noraml = (A, B, C) D = -C*len;
	cv::Vec3f normal(0.0, 0.0, 1.0);
	if(vertTilt)
		rotation(normal, -lat, 0.0, lon);
	else
		rotation(normal, 0.0, lat, lon);
	double D = -normal[2]*distance;
	
	cv::Mat perImg = cv::Mat::zeros(img.rows, img.cols, img.type());
	cv::Vec3f center(0.0, 0.0, distance);
	for(int itRows=0; itRows<perImg.rows; itRows++) {
		for(int itCols=0; itCols<perImg.cols; itCols++) {
			double a = ((double)itCols-hfImgWidth) / lenFocal / zoom;
			double b = ((double)itRows-hfImgHeight) / lenFocal / zoom;
			double c = 1.0;
			double k = -D / (normal[0]*a + normal[1]*b + normal[2]*c);
			cv::Vec3f point(k*a, k*b, k*c);
			//(a,b,c) = point/k;
			point = point - center;
			//point+center
			if(vertTilt)
				rotation(point, lat, 0.0, 0.0);
			else
				rotation(point, 0.0, -lat, 0.0);
			rotation(point, 0.0, 0.0, -lon);
			point = point + center;
			//point-center
			double ixPixel = point[0] * hfImgWidth + hfImgWidth;
			double iyPixel = point[1] * hfImgWidth + hfImgHeight;
			if(ixPixel>=0 && ixPixel<img.cols-1 && iyPixel>=0 && iyPixel<img.rows-1) {
				int ipx_f(ixPixel), ipx_c(ixPixel == img.cols-1 ? ixPixel : ixPixel + 1);
				int ipy_f(iyPixel), ipy_c(iyPixel == img.rows-1 ? iyPixel : iyPixel + 1);
				float ipx_wc(ixPixel - ipx_f), ipx_wf(ipx_c - ixPixel);
				float ipy_wc(iyPixel - ipy_f), ipy_wf(ipy_c - iyPixel);
				perImg.at<cv::Vec3b>(itRows, itCols) =
				img.at<cv::Vec3b>(ipy_f, ipx_f)*ipx_wf*ipy_wf +
				img.at<cv::Vec3b>(ipy_f, ipx_c)*ipx_wc*ipy_wf +
				img.at<cv::Vec3b>(ipy_c, ipx_f)*ipx_wf*ipy_wc +
				img.at<cv::Vec3b>(ipy_c, ipx_c)*ipx_wc*ipy_wc;
			}
		}
	}
	return perImg;
}

cv::Mat MatchingModel::getHomographyFromFront(double distance, double zoom, const cv::Mat img, const double lat, const double lon, const int vertTilt) {
	
	const double lenFocal = img.cols/2;
	const int hfImgWidth = img.cols/2, hfImgHeight = img.rows/2;
	
	//Find the surface equation
	// A * x + B * y + C * c + D = 0
	// noraml = (A, B, C) D = -C*len;
	cv::Vec3f normal(0.0, 0.0, 1.0);
	cv::Vec3f center(0.0, 0.0, distance);
	double D = -normal[2]*distance;
	
	std::vector<cv::Vec2f> corners, originCorners, tfCorners;
	corners.push_back(cv::Vec2f(0, 0));
	corners.push_back(cv::Vec2f(img.cols-1, 0));
	corners.push_back(cv::Vec2f(img.cols-1, img.rows-1));
	corners.push_back(cv::Vec2f(0, img.rows-1));
	
	for(int it=0; it<(int)corners.size(); it++) {
		cv::Vec3f point((corners[it][0] - hfImgWidth)/hfImgWidth, (corners[it][1] - hfImgHeight)/hfImgWidth, distance);
		double a = point[0]/point[2], b = point[1]/point[2];
		originCorners.push_back(cv::Vec2f(a*lenFocal*zoom+hfImgWidth, b*lenFocal*zoom+hfImgHeight));
	}
	
	normal = cv::Vec3f(0.0, 0.0, 1.0);
	if(vertTilt)
		MatchingModel::rotation(normal, -lat, 0.0, lon);
	else
		MatchingModel::rotation(normal, 0.0, lat, lon);
	D = -normal[2]*distance;
	
	for(int it=0; it<(int)corners.size(); it++) {
		cv::Vec3f point((corners[it][0] - hfImgWidth)/hfImgWidth, (corners[it][1] - hfImgHeight)/hfImgWidth, 0);
		if(vertTilt)
			MatchingModel::rotation(point, -lat, 0.0, lon);
		else
			MatchingModel::rotation(point, 0.0, lat, lon);
		point = point + center;
		double a = point[0]/point[2], b = point[1]/point[2];
		tfCorners.push_back(cv::Vec2f(a*lenFocal*zoom+hfImgWidth, b*lenFocal*zoom+hfImgHeight));
	}
	cv::Mat H = cv::findHomography(originCorners, tfCorners, CV_RANSAC);
	return H;
}

cv::Mat MatchingModel::getHomographyFromTilted(double distance, double zoom, const cv::Mat img, const double lat, const double lon, const int vertTilt) {
	
	const double lenFocal = img.cols/2;
	const int hfImgWidth = img.cols/2, hfImgHeight = img.rows/2;
	
	//Find the surface equation
	// A * x + B * y + C * c + D = 0
	// noraml = (A, B, C) D = -C*len;
	cv::Vec3f normal(0.0, 0.0, 1.0);
	cv::Vec3f center(0.0, 0.0, distance);
	if(vertTilt)
		MatchingModel::rotation(normal, -lat, 0.0, 0.0);
	else
		MatchingModel::rotation(normal, 0.0, lat, 0.0);
	double D = -normal[2]*distance;
	
	std::vector<cv::Vec2f> corners, originCorners, tfCorners;
	corners.push_back(cv::Vec2f(0, 0));
	corners.push_back(cv::Vec2f(img.cols-1, 0));
	corners.push_back(cv::Vec2f(img.cols-1, img.rows-1));
	corners.push_back(cv::Vec2f(0, img.rows-1));
	
	for(int it=0; it<(int)corners.size(); it++) {
		cv::Vec3f point((corners[it][0] - hfImgWidth)/hfImgWidth, (corners[it][1] - hfImgHeight)/hfImgWidth, 0);
		if(vertTilt)
			MatchingModel::rotation(point, -lat, 0.0, 0.0);
		else
			MatchingModel::rotation(point, 0.0, lat, 0.0);
		point = point + center;
		double a = point[0]/point[2], b = point[1]/point[2];
		originCorners.push_back(cv::Vec2f(a*lenFocal*zoom+hfImgWidth, b*lenFocal*zoom+hfImgHeight));
	}
	
	normal = cv::Vec3f(0.0, 0.0, 1.0);
	if(vertTilt)
		MatchingModel::rotation(normal, -lat, 0.0, lon);
	else
		MatchingModel::rotation(normal, 0.0, lat, lon);
	D = -normal[2]*distance;
	
	for(int it=0; it<(int)corners.size(); it++) {
		cv::Vec3f point((corners[it][0] - hfImgWidth)/hfImgWidth, (corners[it][1] - hfImgHeight)/hfImgWidth, 0);
		if(vertTilt)
			MatchingModel::rotation(point, -lat, 0.0, lon);
		else
			MatchingModel::rotation(point, 0.0, lat, lon);
		point = point + center;
		double a = point[0]/point[2], b = point[1]/point[2];
		tfCorners.push_back(cv::Vec2f(a*lenFocal*zoom+hfImgWidth, b*lenFocal*zoom+hfImgHeight));
	}
	
	cv::Mat H = cv::findHomography(originCorners, tfCorners, CV_RANSAC);
	return H;
}

cv::Mat MatchingModel::affineTransform(const cv::Mat img, const double lat, const double lon, const int vertTilt) {
	
	int xmin, xmax, ymin, ymax;
	int width=img.cols, height=img.rows;
	
	//transformed image size calculation
	double tilt = 1.0/std::cos(lat);
	double cos = std::cos(-lon);
	double sin = std::sin(-lon);
	
	xmin = xmax = ymin = ymax = 0;
	bound(width - 1, 0, cos, sin, xmin, xmax, ymin, ymax); // bound for x
	bound(0, height - 1, cos, sin, xmin, xmax, ymin, ymax); // bound for y
	bound(width - 1, height - 1, cos, sin, xmin, xmax, ymin, ymax); // bound for x, y
	width = xmax - xmin + 1;
	height = ymax - ymin + 1;
	
	cv::Mat affineImg = cv::Mat::zeros(height, width, img.type());
	double bg = 0.0;
	double dbx, dby;
	int xf, yf;
	double aff, afc, acf, acc;
	float wxc, wyc;
	int bound_x_valid_0, bound_x_valid_1;
	int bound_y_valid_0, bound_y_valid_1;
	for (int ty = ymin; ty <= ymax; ty++) {
		for (int tx = xmin; tx <= xmax; tx++) {
			//Inverse scaling with tilt & roattion with angle
			dbx = cos * (double) tx - sin * tilt * (double)ty;
			dby = sin * (double) tx + cos * tilt * (double)ty;
			xf = std::floor(dbx);
			yf = std::floor(dby);
			wxc = dbx - (double)xf;
			wyc = dby - (double)yf;
			
			bound_x_valid_0 = (xf >= 0 && xf < img.cols);
			bound_x_valid_1 = (xf + 1 >= 0 && xf + 1 < img.cols);
			bound_y_valid_0 = (yf >= 0 && yf < img.rows);
			bound_y_valid_1 = (yf + 1 >= 0 && yf + 1 < img.rows);
			
			aff = (bound_x_valid_0 && bound_y_valid_0 ? img.at<double>(yf, xf) : bg);
			afc = (bound_x_valid_0 && bound_y_valid_1 ? img.at<double>(yf + 1, xf) : bg);
			acf = (bound_x_valid_1 && bound_y_valid_0 ? img.at<double>(yf, xf + 1) : bg);
			acc = (bound_x_valid_1 && bound_y_valid_1 ? img.at<double>(yf + 1, xf + 1) : bg);
			
			affineImg.at<double>(ty - ymin, tx - xmin) = (1.0 - wyc) * ((1.0 - wxc)*aff + wxc*acf) +
			wyc * ((1.0 - wxc)*afc + wxc*acc);
		}
	}
	return affineImg;
}

cv::Mat MatchingModel::affineTransform(const cv::Mat img, const double sx, const double sy, const double g) {
	
	const int width = sx*img.cols + std::abs(g)*sy*img.rows;
	const int stCol = std::min(0.0, g*sy*img.rows);
	const int edCol = sx*img.cols + std::max(0.0, g*sy*img.rows);
	const int height = sy*img.rows;
	const double GaussianSigma = std::log2(1.0/(sx*sy));
	int xf, yf;
	double dbx, dby, wxc, wyc;
	cv::Vec3b aff, afc, acf, acc, bg=0;
	int bound_x_valid_0, bound_x_valid_1;
	int bound_y_valid_0, bound_y_valid_1;
	cv::Mat blurredImg;
	if(sx*sy<1.0)
		cv::GaussianBlur(img, blurredImg, cv::Size(), GaussianSigma);
	else
		blurredImg = img.clone();
	cv::Mat affineImg = cv::Mat::zeros(height, width, img.type());
	for(int itRows=0; itRows<height; itRows++) {
		for(int itCols=stCol; itCols<edCol; itCols++) {
			//Inverse scaling with tilt & roattion with angle
			dbx = ((double)itCols - g*itRows)/sx;
			dby = (double)itRows/sy;
			xf = std::floor(dbx);
			yf = std::floor(dby);
			wxc = dbx - (double)xf;
			wyc = dby - (double)yf;
			
			bound_x_valid_0 = (xf >= 0) && (xf < img.cols);
			bound_x_valid_1 = (xf + 1 >= 0) && (xf + 1 < img.cols);
			bound_y_valid_0 = (yf >= 0) && (yf < img.rows);
			bound_y_valid_1 = (yf + 1 >= 0) && (yf + 1 < img.rows);
			
			aff = bound_x_valid_0 && bound_y_valid_0 ? blurredImg.at<cv::Vec3b>(yf, xf) : bg;
			afc = bound_x_valid_0 && bound_y_valid_1 ? blurredImg.at<cv::Vec3b>(yf + 1, xf) : bg;
			acf = bound_x_valid_1 && bound_y_valid_0 ? blurredImg.at<cv::Vec3b>(yf, xf + 1) : bg;
			acc = bound_x_valid_1 && bound_y_valid_1 ? blurredImg.at<cv::Vec3b>(yf + 1, xf + 1) : bg;
			
			affineImg.at<cv::Vec3b>(itRows, itCols-stCol) = (1.0-wyc) * ((1.0-wxc)*aff + wxc*acf) + wyc*((1.0-wxc)*afc + wxc*acc);
		}
	}
//	cv::imshow("test", affineImg);
//	cv::waitKey();
	return affineImg;
}

cv::Mat MatchingModel::getHomographyFromFront(const cv::Mat img, const double sx, const double sy, const double g) {
	
	const double col_ofs = std::min(0.0, g*sy*img.rows);
	std::vector<cv::Vec2f> originCorners, tfCorners;
	originCorners.push_back(cv::Vec2f(0, 0));
	originCorners.push_back(cv::Vec2f(img.cols-1, 0));
	originCorners.push_back(cv::Vec2f(img.cols-1, img.rows-1));
	originCorners.push_back(cv::Vec2f(0, img.rows-1));
	
	for(int itCnrs=0; itCnrs<(int)originCorners.size(); itCnrs++) {
		double y =sy*originCorners[itCnrs][1];
		double x = sx*originCorners[itCnrs][0] + g*y - col_ofs;
		tfCorners.push_back(cv::Vec2f(x, y));
	}
	cv::Mat H = cv::findHomography(originCorners, tfCorners, CV_RANSAC);
	return H;
}

cv::Mat MatchingModel::getHomography(const cv::Mat img, const View &v1, const View &v2) {
	
	std::vector<cv::Vec2f> originCorners, tfCorners1, tfCorners2;
	originCorners.push_back(cv::Vec2f(0, 0));
	originCorners.push_back(cv::Vec2f(img.cols-1, 0));
	originCorners.push_back(cv::Vec2f(img.cols-1, img.rows-1));
	originCorners.push_back(cv::Vec2f(0, img.rows-1));
	
	double col_ofs = std::min(0.0, v1.g*v1.sy*img.rows);
	for(int itCnrs=0; itCnrs<(int)originCorners.size(); itCnrs++) {
		double y =v1.sy*originCorners[itCnrs][1];
		double x = v1.sx*originCorners[itCnrs][0] + v1.g*y - col_ofs;
		tfCorners1.push_back(cv::Vec2f(x, y));
	}
	col_ofs = std::min(0.0, v2.g*v2.sy*img.rows);
	for(int itCnrs=0; itCnrs<(int)originCorners.size(); itCnrs++) {
		double y =v2.sy*originCorners[itCnrs][1];
		double x = v2.sx*originCorners[itCnrs][0] + v2.g*y - col_ofs;
		tfCorners2.push_back(cv::Vec2f(x, y));
	}
	cv::Mat H = cv::findHomography(tfCorners1, tfCorners2, CV_RANSAC);
	return H;
}




















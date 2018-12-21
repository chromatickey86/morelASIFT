//
//  MatchingModel.hpp
//  MatchingModel
//
//  Created by JoohyukYum on 2016. 7. 18..
//  Copyright © 2016년 JoohyukYum. All rights reserved.
//

#ifndef MatchingModel_hpp
#define MatchingModel_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>

namespace MatchingModel {
	//weighting parameter for fitting the matching model and measured matching accuracy
	const double W1 = 1.33;
	const double W2 = 1.71;
    
	class View {
	public:
		View() : sx(0.0), sy(0.0), g(0.0), sr(0.0), tau(0.0) {};
		View(double sx_, double sy_, double g_) : sx(sx_), sy(sy_), g(g_) { calcParam(); };
		View(double lat_, double lon_) {
			double t = 1.0/std::cos(lat_);
			double cos = std::cos(lon_), sin = std::sin(lon_);
			sx = std::sqrt(t*t*cos*cos + sin*sin)/t;
			sy = 1.0/std::sqrt(t*t*cos*cos + sin*sin);
			g = (1.0/t-t)*sin*cos;
			calcParam();
		}
		void calcParam() { calcParam(W1, W2); }
		void calcParam(const double &w1, const double &w2) {
			sr = w1*std::log10(sx/sy);
			tau = std::atan(g/w2);
		}

		//double lat, lon;
		double sx, sy, g;
		double sr, tau;
	};
	
	class InfoMatch {
	public:
		InfoMatch() : repeatability(0.0), matchingScore(0.0), correspondance(0), correctMatches(0),
		nQueryKps(0), nTrainKps(0), idImg(0) {};
		InfoMatch(int id, View &vwQ, View &vwT, double repeat, double corres, double matches, double mtscore,
					 double nQkps, double nTkps) : idImg(id), vwQuery(vwQ), vwTrain(vwT), repeatability(repeat),
		matchingScore(mtscore), correspondance(corres), correctMatches(matches), nQueryKps(nQkps), nTrainKps(nTkps) {};
		
		double repeatability, matchingScore;
		int correspondance, correctMatches, nQueryKps, nTrainKps;
		View vwQuery, vwTrain;
		int idImg;
	};
	
	//Build//
	//0. Display
	void display(const std::string &nmWnd, const std::vector<InfoMatch> &infoMatches);
	void display(const std::string &nmWnd, const std::vector<View> &views);
	void displayHemisphere(const std::string &nmWnd, const int lat, const int lon);
	cv::Scalar colorCoding(const double &val, const double &min, const double &max);
	
	//1. ParameterCalcuation
	double calc_K(const std::vector<InfoMatch> &infoMatches);
	double calc_K(const std::vector<InfoMatch> &infoMatches, std::vector<double> &history);
	void calcModifiedViews(std::vector<View> &views, const int lat, const int lon);
    
    
	//2. ViewSelection
	int buildOriginalViews(std::vector<View> &views);
	int buildOriginalViews(std::vector<View> &views, int nTilts, float minTilt, double srTilt);
	int buildUniformDistrViews(std::vector<View> &views);
	int buildProposedViews(std::vector<View> &views);
	
	//Data//
	//1. TestImageGen
	void imageSynthesis(const cv::Mat src, const int transType);
	enum { AFFINE_TRANSFORM, PERSPECTIVE_TRANSFORM };
	//2. File I/O
	int wrtMat2csv(const std::string &fnm, const cv::Mat mat);
	cv::Mat readXml(const std::string &fnm); //from Matching util.(matlab)
	int wrtXml(const std::string &fnm, cv::Mat mat);
	void readInfoMatches(std::vector<InfoMatch> &infoMatches);
	void wrtInfoMatches2csv(const std::string &fnm, const std::vector<InfoMatch> &infoMatches);
	void wrtViews2xml(const std::string & fnm, const std::vector<View> &views);
	void readViewsFromXml(const std::string & fnm, std::vector<View> &views);
	
	//Transform//
	void bound(int x, int y, double cos, double sin, int &xmin, int &xmax, int &ymin, int &ymax);
	void rotation(cv::Vec3f &vec, const double rx, const double ry, const double rz);
	cv::Mat perspectiveTransform(double distance, double zoom, const cv::Mat img, const double lat, const double lon, const int vertTilt);
	cv::Mat getHomographyFromFront(double distance, double zoom, const cv::Mat img, const double lat, const double lon, const int vertTilt);
	cv::Mat getHomographyFromFront(const cv::Mat img, const double lat, const double lon, const int vertTilt);
	cv::Mat getHomographyFromFront(const cv::Mat img, const double sx, const double sy, const double g);
	cv::Mat getHomography(const cv::Mat img, const View &v1, const View &v2);
	cv::Mat getHomographyFromTilted(double distance, double zoom, const cv::Mat img, const double lat, const double lon, const int vertTilt);
	cv::Mat affineTransform(const cv::Mat img, const double lat, const double lon, const int vertTilt);
	cv::Mat affineTransform(const cv::Mat img, const double sx, const double sy, const double g);
}


#endif /* MatchingModel_hpp */

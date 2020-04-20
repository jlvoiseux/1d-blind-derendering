#include <opencv2/core/core.hpp>
#include <iostream>
#include <math.h>

static void meshgrid(const cv::Mat &xgv, const cv::Mat &ygv, cv::Mat &X, cv::Mat &Y)
{
	cv::repeat(xgv.reshape(1, 1), ygv.total(), 1, X);
	cv::repeat(ygv.reshape(1, 1).t(), 1, xgv.total(), Y);
}

// helper function (maybe that goes somehow easier)
static void meshgridTest(const cv::Range &xgv, const cv::Range &ygv, cv::Mat &X, cv::Mat &Y)
{
	std::vector<int> t_x, t_y;
	for (int i = xgv.start; i <= xgv.end; i++) t_x.push_back(i);
	for (int i = ygv.start; i <= ygv.end; i++) t_y.push_back(i);
	meshgrid(cv::Mat(t_x), cv::Mat(t_y), X, Y);
}

Mat mSin(cv::Mat inputMat)
{
	Mat temp1, temp2, temp3, temp4;
	cv::pow(inputMat, 3, temp1);
	cv::pow(inputMat, 5, temp2);
	cv::pow(inputMat, 7, temp3);
	cv::pow(inputMat, 9, temp4);
	return inputMat - temp1 / 6.0f + temp2 / 120.0f - temp3 / 5040.0f + temp4 / 362880.0f;
}

Mat mCos(cv::Mat inputMat)
{
	return mSin(M_PI/2 - inputMat);
}
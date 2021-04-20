#include <opencv2/opencv.hpp>

int main()
{
	cv::Mat img = cv::imread("res/data.jpg");
	
	
	cv::Mat out_image;
	
	cv::cvtColor(img, out_image, cv::COLOR_BGR2GRAY);

	cv::imshow("w1", img);
	cv::imshow("w2", out_image);
	
	cv::waitKey();
	cv::destroyAllWindows();

}


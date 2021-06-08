#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/saturate.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <cmath>

#include "feature.h"
#include "processing.h"
#include "color.h"


int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cout << "usage: " << argv[0] << " path_to_input_image path_to_output_image" << std::endl;
        return 1;
    }

    cv::Mat org_image = cv::imread(argv[1]);
    auto pyrGreen = processing::Pyramid(org_image, color::GREEN);
    std::vector<cv::Rect2i> resultsGreen = pyrGreen.analyze();

    auto pyrWhite = processing::Pyramid(org_image, color::WHITE);
    std::vector<cv::Rect2i> resultsWhite = pyrWhite.analyze();

    auto pyrWhiteWider = processing::Pyramid(org_image, color::WHITE_WIDER);
    std::vector<cv::Rect2i> resultsWhiteWider = pyrWhiteWider.analyze();


    for (auto rect : resultsGreen) {
        cv::rectangle(org_image, rect, cv::Scalar(50), 4);
    }

    for (auto rect : resultsWhite) {
        cv::rectangle(org_image, rect, cv::Scalar(50));
    }

    for (auto rect : resultsWhiteWider) {
        cv::rectangle(org_image, rect, cv::Scalar(50));
    }

    cv::imshow("result", org_image);
    cv::waitKey(-1);

    return 0;
}

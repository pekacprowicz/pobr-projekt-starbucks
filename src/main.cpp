#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/saturate.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>

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

    std::vector<cv::Rect2i> results;
    results.insert(std::begin(results), std::begin(resultsGreen), std::end(resultsGreen));
    results.insert(std::begin(results), std::begin(resultsWhite), std::end(resultsWhite));
    results.insert(std::begin(results), std::begin(resultsWhiteWider), std::end(resultsWhiteWider));
    results = structures::concatOverlappingRectanglesSmallestArea(results);

    for (auto rect : results) {
        cv::rectangle(org_image, rect, cv::Scalar(21, 37, 213), 4);
    }

    cv::imwrite(argv[2], org_image);

    return 0;
}

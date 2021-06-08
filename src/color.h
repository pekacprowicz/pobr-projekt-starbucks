#ifndef POBR_PROJEKT_COLOR_H
#define POBR_PROJEKT_COLOR_H

namespace color {

    enum PatternColor {
        GREEN,
        WHITE,
        WHITE_WIDER
    };

    bool inRange(cv::Vec3b &pixel, PatternColor searchedColor) {
        bool H_condition, S_condition, V_condition;
        switch (searchedColor) {
            case GREEN:
                H_condition = pixel[0] > 55 && pixel[0] < 100;
                S_condition = pixel[1] > 60;
                V_condition = pixel[2] > 40 && pixel[2] < 223;
                break;
            case WHITE:
                H_condition = pixel[0] > 0 && pixel[0] < 179;
                S_condition = pixel[1] < 60;
                V_condition = pixel[2] > 210;
                break;
            case WHITE_WIDER:
                H_condition = pixel[0] > 0 && pixel[0] < 179;
                S_condition = pixel[1] < 100;
                V_condition = pixel[2] > 175;
                break;
        }
        return H_condition && S_condition && V_condition;
    }

    bool inRange(cv::Vec3b &pixel, uchar low, uchar high) {
        return pixel[2] > low && pixel[2] < high;
    }

    cv::Mat_<uchar> binaryDebug(const cv::Mat_<uchar> &image) {
        cv::Mat_<uchar> debug = cv::Mat_<uchar>(image.rows, image.cols);
        image.forEach([&](uchar &pixel, const int position[])-> void {
            debug.at<uchar>(position[0], position[1]) = pixel == 1 ? 255 : 0;
        });
//        cv::imshow("binaryDebug", debug);
//        cv::waitKey(-1);
        return debug;
    }

    cv::Mat_<uchar> createBinaryImage(const cv::Mat &image, PatternColor searchedColor) {
        cv::Mat hsv_image = cv::Mat_<uchar>(image.rows, image.cols);
        cv::cvtColor(image, hsv_image, cv::COLOR_BGR2HSV);
        cv::Mat_<uchar> binaryImage = cv::Mat_<uchar>(image.rows, image.cols);

        hsv_image.forEach<cv::Vec3b>([&](auto &pixel, const int position[])-> void {
            binaryImage.at<uchar>(position[0], position[1]) = inRange(pixel, searchedColor) ? 1 : 0;
        });

        return binaryImage;
    }

    void multipleThresholdValue(const cv::Mat &image, std::vector<cv::Mat_<uchar>>& result, PatternColor searchedColor) {
        cv::Mat hsv_image = cv::Mat_<uchar>(image.rows, image.cols);
        cv::cvtColor(image, hsv_image, cv::COLOR_BGR2HSV);

        if (image.rows <= 0 || image.cols <= 0) return;
        cv::Mat_<uchar> thresholdImage = cv::Mat_<uchar>(image.rows, image.cols);
        for (int i = 0; i < 13; i++) {
            hsv_image.forEach<cv::Vec3b>([&](auto &pixel, const int position[])-> void {
                if (inRange(pixel, searchedColor))
                    thresholdImage.at<uchar>(position[0], position[1]) = inRange(pixel, 16 * i, 16 * i + 47) ? 1 : 0;
                else
                    thresholdImage.at<uchar>(position[0], position[1]) = 0;
            });
            result.push_back(thresholdImage.clone());
        }
        return;
    }
}

#endif

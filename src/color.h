#ifndef POBR_PROJEKT_COLOR_H
#define POBR_PROJEKT_COLOR_H

namespace color {

    enum PatternColor {
        GREEN,
        WHITE,
        WHITE_WIDER
    };

    void bgrToHsv(const cv::Mat& source, cv::Mat& destination) {
        source.forEach<cv::Vec3b>([&](auto &pixel, const int position[])-> void {
            double Rp = pixel[2] / 255.0;
            double Gp = pixel[1] / 255.0;
            double Bp = pixel[0] / 255.0;
            double Cmax = std::max(std::max(Rp, Gp), Bp);
            double Cmin = std::min(std::min(Rp, Gp), Bp);
            double delta = Cmax - Cmin;

            if (delta == 0) {
                destination.at<cv::Vec3b>(position[0], position[1])[0] = 0;
            } else if (Cmax == Rp) {
                destination.at<cv::Vec3b>(position[0], position[1])[0] = std::lround((60 * std::fmod(((Gp - Bp) / delta), 6.0)) / 2);
            } else if (Cmax == Gp) {
                destination.at<cv::Vec3b>(position[0], position[1])[0] = std::lround((60 * (((Bp - Rp) / delta) + 2.0)) / 2);
            } else if (Cmax == Bp) {
                destination.at<cv::Vec3b>(position[0], position[1])[0] = std::lround((60 * (((Rp - Gp) / delta) + 4.0)) / 2);
            }

            if (Cmax == 0) {
                destination.at<cv::Vec3b>(position[0], position[1])[1] = 0;
            } else {
                destination.at<cv::Vec3b>(position[0], position[1])[1] = std::lround((delta / Cmax) * 255);
            }

            destination.at<cv::Vec3b>(position[0], position[1])[2] = std::lround(Cmax * 255);
        });
    }

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
        return debug;
    }

    cv::Mat_<uchar> createBinaryImage(const cv::Mat &image, PatternColor searchedColor) {
        cv::Mat hsv_image = cv::Mat_<cv::Vec3b>(image.rows, image.cols);
        bgrToHsv(image, hsv_image);
        cv::Mat_<uchar> binaryImage = cv::Mat_<uchar>(image.rows, image.cols);

        hsv_image.forEach<cv::Vec3b>([&](auto &pixel, const int position[])-> void {
            binaryImage.at<uchar>(position[0], position[1]) = inRange(pixel, searchedColor) ? 1 : 0;
        });

        return binaryImage;
    }

    void multipleThresholdValue(const cv::Mat &image, std::vector<cv::Mat_<uchar>>& result, PatternColor searchedColor) {
        cv::Mat hsv_image = cv::Mat_<cv::Vec3b>(image.rows, image.cols);
        bgrToHsv(image, hsv_image);

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

#ifndef POBR_PROJEKT_FEATURE_H
#define POBR_PROJEKT_FEATURE_H

#include <cmath>

#define _PI 3.1415926535897932384626433832795

namespace feature {
    struct BasicGeometricFeatures {
        int L;
        int S;
        cv::Mat &_image;

        BasicGeometricFeatures(cv::Mat &image) : _image(image) {
            auto L_S = getPerimeterAndArea(_image);
            L = std::get<0>(L_S);
            S = std::get<1>(L_S);
        }

        BasicGeometricFeatures operator=(const BasicGeometricFeatures &right) {
            this->L = right.L;
            this->S = right.S;
            this->_image = right._image;
        }

        std::tuple<int, int> getPerimeterAndArea(cv::Mat& image) {
            int s = 0;
            int l = 0;
            for (int i = 0; i < image.rows; ++i)
                for (int j = 0; j < image.cols; ++j) {
                    if (image.at<uchar>(i, j) == 1) {
                        s++;
                        if (i != 0) {
                            if (j != 0 && image.at<uchar>(i - 1, j - 1) == 0) { l++; continue; }
                            if (image.at<uchar>(i - 1, j) == 0) { l++; continue; }
                            if (j != image.cols - 1 && image.at<uchar>(i - 1, j + 1) == 0) { l++; continue; }
                        }
                        if (j != 0) {
                            if (image.at<uchar>(i, j - 1) == 0) { l++; continue; }
                            if (i != image.rows - 1 && image.at<uchar>(i + 1, j - 1) == 0) { l++; continue; }
                        }
                        if (j != image.cols - 1) {
                            if (image.at<uchar>(i, j + 1) == 0) { l++; continue; }
                            if (i != image.rows - 1 && image.at<uchar>(i + 1, j + 1) == 0) { l++; continue; }
                        }
                        if (i != image.rows - 1 && image.at<uchar>(i + 1, j) == 0) { l++; continue; }
                    }
                }
            return std::make_tuple(l, s);
        }
    };

    struct Moments {
        long double m00, m10, m20, m30, m01, m02, m03, m11, m12, m21;
        double M00, M01, M10, M11, M20, M02, M21, M12, M30, M03;
        double M1, M2, M3, M4, M5, M6, M7, M8, M9;
        double i, j;
        cv::Mat &_image;

        Moments(cv::Mat &image) : _image(image){
            m00 = get_m_p_q(0, 0);
            m10 = get_m_p_q(1, 0);
            m20 = get_m_p_q(2, 0);
            m30 = get_m_p_q(3, 0);
            m01 = get_m_p_q(0, 1);
            m02 = get_m_p_q(0, 2);
            m03 = get_m_p_q(0, 3);
            m11 = get_m_p_q(1, 1);
            m12 = get_m_p_q(1, 2);
            m21 = get_m_p_q(2, 1);

            i = m10 / m00;
            j = m01 / m00;

            M00 = m00;
            M01 = m01 - (m01 / m00) * m00;
            M10 = m10 - (m10 / m00) * m00;
            M11 = m11 - m10 * m01 / m00;
            M20 = m20 - std::pow(m10, 2) / m00;
            M02 = m02 - std::pow(m01, 2) / m00;
            M21 = m21 - 2 * m11 * i - m20 * j + 2 * m01 * std::pow(i, 2);
            M12 = m12 - 2 * m11 * j - m02 * i + 2 * m10 * std::pow(j, 2);
            M30 = m30 - 3 * m20 * i + 2 * m10 * std::pow(i, 2);
            M03 = m03 - 3 * m02 * j + 2 * m01 * std::pow(j, 2);

            M1 = (M20 + M02) / std::pow(m00, 2);
            M2 = (std::pow(M20 - M02, 2) + 4 * std::pow(M11, 2)) / std::pow(m00, 4);
            M3 = (std::pow(M30 - 3 * M12, 2) + std::pow(3 * M21 - M03, 2)) / std::pow(m00, 5);
            M4 = (std::pow(M30 - M12, 2) + std::pow(M21 + M03, 2)) / std::pow(m00, 5);
            M5 = ((M30 - 3 * M12) * (M30 + M12) * (std::pow(M30 + M12, 2) - 3 * std::pow(M21 + M03, 2)) +
                  (3 * M21 - M03) * (M21 + M03) * (3 * std::pow(M30 + M12, 2) - std::pow(M21 + M03, 2))) / std::pow(m00, 10);
            M6 = ((M20 - M02) * (std::pow(M30 + M12, 2) - std::pow(M21 + M03, 2)) + 4 * M11 * (M30 + M12) * (M21 + M03)) / std::pow(m00, 7);
            M7 = (M20 * M02 - std::pow(M11, 2)) / std::pow(m00, 4);
        }

        long double get_m_p_q(int p, int q) {
            long double m = 0;
            for (int i = 0; i < _image.rows; ++i) {
                for (int j = 0; j < _image.cols; ++j) {
                    m += std::pow(j, p) * std::pow(i, q) * _image.at<uchar>(i, j);
                }
            }
            return m;
        }

        static long double get_m_p_q(cv::Mat image, int p, int q) {
            long double m = 0;
            for (int i = 0; i < image.rows; ++i) {
                for (int j = 0; j < image.cols; ++j) {
                    m += std::pow(j, p) * std::pow(i, q) * image.at<uchar>(i, j);
                }
            }
            return m;
        }
    };

    struct Centroid {
        double i, j;
        cv::Mat &_image;

        Centroid(cv::Mat &image) : _image(image) {
            long double m00, m01, m10;
            m00 = Moments::get_m_p_q(_image, 0, 0);
            m01 = Moments::get_m_p_q(_image, 0, 1);
            m10 = Moments::get_m_p_q(_image, 1, 0);

            i = m10 / m00;
            j = m01 / m00;
        }

        Centroid operator=(const Centroid &right) {
            this->i = right.i;
            this->j = right.j;
            this->_image = right._image;
        }

        double distance(cv::Point point) {
            return std::sqrt(std::pow(point.x - i, 2) + std::pow(point.y - j, 2));
        }
    };

    struct Coefficients {
        double W4;
        double W9;
        cv::Mat &_image;
        Centroid _centroid;
        BasicGeometricFeatures _features;

        Coefficients(cv::Mat &image) : _image(image), _centroid(Centroid(_image)), _features(BasicGeometricFeatures(_image)) {
            W4 = getW4();
            W9 = getW9();
        }

        double getW4() {
            long double sum = 0;
            _image.forEach<uchar>([&](uchar &pixel, const int position[])-> void {
                sum += std::pow(_centroid.distance(cv::Point(position[0], position[1])), 2);
            });
            return _features.S / std::sqrt(2 * _PI * sum);
        }

        double getW9() {
            return (2 * std::sqrt(_PI * _features.S)) / _features.L;
        }
    };
}

#endif 

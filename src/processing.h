#ifndef POBR_PROJEKT_PROCESSING_H
#define POBR_PROJEKT_PROCESSING_H

#include "color.h"
#include "structures.h"

namespace processing {

    int sumOfRow(const cv::Mat_<uchar>& image, int row) {
//        std::cout << "Row: " << row << std::endl;
//        cv::imshow("sumOfRow", color::binaryDebug(image));
//        cv::waitKey(-1);
        for (int col = 0; col < image.cols; col++) {
            if (image.at<uchar>(row, col) == 1) return 1;
        }
        return 0;
    }

    int sumOfCol(const cv::Mat_<uchar>& image, int col) {
//        std::cout << "Col: " << col << std::endl;
//        cv::imshow("sumOfCol", color::binaryDebug(image));
//        cv::waitKey(-1);
        for (int row = 0; row < image.rows; row++) {
            uchar val = image.at<uchar>(row, col);
            if (val == 1) return 1;
        }
        return 0;
    }

    double distance(cv::Point pointA, cv::Point pointB) {
        return std::sqrt(std::pow(pointA.x - pointB.x, 2) + std::pow(pointA.y - pointB.y, 2));
    }

    cv::Mat_<cv::Vec3b> expand(const cv::Mat &image, int rows, int cols) {
        int _rows = image.rows + rows;
        int _cols = image.cols + cols;
        //float data[_rows * _cols * 3];
//        cv::Vec3b data[_rows * _cols];
        cv::Mat_<cv::Vec3b>  _I = cv::Mat(_rows, _cols, CV_8UC3);

        for (int i = 0; i < (rows/2); i++)
            for (int j = 0; j < (cols/2); j++)
                for (int ch = 0; ch < 3; ch++)
                    _I(i, j)[ch] = 0;
//                    data[getDataIndex(i, j, _cols)][ch] = 0;

        for (int i = _I.rows - (rows/2); i < _I.rows; i++)
            for (int j = _I.cols - (cols/2); j < _I.cols; j++)
                for (int ch = 0; ch < 3; ch++)
                    _I(i, j)[ch] = 0;
//                    data[getDataIndex(i, j, _cols)][ch] = 0;

        for (int i = 0; i < image.rows; i++)
            for (int j = 0; j < image.cols; j++)
                for (int ch = 0; ch < 3; ch++) {
                    uchar val = image.at<cv::Vec3b>(i, j)[ch];
                    _I(i + (rows / 2), j + (cols / 2))[ch] = val;
                }

        return _I;
    }

    cv::Mat filter(const cv::Mat &image, const cv::Mat_<int> &kernel, double scalar = 1) {
        cv::Mat res(image.rows, image.cols, CV_32FC3);
        int kernelSize = kernel.cols;
        switch (image.channels()) {
            case 3:
                cv::Mat_<cv::Vec3f> _I = expand(image, kernelSize - 1, kernelSize - 1);
                cv::Mat_<int> _ker = kernel;
                float sum;
                for (int i = 2; i < _I.rows - (kernelSize / 2); i++) {
                    for (int j = 2; j < _I.cols - (kernelSize / 2); j++) {
                        for (int ch = 0; ch < 3; ch++) {
                            sum = 0;
                            for (int k = 0; k < 5; k++)
                                for (int m = 0; m < 5; m++) {
                                    sum = sum + _I.at<cv::Vec3f>(i - (kernelSize / 2) + k,
                                                                                j - (kernelSize / 2) + m)[ch] * _ker(k, m);
                                }
                            res.at<cv::Vec3f>(i - (kernelSize/2), j - (kernelSize/2))[ch] = sum * scalar;
//                            data[getDataIndex(i - (kernelSize/2), j - (kernelSize/2), image.cols)][ch] = sum * scalar;
                        }
                    }
                }
                break;
        }
        return res;
    }

    cv::Mat gaussianFilter(const cv::Mat &image) {
        float kernel_data[25] = {1, 4, 6, 4, 1, 4, 16, 24, 16, 4, 6, 24, 36, 24, 6, 4, 16, 24, 16, 4, 1, 4, 6, 4, 1};
        cv::Mat_<float> kernel = cv::Mat_<float>(5, 5, kernel_data) * 1.0/16.0;
        return filter(image, kernel, 1.0/16.0);
    }

    cv::Mat downSize(const cv::Mat &image) {
        cv::Mat res(image.rows / 2, image.cols / 2, CV_32FC3);
        switch (image.channels()) {
            case 3:
                for (int i = 0; i < image.rows; i++) {
                    if (i % 2 == 0) continue;
                    for (int j = 0; j < image.cols; j++) {
                        if (j % 2 == 0) continue;
                        res.at<cv::Vec3f>(i / 2, j / 2) = image.at<cv::Vec3b>(i, j);
                    }
                }
                break;
        }
        return res;
    }

    template<typename MatType>
    struct WindowStruct {
        cv::Mat_<MatType> window;
        cv::Range rangeRows;
        cv::Range rangeCols;
    };

    template<typename MatType>
    void createWindows(const cv::Mat &image, std::vector<WindowStruct<MatType>> &_windows) {
        int windowSize = 100;
        int stepSize = 8;
        int _cols = (int) std::ceil(image.cols / (double) windowSize);
        int _rows = (int) std::ceil(image.rows / (double) windowSize);

        cv::Mat_<MatType> _window;
        int colsRangeBegin, colsRangeEnd, rowsRangeBegin, rowsRangeEnd;
        int col = 0;
        do {
            colsRangeEnd = std::min(col * stepSize + windowSize, image.cols);
            colsRangeBegin = std::max(colsRangeEnd != image.cols ? col * stepSize : (col - 1) * stepSize, 0);
            int row = 0;
            do {
                rowsRangeEnd = std::min(row * stepSize + windowSize, image.rows);
                rowsRangeBegin = std::max(rowsRangeEnd != image.rows ? row * stepSize : (row - 1) * stepSize, 0);
                _window = image(cv::Range(rowsRangeBegin, rowsRangeEnd),
                                cv::Range(colsRangeBegin, colsRangeEnd)
                );
//                cv::imshow("windowbindbg", color::binaryDebug(_window));
//                cv::waitKey(-1);
                _windows.push_back(WindowStruct<uchar>{_window, cv::Range(rowsRangeBegin, rowsRangeEnd), cv::Range(colsRangeBegin, colsRangeEnd)});
                row++;
            } while (row * stepSize + windowSize < image.rows);
            col++;
        } while (col * stepSize + windowSize < image.cols);

//        for (int col = 0; col < _cols; col++) {
//            colsRangeEnd = (col + 2) * windowSize < image.cols ? (col + 2) * windowSize : image.cols;
//            colsRangeBegin = std::max(colsRangeEnd != image.cols ? col * windowSize : (col - 1) * windowSize, 0);
//            for (int row = 0; row < _rows; row++) {
//                rowsRangeEnd = (row + 2) * windowSize < image.rows ? (row + 2) * windowSize : image.rows;
//                rowsRangeBegin = std::max(rowsRangeEnd != image.rows ? row * windowSize : (row - 1) * windowSize, 0);
//                _window = image(cv::Range(rowsRangeBegin, rowsRangeEnd),
//                                cv::Range(colsRangeBegin, colsRangeEnd)
//                );
//                cv::imshow("windowbindbg", color::binaryDebug(_window));
//                cv::waitKey(-1);
//                _windows.push_back(WindowStruct<uchar>{_window, cv::Range(rowsRangeBegin, rowsRangeEnd), cv::Range(colsRangeBegin, colsRangeEnd)});
//            }
//        }
    }

    cv::Rect2i getBoundingBox(const cv::Mat_<uchar> &image) {
        int x_min, x_max, y_min, y_max;
        x_min = image.cols - 1; x_max = 0; y_min = image.rows - 1; y_max = 0;
        image.forEach([&](uchar &pixel, const int position[]) -> void {
            if (pixel == 1) {
                x_min = position[1] < x_min ? position[1] : x_min;
                x_max = position[1] > x_max ? position[1] : x_max;
                y_min = position[0] < y_min ? position[0] : y_min;
                y_max = position[0] > y_max ? position[0] : y_max;
            }
        });
        if (x_max - x_min < 0 || y_max - y_min < 0) return cv::Rect2i(0, 0, 0, 0);
        return cv::Rect2i(std::max(x_min - 1, 0), std::max(y_min - 1, 0),
                          std::min(x_max - x_min + 2, image.cols), std::min(y_max - y_min + 2, image.rows)
                          );
    }

    cv::Rect2i adjustBoundingBox(const cv::Mat_<uchar> &image, cv::Rect2i &bBox) {
        bool shouldMoveLeftToLeft = false, shouldMoveLeftToRight = true;
        for (int i = bBox.y; i < bBox.y + bBox.height + 1; i++) {
            if (image.at<uchar>(i, bBox.x - 1) == 1)
                shouldMoveLeftToLeft = true;
            if (image.at<uchar>(i, bBox.x + 1) == 1)
                shouldMoveLeftToRight = false;
        }

        while (shouldMoveLeftToLeft) {
            shouldMoveLeftToLeft = false;
            if (bBox.x > 0) {bBox.x--; bBox.width++;}
            else break;
            for (int i = bBox.y; i < bBox.y + bBox.height + 1; i++)
                if (image.at<uchar>(i, bBox.x - 1) == 1)
                    shouldMoveLeftToLeft = true;
        }

        while (shouldMoveLeftToRight) {
            if (bBox.width > 0) {bBox.x++; bBox.width--;}
            else break;
            for (int i = bBox.y; i < bBox.y + bBox.height + 1; i++)
                if (image.at<uchar>(i, bBox.x + 1) == 1)
                    shouldMoveLeftToRight = false;
        }

        bool shouldMoveRightToRight = false, shouldMoveRightToLeft = true;
        for (int i = bBox.y; i < bBox.y + bBox.height + 1; i++) {
            if (image.at<uchar>(i, bBox.x + bBox.width + 1) == 1)
                shouldMoveRightToRight = true;
            if (image.at<uchar>(i, bBox.x + bBox.width - 1) == 1)
                shouldMoveRightToLeft = false;
        }

        while (shouldMoveRightToRight) {
            shouldMoveRightToRight = false;
            if (bBox.x + bBox.width < image.cols) bBox.width++;
            else break;
            for (int i = bBox.y; i < bBox.y + bBox.height + 1; i++)
                if (image.at<uchar>(i, bBox.x + bBox.width + 1) == 1)
                    shouldMoveRightToRight = true;
        }

        while (shouldMoveRightToLeft) {
            if (bBox.width > 0) {bBox.width--;}
            else break;
            for (int i = bBox.y; i < bBox.y + bBox.height + 1; i++)
                if (image.at<uchar>(i, bBox.x + bBox.width - 1) == 1)
                    shouldMoveRightToLeft = false;
        }

        bool shouldMoveUpToUp = false, shouldMoveUpToDown = true;
        for (int i = bBox.x; i < bBox.x + bBox.width + 1; i++) {
            if (image.at<uchar>(bBox.y - 1, i) == 1)
                shouldMoveUpToUp = true;
            if (image.at<uchar>(bBox.y + 1, i) == 1)
                shouldMoveUpToDown = false;
        }

        while (shouldMoveUpToUp) {
            shouldMoveUpToUp = false;
            if (bBox.y > 0) {bBox.y--; bBox.height++;}
            else break;
            for (int i = bBox.x; i < bBox.x + bBox.width + 1; i++)
                if (image.at<uchar>(bBox.y - 1, i) == 1)
                    shouldMoveUpToUp = true;
        }

        while (shouldMoveUpToDown) {
            if (bBox.height > 0) {bBox.y++; bBox.height--;}
            else break;
            for (int i = bBox.x; i < bBox.x + bBox.width + 1; i++)
                if (image.at<uchar>(bBox.y + 1, i) == 1)
                    shouldMoveUpToDown = false;
        }

        bool shouldMoveDownToDown = false, shouldMoveDownToUp = true;
        for (int i = bBox.x; i < bBox.x + bBox.width + 1; i++) {
            if (image.at<uchar>(bBox.y + bBox.height + 1, i) == 1)
                shouldMoveDownToDown = true;
            if (image.at<uchar>(bBox.y + bBox.height - 1, i) == 1)
                shouldMoveDownToUp = false;
        }

        while (shouldMoveDownToDown) {
            shouldMoveDownToDown = false;
            if (bBox.y + bBox.height < image.rows) bBox.height++;
            else break;
            for (int i = bBox.x; i < bBox.x + bBox.width + 1; i++)
                if (image.at<uchar>(bBox.y + bBox.height + 1, i) == 1)
                    shouldMoveDownToDown = true;
        }

        while (shouldMoveDownToUp) {
            if (bBox.height > 0) bBox.height--;
            else break;
            for (int i = bBox.x; i < bBox.x + bBox.width + 1; i++)
                if (image.at<uchar>(bBox.y + bBox.height - 1, i) == 1)
                    shouldMoveDownToUp = false;
        }
    }

    void expandBoundingBox(const cv::Mat_<uchar> &image, cv::Rect2i &bBox, feature::Centroid centroid) {
        cv::Point centroidPoint = cv::Point(bBox.x + (int) centroid.i, bBox.y + (int) centroid.j);
        cv::Point minDistancePoint = centroidPoint;
//        image.forEach([&](auto &pixel, const int position[])-> void {
//            if (position[0] > bBox.x && position[0] < bBox.x + bBox.width && position[1] > bBox.y && position[1] < bBox.y + bBox.height) {
//                if (pixel == 1) {
//                    tempPoint = cv::Point(position[0], position[1]);
//                    if (distance(tempPoint, centroidPoint) < minDistance) {
//                        minDistance = distance(tempPoint, centroidPoint);
//                        minDistancePoint = tempPoint;
//                    }
//                }
//            }
//        });
        int step = 1;
        while (true) {
            minDistancePoint.x += step;
            if (image.at<uchar>(minDistancePoint.y, minDistancePoint.x) == 1) { break; }
            minDistancePoint.y += step;
            if (image.at<uchar>(minDistancePoint.y, minDistancePoint.x) == 1) { break; }
            step++;
            minDistancePoint.x -= step;
            if (image.at<uchar>(minDistancePoint.y, minDistancePoint.x) == 1) { break; }
            minDistancePoint.y -= step;
            if (image.at<uchar>(minDistancePoint.y, minDistancePoint.x) == 1) { break; }
            step++;
        }
        cv::Rect2i minimumRect = cv::Rect2i(minDistancePoint.x, minDistancePoint.y, 1, 1);

        cv::Mat dbg = color::binaryDebug(image);
        cv::circle(dbg, centroidPoint, 2, cv::Scalar(125), 2);
        cv::circle(dbg, minDistancePoint, 2, cv::Scalar(150, 150, 150), 2);

        if (minDistancePoint.x < 0 || minDistancePoint.y < 0) return;

        cv::Rect2i tempMinRect = cv::Rect2i(minimumRect.x - 1, minimumRect.y - 1, minimumRect.width + 2, minimumRect.height + 2);

        for (int i = 0; i < 3; i++) {
            while (minimumRect.x > bBox.x && sumOfCol(image(tempMinRect), 0) > 0) {
//            cv::Mat _dbg = dbg.clone();
//            cv::rectangle(_dbg, minimumRect, cv::Scalar(69), 2);
//            cv::imshow("sumOfRow", _dbg);
//            cv::waitKey(-1);
                minimumRect.x--;
                minimumRect.width++;
                tempMinRect = cv::Rect2i(minimumRect.x - 1, minimumRect.y - 1, minimumRect.width + 2, minimumRect.height + 2);
            }

            while (minimumRect.x + minimumRect.width < bBox.x + bBox.width - 1 && sumOfCol(image(tempMinRect), tempMinRect.width - 1) > 0) {
                minimumRect.width++;
                tempMinRect = cv::Rect2i(minimumRect.x - 1, minimumRect.y - 1, minimumRect.width + 2, minimumRect.height + 2);
            }

            while (minimumRect.y > bBox.y && sumOfRow(image(tempMinRect), 0) > 0) {
                minimumRect.y--;
                minimumRect.height++;
                tempMinRect = cv::Rect2i(minimumRect.x - 1, minimumRect.y - 1, minimumRect.width + 2, minimumRect.height + 2);
            }

            while (minimumRect.y + minimumRect.height < bBox.y + bBox.height - 1 && sumOfRow(image(tempMinRect), tempMinRect.height - 1) > 0) {
                minimumRect.height++;
                tempMinRect = cv::Rect2i(minimumRect.x - 1, minimumRect.y - 1, minimumRect.width + 2, minimumRect.height + 2);
            }
            tempMinRect = cv::Rect2i(minimumRect.x - 1, minimumRect.y - 1, minimumRect.width + 2, minimumRect.height + 2);
        }

        minimumRect = tempMinRect.size().area() > 0.1 * bBox.size().area() ? tempMinRect : bBox;
        cv::rectangle(dbg, minimumRect, cv::Scalar(69), 2);
        cv::rectangle(dbg, bBox, cv::Scalar(47), 2);

//        cv::imshow("expandBB", dbg);
//        cv::waitKey(-1);

        bBox = minimumRect;
    }

    cv::Mat erosion(const cv::Mat &image, int size = 1) {
        cv::Mat output = image.clone();

        cv::Rect2i mat_rect = cv::Rect2i(cv::Point(0, 0), image.size());

        output.forEach<uchar>([&](auto &pixel, const int position[]) -> void {
            cv::Rect2i sub_rect = cv::Rect2i(
                    position[1] - size, position[0] - size,
                    2 * size + 1, 2 * size + 1
            ) & mat_rect;

            cv::Mat sub_image = image(sub_rect);

            sub_image.forEach<uchar>([&](uchar &sub_pixel, const int sub_position[]) -> void {
                pixel = std::min(pixel, sub_pixel);
            });
        });

        return output;
    }

    cv::Mat dilation(const cv::Mat &image, int size = 1) {
        cv::Mat output = image.clone();

        cv::Rect2i mat_rect = cv::Rect2i(cv::Point(0, 0), image.size());

        output.forEach<uchar>([&](auto &pixel, const int position[]) -> void {
            cv::Rect2i sub_rect = cv::Rect2i(
                    position[1] - size, position[0] - size,
                    2 * size + 1, 2 * size + 1
            ) & mat_rect;

            cv::Mat sub_image = image(sub_rect);

            sub_image.forEach<uchar>([&](uchar &sub_pixel, const int sub_position[]) -> void {
                pixel = std::max(pixel, sub_pixel);
            });
        });

        return output;
    }

    cv::Mat opening(const cv::Mat &matrix, int size = 1) {
        return dilation(erosion(matrix, size), size);
    }

    cv::Mat closing(const cv::Mat &matrix, int size = 1) {
        return erosion(dilation(matrix, size), size);
    }

    cv::Point getGeometricCenter(cv::Mat image) {
        int x_min, x_max, y_min, y_max;
        x_min = image.rows - 1; x_max = 0; y_min = image.cols - 1; y_max = 0;
        image.forEach<int>([&](auto &pixel, const int position[])-> void {
            if (pixel == 1) {
                x_min = position[0] < x_min ? position[0] : x_min;
                x_max = position[0] > x_max ? position[0] : x_max;
                y_min = position[1] < y_min ? position[1] : y_min;
                y_max = position[1] > y_max ? position[1] : y_max;
            }
        });
        return cv::Point((int) ((x_max + x_min) / 2), (int) ((y_max + y_min)) / 2);
    }

    struct Pyramid {
        struct Element {
            cv::Mat_<uchar> levelImage;
            int level;
        };
        std::vector<Element> _images;
        cv::Mat_<cv::Vec3b> _baseImage;
        cv::Mat_<uchar> _baseImageBinary;
        int sizeThreshold = 150;
        color::PatternColor _pattern;


        Pyramid(const cv::Mat_<cv::Vec3b> &image, color::PatternColor pattern = color::GREEN) : _baseImage(image), _pattern(pattern) {
            cv::Mat_<cv::Vec3b> tempImage = _baseImage;
            _baseImageBinary = color::createBinaryImage(_baseImage, _pattern);
            int level = 0;
            do {
                _images.push_back(Element{opening(color::createBinaryImage(tempImage, _pattern)), level++});
                tempImage = downSize(gaussianFilter(tempImage));
            } while (tempImage.rows > sizeThreshold && tempImage.cols > sizeThreshold);
        }

        void showEveryLevel() {
            for (const auto &imgEl : _images) {
                cv::imshow("level:" + std::to_string(imgEl.level), color::binaryDebug(imgEl.levelImage));
                cv::waitKey(-1);
            }
        }

        std::vector<cv::Rect2i> analyze() {
            std::vector<structures::RangesStruct> ROIs = getROIs();
            std::vector<cv::Rect2i> results;

            std::vector<cv::Mat_<uchar>> multipleThreshold;
            cv::Mat_<uchar> _thImg;
            cv::Rect2i rect;
            cv::Rect2i rectResult;
            int i = 0;
            int roiID;
            switch (_pattern) {
                case color::GREEN:
                    for (auto const &roi : ROIs) {
                        cv::Rect2i roiRect = structures::createRectFromRanges(roi.rangeRows, roi.rangeCols);
                        std::cout << " ::: rect: " << roiRect << std::endl;
                        processing::adjustBoundingBox(_baseImageBinary, roiRect);
                        _thImg = closing(_baseImageBinary(roiRect));
                        feature::BasicGeometricFeatures features(_thImg);
                        std::cout << "S: " << features.S << " ::: L: " << features.L << std::endl;
                        std::cout << "S/L: " << (double) features.S / (double) features.L << std::endl;
//                        cv::imshow("multipleThreshold", color::binaryDebug(_thImg));
//                        cv::waitKey(-1);
                        double s_l = (double) features.S / (double) features.L;
                        if (s_l > 2.1 && s_l < 2.8) {
                            feature::Centroid centroid(_thImg);
                            cv::Point geometricCenter = getGeometricCenter(_thImg);
                            feature::Coefficients thCoeff(_thImg);
                            std::cout << "W9: " << thCoeff.W9 << std::endl;
                            if (thCoeff.W9 > 0.084 && thCoeff.W9 < 0.16) {
                                std::cout << "distance/area: " << centroid.distance(geometricCenter) / roiRect.size().area() << std::endl;

                                if (centroid.distance(geometricCenter) / roiRect.size().area() < 0.001) {
//                                cv::imshow("greenfinal", _baseImage(roiRect));
//                                cv::waitKey(-1);
                                    expandBoundingBox(_baseImageBinary, rectResult, centroid);
                                    results.push_back(roiRect);
                                } else {
//                                    cv::imshow("greenfinal_before", _baseImage(roiRect));
//                                    cv::waitKey(-1);
                                    expandBoundingBox(_baseImageBinary, roiRect, centroid);
//                                    cv::imshow("greenfinal_after", _baseImage(roiRect));
//                                    cv::waitKey(-1);
                                    _thImg = closing(_baseImageBinary(roiRect));
                                    centroid = feature::Centroid(_thImg);
                                    geometricCenter = getGeometricCenter(_thImg);
                                    feature::Coefficients thCoeff = feature::Coefficients(_thImg);
                                    std::cout << "distance/area again: " << centroid.distance(geometricCenter) / roiRect.size().area() << std::endl;
                                    if (centroid.distance(geometricCenter) / roiRect.size().area() < 0.001 && thCoeff.W9 > 0.084 && thCoeff.W9 < 0.16) {
                                        results.push_back(roiRect);
                                    }
                                }
                            }
                        }

                        color::multipleThresholdValue(_baseImage(roiRect), multipleThreshold, _pattern);
                    }

                    for (const cv::Mat_<uchar>& thImg : multipleThreshold) {
                        roiID = i / 13;
                        rect = processing::getBoundingBox(thImg);
                        if (rect.size().area() == 0) continue;
                        _thImg = closing(thImg(rect));
                        feature::BasicGeometricFeatures features(_thImg);
                        std::cout << "S: " << features.S << " ::: L: " << features.L << std::endl;
                        std::cout << "S/L: " << (double) features.S / (double) features.L << std::endl;

                        double s_l = (double) features.S / (double) features.L;
                        if (s_l > 2.1 && s_l < 2.8) {
                            feature::Coefficients thCoeff(_thImg);

                            std::cout << "W9: " << thCoeff.W9 << std::endl;

                            if (thCoeff.W9 > 0.08 && thCoeff.W9 < 0.16) {
                                rectResult = structures::createRectFromRanges(ROIs.at(roiID).rangeRows, ROIs.at(roiID).rangeCols);
                                _thImg = _baseImageBinary(rectResult);
                                feature::Centroid centroid(_thImg);
                                cv::Point geometricCenter = getGeometricCenter(_thImg);
                                std::cout << "distance/area: " << centroid.distance(geometricCenter) / rectResult.size().area() << std::endl;
                                if (centroid.distance(geometricCenter) / rectResult.size().area() < 0.001) {
//                                    cv::imshow("multipleThreshold_before", color::binaryDebug(_baseImageBinary(rectResult)));
//                                    cv::waitKey(-1);
                                    expandBoundingBox(_baseImageBinary, rectResult, centroid);

//                                    cv::imshow("multipleThreshold_after", color::binaryDebug(_baseImageBinary(rectResult)));
//                                    cv::waitKey(-1);
                                    results.push_back(rectResult);
                                }
//                                else {
//                                    expandBoundingBox(_baseImageBinary, rectResult, centroid);
//                                    _thImg = closing(_baseImageBinary(rectResult));
//                                    centroid = feature::Centroid(_thImg);
//                                    geometricCenter = getGeometricCenter(_thImg);
//                                    std::cout << "distance/area again: " << centroid.distance(geometricCenter) / rectResult.size().area() << std::endl;
//                                    if (centroid.distance(geometricCenter) / rectResult.size().area() < 0.001)
//                                        results.push_back(rectResult);
//                                }
//                                cv::imshow("multipleThreshold", _baseImage(rectResult));
//                                cv::waitKey(-1);
                            }
                        }
                        i++;
                    }
                    break;

                default: //color::WHITE:
                    for (auto const &roi : ROIs) {
                        cv::Rect2i roiRect = structures::createRectFromRanges(roi.rangeRows, roi.rangeCols);
                        std::cout << " ::: rect: " << roiRect << std::endl;
                        processing::adjustBoundingBox(_baseImageBinary, roiRect);
                        _thImg = closing(_baseImageBinary(roiRect));
                        feature::BasicGeometricFeatures features(_thImg);
                        std::cout << "S: " << features.S << " ::: L: " << features.L << std::endl;
                        std::cout << "S/L: " << (double) features.S / (double) features.L << std::endl;
                        double s_l = (double) features.S / (double) features.L;
//                                                        cv::imshow("whiteFinalAnalysis", _baseImage(roiRect));
//                        cv::waitKey(-1);
                        if (s_l > 1.8 && s_l < 2.6) {
                            feature::Coefficients thCoeff(_thImg);
                            std::cout << "W9: " << thCoeff.W9 << std::endl;
                            if (thCoeff.W9 > 0.084 && thCoeff.W9 < 0.16) {
                                feature::Centroid centroid(_thImg);
                                cv::Point geometricCenter = getGeometricCenter(_thImg);
                                std::cout << "distance/area: " << centroid.distance(geometricCenter) / roiRect.size().area() << std::endl;
                                                        cv::imshow("whiteFinalAnalysis", _baseImage(roiRect));
                        cv::waitKey(-1);
                                if (centroid.distance(geometricCenter) / roiRect.size().area() < 0.001) {
//                                   cv::imshow("before", _baseImage(roiRect));
//                                    cv::waitKey(-1);
                                    expandBoundingBox(_baseImageBinary, roiRect, centroid);
//                                    cv::imshow("after", _baseImage(roiRect));
//                                    cv::waitKey(-1);
                                    results.push_back(roiRect);
                                } else {
                                    expandBoundingBox(_baseImageBinary, roiRect, centroid);
                                    _thImg = closing(_baseImageBinary(roiRect));
                                    centroid = feature::Centroid(_thImg);
                                    geometricCenter = getGeometricCenter(_thImg);
                                    std::cout << "distance/area again: " << centroid.distance(geometricCenter) / roiRect.size().area() << std::endl;
                                    if (centroid.distance(geometricCenter) / roiRect.size().area() < 0.001)
                                        results.push_back(roiRect);
                                }
                            }
                        }
                    }
                    break;
            }
            return structures::concatOverlappingRectangles(results);
        }

        std::vector<structures::RangesStruct> getROIs() {
            std::vector<structures::RangesStruct> possibleResults;
            for (const auto &level : _images) {
                cv::Mat_<uchar> _img = level.levelImage;
                std::vector<WindowStruct<uchar>> _windows;
                createWindows<uchar>(_img, _windows);

                for (const auto &winStruct: _windows) {
                    cv::Rect2i bBox = getBoundingBox(winStruct.window);
                    int area = bBox.size().area();
                    if (area <= 0) continue;
                    cv::Mat_<uchar> toAnalyze = cv::Mat_<uchar>(winStruct.window, bBox);
                    if (toAnalyze.rows == 0) continue;
//                        std::cout << " ::: M7: " << im.M7 <<  " ::: W4: " << coeff.W4 << std::endl;
//                        std::cout << "M1: " << im.M1 << " ::: M2: " << im.M2 << " ::: M3: " << im.M3 << " ::: M4: " << im.M4 << " ::: M5: " << im.M5 << " ::: M6: " << im.M6 << " ::: M7: " << im.M7 << std::endl;
                    if (couldBeROI(toAnalyze)) {
//                        std::cout << " ::: M7: " << im.M7 <<  " ::: W4: " << coeff.W4 << std::endl;
//                        std::cout << "M1: " << im.M1 << " ::: M2: " << im.M2 << " ::: M3: " << im.M3 << " ::: M4: " << im.M4 << " ::: M5: " << im.M5 << " ::: M6: " << im.M6 << " ::: M7: " << im.M7 << std::endl;
                        if (_pattern == color::WHITE) {
//                            cv::imshow("window_baseImage", _baseImage(cv::Range((winStruct.rangeRows.start + bBox.y) * (level.level + 1),
//                                                                      (winStruct.rangeRows.start + bBox.y + bBox.height) * (level.level + 1)), cv::Range((winStruct.rangeCols.start + bBox.x) * (level.level + 1),
//                                                                                                                                                         (winStruct.rangeCols.start + bBox.x + bBox.width) * (level.level + 1))));
//                            cv::waitKey(-1);

                        }
                        cv::Mat dbg = color::binaryDebug(toAnalyze);
                        possibleResults.push_back(structures::RangesStruct{
                                cv::Range((winStruct.rangeRows.start + bBox.y) * (level.level + 1),
                                          (winStruct.rangeRows.start + bBox.y + bBox.height) * (level.level + 1)),
                                cv::Range((winStruct.rangeCols.start + bBox.x) * (level.level + 1),
                                          (winStruct.rangeCols.start + bBox.x + bBox.width) * (level.level + 1))
                        });

                    }
//
                }
            }
            std::vector<structures::RangesStruct> result = structures::concatOverlappingRanges(possibleResults);
            return result;
        }

        bool couldBeROI(cv::Mat_<uchar>& image) {
            feature::Moments im = feature::Moments(image);
            feature::Coefficients coeff = feature::Coefficients(image);
            switch (_pattern) {
                case color::GREEN:
//                    if (coeff.W4 > 0.1 && coeff.W4 < 0.8 && im.M7 > 0.015 && im.M7 < 0.08) {
//                    std::cout << " ::: M7: " << im.M7 <<  " ::: W4: " << coeff.W4 << std::endl;
//                    std::cout << "M1: " << im.M1 << " ::: M2: " << im.M2 << " ::: M3: " << im.M3 << " ::: M4: " << im.M4 << " ::: M5: " << im.M5 << " ::: M6: " << im.M6 << " ::: M7: " << im.M7 << std::endl;
//                    }
                    return coeff.W4 > 0.1 && coeff.W4 < 0.8 && im.M7 > 0.015 && im.M7 < 0.08/* && im.M3 < 0.01*/;
                case color::WHITE:
                case color::WHITE_WIDER:
//                    if (coeff.W4 > 0.2 && coeff.W4 < 0.4 && im.M7 > 0.04 && im.M7 < 0.07 && im.M3 < 0.1) {
//                    std::cout << " ::: M7: " << im.M7 <<  " ::: W4: " << coeff.W4 << std::endl;
//                    std::cout << "M1: " << im.M1 << " ::: M2: " << im.M2 << " ::: M3: " << im.M3 << " ::: M4: " << im.M4 << " ::: M5: " << im.M5 << " ::: M6: " << im.M6 << " ::: M7: " << im.M7 << std::endl;
//                    }
                    return coeff.W4 > 0.2 && coeff.W4 < 0.4 && im.M7 > 0.04 && im.M7 < 0.07 && im.M3 < 0.1;
            }
        }
    };


}

#endif

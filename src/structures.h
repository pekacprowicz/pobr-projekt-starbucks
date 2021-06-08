#ifndef POBR_PROJEKT_STRUCTURES_H
#define POBR_PROJEKT_STRUCTURES_H

namespace structures {

    struct RangesStruct {
        cv::Range rangeRows;
        cv::Range rangeCols;

        RangesStruct& operator=(const RangesStruct &right) {
            this->rangeCols = right.rangeCols;
            this->rangeRows = right.rangeRows;
        }
    };

    cv::Rect2i createRectFromRanges(cv::Range rowRange, cv::Range colRange) {
        return cv::Rect2i(colRange.start, rowRange.start, colRange.size(), rowRange.size());
    }

    std::vector<RangesStruct> concatOverlappingRanges(std::vector<RangesStruct> ranges) {
        std::vector<RangesStruct> newRanges;
//        std::vector<RangesStruct> rangesToBeAddedAsNew;
        bool overlaps;
        for (auto const& range: ranges) {
//            rangesToBeAddedAsNew.clear();
            if (newRanges.empty())
                newRanges.push_back(range);
            else {
                overlaps = false;
                for (auto & checkedRange : newRanges) {
                    cv::Rect2i _rect = createRectFromRanges(range.rangeRows, range.rangeCols);
                    cv::Rect2i _checkedRect = createRectFromRanges(checkedRange.rangeRows, checkedRange.rangeCols);
                    if ((_rect & _checkedRect).area() > 0) {
                        overlaps = true;
                        int rangeStartRow = std::min(range.rangeRows.start, checkedRange.rangeRows.start);
                        int rangeEndRow = std::max(range.rangeRows.end, checkedRange.rangeRows.end);
                        int rangeStartCol = std::min(range.rangeCols.start, checkedRange.rangeCols.start);
                        int rangeEndCol = std::max(range.rangeCols.end, checkedRange.rangeCols.end);
                        checkedRange = RangesStruct{
                            cv::Range(
                                    rangeStartRow,
                                    rangeEndRow),
                            cv::Range(
                                    rangeStartCol,
                                    rangeEndCol)
                        };
                    }
                }
                if (!overlaps)
                    newRanges.push_back(range);
            }
        }
        return newRanges;
    }

    std::vector<cv::Rect2i> concatOverlappingRectangles(std::vector<cv::Rect2i> rectangles) {
        std::vector<cv::Rect2i> newRectangles;
        bool overlaps;
        for (auto const& rectangle: rectangles) {
//            rangesToBeAddedAsNew.clear();
            if (newRectangles.empty())
                newRectangles.push_back(rectangle);
            else {
                overlaps = false;
                for (auto & checkedRectangle : newRectangles) {
                    if ((rectangle & checkedRectangle).area() > 0) {
                        overlaps = true;
                        int x = std::min(rectangle.x, checkedRectangle.x);
                        int width = std::max(rectangle.x + rectangle.width, checkedRectangle.x + checkedRectangle.width) - x;
                        int y = std::min(rectangle.y, checkedRectangle.y);
                        int height = std::max(rectangle.y + rectangle.height, checkedRectangle.y + checkedRectangle.height) - y;
                        checkedRectangle = cv::Rect2i(x, y, width, height);
                    }
                }
                if (!overlaps)
                    newRectangles.push_back(rectangle);
            }
        }
        return newRectangles;
    }
}

#endif

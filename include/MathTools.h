#ifndef MATHTOOLS_H
#define MATHTOOLS_H

#include <vector>

class MathTools{
public:

    static double L2Norm
    (
        const std::vector<double>& res
    );

    static std::vector<double> vectorSub
    (
        const std::vector<double>& v1,
        const std::vector<double>& v2
    );

    static std::vector<double> vectorAdd
    (
        const std::vector<double>& v1,
        const std::vector<double>& v2
    );

};
#endif
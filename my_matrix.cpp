/*
 *  Copyright 2024 Christopher Vasquez
 */
#include "include/GMatrix.h"
#include <iostream>
#include <cmath>

GMatrix::GMatrix() : GMatrix(1, 0, 0, 0, 1, 0){}

GMatrix GMatrix::Translate(float tx, float ty){
    return GMatrix(1, 0, tx, 0, 1, ty);
}

GMatrix GMatrix::Scale(float sx, float sy){
    return GMatrix(sx, 0, 0, 0, sy, 0);  
}

GMatrix GMatrix::Rotate(float radians){
    const float s = sinf(radians);
    const float c = cosf(radians);
    return GMatrix(c, -s, 0, s, c, 0);
}

GMatrix GMatrix::Concat(const GMatrix& a, const GMatrix& b){
    return GMatrix((a[0] * b[0]) + (a[2] * b[1]),(a[0] * b[2]) + (a[2] * b[3]),(a[0] * b[4]) + (a[2] * b[5]) + a[4],(a[1] * b[0]) + (a[3] * b[1]),(a[1] * b[2]) + (a[3] * b[3]),(a[1] * b[4]) + (a[3] * b[5]) + a[5]);
}

static float dcross(double a, double b, double c, double d){
    return static_cast<float>(a * b - c * d);
}

nonstd::optional<GMatrix> GMatrix::invert() const{
    float det = dcross(fMat[0], fMat[3], fMat[1], fMat[2]);
    if(0 == det){
        return {};
    }
    float idet = 1/det;

    float a =  fMat[3] * idet;
    float c = -fMat[2] * idet;
    float e = dcross(fMat[2], fMat[5], fMat[3], fMat[4]) * idet;
    float b = -fMat[1] * idet;
    float d =  fMat[0] * idet;
    float f = dcross(fMat[1], fMat[4], fMat[0], fMat[5]) * idet;

    return GMatrix(a, c, e, b, d, f);
    

}

void GMatrix::mapPoints(GPoint dst[], const GPoint src[], int count) const{ 
    const auto e0 = this->e0();
    const auto e1 = this->e1();
    const auto origin = this->origin();
    for (int i = 0; i < count; ++i) {
        dst[i] = e0 * src[i].x + e1 * src[i].y + origin;
    }
}
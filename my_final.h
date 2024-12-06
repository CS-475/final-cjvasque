/*
 *  Copyright 2024 Christopher Vasquez
 */

#include "include/GRect.h"
#include "include/GColor.h"
#include "include/GPathBuilder.h"
#include "GEdge.h"
#include "include/GBitmap.h"
#include <deque>
#include "include/GFinal.h"

class MyFinal : public GFinal{
public:
    MyFinal(){}

    std::shared_ptr<GShader> createVoronoiShader(const GPoint points[],
                                                         const GColor colors[],
                                                         int count) override;
    
    std::shared_ptr<GShader> createLinearPosGradient(GPoint p0, GPoint p1,
                                                             const GColor colors[],
                                                             const float pos[],
                                                             int count) override;

private:
};
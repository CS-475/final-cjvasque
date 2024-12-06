#ifndef GEdge_DEFINED
#define GEdge_DEFINED

#include "include/GPoint.h"
#include <algorithm>
#include <iostream>

struct GEdge{
        float m, b, curr_X;
        int top, bottom, w;
    };

GEdge makeEdge(GPoint p1, GPoint p2){
    GEdge e;
    e.m = (p2.x - p1.x)/(p2.y - p1.y);
    e.b = p1.x - (p1.y * e.m);
    e.top = GRoundToInt(std::min(p1.y,p2.y));
    e.bottom = GRoundToInt(std::max(p1.y,p2.y));
    e.w = 0;
    e.curr_X = 0.0f;
    return e;
}

int eval(GEdge e, int y){
    return GRoundToInt((float)((y + 0.5)*e.m) + e.b);
}

bool valid(GEdge e) {
            return e.top < e.bottom;
        }

bool isValid(GEdge e, float y){
    return (float)e.top <= y && y <= (float)e.bottom;
}

bool lessThan(const GEdge& a, const GEdge& b){
    if(a.top == b.top){
        int x0 = eval(a, a.top);
        int x1 = eval(b, b.top);
        return x0 < x1;
    }
    else{
        return a.top < b.top;
    }
}

#endif
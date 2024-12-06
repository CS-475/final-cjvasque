/*
 *  Copyright 2024 Christopher Vasquez
 */
#include "include/GPath.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <deque>

GRect GPath::bounds() const{
    float minX = 0;
    float minY = 0;
    float maxX = 0;
    float maxY = 0;

    if(fPts.empty()){
        return GRect::LTRB(minX, minY, maxX, maxY);
    }

    GPoint top_left;
    GPoint bottom_right;

    GPoint pts[GPath::kMaxNextPoints];
    GPath::Edger edg(*this);
    std::deque<GPoint> ptList;

    while(auto v = edg.next(pts)){
        switch(v.value()){
            case GPathVerb::kLine:
            minX = std::min(pts[0].x, pts[1].x);
            maxX = std::max(pts[0].x, pts[1].x);
            minY = std::min(pts[0].y, pts[1].y);
            maxY = std::max(pts[0].y, pts[1].y);
            top_left = {minX, minY};
            bottom_right = {maxX, maxY};
            ptList.push_front(top_left);
            ptList.push_front(bottom_right);
            break;

            case GPathVerb::kQuad:
            minX = std::min(pts[0].x,pts[2].x);
            minY = std::min(pts[0].y,pts[2].y);
            maxX = std::max(pts[0].x,pts[2].x);
            maxY = std::max(pts[0].y,pts[2].y);
            if(pts[1].x < minX || pts[1].x > maxX || pts[1].y < minY || pts[1].y > maxY){
                float t_x = std::clamp((pts[0].x - pts[1].x)/(pts[0].x - (2*pts[1].x) + pts[2].x), 0.0f, 1.0f);
                float t_y = std::clamp((pts[0].y - pts[1].y)/(pts[0].y - (2*pts[1].y) + pts[2].y), 0.0f, 1.0f);

                float s_x = (1 - t_x);
                float s_y = (1 - t_y);

                float q_x = s_x*s_x*pts[0].x + 2*s_x*t_x*pts[1].x + t_x*t_x*pts[2].x;
                float q_y = s_y*s_y*pts[0].y + 2*s_y*t_y*pts[1].y + t_y*t_y*pts[2].y;

                minX = std::min(minX, q_x);
                minY = std::min(minY, q_y);
                maxX = std::max(maxX, q_x);
                maxY = std::max(maxY, q_y);

                top_left = {minX, minY};
                bottom_right = {maxX, maxY};
                ptList.push_front(top_left);
                ptList.push_front(bottom_right);

            }
            else{
                top_left = {minX, minY};
                bottom_right = {maxX, maxY};
                ptList.push_front(top_left);
                ptList.push_front(bottom_right);
            }
            break;

            case GPathVerb::kCubic:
            minX = std::min(pts[0].x,pts[3].x);
            minY = std::min(pts[0].y,pts[3].y);
            maxX = std::max(pts[0].x,pts[3].x);
            maxY = std::max(pts[0].y,pts[3].y);

            float c_x = -1*pts[0].x + pts[1].x;
            float c_y = -1*pts[0].y + pts[1].y;

            float b_x = pts[0].x - 2*pts[1].x + pts[2].x;
            float b_y = pts[0].y - 2*pts[1].y + pts[2].y;

            float a_x = -1*pts[0].x + 3*pts[1].x - 3*pts[2].x + pts[3].x;
            float a_y = -1*pts[0].y + 3*pts[1].y - 3*pts[2].y + pts[3].y;

            float h_x = (b_x*b_x)  - (a_x*c_x);
            float h_y = (b_y*b_y)  - (a_y*c_y);

            if(h_x > 0){
                h_x = (float)sqrt(h_x);
                float t = (-1*b_x - h_x)/a_x;
                if(a_x == 0){
                    t = 0.5f;
                    float s = 1 - t;
                    float q = s*s*s*pts[0].x + 3*s*s*t*pts[1].x + 3*s*t*t*pts[2].x + t*t*t*pts[3].x;
                    minX = std::min(minX,q);
                    maxX = std::max(maxX,q);
                }
                else if(t > 0 && t < 1){
                    float s = 1 - t;
                    float q = s*s*s*pts[0].x + 3*s*s*t*pts[1].x + 3*s*t*t*pts[2].x + t*t*t*pts[3].x;
                    minX = std::min(minX,q);
                    maxX = std::max(maxX,q);
                }
                t = (-1*b_x + h_x)/a_x;
                if(a_x = 0){
                    t = 0.5f;
                    float s = 1 - t;
                    float q = s*s*s*pts[0].x + 3*s*s*t*pts[1].x + 3*s*t*t*pts[2].x + t*t*t*pts[3].x;
                    minX = std::min(minX,q);
                    maxX = std::max(maxX,q);
                }
                else if(t > 0 && t < 1){
                    float s = 1 - t;
                    float q = s*s*s*pts[0].x + 3*s*s*t*pts[1].x + 3*s*t*t*pts[2].x + t*t*t*pts[3].x;
                    minX = std::min(minX,q);
                    maxX = std::max(maxX,q);
                }
            }

            if(h_y > 0){
                h_y = (float)sqrt(h_y);
                float t = (-1*b_y - h_y)/a_y;
                if(a_y == 0){
                    t = 0.5f;
                    float s = 1 - t;
                    float q = s*s*s*pts[0].y + 3*s*s*t*pts[1].y + 3*s*t*t*pts[2].y + t*t*t*pts[3].y;
                    minY = std::min(minY,q);
                    maxY = std::max(maxY,q);
                }
                else if(t > 0 && t < 1){
                    float s = 1 - t;
                    float q = s*s*s*pts[0].y + 3*s*s*t*pts[1].y + 3*s*t*t*pts[2].y + t*t*t*pts[3].y;
                    minY = std::min(minY,q);
                    maxY = std::max(maxY,q);
                }
                t = (-1*b_y + h_y)/a_y;
                if(a_y == 0){
                    t = 0.5f;
                    float s = 1 - t;
                    float q = s*s*s*pts[0].y + 3*s*s*t*pts[1].y + 3*s*t*t*pts[2].y + t*t*t*pts[3].y;
                    minY = std::min(minY,q);
                    maxY = std::max(maxY,q);
                }
                else if(t > 0 && t < 1){
                    float s = 1 - t;
                    float q = s*s*s*pts[0].y + 3*s*s*t*pts[1].y + 3*s*t*t*pts[2].y + t*t*t*pts[3].y;
                    minY = std::min(minY,q);
                    maxY = std::max(maxY,q);
                }
            }

            top_left = {minX, minY};
            bottom_right = {maxX, maxY};
            ptList.push_front(top_left);
            ptList.push_front(bottom_right);
            break;
        }
    }


    if(!ptList.empty()){
        for(int i = 0; i < ptList.size(); i++){
            minX = std::min(minX, ptList.at(i).x);
            maxX = std::max(maxX, ptList.at(i).x);
            minY = std::min(minY, ptList.at(i).y);
            maxY = std::max(maxY, ptList.at(i).y);
        }
    }

    return GRect::LTRB(minX, minY, maxX, maxY);
}

void GPath::ChopQuadAt(const GPoint src[3], GPoint dst[5], float t){
    
    float f_tx = (1-t) * src[0].x + t * src[1].x;
    float f_ty = (1-t) * src[0].y + t * src[1].y;
    float g_tx = (1-t) * src[1].x + t * src[2].x;
    float g_ty = (1-t) * src[1].y + t * src[2].y;

    GPoint f; f.x = f_tx; f.y = f_ty;
    GPoint g; g.x = g_tx; g.y = g_ty;
    
    float h_tx = (1-t) * f.x + t * g.x;
    float h_ty = (1-t) * f.y + t * g.y;
    GPoint h; h.x = h_tx; h.y = h_ty;

    dst[0] = src[0];
    dst[1] = f;
    dst[2] = h;
    dst[3] = g;
    dst[4] = src[2];
}

void GPath::ChopCubicAt(const GPoint src[4], GPoint dst[7], float t){
    
    float ab_tx = (1-t) * src[0].x + t * src[1].x;
    float ab_ty = (1-t) * src[0].y + t * src[1].y;
    GPoint ab; ab.x = ab_tx; ab.y = ab_ty;

    float bc_tx = (1-t) * src[1].x + t * src[2].x;
    float bc_ty = (1-t) * src[1].y + t * src[2].y;
    GPoint bc; bc.x = bc_tx; bc.y = bc_ty;

    float cd_tx = (1-t) * src[2].x + t * src[3].x;
    float cd_ty = (1-t) * src[2].y + t * src[3].y;
    GPoint cd; cd.x = cd_tx; cd.y = cd_ty;

    float ab_bc_tx = (1-t) * ab.x + t * bc.x;
    float ab_bc_ty = (1-t) * ab.y + t * bc.y;
    GPoint ab_bc; ab_bc.x = ab_bc_tx; ab_bc.y = ab_bc_ty;

    float bc_cd_tx = (1-t) * bc.x + t * cd.x;
    float bc_cd_ty = (1-t) * bc.y + t * cd.y;
    GPoint bc_cd; bc_cd.x = bc_cd_tx; bc_cd.y = bc_cd_ty;

    float abcd_tx = (1-t) * ab_bc.x + t * bc_cd.x;
    float abcd_ty = (1-t) * ab_bc.y + t * bc_cd.y;
    GPoint abcd; abcd.x = abcd_tx; abcd.y = abcd_ty;

    dst[0] = src[0];
    dst[1] = ab;
    dst[2] = ab_bc;
    dst[3] = abcd;
    dst[4] = bc_cd;
    dst[5] = cd;
    dst[6] = src[3];

}
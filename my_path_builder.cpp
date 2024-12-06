/*
 *  Copyright 2024 Christopher Vasquez
 */
#include "include/GPathBuilder.h"
#include <iostream>

void GPathBuilder::addRect(const GRect& rect, GPathDirection d){
    switch(d){
        case GPathDirection::kCW: 
        this->moveTo({rect.left, rect.top});
        this->lineTo({rect.right, rect.top});
        this->lineTo({rect.right, rect.bottom});
        this->lineTo({rect.left, rect.bottom});
        break;
        case GPathDirection::kCCW: 
        this->moveTo({rect.left, rect.top});
        this->lineTo({rect.left, rect.bottom});
        this->lineTo({rect.right, rect.bottom});
        this->lineTo({rect.right, rect.top});
        break;
    }
}

void GPathBuilder::addPolygon(const GPoint pts[], int count){
    this->moveTo(pts[0]);
    for(int i = 1; i < count; i++){
        this->lineTo(pts[i]);
    }
}

void GPathBuilder::addCircle(GPoint center, float radius, GPathDirection d){
    GMatrix mx = GMatrix::Translate(center.x, center.y) * GMatrix::Scale(radius, radius);

    float tan_pi = (float)(sqrt(2) - 1);
    float pi_4 = (float)(sqrt(2)/2);

    GPoint pt_1 = {(mx[0] * (1) + mx[2] * (0) + mx[4]), ((mx[1] * (1) + mx[3] * (0) + mx[5]))};

    GPoint pt_2 = {(mx[0] * (1) + mx[2] * (tan_pi) + mx[4]), ((mx[1] * (1) + mx[3] * (tan_pi) + mx[5]))};

    GPoint pt_3 = {(mx[0] * (pi_4) + mx[2] * (pi_4) + mx[4]), ((mx[1] * (pi_4) + mx[3] * (pi_4) + mx[5]))};

    GPoint pt_4 = {(mx[0] * (tan_pi) + mx[2] * (1) + mx[4]), ((mx[1] * (tan_pi) + mx[3] * (1) + mx[5]))};

    GPoint pt_5 = {(mx[0] * (0) + mx[2] * (1) + mx[4]), ((mx[1] * (0) + mx[3] * (1) + mx[5]))};

    GPoint pt_6 = {(mx[0] * (-tan_pi) + mx[2] * (1) + mx[4]), ((mx[1] * (-tan_pi) + mx[3] * (1) + mx[5]))};

    GPoint pt_7 = {(mx[0] * (-pi_4) + mx[2] * (pi_4) + mx[4]), ((mx[1] * (-pi_4) + mx[3] * (pi_4) + mx[5]))};

    GPoint pt_8 = {(mx[0] * (-1) + mx[2] * (tan_pi) + mx[4]), ((mx[1] * (-1) + mx[3] * (tan_pi) + mx[5]))};

    GPoint pt_9 = {(mx[0] * (-1) + mx[2] * (0) + mx[4]), ((mx[1] * (-1) + mx[3] * (0) + mx[5]))};

    GPoint pt_10 = {(mx[0] * (-1) + mx[2] * (-tan_pi) + mx[4]), ((mx[1] * (-1) + mx[3] * (-tan_pi) + mx[5]))};

    GPoint pt_11 = {(mx[0] * (-pi_4) + mx[2] * (-pi_4) + mx[4]), ((mx[1] * (-pi_4) + mx[3] * (-pi_4) + mx[5]))};

    GPoint pt_12 = {(mx[0] * (-tan_pi) + mx[2] * (-1) + mx[4]), ((mx[1] * (-tan_pi) + mx[3] * (-1) + mx[5]))};

    GPoint pt_13 = {(mx[0] * (0) + mx[2] * (-1) + mx[4]), ((mx[1] * (0) + mx[3] * (-1) + mx[5]))};

    GPoint pt_14 = {(mx[0] * (tan_pi) + mx[2] * (-1) + mx[4]), ((mx[1] * (tan_pi) + mx[3] * (-1) + mx[5]))};

    GPoint pt_15 = {(mx[0] * (pi_4) + mx[2] * (-pi_4) + mx[4]), ((mx[1] * (pi_4) + mx[3] * (-pi_4) + mx[5]))};

    GPoint pt_16 = {(mx[0] * (1) + mx[2] * (-tan_pi) + mx[4]), ((mx[1] * (1) + mx[3] * (-tan_pi) + mx[5]))};

    switch(d){
        case GPathDirection::kCCW:{
            this->moveTo(pt_1);
            this->quadTo(pt_16, pt_15);
            this->quadTo(pt_14, pt_13);
            this->quadTo(pt_12, pt_11);
            this->quadTo(pt_10, pt_9);
            this->quadTo(pt_8, pt_7);
            this->quadTo(pt_6, pt_5);
            this->quadTo(pt_4, pt_3);
            this->quadTo(pt_2, pt_1);
            break;
        }
        

        case GPathDirection::kCW:{
            this->moveTo(pt_1);
            this->quadTo(pt_2, pt_3);
            this->quadTo(pt_4, pt_5);
            this->quadTo(pt_6, pt_7);
            this->quadTo(pt_8, pt_9);
            this->quadTo(pt_10, pt_11);
            this->quadTo(pt_12, pt_13);
            this->quadTo(pt_14, pt_15);
            this->quadTo(pt_16, pt_1);
            break;
        }
        
    }


}
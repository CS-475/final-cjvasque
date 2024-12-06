/*
 *  Copyright 2024 Christopher Vasquez
 */

#ifndef _g_starter_canvas_h_
#define _g_starter_canvas_h_

#include "include/GCanvas.h"
#include "include/GRect.h"
#include "include/GColor.h"
#include "include/GPath.h"
#include "include/GPathBuilder.h"
#include "GEdge.h"
#include "include/GShader.h"
#include "include/GBitmap.h"
#include <deque>

class MyCanvas : public GCanvas {
public:
    MyCanvas(const GBitmap& device) : fDevice(device) {
        ctm = GMatrix();
        ctmList.push_front(ctm);
        haveSaved = false;
    }

    void save();

    void restore();

    void concat(const GMatrix& matrix);

    void clear(const GColor& color) override;

    void translate(float x, float y);

    void scale(float x, float y);

    void rotate(float radians);

    void fillRect(const GRect& rect, const GColor& color);

    void drawRect(const GRect& rect, const GPaint& paint) override;

    void drawConvexPolygon(const GPoint pts[] , int count, const GPaint&) override;

    void drawPath(const GPath& path, const GPaint& paint) override;

    void drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[], int count, const int indices[], const GPaint& paint) override;

    void drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4], int level, const GPaint& paint) override;

    void drawTriangleWithColor(const GPoint verts[3], const GColor colors[3], const GPaint& paint);

    void drawTriangleWithTex(const GPoint verts[3], const GPoint tex[3], const GPaint& paint);

    void drawTriangleWithBoth(const GPoint verts[3], const GColor colors[3], const GPoint tex[3], const GPaint& paint);
    
    void blit(int left, int y, int width, const GPaint& paint);

    //void remove_edge(std::deque<GEdge>& deque, const GEdge& value);

private:
    // Note: we store a copy of the bitmap
    const GBitmap fDevice;

    // Add whatever other fields you need
    GMatrix ctm;

    std::deque<GMatrix> ctmList;

    std::shared_ptr<GShader> sh;

    int evaluate(GEdge a, float y);

    float evaluate_precise(GEdge a, float y);

    GColor lerp_color(const GColor& a, const GColor& b, float t);

    GPoint lerp_line(const GPoint& a, const GPoint& b, float t);

    GPoint lerp_quad(const GPoint& a, const GPoint& b, const GPoint& c, float t);

    GPoint lerp_cubic(const GPoint& a, const GPoint& b, const GPoint& c, const GPoint& d, float t);

    bool haveSaved;

    GPixel blend(GPixel& src, GPixel& dst, GBlendMode mode);

    GPixel clear_mode(GPixel& src, GPixel& dst);

    GPixel src_mode(GPixel& src, GPixel& dst);

    GPixel dst_mode(GPixel& src, GPixel& dst);

    GPixel srcOver_mode(GPixel& src, GPixel& dst);

    GPixel dstOver_mode(GPixel& src, GPixel& dst);

    GPixel srcIn_mode(GPixel& src, GPixel& dst);

    GPixel dstIn_mode(GPixel& src, GPixel& dst);

    GPixel srcOut_mode(GPixel& src, GPixel& dst);

    GPixel dstOut_mode(GPixel& src, GPixel& dst);

    GPixel srcAtop_mode(GPixel& src, GPixel& dst);

    GPixel dstAtop_mode(GPixel& src, GPixel& dst);

    GPixel xor_mode(GPixel& src, GPixel& dst);

};

#endif
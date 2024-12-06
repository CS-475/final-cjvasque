/*
 *  Copyright 2024 Christopher Vasquez
 */

#include "my_final.h"
#include "my_header.h"
#include <cmath>
#include <optional>
#include <algorithm>
#include <iostream>
#include <list>
#include <deque>
#include <vector>
#include <stdexcept>
#include <array>

class MyLinearPosGradient : public GShader{
    public:
    MyLinearPosGradient(GPoint p0, GPoint p1, const GColor colors[], const float pos[], int count) : p0(p0), p1(p1), count(count), ctm(GMatrix()){
        float dx = p1.x - p0.x;
        float dy = p1.y - p0.y;
        lgm = GMatrix(dx, -dy, p0.x, dy, dx, p0.y);
        k = count - 1;
        c = new GColor[count];
        std::copy(colors, colors + count, c);
        p = new float[count];
        std::copy(pos, pos + count, p);
    }

    bool isOpaque(){
        for(int i = 0; i < count; i++){
            if((c+i)->a != 1.0f){
                return false;
            }
        }
        return true;
    }
    bool setContext(const GMatrix& ctm){
        if(!ctm.invert().has_value()){
            return false;
        }
        else{
            this->ctm = ctm;
            imx = lgm.invert().value() * ctm.invert().value();
            return true;
        }
    }

    void shadeRow(int x, int y, int count, GPixel row[]){
        GColor color;
        for(int i = 0; i < count; i++){
            float px = (float)(imx.value()[0] * (x + i + 0.5) + imx.value()[2] * (y + 0.5) + imx.value()[4]);
            if(px < 0.0f){
                px = 0.0f;
            }
            if(px > 1.0f){
                px = 1.0f;
            }
            if((px == p[0] || this->count == 1)){
                row[i] = GPixel_PackARGB(GRoundToInt(c[0].a*255.0f),GRoundToInt(c[0].r*(c[0].a)*255.0f),GRoundToInt(c[0].g*(c[0].a)*255.0f),GRoundToInt(c[0].b*(c[0].a)*255.0f));
            }
            else if(px == p[k]){
                row[i] = GPixel_PackARGB(GRoundToInt(c[k].a*255.0f),GRoundToInt(c[k].r*(c[k].a)*255.0f),GRoundToInt(c[k].g*(c[k].a)*255.0f),GRoundToInt(c[k].b*(c[k].a)*255.0f));
            }
            else{
                float start_pos = p[0];
                float end_pos = p[k];
                float t = 0;
                int start_index = 0;
                int end_index = k;
                for(int i = 0; i < count; i++){
                    if(px <= p[i]){
                        end_pos = p[i];
                        start_pos = p[i - 1];
                        start_index = i - 1;
                        end_index = i;
                        t = (px - start_pos)/(-start_pos + end_pos);
                        break;
                    }
                }
                color.a = c[start_index].a + (c[end_index].a - c[start_index].a) * t;
                color.r = c[start_index].r + (c[end_index].r - c[start_index].r) * t;
                color.g = c[start_index].g + (c[end_index].g - c[start_index].g) * t;
                color.b = c[start_index].b + (c[end_index].b - c[start_index].b) * t;
                row[i] = GPixel_PackARGB(GRoundToInt(color.a*255.0f),GRoundToInt(color.r*color.a*255.0f),GRoundToInt(color.g*color.a*255.0f),GRoundToInt(color.b*color.a*255.0f));

            }
        }
    }

    private:
    GPoint p0;
    GPoint p1;
    int count;
    GMatrix ctm;
    GMatrix lgm;
    GMatrix scmx;
    nonstd::optional<GMatrix> imx;
    int k;
    GColor* c;
    float* p;
};
class MyVoronoiShader : public GShader{
    public:
    MyVoronoiShader(const GPoint points[], const GColor colors[], int count) : count(count){
        c = new GColor[count];
        std::copy(colors, colors + count, c);
        p = new GPoint[count];
        std::copy(points, points + count, p);
    }

    bool isOpaque(){
        for(int i = 0; i < count; i++){
            if((c+i)->a != 1.0f){
                    return false;
            }
        }
        return true;
    }

    bool setContext(const GMatrix& ctm){
        if(!ctm.invert().has_value()){
                return false;
            }
            else{
                this->ctm = ctm;
                imx = ctm.invert().value();
                return true;
            }
    }

    void shadeRow(int x, int y, int count, GPixel row[]){
        for(int i = 0; i < count; i++){
            GPoint closest_point = p[0];
            int closest_index = 0;
            float px = (float)(imx.value()[0] * (x + i + 0.5) + imx.value()[2] * (y + 0.5) + imx.value()[4]);
            float py = (float)(imx.value()[1] * (x + i + 0.5) + imx.value()[3] * (y + 0.5) + imx.value()[5]);
            float min_distance = (float)sqrt(pow((px - p[0].x), 2) + pow((py - p[0].y), 2));
            for(int j = 0; j < this->count; j++){
                float temp_distance = min_distance;
                min_distance = std::min(min_distance, (float)sqrt(pow((px - p[j].x), 2) + pow((py - p[j].y), 2)));
                if(temp_distance != min_distance){
                    closest_index = j;
                }
            }
            row[i] = GPixel_PackARGB(GRoundToInt(c[closest_index].a*255.0f),GRoundToInt(c[closest_index].r*c[closest_index].a*255.0f),GRoundToInt(c[closest_index].g*c[closest_index].a*255.0f),GRoundToInt(c[closest_index].b*c[closest_index].a*255.0f));
        }
    }

    private:
    GMatrix ctm;
    nonstd::optional<GMatrix> imx;
    GColor* c;
    GPoint* p;
    int count;

};

class MyTriColorShader : public GShader{
    public:
    MyTriColorShader(GPoint p0, GPoint p1, GPoint p2, 
     const GColor& c0, const GColor& c1, const GColor& c2) : p0(p0), p1(p1), p2(p2), c0(c0), c1(c1), c2(c2){
        float u_x = p1.x - p0.x;
        float u_y = p1.y - p0.y;
        float v_x = p2.x = p0.x;
        float v_y = p2.y - p0.y;
        m = GMatrix(u_x, v_x, p0.x, u_y, v_y, p0.y);
    }

    bool isOpaque(){
        return c0.a == 1 && c1.a == 1 && c2.a == 1;
    }

    bool setContext(const GMatrix& ctm){
        if(!ctm.invert().has_value()){
                return false;
            }
            else{
                this->ctm = ctm;
                imx = m.invert().value() * ctm.invert().value();
                return true;
            }
    }

    void shadeRow(int x, int y, int count, GPixel row[]){
        GColor color;
        for(int i = 0; i < count; i++){
            //std::cout << "c0.a: " << c0.a << std::endl;
            //std::cout << "c0.r: " << c0.r << std::endl;
            //std::cout << "c0.g: " << c0.g << std::endl;
            //std::cout << "c0.b: " << c0.b << std::endl;
            //std::cout << "c1.a: " << c1.a << std::endl;
            //std::cout << "c1.r: " << c1.r << std::endl;
            //std::cout << "c1.g: " << c1.g << std::endl;
            //std::cout << "c1.b: " << c1.b << std::endl;
            //std::cout << "c2.a: " << c2.a << std::endl;
            //std::cout << "c2.r: " << c2.r << std::endl;
            //std::cout << "c2.g: " << c2.g << std::endl;
            //std::cout << "c2.b: " << c2.b << std::endl;
            float px = (float)(imx.value()[0] * (x + i + 0.5) + imx.value()[2] * (y + 0.5) + imx.value()[4]);
            float py = (float)(imx.value()[1] * (x + i + 0.5) + imx.value()[3] * (y + 0.5) + imx.value()[5]);
            if(px < 0){
                px = 0;
            }
            if(px > 1){
                px = 1;
            }
            if(py < 0){
                py = 0;
            }
           if(py > 1){
                py = 1;
            }
            //std::cout << "px: " << px << std::endl;
            //std::cout << "py: " << py << std::endl;

            color.a = (abs)(px * c1.a + py * c2.a + (1 - px - py) * c0.a);
            color.r = (abs)(px * c1.r + py * c2.r + (1 - px - py) * c0.r);
            color.g = (abs)(px * c1.g + py * c2.g + (1 - px - py) * c0.g);
            color.b = (abs)(px * c1.b + py * c2.b + (1 - px - py) * c0.b);
            //std::cout << "color a value: " << color.a << std::endl;
            //std::cout << "color r value: " << color.r << std::endl;
            //std::cout << "color g value: " << color.g << std::endl;
            //std::cout << "color b value: " << color.b << std::endl;

            row[i] = GPixel_PackARGB(GRoundToInt(color.a*255.0f),GRoundToInt(color.r*color.a*255.0f),GRoundToInt(color.g*color.a*255.0f),GRoundToInt(color.b*color.a*255.0f));
        }
    }

    private:
    GPoint p0;
    GPoint p1;
    GPoint p2;
    GColor c0;
    GColor c1;
    GColor c2;
    GMatrix m;
    GMatrix ctm;
    nonstd::optional<GMatrix> imx;
};

class ProxyShader : public GShader {
public:
    ProxyShader(GShader* shader, const GMatrix& extraTransform)
        : fRealShader(shader), fExtraTransform(extraTransform) {}

    bool isOpaque() override { return fRealShader->isOpaque(); }

    bool setContext(const GMatrix& ctm) override {
        return fRealShader->setContext(ctm * fExtraTransform);
    }
    
    void shadeRow(int x, int y, int count, GPixel row[]) override {
        fRealShader->shadeRow(x, y, count, row);
    }
    private:
    GShader* fRealShader;
    GMatrix  fExtraTransform;
};

class ComposShader : public GShader {
public:
    ComposShader(std::shared_ptr<GShader> shader0, std::shared_ptr<GShader> shader1) : fShader_0(shader0), fShader_1(shader1){}

    bool isOpaque() override{
        return (fShader_0->isOpaque() && fShader_1->isOpaque());
    }

    bool setContext(const GMatrix& ctm) override{
        return (fShader_0->setContext(ctm) && fShader_1->setContext(ctm));
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override{
        GPixel s0[count];
        GPixel s1[count];
        fShader_0->shadeRow(x, y, count, s0);
        fShader_1->shadeRow(x, y, count, s1);
        for(int i = 0; i < count; i++){
            //std::cout << "count: " << count << std::endl;
            //std::cout << "s0[" << i << "].a: " << GPixel_GetA(s0[i]) << std::endl;
            //std::cout << "s0[" << i << "].r: " << GPixel_GetR(s0[i]) << std::endl;
            //std::cout << "s0[" << i << "].g: " << GPixel_GetG(s0[i]) << std::endl;
            //std::cout << "s0[" << i << "].b: " << GPixel_GetB(s0[i]) << std::endl;
            //std::cout << "s1[" << i << "].a: " << GPixel_GetA(s1[i]) << std::endl;
            //std::cout << "s1[" << i << "].r: " << GPixel_GetR(s1[i]) << std::endl;
            //std::cout << "s1[" << i << "].g: " << GPixel_GetG(s1[i]) << std::endl;
            //std::cout << "s1[" << i << "].b: " << GPixel_GetB(s1[i]) << std::endl;

            int resultA = (int)((GPixel_GetA(s0[i]) * GPixel_GetA(s1[i]))/255.0f);
            int resultR = (int)((GPixel_GetR(s0[i]) * GPixel_GetR(s1[i]))/255.0f);
            int resultG = (int)((GPixel_GetG(s0[i]) * GPixel_GetG(s1[i]))/255.0f);
            int resultB = (int)((GPixel_GetB(s0[i]) * GPixel_GetB(s1[i]))/255.0f);
            //std::cout << "resultA: " << resultA << std::endl;
            //std::cout << "resultR: " << resultR << std::endl;
            //std::cout << "resultG: " << resultG << std::endl;
            //std::cout << "resultB: " << resultB << std::endl;
            row[i] = GPixel_PackARGB(resultA, resultR, resultG, resultB);
        }

    }

    private:
    std::shared_ptr<GShader> fShader_0;
    std::shared_ptr<GShader> fShader_1;

};


std::shared_ptr<GShader> GCreateMyTriColorShader(GPoint p0, GPoint p1, GPoint p2, const GColor& c0, const GColor& c1, const GColor& c2){
    return std::make_shared<MyTriColorShader>(p0, p1, p2, c0, c1, c2);
}

std::shared_ptr<GShader> GCreateProxyShader(GShader* shader, const GMatrix& extraTransform){
    return std::make_shared<ProxyShader>(shader, extraTransform);
}

std::shared_ptr<GShader> GCreateComposShader(std::shared_ptr<GShader> shader0, std::shared_ptr<GShader> shader1){
    return std::make_shared<ComposShader>(shader0, shader1);
}

std::shared_ptr<GShader> MyFinal::createVoronoiShader(const GPoint points[], const GColor colors[], int count){
    return std::make_shared<MyVoronoiShader>(points, colors, count);
}

std::shared_ptr<GShader> MyFinal::createLinearPosGradient(GPoint p0, GPoint p1, const GColor colors[], const float pos[], int count){
    return std::make_shared<MyLinearPosGradient>(p0, p1, colors, pos, count);
}

void MyCanvas::save(){
    this->ctmList.push_front(this->ctm);
    haveSaved = true;
}

void MyCanvas::restore(){
    if(!haveSaved){
        throw std::runtime_error("No previous call to save");
    }
    else{
        if(this->ctmList.size() > 1){
            this->ctmList.pop_front();
            this->ctm = this->ctmList.front();
        }
    }
}

void MyCanvas::concat(const GMatrix& matrix){
    this->ctm = this->ctm * matrix;
    this->ctmList.front() = this->ctm;
}

void MyCanvas::clear(const GColor& color) {
    // your code here
    GPixel pixel_color = GPixel_PackARGB(GRoundToInt(color.a*255.0f),GRoundToInt(color.r*color.a*255.0f),GRoundToInt(color.g*color.a*255.0f),GRoundToInt(color.b*color.a*255.0f));
    int width = fDevice.width();
    int height = fDevice.height();

    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
        *fDevice.getAddr(x,y) = pixel_color;
        }
  }
}

void MyCanvas::fillRect(const GRect& rect, const GColor& color) {
    // your code here
    this->drawRect(rect, GPaint(color));
}

void MyCanvas::drawRect(const GRect& rect, const GPaint& paint){
    GPoint pts[4];

    pts[0].x = rect.x(); pts[0].y = rect.y();
    pts[1].x = rect.x(); pts[1].y =  rect.y() + rect.height();
    pts[2].x = rect.x() + rect.width(); pts[2].y = rect.y() + rect.height();
    pts[3].x = rect.x() + rect.width(); pts[3].y = rect.y();
    
    drawConvexPolygon(pts, 4, paint);
}

void MyCanvas::drawConvexPolygon(const GPoint pts[], int count, const GPaint& paint){
    GPoint dst[count];
    if(ctm.invert().has_value()){
        ctm.mapPoints(dst,pts,count);
    }
    else{return;}
    std::deque<GEdge> edgeList;
    int edgeCount = 0;
    for(int i = 0; i < count; i++){
        int k = (i + 1) % count;
        GEdge e;
        e = makeEdge(dst[i],dst[k]);
        if(valid(e)){
            edgeList.push_front(e);
            edgeCount++;
        }
    }
    std::sort(edgeList.begin(),edgeList.end(),lessThan);
    if(!edgeList.empty()){
    int top = edgeList.at(0).top;
    int bottom = edgeList.at(edgeCount-1).bottom;
    
    GEdge a = edgeList.at(0);
    GEdge b = edgeList.at(1);
    edgeList.pop_front();
    edgeList.pop_front();
    
   for(int y = top; y < bottom; y++){
        if(y > a.bottom && y > b.bottom && !edgeList.empty()){
            a = edgeList.at(0);
            edgeList.pop_front();
            if(!edgeList.empty()){
                b = edgeList.at(0);
                edgeList.pop_front();
            }
        }
        else if(y > a.bottom && y < b.bottom && !edgeList.empty()){
            a = edgeList.at(0);
            edgeList.pop_front();
        }
        else if(y < a.bottom && y > b.bottom && !edgeList.empty()){
            b = edgeList.at(0);
            edgeList.pop_front();
        }
        float yMid = float(y + 0.5);
        int xA = evaluate(a, yMid);
        int xB = evaluate(b, yMid);
        blit(xA, y, xB, paint);
    }
    }
}

void MyCanvas::drawPath(const GPath& path, const GPaint& paint){

    std::shared_ptr<GPath> copy = path.transform(ctm);
    std::deque<GEdge> edgeList;
    float tolerance = 0.25f;
    
    int yMin;
    int yMax;

    GPoint pts[GPath::kMaxNextPoints];
    GPath::Edger edg(*copy);

    while (auto v = edg.next(pts)) {

        switch (v.value()) {
            case GPathVerb::kLine:
            {
                if (GRoundToInt(pts[0].y) != GRoundToInt(pts[1].y)) {
                    GEdge e = makeEdge(pts[0], pts[1]);
                    if(GRoundToInt(pts[0].y) > GRoundToInt(pts[1].y)){
                        e.w = 1;
                    }
                    else{
                        e.w = -1;
                    }
                    edgeList.push_front(e);
                }
                break;
            }

            case GPathVerb::kQuad:
            {
                float e_x = (pts[0].x - (2 * pts[1].x) + pts[2].x)/4;
                float e_y = (pts[0].y - (2 * pts[1].y) + pts[2].y)/4;
                float e_m = sqrt((e_x * e_x) + (e_y * e_y));
                int num_segs = (int)std::ceil(sqrt(e_m/tolerance));
               
                for(int i = 0; i < num_segs; i++){
                    GPoint a = lerp_quad(pts[0], pts[1], pts[2], (float)i/num_segs);
                    GPoint b = lerp_quad(pts[0], pts[1], pts[2], (float)(i+1)/num_segs);

                    if (GRoundToInt(a.y) != GRoundToInt(b.y)) {
                    GEdge e = makeEdge(a, b);
                    if(GRoundToInt(a.y) > GRoundToInt(b.y)){
                        e.w = 1;
                    }
                    else{
                        e.w = -1;
                    }
                    edgeList.push_front(e);
                    }
                }
                break;
            }

            case GPathVerb::kCubic:
            {
                float e0_x = pts[0].x - (2 * pts[1].x) + pts[2].x;
                float e0_y = pts[0].y - (2 * pts[1].y) + pts[2].y;

                float e1_x = pts[1].x - (2 * pts[2].x) + pts[3].x;
                float e1_y = pts[1].y - (2 * pts[2].y) + pts[3].y;

                float e_x = std::max(abs(e0_x), abs(e1_x));
                float e_y = std::max(abs(e0_y), abs(e1_y));
                float e_m = sqrt((e_x * e_x) + (e_y * e_y));
                int num_segs = (int)std::ceil(sqrt((e_m * 3)/(4 * tolerance)));
               
                for(int i = 0; i < num_segs; i++){
                    GPoint a = lerp_cubic(pts[0], pts[1], pts[2], pts[3], (float)i/num_segs);
                    GPoint b = lerp_cubic(pts[0], pts[1], pts[2], pts[3], (float)(i+1)/num_segs);

                    if (GRoundToInt(a.y) != GRoundToInt(b.y)) {
                    GEdge e = makeEdge(a, b);
                    if(GRoundToInt(a.y) > GRoundToInt(b.y)){
                        e.w = 1;
                    }
                    else{
                        e.w = -1;
                    }
                    edgeList.push_front(e);
                    }
                }
                break;
            }
        }
    }

    if(!edgeList.empty()){

    std::sort(edgeList.begin(),edgeList.end(),lessThan);
    yMin = edgeList.at(0).top;
    yMax = edgeList.at(0).bottom;

       for(int i = 0; i < edgeList.size(); i++){
            yMin = std::min(yMin, edgeList.at(i).top);
            yMax = std::max(yMax, edgeList.at(i).bottom);
       }
    
    for(int y = yMin; y < yMax; y++){
        size_t i = 0;
        int w = 0;
        int L;
        while (i < edgeList.size() && isValid(edgeList.at(i), (float)(y + 0.5))) {
            int x = evaluate(edgeList.at(i), (float)(y + 0.5));
            if (w == 0) {
            L = x;
            }
            w += edgeList.at(i).w;  // +1 or -1
            if (w == 0) {
		    int R = x;
            blit(L, y, R, paint);
            }
        if (isValid(edgeList.at(i), (float)(y + 1.5))) {
            // if you track “currX”, bump it by slope here
            i += 1;
        } else {
        	// we’re done with this edge
            edgeList.erase(edgeList.begin()+i);
        }
        }
        assert(w == 0);
    // now i is the number of remaining valid edges

    // account for any new edges that will be valid for next y
        while (i < edgeList.size() && isValid(edgeList.at(i), (float)(y + 1.5))) {
            i += 1;
        }
    // now i also includes the number of edges that will be valid for next y
        auto comp = [=](const GEdge& a, const GEdge& b){
            return evaluate_precise(a, (float)(y + 1.5)) < evaluate_precise(b, (float)(y + 1.5));
        };
        std::sort(edgeList.begin(),edgeList.begin() + i,comp);
        
    }
    }
}

void MyCanvas::drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[], int count, const int indices[], const GPaint& paint){
    int n = 0;
    for(int i = 0; i < count; i++){
        GPoint points[3] = {verts[indices[n+0]], verts[indices[n+1]], verts[indices[n+2]]};

        if(colors != nullptr && texs == nullptr){
            GColor color[3] = {colors[indices[n+0]], colors[indices[n+1]], colors[indices[n+2]]};
            this->drawTriangleWithColor(points, color, paint);
        }

        else if(colors == nullptr && texs != nullptr){
            GPoint tex[3] = {texs[indices[n+0]], texs[indices[n+1]], texs[indices[n+2]]};
            this->drawTriangleWithTex(points, tex, paint);
        }

        else if(colors != nullptr && texs != nullptr){
            GColor color[3] = {colors[indices[n+0]], colors[indices[n+1]], colors[indices[n+2]]};
            GPoint tex[3] = {texs[indices[n+0]], texs[indices[n+1]], texs[indices[n+2]]};
            this->drawTriangleWithBoth(points, color, tex, paint);
        }      
        n += 3;
    }
}

void MyCanvas::drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4], int level, const GPaint& paint){
    int num_columns = level + 1;
    int num_rows = num_columns;
  
    if(colors != nullptr && texs == nullptr){
    for(int i = 0; i < num_rows; i++){
        for(int j = 0; j < num_columns; j++){

                GPoint a = this->lerp_line(this->lerp_line(verts[0], verts[3], (float)i/num_rows), this->lerp_line(verts[1], verts[2], (float)i/num_rows), (float)j/num_columns);
                GPoint b = this->lerp_line(this->lerp_line(verts[0], verts[3], (float)i/num_rows), this->lerp_line(verts[1], verts[2], (float)i/num_rows), (float)(j + 1)/num_columns);
                GPoint c = this->lerp_line(this->lerp_line(verts[0], verts[3], (float)(i + 1)/num_rows), this->lerp_line(verts[1], verts[2], (float)(i + 1)/num_rows), (float)(j + 1)/num_columns);
                GPoint d = this->lerp_line(this->lerp_line(verts[0], verts[3], (float)(i + 1)/num_rows), this->lerp_line(verts[1], verts[2], (float)(i + 1)/num_rows), (float)j/num_columns);

                GColor c_a = this->lerp_color(this->lerp_color(colors[0], colors[3], (float)i/num_rows), this->lerp_color(colors[1], colors[2], (float)i/num_rows), (float)j/num_columns);
                GColor c_b = this->lerp_color(this->lerp_color(colors[0], colors[3], (float)i/num_rows), this->lerp_color(colors[1], colors[2], (float)i/num_rows), (float)(j + 1)/num_columns);
                GColor c_c = this->lerp_color(this->lerp_color(colors[0], colors[3], (float)(i + 1)/num_rows), this->lerp_color(colors[1], colors[2], (float)(i + 1)/num_rows), (float)(j + 1)/num_columns);
                GColor c_d = this->lerp_color(this->lerp_color(colors[0], colors[3], (float)(i + 1)/num_rows), this->lerp_color(colors[1], colors[2], (float)(i + 1)/num_rows), (float)j/num_columns);

                GPoint poly1[3] = {a, b, d};
                GPoint poly2[3] = {c, b, d};
                GColor colors1[3] = {c_a, c_b, c_d};
                GColor colors2[3] = {c_c, c_b, c_d};
                this->drawTriangleWithColor(poly1, colors1, paint);
                this->drawTriangleWithColor(poly2, colors2, paint);               
        }
    }
    }

    else if(colors == nullptr && texs != nullptr){
    for(int i = 0; i < num_rows; i++){
        for(int j = 0; j < num_columns; j++){

                GPoint a = this->lerp_line(this->lerp_line(verts[0], verts[3], (float)i/num_rows), this->lerp_line(verts[1], verts[2], (float)i/num_rows), (float)j/num_columns);
                GPoint b = this->lerp_line(this->lerp_line(verts[0], verts[3], (float)i/num_rows), this->lerp_line(verts[1], verts[2], (float)i/num_rows), (float)(j + 1)/num_columns);
                GPoint c = this->lerp_line(this->lerp_line(verts[0], verts[3], (float)(i + 1)/num_rows), this->lerp_line(verts[1], verts[2], (float)(i + 1)/num_rows), (float)(j + 1)/num_columns);
                GPoint d = this->lerp_line(this->lerp_line(verts[0], verts[3], (float)(i + 1)/num_rows), this->lerp_line(verts[1], verts[2], (float)(i + 1)/num_rows), (float)j/num_columns);


                GPoint tex_a = this->lerp_line(this->lerp_line(texs[0], texs[3], (float)i/num_rows), this->lerp_line(texs[1], texs[2], (float)i/num_rows), (float)j/num_columns);
                GPoint tex_b = this->lerp_line(this->lerp_line(texs[0], texs[3], (float)i/num_rows), this->lerp_line(texs[1], texs[2], (float)i/num_rows), (float)(j + 1)/num_columns);
                GPoint tex_c = this->lerp_line(this->lerp_line(texs[0], texs[3], (float)(i + 1)/num_rows), this->lerp_line(texs[1], texs[2], (float)(i + 1)/num_rows), (float)(j + 1)/num_columns);
                GPoint tex_d = this->lerp_line(this->lerp_line(texs[0], texs[3], (float)(i + 1)/num_rows), this->lerp_line(texs[1], texs[2], (float)(i + 1)/num_rows), (float)j/num_columns);
                

                GPoint poly1[3] = {a, b, d};
                GPoint tex1[3] = {tex_a, tex_b, tex_d};
                GPoint poly2[3] = {c, b, d};
                GPoint tex2[3] = {tex_c, tex_b, tex_d};
    
                this->drawTriangleWithTex(poly1, tex1, paint);
                this->drawTriangleWithTex(poly2, tex2, paint);               
        }
    }
    }

    else if(colors != nullptr && texs != nullptr){
    for(int i = 0; i < num_rows; i++){
        for(int j = 0; j < num_columns; j++){

                GPoint a = this->lerp_line(this->lerp_line(verts[0], verts[3], (float)i/num_rows), this->lerp_line(verts[1], verts[2], (float)i/num_rows), (float)j/num_columns);
                GPoint b = this->lerp_line(this->lerp_line(verts[0], verts[3], (float)i/num_rows), this->lerp_line(verts[1], verts[2], (float)i/num_rows), (float)(j + 1)/num_columns);
                GPoint c = this->lerp_line(this->lerp_line(verts[0], verts[3], (float)(i + 1)/num_rows), this->lerp_line(verts[1], verts[2], (float)(i + 1)/num_rows), (float)(j + 1)/num_columns);
                GPoint d = this->lerp_line(this->lerp_line(verts[0], verts[3], (float)(i + 1)/num_rows), this->lerp_line(verts[1], verts[2], (float)(i + 1)/num_rows), (float)j/num_columns);

                GColor c_a = this->lerp_color(this->lerp_color(colors[0], colors[3], (float)i/num_rows), this->lerp_color(colors[1], colors[2], (float)i/num_rows), (float)j/num_columns);
                GColor c_b = this->lerp_color(this->lerp_color(colors[0], colors[3], (float)i/num_rows), this->lerp_color(colors[1], colors[2], (float)i/num_rows), (float)(j + 1)/num_columns);
                GColor c_c = this->lerp_color(this->lerp_color(colors[0], colors[3], (float)(i + 1)/num_rows), this->lerp_color(colors[1], colors[2], (float)(i + 1)/num_rows), (float)(j + 1)/num_columns);
                GColor c_d = this->lerp_color(this->lerp_color(colors[0], colors[3], (float)(i + 1)/num_rows), this->lerp_color(colors[1], colors[2], (float)(i + 1)/num_rows), (float)j/num_columns);

                GPoint tex_a = this->lerp_line(this->lerp_line(texs[0], texs[3], (float)i/num_rows), this->lerp_line(texs[1], texs[2], (float)i/num_rows), (float)j/num_columns);
                GPoint tex_b = this->lerp_line(this->lerp_line(texs[0], texs[3], (float)i/num_rows), this->lerp_line(texs[1], texs[2], (float)i/num_rows), (float)(j + 1)/num_columns);
                GPoint tex_c = this->lerp_line(this->lerp_line(texs[0], texs[3], (float)(i + 1)/num_rows), this->lerp_line(texs[1], texs[2], (float)(i + 1)/num_rows), (float)(j + 1)/num_columns);
                GPoint tex_d = this->lerp_line(this->lerp_line(texs[0], texs[3], (float)(i + 1)/num_rows), this->lerp_line(texs[1], texs[2], (float)(i + 1)/num_rows), (float)j/num_columns);
                
                GPoint poly1[3] = {a, b, d};
                GPoint poly2[3] = {c, b, d};
                GColor colors1[3] = {c_a, c_b, c_d};
                GColor colors2[3] = {c_c, c_b, c_d};
                GPoint tex1[3] = {tex_a, tex_b, tex_d};
                GPoint tex2[3] = {tex_c, tex_b, tex_d};
                this->drawTriangleWithBoth(poly1, colors1, tex1, paint);
                this->drawTriangleWithBoth(poly2, colors2, tex2, paint);               
        }
    }
    }
}

void MyCanvas::drawTriangleWithColor(const GPoint verts[3], const GColor colors[3], const GPaint& paint){
    GPaint p = paint;
    std::shared_ptr<GShader> sh = GCreateMyTriColorShader(verts[0], verts[1], verts[2], colors[0], colors[1], colors[2]);
    p.setShader(sh);
    GPoint poly[] = {verts[0], verts[1], verts[2]};
    this->drawConvexPolygon(poly, 3, p);
}

void MyCanvas::drawTriangleWithTex(const GPoint verts[3], const GPoint tex[3], const GPaint& paint){
    GPaint p = paint;
    float tu_x = tex[1].x - tex[0].x;
    float tu_y = tex[1].y - tex[0].y;

    float tv_x = tex[2].x - tex[0].x;
    float tv_y = tex[2].y - tex[0].y;

    float pu_x = verts[1].x - verts[0].x;
    float pu_y = verts[1].y - verts[0].y;

    float pv_x = verts[2].x - verts[0].x;
    float pv_y = verts[2].y - verts[0].y;

    GMatrix t = GMatrix(tu_x, tv_x, tex[0].x, tu_y, tv_y, tex[0].y);
    GMatrix pt = GMatrix(pu_x, pv_x, verts[0].x, pu_y, pv_y, verts[0].y);
    GMatrix imx = pt * t.invert().value();
    std::shared_ptr<GShader> sh = GCreateProxyShader(paint.peekShader(), imx);
    p.setShader(sh);
    GPoint poly[] = {verts[0], verts[1], verts[2]};
    this->drawConvexPolygon(poly, 3, p);
}

void MyCanvas::drawTriangleWithBoth(const GPoint verts[3], const GColor colors[3], const GPoint tex[3], const GPaint& paint){
    GPaint p = paint;
    float tu_x = tex[1].x - tex[0].x;
    float tu_y = tex[1].y - tex[0].y;

    float tv_x = tex[2].x - tex[0].x;
    float tv_y = tex[2].y - tex[0].y;

    float pu_x = verts[1].x - verts[0].x;
    float pu_y = verts[1].y - verts[0].y;

    float pv_x = verts[2].x - verts[0].x;
    float pv_y = verts[2].y - verts[0].y;

    GMatrix t = GMatrix(tu_x, tv_x, tex[0].x, tu_y, tv_y, tex[0].y);
    GMatrix pt = GMatrix(pu_x, pv_x, verts[0].x, pu_y, pv_y, verts[0].y);
    GMatrix imx = pt * t.invert().value();

    std::shared_ptr<GShader> triShader = GCreateMyTriColorShader(verts[0], verts[1], verts[2], colors[0], colors[1], colors[2]);
    std::shared_ptr<GShader> texShader = GCreateProxyShader(paint.peekShader(), imx);
    std::shared_ptr<GShader> sh = GCreateComposShader(triShader, texShader);
    p.setShader(sh);
    GPoint poly[] = {verts[0], verts[1], verts[2]};
    this->drawConvexPolygon(poly, 3, p);
}


void MyCanvas::blit(int left, int y, int right, const GPaint& paint){
    if(right < left){
        int t = right;
        right = left;
        left = t;
    }
    if((y < 0) || (y >= fDevice.height()) || (left < 0 && right < 0) || (left >= fDevice.width() && right >= fDevice.width())){
        return;
    }
    if(left < 0 && right >= fDevice.width()){
        left = 0;
        right = fDevice.width() - 1;
    }
    if(left < 0){
        left = 0;
    }
    if(right >= fDevice.width()){
        right = fDevice.width() - 1;
    }
    if(paint.peekShader() != nullptr){
        int width = (right - left + 1);
        GPixel shade_addr[width];
        GPixel *row_addr = nullptr;
        if(paint.peekShader()->setContext(this->ctm)){
            paint.peekShader()->shadeRow(left, y, width, shade_addr);
            row_addr = fDevice.getAddr(left, y);
            for (int x = 0; x < width; x++){
                GPixel dst_pixel = row_addr[x];
                GPixel src_pixel = shade_addr[x];
                row_addr[x] = blend(src_pixel, dst_pixel, paint.getBlendMode());
            }
        }
    }
    else{
        GPixel src_pixel = GPixel_PackARGB(GRoundToInt(paint.getColor().a*255.0f), GRoundToInt(paint.getColor().r*paint.getColor().a*255.0f), GRoundToInt(paint.getColor().g*paint.getColor().a*255.0f), GRoundToInt(paint.getColor().b*paint.getColor().a*255.0f));
        GPixel *row_addr = nullptr;   
        row_addr = fDevice.getAddr(0,y);
        for (int x = left; x < right; x++){
            GPixel dst_pixel = row_addr[x];
            row_addr[x] = blend(src_pixel, dst_pixel, paint.getBlendMode());
        }
    }
    
}

GColor MyCanvas::lerp_color(const GColor& a, const GColor& b, float t){
    GColor color;
    color.a = a.a + (b.a - a.a) * t;
    color.r = a.r + (b.r - a.r) * t;
    color.g = a.g + (b.g - a.g) * t;
    color.b = a.b + (b.b - a.b) * t;
    return color;
}

GPoint MyCanvas::lerp_line(const GPoint& a, const GPoint& b, float t){
    float s = 1 - t;
    float q_x = (s * a.x) + (t * b.x);
    float q_y = (s * a.y) + (t * b.y);
    GPoint g; g.x = q_x; g.y = q_y;
    return g;
}

GPoint MyCanvas::lerp_quad(const GPoint& a, const GPoint& b, const GPoint& c, float t){
    float s = 1 - t;
    float q_x = (a.x * s * s) + (2 * b.x * t * s) + (c.x * t * t);
    float q_y = (a.y * s * s) + (2 * b.y * t * s) + (c.y * t * t);
    GPoint g; g.x = q_x; g.y = q_y;
    return g;
}

GPoint MyCanvas::lerp_cubic(const GPoint& a, const GPoint& b, const GPoint& c, const GPoint& d, float t){
    float s = 1 - t;
    float q_x = (a.x * s * s * s) + (3 * b.x * t * s * s) + (3 * c.x * s * t * t) + (d.x * t * t * t);
    float q_y = (a.y * s * s * s) + (3 * b.y * t * s * s) + (3 * c.y * s * t * t) + (d.y * t * t * t);
    GPoint g; g.x = q_x; g.y = q_y;
    return g;
}

int MyCanvas::evaluate(GEdge a, float y){
    return GRoundToInt((y*a.m) + a.b);
}

float MyCanvas::evaluate_precise(GEdge a, float y){
    return (y*a.m) + a.b;
}

GPixel MyCanvas::blend(GPixel& src, GPixel& dst, GBlendMode mode){
    switch(mode){
        case GBlendMode::kClear: return clear_mode(src, dst);

        case GBlendMode::kSrc: return src_mode(src, dst);

        case GBlendMode::kDst: return dst_mode(src, dst); 

        case GBlendMode::kSrcOver: return srcOver_mode(src, dst);

        case GBlendMode::kDstOver: return dstOver_mode(src, dst); 

        case GBlendMode::kSrcIn: return srcIn_mode(src, dst); 

        case GBlendMode::kDstIn: return dstIn_mode(src, dst); 

        case GBlendMode::kSrcOut: return srcOut_mode(src, dst);

        case GBlendMode::kDstOut: return dstOut_mode(src, dst); 

        case GBlendMode::kSrcATop: return srcAtop_mode(src, dst); 

        case GBlendMode::kDstATop: return dstAtop_mode(src, dst); 
        
        case GBlendMode::kXor: return xor_mode(src, dst);

        default: return GPixel_PackARGB(0,0,0,0); 
    }
    
}

GPixel MyCanvas::clear_mode(GPixel& src, GPixel& dst){
    return GPixel_PackARGB(0,0,0,0);
}

GPixel MyCanvas::src_mode(GPixel& src, GPixel& dst){
    return src;
}

GPixel MyCanvas::dst_mode(GPixel& src, GPixel& dst){
    return dst;
}

GPixel MyCanvas::srcOver_mode(GPixel& src, GPixel& dst){
    float srcA = GPixel_GetA(src)/255.0f;
    float dstA = GPixel_GetA(dst)/255.0f;
    int resultA = (int)(255.0f * (srcA + dstA * (1.0f - srcA)));
    int resultR = (int)(GPixel_GetR(src) + GPixel_GetR(dst) * (1.0f - srcA));
    int resultG = (int)(GPixel_GetG(src) + GPixel_GetG(dst) * (1.0f - srcA));
    int resultB = (int)(GPixel_GetB(src) + GPixel_GetB(dst) * (1.0f - srcA));
    return GPixel_PackARGB(resultA,resultR,resultG,resultB);
}

GPixel MyCanvas::dstOver_mode(GPixel& src, GPixel& dst){
    float srcA = GPixel_GetA(src)/255.0f;
    float dstA = GPixel_GetA(dst)/255.0f;
    int resultA = (int)(255.0f * (dstA + srcA * (1.0f - dstA)));
    int resultR = (int)(GPixel_GetR(dst) + GPixel_GetR(src) * (1.0f - dstA));
    int resultG = (int)(GPixel_GetG(dst) + GPixel_GetG(src) * (1.0f - dstA));
    int resultB = (int)(GPixel_GetB(dst) + GPixel_GetB(src) * (1.0f - dstA));
    return GPixel_PackARGB(resultA,resultR,resultG,resultB);
}

GPixel MyCanvas::srcIn_mode(GPixel& src, GPixel& dst){
    float srcA = GPixel_GetA(src)/255.0f;
    float dstA = GPixel_GetA(dst)/255.0f;
    int resultA = (int)(255.0f * (dstA * srcA));
    int resultR = (int)(dstA * GPixel_GetR(src));
    int resultG = (int)(dstA * GPixel_GetG(src));
    int resultB = (int)(dstA * GPixel_GetB(src));
    return GPixel_PackARGB(resultA,resultR,resultG,resultB);
}

GPixel MyCanvas::dstIn_mode(GPixel& src, GPixel& dst){
    float srcA = GPixel_GetA(src)/255.0f;
    float dstA = GPixel_GetA(dst)/255.0f;
    int resultA = (int)(255.0f * (srcA * dstA));
    int resultR = (int)(srcA * GPixel_GetR(dst));
    int resultG = (int)(srcA * GPixel_GetG(dst));
    int resultB = (int)(srcA * GPixel_GetB(dst));
    return GPixel_PackARGB(resultA,resultR,resultG,resultB);
}

GPixel MyCanvas::srcOut_mode(GPixel& src, GPixel& dst){
    float srcA = GPixel_GetA(src)/255.0f;
    float dstA = GPixel_GetA(dst)/255.0f;
    int resultA = (int)(255.0f * ((1.0f - dstA) * srcA));
    int resultR = (int)((1.0f - dstA) * GPixel_GetR(src));
    int resultG = (int)((1.0f - dstA) * GPixel_GetG(src));
    int resultB = (int)((1.0f - dstA) * GPixel_GetB(src));
    return GPixel_PackARGB(resultA,resultR,resultG,resultB);
}

GPixel MyCanvas::dstOut_mode(GPixel& src, GPixel& dst){
    float srcA = GPixel_GetA(src)/255.0f;
    float dstA = GPixel_GetA(dst)/255.0f;
    int resultA = (int)(255.0f * ((1.0f - srcA) * dstA));
    int resultR = (int)((1.0f - srcA) * GPixel_GetR(dst));
    int resultG = (int)((1.0f - srcA) * GPixel_GetG(dst));
    int resultB = (int)((1.0f - srcA) * GPixel_GetB(dst));
    return GPixel_PackARGB(resultA,resultR,resultG,resultB);
}

GPixel MyCanvas::srcAtop_mode(GPixel& src, GPixel& dst){
    float srcA = GPixel_GetA(src)/255.0f;
    float dstA = GPixel_GetA(dst)/255.0f;
    int resultA = (int)(255.0f * (dstA * srcA + (1.0f - srcA) * dstA));
    int resultR = (int)(dstA * GPixel_GetR(src) + (1.0f - srcA) * GPixel_GetR(dst));
    int resultG = (int)(dstA * GPixel_GetG(src) + (1.0f - srcA) * GPixel_GetG(dst));
    int resultB = (int)(dstA * GPixel_GetB(src) + (1.0f - srcA) * GPixel_GetB(dst));
    return GPixel_PackARGB(resultA,resultR,resultG,resultB);
}

GPixel MyCanvas::dstAtop_mode(GPixel& src, GPixel& dst){
    float srcA = GPixel_GetA(src)/255.0f;
    float dstA = GPixel_GetA(dst)/255.0f;
    int resultA = (int)(255.0f * (srcA * dstA+ (1 - dstA) * srcA));
    int resultR = (int)(srcA * GPixel_GetR(dst) + (1.0f - dstA) * GPixel_GetR(src));
    int resultG = (int)(srcA * GPixel_GetG(dst) + (1.0f - dstA) * GPixel_GetG(src));
    int resultB = (int)(srcA * GPixel_GetB(dst) + (1.0f - dstA) * GPixel_GetB(src));
    return GPixel_PackARGB(resultA,resultR,resultG,resultB);
}

GPixel MyCanvas::xor_mode(GPixel& src, GPixel& dst){
    float srcA = GPixel_GetA(src)/255.0f;
    float dstA = GPixel_GetA(dst)/255.0f;
    int resultA = (int)(255.0f * ((1.0f- srcA) * dstA + (1.0f - dstA) * srcA));
    int resultR = (int)((1.0f - srcA) * GPixel_GetR(dst) + (1.0f - dstA) * GPixel_GetR(src));
    int resultG = (int)((1.0f - srcA) * GPixel_GetG(dst) + (1.0f - dstA) * GPixel_GetG(src));
    int resultB = (int)((1.0f - srcA) * GPixel_GetB(dst) + (1.0f - dstA) * GPixel_GetB(src));
    return GPixel_PackARGB(resultA,resultR,resultG,resultB);
}

void MyCanvas::translate(float x, float y) {
        this->concat(GMatrix::Translate(x, y));
    }

void MyCanvas::scale(float x, float y) {
        this->concat(GMatrix::Scale(x, y));
    }

void MyCanvas::rotate(float radians) {
        this->concat(GMatrix::Rotate(radians));
    }

std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& device) {
    return std::unique_ptr<GCanvas>(new MyCanvas(device));
}

std::unique_ptr<GFinal> GCreateFinal(){
    return std::make_unique<MyFinal>();
}

std::string GDrawSomething(GCanvas* canvas, GISize dim) {
    // as fancy as you like
    // ...
    // canvas->clear(...);
    // canvas->fillRect(...);
    GColor blue; blue.r = 0; blue.b = 1; blue.g = 0; blue.a = 1;
    GColor yellow; yellow.r = 1; yellow.b = 0; yellow.g = 1; yellow.a = 1;
    GColor green; green.r = 0; green.b = 0; green.g = 1; green.a = 1;
    GPathBuilder bu;
    bu.moveTo(220.0f, 120.0f);
    bu.lineTo(39.0983f, 178.779f);
    bu.lineTo(150.902f, 24.8943f);
    bu.lineTo(150.902f, 215.106f);
    bu.lineTo(39.0983f, 61.2215f);
    //bu.reset();
    //GRect r; r.left = 100; r.right = 200; r.top = 50; r.bottom = 300;
    //bu.addRect(r, GPathDirection::kCW);
    std::shared_ptr<GPath> path = bu.detach();
    canvas->clear(blue);
    canvas->drawPath(*path, GPaint(green));
    return "green star";
}
/*
 *  Copyright 2024 Christopher Vasquez
 */

#include "include/GShader.h"
#include <optional>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>
#include "include/GMatrix.h"
#include "include/GBitmap.h"
class GMatrix;
class GBitmap;

class bmShader : public GShader{
        public:
        bmShader(const GBitmap& bm, const GMatrix& localMatrix, GTileMode mode) : bm(bm), mx(localMatrix), ctm(GMatrix()), mode(mode){}

        bool isOpaque(){
            return bm.isOpaque();
        }
        bool setContext(const GMatrix& ctm){
            if(!ctm.invert().has_value()){
                return false;
            }
            else{
                this->ctm = ctm;
                imx = mx.invert().value() * ctm.invert().value();
                return true;
            }
        }
        void shadeRow(int x, int y, int count, GPixel row[]){
        switch(mode){
            case GTileMode::kClamp:
            {
                for(int i = 0; i < count; i++){
                float px = (float)(imx.value()[0] * (x + i + 0.5) + imx.value()[2] * (y + 0.5) + imx.value()[4]);
                float py = (float)(imx.value()[1] * (x + i + 0.5) + imx.value()[3] * (y + 0.5) + imx.value()[5]);
                int ix = (int)std::floor(px);                
                int iy = (int)std::floor(py);
                if(ix < 0){
                    ix = 0;
                }
                if(ix >= bm.width()){
                    ix = bm.width()-1;
                }
                if(iy < 0){
                    iy = 0;
                }
                if(iy >= bm.height()){
                    iy = bm.height()-1;
                }
                row[i] = *bm.getAddr(ix, iy);
            }
            break;
            }

            case GTileMode::kRepeat:{
                for(int i = 0; i < count; i++){
                float px = (float)(imx.value()[0] * (x + i + 0.5) + imx.value()[2] * (y + 0.5) + imx.value()[4]);
                float py = (float)(imx.value()[1] * (x + i + 0.5) + imx.value()[3] * (y + 0.5) + imx.value()[5]);
                int ix = (int)std::floor(px);                
                int iy = (int)std::floor(py);
                ix = abs(ix % bm.width());
                iy = abs(iy % bm.height());
                row[i] = *bm.getAddr(ix, iy);
            }
            break;
            }

            case GTileMode::kMirror:{
                
                for(int i = 0; i < count; i++){
                float px = (float)(imx.value()[0] * (x + i + 0.5) + imx.value()[2] * (y + 0.5) + imx.value()[4]);
                float py = (float)(imx.value()[1] * (x + i + 0.5) + imx.value()[3] * (y + 0.5) + imx.value()[5]);
        
                if(((abs)((int)std::floor(px/bm.width())) % 2) == 0 && ((abs)((int)std::floor(py/bm.height())) % 2) == 0){
                    int ix = (abs)((int)std::floor(px) % bm.width());
                    int iy = (abs)((int)std::floor(py) % bm.height());
                    row[i] = *bm.getAddr(ix, iy);
                }

                else if(((abs)((int)std::floor(px/bm.width())) % 2) == 1 && ((abs)((int)std::floor(py/bm.height())) % 2) == 0){
                    int ix = (bm.width() - 1) - ((abs)((int)std::floor(px) % bm.width()));
                    int iy = (abs)((int)std::floor(py) % bm.height());
                    row[i] = *bm.getAddr(ix, iy);
                }

                else if(((abs)((int)std::floor(px/bm.width())) % 2) == 0 && ((abs)((int)std::floor(py/bm.height())) % 2) == 1){
                    int ix = (abs)((int)std::floor(px) % bm.width());
                    int iy = (bm.height() - 1) - ((abs)((int)std::floor(py) % bm.height()));
                    row[i] = *bm.getAddr(ix, iy);
                }

                else if(((abs)((int)std::floor(px/bm.width())) % 2) == 1 && ((abs)((int)std::floor(py/bm.height())) % 2) == 1){
                    int ix = (bm.width() - 1) - ((abs)((int)std::floor(px) % bm.width()));
                    int iy = (bm.height() - 1) - ((abs)((int)std::floor(py) % bm.height()));
                    row[i] = *bm.getAddr(ix, iy);
                }

            }
            break;
            }
        
        }
        }
        private:
        GBitmap bm;
        GMatrix mx;
        GMatrix ctm;
        nonstd::optional<GMatrix> imx;
        GTileMode mode;      
    };

    class lgShader : public GShader{
        public:
        lgShader(GPoint p0, GPoint p1, const GColor colors[], int count, GTileMode mode) : p0(p0), p1(p1), count(count), ctm(GMatrix()), mode(mode){
            float dx = p1.x - p0.x;
            float dy = p1.y - p0.y;
            lgm = GMatrix(dx, -dy, p0.x, dy, dx, p0.y);
            scmx = GMatrix::Scale((float)(count-1), (float)(count-1));
            k = count - 1;
            c = new GColor[count];
            std::copy(colors, colors + count, c);
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
                imx = scmx * lgm.invert().value() * ctm.invert().value();
                return true;
            }
        }
        void shadeRow(int x, int y, int count, GPixel row[]){
            switch(mode){
                case GTileMode::kClamp:
                {
                GColor color;
                for(int i = 0; i < count; i++){
                float px = (float)(imx.value()[0] * (x + i + 0.5) + imx.value()[2] * (y + 0.5) + imx.value()[4]);
                if(px < 0.0f){
                    px = 0.0f;
                }
                if(px > (float)k){
                    px = (float)k;
                }

                int ix = (int)(std::floor(px));                
                
                float t = (float)(px - ix);

                if((px == 0.0f || this->count == 1)){
                    row[i] = GPixel_PackARGB(GRoundToInt(c[0].a*255.0f),GRoundToInt(c[0].r*(c[0].a)*255.0f),GRoundToInt(c[0].g*(c[0].a)*255.0f),GRoundToInt(c[0].b*(c[0].a)*255.0f));
                }
                else if(px == (float)k){
                    row[i] = GPixel_PackARGB(GRoundToInt(c[k].a*255.0f),GRoundToInt(c[k].r*(c[k].a)*255.0f),GRoundToInt(c[k].g*(c[k].a)*255.0f),GRoundToInt(c[k].b*(c[k].a)*255.0f));
                }
                else{
                    color.a = c[ix].a + (c[ix + 1].a - c[ix].a) * t;
                    color.r = c[ix].r + (c[ix + 1].r - c[ix].r) * t;
                    color.g = c[ix].g + (c[ix + 1].g - c[ix].g) * t;
                    color.b = c[ix].b + (c[ix + 1].b - c[ix].b) * t;
                    row[i] = GPixel_PackARGB(GRoundToInt(color.a*255.0f),GRoundToInt(color.r*color.a*255.0f),GRoundToInt(color.g*color.a*255.0f),GRoundToInt(color.b*color.a*255.0f));
                }
                }
                break;
                }
                

                case GTileMode::kRepeat:{
                GColor color;
                for(int i = 0; i < count; i++){
                float px = (float)(imx.value()[0] * (x + i + 0.5) + imx.value()[2] * (y + 0.5) + imx.value()[4]);

                int ix = (abs)((int)std::floor(px) % k);                
                
                float t = (abs)(px - std::floor(px));

                color.a = c[ix].a + (c[ix + 1].a - c[ix].a) * t;
                color.r = c[ix].r + (c[ix + 1].r - c[ix].r) * t;
                color.g = c[ix].g + (c[ix + 1].g - c[ix].g) * t;
                color.b = c[ix].b + (c[ix + 1].b - c[ix].b) * t;
                row[i] = GPixel_PackARGB(GRoundToInt(color.a*255.0f),GRoundToInt(color.r*color.a*255.0f),GRoundToInt(color.g*color.a*255.0f),GRoundToInt(color.b*color.a*255.0f));
                
                }
                break;
                }
            

                case GTileMode::kMirror:{
                GColor color;
                for(int i = 0; i < count; i++){
                float px = (float)(imx.value()[0] * (x + i + 0.5) + imx.value()[2] * (y + 0.5) + imx.value()[4]);

                if((((abs)((int)std::floor(px/k)) % 2)) == 0){
                    int ix = (abs)((int)std::floor(px) % k);                
                
                    float t = (abs)(px - std::floor(px));

        
                    color.a = c[ix].a + (c[ix + 1].a - c[ix].a) * t;
                    color.r = c[ix].r + (c[ix + 1].r - c[ix].r) * t;
                    color.g = c[ix].g + (c[ix + 1].g - c[ix].g) * t;
                    color.b = c[ix].b + (c[ix + 1].b - c[ix].b) * t;
                    row[i] = GPixel_PackARGB(GRoundToInt(color.a*255.0f),GRoundToInt(color.r*color.a*255.0f),GRoundToInt(color.g*color.a*255.0f),GRoundToInt(color.b*color.a*255.0f));
                }
                else if(((abs)((int)(std::floor(px/k)) % 2)) == 1){
                    int ix = (k - 1) - (abs)((int)std::floor(px) % k);                
                
                    float t = 1 - (abs)(px - std::floor(px));

                    color.a = c[ix].a + (c[ix + 1].a - c[ix].a) * t;
                    color.r = c[ix].r + (c[ix + 1].r - c[ix].r) * t;
                    color.g = c[ix].g + (c[ix + 1].g - c[ix].g) * t;
                    color.b = c[ix].b + (c[ix + 1].b - c[ix].b) * t;
                    row[i] = GPixel_PackARGB(GRoundToInt(color.a*255.0f),GRoundToInt(color.r*color.a*255.0f),GRoundToInt(color.g*color.a*255.0f),GRoundToInt(color.b*color.a*255.0f));
                }
                }
                break;
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
        GTileMode mode;
    };


std::shared_ptr<GShader> GCreateBitmapShader(const GBitmap& bm, const GMatrix& localMatrix, GTileMode mode){
    if(bm.pixels() == NULL || !localMatrix.invert().has_value()){
        return nullptr;
    }
    else{
        return std::make_shared<bmShader>(bm, localMatrix, mode);
    }
}

std::shared_ptr<GShader> GCreateLinearGradient(GPoint p0, GPoint p1, const GColor colors[], int count, GTileMode mode){
    if(count < 1){
        return nullptr;
    }
    else{
        return std::make_shared<lgShader>(p0, p1, colors, count, mode);
    }
}
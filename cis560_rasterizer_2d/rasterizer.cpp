#include "rasterizer.h"
#include <iostream>
#include <vec4.h>

float min(float a, float b, float c){
    if(a<=b && a<=c) return a;
    if(b<=a && b<=c) return b;
    return c;
}
float max(float a, float b, float c){
    if(a>=b && a>=c) return a;
    if(b>=a && b>=c) return b;
    return c;
}

void BoundingBox::Computing(vec4 v1, vec4 v2, vec4 v3){
    Xl = min(v1[0], v2[0], v3[0]);
    Xr = max(v1[0], v2[0], v3[0]);
    Yl = min(v1[1], v2[1], v3[1]);
    Yr = max(v1[1], v2[1], v3[1]);
}


Rasterizer::Rasterizer(const std::vector<Polygon> &polygons)
    : polygons(polygons)
{}

vec4 Barycentric_Interpolation(vec4 p, Vertex S1, Vertex S2, Vertex S3){
    vec4 p1(S1.pos[0], S1.pos[1], 0, 0), p2(S2.pos[0], S2.pos[1], 0, 0), p3(S3.pos[0], S3.pos[1], 0, 0);
    vec4 c, c1 = S1.color, c2 = S2.color, c3 = S3.color;
    float A, A1, A2, A3;
    A = 0.5 * length(cross((p1 - p2),(p3 - p2)));
    A1 = 0.5 * length(cross((p3 - p),(p2 - p)));
    A2 = 0.5 * length(cross((p1 - p),(p3 - p)));
    A3 = 0.5 * length(cross((p2 - p),(p1 - p)));
    c = c1 * A1 / A + c2 * A2 / A + c3 * A3 / A;
    return c;
}

static inline Comparator(Polygon A, Polygon B){
    return A.Depth() < B.Depth();
}

QImage Rasterizer::RenderScene()
{
    //TODO: Complete the various components of code that make up this function
    //It should return the rasterized image of the current scene

    //Notice that even though it should return a QImage this function compiles.
    //Remember, C++'s compiler is dumb (though it will at least warn you that
    //the function doesn't return anything).
    //BEWARE! If you use a function that is missing its return statement,
    //it will return garbage memory!
    std::cout<<polygons.size()<<"\n";
    QImage img(1024,1024,QImage::Format_RGB32);
    BoundingBox B_Box;
    std::sort(polygons.begin(), polygons.end(), Comparator);
    for(unsigned int i = 0; i < polygons.size(); i++){
        for(unsigned int j =0; j < polygons[i].Triangle_Size(); j++){
            Triangle Tr = polygons[i].TriAt(j);
            Vertex V1, V2, V3;
            V1 = polygons[i].VertAt(Tr.indices[0]);
            V2 = polygons[i].VertAt(Tr.indices[1]);
            V3 = polygons[i].VertAt(Tr.indices[2]);
            B_Box.Computing(V1.pos, V2.pos, V3.pos);
            Line_Segment L1(V1,V2),L2(V2,V3),L3(V3,V1);
            //for each row, compute the Vl, Vr, the color should be paint into Vl -> Vr
            for(int Yrow = B_Box.Yl; Yrow <= B_Box.Yr; Yrow++){
                vec4 p0(0,0,0,0), p1(B_Box.Xr+1,B_Box.Yr+1,0,0);
                Vertex Va(p1,p0,p0,p0), Vb(p1,p0,p0,p0),Vl,Vr;
                bool flagL = false, flagR = false;
                L1.IsIntersection(Yrow, Va, Vb);
                if(Vb.pos[0]<=B_Box.Xr){ // the L1 is parrell with Line Yrow
                    Vl = Va;
                    Vr = Vb;
                    flagL = flagR = true;
                    //Paint By Row
                    for(int k = Vl.pos[0]; k <= Vr.pos[0]; k ++){
                        vec4 p(k, Yrow, 0, 1);
                        vec4 c = Barycentric_Interpolation(p, V1, V2, V3);
                        QRgb value = qRgb(c[0], c[1], c[2]);
                        img.setPixel(k, Yrow, value);
                    }
                    continue;
                }
                if(Va.pos[0]>=B_Box.Xl && Va.pos[0] <=B_Box.Xr){
                    flagL = flagR = true;
                    Vl = Vr = Va;
                }
                L2.IsIntersection(Yrow, Va, Vb);
                if(Vb.pos[0]<=B_Box.Xr){ // the L2 is parrell with Line Yrow
                    Vl = Va;
                    Vr = Vb;
                    flagL = flagR = true;
                    for(int k = Vl.pos[0]; k <= Vr.pos[0]; k ++){
                        vec4 p(k, Yrow, 0, 1);
                        vec4 c = Barycentric_Interpolation(p, V1, V2, V3);
                        QRgb value = qRgb(c[0], c[1], c[2]);
                        img.setPixel(k, Yrow, value);
                    }
                    continue;
                }
                if(Va.pos[0]>=B_Box.Xl && Va.pos[0]<=B_Box.Xr){
                    if(flagL && flagR){
                        if(Va.pos[0]<=Vl.pos[0]) Vl = Va;
                        if(Va.pos[0]>=Vr.pos[0]) Vr = Va;
                    }
                    else{
                        flagL = flagR = true;
                        Vl = Vr = Va;
                    }
                }
                L3.IsIntersection(Yrow, Va, Vb);
                if(Vb.pos[0]<=B_Box.Xr){ // the L3 is parrell with Line Yrow
                    Vl = Va;
                    Vr = Vb;
                    flagL = flagR = true;
                    for(int k = Vl.pos[0]; k <= Vr.pos[0]; k ++){
                        vec4 p(k, Yrow, 0, 1);
                        vec4 c = Barycentric_Interpolation(p, V1, V2, V3);
                        QRgb value = qRgb(c[0], c[1], c[2]);
                        img.setPixel(k, Yrow, value);
                    }
                    continue;
                }
                if(Va.pos[0]>=B_Box.Xl && Va.pos[0]<=B_Box.Xr){
                    if(flagL && flagR){
                        if(Va.pos[0]<=Vl.pos[0]) Vl = Va;
                        if(Va.pos[0]>=Vr.pos[0]) Vr = Va;
                    }
                    else{
                        flagL = flagR = true;
                        Vl = Vr = Va;
                    }
                }
                for(int k = Vl.pos[0]; k <= Vr.pos[0]; k ++){
                    vec4 p(k, Yrow, 0, 1);
                    vec4 c = Barycentric_Interpolation(p, V1, V2, V3);
                    QRgb value = qRgb(c[0], c[1], c[2]);
                    img.setPixel(k, Yrow, value);
                }
            }
        }
    }
    return img;
}

QImage Rasterizer::RenderLine()
{
    QImage img(1024,1024,QImage::Format_RGB32);
    BoundingBox B_Box;
    std::sort(polygons.begin(), polygons.end(), Comparator);
    for(unsigned int i = 0; i < polygons.size(); i++){
        for(unsigned int j =0; j < polygons[i].Triangle_Size(); j++){
            Triangle Tr = polygons[i].TriAt(j);
            Vertex V1, V2, V3;
            V1 = polygons[i].VertAt(Tr.indices[0]);
            V2 = polygons[i].VertAt(Tr.indices[1]);
            V3 = polygons[i].VertAt(Tr.indices[2]);
            B_Box.Computing(V1.pos, V2.pos, V3.pos);
            Line_Segment L1(V1,V2),L2(V2,V3),L3(V3,V1);
            //for each row, compute the Vl, Vr, the color should be paint into Vl -> Vr
            for(int Yrow = B_Box.Yl; Yrow <= B_Box.Yr; Yrow++){
                vec4 p0(0,0,0,0), p1(B_Box.Xr+1,B_Box.Yr+1,0,0);
                Vertex Va(p1,p0,p0,p0), Vb(p1,p0,p0,p0),Vl,Vr;
                bool flagL = false, flagR = false;
                L1.IsIntersection(Yrow, Va, Vb);
                if(Vb.pos[0]<=B_Box.Xr){ // the L1 is parrell with Line Yrow
                    Vl = Va;
                    Vr = Vb;
                    flagL = flagR = true;
                    //Paint By Row
                    for(int k = Vl.pos[0]; k <= Vr.pos[0]; k ++){
                        vec4 p(k, Yrow, 0, 1);
                        vec4 c = Barycentric_Interpolation(p, V1, V2, V3);
                        QRgb value = qRgb(c[0], c[1], c[2]);
                        img.setPixel(k, Yrow, value);
                    }
                    continue;
                }
                if(Va.pos[0]>=B_Box.Xl && Va.pos[0] <=B_Box.Xr){
                    flagL = flagR = true;
                    Vl = Vr = Va;
                }
                L2.IsIntersection(Yrow, Va, Vb);
                if(Vb.pos[0]<=B_Box.Xr){ // the L2 is parrell with Line Yrow
                    Vl = Va;
                    Vr = Vb;
                    flagL = flagR = true;
                    for(int k = Vl.pos[0]; k <= Vr.pos[0]; k ++){
                        vec4 p(k, Yrow, 0, 1);
                        vec4 c = Barycentric_Interpolation(p, V1, V2, V3);
                        QRgb value = qRgb(c[0], c[1], c[2]);
                        img.setPixel(k, Yrow, value);
                    }
                    continue;
                }
                if(Va.pos[0]>=B_Box.Xl && Va.pos[0]<=B_Box.Xr){
                    if(flagL && flagR){
                        if(Va.pos[0]<=Vl.pos[0]) Vl = Va;
                        if(Va.pos[0]>=Vr.pos[0]) Vr = Va;
                    }
                    else{
                        flagL = flagR = true;
                        Vl = Vr = Va;
                    }
                }
                L3.IsIntersection(Yrow, Va, Vb);
                if(Vb.pos[0]<=B_Box.Xr){ // the L3 is parrell with Line Yrow
                    Vl = Va;
                    Vr = Vb;
                    flagL = flagR = true;
                    for(int k = Vl.pos[0]; k <= Vr.pos[0]; k ++){
                        vec4 p(k, Yrow, 0, 1);
                        vec4 c = Barycentric_Interpolation(p, V1, V2, V3);
                        QRgb value = qRgb(c[0], c[1], c[2]);
                        img.setPixel(k, Yrow, value);
                    }
                    continue;
                }
                if(Va.pos[0]>=B_Box.Xl && Va.pos[0]<=B_Box.Xr){
                    if(flagL && flagR){
                        if(Va.pos[0]<=Vl.pos[0]) Vl = Va;
                        if(Va.pos[0]>=Vr.pos[0]) Vr = Va;
                    }
                    else{
                        flagL = flagR = true;
                        Vl = Vr = Va;
                    }
                }
                vec4 p(Vl.pos[0], Yrow, 0, 1);
                vec4 c = Barycentric_Interpolation(p, V1, V2, V3);
                QRgb value = qRgb(c[0], c[1], c[2]);
                img.setPixel(Vl.pos[0], Yrow, value);
                p[0] = Vr.pos[0];
                c = Barycentric_Interpolation(p, V1, V2, V3);
                value = qRgb(c[0], c[1], c[2]);
                img.setPixel(Vr.pos[0], Yrow, value);

            }
        }
    }
    return img;
}

QImage Rasterizer::AntiAliasing(float AA)
{
    //TODO: Complete the various components of code that make up this function
    //It should return the rasterized image of the current scene

    //Notice that even though it should return a QImage this function compiles.
    //Remember, C++'s compiler is dumb (though it will at least warn you that
    //the function doesn't return anything).
    //BEWARE! If you use a function that is missing its return statement,
    //it will return garbage memory!
    QImage img(512,512,QImage::Format_RGB32);
    QImage Scaleimg(512*AA,512*AA,QImage::Format_RGB32);
    BoundingBox B_Box;
    std::sort(polygons.begin(), polygons.end(), Comparator);
    for(unsigned int i = 0; i < polygons.size(); i++){
        for(unsigned int j =0; j < polygons[i].Triangle_Size(); j++){
            Triangle Tr = polygons[i].TriAt(j);
            Vertex V1, V2, V3;
            V1 = polygons[i].VertAt(Tr.indices[0]);
            V2 = polygons[i].VertAt(Tr.indices[1]);
            V3 = polygons[i].VertAt(Tr.indices[2]);
            V1.pos = V1.pos * AA;
            V2.pos = V2.pos * AA;
            V3.pos = V3.pos * AA;
            B_Box.Computing(V1.pos, V2.pos, V3.pos);
            Line_Segment L1(V1,V2),L2(V2,V3),L3(V3,V1);
            //for each row, compute the Vl, Vr, the color should be paint into Vl -> Vr
            for(int Yrow = B_Box.Yl; Yrow <= B_Box.Yr; Yrow++){
                vec4 p0(0,0,0,0), p1(B_Box.Xr+1,B_Box.Yr+1,0,0);
                Vertex Va(p1,p0,p0,p0), Vb(p1,p0,p0,p0),Vl,Vr;
                bool flagL = false, flagR = false;
                L1.IsIntersection(Yrow, Va, Vb);
                if(Vb.pos[0]<=B_Box.Xr){ // the L1 is parrell with Line Yrow
                    Vl = Va;
                    Vr = Vb;
                    flagL = flagR = true;
                    //Paint By Row
                    for(int k = Vl.pos[0]; k <= Vr.pos[0]; k ++){
                        vec4 p(k, Yrow, 0, 1);
                        vec4 c = Barycentric_Interpolation(p, V1, V2, V3);
                        QRgb value = qRgb(c[0], c[1], c[2]);
                        Scaleimg.setPixel(k, Yrow, value);
                    }
                    continue;
                }
                if(Va.pos[0]>=B_Box.Xl && Va.pos[0] <=B_Box.Xr){
                    flagL = flagR = true;
                    Vl = Vr = Va;
                }
                L2.IsIntersection(Yrow, Va, Vb);
                if(Vb.pos[0]<=B_Box.Xr){ // the L2 is parrell with Line Yrow
                    Vl = Va;
                    Vr = Vb;
                    flagL = flagR = true;
                    for(int k = Vl.pos[0]; k <= Vr.pos[0]; k ++){
                        vec4 p(k, Yrow, 0, 1);
                        vec4 c = Barycentric_Interpolation(p, V1, V2, V3);
                        QRgb value = qRgb(c[0], c[1], c[2]);
                        Scaleimg.setPixel(k, Yrow, value);
                    }
                    continue;
                }
                if(Va.pos[0]>=B_Box.Xl && Va.pos[0]<=B_Box.Xr){
                    if(flagL && flagR){
                        if(Va.pos[0]<=Vl.pos[0]) Vl = Va;
                        if(Va.pos[0]>=Vr.pos[0]) Vr = Va;
                    }
                    else{
                        flagL = flagR = true;
                        Vl = Vr = Va;
                    }
                }
                L3.IsIntersection(Yrow, Va, Vb);
                if(Vb.pos[0]<=B_Box.Xr){ // the L3 is parrell with Line Yrow
                    Vl = Va;
                    Vr = Vb;
                    flagL = flagR = true;
                    for(int k = Vl.pos[0]; k <= Vr.pos[0]; k ++){
                        vec4 p(k, Yrow, 0, 1);
                        vec4 c = Barycentric_Interpolation(p, V1, V2, V3);
                        QRgb value = qRgb(c[0], c[1], c[2]);
                        Scaleimg.setPixel(k, Yrow, value);
                    }
                    continue;
                }
                if(Va.pos[0]>=B_Box.Xl && Va.pos[0]<=B_Box.Xr){
                    if(flagL && flagR){
                        if(Va.pos[0]<=Vl.pos[0]) Vl = Va;
                        if(Va.pos[0]>=Vr.pos[0]) Vr = Va;
                    }
                    else{
                        flagL = flagR = true;
                        Vl = Vr = Va;
                    }
                }
                for(int k = Vl.pos[0]; k <= Vr.pos[0]; k ++){
                    vec4 p(k, Yrow, 0, 1);
                    vec4 c = Barycentric_Interpolation(p, V1, V2, V3);
                    QRgb value = qRgb(c[0], c[1], c[2]);
                    Scaleimg.setPixel(k, Yrow, value);
                }
            }
        }
    }
    for(int m=0; m<512*AA; m++)
        for(int n=0; n<512*AA; n++){
            img.setPixel(m/AA, n/AA, Scaleimg.pixel(m,n));
        }
    return img;
}

void Rasterizer::ClearScene()
{
    polygons.clear();
}

Line_Segment::Line_Segment(Vertex leftend, Vertex rightend){
    endpoints[0] = leftend;
    endpoints[1] = rightend;
    Y = leftend.pos[1] - rightend.pos[1];
    X = leftend.pos[0] - rightend.pos[0];
}

bool Line_Segment::IsIntersection(int Yr, Vertex &V1, Vertex &V2){
    if((Yr > endpoints[0].pos[1] && Yr > endpoints[1].pos[1]) || (Yr < endpoints[0].pos[1] && Yr < endpoints[1].pos[1]))
        return false;
    if(Y == 0){
        if(endpoints[0].pos[0] <= endpoints[1].pos[0]){
            V1 = endpoints[0];
            V2 = endpoints[1];
        }
        else{
            V1 = endpoints[1];
            V2 = endpoints[0];
        }
        return true;
    }
    if(X == 0){
        V1.pos[0] = endpoints[0].pos[0];
        V1.pos[1] = Yr;
        V1.color = endpoints[0].color;
        return true;
    }
    float slope = float(Y)/float(X);
    V1.pos[0] = endpoints[0].pos[0] + (Yr - endpoints[0].pos[1])/slope;
    V1.pos[1] = Yr;
    V1.color = endpoints[0].color;
    return true;
}



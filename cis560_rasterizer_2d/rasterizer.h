#pragma once
#include <polygon.h>
#include <QImage>

struct BoundingBox
{
    int Xl,Xr,Yl,Yr;
    BoundingBox(){Xl =0; Xr =0; Yl =0; Yr = 0;}
    void Computing(vec4 v1, vec4 v2, vec4 v3);
};

class Rasterizer
{
private:
    //This is the set of Polygons loaded from a JSON scene file
    std::vector<Polygon> polygons;
public:
    Rasterizer(const std::vector<Polygon> &polygons);
    QImage RenderScene();
    void ClearScene();
    QImage RenderLine();
    QImage AntiAliasing(float AA);
};

class Line_Segment
{
private:
    Vertex endpoints[2];
    int Y;
    int X;
public:
    Line_Segment(Vertex leftend, Vertex rightend);
    bool IsIntersection(int Yr, Vertex &V1, Vertex &V2);
};

vec4 Barycentric_Interpolation(vec4 p, Vertex S1, Vertex S2, Vertex S3);


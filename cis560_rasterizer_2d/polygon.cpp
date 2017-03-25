#include "polygon.h"
#include <math.h>

void Polygon::Triangulate()
{
    //TODO: Populate list of triangles
    Triangle Trag;
    int j = 0;
    while(!verts[j].Using) j++;
    if(j >= verts.size()-2) return;
    for(unsigned int i = j+1; i < verts.size()-1; i++){
        Trag.indices[0] = j;
        while(!verts[i].Using) i++;
        if(i >= verts.size() - 1) return;
        Trag.indices[1] = i;
        while(!verts[i+1].Using) i++;
        if(i >= verts.size()) return;
        Trag.indices[2] = i+1;
        Trag.depth = verts[i].pos[2];
        tris.push_back(Trag);
    }
}




//Creates a polygon from the input list of vertex positions and colors
Polygon::Polygon(const QString &name, const std::vector<vec4> &pos, const std::vector<vec4> &col)
    : tris(), name(name), texture(nullptr)
{
    for(unsigned int i = 0; i < pos.size(); i++)
    {
        verts.push_back(Vertex(pos[i], col[i], vec4(), vec4()));
    }
    Triangulate();
    ConcaveToConvex();
}

//Creates a regular polygon with a number of sides indicated by the "sides" input integer.
//All of its vertices are of color "color", and the polygon is centered at "pos".
//It is rotated about its center by "rot" degrees, and is scaled from its center by "scale" units
Polygon::Polygon(const QString &name, int sides, vec4 color, vec4 pos, float rot, vec4 scale)
    :tris(), name(name), texture(nullptr)
{
    vec4 v(0,1,0,1);
    float angle = 360.f / sides;
    for(int i = 0; i < sides; i++)
    {
        vec4 vert_pos = (mat4::translate(pos[0], pos[1], pos[2]) * mat4::rotate(rot, 0, 0, 1)
                * mat4::scale(scale[0], scale[1], scale[2]) * mat4::rotate(i * angle, 0,0,1) * v);
        verts.push_back(Vertex(vert_pos, color, vec4(), vec4()));
    }
    for(int i = 0; i < verts.size(); i ++){
        verts[i].pos[0] = int(verts[i].pos[0]);
        verts[i].pos[1] = int(verts[i].pos[1]);
        verts[i].pos[2] = int(verts[i].pos[2]);
    }
    Triangulate();
}

Polygon::Polygon(const QString &name)
    : tris(), name(name)
{}

Polygon::~Polygon()
{
    delete texture;
}

void Polygon::SetTexture(QImage *i)
{
    texture = i;
}

void Polygon::AddTriangle(const Triangle &t)
{
    tris.push_back(t);
}

void Polygon::AddVertex(const Vertex &v)
{
    verts.push_back(v);
}

void Polygon::ClearTriangles()
{
    tris.clear();
}

Triangle& Polygon::TriAt(unsigned int i)
{
    return tris[i];
}

Triangle Polygon::TriAt(unsigned int i) const
{
    return tris[i];
}

Vertex &Polygon::VertAt(unsigned int i)
{
    return verts[i];
}

Vertex Polygon::VertAt(unsigned int i) const
{
    return verts[i];
}

bool Polygon::IsConvex(Vertex &Vnum){
    std::vector<Vertex> vts;
    for(int j = 0; j < verts.size(); j++){
        if(verts[j].Using == true) vts.push_back(verts[j]);
    }
    int n = vts.size();
    for(int i = 0; i < n; i++){
        vec4 L1, L2;
        if(i == 0){
            L1 = vts[i].pos - vts[n-1].pos;
            L2 = vts[i].pos - vts[i+1].pos;
        }
        else if(i == n-1){
            L1 = vts[i].pos - vts[i-1].pos;
            L2 = vts[i].pos - vts[0].pos;
        }
        else{
            L1 = vts[i].pos - vts[i-1].pos;
            L2 = vts[i].pos - vts[i+1].pos;
        }
        if(cross(L1, L2)[2] < 0){
            Vnum = vts[i];
            return false;
        }
    }

    return true;
}

void Polygon::ConcaveToConvex(){
    Vertex Vnum;
    tris.clear();
    while(!IsConvex(Vnum)){
        int num = -1, i;
        for(i = 0; i < verts.size(); i++){
            if(Vnum.pos == verts[i].pos) num = i;
        }
        Triangle T1;
        T1.indices[0] = num;
        while(!verts[num+1].Using) num++;
        if(num >= verts.size()-2) break;
        T1.indices[1] = num+1;
        T1.indices[2] = num+2;
        T1.depth = verts[num].pos[2];
        tris.push_back(T1);
        verts[num+1].Using = false;
    }
    Triangulate();
}

float Polygon::Depth(){
    return verts[0].pos[2];
}


#pragma once
#include <vec4.h>
#include <mat4.h>
#include <vector>
#include <QString>
#include <QImage>

struct Vertex
{
    vec4 pos;    // The position of the vertex. In hw02, this is in pixel space.
    vec4 color;  // The color of the vertex. X corresponds to Red, Y corresponds to Green, and Z corresponds to Blue.
    vec4 normal; // The surface normal of the vertex (not yet used)
    vec4 uv;     // The texture coordinates of the vertex (not yet used)

    bool Using;

    Vertex(){vec4 v; pos = color = normal = uv = v; Using = true;}
    Vertex(vec4 p, vec4 c, vec4 n, vec4 u)
        : pos(p), color(c), normal(n), uv(u), Using(true)
    {}
};

struct Triangle
{
    //The indices of the Vertices that make up this triangle
    //The indices correspond to the std::vector of Vertices stored in the Polygon
    //which stores this Triangle
    Triangle(){indices[0] = indices[1] = indices[2] = depth = 0;}
    Triangle(int i, int j, int k, int d){indices[0] = i;indices[1] = j;indices[2] = k; depth = d;}
    unsigned int indices[3];
    unsigned int depth;
};

class Polygon
{
private:
    //TODO: Populate this list of triangles in Triangulate()
    std::vector<Triangle> tris;
    //The list of Vertices that define this polygon. This is already filled by the Polygon constructor.
    std::vector<Vertex> verts;
    //The name of this polygon, primarily to help you debug
    QString name;
    //The image that can be read to determine pixel color when used in conjunction with UV coordinates
    QImage* texture;

public:
    Polygon(const QString &name, const std::vector<vec4> &pos, const std::vector<vec4> &col);
    Polygon(const QString &name, int sides, vec4 color, vec4 pos, float rot, vec4 scale);
    Polygon(const QString &name);
    ~Polygon();

    //TODO: Complete the body of Triangulate() in polygon.cpp
    void Triangulate();
    //Copies the input QImage into this Polygon's texture
    void SetTexture(QImage*);

    //Various getter, setter, and adder functions
    void AddVertex(const Vertex&);
    void AddTriangle(const Triangle&);
    void ClearTriangles();
    int Vertex_Size() {return verts.size();}
    int Triangle_Size() { return tris.size();}

    bool IsConvex(Vertex &Vnum);
    void ConcaveToConvex();
    float Depth();

    Triangle& TriAt(unsigned int i);
    Triangle TriAt(unsigned int i) const;
    Vertex& VertAt(unsigned int i);
    Vertex VertAt(unsigned int i) const;
};



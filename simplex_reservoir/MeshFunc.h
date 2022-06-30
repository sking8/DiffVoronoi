//////////////////////////////////////////////////////////////////////////
// Mesh functions
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __AuxMeshFunc_h__
#define __AuxMeshFunc_h__
#include "Common.h"
template<int> class Box;

namespace MeshFunc{
    ////Normals
    Vector2 Normal(const Vector2& v1,const Vector2& v2);
    Vector1 Normal(const Vector2& p1,const Vector2& p2,const Vector2& p3);
    Vector3 Normal(const Vector3& v1,const Vector3& v2,const Vector3& v3);
    template<int d> Vector<real,d> Normal(const ArrayF<Vector<real,d>,d>& v);
    template<int d,int e_d> void Update_Normals(const Array<Vector<real,d> >& vertices,const Array<Vector<int,e_d> >& elements,Array<Vector<real,d> >& normals);
    template<class T_MESH> void Update_Normals(T_MESH& mesh)
	{if(mesh.normals==nullptr)mesh.normals.reset(new Array<Vector<real,T_MESH::Dim()> >());Update_Normals((*mesh.vertices),mesh.elements,(*mesh.normals));}
    template<int d,int e_d> Vector<real,d> Element_Normal(const Array<Vector<real,d> >& vertices,const Array<Vector<int,e_d> >& elements,const int i);
    
	////Area-weighted normals
    Vector2 Area_Weighted_Normal(const Vector2& v1,const Vector2& v2);
    Vector1 Area_Weighted_Normal(const Vector2& p1,const Vector2& p2,const Vector2& p3);
    Vector3 Area_Weighted_Normal(const Vector3& v1,const Vector3& v2,const Vector3& v3);
    template<int d,int e_d> Vector<real,d> Element_Area_Weighted_Normal(const Array<Vector<real,d> >& vertices,const Array<Vector<int,e_d> >& elements,const int i);

	////Element size and center
    real Tetrahedron_Volume(const Vector3& v0,const Vector3& v1,const Vector3& v2,const Vector3& v3);
	real Tetrahedron_Volume(const Vector2&,const Vector2&,const Vector2&,const Vector2&);	////TEMP: for template compilation only, should never be called
    real Tetrahedron_Volume(const ArrayF<Vector3,4>& tet);
    
	real Triangle_Area(const Vector2& v0,const Vector2& v1,const Vector2& v2);
    real Triangle_Area(const ArrayF<Vector2,3>& tri);
    real Triangle_Area(const Vector3& v0,const Vector3& v1,const Vector3& v2);
    real Triangle_Area(const ArrayF<Vector3,3>& tri);

	real Simplex_Size(const ArrayF<Vector2,2>& seg);
	real Simplex_Size(const ArrayF<Vector3,2>& seg);
    real Simplex_Size(const ArrayF<Vector3,4>& tet);
    real Simplex_Size(const ArrayF<Vector2,3>& tri);

	template<int d> real Element_Size(Array<Vector<real,d> >& v);
    template<int d,int e_d> real Element_Size(const Array<Vector<real,d> >& vertices,const Array<Vector<int,e_d> >& elements,const int i);
	////Center
	template<int d> Vector<real,d> Simplex_Center(const ArrayF<Vector<real,d>,d+1>& v);
	template<int d> Vector<real, d> Element_Center(const Array<Vector<real,d> >& v);
    template<int d,int e_d> Vector<real,d> Element_Center(const Array<Vector<real,d> >& vertices,const Array<Vector<int,e_d> >& elements,const int i);
	////Angle
	template<int d> void Triangle_Angles(const ArrayF<Vector<real,d>,3>& tri,ArrayF<real,3>& angles);
	template<int d> real Three_Point_Angle(const Vector<real,d>& a,const Vector<real,d>& b,const Vector<real,d>& c);
    ////Element edges
    template<class T_ARRAY> int Element_Edges(const Vector2i& v,T_ARRAY& edges);
    template<class T_ARRAY> int Element_Edges(const Vector3i& v,T_ARRAY& edges);
    template<class T_ARRAY> int Element_Edges(const Vector4i& v, T_ARRAY& edges);
    template<class T_ARRAY> int Element_Edges(const Vector4i& v,const int ds,T_ARRAY& edges);
    template<class T_ARRAY> int Quad_Edges(const Vector4i& v,T_ARRAY& edges);
	////Element faces
    template<class T_ARRAY> int Element_Faces(const Vector3i& v,T_ARRAY& faces);
    template<class T_ARRAY> int Element_Faces(const Vector4i& v,T_ARRAY& faces);
    template<class T_ARRAY> int Element_Faces(const Vector4i& v,const int ds,T_ARRAY& faces);
	////Mesh conversion (volume to surface)
    void Volumetric_Mesh_Surface(const Array<Vector4i>& vol_elements,Array<Vector3i>& surf_elements);       ////Compute a surface mesh from a tet mesh
	void Volumetric_Mesh_Surface(const Array<Vector3i>& vol_elements,Array<Vector2i>& surf_elements);		////Compute a segment mesh from a tri mesh
	////Edge length
	template<int d> real Average_Edge_Length(const Array<Vector<real,d> >& v,const Array<Vector<int,3> >& tri);
	template<int d> real Average_Edge_Length(const Array<Vector<real,d> >& v,const Array<Vector<int,2> >& seg);

	////Barycentric coordinates
	Vector2 Barycentric_Coord(const Vector2& p0,const Vector2& p1,const Vector2& p2,const Vector2& p);
	Vector3 Barycentric_Coord(const Vector3& p0,const Vector3& p1,const Vector3& p2,const Vector3& p3,const Vector3& p);	////TOIMPL
	template<int d> Vector<real,d> Barycentric_Coord(const ArrayF<Vector<real,d>,d+1>& vtx,const Vector<real,d>& p);
	template<class T,int d> T Barycentric_Interpolation(const ArrayF<T,d+1>& values,const Vector<real,d>& coord);
	template<int d> bool Inside(const ArrayF<Vector<real,d>,d+1>& vtx,const Vector<real,d>& p);

	////Mesh transformation
	template<int d> Box<d> Bounding_Box(const Array<Vector<real,d> >& vertices);
	template<int d> Box<d> Bounding_Box(const Vector<real,d>* vertices,int vn);
};
#endif
//////////////////////////////////////////////////////////////////////////
// Mesh functions
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <numeric>
#include <iostream>
#include "GeometryPrimitives.h"
#include "Constants.h"
#include "Hashtable.h"
#include "Grid.h"
#include "MeshFunc.h"
//#ifdef USE_TRI2D
//#include "Triangulation2D.h"	
//#endif

namespace MeshFunc{
    using namespace AuxFunc;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Normals
    Vector2 Normal(const Vector2& v1,const Vector2& v2)
	{Vector2 v12=(v2-v1).normalized();return Vector2(v12[1],-v12[0]);}
    
	Vector1 Normal(const Vector2& p1,const Vector2& p2,const Vector2& p3)
	{Vector2 p12=p2-p1;Vector2 p13=p3-p1;return Cross(p12,p13).normalized();}
    
	Vector3 Normal(const Vector3& p1,const Vector3& p2,const Vector3& p3)
	{return (p2-p1).cross(p3-p1).normalized();}

    template<> Vector<real,2> Normal<2>(const ArrayF<Vector<real,2>,2>& v){return Normal(v[0],v[1]);}
    template<> Vector<real,3> Normal<3>(const ArrayF<Vector<real,3>,3>& v){return Normal(v[0],v[1],v[2]);}

    template<> void Update_Normals<2,2>(const Array<Vector2>& vertices,const Array<Vector2i>& elements,Array<Vector2>& normals)
    {
        normals.resize(vertices.size(),Vector2::Zero());
        for(const auto& v:elements){Vector2 n=Normal(vertices[v[0]],vertices[v[1]]);for(int j=0;j<2;j++){normals[v[j]]+=n;}}
        for(auto& n:normals){n.normalize();}
    }
    template<> void Update_Normals<3,3>(const Array<Vector3>& vertices,const Array<Vector3i>& elements,Array<Vector3>& normals)
    {
        normals.resize(vertices.size(),Vector3::Zero());
        for(const auto& v:elements){Vector3 n=Normal(vertices[v[0]],vertices[v[1]],vertices[v[2]]);for(int j=0;j<3;j++){normals[v[j]]+=n;}}
        for(auto& n:normals){n.normalize();}
    }

    template<> Vector2 Element_Normal<2,2>(const Array<Vector2>& vertices,const Array<Vector2i>& elements,const int i)
    {return Normal(vertices[elements[i][0]],vertices[elements[i][1]]);}
    
	template<> Vector3 Element_Normal<3,3>(const Array<Vector3>& vertices,const Array<Vector3i>& elements,const int i)
    {return Normal(vertices[elements[i][0]],vertices[elements[i][1]],vertices[elements[i][2]]);}

    template<> Vector2 Element_Area_Weighted_Normal<2,2>(const Array<Vector2>& vertices,const Array<Vector2i>& elements,const int i)
    {return Area_Weighted_Normal(vertices[elements[i][0]],vertices[elements[i][1]]);}
    
	template<> Vector3 Element_Area_Weighted_Normal<3,3>(const Array<Vector3>& vertices,const Array<Vector3i>& elements,const int i)
    {return Area_Weighted_Normal(vertices[elements[i][0]],vertices[elements[i][1]],vertices[elements[i][2]]);}

	Vector2 Area_Weighted_Normal(const Vector2& v1,const Vector2& v2)
	{Vector2 v12=(v2-v1).normalized();return Vector2(v12[1],-v12[0]);}

	Vector1 Area_Weighted_Normal(const Vector2& p1,const Vector2& p2,const Vector2& p3)
	{Vector2 p12=p2-p1;Vector2 p13=p3-p1;return Cross(p12,p13);}

	Vector3 Area_Weighted_Normal(const Vector3& p1,const Vector3& p2,const Vector3& p3)
	{return (p2-p1).cross(p3-p1);}

	////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Element size and center
    real Tetrahedron_Volume(const Vector3& v0,const Vector3& v1,const Vector3& v2,const Vector3& v3){return (real)one_sixth*abs((v1-v0).cross(v2-v0).dot(v3-v0));}
    real Tetrahedron_Volume(const ArrayF<Vector3,4>& tet){return Tetrahedron_Volume(tet[0],tet[1],tet[2],tet[3]);}
	real Tetrahedron_Volume(const Vector2&,const Vector2&,const Vector2&,const Vector2&){return (real)0;}

	real Triangle_Area(const Vector2& v0,const Vector2& v1,const Vector2& v2){Vector2 v01=v1-v0;Vector2 v02=v2-v0;return (real).5*abs(Cross(v01,v02)[0]);}
    real Triangle_Area(const ArrayF<Vector2,3>& tri){return Triangle_Area(tri[0],tri[1],tri[2]);}
    real Triangle_Area(const Vector3& v0,const Vector3& v1,const Vector3& v2){return (real).5*((v1-v0).cross(v2-v0)).norm();}
    real Triangle_Area(const ArrayF<Vector3,3>& tri){return Triangle_Area(tri[0],tri[1],tri[2]);}
	
	real Simplex_Size(const ArrayF<Vector2,2>& seg){return (seg[0]-seg[1]).norm();}
	real Simplex_Size(const ArrayF<Vector3,2>& seg){return (seg[0]-seg[1]).norm();}
	real Simplex_Size(const ArrayF<Vector2,3>& tri){return Triangle_Area(tri);}
	real Simplex_Size(const ArrayF<Vector3,4>& tet){return Tetrahedron_Volume(tet);}

    template<> real Element_Size<2>(Array<Vector2>& v){return (v[1]-v[0]).norm();}
    template<> real Element_Size<3>(Array<Vector3>& v){return (real).5*((v[1]-v[0]).cross(v[2]-v[0])).norm();}

    template<> real Element_Size<2,2>(const Array<Vector2>& vertices,const Array<Vector2i>& elements,const int i)
    {return (vertices[elements[i][1]]-vertices[elements[i][0]]).norm();}
    template<> real Element_Size<3,3>(const Array<Vector3>& vertices,const Array<Vector3i>& elements,const int i)
    {return (real).5*((vertices[elements[i][1]]-vertices[elements[i][0]]).cross(vertices[elements[i][2]]-vertices[elements[i][0]])).norm();}

	template<int d> Vector<real,d> Simplex_Center(const ArrayF<Vector<real,d>,d+1>& v)
	{
		Vector<real,d> center=Vector<real,d>::Zero();
		for(int i=0;i<v.size();i++){center+=v[i];}center/=(real)(v.size());return center;		
	}
	template Vector2 Simplex_Center<2>(const ArrayF<Vector2,3>&);
	template Vector3 Simplex_Center<3>(const ArrayF<Vector3,4>&);

    template<> Vector<real,2> Element_Center<2>(const Array<Vector2>& v){return (real).5*(v[0]+v[1]);}
    template<> Vector<real,3> Element_Center<3>(const Array<Vector3>& v){return (real)one_third*(v[0]+v[1]+v[2]);}

    template<int d,int e_d> Vector<real,d> Element_Center(const Array<Vector<real,d> >& vertices,const Array<Vector<int,e_d> >& elements,const int i)
    {Vector<real,d> center=Vector<real,d>::Zero();for(int j=0;j<e_d;j++)center+=vertices[elements[i][j]];center/=(real)e_d;return center;}

#define comma ,
#define Inst_Helper(d,e_d) \
template Vector<real,d> Element_Center<d,e_d>(const Array<Vector<real,d> >&,const Array<Vector<int,e_d> >&,const int)
Inst_Helper(2,2);Inst_Helper(2,3);Inst_Helper(3,3);Inst_Helper(3,4);
#undef Inst_Helper
#undef comma

	template<int d> void Triangle_Angles(const ArrayF<Vector<real,d>,3>& tri,ArrayF<real,3>& angles)
	{
		Vector<real,d> v01=(tri[1]-tri[0]).normalized();Vector<real,d> v02=(tri[2]-tri[0]).normalized();angles[0]=Angle_Between(v01,v02);
		Vector<real,d> v12=(tri[2]-tri[1]).normalized();angles[1]=Angle_Between(v12,-v01);
		angles[2]=pi-angles[0]-angles[1];	
	}
	template void Triangle_Angles<2>(const ArrayF<Vector2,3>&,ArrayF<real,3>&);
	template void Triangle_Angles<3>(const ArrayF<Vector3,3>&,ArrayF<real,3>&);
	
	template<int d> real Three_Point_Angle(const Vector<real,d>& a,const Vector<real,d>& b,const Vector<real,d>& c)
	{
		Vector<real,d> ba=(a-b).normalized();Vector<real,d> bc=(c-b).normalized();return Angle_Between(bc,ba);
	}
	template real Three_Point_Angle<2>(const Vector2&,const Vector2&,const Vector2&);
	template real Three_Point_Angle<3>(const Vector3&,const Vector3&,const Vector3&);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Element edges
    template<class T_ARRAY> int Element_Edges(const Vector2i& v,T_ARRAY& edges)
    {edges[0]=v;return 1;}
    template<class T_ARRAY> int Element_Edges(const Vector3i& v,T_ARRAY& edges)
    {edges[0]=Vector2i(v[0],v[1]);edges[1]=Vector2i(v[1],v[2]);edges[2]=Vector2i(v[2],v[0]);return 3;}
    template<class T_ARRAY> int Element_Edges(const Vector4i& v,T_ARRAY& edges)
    {edges[0]=Vector2i(v[0],v[1]);edges[1]=Vector2i(v[0],v[2]);edges[2]=Vector2i(v[0],v[3]);edges[3]=Vector2i(v[1],v[2]);edges[4]=Vector2i(v[2],v[3]);edges[5]=Vector2i(v[3],v[1]);return 6;}
    template<class T_ARRAY> int Element_Edges(const Vector4i& v,const int ds,T_ARRAY& edges)
    {switch(ds){case 2:return Element_Edges(Vector2i(v[0],v[1]),edges);case 3:return Element_Edges(Vector3i(v[0],v[1],v[2]),edges);case 4:return Element_Edges(v,edges);default:return 0;}}

    template<class T_ARRAY> int Quad_Edges(const Vector4i& v,T_ARRAY& edges)
    {edges[0]=Vector2i(v[0],v[1]);edges[1]=Vector2i(v[1],v[2]);edges[2]=Vector2i(v[2],v[3]);edges[3]=Vector2i(v[3],v[0]);return 4;}

    template<class T_ARRAY> int Element_Faces(const Vector3i& v,T_ARRAY& faces)
    {faces[0]=v;return 1;}
    template<class T_ARRAY> int Element_Faces(const Vector4i& v,T_ARRAY& faces)
    {faces[0]=Vector3i(v[0],v[1],v[3]);faces[1]=Vector3i(v[1],v[2],v[3]);faces[2]=Vector3i(v[2],v[0],v[3]);faces[3]=Vector3i(v[1],v[0],v[2]);return 4;}
    template<class T_ARRAY> int Element_Faces(const Vector4i& v,const int ds,T_ARRAY& faces)
    {switch(ds){case 3:return Element_Faces(Vector3i(v[0],v[1],v[2]),faces);case 4:return Element_Faces(v,faces);default:return 0;}}

#define comma ,
#define Inst_Helper(T_VEC,T_ARRAY) \
template int Element_Edges<T_ARRAY >(const T_VEC& v,T_ARRAY& edges)
    Inst_Helper(Vector2i,Array<Vector2i>);Inst_Helper(Vector2i,ArrayF<Vector2i comma 1>);Inst_Helper(Vector2i,ArrayF<Vector2i comma 6>);
    Inst_Helper(Vector3i,Array<Vector2i>);Inst_Helper(Vector3i,ArrayF<Vector2i comma 3>);Inst_Helper(Vector3i,ArrayF<Vector2i comma 6>);
    Inst_Helper(Vector4i,Array<Vector2i>);Inst_Helper(Vector4i,ArrayF<Vector2i comma 6>);
#undef Inst_Helper
#define Inst_Helper(T_VEC,T_ARRAY) \
template int Element_Edges<T_ARRAY >(const T_VEC& v,const int ds,T_ARRAY& edges)
    Inst_Helper(Vector4i,Array<Vector2i>);Inst_Helper(Vector4i,ArrayF<Vector2i comma 6>);
#undef Inst_Helper
#define Inst_Helper(T_VEC,T_ARRAY) \
template int Quad_Edges<T_ARRAY >(const T_VEC& v,T_ARRAY& edges)
    Inst_Helper(Vector4i,Array<Vector2i>);Inst_Helper(Vector4i,ArrayF<Vector2i comma 4>);
#undef Inst_Helper
#define Inst_Helper(T_VEC,T_ARRAY) \
template int Element_Faces<T_ARRAY >(const T_VEC& v,T_ARRAY& faces)
    Inst_Helper(Vector3i,Array<Vector3i>);Inst_Helper(Vector3i,ArrayF<Vector3i comma 1>);Inst_Helper(Vector3i,ArrayF<Vector3i comma 4>);
    Inst_Helper(Vector4i,Array<Vector3i>);Inst_Helper(Vector4i,ArrayF<Vector3i comma 4>);
#undef Inst_Helper
#define Inst_Helper(T_VEC,T_ARRAY) \
template int Element_Faces<T_ARRAY >(const T_VEC& v,const int ds,T_ARRAY& faces)
    Inst_Helper(Vector4i,Array<Vector3i>);Inst_Helper(Vector4i,ArrayF<Vector3i comma 4>);
#undef Inst_Helper
#undef comma

	//////////////////////////////////////////////////////////////////////////
	////Mesh edges

	template<int d> real Average_Edge_Length(const Array<Vector<real,d> >& X,const Array<Vector<int,3> >& tri)
	{
		real length=(real)0;int n=0;
		for(auto i=0;i<tri.size();i++){const Vector3i& v=tri[i];
			length+=((X[v[0]]-X[v[1]]).norm()+(X[v[1]]-X[v[2]]).norm()+(X[v[2]]-X[v[0]]).norm());n+=3;}		
		if(n>0){return length/(real)n;}return length;
	}
	template real Average_Edge_Length<2>(const Array<Vector<real,2> >&,const Array<Vector3i>&);
	template real Average_Edge_Length<3>(const Array<Vector<real,3> >&,const Array<Vector3i>&);
	
	template<int d> real Average_Edge_Length(const Array<Vector<real,d> >& X,const Array<Vector<int,2> >& seg)
	{
		real length=(real)0;int n=(int)seg.size();
		for(auto i=0;i<n;i++){const Vector2i& v=seg[i];
			length+=((X[v[0]]-X[v[1]]).norm());}		
		if(n>0)return length/(real)n;
		return length;
	}
	template real Average_Edge_Length<2>(const Array<Vector<real,2> >&,const Array<Vector2i>&);
	template real Average_Edge_Length<3>(const Array<Vector<real,3> >&,const Array<Vector2i>&);
	
    //////////////////////////////////////////////////////////////////////////
    ////Mesh converter
    void Volumetric_Mesh_Surface(const Array<Vector4i>& vol_elements,Array<Vector3i>& surf_elements)
    {
        ArrayF<Vector3i,4> faces;Hashtable<Vector3i,int> boundary_face_hashtable;
        for(int i=0;i<(int)vol_elements.size();i++){const Vector4i& vtx=vol_elements[i];
            Element_Faces(vtx,faces);for(int i=0;i<4;i++){const Vector3i& key=Sorted(faces[i]);
                auto find=boundary_face_hashtable.find(key);
                if(find!=boundary_face_hashtable.end()){boundary_face_hashtable.erase(key);}
                else{int sign=Is_Same_Order(key,faces[i])?1:-1;boundary_face_hashtable.insert(std::make_pair(key,sign));}}}
        for(const auto& iter:boundary_face_hashtable){
	    Vector3i tri=iter.first;
            if(iter.second<0)tri=Reversed(iter.first);
	    surf_elements.push_back(tri);
	}
    }

	void Volumetric_Mesh_Surface(const Array<Vector3i>& vol_elements,Array<Vector2i>& surf_elements)
	{
		ArrayF<Vector2i,3> edges;Hashtable<Vector2i,int> boundary_edge_hashtable;
		for(int i=0;i<(int)vol_elements.size();i++){const Vector3i& vtx=vol_elements[i];Element_Edges(vtx,edges);
			for(int i=0;i<3;i++){const Vector2i& key=Sorted(edges[i]);
				auto find=boundary_edge_hashtable.find(key); 
				if(find!=boundary_edge_hashtable.end()){boundary_edge_hashtable.erase(key);}
				else{int sign=(key==edges[i])?1:-1;boundary_edge_hashtable.insert(std::make_pair(key,sign));}}}
		for(const auto& iter:boundary_edge_hashtable){Vector2i edge=iter.first;
			if(iter.second<0)edge=Reversed(edge);
			surf_elements.push_back(edge);}
	}

	Vector2 Barycentric_Coord(const Vector2& p0,const Vector2& p1,const Vector2& p2,const Vector2& p)
	{
		Matrix2 m;m<<p0[0]-p[0],p1[0]-p2[0],p0[1]-p2[1],p1[1]-p2[1];return m.inverse()*(p-p2);
	}

	Vector3 Barycentric_Coord(const Vector3& p0,const Vector3& p1,const Vector3& p2,const Vector3& p3,const Vector3& p)
	{
		////TOIMPL
		return Vector3::Zero();
	}

	template<> Vector2 Barycentric_Coord<2>(const ArrayF<Vector2,3>& vtx,const Vector2& p)
	{
		return Barycentric_Coord(vtx[0],vtx[1],vtx[2],p);
	}

	template<> Vector3 Barycentric_Coord<3>(const ArrayF<Vector3,4>& vtx,const Vector3& p)
	{
		return Barycentric_Coord(vtx[0],vtx[1],vtx[2],vtx[3],p);
	}

	template<class T,int d> T Barycentric_Interpolation(const ArrayF<T,d+1>& values,const Vector<real,d>& coord)
	{
		T v=Zero<T>();real c=(real)1;
		for(int i=0;i<d;i++){c-=coord[i];v+=values[i]*coord[i];}v+=values[d]*c;return v;
	}
	#define Inst_Helper(T,d)	\
	template T Barycentric_Interpolation<T,d>(const ArrayF<T,d+1>&,const Vector<real,d>&);
	Inst_Helper(real,2);Inst_Helper(real,3);Inst_Helper(Vector2,2);Inst_Helper(Vector3,3);
	#undef Inst_Helper
	
	template<int d> bool Inside(const ArrayF<Vector<real,d>,d+1>& vtx,const Vector<real,d>& p)
	{
		Vector<real,d> bc=Barycentric_Coord<d>(vtx,p);real c=1;for(int i=0;i<d;i++){if(bc[i]<0||bc[i]>1)return false;c-=bc[i];}if(c<0||c>1)return false;return true;
	}
	#define Inst_Helper(d)	\
	template bool Inside<d>(const ArrayF<Vector<real,d>,d+1>&,const Vector<real,d>&);
	Inst_Helper(2);Inst_Helper(3);
	#undef Inst_Helper
	
	//////////////////////////////////////////////////////////////////////////
	////Mesh transformation
	template<int d> Box<d> Bounding_Box(const Array<Vector<real,d> >& vertices)
	{
		Box<d> box=Box<d>::Infi_Min();
		for(auto& v:vertices){box.min_corner=box.min_corner.cwiseMin(v);box.max_corner=box.max_corner.cwiseMax(v);}return box;
	}
	template Box<2> Bounding_Box<2>(const Array<Vector2>&);
	template Box<3> Bounding_Box<3>(const Array<Vector3>&);

	template<int d> Box<d> Bounding_Box(const Vector<real,d>* vertices,int vn)
	{
		Box<d> box=Box<d>::Infi_Min();
		for(int i=0;i<vn;i++){box.min_corner=box.min_corner.cwiseMin(vertices[i]);box.max_corner=box.max_corner.cwiseMax(vertices[i]);}return box;
	}
	template Box<2> Bounding_Box<2>(const Vector2*,int);
	template Box<3> Bounding_Box<3>(const Vector3*, int);
};

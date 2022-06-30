//////////////////////////////////////////////////////////////////////////
// Soft body mass spring
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __SoftBodyMassSpring_h__
#define __SoftBodyMassSpring_h__
#include "Hashtable.h"
#include "Particles.h"
#include "SparseFunc.h"

template<int d> class SoftBodyMassSpring
{Typedef_VectorDii(d);Typedef_MatrixD(d);
public:
	real ks_0=(real)1e5;
	real kd_0=(real)1e3;
	real ks_bending_0=(real)1e4;
	real kd_bending_0=(real)1e1;
protected:
	Particles<d>* particles_ptr=nullptr;
public:	
	Particles<d>& particles;
	bool own_data=false;

	Array<Vector2i> springs;
	Array<Vector4i> dihedrals;
	Array<real> rest_length;
	Array<real> rest_angle; //// equals to sin((pi-theta)/2) pi-theta is the dihedra
	Array<Array<int> > particle_springs;		////incident springs of each particle, to accelerate the spring calculations
	bool use_bending_spring_line=false;			////whether add bending springs for 1D geometry
	bool use_bending_spring_surface=false;		////whether add bending springs for 2D geometry
	bool use_bending_force=false;

	int p_num=0;
	int s_num=0;

	bool use_varying_k=false;
	Array<real> ks;
	Array<real> kd;
	Array<real> ks_bending;
	Array<real> kd_bending;
	Array<VectorD> ext_force;

	Hashtable<int,VectorD> psi_D_values;		////displacements
	Hashtable<int,VectorD> psi_N_values;		////forces
	VectorD g=VectorD::Unit(1)*(real)-9.8;
	bool use_body_force=true;
	
	bool use_implicit=true;
	bool use_newmark=false;
	bool use_newmark2=false;
	SparseMatrixT K;
	VectorX u,b;

	SoftBodyMassSpring(Particles<d>* _p=nullptr)
		:particles_ptr(_p==nullptr?new Particles<d>():_p),particles(*particles_ptr){if(_p==nullptr)own_data=true;}
	~SoftBodyMassSpring(){if(own_data)delete particles_ptr;}

	virtual void Initialize()	////assuming springs and particles are initialized
	{
		s_num=(int)springs.size();
		p_num=particles.Size();
		
		if(!use_varying_k){
			ks.resize(s_num,ks_0);
			kd.resize(s_num,kd_0);}

		if (use_bending_spring_line) {
			Construct_Bending_Spring_Line();} //// extend springs, ks, kd and modify s_num
		else if (use_bending_spring_surface) {
			Construct_Bending_Spring_Surface();} //// extend springs, ks, kd and modify s_num

		rest_length.resize(s_num);
		#pragma omp parallel for
		for(int i=0;i<s_num;i++){const Vector2i& s=springs[i];
			rest_length[i]=(particles.X(s[0])-particles.X(s[1])).norm();}
		
		ext_force.resize(p_num,VectorD::Zero());

		particle_springs.resize(p_num);
		for(int i=0;i<s_num;i++){
			const Vector2i& s=springs[i];
			particle_springs[s[0]].push_back(i);
			particle_springs[s[1]].push_back(i);}

		if(use_bending_force){
			if constexpr (d==2) {
				std::cout<<"dihedral constraint doesn't supported in 2D"<<std::endl;
				return;
			} else {
				for(int i=0;i<s_num;i++){
					const Vector2i& s=springs[i];
					int u3=s[0],u4=s[1];
					Hashset<int> u3_neighbors;
					Array<int> public_neighbors;
					for (int j=0;j<particle_springs[u3].size();j++) {
						int another_node=springs[particle_springs[u3][j]][0];
						if (another_node==u3) another_node=springs[particle_springs[u3][j]][1];
						u3_neighbors.insert(another_node);
					}
					for (int j=0;j<particle_springs[u4].size();j++) {
						int another_node=springs[particle_springs[u4][j]][0];
						if (another_node==u4) another_node=springs[particle_springs[u4][j]][1];
						if(u3_neighbors.find(another_node)!=u3_neighbors.end()){
							public_neighbors.push_back(another_node);}
					}
					if(public_neighbors.size()==2){
						dihedrals.push_back(Vector4i(public_neighbors[0],public_neighbors[1],u3,u4));}
				}

				rest_angle.resize(dihedrals.size());
				for (int i=0;i<dihedrals.size();i++){
					const Vector4i& dih=dihedrals[i];
					VectorD x1=particles.X(dih[0]),x2=particles.X(dih[1]),x3=particles.X(dih[2]),x4=particles.X(dih[3]);
					VectorD n1=(x1-x3).cross(x1-x4),n2=(x2-x4).cross(x2-x3),e=x4-x3;
					n1/=n1.norm();n2/=n2.norm();
					real angle=sqrt(((real)1.0-n1.dot(n2))/2);
					if (n1.cross(n2).dot(e)<0) angle=-angle;
					rest_angle[i]=angle;
				}

				if(!use_varying_k){
					ks_bending.resize(dihedrals.size(),ks_bending_0);
					kd_bending.resize(dihedrals.size(),kd_bending_0);}}
		}

		if(use_implicit || use_newmark || use_newmark2) Initialize_Implicit_K();
	}

	virtual void Advance(const real dt)
	{
		if (use_newmark2) {
			real cfl=(real).1;int n = (int)(1./cfl); ////TOFIX: fix cfl  
			for (int i = 0; i < n; i++)
				Advance_Newmark2(dt*cfl);
		}
		else if (use_newmark) {
			real cfl=(real).01;int n = (int)(1./cfl); ////TOFIX: fix cfl  
			for (int i = 0; i < n; i++)
				Advance_Newmark(dt*cfl);
		}
		else if(!use_implicit){
			real cfl=(real).005;int n=(int)(1./cfl);	////TOFIX: fix cfl
			for(int i=0;i<n;i++)Advance_Midpoint(dt*cfl);}
		else Advance_Implicit(dt);

		//Advance_Explicit_Euler(dt);
		//Advance_Midpoint(dt);
	}

	//// add springs between the other vertices of every two edges which contains one same vertex.
	void Construct_Bending_Spring_Line() {
		Hashset<Vector2i> spring_hashset; 
		for (int i = 0; i < s_num; i++) {
			spring_hashset.insert(Sorted(springs[i]));}

		// too brute force
		for (int i = 0; i < p_num; i++) {
			Array<int> adjacent_vertices;
			for (int j = 0; j < p_num; j++) {
				Vector2i spring(i, j);spring = Sorted(spring);
				if (spring_hashset.find(spring) != spring_hashset.end()) {
					adjacent_vertices.push_back(j);
				}
			}
			if (adjacent_vertices.size()==2) {
				Vector2i spring(adjacent_vertices[0], adjacent_vertices[1]);
				if (spring_hashset.find(Sorted(spring)) == spring_hashset.end()) {
					s_num++;springs.push_back(spring);kd.push_back(kd_bending_0);ks.push_back(ks_bending_0);}
			}
		}
	}

	//// add springs between the corner vertices of every two adjacent triangles 
	void Construct_Bending_Spring_Surface() {
		if constexpr (d==2){
			return;
		}
		Hashset<Vector2i> spring_hashset;
		for (int i = 0; i < s_num; i++) {
			spring_hashset.insert(Sorted(springs[i]));}

		//// too brute force
		for (int s_i = 0; s_i < s_num; s_i++) {
			Array<int> adjacent_vertices;
			int i = springs[s_i][0], j = springs[s_i][1];
			for (int k = 0; k < p_num; k++) {
				Vector2i e1(i, k);Vector2i e2(j, k);
				if (spring_hashset.find(Sorted(e1)) != spring_hashset.end() && 
					spring_hashset.find(Sorted(e2)) != spring_hashset.end()) {
					adjacent_vertices.push_back(k);
				}
			}

			if (adjacent_vertices.size() == 2) {
				Vector2i spring(adjacent_vertices[0], adjacent_vertices[1]);
				if (spring_hashset.find(Sorted(spring)) == spring_hashset.end()) {
					s_num++;springs.push_back(spring);kd.push_back(kd_bending_0);ks.push_back(ks_bending_0);}
			}
		}
	}

	void Update_Forces(const Array<VectorD>& X,const Array<VectorD>& V,Array<VectorD>& F)
	{
		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			F[i]=ext_force[i];
			if(use_body_force)F[i]+=particles.M(i)*g;
			int s_n=particle_springs[i].size();
			for(int j=0;j<s_n;j++){
				int s=particle_springs[i][j];
				VectorD f=Spring_Force(X,V,s);
				F[i]+=(springs[s][0]==i?(real)1:(real)-1)*f;}}

		////sequential update force
		//for(int i=0;i<(int)springs.size();i++){
		//	VectorD f=Spring_Force(X,V,i);
		//	F[springs[i][0]]+=f;F[springs[i][1]]-=f;}
		//if(use_body_force)for(int i=0;i<p_num;i++){
		//	F[i]+=particles.M(i)*g;}

		//// Reference: Simulation of Clothing with Folds and Wrinkles
		if (use_bending_force){
			if constexpr (d==2) {
				std::cout<<"dihedral constraint doesn't supported in 2D"<<std::endl;
			}
			else {
				int d_num=dihedrals.size();
				for(int d_i=0;d_i<d_num;d_i++){
					const Array<VectorD>& forces=Bending_Force(X,V,d_i);
					const Vector4i dih=dihedrals[d_i];
					for(int i=0;i<4;i++) {
						F[dih[i]]+=forces[i];
					}
				}
			}
		}

		for(auto p:psi_N_values){int idx=p.first;const VectorD& f=p.second;F[idx]+=f;}
	}

	VectorD Spring_Force(const Array<VectorD>& X,const Array<VectorD>& V,const int s_idx)	////return f_ij=f_s+f+d
	{
		int i=springs[s_idx][0];int j=springs[s_idx][1];
		VectorD dir=X[j]-X[i];
		real length=dir.norm();
		dir/=length;
		VectorD rel_dir=V[j]-V[i];
		VectorD f_s=ks[s_idx]*(length-rest_length[s_idx])*dir;
		VectorD f_d=kd[s_idx]*rel_dir.dot(dir)*dir;
		return f_s+f_d;
	}

	static VectorD Spring_Force(const VectorD& Xi,const VectorD& Xj,const VectorD& Vi,const VectorD& Vj,const real rest_length,const real ks,const real kd)	////return f_ij=f_s+f+d
	{
		VectorD dir=Xj-Xi;
		real length=dir.norm();
		dir/=length;
		VectorD rel_dir=Vj-Vi;
		VectorD f_s=ks*(length-rest_length)*dir;
		VectorD f_d=kd*rel_dir.dot(dir)*dir;
		return f_s+f_d;
	}

	Array<VectorD> Bending_Force(const Array<VectorD>& X,const Array<VectorD>& V,const int d_idx)
	{
		Array<VectorD> forces;
		const Vector4i& dih=dihedrals[d_idx];
		VectorD x1=X[dih[0]],x2=X[dih[1]],x3=X[dih[2]],x4=X[dih[3]];
		VectorD N1=(x1-x3).cross(x1-x4),N2=(x2-x4).cross(x2-x3),E=x4-x3;
		real N1_sqr=N1.dot(N1),N2_sqr=N2.dot(N2),E_norm=E.norm();
		Array<VectorD> u;
		u.push_back(E_norm*N1/N1_sqr);
		u.push_back(E_norm*N2/N2_sqr);
		u.push_back((x1-x4).dot(E)/E_norm*N1/N1_sqr+(x2-x4).dot(E)/E_norm*N2/N2_sqr);
		u.push_back(-(x1-x3).dot(E)/E_norm*N1/N1_sqr-(x2-x3).dot(E)/E_norm*N2/N2_sqr);

		N1/=N1.norm();N2/=N2.norm();
		real sin_angle=sqrt(std::max(((real)1-N1.dot(N2))/2, (real)0.0));
		if(N1.cross(N2).dot(E)<0) sin_angle=-sin_angle;
		real k_para=ks_bending[d_idx]*E.dot(E)/(N1.norm()+N2.norm())*(sin_angle-rest_angle[d_idx]);
		for (int i=0;i<4;i++){
			forces.push_back(k_para*u[i]);
		}
		return forces;
	}

	void Set_Psi_D(const int p,const VectorD v=VectorD::Zero()){psi_D_values[p]=v;}
	bool Is_Psi_D(const int p){return psi_D_values.find(p)!=psi_D_values.end();}
	void Set_Psi_N(const int p,const VectorD f){psi_N_values[p]=f;}
	bool Is_Psi_N(const int p){return psi_N_values.find(p)!=psi_N_values.end();}

	void Enforce_Boundary_Conditions(Array<VectorD>& V,Array<VectorD>& F)
	{
		for(auto p:psi_D_values){int idx=p.first;const VectorD& v=p.second;V[idx]=v;F[idx]=VectorD::Zero();}
	}

	virtual real CFL(){return (real).002;}

	void Advance_Explicit_Euler(const real dt)
	{
		Update_Forces(particles.XRef(),particles.VRef(),particles.FRef());
		Enforce_Boundary_Conditions(particles.VRef(),particles.FRef());

		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			particles.V(i)+=particles.F(i)/particles.M(i)*dt;
			particles.X(i)+=particles.V(i)*dt;}	
	}

	void Advance_Midpoint(const real dt)
	{
		Update_Forces(particles.XRef(),particles.VRef(),particles.FRef());	////update f_n and store it in particles.F
		Enforce_Boundary_Conditions(particles.VRef(),particles.FRef());
		Array<VectorD> v_half(p_num,VectorD::Zero());
		Array<VectorD> x_half(p_num,VectorD::Zero());

		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			x_half[i]=particles.X(i)+dt*(real).5*particles.V(i);	////x_half=x_n+half_dt*v_n
			v_half[i]=particles.V(i)+dt*(real).5*particles.F(i)/particles.M(i);}	////v_half=v_n+half_dt*f_n/m

		Update_Forces(x_half,v_half,particles.FRef());			////update f_half and store it in particles.F
		Enforce_Boundary_Conditions(v_half,particles.FRef());

		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			particles.X(i)+=dt*v_half[i];						////x_{n+1}=x_n+dt*v_half
			particles.V(i)+=dt*particles.F(i)/particles.M(i);}	////v_{n+1}=v_n+dt*f_half/m
	}

	////Implicit time integration
	void Advance_Implicit(const real dt)
	{
		Update_Forces(particles.XRef(),particles.VRef(),particles.FRef());
		Enforce_Boundary_Conditions(particles.VRef(),particles.FRef());

		Update_Implicit_K(particles.XRef(),particles.VRef(),particles.FRef(),dt);

		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			VectorD v=particles.V(i);
			for(int j=0;j<d;j++)u[i*d+j]=v[j];}	////set initial guess to be the velocity from the last time step
		SparseSolver::Conjugate_Gradient(K,u,b);

		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			VectorD v;for(int j=0;j<d;j++)v[j]=u[i*d+j];
			particles.V(i)=v;
			particles.X(i)+=particles.V(i)*dt;}
	}

	//// Mixed Explicit/Implicit Time Integration(nemark scheme) 
	//// Reference: Simulation of Clothing with Folds and Wrinkles 
	void Advance_Newmark(const real dt)
	{
		Update_Forces(particles.XRef(),particles.VRef(),particles.FRef());	////update f_n and store it in particles.F
		Enforce_Boundary_Conditions(particles.VRef(),particles.FRef());
		Array<VectorD> v_half_explicit(p_num,VectorD::Zero());

		Update_Forces(particles.XRef(), particles.VRef(), particles.FRef());
		Enforce_Boundary_Conditions(particles.VRef(), particles.FRef());

		for (int i = 0; i < p_num; i++) {
			v_half_explicit[i] = particles.V(i) + dt*(real).5*particles.F(i)/particles.M(i);
			particles.X(i) += v_half_explicit[i] * dt;
		}

		Update_Forces(particles.XRef(), v_half_explicit, particles.FRef());			////update f_half and store it in particles.F
		Enforce_Boundary_Conditions(v_half_explicit, particles.FRef());

		Update_Newmark_K(particles.XRef(), v_half_explicit, particles.FRef(), dt*(real).5);

		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			VectorD v=v_half_explicit[i];
			for(int j=0;j<d;j++)u[i*d+j]=v[j];}	////set initial guess to be the velocity from the last half time step
		SparseSolver::Conjugate_Gradient(K,u,b);

		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			VectorD v;for(int j=0;j<d;j++)v[j]=u[i*d+j];
			particles.V(i)=v;}
	}

	// Modified newmark scheme in the same paper
	void Advance_Newmark2(const real dt) {
		Array<VectorD> v_half(p_num, VectorD::Zero());
		Array<VectorD> x_old(p_num, VectorD::Zero());

		for (int i = 0; i < p_num; i++) {
			x_old[i] = particles.X(i);
		}

		Update_Forces(particles.XRef(), particles.VRef(), particles.FRef());
		Enforce_Boundary_Conditions(particles.VRef(), particles.FRef());

		Update_Newmark_K(particles.XRef(), particles.VRef(), particles.FRef(), dt*(real).5);

		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			VectorD v=particles.V(i);
			for(int j=0;j<d;j++)u[i*d+j]=v[j];}	////set initial guess to be the velocity from the last time step
		SparseSolver::Conjugate_Gradient(K,u,b);

		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			VectorD v;for(int j=0;j<d;j++)v[j]=u[i*d+j];
			v_half[i]=v;} ////\tilde{v_half}=v_n+half_dt*a(x_n, \tilde{v_half})
		
		#pragma omp parallel for
		for (int i = 0; i < p_num; i++) {
			particles.X(i) += v_half[i] * dt; ////x_{n+1}=x_n+dt*\tilde{v_half}
		}

		Update_Forces(x_old, particles.VRef(), particles.FRef());
		Enforce_Boundary_Conditions(particles.VRef(), particles.FRef());
		for (int i = 0; i < p_num; i++) {
			particles.V(i) += dt * (real).5 * particles.F(i)/particles.M(i); ////v_half = v_n + half_dt*a(x_n, v_n)
		}

		Update_Forces(particles.XRef(), particles.VRef(), particles.FRef());
		Enforce_Boundary_Conditions(particles.VRef(), particles.FRef());

		Update_Newmark_K(particles.XRef(), particles.VRef(), particles.FRef(), dt*(real).5);

		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			VectorD v=particles.V(i);
			for(int j=0;j<d;j++)u[i*d+j]=v[j];}	////set initial guess to be the v_half
		SparseSolver::Conjugate_Gradient(K,u,b);

		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			VectorD v;for(int j=0;j<d;j++)v[j]=u[i*d+j];
			particles.V(i)=v;} ////v_{n+1}=v_half+half_dt*a(x_{n+1},v_{n+1}) 
	}

	void Initialize_Implicit_K()
	{
		////this function is not parallelized
		int n=d*p_num;
		K.resize(n,n);u.resize(n);u.fill((real)0);b.resize(n);b.fill((real)0);
		Array<TripletT> elements;
		for(int s=0;s<s_num;s++){int i=springs[s][0];int j=springs[s][1];
			Add_Triplet_Helper(i,i,elements);
			Add_Triplet_Helper(i,j,elements);
			Add_Triplet_Helper(j,i,elements);
			Add_Triplet_Helper(j,j,elements);}
		K.setFromTriplets(elements.begin(),elements.end());
		K.makeCompressed();	
	}

	void Update_Implicit_K(const Array<VectorD>& X,const Array<VectorD>& V,const Array<VectorD>& F,const real dt)
	{
		////Clear K
		SparseFunc::Set_Value(K,(real)0);
		
		////Set K diagonal blocks and rhs
		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			for(int ii=0;ii<d;ii++){K.coeffRef(i*d+ii,i*d+ii)+=particles.M(i);}
			VectorD rhs=particles.M(i)*V[i]+dt*F[i];Set_Block(b,i,rhs);}

		////sequential update implicit K
		//////Set K blocks with spring and damping Jacobians
		//MatrixD Ks,Kd;real dt_sq=pow(dt,2);
		//for(int s=0;s<s_num;s++){int i=springs[s][0];int j=springs[s][1];
		//	Compute_Ks_Block(X,s,Ks);
		//	Ks*=-dt_sq;
		//	Add_Block_Helper(K,i,j,Ks);
		//	Compute_Kd_Block(X,s,Kd);
		//	Kd*=-dt;
		//	Add_Block_Helper(K,i,j,Kd);

		//	VectorD rhs_d_i=Kd*(V[i]-V[j]);
		//	if(!Is_Psi_D(i))Add_Block(b,i,rhs_d_i);
		//	if(!Is_Psi_D(j))Add_Block(b,j,-rhs_d_i);}

		////parallelized assembling K
		real dt_sq=pow(dt,2);
		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			MatrixD Ks,Kd;
			int s_n=particle_springs[i].size();
			for(int j=0;j<s_n;j++){
				int s=particle_springs[i][j];
				Compute_Ks_Block(X,s,Ks);
				Ks*=-dt_sq;
				Compute_Kd_Block(X,s,Kd);
				Kd*=-dt;

				int si=springs[s][0];int sj=springs[s][1];
				if(i==si){
					Add_Block_Helper_Row_I(K,si,sj,Ks);
					Add_Block_Helper_Row_I(K,si,sj,Kd);}
				else{	////i==sj
					Add_Block_Helper_Row_J(K,si,sj,Ks);
					Add_Block_Helper_Row_J(K,si,sj,Kd);}
				
				VectorD rhs_d_i=Kd*(V[si]-V[sj]);
				if(!Is_Psi_D(si)&&i==si)Add_Block(b,si,rhs_d_i);
				if(!Is_Psi_D(sj)&&i==sj)Add_Block(b,sj,-rhs_d_i);}}
	}

	void Update_Newmark_K(const Array<VectorD>& X,const Array<VectorD>& V,const Array<VectorD>& F,const real dt) {
		////Clear K
		SparseFunc::Set_Value(K,(real)0);
		
		////Set K diagonal blocks and rhs
		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			for(int ii=0;ii<d;ii++){K.coeffRef(i*d+ii,i*d+ii)+=particles.M(i);}
			VectorD rhs=particles.M(i)*V[i]+dt*F[i];Set_Block(b,i,rhs);}

		////parallelized assembling K
		#pragma omp parallel for
		for(int i=0;i<p_num;i++){
			MatrixD Kd;
			int s_n=particle_springs[i].size();
			for(int j=0;j<s_n;j++){
				int s=particle_springs[i][j];
				Compute_Kd_Block(X,s,Kd);
				Kd*=-dt;

				int si=springs[s][0];int sj=springs[s][1];
				if(i==si){
					Add_Block_Helper_Row_I(K,si,sj,Kd);}
				else{	////i==sj
					Add_Block_Helper_Row_J(K,si,sj,Kd);}
				
				VectorD rhs_d_i=Kd*(V[si]-V[sj]);
				if(!Is_Psi_D(si)&&i==si)Add_Block(b,si,rhs_d_i);
				if(!Is_Psi_D(sj)&&i==sj)Add_Block(b,sj,-rhs_d_i);}}
	}

	////Spring force derivative
	void Compute_Ks_Block(const Array<VectorD>& X,const int s,MatrixD& Ks)
	{
		int i=springs[s][0];int j=springs[s][1];
		VectorD x_ij=X[j]-X[i];
		real length=x_ij.norm();
		real length_0=rest_length[s];
		Ks=ks[s]*((length_0/length-(real)1)*MatrixD::Identity()-(length_0/pow(length,3))*x_ij*x_ij.transpose());
	}

	////Damping force derivative
	void Compute_Kd_Block(const Array<VectorD>& X,const int s,MatrixD& Kd)
	{
		int i=springs[s][0];int j=springs[s][1];
		VectorD n_ij=(X[j]-X[i]).normalized();
		Kd=-kd[s]*n_ij*n_ij.transpose();
	}

protected:
	void Add_Triplet_Helper(const int i,const int j,Array<TripletT>& elements)
	{
		for(int ii=0;ii<d;ii++)for(int jj=0;jj<d;jj++)elements.push_back(TripletT(i*d+ii,j*d+jj,(real)0));
	}

	void Add_Block_Helper(SparseMatrixT& K,const int i,const int j,const MatrixD& Ks)
	{
		SparseFunc::Add_Block<d,MatrixD>(K,i,i,Ks);
		SparseFunc::Add_Block<d,MatrixD>(K,j,j,Ks);
		if(!Is_Psi_D(i)&&!Is_Psi_D(j)){
			SparseFunc::Add_Block<d,MatrixD>(K,i,j,-Ks);
			SparseFunc::Add_Block<d,MatrixD>(K,j,i,-Ks);}
	}

	////parallellized add block helpers
	void Add_Block_Helper_Row_I(SparseMatrixT& K,const int i,const int j,const MatrixD& Ks)
	{
		SparseFunc::Add_Block<d,MatrixD>(K,i,i,Ks);
		if(!Is_Psi_D(i)&&!Is_Psi_D(j)){
			SparseFunc::Add_Block<d,MatrixD>(K,i,j,-Ks);}
	}

	void Add_Block_Helper_Row_J(SparseMatrixT& K,const int i,const int j,const MatrixD& Ks)
	{
		SparseFunc::Add_Block<d,MatrixD>(K,j,j,Ks);
		if(!Is_Psi_D(i)&&!Is_Psi_D(j)){
			SparseFunc::Add_Block<d,MatrixD>(K,j,i,-Ks);}
	}

	void Set_Block(VectorX& b,const int i,const VectorD& bi)
	{for(int ii=0;ii<d;ii++)b[i*d+ii]=bi[ii];}

	void Add_Block(VectorX& b,const int i,const VectorD& bi)
	{for(int ii=0;ii<d;ii++)b[i*d+ii]+=bi[ii];}
};

#endif

/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2012-09-13

  Copyright (C) 2012 Universite de Strasbourg

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file maxwell_DG.hpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \author Philippe Helluy <helluy@math.unistra.fr>
   \author Thomas Strub <strub@math.unistra.fr>
   \author Pierre GERHARD <pierre.gerhard@gmail.com>
   \date 2014-11-01
 */

#ifndef __Maxwell_DG_H
#define __Maxwell_DG_H 1

  
// #include <feel/options.hpp>
// #include <feel/feelalg/backend.hpp>
// #include <feel/feeldiscr/functionspace.hpp>
// #include <feel/feelvf/vf.hpp>
// #include <feel/feeldiscr/mesh.hpp>
// #include <feel/feelfilters/geo.hpp>
// #include <feel/feelfilters/creategmshmesh.hpp>
// #include <feel/feelfilters/exporter.hpp>
// #include <feel/feeldiscr/projector.hpp>

#include <feel/feel.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace method
{
    enum class MethodType
    {
        EULER = 0,
        RK4 = 1
    };
}

template<int Dim,int Order_poly,int Order_geo>
class Maxwell_DG 
:
public Simget
{
    typedef Simget super;
public:
    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef Simplex<Dim,Order_geo> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
     /*
     * Discontinuous space
     */
#if defined( USE_LEGENDRE )
	/*TODO DEFINE LEGENDRE PROPRELY ACTUALLY DO NOT WORK
    typedef bases<Legendre<2,1,Order_poly, double> > d_basis_type; DO NOT WORK*/
#else
   typedef bases<Lagrange<Order_poly,Scalar,Discontinuous,PointSetEquiSpaced> > d_basis_type;
#endif
    typedef FunctionSpace<mesh_type, d_basis_type> d_space_type;
    typedef boost::shared_ptr<d_space_type> d_space_ptrtype;
    typedef typename d_space_type::element_type d_element_type;
       
    
    // typedef bases<d_basis_type,d_basis_type,d_basis_type,d_basis_type,d_basis_type,d_basis_type> w_basis_type;
    // typedef FunctionSpace<mesh_type, w_basis_type> w_space_type;
    // typedef boost::shared_ptr<w_space_type> w_space_ptrtype;
    // typedef typename w_space_type::element_type w_element_type;
    // typedef typename w_element_type:: sub_element<0>::type w_element_0_type;
    // typedef typename w_element_type:: sub_element<1>::type w_element_1_type;
    // typedef typename w_element_type:: sub_element<2>::type w_element_2_type;
    // typedef typename w_element_type:: sub_element<3>::type w_element_3_type;
    // typedef typename w_element_type:: sub_element<4>::type w_element_4_type;
    // typedef typename w_element_type:: sub_element<5>::type w_element_5_type;
  
    /*
     * Continuous space
     */
    typedef bases<Lagrange<Order_poly,Scalar> > c_basis_type;
    typedef FunctionSpace<mesh_type, c_basis_type> c_space_type;
    typedef boost::shared_ptr<c_space_type> c_space_ptrtype;
    typedef typename c_space_type::element_type c_element_type;
    
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    /*
     * Constructor
     */
    Maxwell_DG()
    :
    super(),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize(this->vm()["hsize"].template as<double>() ) ,
    CFL(this->vm()["CFL"].template as<double>() ) ,
	Tfinal(this->vm()["Tfinal"].template as<double>() ),
    geo_filename(this->vm()["gmsh.filename"].template as<std::string>() ),
    shape( this->vm()["shape"].template as<std::string>() ),
    t_m( (method::MethodType)this->vm()["t_m"].template as<int>() ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
    timers()
    {}

    FEELPP_DONT_INLINE
    void init();
    FEELPP_DONT_INLINE
    void run();

private:
    backend_ptrtype M_backend;
    double meshSize;
    double dt ;
    double Tfinal;
	double CFL;
    std::string shape;
    const std::string geo_filename;
    method::MethodType t_m;
    double time;
    
    d_space_ptrtype Xh;
    c_space_ptrtype Ch;

    sparse_matrix_ptrtype D;
    vector_ptrtype RHS_Ex, RHS_Ey, RHS_Bz;
        
    std::map<std::string, std::pair<boost::timer, double> > timers;

    mesh_ptrtype mesh;
    boost::shared_ptr<export_type> exporter;

private : 
	

   void assemble_LHS() ;
    
   void assemble_RHS(  
                                    double time,
                                    d_element_type& Exn,
                                    d_element_type& Eyn,
                                    d_element_type& Bzn
                                    );
                                    
    void solve(d_element_type& dtExn,
                d_element_type& dtEyn,
                d_element_type& dtBzn
                );
                
    void time_integrator(double dt, 
                                d_element_type& Exn,
                                d_element_type& Eyn,
                                d_element_type& Bzn
                                );
                                
                                
                                
    void exportResults(double time, 
                                d_element_type& Exn,
                                d_element_type& Eyn,
                                d_element_type& Bzn
                              );

};

template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::init()
{
    if ( !this->vm().count( "nochdir" ) == false ) //change to true
        Environment::changeRepository( boost::format( "examples/maxwell/%1%/%2%-%3%/P%4%-P%5%/h_%6%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order_poly
                                       % Order_geo
                                       % meshSize);
    
    
    mesh = createGMSHMesh( _mesh=new mesh_type, 
                           _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                           _desc=geo( _filename=geo_filename, _h=meshSize ));
    
    
    export_ptrtype exporter( export_type::New( this->vm(),
                                     ( boost::format( "%1%-%2%" )
                                       % this->about().appName()
                                       % shape ).str() ) );
        
    
    Xh= d_space_type::New( _mesh=mesh, _worldscomm=std::vector<WorldComm>( 1,mesh->worldComm() ),
                                          _extended_doftable=std::vector<bool>( 1,true ) ); //buildExtendedDofTable=true
    

    
    Ch = c_space_type::New( _mesh=mesh ); 
    D = M_backend->newMatrix( _test=Xh, _trial=Xh );
    RHS_Ex = M_backend->newVector( Xh );
    RHS_Ey = M_backend->newVector( Xh );
    RHS_Bz = M_backend->newVector( Xh );
}


template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::run()
{
    this->init();
  
    d_element_type Exn = Xh->element();
    d_element_type Eyn = Xh->element();
    d_element_type Bzn = Xh->element();
   
   d_element_type dtExn = Xh->element();
    d_element_type dtEyn = Xh->element();
    d_element_type dtBzn = Xh->element();
		
    if( Dim == 2)
    {
    
   
	double pi=M_PI;
    
    /*Declare here the exact solution, used for initialising fields
	Duplicated code in assemble_RHS, should be refactorised
	*/
	
	auto Ex_ex_expr = cst(0.);
    auto Ey_ex_expr = cos( 6*pi*(Px() - cst_ref(time)));
    auto Bz_ex_expr = cos( 6*pi*(Px() - cst_ref(time)));
	
	 auto L2ProjDisc = opProjection( _domainSpace=Xh, _imageSpace=Xh, _type=L2 );
	 int counter=0;
	 

    dt= CFL*meshSize/( 2*Order_poly+1 );
    time = 0;
  

   /*Init fields*/
     Exn = L2ProjDisc->project( Ex_ex_expr );
     Eyn = L2ProjDisc->project( Ey_ex_expr );
     Bzn = L2ProjDisc->project( Bz_ex_expr );
	 this->exportResults(time,Exn, Eyn, Bzn);

	this-> assemble_LHS();
	time = time +dt;  
	while (time <= Tfinal){
		counter++;
		this->assemble_RHS(time, Exn, Eyn, Bzn);
		this->solve(dtExn,dtEyn,dtBzn);
		/*Explicit Euler Stage*/
		Exn.add(dt, dtExn);
		Eyn.add(dt, dtEyn);
		Bzn.add(dt, dtBzn);
		time=time+dt;
		/*RK2 stage
		Exn.add( dt/2.0, dtExn);
		Eyn.add( dt/2.0, dtEyn);
		Bzn.add( dt/2.0, dtBzn);
		time=time + dt/2.0;
		this->assemble_RHS(time, Exn, Eyn, Bzn);
		this->solve(dtExn,dtEyn,dtBzn);
		Exn.add( dt, dtExn);
		Eyn.add( dt, dtEyn);
		Bzn.add( dt, dtBzn);
		time = time + dt/2.0;
		*/
  
		if ( counter % 3 == 0 ){
		this->exportResults(time,Exn, Eyn, Bzn);
		}   
		auto error = vf::project( _space=Xh, _expr=idv(Exn)-Ex_ex_expr);
		double maxerror = error.linftyNorm();
		double errorL2 = integrate( elements(mesh), (idv( Eyn )- Ey_ex_expr )*( idv( Eyn )-Ey_ex_expr ) ).evaluate()( 0, 0 );
		if ( Environment::worldComm().isMasterRank() ){
			/*Compute error L2 Linf*/
			 std::cout << "time :"<< time << "s - ||u_error||_2 = " << math::sqrt( errorL2 ) << "\n";
			 std::cout << "max error at time " << time  << " s :" << std::setprecision(16) << maxerror << "\n";
			 }
    }
	} 
}
        
template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::solve(d_element_type& dtExn,
																		  d_element_type& dtEyn,
																		  d_element_type& dtBzn
																		  )
{
        timers["solve"].first.restart();
        LOG(INFO) << "[solve] :: Start solving \n";
 if ( Dim == 2 )
{  
        M_backend->solve(_matrix=D,_solution=dtExn, _rhs=RHS_Ex );
        M_backend->solve(_matrix=D,_solution=dtEyn, _rhs=RHS_Ey );
        M_backend->solve(_matrix=D,_solution=dtBzn, _rhs=RHS_Bz );
 
    if ( Environment::worldComm().isMasterRank() )
    {
          std::cout << "[solve] :: System solved in "<< timers["solve"].first.elapsed() << "s\n";
    }
}
}

template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::assemble_LHS()
{
    timers["LHS"].first.restart();
    auto v=Xh->element();
    auto u=Xh->element();
   
   LOG(INFO) <<  "[assemble LHS] :: Assemblage  Start\n";
   
    auto a = form2( _test=Xh, _trial=Xh, _matrix=D);
    a = integrate( elements( mesh ), idt( u )*id( v ) );
    if ( Environment::worldComm().isMasterRank() )
    {
         std::cout<< "[assemble LHS] :: Assemblage Done in " << timers["LHS"].first.elapsed() << "s\n";
    }
	
	/*Matlab output, crash if matrix size is too big
	D->printMatlab("mass.m");
	*/
}


/*Compute the RHS*/


template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::assemble_RHS(                double time,
                                                                                                        d_element_type& Exn,
                                                                                                        d_element_type& Eyn,
                                                                                                        d_element_type& Bzn
                                                                                                      )
{    

	LOG(INFO) << "[assemble_RHS] :: Starting \n";
    timers["RHS"].first.restart();

 if ( Dim == 2 )
    {
	 /*
	 Define here the exact solution used in boundary condition
	 Should be refactorized 'coz of duplicated code in initialisation
	 */
	 
	double pi=M_PI;
	auto Ex_ex_expr = cst(0.);
    auto Ey_ex_expr = cos( 6*pi*(Px() - cst_ref(time)));
    auto Bz_ex_expr = cos( 6*pi*(Px() - cst_ref(time)));
	
	/*
	Another exact solution
	Bias 2D plane wave parameter theta
	
	double k=2*pi; //=0;
    double theta=0.2;
    double vu=cos( theta );
    double vv=sin( theta );
    auto c=cos( k * ( vu * Px() + vv * Py() - cst_ref( time ) ) + cst( pi/2.0 ) );
    auto Ex_ex_expr = -vv*c;
    auto Ey_ex_expr = vu*c;
    auto Bz_ex_expr = c;
	 
	*/
	/*
	Another exact solution
	2D n,m mode Cavity exact solution
    double c=1.;
	double omega=1;
	double n=1,m=1;
    auto Ex_ex_expr = -cst(c*c*( n*pi*omega ))*cos( m*pi*Px() )*sin( n*pi*Py() )*sin( omega*cst(time) );
	auto Ey_ex_expr =  cst(c*c*( m*pi/omega ))*sin( m*pi*Px() )*cos( n*pi*Py() )*sin( omega*cst(time) );
	auto Bz_ex_expr =  cos( m*pi*Px() )*cos( n*pi*Py() )*cos( omega*cst(time) );
	*/
				 
	auto lf_Ex=form1( _test=Xh, _vector=RHS_Ex,_init=true );
    auto lf_Ey=form1( _test=Xh, _vector=RHS_Ey,_init=true );
    auto lf_Bz=form1( _test=Xh, _vector=RHS_Bz,_init=true );
	
	/*Test function*/
	auto v = Xh->element();

	/*Stiffness part*/
	lf_Ex=integrate(_range=elements(mesh),_expr= id(v)*dyv(Bzn));
    lf_Ey=integrate(_range=elements(mesh),_expr=-id(v)*dxv(Bzn));
    lf_Bz=integrate(_range=elements(mesh),_expr=-id(v)*dxv(Eyn)+id(v)*dyv(Exn));

    /*	Define Right and Left value */
	auto wR = vec( rightfacev( idv( Exn ) ),rightfacev( idv( Eyn ) ),rightfacev( idv( Bzn ) ) );
    auto wL = vec( leftfacev( idv( Exn ) ),leftfacev( idv( Eyn ) ),leftfacev( idv( Bzn ) ) );
	
    /* Upwind Flux : 	Matrices Aini(+) and Aini(-)  for upwind flux 	*/
	auto Anp_1 = vec( 0.5*Ny()*Ny(), -0.5*Nx()*Ny(), -0.5*Ny() );
    auto Anp_2 = vec( -0.5*Nx()*Ny(), 0.5*Nx()*Nx(), 0.5*Nx() );
    auto Anp_3 = vec( -0.5*Ny(), 0.5*Nx(), cst(0.5) );
    auto Anm_1 = vec( -0.5*Ny()*Ny(), 0.5*Nx()*Ny(),-0.5*Ny() );
    auto Anm_2 = vec( 0.5*Nx()*Ny(), -0.5*Nx()*Nx(), 0.5*Nx() );
    auto Anm_3 = vec( -0.5*Ny(), 0.5*Nx(), cst(-0.5) );
		
	
	lf_Ex +=integrate(_range=internalfaces(mesh),_expr=( trans(Anm_1)*(wL-wR))*leftface(id(v) )
																			   + ( trans(Anp_1)*(wL-wR))*rightface(id(v) ) );
	
	lf_Ey +=integrate(_range=internalfaces(mesh),_expr=( trans(Anm_2)*(wL-wR))*leftface(id(v) )
																			   + ( trans(Anp_2)*(wL-wR))*rightface(id(v) ) );
	
	
	lf_Bz +=integrate(_range=internalfaces(mesh),_expr=( trans(Anm_3)*(wL-wR))*leftface(id(v) )
																			   + ( trans(Anp_3)*(wL-wR))*rightface(id(v) ) );
	
	
	
	/* Boundary flux : Imposed weakly using ghost state wL*/
	auto wex=vec(Ex_ex_expr,Ey_ex_expr,Bz_ex_expr);
    auto wMetal = vec( -leftfacev( idv( Exn ) ),-leftfacev( idv( Eyn ) ),leftfacev( idv( Bzn ) ) );
    // auto wSivlerM= vec( leftfacev( idv( Exn ) ),-leftfacev( idv( Eyn ) ),-leftfacev( idv( Bzn ) ) );
	
	lf_Ex += integrate(markedfaces( mesh, "Dirichlet" ),trans(Anm_1)*(wL-wex)*id(v) );
	lf_Ey += integrate(markedfaces( mesh, "Dirichlet" ),trans(Anm_2)*(wL-wex)*id(v) );
	lf_Bz += integrate(markedfaces( mesh, "Dirichlet" ),trans(Anm_3)*(wL-wex)*id(v) );
	
	lf_Ex += integrate(markedfaces( mesh, "Metal" ),trans(Anm_1)*(wL-wMetal)*id(v) );
	lf_Ey += integrate(markedfaces( mesh, "Metal" ),trans(Anm_2)*(wL-wMetal)*id(v) );
	lf_Bz += integrate(markedfaces( mesh, "Metal" ),trans(Anm_3)*(wL-wMetal)*id(v) );

	
	/* Centered Flux	
	auto Aini_1 = vec(cst(0.), cst(0.),-Ny() );
	auto Aini_2 = vec(cst(0.), cst(0.), Nx() );
	auto Aini_3 = vec(-Ny() , Nx() ,cst(0.) );
	
	lf_Ex+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_1)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_Ey+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_2)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_Bz+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_3)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	*/
	
	/* Rusanov fux :  Here Lmax denotes the max velocity of waves, physicaly the speed. Matrices Aini are re-used.
	lf_Ex+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_1)*(wL*(leftface(id(v)) - wR*rightface(id(v))))) - lmax*0.5*(trans(Aini_1)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_Ey+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_2)*(wL*(leftface(id(v)) - wR*rightface(id(v))))) - lmax*0.5*(trans(Aini_2)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_Bz+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_3)*(wL*(leftface(id(v)) - wR*rightface(id(v))))) - lmax*0.5*(trans(Aini_3)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	*/
	}
   
 
	/*
	This part is the begining of the 3D extension
	
	/_!_\ Actually vec() function is NOT overloaded for 6th component !! DO NOT WORK
		
	if ( Dim == 3 )
    {
  
   TODO : Expres the stiffness part 
   
   
   	auto wR  = vec( rightfacev( idv( Exn ) ), rightfacev( idv( Eyn ) ), rightfacev( idv( Ezn ) ), rightfacev( idv( Bxn ) ), rightfacev( idv( Byn ) ), rightfacev( idv( Bzn ) ) );
    auto wL  = vec( leftfacev( idv( Exn ) ), leftfacev( idv( Eyn ) ), leftfacev( idv( Ezn ) ), leftfacev( idv( Bxn ) ), leftfacev( idv( Byn ) ), leftfacev( idv( Bzn ) ) );
      
   Updwind Flux Matrices Aini(+) and Aini(-)
 
	auto Anp_1_3D = vec(  0.5*(Ny()*Ny()+Nz()*Nz()) , -0.5*Nx()*Ny(), -0.5*Ny()*Nz(), 0, 0.5*Nz(), -0.5*Nz() );
    auto Anp_2_3D = vec( -0.5*Nx()*Ny(), 0.5*(Nx()*Nx()+ Nz()*Nz()), -0.5*Ny()*Nz(), -0.5*Nz(),0,0.5*Nx() );
	auto Anp_3_3D = vec(-0.5*Nx()*Nz(), -0.5*Ny()*Nz(), 0.5*(Ny()*Ny()+Nx()*Nx() ),0.5*Ny(),-0.5*Nx(),0);
	auto Anp_4_3D = vec(0, -0.5*Nz(), 0.5*Ny(), 0.5*(Ny()*Ny()+Nz()*Nz() ),-0.5*Nx()*Ny(),-0.5*Nx()*Nz());
	auto Anp_5_3D = vec(0.5*Nz(), 0, -0.5*Nx(),-0.5*Nx()*Ny(),0.5*(Nx()*Nx()+Nz()*Nz() ),-0.5*Ny()*Nz());
	auto Anp_6_3D = vec(-0.5*Ny(), -0.5*Nx(), 0,-0.5*Nx()*Nz(),-0.5*Ny()*Nz(),-0.5*(Ny()*Ny()+Nx()*Nx() ));
      
	auto Anm_1_3D = vec( -0.5*(Ny()*Ny()+Nz()*Nz()) , 0.5*Nx()*Ny(), 0.5*Ny()*Nz(), 0, 0.5*Nz(), -0.5*Nz() );
    auto Anm_2_3D = vec( 0.5*Nx()*Ny(), -0.5*(Nx()*Nx()+ Nz()*Nz()), 0.5*Ny()*Nz(), 0.5*Nz(),0,0.5*Nx() );
	auto Anm_3_3D = vec(0.5*Nx()*Nz(), 0.5*Ny()*Nz(), -0.5*(Ny()*Ny()+Nx()*Nx() ),0.5*Ny(),-0.5*Nx(),0);
	auto Anm_4_3D = vec(0, -0.5*Nz(), 0.5*Ny(), -0.5*(Ny()*Ny()+Nz()*Nz() ),0.5*Nx()*Ny(),0.5*Nx()*Nz());
	auto Anm_5_3D = vec(0.5*Nz(), 0, -0.5*Nx(),0.5*Nx()*Ny(), -0.5*(Nx()*Nx()+Nz()*Nz() ),0.5*Ny()*Nz());
	auto Anm_6_3D = vec(-0.5*Ny(), -0.5*Nx(), 0,0.5*Nx()*Nz(),0.5*Ny()*Nz(),0.5*(Ny()*Ny()+Nx()*Nx() )); 
	
		
	lf_Ex +=integrate(_range=internalfaces(mesh),_expr=( trans(Anm_1)*(wL-wR))*leftface(id(v) )
																			   + ( trans(Anp_1)*(wL-wR))*rightface(id(v) ) );
	
	lf_Ey +=integrate(_range=internalfaces(mesh),_expr=( trans(Anm_2)*(wL-wR))*leftface(id(v) )
																			   + ( trans(Anp_2)*(wL-wR))*rightface(id(v) ) );
	
	
	lf_Ez +=integrate(_range=internalfaces(mesh),_expr=( trans(Anm_3)*(wL-wR))*leftface(id(v) )
																			   + ( trans(Anp_3)*(wL-wR))*rightface(id(v) ) );
	
	
	
	lf_Bx +=integrate(_range=internalfaces(mesh),_expr=( trans(Anm_4)*(wL-wR))*leftface(id(v) )
																			   + ( trans(Anp_4)*(wL-wR))*rightface(id(v) ) );
	
	lf_By +=integrate(_range=internalfaces(mesh),_expr=( trans(Anm_5)*(wL-wR))*leftface(id(v) )
																			   + ( trans(Anp_5)*(wL-wR))*rightface(id(v) ) );
	
	
	lf_Bz +=integrate(_range=internalfaces(mesh),_expr=( trans(Anm_6)*(wL-wR))*leftface(id(v) )
																			   + ( trans(Anp_6)*(wL-wR))*rightface(id(v) ) );
	
	
		
	Centered Flux
	
	auto Aini_1 = vec(cst(0.), cst(0.), cst(0.), cst(0.), Nz(), -Ny() );
	auto Aini_2 = vec(cst(0.), cst(0.), cst(0.), -Nz(), cst(0.), Nx() );
	auto Aini_3 = vec(cst(0.), cst(0.), cst(0.), Ny(), -Nx(), cst(0.) );
	auto Aini_4 = vec(cst(0.), -Nz(), Ny(), cst(0.), cst(0.), cst(0.) );
	auto Aini_5 = vec(Nz(), cst(0.), -Nx(), cst(0.), cst(0.), cst(0.) );
	auto Aini_6 = vec(-Ny(), Nx(), cst(0.), cst(0.), cst(0.), cst(0.) );
	
	
	lf_Ex+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_1)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_Ey+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_2)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_Ez+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_3)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_Bx+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_4)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_By+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_5)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_Bz+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_6)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	
	
	Rusanov fux :  Here Lmax denotes the max velocity of waves, physicaly the speed of light in medium. Matrices Aini are re-used.
	lf_Ex+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_1)*(wL*(leftface(id(v)) - wR*rightface(id(v))))) - lmax*0.5*(trans(Aini_1)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_Ey+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_2)*(wL*(leftface(id(v)) - wR*rightface(id(v))))) - lmax*0.5*(trans(Aini_2)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_Ez+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_3)*(wL*(leftface(id(v)) - wR*rightface(id(v))))) - lmax*0.5*(trans(Aini_3)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_Bx+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_4)*(wL*(leftface(id(v)) - wR*rightface(id(v))))) - lmax*0.5*(trans(Aini_4)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_By+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_5)*(wL*(leftface(id(v)) - wR*rightface(id(v))))) - lmax*0.5*(trans(Aini_5)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
	lf_Bz+=integrate( _range=internalfaces(mesh),_expr=0.5*(trans(Aini_6)*(wL*(leftface(id(v)) - wR*rightface(id(v))))) - lmax*0.5*(trans(Aini_6)*(wL*(leftface(id(v)) - wR*rightface(id(v))))));
			
	
	
	Boundary flux : Imposed weakly using ghost state wL
	auto wex = vec(Ex_ex_expr, Ey_ex_expr, Ez_ex_expr, Bx_ex_expr, By_ex_expr, Bz_ex_expr);
   
   TO DO : Express the wMetal and Silver Muller
    auto wMetal = ;
	
	lf_Ex += integrate(markedfaces( mesh, "Dirichlet" ),trans(Anm_1)*(wL-wex)*id(v) );
	lf_Ey += integrate(markedfaces( mesh, "Dirichlet" ),trans(Anm_2)*(wL-wex)*id(v) );
	lf_Ez += integrate(markedfaces( mesh, "Dirichlet" ),trans(Anm_3)*(wL-wex)*id(v) );
	lf_Bx += integrate(markedfaces( mesh, "Dirichlet" ),trans(Anm_4)*(wL-wex)*id(v) );
	lf_By += integrate(markedfaces( mesh, "Dirichlet" ),trans(Anm_5)*(wL-wex)*id(v) );
	lf_Bz += integrate(markedfaces( mesh, "Dirichlet" ),trans(Anm_6)*(wL-wex)*id(v) );
	
	lf_Ex += integrate(markedfaces( mesh, "Metal" ),trans(Anm_1)*(wL-wMetal)*id(v) );
	lf_Ey += integrate(markedfaces( mesh, "Metal" ),trans(Anm_2)*(wL-wMetal)*id(v) );
	lf_Ez += integrate(markedfaces( mesh, "Metal" ),trans(Anm_3)*(wL-wMetal)*id(v) );
	lf_Bx += integrate(markedfaces( mesh, "Metal" ),trans(Anm_4)*(wL-wMetal)*id(v) );
	lf_By += integrate(markedfaces( mesh, "Metal" ),trans(Anm_5)*(wL-wMetal)*id(v) );
	lf_Bx += integrate(markedfaces( mesh, "Metal" ),trans(Anm_6)*(wL-wMetal)*id(v) );
   }
*/
    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "[assemble_RHS] :: Completed \n";
    }
} 


/*
Time integrator function, Actually not used
*/
template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::time_integrator(double dt, 
                                                                    d_element_type& Exn,
                                                                    d_element_type& Eyn,
                                                                    d_element_type& Bzn
                                                                  )
                                                                  
{

    d_element_type Exk = Xh->element();
    d_element_type Eyk = Xh->element();
    d_element_type Bzk = Xh->element();

    assemble_RHS(Exn, Eyn,  Bzn);
    solve(Exn, Eyn, Bzn);

    Exn.add(dt, Exn);
    Eyn.add(dt, Eyn);
    Bzn.add(dt, Bzn);
    
     
    /*3D RK4 Stage WARNING UNTESTED
    Exk=Exn;
    Eyk=Eyn;
    Ezk=Ezn;
    Bxk=Bxn;
    Byk=Byn;
    Bzk=Bzn;
     
    Exn.add( dt/6.0, Exn );
    Eyn.add( dt/6.0, Eyn ); 
    Ezn.add( dt/6.0, Ezn ); 
    Bxn.add( dt/6.0, Bxn );
    Byn.add( dt/6.0, Byn ); 
    Bzn.add( dt/6.0, Bzn ); 

    Exk.add( dt/2.0, Exk );
    Eyk.add( dt/2.0, Eyk ); 
    Ezk.add( dt/2.0, Ezk );  
    Bxk.add( dt/2.0, Bxk );
    Byk.add( dt/2.0, Byk ); 
    Bzk.add( dt/2.0, Bzk ); 

    assemble_RHS(Exk,  Eyk, Ezk, Bxk, Byk, Bzk);
    solve(Exk,  Eyk, Ezk, Bxk, Byk, Bzk);
  
    Exn.add( dt/3.0, Exk );
    Eyn.add( dt/3.0, Eyk ); 
    Ezn.add( dt/3.0, Ezk ); 
    Bxn.add( dt/3.0, Bxk );
    Byn.add( dt/3.0, Byk ); 
    Bzn.add( dt/3.0, Bzk ); 

    Exk.add( dt/2.0, Exk );
    Eyk.add( dt/2.0, Eyk ); 
    Ezk.add( dt/2.0, Ezk );  
    Bxk.add( dt/2.0, Bxk );
    Byk.add( dt/2.0, Byk ); 
    Bzk.add( dt/2.0, Bzk ); 

    assemble_RHS(Exk,  Eyk, Ezk, Bxk, Byk, Bzk);
    solve(Exk,  Eyk, Ezk, Bxk, Byk, Bzk);  


    Exn.add( dt/3.0, Exk );
    Eyn.add( dt/3.0, Eyk ); 
    Ezn.add( dt/3.0, Ezk ); 
    Bxn.add( dt/3.0, Bxk );
    Byn.add( dt/3.0, Byk ); 
    Bzn.add( dt/3.0, Bzk );   

    Exk.add( dt, Exk );
    Eyk.add( dt, Eyk ); 
    Ezk.add( dt, Ezk );  
    Bxk.add( dt, Bxk );
    Byk.add( dt, Byk ); 
    Bzk.add( dt, Bzk );  

    assemble_RHS(Exk,  Eyk, Ezk, Bxk, Byk, Bzk);
    solve(Exk,  Eyk, Ezk, Bxk, Byk, Bzk);    


    Exn.add( dt/6.0, Exk );
    Eyn.add( dt/6.0, Eyk ); 
    Ezn.add( dt/6.0, Ezk ); 
    Bxn.add( dt/6.0, Bxk );
    Byn.add( dt/6.0, Byk ); 
    Bzn.add( dt/6.0, Bzk ); 
*/

}


/*Exporter function*/

template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::exportResults(double time, 
                                                                    d_element_type& Exn,
                                                                    d_element_type& Eyn,
                                                                    d_element_type& Bzn
                                                                  )
{
    auto L2ProjCon = opProjection( _domainSpace=Ch, _imageSpace=Ch, _type=L2 );

    if ( exporter->doExport()){
        exporter->step( time )->setMesh( mesh );
       if(Dim==2){
        exporter->step( time )->add( "Ex", L2ProjCon->project( idv( Exn ) ) );
        exporter->step( time )->add( "Ey", L2ProjCon->project( idv( Eyn ) ) );
        exporter->step( time )->add( "Bz", L2ProjCon->project( idv( Bzn ) ) );
       }

        exporter->save();
        LOG(INFO) << "exportResults done\n";
    }
}
#endif /* __Maxwell_DG_H */
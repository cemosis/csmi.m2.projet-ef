
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
    typedef bases<Legendre<Order_poly,Scalar> > d_basis_type;
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
    typedef Exporter<mesh_type,Order_geo> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    /*
     * Constructor
     */
    Maxwell_DG()
    :
    super(),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize(this->vm()["hsize"].template as<double>() ) ,
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
    std::string shape;
    const std::string geo_filename;
    method::MethodType t_m;
    double time;
    d_space_ptrtype Xh;
    c_space_ptrtype Ch;
    mesh_ptrtype mesh;
    sparse_matrix_ptrtype D;
    std::map<std::string, std::pair<boost::timer, double> > timers;
    export_ptrtype exporter;
    vector_ptrtype RHS_Ex, RHS_Ey, RHS_Ez, RHS_Bx, RHS_By, RHS_Bz;

private : 
    void assemble_LHS() ;
    void assemble_RHS(  
                                    d_element_type& Ex_ex,
                                    d_element_type& Ey_ex,
                                    d_element_type& Ez_ex,
                                    d_element_type& Bx_ex,
                                    d_element_type& By_ex,
                                    d_element_type& Bz_ex,
                                    d_element_type& Exn,
                                    d_element_type& Eyn,
                                    d_element_type& Ezn,
                                    d_element_type& Bxn,
                                    d_element_type& Byn,
                                    d_element_type& Bzn
                                    );
                                    
    void solve(d_element_type& dtExn,
                d_element_type& dtEyn,
                d_element_type& dtEzn,
                d_element_type& dtBxn,
                d_element_type& dtByn,
                d_element_type& dtBzn
                );
                
    void time_integrator(double dt, 
                                d_element_type& Exn,
                                d_element_type& Eyn,
                                d_element_type& Ezn,
                                d_element_type& Bxn,
                                d_element_type& Byn,
                                d_element_type& Bzn
                                );
                                
                                
                                
    void exportResults(double time, 
                                d_element_type& Exn,
                                d_element_type& Eyn,
                                d_element_type& Ezn,
                                d_element_type& Bxn,
                                d_element_type& Byn,
                                d_element_type& Bzn
                              );

};

template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::init()
{
    if ( !this->vm().count( "nochdir" ) == false )
        Environment::changeRepository( boost::format( "examples/maxwell/%1%/%2%-%3%/P%4%-P%5%/h_%6%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order_poly
                                       % Order_geo
                                       % meshSize);
    
    mesh = createGMSHMesh( _mesh=new Mesh<Simplex<Dim,Order_geo>>, 
                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                            _desc=geo( _filename=geo_filename, _h=meshSize ) );
    
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
    RHS_Ez = M_backend->newVector( Xh );
    RHS_Bx = M_backend->newVector( Xh );
    RHS_By = M_backend->newVector( Xh );
    RHS_Bz = M_backend->newVector( Xh );
}


template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::run()
{
    this->init();
  
 
    d_element_type Exn = Xh->element();
    d_element_type Eyn = Xh->element();
    d_element_type Ezn = Xh->element();
    d_element_type Bxn = Xh->element();
    d_element_type Byn = Xh->element();
    d_element_type Bzn = Xh->element();
    
    d_element_type Ex_ex = Xh->element();
    d_element_type Ey_ex = Xh->element();
    d_element_type Ez_ex = Xh->element();
    d_element_type Bx_ex = Xh->element();
    d_element_type By_ex = Xh->element();
    d_element_type Bz_ex = Xh->element();
    
    
    
    
    auto dtExn =  Xh->element();
    if( Dim == 2)
    {
    double pi=M_PI;
    double k=2*pi; //=0;
    double theta=0.2;
    double vu=cos( theta );
    double vv=sin( theta );
    auto c=cos( k * ( vu * Px() + vv * Py() - cst_ref( time ) ) + cst( pi/2.0 ) );
    auto Ex_ex_expr = -vv*c;
    auto Ey_ex_expr = vu*c;
    auto Bz_ex_expr = c;
    auto L2ProjDisc = opProjection( _domainSpace=Xh, _imageSpace=Xh, _type=L2 );
        
    

    Ezn = L2ProjDisc->project( Ex_ex_expr );
    Bxn = L2ProjDisc->project( Ey_ex_expr );
    Byn = L2ProjDisc->project( Bz_ex_expr );

    dt= 0.1*meshSize/( 2*Order_poly+1 );

    this-> assemble_LHS();
    
    while (time <= Tfinal)
    {
    
       //this->time_integrator(dt, Exn, Eyn, Ezn, Bxn, Byn, Bzn );
        // exportResults(time, Exn, Eyn, Ezn, Bxn, Byn, Bzn );     
                  
    Ex_ex = L2ProjDisc->project( Ex_ex_expr );
    Ey_ex = L2ProjDisc->project( Ey_ex_expr );
    Bz_ex = L2ProjDisc->project( Bz_ex_expr );
     
    //this->exportResults(time, Exn, Eyn, Ezn, Bxn, Byn, Bzn );     
    this->assemble_RHS(Ex_ex,Ey_ex,Ez_ex,Bx_ex,By_ex,Bz_ex,Exn, Eyn, Ezn, Bxn, Byn, Bzn);
    this->solve(Exn, Eyn, Ezn, Bxn, Byn, Bzn);
    Exn.add(dt, Exn);
    Eyn.add(dt, Eyn);
    Bzn.add(dt, Bzn);
   
    if ( exporter->doExport()){
    LOG(INFO) << "exportResults done\n";
    this->exportResults(time,Exn, Eyn, Ezn, Bxn, Byn, Bzn);
    }
    
   
    auto error = vf::project( _space=Xh, _expr=idv(Exn)-Ex_ex_expr);
    double maxerror = error.linftyNorm();
 
        if ( Environment::worldComm().isMasterRank() ){
         std::cout << "max error at time " << time  << " s :" << std::setprecision(16) << maxerror << "\n";
         }
        
        time = time + dt;
        
        
    }


    // this->assemble_RHS(Exn,  Eyn, Ezn, Bxn, Byn, Bzn);
    // this->solve(Exn,  Eyn, Ezn, Bxn, Byn, Bzn);
  } 
}
        
template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::solve( d_element_type& dtExn,
                                                                                      d_element_type& dtEyn,
                                                                                      d_element_type& dtEzn,
                                                                                      d_element_type& dtBxn,
                                                                                      d_element_type& dtByn,
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
}


template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::assemble_RHS(                
                                                                                                        d_element_type& Ex_ex,
                                                                                                        d_element_type& Ey_ex,
                                                                                                        d_element_type& Ez_ex,
                                                                                                        d_element_type& Bx_ex,
                                                                                                        d_element_type& By_ex,
                                                                                                        d_element_type& Bz_ex,
                                                                                                        d_element_type& Exn,
                                                                                                        d_element_type& Eyn,
                                                                                                        d_element_type& Ezn,
                                                                                                        d_element_type& Bxn,
                                                                                                        d_element_type& Byn,
                                                                                                        d_element_type& Bzn
                                                                                                      )
{    
 if ( Dim == 2 )
    {
    auto w = vec( idv( Exn ),idv( Eyn ),idv( Bzn ) );
    auto wex= vec( idv( Ex_ex ),idv( Ey_ex ),idv( Bz_ex ) );
    auto wR = vec( rightfacev( idv( Exn ) ),rightfacev( idv( Eyn ) ),rightfacev( idv( Bzn ) ) );

    auto wL = vec( leftfacev( idv( Exn ) ),leftfacev( idv( Eyn ) ),leftfacev( idv( Bzn ) ) );
    //On fait rebondir le flux
    //SIGMA W_L
    auto wMetal = vec( -leftfacev( idv( Exn ) ),-leftfacev( idv( Eyn ) ),leftfacev( idv( Bzn ) ) );
    auto wSivlerM= vec( leftfacev( idv( Exn ) ),-leftfacev( idv( Eyn ) ),-leftfacev( idv( Bzn ) ) );
      // FLux matrix GHOST CELL
    auto Anp_1 = vec( 0.5*Ny()*Ny(), -0.5*Nx()*Ny(), -0.5*Ny() );
    auto Anp_2 = vec( -0.5*Nx()*Ny(), 0.5*Nx()*Nx(), 0.5*Nx() );
    auto Anp_3 = vec( -0.5*Ny(), 0.5*Nx(), cst(0.5) );
    auto Anm_1 = vec( -0.5*Ny()*Ny(), 0.5*Nx()*Ny(),-0.5*Ny() );
    auto Anm_2 = vec( 0.5*Nx()*Ny(), -0.5*Nx()*Nx(), 0.5*Nx() );
    auto Anm_3 = vec( -0.5*Ny(), 0.5*Nx(), cst(-0.5) );
    auto v = Xh->element();
         
    auto lf_Ex=form1( _test=Xh, _vector=RHS_Ex );
    auto lf_Ey=form1( _test=Xh, _vector=RHS_Ey );
    // auto lf_Ez=form1( _test=Xh, _vector=RHS_Ez );
    // auto lf_Bx=form1( _test=Xh, _vector=RHS_Bx );
    // auto lf_By=form1( _test=Xh, _vector=RHS_By );
    auto lf_Bz=form1( _test=Xh, _vector=RHS_Bz );
     
    LOG(INFO) << "[assemble_RHS] :: Starting \n";
    timers["RHS"].first.restart();
    
    
    
    //First Term Decentred flux should be ok
    lf_Ex = integrate(  _range=elements( mesh ), 
    _expr=+id( v )*dyv( Bzn ) );

    lf_Ey = integrate( _range=elements( mesh ),  
    _expr=+id( v )*dxv( Bzn ) ); //sign change - -> +

    lf_Bz = integrate(_range=elements( mesh ), 
    _expr=-id( v )*dxv( Eyn ) + id( v )*dyv( Exn ) );
    
    LOG(INFO) << "[assemble_RHS] :: TM MODE - First RHS Term done \n";
    
    //Second term Decentred Flux
    lf_Ex += integrate( _range=internalfaces( mesh ),
    _expr= ( trans( Anm_1 )*( wL-wR ) )*leftface( id( v ) )
    + ( trans( Anp_1  )*( wL-wR ) )*rightface( id( v ) ) );

    lf_Ey += integrate( _range=internalfaces( mesh ),
    _expr= ( trans( Anm_2 )*( wL-wR ) )*leftface( id( v ) )
    + ( trans( Anp_2 )*( wL-wR ) )*rightface( id( v ) ) );

    lf_Bz += integrate( _range=internalfaces( mesh ),
    _expr= ( trans( Anm_3 )*( wL-wR ) )*leftface( id( v ) ) 
    + ( trans( Anp_3 )*( wL-wR ) )*rightface( id( v ) ) );    

    LOG(INFO) << "[assemble_RHS] :: TM MODE - Second RHS Term done \n";
    
    lf_Ex += integrate(  _range=markedfaces( mesh, "Metal" ), 
    _expr=( trans( Anm_1 )*( wL-wMetal ) )*leftface( id( v ) ) );
    
    lf_Ey += integrate( _range=markedfaces( mesh, "Metal" ),   
    _expr=( trans( Anm_2 )*( wL-wMetal ) )*leftface( id( v ) ) );                                    
          
    lf_Bz += integrate( _range=markedfaces( mesh, "Metal" ), 
    _expr=( trans( Anm_3 )*( wL-wMetal ) )*leftface( id( v ) ) );

    LOG(INFO) << "[assemble_RHS] :: TM MODE - Metalic Boundary RHS Term done \n";
    
    
    //TCHECK IF GOOD.
    lf_Ex += integrate( _range=markedfaces( mesh, "Dirichlet" ), 
    _expr=-( trans( Anm_1 )*( wL-wex ) )*leftface( id( v ) ));

    lf_Ey += integrate(_range=markedfaces( mesh, "Dirichlet" ), 
    _expr=-( trans( Anm_2 )*( wL-wex ) )*leftface( id( v ) ) );

    lf_Bz += integrate(_range=markedfaces( mesh, "Dirichlet" ), 
    _expr=-( trans( Anm_3 )*( wL-wex ) )*leftface( id( v ) ) );                               
                                    
    LOG(INFO) << "[assemble_RHS] :: TM MODE - Dirichlet Boundary RHS Term done \n";
        
   }
    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "[assemble_RHS] :: Completed \n";
    }
} 


template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::time_integrator(double dt, 
                                                                    d_element_type& Exn,
                                                                    d_element_type& Eyn,
                                                                    d_element_type& Ezn,
                                                                    d_element_type& Bxn,
                                                                    d_element_type& Byn,
                                                                    d_element_type& Bzn
                                                                  )
                                                                  
{

    d_element_type Exk = Xh->element();
    d_element_type Eyk = Xh->element();
    d_element_type Ezk = Xh->element();
    d_element_type Bxk = Xh->element();
    d_element_type Byk = Xh->element();
    d_element_type Bzk = Xh->element();

    assemble_RHS(Exn, Eyn, Ezn, Bxn, Byn, Bzn);
    solve(Exn, Eyn, Ezn, Bxn, Byn, Bzn);

    Exn.add(dt, Exn);
    Eyn.add(dt, Eyn);
    Bzn.add(dt, Bzn);
    
    }
    
    
    
    
    /*
    //Ex,Ey ,etc... stockent le stage 1

    Exk=Exn;
    Eyk=Eyn;
    // Ezk=Ezn;
    // Bxk=Bxn;
    // Byk=Byn;
    Bzk=Bzn;
     
     //On actualise EX au fur et a mesure
    Exn.add( dt/6.0, Exn );
    Eyn.add( dt/6.0, Eyn ); 
    // Ezn.add( dt/6.0, Ezn ); 
    // Bxn.add( dt/6.0, Bxn );
    // Byn.add( dt/6.0, Byn ); 
    Bzn.add( dt/6.0, Bzn ); 

    Exk.add( dt/2.0, Exk );
    Eyk.add( dt/2.0, Eyk ); 
    // Ezk.add( dt/2.0, Ezk );  
    // Bxk.add( dt/2.0, Bxk );
    // Byk.add( dt/2.0, Byk ); 
    Bzk.add( dt/2.0, Bzk ); 

    assemble_RHS(Exk,  Eyk, Ezk, Bxk, Byk, Bzk);
    solve(Exk,  Eyk, Ezk, Bxk, Byk, Bzk);

    //Exk,Eyk ,etc... stockent le stage 2
    
    Exn.add( dt/3.0, Exk );
    Eyn.add( dt/3.0, Eyk ); 
    // Ezn.add( dt/3.0, Ezk ); 
    // Bxn.add( dt/3.0, Bxk );
    // Byn.add( dt/3.0, Byk ); 
    Bzn.add( dt/3.0, Bzk ); 

    Exk.add( dt/2.0, Exk );
    Eyk.add( dt/2.0, Eyk ); 
    // Ezk.add( dt/2.0, Ezk );  
    // Bxk.add( dt/2.0, Bxk );
    // Byk.add( dt/2.0, Byk ); 
    Bzk.add( dt/2.0, Bzk ); 

    assemble_RHS(Exk,  Eyk, Ezk, Bxk, Byk, Bzk);
    solve(Exk,  Eyk, Ezk, Bxk, Byk, Bzk);  
    //Exk,Eyk ,etc... stockent le stage 3

    Exn.add( dt/3.0, Exk );
    Eyn.add( dt/3.0, Eyk ); 
    // Ezn.add( dt/3.0, Ezk ); 
    // Bxn.add( dt/3.0, Bxk );
    // Byn.add( dt/3.0, Byk ); 
    Bzn.add( dt/3.0, Bzk );   

    Exk.add( dt, Exk );
    Eyk.add( dt, Eyk ); 
    // Ezk.add( dt, Ezk );  
    // Bxk.add( dt, Bxk );
    // Byk.add( dt, Byk ); 
    Bzk.add( dt, Bzk );  

    assemble_RHS(Exk,  Eyk, Ezk, Bxk, Byk, Bzk);
    solve(Exk,  Eyk, Ezk, Bxk, Byk, Bzk);    
    //Exk,Eyk ,etc... stockent le stage 4

    Exn.add( dt/6.0, Exk );
    Eyn.add( dt/6.0, Eyk ); 
    // Ezn.add( dt/6.0, Ezk ); 
    // Bxn.add( dt/6.0, Bxk );
    // Byn.add( dt/6.0, Byk ); 
    Bzn.add( dt/6.0, Bzk ); 

}
*/


template<int Dim,int Order_poly,int Order_geo>
void 
Maxwell_DG<Dim, Order_poly, Order_geo>::exportResults(double time, 
                                                                    d_element_type& Exn,
                                                                    d_element_type& Eyn,
                                                                    d_element_type& Ezn,
                                                                    d_element_type& Bxn,
                                                                    d_element_type& Byn,
                                                                    d_element_type& Bzn
                                                                  )
{
    auto L2ProjCon = opProjection( _domainSpace=Ch, _imageSpace=Ch, _type=L2 );
    
    
    exporter->step( time )->setMesh( mesh );
    
    c_element_type Exc = Ch->element();
    c_element_type Eyc = Ch->element();
    c_element_type Bzc = Ch->element();
       
    Exc = L2ProjCon->project( idv( Exn ) );
    Eyc = L2ProjCon->project( idv( Eyn ) );
    Bzc = L2ProjCon->project( idv( Bzn ) );
    
    if(Dim==2){
    exporter->step( time )->add( "Ex", Exc );
    exporter->step( time )->add( "Ey", Eyc );
    exporter->step( time )->add( "Bz", Bzc );
    exporter->save();
    }
}























#endif /* __Maxwell_DG_H */



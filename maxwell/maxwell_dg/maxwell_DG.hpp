
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
    timers()
    {}

    FEELPP_DONT_INLINE
    void init();
    FEELPP_DONT_INLINE
    void run();

private:
    backend_ptrtype M_backend;
    double meshSize;
    double timeStep ;
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
    void assemble_RHS(  d_element_type &Exn,
                            d_element_type &Eyn, 
                            d_element_type &Ezn, 
                            d_element_type &Bxn,
                            d_element_type &Byn,
                            d_element_type &Bzn
                            );
    void solve(d_element_type& dtExn,
                d_element_type& dtEyn,
                d_element_type& dtEzn,
                d_element_type& dtBxn,
                d_element_type& dtByn,
                d_element_type& dtBzn
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
    auto dtExn =  Xh->element();
    if( Dim == 2)
    {
    double pi=M_PI;
    double k=2*pi; //=0;
    double theta=0;
    double vu=cos( theta );
    double vv=sin( theta );
    auto c=cos( k * ( vu * Px() + vv * Py() - cst_ref( time ) ) + cst( pi/2.0 ) );
    auto Ex_exact = -vv*c;
    auto Ey_exact = vu*c;
    auto Bz_exact = c;
    auto L2ProjDisc = opProjection( _domainSpace=Xh, _imageSpace=Xh, _type=L2 );
          
    Exn = L2ProjDisc->project( Ex_exact );
    Eyn = L2ProjDisc->project( Ey_exact );
    Bzn = L2ProjDisc->project( Bz_exact );
    this-> assemble_LHS();
    this->assemble_RHS(Exn,  Eyn, Ezn, Bxn, Byn, Bzn);
    this->solve(Exn,  Eyn, Ezn, Bxn, Byn, Bzn);
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
Maxwell_DG<Dim, Order_poly, Order_geo>::assemble_RHS(  d_element_type& Exn,
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
    auto wR = vec( rightfacev( idv( Exn ) ),rightfacev( idv( Eyn ) ),rightfacev( idv( Bzn ) ) );
    auto wL = vec( leftfacev( idv( Exn ) ),leftfacev( idv( Eyn ) ),leftfacev( idv( Bzn ) ) );
    auto wMetal = vec( -leftfacev( idv( Exn ) ),-leftfacev( idv( Eyn ) ),leftfacev( idv( Bzn ) ) );
      // FLux matrix
    auto Anp_1 = vec( +Ny() * Ny() / 2., -Nx() * Ny() / 2., -Ny() / 2. );
    auto Anp_2 = vec( -Nx() * Ny() / 2., Nx() * Nx() / 2., Nx() / 2. );
    auto Anp_3 = vec( -Ny() / 2., Nx() / 2., cst( 1. / 2. ) );
    auto Anm_1 = vec(  -Ny() * Ny() / 2., Nx() * Ny() / 2., -Ny() / 2. );
    auto Anm_2 = vec( Nx() * Ny() / 2., -Nx() * Nx() / 2., Nx() / 2. );
    auto Anm_3 = vec( -Ny() / 2., Nx() / 2., cst( -1. / 2. ) );
    auto v = Xh->element();
         
    auto lf_Ex=form1( _test=Xh, _vector=RHS_Ex );
    auto lf_Ey=form1( _test=Xh, _vector=RHS_Ey );
    // auto lf_Ez=form1( _test=Xh, _vector=RHS_Ez );
    // auto lf_Bx=form1( _test=Xh, _vector=RHS_Bx );
    // auto lf_By=form1( _test=Xh, _vector=RHS_By );
    auto lf_Bz=form1( _test=Xh, _vector=RHS_Bz );
        
    LOG(INFO) << "[assemble_RHS] :: Starting \n";
    timers["RHS"].first.restart();
    lf_Ex = integrate(  _range=elements( mesh ), 
                          _expr=id( v)*dyv( Bzn ) );
                              

    lf_Ex += integrate( _range=internalfaces( mesh ),
                                    _expr= ( trans( Anm_1 )*( wL-wR ) )*leftface( id( v ) )
                                                      + ( trans( Anp_1 )*( wL-wR ) )*rightface( id( v ) ) );

    lf_Ex += integrate(  _range=markedfaces( mesh, "Metal" ), 
                                        _expr=( trans( Anm_1 )*( wL-wMetal ) )*id( v ) );

    lf_Ex += integrate( _range=markedfaces( mesh, "Dirichlet" ), 
                                     _expr=( trans( Anm_1 )*( wL-w ) )*id( v ) );

    LOG(INFO) << "[assemble_RHS] :: Assemblage RHS -TM- Ex done \n";


    lf_Ey = integrate( _range=elements( mesh ),  
                                    _expr=-id( v)*dxv( Bzn ) );

    lf_Ey += integrate( _range=internalfaces( mesh ),
                                    _expr= ( trans( Anm_2 )*( wL-wR ) )*leftface( id( v ) )
                                                    + ( trans( Anp_2 )*( wL-wR ) )*rightface( id( v ) ) );

    lf_Ey += integrate( _range=markedfaces( mesh, "Metal" ),   
                                    _expr=( trans( Anm_2 )*( wL-wMetal ) )*id( v ) );

    lf_Ey += integrate(_range=markedfaces( mesh, "Dirichlet" ), 
                                    _expr=( trans( Anm_2 )*( wL-w ) )*id( v ) );


    LOG(INFO) << "[assemble_RHS] :: Assemblage RHS -TM- Ey done \n";

    lf_Bz = integrate(_range=elements( mesh ), 
                                _expr=-id( v)*dxv( Eyn ) + id( v)*dyv( Exn ) );

    lf_Bz += integrate( _range=internalfaces( mesh ),
                                    _expr= ( trans( Anm_3 )*( wL-wR ) )*leftface( id( v ) )
                                                    + ( trans( Anp_3 )*( wL-wR ) )*rightface( id( v ) ) );

    lf_Bz += integrate( _range=markedfaces( mesh, "Metal" ), 
                                    _expr=( trans( Anm_3 )*( wL-wMetal ) )*id( v ) );
    lf_Bz += integrate(_range=markedfaces( mesh, "Dirichlet" ), 
                                    _expr=( trans( Anm_3 )*( wL-w ) )*id( v ) );
        
   }
    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "[assemble_RHS] :: Completed ";
    }
} 
#endif /* __Maxwell_DG_H */



/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
      Date: 2011-06-01

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <feel/feel.hpp>
int main( int argc, char** argv ) {
  using namespace Feel;
  po::options_description options ( "Allowed Option") ;
  options.add_options()
  ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
  ( "Tfinal", po::value<double>()->default_value( 1 ), "final time" )
	( "gmsh.filename", po::value<std::string>()->default_value("diode-simplex.geo"), "geo filename" );
  Environment env(  _argc=argc, _argv=argv,
                                 _desc=options,
                                 _about=about(_name="maxwell",
                                 _author="Feel++ Consortium",
                                 _email="feelpp-devel@feelpp.org"));
 

  double Tfinal = option(_name="Tfinal").as<double>();
  const std::string geo_filename = option( _name="gmsh.filename" ).as<std::string>();
 // Load simplex geo
 // TODO :: Check if simplex of hypercube and load appropriate geo
  auto mesh = createGMSHMesh( _mesh=new Mesh<Simplex<2>>,
                                                        _desc=geo( _filename=geo_filename,
                                                        _h=option(_name="hsize").as<double>() ) );

                                                
  // auto mesh = unitSquare();                                              
                                                
  std::cout << "Geo and mesh Generated " << "\n";
  

 
  //TODO :: Add Exact solution, and init field with her
  double dt=0;
  auto Xh = Pdh<4>(mesh,true);
  auto Ch = Pch<1>(mesh);
  auto Exn = Xh->element();
  auto Eyn = Xh->element();
  auto Bzn = Xh->element();
  auto dtExn = Xh->element();
  auto dtEyn = Xh->element();
  auto dtBzn = Xh->element();
  auto v = Xh->element();
  auto u =  Xh->element();
  auto w = vec( idv( Exn ),idv( Eyn ),idv( Bzn ) );
  auto wR = vec( rightfacev( idv( Exn ) ),rightfacev( idv( Eyn ) ),rightfacev( idv( Bzn ) ) );
  auto wL = vec( leftfacev( idv( Exn ) ),leftfacev( idv( Eyn ) ),leftfacev( idv( Bzn ) ) );
  auto wMetal = vec( -leftfacev( idv( Exn ) ),-leftfacev( idv( Eyn ) ),leftfacev( idv( Bzn ) ) );
  
  auto a = form2( _test=Xh, _trial=Xh );
  auto lf_Ex=form1(_test=Xh);
  auto lf_Ey=form1(_test=Xh);
  auto lf_Bz=form1(_test=Xh);
  
  auto Anp_1 = vec( +Ny() * Ny() / 2., -Nx() * Ny() / 2., -Ny() / 2. );
  auto Anp_2 = vec( -Nx() * Ny() / 2., Nx() * Nx() / 2., Nx() / 2. );
  auto Anp_3 = vec( -Ny() / 2., Nx() / 2., cst( 1. / 2. ) );
  auto Anm_1 = vec(  -Ny() * Ny() / 2., Nx() * Ny() / 2., -Ny() / 2. );
  auto Anm_2 = vec( Nx() * Ny() / 2., -Nx() * Nx() / 2., Nx() / 2. );
  auto Anm_3 = vec( -Ny() / 2., Nx() / 2., cst( -1. / 2. ) );
  
  auto L2ProjC = opProjection( _domainSpace=Ch, _imageSpace=Ch, _type=L2 );
  auto L2ProjD = opProjection( _domainSpace=Xh, _imageSpace=Xh, _type=L2 );
  Exn = L2ProjD->project( cst(0.) );
  Eyn = L2ProjD->project( cst(0.) );
  Bzn = L2ProjD->project( cst(0.) );
    
  //Bilinear form (the trial function u is only here for pretty naming)
  a = integrate( elements( mesh ), idt( u )*id( v ) );   
                                  
  //Linear form for Ex,Ey,Bz
  lf_Ex = integrate(  _range=elements( mesh ), 
                              _expr=id( v)*dyv( Bzn ) );
  lf_Ex += integrate( _range=internalfaces( mesh ),
                               _expr= ( trans( Anm_1 )*( wL-wR ) )*leftface( id( v ) )
                                            + ( trans( Anp_1 )*( wL-wR ) )*rightface( id( v ) ) );
                                            
  lf_Ex += integrate(  _range=markedfaces( mesh, "Metal" ), 
                                _expr=( trans( Anm_1 )*( wL-wMetal ) )*id( v ) );
                                
  lf_Ex += integrate( _range=markedfaces( mesh, "Dirichlet" ), 
                               _expr=( trans( Anm_1 )*( wL-w ) )*id( v ) );
  
  lf_Ey = integrate( _range=elements( mesh ),  
                             _expr=-id( v)*dxv( Bzn ) );
  lf_Ey += integrate( _range=internalfaces( mesh ),
                               _expr= ( trans( Anm_2 )*( wL-wR ) )*leftface( id( v ) )
                                            + ( trans( Anp_2 )*( wL-wR ) )*rightface( id( v ) ) );
  
  lf_Ey += integrate( _range=markedfaces( mesh, "Metal" ),   
                              _expr=( trans( Anm_2 )*( wL-wMetal ) )*id( v ) );
                              
  lf_Ey += integrate(_range=markedfaces( mesh, "Dirichlet" ), 
                              _expr=( trans( Anm_2 )*( wL-w ) )*id( v ) );
  
  lf_Bz = integrate(_range=elements( mesh ), 
                            _expr=-id( v)*dxv( Eyn ) + id( v)*dyv( Exn ) );
  lf_Bz += integrate( _range=internalfaces( mesh ),
                               _expr= ( trans( Anm_3 )*( wL-wR ) )*leftface( id( v ) )
                                            + ( trans( Anp_3 )*( wL-wR ) )*rightface( id( v ) ) );
  
  lf_Bz += integrate( _range=markedfaces( mesh, "Metal" ), 
                               _expr=( trans( Anm_3 )*( wL-wMetal ) )*id( v ) );
  lf_Bz += integrate(_range=markedfaces( mesh, "Dirichlet" ), 
                              _expr=( trans( Anm_3 )*( wL-w ) )*id( v ) );

   std::cout << "Assemblage done " << "\n";
   a.solve( _solution=dtExn, _rhs=lf_Ex );
   a.solve( _solution=dtEyn, _rhs=lf_Ey );
   a.solve( _solution=dtBzn, _rhs=lf_Bz );
   std::cout << "Solve done " << "\n"; 
  //TODO :: Integrate petsc TS solver
  //                What about CFL ?   
 //Very simple explicit Euler scheme for time integration 
  Exn = Exn + dt*dtExn;
  Eyn = Eyn + dt*dtEyn;
  Bzn = Bzn + dt*dtExn;
  std::cout << "Time marching Done " << "\n"; 
    
}


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
  ( "verbose", po::value<bool>()->default_value( 1 ), "verbose output" )
  ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
  ( "dt", Feel::po::value<double>()->default_value( 0.01 ), "timestep value" )
  ( "Tfinal", Feel::po::value<double>()->default_value( 1 ), "final time" )
  ;
   Environment env(  _argc=argc, _argv=argv,
                                  _desc=options,
                                  _about=about(_name="maxwell",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
 
  //Only here for compilation test
  auto mesh = unitSquare();
  double dt;
  auto Xh = Pdh<4>(mesh);
  auto Exn = Xh->element();
  auto Eyn = Xh->element();
  auto Bzn = Xh->element();
  auto dtExn = Xh->element();
  auto dtEyn = Xh->element();
  auto dtBzn = Xh->element();
  
  auto Exnp1 = Xh->element();
  auto Eynp1 = Xh->element();
  auto Bznp1 = Xh->element();
  auto v = Xh->element();
  auto u =  Xh->element();
  auto w = vec( idv( Exn ),idv( Eyn ),idv( Bzn ) );
  auto wR = vec( rightfacev( idv( Exn ) ),rightfacev( idv( Eyn ) ),rightfacev( idv( Bzn ) ) );
  auto wL = vec( leftfacev( idv( Exn ) ),leftfacev( idv( Eyn ) ),leftfacev( idv( Bzn ) ) );
  auto wMetal = vec( -leftfacev( idv( Exn ) ),-leftfacev( idv( Eyn ) ),leftfacev( idv( Bzn ) ) );
  
  auto a = form2( _test=Xh, _trial=Xh );
  a = integrate( elements( mesh ), idt( u )*id( v ) );
    
  auto lf_Ex=form1(_test=Xh);
  auto lf_Ey=form1(_test=Xh);
  auto lf_Bz=form1(_test=Xh);
    
  auto Anp_1 = vec( +Ny() * Ny() / 2., -Nx() * Ny() / 2., -Ny() / 2. );
  auto Anp_2 = vec( -Nx() * Ny() / 2., Nx() * Nx() / 2., Nx() / 2. );
  auto Anp_3 = vec( -Ny() / 2., Nx() / 2., cst( 1. / 2. ) );
  auto Anm_1 = vec(  -Ny() * Ny() / 2., Nx() * Ny() / 2., -Ny() / 2. );
  auto Anm_2 = vec( Nx() * Ny() / 2., -Nx() * Nx() / 2., Nx() / 2. );
  auto Anm_3 = vec( -Ny() / 2., Nx() / 2., cst( -1. / 2. ) );

  lf_Ex = integrate(  _range=elements( mesh ), 
                              _expr=id( v)*dyv( Bzn ) );
  lf_Ex += integrate( _range=internalfaces( mesh ),
                               _expr= ( trans( Anm_1 )*( wL-wR ) )*leftface( id( v) )
                                            + ( trans( Anp_1 )*( wL-wR ) )*rightface( id( v) ) );
                                            
  lf_Ex += integrate(  _range=markedfaces( mesh, "Metal" ), 
                                _expr=( trans( Anm_1 )*( wL-wMetal ) )*id( v) );
                                
  lf_Ex += integrate( _range=markedfaces( mesh, "Dirichlet" ), 
                               _expr=( trans( Anm_1 )*( wL-w ) )*id( v) );
  
  lf_Ey = integrate( _range=elements( mesh ),  
                             _expr=-id( v)*dxv( Bzn ) );
  lf_Ey += integrate( _range=internalfaces( mesh ),
                               _expr= ( trans( Anm_2 )*( wL-wR ) )*leftface( id( v) )
                                            + ( trans( Anp_2 )*( wL-wR ) )*rightface( id( v) ) );
  
  lf_Ey += integrate( _range=markedfaces( mesh, "Metal" ),   
                              _expr=( trans( Anm_2 )*( wL-wMetal ) )*id( v) );
                              
  lf_Ey += integrate(_range=markedfaces( mesh, "Dirichlet" ), 
                              _expr=( trans( Anm_2 )*( wL-w ) )*id( v) );
  
  lf_Bz = integrate(_range=elements( mesh ), 
                            _expr=-id( v)*dxv( Eyn ) + id( v)*dyv( Exn ) );
  lf_Bz += integrate( _range=internalfaces( mesh ),
                               _expr= ( trans( Anm_3 )*( wL-wR ) )*leftface( id( v) )
                                            + ( trans( Anp_3 )*( wL-wR ) )*rightface( id( v) ) );
  
  lf_Bz += integrate( _range=markedfaces( mesh, "Metal" ), 
                               _expr=( trans( Anm_3 )*( wL-wMetal ) )*id( v) );
  lf_Bz += integrate(_range=markedfaces( mesh, "Dirichlet" ), 
                              _expr=( trans( Anm_3 )*( wL-w ) )*id( v) );

   a.solve( _solution=dtExn, _rhs=lf_Ex );
   a.solve( _solution=dtEyn, _rhs=lf_Ey );
   a.solve( _solution=dtBzn, _rhs=lf_Bz );
    
  Exnp1 = Exn + dt*dtExn;
  Eynp1 = Eyn + dt*dtEyn;
  Bznp1 = Bzn + dt*dtExn;
  
    
}


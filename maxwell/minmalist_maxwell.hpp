/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
This file is part of the Feel library
Author(s): Samuel Quinodoz
Christophe Prud'homme <christophe.prudhomme@feelpp.org>
Date: 2009-02-25
Copyright (C) 2007 Samuel Quinodoz
Copyright (C) 2009 Université Joseph Fourier (Grenoble I)
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3.0 of the License, or (at your option) any later version.
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/
#ifndef __Maxwell_DG
#define __Maxwell_DG 1

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;
#if !defined( DG_DIM )
#define CONVECTION_DIM 1
#endif
#if !defined( DG_GEO_ORDER )
#define DG_GEO_ORDER 1
#endif
#if !defined( DG_POLY_ORDER )
#define DG_POLY_ORDER 4
#endif
inline
po::options_description
makeOptions()
{
    po::options_description maxwell_options( "Maxwell_DG options" );
   maxwell_options.add_options()
  ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
  ( "Tfinal", po::value<double>()->default_value( 1 ), "final time" )
	( "gmsh.filename", po::value<std::string>()->default_value("diode-simplex.geo"), "geo filename" );
    ;
    return maxwell_options.add( Feel::feel_options() );
}

class Maxwell_DG : public Application
{
typedef Application super;
public:
  static const int Dim = DG_DIM;
  static const int Order_geo = DG_GEO_ORDER;
  static const int Order_poly = DG_POLY_ORDER;
  //Algebra
  typedef double value_type;
  typedef Backend<value_type> backend_type;
  typedef boost::shared_ptr<backend_type> backend_ptrtype;
  typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
  typedef backend_type::vector_ptrtype vector_ptrtype;
 
 //Mesh What for Hypercube ?
  typedef Simplex<Dim,Order_geo> entity_type;
  typedef Mesh<entity_type> mesh_type;
  typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

  //Discontinous
   #if defined( USE_LEGENDRE )
   typedef bases<Legendre<Order_poly,Scalar> > d_basis_type;
  #else
  typedef bases<Lagrange<Order_poly,Scalar,Discontinuous> > d_basis_type;
  #endif
  typedef FunctionSpace<mesh_type, d_basis_type> d_space_type;
  typedef boost::shared_ptr<d_space_type> d_space_ptrtype;
  typedef d_space_type::element_type d_element_type;

  //Continuous
  typedef bases<Lagrange<Order_poly,Scalar> > c_basis_type;
  typedef FunctionSpace<mesh_type, c_basis_type> c_space_type;
  typedef boost::shared_ptr<c_space_type> c_space_ptrtype;
  typedef c_space_type::element_type c_element_type;
 
  //Exporter
  typedef Exporter<mesh_type,Order_geo> export_type;
  typedef boost::shared_ptr<export_type> export_ptrtype;

 
   Maxwell_DG( int argc , char** argv , AboutData const& , po::options_description const& );
private:
    sparse_matrix_ptrtype D;
    vector_ptrtype Ex, Ey, Ez, Bx, By, Bz;
    mesh_ptrtype mesh;
    space_ptrtype Xh;
}
#endif

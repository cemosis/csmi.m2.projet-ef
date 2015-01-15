
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
   \file maxwell_dg_main.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \author Philippe Helluy <helluy@math.unistra.fr>
   \author Thomas Strub <strub@math.unistra.fr>
   \author Pierre GERHARD <pierre.gerhard@gmail.com>
   \date 2014-11-01
 */
#include "maxwell_DG.hpp"

inline
po::options_description
makeOptions()
{
    po::options_description maxwell_options( "Maxwell_DG options" );
    maxwell_options.add_options()
    ( "hsize", po::value<double>()->default_value( 0.25 ), "mesh size" )
    ( "CFL", po::value<double>()->default_value( 0.1 ), "CFL" )
    ("shape", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the domain (either simplex or hypercube)")
    ( "Tfinal", po::value<double>()->default_value( 1 ), "final time" )
    ( "gmsh.filename", po::value<std::string>()->default_value("circle2.geo"), "geo filename" )
    ( "t_m", po::value<int>()->default_value( 0 ), "Euler(=0), RK4(=1) method" )
    ;
    return maxwell_options.add( Feel::feel_options() );
}

inline AboutData
makeAbout()
{
    AboutData about( "Maxwell_DG",
                                "diode",
                                "0.2", 			
                                "Resolve 2D Maxwell equation with Discontinuous Galerkin",
                                AboutData::License_GPL ,	
                                "Copyright (c) 2014 Universite de Strasbourg" );

    about.addAuthor( "Christophe Prud'homme","Maintainer","prudhomme@unistra.fr" ,"" );
    about.addAuthor( "Philippe Helluy", "Developer", "helluy@math.unistra.fr", "" );
    about.addAuthor( "Thomas Strub", "Developer", "strub@math.unistra.fr", "" );
    about.addAuthor( "Pierre GERHARD","Developper","pierre.gerhard@gmail.com" ,"" );
    return about;
};


int
main( int argc, char** argv )
{
using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                                  _desc=makeOptions(),
                                  _about=makeAbout()  );
                                  
    Application app; 
    
     if ( app.vm().count( "help" ) ){
    LOG(INFO) << app.optionsDescription() << "\n";
    return 0;
    }
    app.add( new Maxwell_DG<2,2,1>() );
    app.run();
}

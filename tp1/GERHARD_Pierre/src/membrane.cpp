#include <feel/feel.hpp>
int main( int argc, char** argv ) {
  using namespace Feel;
  po::options_description options ( "Membrane options ") ;
  options.add_options()
  ( "T", po::value<double>()->default_value(1.0), "Tension" )
  ( "theta", po::value<double>()->default_value(0.2), "theta" )
  // ( "y0", po::value<double>()->default_value(0.0), "y-coordinate" )
  ( "sigma", po::value<double>()->default_value(0.02), "sigma" )
  ( "R", po::value<double>()->default_value(0.3), "Radius" )
  ( "A", po::value<double>()->default_value(1.0), "Amplitude" );
   Environment env( _argc=argc, _argv=argv,
                              _desc=options,
                              _about=about(_name="membrane",
                              _author="Feel++ Consortium",
                              _email="feelpp-devel@feelpp.org"));
                            
  double T = option(_name="T").as<double>();
  // double x0 = option(_name="x0").as<double>(); 
  double theta = option(_name="theta").as<double>(); 
  // double y0  = option(_name="y0").as<double>();
  double sigma = option(_name="sigma").as<double>();
  double R = option(_name="R").as<double>();
  double A = option(_name="A").as<double>();
  double x0 = 0.5*R*cos(theta);
  double y0 = 0.5*R*sin(theta);
  auto mesh = unitCircle();
  auto Xh = Pch<2>(mesh);
  auto w = Xh->element();
  auto v = Xh->element();
  auto d = Xh->element();
  auto we = expr("1-x*x-y*y:x:y");
  auto f=expr("4*exp(-0.5*(pow((R*x-x0)/sigma,2)) - 0.5*(pow((R*y-y0)/sigma,2))):x:y:R:sigma:x0:y0"); 
  //linear and bilinear form
  auto l = form1(_test=Xh);
  l = integrate(_range=elements(mesh),
                     _expr= f*id(v) );
  auto a = form2(_test=Xh,_trial=Xh);
  a = integrate(_range=elements(mesh),
                      _expr=gradt(w)*trans(grad(v)) );
  a+=on(_range=boundaryfaces(mesh),
             _rhs=l, 
             _element=w, 
             _expr=cst(0.));

  a.solve( _solution=w, _rhs=l );
  d=(A*R*R)/(8*pi*sigma*T)*w;
  std::cout << "|w-we|= "
                << normL2(_range=elements(mesh),_expr=(idv(w)-we) )
                << "\n";

  std::cout << "max(w) = "
                << normLinf( _range=elements(mesh), _expr=idv(w), _pset=_Q<5>() )
                << "\n";
                
  std::cout << "max(d) = "
               << normLinf( _range=elements(mesh), _expr=idv(d), _pset=_Q<5>() )
               << "\n";          
             
  auto e = exporter( _mesh=mesh ); 
  e->add("w", w );
  e->add("f",vf::project(_space=Xh, _range=elements(mesh), _expr=f));
  e->add("d", d);
  // v.on(_range=elements(mesh),_expr=we );
  // e->add( "we", v );
  e->save();
}

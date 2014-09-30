#include <feel/feel.hpp>
int main( int argc, char** argv ) {
  using namespace Feel;
  Environment env( _argc=argc, _argv=argv,
                   _about=about(_name="lap2",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));
  auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
  // $\mathbb{P}_{2}$ finite element space on
  // triangular elements
  auto Xh = Pch<2>(mesh);
  auto u = Xh->element();
  auto v = Xh->element();
  auto ue = expr("1+x*x+2*y*y:x:y");
  auto f = expr("-6:x:y");
  // $\int_\Omega f v$
  auto l = form1(_test=Xh);
  l= integrate( elements(mesh), f*id(v) );
  // $\int_\Omega \nabla u \cdot \nabla v$
  auto a = form2(_test=Xh,_trial=Xh);
  a = integrate( _range=elements(mesh),
                 _expr=gradt(u)*trans(grad(v)) );
  a+=on( boundaryfaces(mesh), u, l, ue );
  // solve a( u, v ) = l( v )
  a.solve( _solution=u, _rhs=l );

  std::cout << "|u-ue|=$"
            << normL2(_range=elements(mesh),
                      _expr=(idv(u)-ue) )
            << "\n";

  auto e = exporter( _mesh=mesh ); 
  e->add( "u", u );
  // v interpolate ue
  v.on(_range=elements(mesh),_expr=ue );
  e->add( "ue", v );
  e->save();
}

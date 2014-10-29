#include<feel/feel.hpp>

int main(int argc,char **argv) {

	using namespace Feel;
	using namespace Feel::vf;
	/*v
	   Setting option for the membrane application
	 */
	po::options_description myoptions("membrane options");
	myoptions.add_options()
		("T",po::value<double>()->default_value(1.0),"Tension") // double >() - > default_value ( 1.0 ) , " Tension " )
		("x0",po::value<double>()->default_value(0.0),"x-coordinate") 
		("y0",po::value<double>()->default_value(0.0),"y-coordinate") 
		("sigma",po::value<double>()->default_value(0.02),"sigma")
		("R",po::value<double>()->default_value(0.3),"Radius") 
		("theta",po::value<double>()->default_value(0.0),"angle") 
		("A",po::value<double>()->default_value(1.0),"Amplitude");
	Environment env ( _argc = argc , _argv = argv ,
			_desc =myoptions  ,
			_about=about(_name="membrane",
				     _author="Feel++ Consortium",
				_email="feelpp-devel@feelpp.org"));
         /*
           Create mesh over unit circle */
         auto mesh = unitCircle();
       
        /* Functions spaces */
        auto Vh = Pch<2>(mesh) ;
       
        /* Element of functions spaces */
        auto u = Vh->element();
        auto v = Vh->element();
        auto F = Vh->element();
        auto deflection = Vh->element();

//f=exp(-0.5*(pow((R*x-x0)/sigma,2)))*(-4)*exp(-0.5*(pow((R*y-y0)/sigma,2))):x:y:T:sigma:x0:y0
         /* 
          compute options */
       //  auto         f = expr( soption(_name ="functions.f"),"f");

         double       T = doption(_name = "T");
	     double      x0 = doption(_name = "x0");
	     double      y0 = doption(_name = "y0");
	     double   sigma = doption(_name = "sigma");
	     double      R  = doption(_name = "R");
	     double  theta  = doption(_name = "theta");
	     double       A = doption(_name  = "A"); 
                    x0  = 0.6*R*cos(theta);
                     y0 = 0.6*R*sin(theta);
         auto        f  =   expr("exp(-0.5*(pow((R*x-x0)/sigma,2)))*(-4)*exp(-0.5*(pow((R*y-y0)/sigma,2))):x:y:x0:y0:R:sigma");
         //  std::cout<<"f = "<<f<<"\n";
         F = project(_space = Vh, _range =elements(mesh),_expr =f);
	/* linear form */
      	auto l = form1(_test= Vh);
             l = integrate (_range = elements(mesh) ,_expr = f*id(v));
       /* bilinear form */
       auto  a = form2(_trial =Vh,_test =Vh);
             a = integrate(_range =elements(mesh) ,_expr=T*(gradt(u)*trans(grad(v)))) ;
             a += on(_range = boundaryfaces(mesh),_rhs = l,_element =u,_expr = cst(0.));
             a.solve(_rhs=l,_solution=u);
     auto  dfl = A*R*R*idv(u)/(8*pi*sigma*T);
     auto  maximum = max(idv(u),dfl);
/*     deflection = project(_space =Vh,_range =elements(mesh),_expr =maximum); */
       deflection =vf:: project(_space = Vh, _range =elements(mesh),_expr =maximum);
       double L2Error = normL2(_range = elements(mesh),_expr=(maximum));
       std:: cout << "||Erreur||_2  " << L2Error<<"\n";
	/* add exporter */
        auto e = exporter(_mesh =mesh);
	     e->add("u",u);
         e->add("F",F);
        e->add("dfl",deflection);
         e->save();
	return 0;


}



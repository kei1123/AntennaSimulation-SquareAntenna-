//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
//			•¡‘f”ŒvZ—pƒwƒbƒ_[@@complex.h
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////
//				ƒNƒ‰ƒX
//////////////////////////////////////////////////////////////////////////////////

class complex{									//complexƒNƒ‰ƒX
public:
	double re;									//realiÀ”•”j
	double im;									//imaginaryi‹•”•”j
};


//////////////////////////////////////////////////////////////////////////////////
//				•¡‘f”‚ÌŠeí‰‰Z
//////////////////////////////////////////////////////////////////////////////////
complex Complex(double x,double y)				//À”‚©‚ç•¡‘f”‚ğì‚é
{
	complex z;
	z.re = x;
	z.im = y;
	return z;
}

complex operator+(complex x,complex y)			//•¡‘f”‚Æ•¡‘f”‚Ì‰ÁZ
{
	complex z;
	z.re = x.re + y.re;
	z.im = x.im + y.im;
	return z;
}

complex operator+(complex x,double y)			//•¡‘f”‚ÆÀ”‚Ì‰ÁZ(•¡‘f”+À”)
{
	complex z;
	z.re = y + x.re;
	z.im = x.im;
	return z;
}

complex operator+(double x,complex y)			//•¡‘f”‚ÆÀ”‚Ì‰ÁZ(À”+•¡‘f”)
{
	complex z;
	z.re = x + y.re ;
	z.im = y.im;
	return z;
}

complex operator-(complex x,complex y)			//•¡‘f”‚Æ•¡‘f”‚ÌŒ¸Z
{
	complex z;
	z.re = x.re - y.re;
	z.im = x.im - y.im;
	return z;
}

complex operator-(complex x,double y)			//•¡‘f”‚ÆÀ”  ‚ÌŒ¸Z(•¡‘f”-À”)
{
	complex z;
	z.re = x.re -y;
	z.im = x.im;
	return z;
}

complex operator-(double x ,complex y)			//•¡‘f”‚ÆÀ”  ‚ÌŒ¸Z(À”-•¡‘f”)
{   
	complex z;
	z.re = x - y.re;
	z.im = 0.0 - y.im;
	return z;
}

complex operator*(complex x,complex y)			//•¡‘f”‚Æ•¡‘f”‚ÌæZ
{
	complex z;
	z.re = x.re*y.re - x.im*y.im;
	z.im = x.re*y.im + x.im*y.re;
	return z;
}

complex operator*(complex x,double y)			//•¡‘f”‚ÆÀ”‚ÌæZ(•¡‘f”*À”)
{
	complex z;
	z.re = y * x.re;
	z.im = y * x.im;
	return z;
}

complex operator*(double x,complex y)			//•¡‘f”‚ÆÀ”‚ÌæZ(À”*•¡‘f”)
{
	complex z;
	z.re = x * y.re;
	z.im = x * y.im;
	return z;
}

complex operator/(complex x,complex y)			//•¡‘f”‚Æ•¡‘f”‚ÌœZ
{                                             
	complex z;
	z.re = (x.re*y.re + x.im*y.im)/(y.re*y.re + y.im*y.im);
	z.im = (x.im*y.re - x.re*y.im)/(y.re*y.re + y.im*y.im);
	return z;
}

complex operator/(complex x,double y)			//•¡‘f”‚ÆÀ”‚ÌœZ(•¡‘f”€À”)
{
	complex z;
	z.re = x.re / y;
	z.im = x.im / y;
	return z;
}

complex operator/(double x,complex y)			//•¡‘f”‚ÆÀ”‚ÌœZ(À”€•¡‘f”)
{
	complex z;
	z.re = (x*y.re)/(y.re*y.re + y.im*y.im);
	z.im = (-1.0* x*y.im)/(y.re*y.re + y.im*y.im);
	return z;
}

complex Exp(complex x)							//•¡‘f”‚Ìw”ŠÖ”
{												//exp(jb)=cos(b)+jsin(b)
	complex z;                             
	z.re = exp(x.re)*cos(x.im);
	z.im = exp(x.re)*sin(x.im);
	return z;
}

double Abs(complex x)							//•¡‘f”‚Ìâ‘Î’l
{
	return sqrt(x.re*x.re + x.im*x.im);
}
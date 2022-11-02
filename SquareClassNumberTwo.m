//This is an implementation of the algorithm applied in Remark 2.5

Is_in_Z:=function(phi,Aut,L)
    t:=true;
	for k in [1..#Aut] do
		psi:=Aut[k];
		if psi(L[4]) eq L[4] then
		    for l in L do t:=t and (phi(psi(l)) eq psi(phi(l))); end for;
		end if;
	end for;    
    return t;
end function;    
    
Is_in_Gplus:=function(phi,j)
    return j eq phi(j);
end function;

has_order_2:=function(phi,L)
    t:=true;
    for l in L do t:=t and (phi(phi(l)) eq l); end for;
    return t; 
end function;

GenTwoCurves:=procedure(d);
	K<j>:=NumberField(HilbertClassPolynomial(-d));
	E:=EllipticCurveFromjInvariant(j);
	print(d);
	WE:=WeierstrassModel(E);
	a:=Coefficients(WE)[4];
	b:=Coefficients(WE)[5];

	P<t>:=PolynomialRing(K);
	H<sd>:=NumberField(t^2+d);
	PP<x>:=PolynomialRing(H);
	H2:=SplittingField(x^3+a*x+b);
	Q<z>:=PolynomialRing(H2);

	Arr:=Roots(z^3+a*z+b);
	alpha1:=Arr[1][1];
	alpha2:=Arr[2][1];
	alpha3:=Arr[3][1];
	
	Aut:=Automorphisms(H2);

	Sigmas:=[];
	for i in [1..#Aut] do
		phi:=Aut[i];
		if not(Is_in_Gplus(phi,j)) and Is_in_Z(phi,Aut,[alpha1,alpha2,alpha3,j,sd]) and has_order_2(phi,[alpha1,alpha2,alpha3,j,sd]) then Sigmas:=Append(Sigmas,phi); end if;
	end for;
	print "Elements in G-Gplus commuting with Gplus of order 2:", #Sigmas;
	
	for phi in Sigmas do

		beta1:=phi(alpha1);
		beta2:=phi(alpha2);
		beta3:=phi(alpha3);

		a1 := (alpha3-alpha2)^2/(beta3-beta2) + (alpha2-alpha1)^2/(beta2-beta1) +  (alpha1-alpha3)^2/(beta1-beta3);
		b1 := (beta3-beta2)^2/(alpha3-alpha2) + (beta2-beta1)^2/(alpha2-alpha1) +  (beta1-beta3)^2/(alpha1-alpha3);

		a2 := alpha1*(beta3-beta2) + alpha2*(beta1-beta3) + alpha3*(beta2-beta1);
		b2 := beta1*(alpha3-alpha2) + beta2*(alpha1-alpha3) + beta3*(alpha2-alpha1);

		A:= Discriminant((z-beta1)*(z-beta2)*(z-beta3))*a1/a2;
		B:= Discriminant((z-alpha1)*(z-alpha2)*(z-alpha3))*b1/b2;

		h:= - (A*(alpha2-alpha1)*(alpha1-alpha3)*z^2 + B*(beta2-beta1)*(beta1-beta3))*(A*(alpha3-alpha2)*(alpha2-alpha1)*z^2 + B*(beta3-beta2)*(beta2-beta1))*(A*(alpha1-alpha3)*(alpha3-alpha2)*z^2 + B*(beta1-beta3)*(beta3-beta2));
		c6:=Coefficients(Q!h)[7];
		c4:=Coefficients(Q!h)[5];
		c2:=Coefficients(Q!h)[3];
		c0:=Coefficients(Q!h)[1];
	
		c0:=K!(c0/c6);
		c2:=K!(c2/c6);
		c4:=K!(c4/c6);
	
		h:=P!(x^6+c4*x^4+c2*x^2+c0);
	
		J2:=IgusaInvariants(h)[1];
		J4:=IgusaInvariants(h)[2];
		J6:=IgusaInvariants(h)[3];
		J8:=IgusaInvariants(h)[4];
		J10:=IgusaInvariants(h)[5];

		g1 := J2^5/J10;
		g2:= J2^3*J4/J10;
		g3:= J2^2*J6/J10;

C := ReducedWamelenModel(IntegralModel(HyperellipticCurveFromG2Invariants([Rationals()!g1,Rationals()!g2,Rationals()!g3])));
		print C;
	end for;
end procedure;


D2:=[15,20,24,35,40,51,52,88,91,115,123,148,187,232,235,267,403,427];

for d in D2 do
	GenTwoCurves(d);
end for;


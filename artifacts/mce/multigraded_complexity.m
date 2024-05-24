load "../trilinear_XL_estimator.m";

precision := 34;

P<a> := PowerSeriesRing(Integers(), precision);
P<b> := PowerSeriesRing(P, precision);
P<c> := PowerSeriesRing(P, precision);

// d_solv old parameters
for n in [14,22,30] do
	m := n;
	k := n;

	ms := 1 / ((1-a)^(n) * (1-b)^(m) * (1-c)^(k));
	hs := ms * (1 - a*b)^k * (1 - a*c)^m * (1 - b*c)^n * (1-a*b*c)^(-2);

	monomials, d_solv := TrilinearXLComplexity(hs, ms, precision);
    print "d_solv", n, m, k, ChangePrecision(Log(2, n*m*monomials^2),4);//, "d_solv:", d_solv;
end for;

// d_solv new parameters
for n in [26,35,45] do
	m := n-1;
	k := n-1;

	ms := 1 / ((1-a)^(n) * (1-b)^(m) * (1-c)^(k));
	hs := ms * (1 - a*b)^k * (1 - a*c)^m * (1 - b*c)^n * (1-a*b*c)^(-2);

	monomials, d_solv := TrilinearXLComplexity(hs, ms, precision);
    print "d_solv", n, m, k, ChangePrecision(Log(2, n*m*monomials^2),4);//, "d_solv:", d_solv;
end for;

// d_ff old parameters
for n in [14,22,30] do
	m := n;
	k := n;

	ms := 1 / ((1-a)^(n-1) * (1-b)^(m-1) * (1-c)^(k-1));
	hs := ms * (1 - a*b)^k * (1 - a*c)^m * (1 - b*c)^n * (1-a*b*c)^(-2);

	monomials, d_ff := TrilinearXLComplexity(hs, ms, precision);
    print "d_ff  ", n, m, k, ChangePrecision(Log(2, n*m*monomials^2),4);//, "d_ff:  ", d_ff;
end for;

// d_ff new parameters
for n in [26,35,45] do
	m := n-1;
	k := n-1;

	ms := 1 / ((1-a)^(n-1) * (1-b)^(m-1) * (1-c)^(k-1));
	hs := ms * (1 - a*b)^k * (1 - a*c)^m * (1 - b*c)^n * (1-a*b*c)^(-2);

	monomials, d_ff := TrilinearXLComplexity(hs, ms, precision);
    print "d_ff  ", n, m, k, ChangePrecision(Log(2, n*m*monomials^2),4);//, "d_ff:  ", d_ff;
end for;
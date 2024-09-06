load "../trilinear_XL_estimator.m";

precision := 18;
P<t> := PowerSeriesRing(Integers(), precision);
P<s> := PowerSeriesRing(P, precision);
P<u> := PowerSeriesRing(P, precision);

for n in [13, 20, 25] do

    density := (n-2)*(n-3) + 1;

	ms := 1 / ((1-t)^(n-2) * (1-s)^(n-2) * (1-u)^(n-2));
	hs := ms * ((1-t*s)*(1-u*s)*(1-t*u))^n * ((1-s*u*t)^2 * (1-s*s*u) * (1-s*s*t) * (1-u*u*t) * (1-u*u*s) * (1-t*t*u) * (1-t*t*s))^(-1);

	monomials, d_solv := TrilinearXLComplexity(hs, ms, precision);
    print n, "d_solv", ChangePrecision(Log(2, density*monomials^2),4);//, "d_solv:", d_solv;

    // ==================

    ms := 1 / ((1-t)^(n-3) * (1-s)^(n-3) * (1-u)^(n-3));
	hs := ms * ((1-t*s)*(1-u*s)*(1-t*u))^n * ((1-s*u*t)^2 * (1-s*s*u) * (1-s*s*t) * (1-u*u*t) * (1-u*u*s) * (1-t*t*u) * (1-t*t*s))^(-1);

	monomials, d_ff := TrilinearXLComplexity(hs, ms, precision);
    print n, "d_ff  ", ChangePrecision(Log(2, density*monomials^2),4);//, "d_ff  :", d_ff;
end for;

exit;
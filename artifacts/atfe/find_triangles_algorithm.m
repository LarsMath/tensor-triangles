function FindTriangles(phi : verbosity:=0)
    F := Parent(phi[1][1][1]);
    n := #phi;

    R<[vars]> := PolynomialRing(F, 3*(n-3), "grevlex");
    x := [1,0,0] cat vars[1..n-3];
    y := [0,1,0] cat vars[n-2..2*n-6];
    z := [0,0,1] cat vars[2*n-5..3*n-9];

    system :=   [&+[phi[i][j][l] * x[i] * y[j] : i in [1..n], j in [1..n]] : l in [1..n]];
    system cat:=[&+[phi[i][j][l] * y[i] * z[j] : i in [1..n], j in [1..n]] : l in [1..n]];
    system cat:=[&+[phi[i][j][l] * z[i] * x[j] : i in [1..n], j in [1..n]] : l in [1..n]];

    SetVerbose("Groebner", verbosity);
    solutions := Variety(Ideal(system));
    SetVerbose("Groebner", 0);
    return solutions;
end function;

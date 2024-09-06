function FindTriangles(C : verbosity:=0)
    F := Parent(C[1][1][1]);
    n := #C;
    m := #C[1];
    k := #C[1][1];

    R<[vars]> := PolynomialRing(F, (n-1) + (m-1) + (k-1), "grevlex");
    x := [1] cat vars[1..(n-1)];
    y := [1] cat vars[(n-1)+1..(n-1)+(m-1)];
    z := [1] cat vars[(n-1)+(m-1)+1..(n-1)+(m-1)+(k-1)];

    system :=   [&+[C[i][j][l] * x[i] * y[j] : i in [1..n], j in [1..m]] : l in [1..k]];
    system cat:=[&+[C[i][j][l] * y[j] * z[l] : j in [1..m], l in [1..k]] : i in [1..n]];
    system cat:=[&+[C[i][j][l] * z[l] * x[i] : l in [1..k], i in [1..n]] : j in [1..m]];

    SetVerbose("Groebner", verbosity);
    solutions := Variety(Ideal(system));
    SetVerbose("Groebner", 0);
    return solutions;
end function;

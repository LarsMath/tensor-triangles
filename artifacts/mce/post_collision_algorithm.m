load "../tools.m";

// Pre-condition: C and D have a triangle in standard position
function PostCollision(C, D : verbosity:=0)

    F := Parent(C[1][1][1]);
    n := #C;
    m := #C[1];
    k := #C[1][1];

    for i in [1..n] do assert C[i][1][1] eq 0 and D[i][1][1] eq 0; end for;
    for j in [1..m] do assert C[1][j][1] eq 0 and D[1][j][1] eq 0; end for;
    for l in [1..k] do assert C[1][1][l] eq 0 and D[1][1][l] eq 0; end for;

    // ========================================================

    R<[x]> := PolynomialRing(F, 2 * (1 + n*(n-1) + m*(m-1) + k*(k-1)), "grevlex");

    // Build variable matrices
    A := HorizontalJoin(Matrix(R, n, 1, [[i eq 1 select x[2*(n*(n-1) + m*(m-1) + k*(k-1))+1] else 0] : i in [1..n]]), Matrix(R, n, n-1, x[1..n*(n-1)]));
    Ai := HorizontalJoin(Matrix(R, n, 1, [[i eq 1 select x[2*(n*(n-1) + m*(m-1) + k*(k-1))+2] else 0] : i in [1..n]]), Matrix(R, n, n-1, x[n*(n-1)+1..2*n*(n-1)]));
    B := HorizontalJoin(Matrix(R, m, 1, [[i eq 1 select 1 else 0] : i in [1..m]]), Matrix(R, m, m-1, x[2*n*(n-1)+1..2*n*(n-1)+m*(m-1)]));
    Bi := HorizontalJoin(Matrix(R, m, 1, [[i eq 1 select 1 else 0] : i in [1..m]]), Matrix(R, m, m-1, x[2*n*(n-1)+m*(m-1)+1..2*n*(n-1)+2*m*(m-1)]));
    T := HorizontalJoin(Matrix(R, k, 1, [[i eq 1 select 1 else 0] : i in [1..k]]), Matrix(R, k, k-1, x[2*n*(n-1)+2*m*(m-1)+1..2*n*(n-1)+2*m*(m-1)+k*(k-1)]));
    Ti := HorizontalJoin(Matrix(R, k, 1, [[i eq 1 select 1 else 0] : i in [1..k]]), Matrix(R, k, k-1, x[2*n*(n-1)+2*m*(m-1)+k*(k-1)+1..2*n*(n-1)+2*m*(m-1)+2*k*(k-1)]));

    //=========================== Build system ========================================

    In := ScalarMatrix(R, n, 1);
    Im := ScalarMatrix(R, m, 1);
    Ik := ScalarMatrix(R, k, 1);

    phi_AB := ApplyGroupAction3TI(A, B, Ik, C, R);
    phi_AT := ApplyGroupAction3TI(A, Im, T, C, R);
    phi_BT := ApplyGroupAction3TI(In, B, T, C, R);
    phi_A := ApplyGroupAction3TI(A, Im, Ik, C, R);
    phi_B := ApplyGroupAction3TI(In,  B, Ik, C, R);
    phi_T := ApplyGroupAction3TI(In, Im, T, C, R);

    psi_AB := ApplyGroupAction3TI(Ai, Bi, Ik, D, R);
    psi_AT := ApplyGroupAction3TI(Ai, Im, Ti, D, R);
    psi_BT := ApplyGroupAction3TI(In, Bi, Ti, D, R);
    psi_A := ApplyGroupAction3TI(Ai, Im, Ik, D, R);
    psi_B := ApplyGroupAction3TI(In, Bi, Ik, D, R);
    psi_T := ApplyGroupAction3TI(In, Im, Ti, D, R);

    AAi := Matrix(A) * Matrix(Ai);
    BBi := Matrix(B) * Matrix(Bi);
    TTi := Matrix(T) * Matrix(Ti);
    AiA := Matrix(Ai) * Matrix(A);
    BiB := Matrix(Bi) * Matrix(B);
    TiT := Matrix(Ti) * Matrix(T);

    system := [];
    system := system cat [phi_A[i][j][l] - psi_BT[i][j][l] : i in [1..n], j in [1..m], l in [1..k]];
    system := system cat [phi_B[i][j][l] - psi_AT[i][j][l] : i in [1..n], j in [1..m], l in [1..k]];
    system := system cat [phi_T[i][j][l] - psi_AB[i][j][l] : i in [1..n], j in [1..m], l in [1..k]];
    system := system cat [phi_AB[i][j][l] - psi_T[i][j][l] : i in [1..n], j in [1..m], l in [1..k]];
    system := system cat [phi_AT[i][j][l] - psi_B[i][j][l] : i in [1..n], j in [1..m], l in [1..k]];
    system := system cat [phi_BT[i][j][l] - psi_A[i][j][l] : i in [1..n], j in [1..m], l in [1..k]];

    system := system cat [AAi[i][j] - In[i][j] : i in [1..n], j in [1..n]];
    system := system cat [BBi[i][j] - Im[i][j] : i in [1..m], j in [1..m]];
    system := system cat [TTi[i][j] - Ik[i][j] : i in [1..k], j in [1..k]];
    system := system cat [AiA[i][j] - In[i][j] : i in [1..n], j in [1..n]];
    system := system cat [BiB[i][j] - Im[i][j] : i in [1..m], j in [1..m]];
    system := system cat [TiT[i][j] - Ik[i][j] : i in [1..k], j in [1..k]];

    // print "linears", #Reduce([f : f in system | Degree(f) eq 1]);
    // print "quadratics", #Reduce([f : f in system | Degree(f) eq 2]); // Not necessary but gives some statistics

    SetVerbose("Groebner", verbosity);
    solutions := Variety(Ideal(system));
    SetVerbose("Groebner", 0);

    return solutions;
end function;

load "../tools.m";

// Pre-condition: phi and psi have a triangle in standard position
function PostCollision(phi, psi : wtA12:=3, wtA22:=3, verbosity:=0)

    F := Parent(phi[1][1][1]);
    n := #phi;

    for i in [1..n] do 
        assert phi[1][2][i] eq 0 and phi[2][3][i] eq 0 and phi[3][1][i] eq 0;
        assert psi[1][2][i] eq 0 and psi[2][3][i] eq 0 and psi[3][1][i] eq 0;
    end for;

    // ================================================================

    R<[x]> := PolynomialRing(F, [1 : _ in [1..9]] cat [wtA12 : _ in [1..(n-3)*3]] cat [wtA22 : _ in [1..(n-3)^2]] cat [1 : _ in [1..9]] cat [wtA12 : _ in [1..(n-3)*3]] cat [wtA22 : _ in [1..(n-3)^2]]);

    // Note: Bx representes the inverse of A
    Ax := Matrix([(i le 3 select [x[j+3*(i-1)] : j in [1..3]] else [0,0,0]) cat [x[j + (n-3)*(i-1) + 9] : j in [1..n-3]] : i in [1..n]]);
    Bx := Matrix([(i le 3 select [x[j+3*(i-1) + 9 + n*(n-3)] : j in [1..3]] else [0,0,0]) cat [x[j + (n-3)*(i-1) + 9 + n*(n-3) + 9] : j in [1..n-3]] : i in [1..n]]);

    //=========================== Build system ========================================

    I := ScalarMatrix(R, n, 1);

    phi_AA := ApplyGroupAction3TI(Ax, Ax, I, phi, R);
    phi_A := ApplyGroupAction3TI(Ax, I, I, phi, R);

    psi_BB := ApplyGroupAction3TI(I, Bx, Bx, psi, R);
    psi_B := ApplyGroupAction3TI(I, I, Bx, psi, R);

    AB := Matrix(Ax) * Matrix(Bx);
    BA := Matrix(Bx) * Matrix(Ax);

    system := Eltseq(AB - I) cat Eltseq(BA - I);
    system := system cat [phi_A[i][j][l] - psi_BB[i][j][l] : i in [1..n], j in [1..n], l in [1..n]];
    system := system cat [phi_AA[i][j][l] - psi_B[i][j][l] : i in [1..n], j in [1..n], l in [1..n]];

    SetVerbose("Groebner", verbosity);
    solutions := Variety(Ideal(GroebnerBasis(system, 6)));
    SetVerbose("Groebner", 0);

    return solutions;
end function;

function MoveTriangleGeneralPosition(phi, x, y, z)
    F := Parent(phi[1][1][1]);
    n := #phi;

    for l in [1..n] do
        assert &+[phi[i][j][l] * x[i] * y[j] : i in [1..n], j in [1..n]] eq 0;
        assert &+[phi[i][j][l] * y[i] * z[j] : i in [1..n], j in [1..n]] eq 0;
        assert &+[phi[i][j][l] * z[i] * x[j] : i in [1..n], j in [1..n]] eq 0;
    end for;

    while true do
        A := Transpose(Matrix(F, n, n, [x, y, z] cat [[Random(F) : _ in [1..n]] : _ in [1..n-3]]));
        if Rank(A) eq n then
            return ApplyGroupAction3TI(A, A, A, phi, F), A;
        end if;
    end while;
end function;

function GenerateMonomials(vars, degree)
    monomials := [{}, {1}];
    for i in [1..degree] do
        monomials := monomials cat [{m * x : m in monomials[i+1], x in vars}];
    end for;
    return monomials;
end function;

function ApplyGroupAction3TI(A, B, T, C, R)
    n, m, k := Explode(<Nrows(A), Nrows(B), Nrows(T)>);
    C_slices := [Matrix(R, m, k, [[C[i][j][l] : l in [1..k]] : j in [1..m]]) : i in [1..n]];
    return [&+[A[i][ii] * (Transpose(B) * C_slices[i] * T) : i in [1..n]] : ii in [1..n]];
end function;

function Antisymmetrize(phi)
    return [[[phi[i][j][k] + phi[j][k][i] + phi[k][i][j] - phi[i][k][j] - phi[k][j][i] - phi[j][i][k] : i in [1..#phi]] : j in [1..#phi]] : k in [1..#phi]];
end function;

function GenerateATFEKeyPair(F, n)
    phi := Antisymmetrize([[[Random(F) : _ in [1..n]] : _ in [1..n]] : _ in [1..n]]);
    A := RandomMatrix(F, n, n); while Rank(A) ne n do A := RandomMatrix(F, n, n); end while;
    psi := ApplyGroupAction3TI(A, A, A, phi, F);
    return phi, psi;
end function;

function GenerateATFEKeyPairWithPlantedTriangle(F, n)
    phi_0 := Antisymmetrize([[[Random(F) : _ in [1..n]] : _ in [1..n]] : _ in [1..n]]);
    for i in [1..3] do
        for j in [1..3] do
            for k in [1..n] do
                phi_0[i][j][k] := 0;
                phi_0[k][i][j] := 0;
                phi_0[j][k][i] := 0;
            end for;
        end for;
    end for;

    A_0 := RandomMatrix(F, n, n); while Rank(A_0) ne n do A_0 := RandomMatrix(F, n, n); end while;
    B_0 := RandomMatrix(F, n, n); while Rank(B_0) ne n do B_0 := RandomMatrix(F, n, n); end while;

    phi := ApplyGroupAction3TI(A_0, A_0, A_0, phi_0, F);
    psi := ApplyGroupAction3TI(B_0, B_0, B_0, phi_0, F);
    return phi, psi;
end function;

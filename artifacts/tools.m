function GenerateMonomials(vars, degree)
    monomials := [{}, {1}];
    for i in [1..degree] do
        monomials := monomials cat [{m * x : m in monomials[i+1], x in vars}];
    end for;
    return monomials;
end function;

function ApplyGroupAction3TI(A, B, C, phi, R)
    n, m, k := Explode(<Nrows(A), Nrows(B), Nrows(C)>);
    phi_m := [Matrix(R, m, k, [[phi[i][j][jj] : jj in [1..k]] : j in [1..m]]) : i in [1..n]];
    return [&+[A[i][ii] * (Transpose(B) * phi_m[i] * C) : i in [1..n]] : ii in [1..n]];
end function;

function Antisymmetrize(phi)
    return [[[phi[i][j][k] + phi[j][k][i] + phi[k][i][j] - phi[i][k][j] - phi[k][j][i] - phi[j][i][k] : i in [1..#phi]] : j in [1..#phi]] : k in [1..#phi]];
end function;

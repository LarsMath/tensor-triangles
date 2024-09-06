load "../tools.m";
load "./find_triangles_algorithm.m";
load "./post_collision_algorithm.m";

n := 13;
q := 2^32 - 5;

printf "Solve ATFE with a planted (but unknown) triangle (n: %o, q: %o)\n", n, q;

F := GF(q);

// =================================== Construct an instance with triangle =======================================

phi, psi := GenerateATFEKeyPairWithPlantedTriangle(F, n);

// ================================= Find both triangles ==================================================

while true do

    // Rerandomize
    A := RandomMatrix(F, n, n); while Rank(A) ne n do A := RandomMatrix(F, n, n); end while;
    B := RandomMatrix(F, n, n); while Rank(B) ne n do B := RandomMatrix(F, n, n); end while;

    phii := ApplyGroupAction3TI(A, A, A, phi, F);
    psii := ApplyGroupAction3TI(B, B, B, psi, F);

    solutions_phii := FindTriangles(phii);
    solutions_psii := FindTriangles(psii);

    if #solutions_phii ge 2 or #solutions_psii ge 2 then
        print "Found multiple triangles in the same ATF, currently not implemented stopping...";
        exit;
    end if;

    if #solutions_phii eq 1 and #solutions_psii eq 1 then break; end if;
end while;

// ===============================Extract triangles ==============================

x_phi := [1, 0, 0] cat [solutions_phii[1][i] : i in [1..n-3]];
y_phi := [0, 1, 0] cat [solutions_phii[1][i] : i in [n-2..2*n-6]];
z_phi := [0, 0, 1] cat [solutions_phii[1][i] : i in [2*n-5..3*n-9]];

x_psi := [1, 0, 0] cat [solutions_psii[1][i] : i in [1..n-3]];
y_psi := [0, 1, 0] cat [solutions_psii[1][i] : i in [n-2..2*n-6]];
z_psi := [0, 0, 1] cat [solutions_psii[1][i] : i in [2*n-5..3*n-9]];

// ===============================Transform instances===================================

phi_e, A_e := MoveTriangleGeneralPosition(phii, x_phi, y_phi, z_phi);
psi_e, B_e := MoveTriangleGeneralPosition(psii, x_psi, y_psi, z_psi);

// ========================= Post collision =================================

found := PostCollision(phi_e, psi_e );

assert #found ne 0;

s := found[1];
A_f := Matrix([(i le 3 select [s[j+3*(i-1)] : j in [1..3]] else [0,0,0]) cat [s[j + (n-3)*(i-1) + 9] : j in [1..n-3]] : i in [1..n]]);

A_tot := A * A_e * A_f * B_e^(-1) * B^(-1);

print "Found correct isometry:", ApplyGroupAction3TI(A_tot, A_tot, A_tot, phi, F) eq psi;

exit;

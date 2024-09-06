load "../tools.m";
load "./find_triangles_algorithm.m";
load "./post_collision_algorithm.m";

n := 10;
q := 29;

printf "Search for ATFE pair that has a triangle and solve it (n: %o, q: %o)\n", n, q;

F := GF(q);

// ================================= Find both triangles ==================================================

count := 0;
while true do
    count +:= 1;

    // Generate problem
    phi, psi := GenerateATFEKeyPair(F, n);

    // See if it contains a triangle
    solutions_phi := FindTriangles(phi);

    if #solutions_phi eq 1 then
        solutions_psi := FindTriangles(psi);
        if #solutions_psi eq 1 then
            break;
        end if;
    end if;
end while;

printf "Found an instance with a unique triangle in %o tries\n", count;

// ===============================Extract triangles ==============================

x_phi := [1, 0, 0] cat [solutions_phi[1][i] : i in [1..n-3]];
y_phi := [0, 1, 0] cat [solutions_phi[1][i] : i in [n-2..2*n-6]];
z_phi := [0, 0, 1] cat [solutions_phi[1][i] : i in [2*n-5..3*n-9]];

x_psi := [1, 0, 0] cat [solutions_psi[1][i] : i in [1..n-3]];
y_psi := [0, 1, 0] cat [solutions_psi[1][i] : i in [n-2..2*n-6]];
z_psi := [0, 0, 1] cat [solutions_psi[1][i] : i in [2*n-5..3*n-9]];

// ===============================Transform instances===================================

phi_e, A := MoveTriangleGeneralPosition(phi, x_phi, y_phi, z_phi);
psi_e, B := MoveTriangleGeneralPosition(psi, x_psi, y_psi, z_psi);

// ========================= Post collision =================================

found := PostCollision(phi_e, psi_e );

assert #found ne 0;

s := found[1];
A_f := Matrix([(i le 3 select [s[j+3*(i-1)] : j in [1..3]] else [0,0,0]) cat [s[j + (n-3)*(i-1) + 9] : j in [1..n-3]] : i in [1..n]]);

A_tot := A * A_f * B^(-1);

print "Found correct isometry:", ApplyGroupAction3TI(A_tot, A_tot, A_tot, phi, F) eq psi;

exit;


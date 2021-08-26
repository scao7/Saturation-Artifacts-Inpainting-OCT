function k = lamb2k_v2(A, Chirp)

cd1 = Chirp;

k = spline(cd1, A(1:1:2048), 1:2048);

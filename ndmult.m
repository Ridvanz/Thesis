function M = ndmult (A,B,dim)
  dA = dim(1);
  dB = dim(2);

  %reshape A into 2d
  sA = size (A);
  nA = length (sA);
  perA = [1:(dA-1) (dA+1):(nA-1) nA dA](1:nA);
  Ap = permute (A, perA);
  Ap = reshape (Ap, prod (sA(perA(1:end-1))), sA(perA(end)));

  %reshape B into 2d
  sB = size (B);
  nB = length (sB);
  perB = [dB 1:(dB-1) (dB+1):(nB-1) nB](1:nB);
  Bp = permute (B, perB);
  Bp = reshape (Bp, sB(perB(1)), prod (sB(perB(2:end))));

  %multiply
  M = Ap * Bp;

  % reshape back to original format
  s = [sA(perA(1:end-1)) sB(perB(2:end))];
  M = squeeze (reshape (M, s));
endfunction
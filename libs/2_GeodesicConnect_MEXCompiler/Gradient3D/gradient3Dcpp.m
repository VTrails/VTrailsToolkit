function [Ux,Uy,Uz] = gradient3Dcpp(U)

assert(isreal(U));

[Ux,Uy,Uz] = gradient(U);
# -----------------------------------------------------------------------------
# Context (from docs/main.tex): Massive scalar perturbations on Schwarzschild
# Matrices assembling routines for QNMs computations
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology,
#         Abu Dhabi, UAE)
# -----------------------------------------------------------------------------

MatrixAssembler := proc (
  d   ::integer, #    d : number of digits used in computations
  n   ::integer, #    n : number of Tchebyshev modes
  mu  ::numeric, #   mu : mass parameter
  L   ::numeric, #    L : angular momentum
  path::string   # path : string containing the path where we save the assembled matrices
  )
  local L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, L30::function, L31::function, L32::function, L40::function, L41::function, L42::function, L50::function, L51::function, L52::function, L60::function, L61::function, L62::function, L70::function, L71::function, L72::function, M0::Matrix, M1::Matrix, M2::Matrix, M3::Matrix, M4::Matrix, M5::Matrix, M6::Matrix, M7::Matrix, i::integer, j::integer, xi::numeric, path2::string, nstr::string, lambda::numeric, l00, l01, l02, l10, l11, l12, l20, l21, l22, l30, l31, l32, l40, l41, l42, l50, l51, l52, l60, l61, l62, l70, l71, l72, Tp, Tc, Tn, dTp, dTc, dTn, d2Tp, d2Tc, d2Tn:
  with(LinearAlgebra):
  Digits := d:
  # Definition of the differential operator coefficients Lij(y):
  # - L0j correspond to terms independent of the spectral parameter (e.g., Omega^0)
  # - i*L1j correspond to terms linear in the spectral parameter (Omega^1)
  # - L2j correspond to terms linear in the spectral parameter (Omega^2)
  # - i*L3j correspond to terms linear in the spectral parameter (Omega^3)
  # - L4j correspond to terms linear in the spectral parameter (Omega^4)
  # - i*L5j correspond to terms linear in the spectral parameter (Omega^5)
  # - L6j correspond to terms linear in the spectral parameter (Omega^6)
  # - i*L7j correspond to terms linear in the spectral parameter (Omega^7)
  # Mapping to the continuous form at a point y:
  #   L00*Phi + L01*Phi' + L02*Phi'' -> contributes to M0
  #   L10*Phi + L11*Phi' + L12*Phi'' -> contributes to M1
  #   L20*Phi + L21*Phi' + L22*Phi'' -> contributes to M2
  #   L30*Phi + L31*Phi' + L32*Phi'' -> contributes to M3
  #   L40*Phi + L41*Phi' + L42*Phi'' -> contributes to M4
  #   L50*Phi + L51*Phi' + L52*Phi'' -> contributes to M5
  #   L60*Phi + L61*Phi' + L62*Phi'' -> contributes to M6
  #   L70*Phi + L71*Phi' + L72*Phi'' -> contributes to M7
  # Note: In `docs/main.tex`, for equation (SE10) and (TSCH), specific S2,S1,S0
  #       yield particular Lij. The definitions below encode a specific model
  #       instance consistent with that framework.
  lambda := L*(L+1):
  L00 := y -> 0:
  L01 := y -> -8*mu^6:
  L02 := y -> 0:
  L10 := y -> mu^4*(1 + 2*lambda - y):
  L11 := y -> -mu^4*(3*y^2 - 2*y - 1):
  L12 := y -> -mu^4*(y + 1)*(1 - y)^2:
  L20 := y -> 4*mu^4:
  L21 := y -> 8*mu^4*(y + 2):
  L22 := y -> 0:
  L30 := y -> -2*mu^2*(1 + 2*lambda - y):
  L31 := y -> 2*mu^2*(3*y^2 - 2*y - 1):
  L32 := y -> 2*mu^2*(y + 1)*(1 - y)^2:
  L40 := y -> 4*mu^2*(y - 2):
  L41 := y -> 4*mu^2*(y^2 - 4*y - 3):
  L42 := y -> 0:
  L50 := y -> 16*mu^2 + 2*lambda - y + 1:
  L51 := y -> -3*y^2 + 2*y + 1:
  L52 := y -> -(y + 1)*(1 - y)^2:
  L60 := y -> 4*(1 - y):
  L61 := y -> -4*(y^2 - 2*y - 1):
  L62 := y -> 0:
  L70 := y -> 4*(-3 + y):
  L71 := y -> 0:
  L72 := y -> 0:
  # Assemble matrices using Chebyshev three-term recurrence (purely numerical).
  # M_k[i,j] = L_{k,0}(xi)*T_{j-1}(xi) + L_{k,1}(xi)*T'_{j-1}(xi) + L_{k,2}(xi)*T''_{j-1}(xi)
  # where T_j, T'_j, T''_j are computed via:
  #   T_0=1, T_1=x, T_{j+1} = 2x*T_j - T_{j-1}
  #   T'_0=0, T'_1=1, T'_{j+1} = 2*T_j + 2x*T'_j - T'_{j-1}
  #   T''_0=0, T''_1=0, T''_{j+1} = 4*T'_j + 2x*T''_j - T''_{j-1}
  M0 := Matrix(n):
  M1 := Matrix(n):
  M2 := Matrix(n):
  M3 := Matrix(n):
  M4 := Matrix(n):
  M5 := Matrix(n):
  M6 := Matrix(n):
  M7 := Matrix(n):
  for i from 1 to n do
    xi := evalf(cos((2.0*i-1.0)*Pi/(2.0*n))):
    # Evaluate operator coefficients at this collocation point:
    l00 := evalf(L00(xi)): l01 := evalf(L01(xi)): l02 := evalf(L02(xi)):
    l10 := evalf(L10(xi)): l11 := evalf(L11(xi)): l12 := evalf(L12(xi)):
    l20 := evalf(L20(xi)): l21 := evalf(L21(xi)): l22 := evalf(L22(xi)):
    l30 := evalf(L30(xi)): l31 := evalf(L31(xi)): l32 := evalf(L32(xi)):
    l40 := evalf(L40(xi)): l41 := evalf(L41(xi)): l42 := evalf(L42(xi)):
    l50 := evalf(L50(xi)): l51 := evalf(L51(xi)): l52 := evalf(L52(xi)):
    l60 := evalf(L60(xi)): l61 := evalf(L61(xi)): l62 := evalf(L62(xi)):
    l70 := evalf(L70(xi)): l71 := evalf(L71(xi)): l72 := evalf(L72(xi)):
    # Compute T_j(xi), T'_j(xi), T''_j(xi) and assemble row i:
    for j from 0 to n-1 do
      if j = 0 then
        Tc := evalf(1): dTc := evalf(0): d2Tc := evalf(0):
      elif j = 1 then
        Tp := Tc: dTp := dTc: d2Tp := d2Tc:
        Tc := xi: dTc := evalf(1): d2Tc := evalf(0):
      else
        Tn   := evalf(2*xi*Tc - Tp):
        dTn  := evalf(2*Tc + 2*xi*dTc - dTp):
        d2Tn := evalf(4*dTc + 2*xi*d2Tc - d2Tp):
        Tp := Tc: Tc := Tn:
        dTp := dTc: dTc := dTn:
        d2Tp := d2Tc: d2Tc := d2Tn:
      end if:
      M0[i,j+1] := l00*Tc + l01*dTc + l02*d2Tc:
      M1[i,j+1] := l10*Tc + l11*dTc + l12*d2Tc:
      M2[i,j+1] := l20*Tc + l21*dTc + l22*d2Tc:
      M3[i,j+1] := l30*Tc + l31*dTc + l32*d2Tc:
      M4[i,j+1] := l40*Tc + l41*dTc + l42*d2Tc:
      M5[i,j+1] := l50*Tc + l51*dTc + l52*d2Tc:
      M6[i,j+1] := l60*Tc + l61*dTc + l62*d2Tc:
      M7[i,j+1] := l70*Tc + l71*dTc + l72*d2Tc:
    end do:
  end do:
  # We finally export the data from Maple and save in files:
  path2 := cat(path, "/assemble/"):
  nstr := convert(n, string);
  ExportMatrix(cat(path2, "M0_", nstr, ".mat"), M0, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path2, "M1_", nstr, ".mat"), M1, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path2, "M2_", nstr, ".mat"), M2, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path2, "M3_", nstr, ".mat"), M3, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path2, "M4_", nstr, ".mat"), M4, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path2, "M5_", nstr, ".mat"), M5, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path2, "M6_", nstr, ".mat"), M6, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path2, "M7_", nstr, ".mat"), M7, target=MATLAB, mode=ascii):
end proc:
% Polynomial eigenvalue solver for massive scalar spectral roots on Schwarzschild
% Sheet-resolved output revision for DR13837.
%
% Authors: Davide Batic (Khalifa University of Science and Technology,
%          Abu Dhabi, UAE)
%          Anna Chrysostomou (LPTHE, Sorbonne Universite, CNRS, Paris, France)
%          Alan S. Cornell (University of Johannesburg, Auckland Park,
%          South Africa)
%          Dr. Denys Dutykh (Khalifa University of Science and Technology,
%          Abu Dhabi, UAE)
% Last modified: 2 August 2026
%
% The original routine mapped the uniformising eigenvalue Lambda to Omega and
% archived only Re(Omega), Im(Omega).  Because Lambda and mu^2/Lambda produce
% the same Omega but opposite k, that two-column output is insufficient for a
% posteriori sheet classification.  This version retains Lambda and k.

close
clear
format longE

addpath('/home/dds/Soft/advanpix/');

list = load('resolutions.txt')';
L = length(list);

% Read mu from shared config.
fid = fopen('params.txt', 'r');
mu = mp(strtrim(fgetl(fid)));
fclose(fid);

for idx = 1:L
  n = list(idx);
  nstr = num2str(n);

  fprintf('Computing eigs for n = %3d ... ', n);
  mp.Digits(n);

  M0 = mp.read(strcat('assemble/M0_', nstr, '.mat'));
  M1 = mp.read(strcat('assemble/M1_', nstr, '.mat'));
  M2 = mp.read(strcat('assemble/M2_', nstr, '.mat'));
  M3 = mp.read(strcat('assemble/M3_', nstr, '.mat'));
  M4 = mp.read(strcat('assemble/M4_', nstr, '.mat'));
  M5 = mp.read(strcat('assemble/M5_', nstr, '.mat'));
  M6 = mp.read(strcat('assemble/M6_', nstr, '.mat'));
  M7 = mp.read(strcat('assemble/M7_', nstr, '.mat'));

  % Seventh-degree polynomial eigenvalue problem in the uniformising variable.
  [V, Lambda, condnum] = polyeig(M0, mp('1i')*M1, M2, mp('1i')*M3, ...
                                 M4, mp('1i')*M5, M6, mp('1i')*M7); %#ok<ASGLU>

  Omega = (Lambda.^2 + mu^2)./(mp('2')*Lambda);
  kappa = (Lambda.^2 - mu^2)./(mp('2')*Lambda);

  % Keep Lambda, Omega, and k aligned when removing singular algebraic roots.
  good = isfinite(Lambda) & isfinite(Omega) & isfinite(kappa);
  Lambda = Lambda(good);
  Omega  = Omega(good);
  kappa  = kappa(good);

  % Preserve the original two-column file for backward compatibility.
  fname = strcat('results/eigs_', nstr, '.dat');
  fid = fopen(fname, 'w');
  for j = 1:length(Omega)
    fprintf(fid, '%s %s\n', num2str(real(Omega(j)), n), ...
                          num2str(imag(Omega(j)), n));
  end
  fclose(fid);

  % New sheet-resolved file.  The final two integer columns are:
  % sheet_flag   = +1 exterior/outgoing sheet (|Lambda|>mu),
  %                -1 reciprocal sheet (|Lambda|<mu), 0 on branch circle;
  % spatial_flag = +1 exponential decay (Im k>0),
  %                -1 exponential growth (Im k<0), 0 neutral.
  fname = strcat('results/eigs_sheet_resolved_', nstr, '.dat');
  fid = fopen(fname, 'w');
  fprintf(fid, ['# ReLambda ImLambda ReOmega ImOmega Rek Imk absLambda ' ...
                'sheet_flag spatial_flag\n']);

  for j = 1:length(Omega)
    sheet_flag = 0;
    if abs(Lambda(j)) > mu
      sheet_flag = 1;
    elseif abs(Lambda(j)) < mu
      sheet_flag = -1;
    end

    spatial_flag = 0;
    if imag(kappa(j)) > 0
      spatial_flag = 1;
    elseif imag(kappa(j)) < 0
      spatial_flag = -1;
    end

    fprintf(fid, '%s %s %s %s %s %s %s %d %d\n', ...
      num2str(real(Lambda(j)), n), num2str(imag(Lambda(j)), n), ...
      num2str(real(Omega(j)), n),  num2str(imag(Omega(j)), n),  ...
      num2str(real(kappa(j)), n),  num2str(imag(kappa(j)), n),  ...
      num2str(abs(Lambda(j)), n), sheet_flag, spatial_flag);
  end
  fclose(fid);

  fprintf('done.\n');
end

% Polynomial eigenvalue solver for massive scalar QNMs on Schwarzschild
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology,
%         Abu Dhabi, UAE)

close
clear
format longE

addpath('/home/dds/Soft/advanpix/');

list = load('resolutions.txt')';
L = length(list);

% Read mu from shared config:
fid = fopen('params.txt', 'r');
mu = mp(strtrim(fgetl(fid)));
fclose(fid);

for idx = 1:L
  n = list(idx);
  nstr = num2str(n);

  fprintf('Computing eigs for n = %3d ... ', n);

  mp.Digits(n);

  m0name = strcat('assemble/M0_', nstr, '.mat');
  m1name = strcat('assemble/M1_', nstr, '.mat');
  m2name = strcat('assemble/M2_', nstr, '.mat');
  m3name = strcat('assemble/M3_', nstr, '.mat');
  m4name = strcat('assemble/M4_', nstr, '.mat');
  m5name = strcat('assemble/M5_', nstr, '.mat');
  m6name = strcat('assemble/M6_', nstr, '.mat');
  m7name = strcat('assemble/M7_', nstr, '.mat');

  M0 = mp.read(m0name);
  M1 = mp.read(m1name);
  M2 = mp.read(m2name);
  M3 = mp.read(m3name);
  M4 = mp.read(m4name);
  M5 = mp.read(m5name);
  M6 = mp.read(m6name);
  M7 = mp.read(m7name);

  % Original polynomial eigenvalue formulation:
  [V, Lambda, k] = polyeig(M0, mp('1i')*M1, M2, mp('1i')*M3, M4, mp('1i')*M5, M6, mp('1i')*M7);

  e = (Lambda.^2 + mu^2)./(mp('2')*Lambda);

  % Discard spurious modes (Inf/NaN from zero or degenerate Lambda):
  e = e(isfinite(e));

  % Save QNMs to platform-agnostic text file (full multiprecision):
  fname = strcat('results/eigs_', nstr, '.dat');
  fid = fopen(fname, 'w');
  for k = 1:length(e)
    fprintf(fid, '%s %s\n', num2str(real(e(k)), n), num2str(imag(e(k)), n));
  end
  fclose(fid);

  fprintf('done.\n');

end % for

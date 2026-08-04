clear all
clc
close all

%% Analytical and Numerical Soln
[C9A, dy9] = Analytical(9);
[C17A, dy17] = Analytical(17);
[C33A, dy33] = Analytical(33);
[C65A, dy65] = Analytical(65);
[C129A, dy129] = Analytical(129);

[C9N] = NumericalFTCS(9);
[C17N] = NumericalFTCS(17);
[C33N] = NumericalFTCS(33);
[C65N] = NumericalFTCS(65);
[C129N] = NumericalFTCS(129);

%% Error Norm

% nx = 9 row
L1_9 = norm(C9N - C9A, 1)*dy9;
L2_9 = norm(C9N - C9A, 2)*sqrt(dy9);
LInf_9 = norm(C9N - C9A, Inf);

% nx = 17 row
L1_17 = norm(C17N - C17A, 1)*dy17;
L2_17 = norm(C17N - C17A, 2)*sqrt(dy17);
LInf_17 = norm(C17N - C17A, Inf);

% nx = 33 row
L1_33 = norm(C33N - C33A, 1)*dy33;
L2_33 = norm(C33N - C33A, 2)*sqrt(dy33);
LInf_33 = norm(C33N - C33A, Inf);

% nx = 65 row
L1_65 = norm(C65N - C65A, 1)*dy65;
L2_65 = norm(C65N - C65A, 2)*sqrt(dy65);
LInf_65 = norm(C65N - C65A, Inf);

% nx = 129 row
L1_129 = norm(C129N - C129A, 1)*dy129;
L2_129 = norm(C129N - C129A, 2)*sqrt(dy129);
LInf_129 = norm(C129N - C129A, Inf);

%% Order of Accuracy

% n(L1) column
nL1_17 = (log10(L1_9/L1_17))/(log10(dy9/dy17));
nL1_33 = (log10(L1_17/L1_33))/(log10(dy17/dy33));
nL1_65 = (log10(L1_33/L1_65))/(log10(dy33/dy65));
nL1_129 = (log10(L1_65/L1_129))/(log10(dy65/dy129));

% n(L2) column
nL2_17 = (log10(L2_9/L2_17))/(log10(dy9/dy17));
nL2_33 = (log10(L2_17/L2_33))/(log10(dy17/dy33));
nL2_65 = (log10(L2_33/L2_65))/(log10(dy33/dy65));
nL2_129 = (log10(L2_65/L2_129))/(log10(dy65/dy129));

% n(LInf) column
nLInf_17 = (log10(LInf_9/LInf_17))/(log10(dy9/dy17));
nLInf_33 = (log10(LInf_17/LInf_33))/(log10(dy17/dy33));
nLInf_65 = (log10(LInf_33/LInf_65))/(log10(dy33/dy65));
nLInf_129 = (log10(LInf_65/LInf_129))/(log10(dy65/dy129));
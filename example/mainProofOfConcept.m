%% Implement spectral + energy/independent bounds
%% Generate data: circulant covariance matrices
%We only need to define first and second moments for each column.
M=30;
N=1e4;

switch k
    case 1
%Case 1: same ESD profile, same energy for all columns

    case 2
%Case 2: same ESD profile, different energy

    case 3
%Case 3: different ESD profile, same energy

    case 4
%Case 4: different ESD profile, different energy (!)

    case 5 %same covariance matrix, same energy(?) but non-circulant matrix. Still using DFT ESD for bound

end

%% First, assume we know the ESD of all columns:

bound1=VAF(F); %Bound 1 through DFT ESD. This is the best possible for circulant covariance matrices
bound1Alt=VAF(eig(Sigma)); %For the case of non-circulant matrix, we can compute the bound
%from the actual eigenvalues of the TRUE covariance matrix, which should
%give a tighter bound. Sigma is the sum of all the NxN covariance matrices

%% Second, assuming independence, get a bound:
bound2=

%% Third, merge bounds:
mixedBound=mergeBounds(bound1,bound2)

%% Fourth, under some independence between bounds assumption, we can do a heuristic improvement:
E=tr(Sigma); %Total expected energy
bestBound=bound1+(E-bound1).*(bound2/E);
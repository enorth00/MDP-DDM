# MDP-DDM
MatLab implementation of a DDM framework built on top of the Method of Difference Potentials.

## Requirements to run:
- Download the package 'chebfun' from the following website: https://www.chebfun.org/download/
- Add 'chebfun' to your matlab path. I added mine permanently, but it can be done on a session-by-session basis.

## Summary:
Main contributing code for Evan North's dissertation, found at: https://repository.lib.ncsu.edu/handle/1840.20/39942

This repository implements a solver for the Helmholtz equation over rectangular domains. These rectangular domains are constructed from square subdomains which are allowed to vary in material properties (i.e. Helmholtz wavenumber). In addition to the uniform square subdomains, "rod" subdomains have been implemented that support a scattering rod of arbitrary radius 'r' centered in the square subdomain. Together, these subdomains allow for the simulation of our main application of interest, the Photonic Crystal Ring Resonator (PCRR). 

## How to Use
TEST_exact.m contains the complete workflow for a working example, and the process should be able to be derived from there. Here's a layout though:
1) Pick a domain/init case and run that script.
2) Set the following parameters:
    a) OPTIONS.GRIDS = each n you want so that a subdomain's discretization has 2^n + 1 nodes.
    b) OPTIONS.FUNC = which test function from u_exact.m you want to use.
    c) OPTIONS.PARALLEL = true/false, are you running in parallel
    d) OPTIONS.REAL = true/false, is your solution real or complex? If real, this saves some computation.
    e) OPTIONS.RADIATION = absorbing boundary condition? 0 for no, 1-3 for implemented Engquist-Majda ABCs.
    f) OPTIONS.FOURIER_AUX = If fourier functions are used along edges, these are their intervals. Not a recommended method, and this shouldn't affect anything when not using them.
3) Set subdomains info with SET_SubdomainInfo
4) Set discretization grids
5) Run the problem_solver() workhorse
6) Compare computed solution to exact solution

## Explanation of inits
A lot of critical parameters are set in the init files. Most of them should be commented, but this should be a master collection of their meanings as well:
- Mc is the number of chebyshev functions (dimension of the expansion) of each side of a subdomain.
- Mf is the number of fourier functions per side (if using the fourier expansion option, not recommended).
- Mr is the number of fourier functions for the rod.
- k1 is the wavenumber of the square itself.
- k2 is the wavenumber of the rod. If there's no rod, this should be set to zero.
- TYPES identifies each unique block you want to use in construction of the domain. index (i,:) indicates all of the info for the ith block. In order, the parameters are [{type}, {basis for functions along each edge}, {[k1, k2]}]
  - type is either 1=empty or 2=rod
  - basis functions are a string of 'c' or 'f' characters. 'cccc' is Chebyshev along all four edges of the square. 'ffcc' is Fourier functions on the left and right edges, and Chebyshev along the top and bottom. Edge order is left-right-bottom-top.
  - Wavenumber is always specified with two values. If the type is 1 (empty, no rod), it is good practice to explicitly set this to zero rather than setting k2 to zero and using it here.
- DOMAIN is a matrix that provides a visual representation of how the domain is constructed from the above defined blocks in TYPES. Use the indices of the defined blocks to define where they sit.
  - Example: Let TYPES(1) be an empty subdomain and TYPES(2) be a basic rod. Then DOMAIN=[1, 1; 1, 2] creates a 2x2 domain where only the bottom right one contains a rod.
- start is the coordinates of the left-right-bottom-top sides of the first subdomain (i.e. DOMAIN(1,1)). So start=[-3, -1, 5, 7] means that DOMAIN(1,1) spans x=-3 to -1 and y=5 to 7.

Start by looking at these set examples. From there, it should be clear what the other examples are creating/accomplishing.
- empty_base.m: Single subdomain with a uniform wavenumber (no rod). 
- rod_base.m: Single subdomain with a scattering rod (two wavenumbers).
- corner_scatter_3x3.m: Bottom-left subdomain is a scatterer, the rest are empty.
- tunnel_3x3.m: 3 uniform subdomains (horizontal) with a row of three rod subdomains above and below (makes a 3x3 square).

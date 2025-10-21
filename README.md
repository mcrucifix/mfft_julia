# mfft_julia

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Embryonic Julia implementation of the Modified Fourier Transform (MFFT).

## Overview

The Modified Fourier Transform (MFFT) was, to the best of the author’s knowledge, first introduced in **Laskar (1988)**, Sect. 3, p. 343.  
The present implementation follows the approach described by **Šidlichovský and Nesvorný (1997)**.

Further developments are available in the [**gtseries**](https://github.com/mcrucifix/gtseries) package, where the MFFT is extended to real numbers and includes the frequency correction algorithm of Šidlichovský & Nesvorný (1997), effectively upgrading Laskar’s MFFT into a *Frequency Modified Fourier Transform (FMFT)*.

## Code provenance

The original C implementation by Šidlichovský and Nesvorný can be found at:  
[https://www2.boulder.swri.edu/~davidn/fmft/fmft.c](https://www2.boulder.swri.edu/~davidn/fmft/fmft.c)

This file is made publicly available by the authors, but no explicit license is provided. Consequently, its direct reuse or redistribution is not permitted.  
The present Julia code does **not** include any portion of that C code; it is a clean reimplementation based on the algorithmic description provided in their publication.  
Best efforts have been made to contact the authors regarding licensing. All algorithmic steps are fully documented in their paper, and the reimplementation could have been produced without consulting the source.

## Disclaimer

Code is provided *as is*, as part of an ongoing research project.

For a stable, documented, and fully tested version in R, see:  
[https://github.com/mcrucifix/gtseries](https://github.com/mcrucifix/gtseries)

## References

- Laskar, J. (1988). *Secular evolution of the solar system over 10 million years*.  
  **Astronomy and Astrophysics**, 198, 341–362.  
- Šidlichovský, M., & Nesvorný, D. (1997). *Frequency Modified Fourier Transform and its Application to Asteroids.*  
  In **The Dynamical Behaviour of our Planetary System** (pp. 138–148). Springer.  
  doi:[10.1007/978-94-011-5510-6_9](https://doi.org/10.1007/978-94-011-5510-6_9)

## License

This project is licensed under the [MIT License](LICENSE).

## Sample code

```julia
using Statistics
include("mfft.jl")

V = 0.02 * randn(1024) + im* 0.0004 * randn(1024) .+ 
             3.0.*(cis(0.142 * t) for t in (0:1023)) .+
             (2.93 + 0.14im) .*(cis(0.192 * t) for t in (0:1023)) .+
             (0.93 + 0.54im) .*(cis(0.242 * t) for t in (0:1023)) .+
             (0.55 + 0.93im) .*(cis(0.441 * t) for t in (0:1023)) 

Amp, freq = mfft(V)
```

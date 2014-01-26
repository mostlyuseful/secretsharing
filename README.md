Secret sharing in Python
========================

Secret sharing via Adi Shamir's method.

There are three methods:

  - Sharing an arithmetic secret (a natural number in the set N_0) resulting in big integer coefficients
  - Sharing a secret (a natural number in the set N_0) in the finite field resulting in big integer coefficients
  - Sharing a secret string split into smaller chunks resulting in coefficients of controllable sizes

To reconstruct the underlying formula Lagrange polynomials are used in an O^2 implementation.

This implementation is mostly done to easily follow the general idea, it is not optimized for speed.

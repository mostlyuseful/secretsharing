# -*- coding: utf-8 -*-

# Copyright (c) 2014, Maurice-Pasca
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Maurice-Pascal Sonnemann nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL MAURICE-PASCAL SONNEMANN BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import random

def error_probability(k):
    '''Returns the error probability that a composite is detected as a prime
    '''
    return 1./(4**k)

def forever(start = 0, step=1):
    i = start
    while True:
        yield i
        i += step

def is_probable_prime(n, k=10):
    '''Returns True if n is probably prime
    n : Number to test
    k : Count of tests, the higher k, the lower the probability of an error in
        determining the primality of n
        
    Returns True if n is a probable prime, False if it is definitely composite
    
    Results are provably correct for n < 341,550,071,728,321.
    Resultes for n >= 341,550,071,728,321 have a error probability p_e(k)=1/(4**k), with k=10 being the default: p_e(10) ≅ 95.367e-06 percent
        
    Applies the Miller-Rabin test, using some carefully chosen bases, see http://primes.utm.edu/prove/prove2_3.html
    '''
    
    def factor2(n):
        d = n-1
        for s in forever():
            D = d >> s
            if D % 2:
                return D,s
    
    assert n >= 0
    if n < 4:
        #  n =      0,     1,    2,     3
        return (False, False, True, True)[n]
    if n % 2 == 0:
        return False
    
    d, s = factor2(n)
    assert(2**s*d == n-1)
    
    def is_maybe_composite(a):
        if pow(a, d, n) == 1:
            return False
        for i in range(s):
            if pow(a, 2**i * d, n) == n-1:
                return False
        return True # n definitely composite
    
    # Try exact proofs before fail-over to probabilistic test
    
    def definitely_prime(a_list):
        return not any(is_maybe_composite(a) for a in a_list)
    
    table = [(1373653, (2,3)),
             (9080191, (31,73)),
             (25326001,(2,3,5)),
             (4759123141,(2,7,61))]
    for limit, a_list in table:
        if n < limit:
            return definitely_prime(a_list)

    if n < 118670087467:
        if n == 3215031751:
            return False
        else:
            return not any(is_maybe_composite(a) for a in (2,3,5,7))
    
    table = [(2152302898747, (2,3,5,7,11)),
             (3474749660383, (2,3,5,7,11,13)),
             (341550071728321,(2,3,5,7,11,13,17))]
    for limit, a_list in table:
        if n < limit:
            return definitely_prime(a_list)
    
    # Probabilistic test follows here:
    for i in range(k):
        a = random.randrange(2,n)
        if is_maybe_composite(a):
            return False
    return True

PRIME_TABLE = [2, 3] + [x for x in range(5, 1000, 2) if is_probable_prime(x)]


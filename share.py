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

from __future__ import division
from miller_rabin import is_probable_prime, forever
import numpy as np
import random
import gmpy

# Helper functions

def generate_prime_bigger_than(x):
    if (x%2)==0:
        start = x+1
    else:
        start = x+2
    for t in forever(start=start, step=2):
        if is_probable_prime(t):
            return t
        
def str_to_dec(s):
    d = 0
    for i,c in enumerate(s):
        d += ord(c) * 256**(len(s)-i-1)
    return d

def dec_to_str(x):
    l = []
    running = gmpy.mpq(x)
    i = 100
    
    if x > 0:
        i = int(np.ceil(np.log10(float(x))/np.log10(256)))+1
    else:
        return ''
    while x//256**i == 0:
        i -= 1
    
    while running > 0:
        q, r = divmod(running, 256**i)
        l.append(chr(int(q)))
        running = r
        i -= 1
    return "".join(l)

class ArithmeticLagrangePolynomial(object):
    '''Slow! Just an ordinary O² implementation!
    See here for more intelligent approaches: 
        Aho, Hopcroft and Ullman: The Design and Analysis of Computer Algorithms, Addison-Wesley, 1974
        Knuth: TAOCP Vol. 2, 1969, page 503ff.
        http://en.wikipedia.org/wiki/Lagrange_polynomial#Barycentric_interpolation
    '''
    
    def __init__(self, xx, yy):
        self.xx = [ gmpy.mpq(x) for x in xx ]
        self.yy = [ gmpy.mpq(y) for y in yy ]
        assert len(self.xx) == len(self.yy)
        
    
    def _base(self, j, x):
        l = 1
        k = len(self.xx)
        x = gmpy.mpq(x)
        for m in range(k):
            if j==m:
                continue
            l *= (x-self.xx[m])/(self.xx[j]-self.xx[m])
        return l
    
    def _at(self, x):
        k = len(self.xx)
        return sum( self.yy[j]*self._base(j,x) for j in range(k) )
    
    def __call__(self, x):
        if hasattr(x, '__iter__'):
            return [self(v) for v in x]
        else:
            return self._at(x)

class ArithmeticSecret(object):
    
    @classmethod
    def polyval(cls, a, x):
        y = gmpy.mpq(0)
        for i, a_i in enumerate(reversed(a)):
            y += a_i*x**i
        return y
    
    @classmethod
    def share(cls, D, k, n):
        # Prepare polynomial
        a = [ random.randint(-D,D) for _ in range(k-1) ] + [D]
        for i in range(n):
            x = i+1
            D_i = gmpy.mpq(cls.polyval(a, x))
            if D_i.denominator == 1:
                D_i = gmpy.mpz(D_i)
            yield (x, D_i)
        
    
    @classmethod
    def join(cls, D_seq):
        D_seq = list(D_seq)
        xx = [D_i[0] for D_i in D_seq]
        yy = [D_i[1] for D_i in D_seq]
        if len(set(xx)) < len(xx):
            raise RuntimeError("D_i's must be unique!")
        L = ArithmeticLagrangePolynomial(xx, yy)
        D = L(0)
        if D.denominator == 1:
            return gmpy.mpz(D)
        else:
            return D
    
class FiniteFieldLagrangePolynomial(object):
    '''Slow! Just an ordinary O² implementation!
    See here for more intelligent approaches: 
        Aho, Hopcroft and Ullman: The Design and Analysis of Computer Algorithms, Addison-Wesley, 1974
        Knuth: TAOCP Vol. 2, 1969, page 503ff.
        http://en.wikipedia.org/wiki/Lagrange_polynomial#Barycentric_interpolation
    '''
    
    def __init__(self, xx, yy, P):
        self.xx = [ gmpy.mpq(x) for x in xx ]
        self.yy = [ gmpy.mpq(y) for y in yy ]
        self.P = P
        assert len(self.xx) == len(self.yy)
        
    
    def _base(self, j, x):
        l = 1
        k = len(self.xx)
        x = gmpy.mpq(x)
        for m in range(k):
            if j==m:
                continue
            a = (x-self.xx[m])
            b = (self.xx[j]-self.xx[m])
            l *= (a/b)
        return l
    
    def _at(self, x):
        k = len(self.xx)
        z = gmpy.mpq(0)
        for j in range(k):
            y = self.yy[j]
            b = self._base(j,x)
            z += (y*b)
            z %= self.P
        return z
    
    def __call__(self, x):
        if hasattr(x, '__iter__'):
            return [self(v) for v in x]
        else:
            return self._at(x)
        
class FiniteFieldSecret(object):
    
    @classmethod
    def polyval(cls, a, p, x):
        y = gmpy.mpq(0)
        for i, a_i in enumerate(reversed(a)):
            y += (a_i * ((x**i) % p)) % p
            y %= p
        return y
    
    @classmethod
    def share(cls, D, k, n, P=None):
        if P is None or not is_probable_prime(P):            
            P = generate_prime_bigger_than(random.randint(max(D, n),100000000*max(D, n)))
        assert P > D
        assert P > n
        assert D >= 0
        if D > 0:
            a = [ random.randint(0, D-1) for _ in range(k-1) ] + [D]
        else:
            a = [ random.randint(0, 1000) for _ in range(k-1) ] + [D]
        shares = []
        for i in range(n):
            x = i+1
            D_i = gmpy.mpq(cls.polyval(a, P, x))
            if D_i.denominator == 1:
                D_i = int(gmpy.mpz(D_i))
            shares.append( (x, D_i) )
        return (P, shares)
    
    @classmethod
    def join(cls, P, D_seq):
        D_seq = list(D_seq)
        xx = [D_i[0] for D_i in D_seq]
        yy = [D_i[1] for D_i in D_seq]
        if len(set(xx)) < len(xx):
            raise RuntimeError("Shares must be unique!")
        L = FiniteFieldLagrangePolynomial(xx, yy, P)
        D = L(0)
        if D.denominator == 1:
            return int(gmpy.mpz(D))
        else:
            return D

class FiniteFieldByteSecret(object):
    
    @classmethod
    def share_small(cls, msg_bytes, k, n, share_length, P=None):
        D_max = 256**share_length
        if P is None or not is_probable_prime(P):            
            P = generate_prime_bigger_than(random.randint(max(D_max, n),100000000*max(D_max, n)))
        assert P > D_max
        assert P > n
        assert 1 <= share_length <= len(msg_bytes)
        all_shares = {}
        for start_idx in range(0, len(msg_bytes), share_length):
            msg_part = msg_bytes[start_idx:start_idx+share_length]
            if not msg_part:
                continue
            D_part = str_to_dec(msg_part)
            shares_part = FiniteFieldSecret.share(D_part, k, n, P=P)[1]
            for i,share in enumerate(shares_part):
                l = all_shares.get(i, [])
                l.append(share)
                all_shares[i] = l
        return share_length, P, len(msg_bytes), all_shares.values()
    
    @classmethod
    def join_small(cls, share_length, P, msg_len, shares):
        
        def pad(s, min_len):
            if len(s) >= min_len:
                return s
            return s+(min_len-len(s))*'\0'
        
        msg = []
        sub_share_count = len(shares[0])
        participant_count = len(shares)
        for share_idx in range(sub_share_count):
            sub_shares = [ shares[participant][share_idx] for participant in range(participant_count) ]
            D_sub = FiniteFieldSecret.join(P, sub_shares)
            msg_sub = pad(dec_to_str(D_sub), share_length)
            msg.append(msg_sub)
        return (''.join(msg))[:msg_len]

if __name__=='__main__':
    msg = "Hello, world! Everything nice and wobbly?"
    assert ArithmeticSecret.join( ArithmeticSecret.share(str_to_dec(msg),3,60)) == str_to_dec(msg)

    P, shares = FiniteFieldSecret.share(str_to_dec(msg),3,60)
    assert FiniteFieldSecret.join(P, shares) == str_to_dec(msg)
    
    share_len, P2, msg_len, shares2 = FiniteFieldByteSecret.share_small(msg, 3, 5, 2)
    msg_decoded = FiniteFieldByteSecret.join_small(share_len, P2, msg_len, shares2)
    print repr(msg_decoded)
    assert msg_decoded == msg
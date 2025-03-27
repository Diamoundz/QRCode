import array
import math
try:
    bytearray
    _bytearray = bytearray
except NameError:
    def _bytearray(obj = 0, encoding = "latin-1"):
        '''Simple pure-python bytearray replacement if not implemented'''
        if isinstance(obj, str):
            obj = [ord(ch) for ch in obj.encode(encoding)]
        elif isinstance(obj, int):
            obj = [0] * obj
        return array.array("B", obj)
try:
    xrange
except NameError:
    xrange = range
class ReedSolomonError(Exception):
    pass
gf_exp = _bytearray([1] * 512)
gf_log = _bytearray(256)
field_charac = int(2**8 - 1)
def rwh_primes1(n):
    ''' Returns a list of primes < n '''
    n_half = int(n/2)
    sieve = [True] * n_half
    for i in xrange(3,int(math.pow(n,0.5))+1,2):
        if sieve[int(i/2)]:
            sieve[int((i*i)/2)::i] = [False] * int((n-i*i-1)/(2*i)+1)
    return array.array('i', [2] + [2*i+1 for i in xrange(1,n_half) if sieve[i]])
def find_prime_polys(generator=2, c_exp=8, fast_primes=False, single=False):
    '''Compute the list of prime polynomials for the given generator and galois field characteristic exponent.'''
    root_charac = 2
    field_charac = int(root_charac**c_exp - 1)
    field_charac_next = int(root_charac**(c_exp+1) - 1)
    if fast_primes:
        prim_candidates = rwh_primes1(field_charac_next)
        prim_candidates = array.array('i', [x for x in prim_candidates if x > field_charac])
    else:
        prim_candidates = array.array('i', xrange(field_charac+2, field_charac_next, root_charac))
    correct_primes = array.array('i', [])
    for prim in prim_candidates:
        seen = _bytearray(field_charac+1)
        conflict = False
        x = 1
        for i in xrange(field_charac):
            x = gf_mult_noLUT(x, generator, prim, field_charac+1)
            if x > field_charac or seen[x] == 1:
                conflict = True
                break
            else:
                seen[x] = 1
        if not conflict:
            correct_primes.append(prim)
            if single: return array.array('i', [prim])
    return correct_primes
def init_tables(prim=0x11d, generator=2, c_exp=8):
    '''Precompute the logarithm and anti-log tables for faster computation later, using the provided primitive polynomial.
    These tables are used for multiplication/division since addition/substraction are simple XOR operations inside GF of characteristic 2.
    The basic idea is quite simple: since b**(log_b(x), log_b(y)) == x * y given any number b (the base or generator of the logarithm), then we can use any number b to precompute logarithm and anti-log (exponentiation) tables to use for multiplying two numbers x and y.
    That's why when we use a different base/generator number, the log and anti-log tables are drastically different, but the resulting computations are the same given any such tables.
    For more infos, see https://en.wikipedia.org/wiki/Finite_field_arithmetic
    '''
    global _bytearray
    if c_exp <= 8:
        _bytearray = bytearray
    else:
        def _bytearray(obj = 0, encoding = "latin-1"):
            '''Fake bytearray replacement, supporting int values above 255'''
            if isinstance(obj, str):
                obj = obj.encode(encoding)
                if isinstance(obj, str):
                    obj = [ord(chr) for chr in obj]
                elif isinstance(obj, bytes):
                    obj = [int(chr) for chr in obj]
                else:
                    raise(ValueError, "Type of object not recognized!")
            elif isinstance(obj, int):
                obj = [0] * obj
            elif isinstance(obj, bytes):
                obj = [int(b) for b in obj]
            return array.array("i", obj)
    global gf_exp, gf_log, field_charac
    field_charac = int(2**c_exp - 1)
    gf_exp = _bytearray(field_charac * 2)
    gf_log = _bytearray(field_charac+1)
    x = 1
    for i in xrange(field_charac):
        gf_exp[i] = x
        gf_log[x] = i
        x = gf_mult_noLUT(x, generator, prim, field_charac+1)
    for i in xrange(field_charac, field_charac * 2):
        gf_exp[i] = gf_exp[i - field_charac]
    return gf_log, gf_exp, field_charac
def gf_add(x, y):
    '''Add two galois field integers'''
    return x ^ y
def gf_sub(x, y):
    '''Subtract two galois field integers'''
    return x ^ y
def gf_neg(x):
    '''Negate one galois field integer (does nothing)'''
    return x
def gf_inverse(x):
    '''Inverse of a galois field integer'''
    return gf_exp[field_charac - gf_log[x]]
def gf_mul(x, y):
    '''Multiply two galois field integers'''
    if x == 0 or y == 0:
        return 0
    return gf_exp[(gf_log[x] + gf_log[y]) % field_charac]
def gf_div(x, y):
    '''Divide x by y galois field integers'''
    if y == 0:
        raise ZeroDivisionError()
    if x == 0:
        return 0
    return gf_exp[(gf_log[x] + field_charac - gf_log[y]) % field_charac]
def gf_pow(x, power):
    '''Power of x galois field integer'''
    return gf_exp[(gf_log[x] * power) % field_charac]
def gf_mult_noLUT_slow(x, y, prim=0):
    '''Multiplication in Galois Fields on-the-fly without using a precomputed look-up table (and thus it's slower) by using the standard carry-less multiplication + modular reduction using an irreducible prime polynomial.'''
    def cl_mult(x,y):
        '''Bitwise carry-less multiplication on integers'''
        z = 0
        i = 0
        while (y>>i) > 0:
            if y & (1<<i):
                z ^= x<<i
            i += 1
        return z
    def bit_length(n):
        '''Compute the position of the most significant bit (1) of an integer. Equivalent to int.bit_length()'''
        bits = 0
        while n >> bits: bits += 1
        return bits
    def cl_div(dividend, divisor=None):
        '''Bitwise carry-less long division on integers and returns the remainder'''
        dl1 = bit_length(dividend)
        dl2 = bit_length(divisor)
        if dl1 < dl2:
            return dividend
        for i in xrange(dl1-dl2,-1,-1):
            if dividend & (1 << i+dl2-1):
                dividend ^= divisor << i
        return dividend
    result = cl_mult(x,y)
    if prim > 0:
        result = cl_div(result, prim)
    return result
def gf_mult_noLUT(x, y, prim=0, field_charac_full=256, carryless=True):
    '''Galois Field integer multiplication on-the-fly without using a look-up table, using Russian Peasant Multiplication algorithm (faster than the standard multiplication + modular reduction). This is still slower than using a look-up table, but is the fastest alternative, and is often used in embedded circuits where storage space is limited (ie, no space for a look-up table).
    If prim is 0 and carryless=False, then the function produces the result for a standard integers multiplication (no carry-less arithmetics nor modular reduction).'''
    r = 0
    while y:
        if y & 1: r = r ^ x if carryless else r + x
        y = y >> 1
        x = x << 1
        if prim > 0 and x & field_charac_full: x = x ^ prim
    return r
def gf_poly_scale(p, x):
    '''Scale a galois field polynomial with a factor x (an integer)'''
    out = _bytearray(len(p))
    for i in range(len(p)):
        out[i] = gf_mul(p[i], x)
    return out
def gf_poly_add(p, q):
    '''Add two galois field polynomials'''
    q_len = len(q)
    r = _bytearray( max(len(p), q_len) )
    r[len(r)-len(p):len(r)] = p
    for i in xrange(q_len):
        r[i + len(r) - q_len] ^= q[i]
    return r
def gf_poly_mul(p, q):
    '''Multiply two polynomials, inside Galois Field (but the procedure is generic). Optimized function by precomputation of log.'''
    r = _bytearray(len(p) + len(q) - 1)
    lp = [gf_log[p[i]] for i in xrange(len(p))]
    for j in xrange(len(q)):
        qj = q[j]
        if qj != 0:
            lq = gf_log[qj]
            for i in xrange(len(p)):
                if p[i] != 0:
                    r[i + j] ^= gf_exp[lp[i] + lq]
    return r
def gf_poly_mul_simple(p, q):
    '''Multiply two polynomials, inside Galois Field'''
    r = _bytearray(len(p) + len(q) - 1)
    for j in xrange(len(q)):
        for i in xrange(len(p)):
            r[i + j] ^= gf_mul(p[i], q[j])
    return r
def gf_poly_neg(poly):
    '''Returns the polynomial with all coefficients negated. In GF(2^p), negation does not change the coefficient, so we return the polynomial as-is.'''
    return poly
def gf_poly_div(dividend, divisor):
    '''Fast polynomial division by using Extended Synthetic Division and optimized for GF(2^p) computations (doesn't work with standard polynomials outside of this galois field).'''
    msg_out = _bytearray(dividend)
    divisor_len = len(divisor)
    for i in xrange(len(dividend) - (divisor_len-1)):
        coef = msg_out[i]
        if coef != 0:
            for j in xrange(1, divisor_len):
                if divisor[j] != 0:
                    msg_out[i + j] ^= gf_mul(divisor[j], coef)
    separator = -(divisor_len-1)
    return msg_out[:separator], msg_out[separator:]
def gf_poly_square(poly):
    '''Linear time implementation of polynomial squaring. For details, see paper: "A fast software implementation for arithmetic operations in GF (2n)". De Win, E., Bosselaers, A., Vandenberghe, S., De Gersem, P., & Vandewalle, J. (1996, January). In Advances in Cryptology - Asiacrypt'96 (pp. 65-76). Springer Berlin Heidelberg.'''
    length = len(poly)
    out = _bytearray(2*length - 1)
    for i in xrange(length-1):
        p = poly[i]
        k = 2*i
        if p != 0:
            out[k] = gf_exp[2*gf_log[p]]
    out[2*length-2] = gf_exp[2*gf_log[poly[length-1]]]
    if out[0] == 0: out[0] = 2*poly[1] - 1
    return out
def gf_poly_eval(poly, x):
    '''Evaluates a polynomial in GF(2^p) given the value for x. This is based on Horner's scheme for maximum efficiency.'''
    y = poly[0]
    for i in xrange(1, len(poly)):
        y = gf_mul(y, x) ^ poly[i]
    return y
def rs_generator_poly(nsym, fcr=0, generator=2):
    '''Generate an irreducible generator polynomial (necessary to encode a message into Reed-Solomon)'''
    g = _bytearray([1])
    for i in xrange(nsym):
        g = gf_poly_mul(g, [1, gf_pow(generator, i+fcr)])
    return g
def rs_generator_poly_all(max_nsym, fcr=0, generator=2):
    '''Generate all irreducible generator polynomials up to max_nsym (usually you can use n, the length of the message+ecc). Very useful to reduce processing time if you want to encode using variable schemes and nsym rates.'''
    g_all = [[1]] * max_nsym
    for nsym in xrange(max_nsym):
        g_all[nsym] = rs_generator_poly(nsym, fcr, generator)
    return g_all
def rs_simple_encode_msg(msg_in, nsym, fcr=0, generator=2, gen=None):
    '''Simple Reed-Solomon encoding (mainly an example for you to understand how it works, because it's slower than the inlined function below)'''
    global field_charac
    if (len(msg_in) + nsym) > field_charac: raise ValueError("Message is too long (%i when max is %i)" % (len(msg_in)+nsym, field_charac))
    if gen is None: gen = rs_generator_poly(nsym, fcr, generator)
    _, remainder = gf_poly_div(msg_in + _bytearray(len(gen)-1), gen)
    msg_out = msg_in + remainder
    return msg_out
def rs_encode_msg(msg_in, nsym, fcr=0, generator=2, gen=None):
    '''Reed-Solomon main encoding function, using polynomial division (Extended Synthetic Division, the fastest algorithm available to my knowledge), better explained at http://research.swtch.com/field'''
    global field_charac
    if (len(msg_in) + nsym) > field_charac: raise ValueError("Message is too long (%i when max is %i)" % (len(msg_in)+nsym, field_charac))
    if gen is None: gen = rs_generator_poly(nsym, fcr, generator)
    msg_in = _bytearray(msg_in)
    msg_out = _bytearray(msg_in) + _bytearray(len(gen)-1)
    lgen = _bytearray([gf_log[gen[j]] for j in xrange(len(gen))])
    msg_in_len = len(msg_in)
    gen_len = len(gen)
    for i in xrange(msg_in_len):
        coef = msg_out[i]
        if coef != 0:
            lcoef = gf_log[coef]
            for j in xrange(1, gen_len):
                msg_out[i + j] ^= gf_exp[lcoef + lgen[j]]
    msg_out[:msg_in_len] = msg_in
    return msg_out
def inverted(msg):
    '''Implements msg[::-1] explicitly to make the library compatible with MicroPython which does not support stepped slices.'''
    return _bytearray(reversed(msg))
def rs_calc_syndromes(msg, nsym, fcr=0, generator=2):
    '''Given the received codeword msg and the number of error correcting symbols (nsym), computes the syndromes polynomial.
    Mathematically, it's essentially equivalent to a Fourrier Transform (Chien search being the inverse).
    '''
    return [0] + [gf_poly_eval(msg, gf_pow(generator, i+fcr)) for i in xrange(nsym)]
def rs_correct_errata(msg_in, synd, err_pos, fcr=0, generator=2):
    '''Forney algorithm, computes the values (error magnitude) to correct the input message.'''
    global field_charac
    msg = _bytearray(msg_in)
    coef_pos = _bytearray(len(msg) - 1 - p for p in err_pos)
    err_loc = rs_find_errata_locator(coef_pos, generator)
    err_eval = inverted(rs_find_error_evaluator(inverted(synd), err_loc, len(err_loc)-1))
    X = _bytearray(len(coef_pos))
    for i in xrange(len(coef_pos)):
        l = field_charac - coef_pos[i]
        X[i] = gf_pow(generator, -l)
    E = _bytearray(len(msg))
    X_len = len(X)
    for i, Xi in enumerate(X):
        Xi_inv = gf_inverse(Xi)
        err_loc_prime = 1
        for j in xrange(X_len):
            if j != i:
                err_loc_prime = gf_mul(err_loc_prime, gf_sub(1, gf_mul(Xi_inv, X[j])))
        if err_loc_prime == 0:
            raise ReedSolomonError("Decoding failed: Forney algorithm could not properly detect where the errors are located (errata locator prime is 0).")
        y = gf_poly_eval(inverted(err_eval), Xi_inv)
        y = gf_mul(gf_pow(Xi, 1-fcr), y)
        magnitude = gf_div(y, err_loc_prime)
        E[err_pos[i]] = magnitude
    msg = gf_poly_add(msg, E)
    return msg
def rs_find_error_locator(synd, nsym, erase_loc=None, erase_count=0):
    '''Find error/errata locator and evaluator polynomials with Berlekamp-Massey algorithm'''
    if erase_loc:
        err_loc = _bytearray(erase_loc)
        old_loc = _bytearray(erase_loc)
    else:
        err_loc = _bytearray([1])
        old_loc = _bytearray([1])
    synd_shift = 0
    if len(synd) > nsym: synd_shift = len(synd) - nsym
    for i in xrange(nsym-erase_count):
        if erase_loc:
            K = erase_count+i+synd_shift
        else:
            K = i+synd_shift
        delta = synd[K]
        for j in xrange(1, len(err_loc)):
            delta ^= gf_mul(err_loc[-(j+1)], synd[K - j])
        old_loc = old_loc + _bytearray([0])
        if delta != 0:
            if len(old_loc) > len(err_loc):
                new_loc = gf_poly_scale(old_loc, delta)
                old_loc = gf_poly_scale(err_loc, gf_inverse(delta))
                err_loc = new_loc
            err_loc = gf_poly_add(err_loc, gf_poly_scale(old_loc, delta))
    for i, x in enumerate(err_loc):
        if x != 0:
            err_loc = err_loc[i:]
            break
    errs = len(err_loc) - 1
    if (errs-erase_count) * 2 + erase_count > nsym:
        raise ReedSolomonError("Too many errors to correct")
    return err_loc
def rs_find_errata_locator(e_pos, generator=2):
    '''Compute the erasures/errors/errata locator polynomial from the erasures/errors/errata positions (the positions must be relative to the x coefficient, eg: "hello worldxxxxxxxxx" is tampered to "h_ll_ worldxxxxxxxxx" with xxxxxxxxx being the ecc of length n-k=9, here the string positions are [1, 4], but the coefficients are reversed since the ecc characters are placed as the first coefficients of the polynomial, thus the coefficients of the erased characters are n-1 - [1, 4] = [18, 15] = erasures_loc to be specified as an argument.'''
    e_loc = [1]
    for i in e_pos:
        e_loc = gf_poly_mul( e_loc, gf_poly_add(_bytearray([1]), [gf_pow(generator, i), 0]) )
    return e_loc
def rs_find_error_evaluator(synd, err_loc, nsym):
    '''Compute the error (or erasures if you supply sigma=erasures locator polynomial, or errata) evaluator polynomial Omega from the syndrome and the error/erasures/errata locator Sigma. Omega is already computed at the same time as Sigma inside the Berlekamp-Massey implemented above, but in case you modify Sigma, you can recompute Omega afterwards using this method, or just ensure that Omega computed by BM is correct given Sigma.'''
    remainder = gf_poly_mul(synd, err_loc)
    remainder = remainder[len(remainder)-(nsym+1):]
    return remainder
def rs_find_errors(err_loc, nmess, generator=2):
    '''Find the roots (ie, where evaluation = zero) of error polynomial by smart bruteforce trial. This is a faster form of chien search, processing only useful coefficients (the ones in the messages) instead of the whole 2^8 range. Besides the speed boost, this also allows to fix a number of issue: correctly decoding when the last ecc byte is corrupted, and accepting messages of length n > 2^8.'''
    err_pos = []
    for i in xrange(nmess):
        if gf_poly_eval(err_loc, gf_pow(generator, i)) == 0:
            err_pos.append(nmess - 1 - i)
    errs = len(err_loc) - 1
    if len(err_pos) != errs:
        raise ReedSolomonError("Too many (or few) errors found by Chien Search for the errata locator polynomial!")
    return _bytearray(err_pos)
def rs_forney_syndromes(synd, pos, nmess, generator=2):
    erase_pos_reversed = _bytearray(nmess-1-p for p in pos)
    fsynd = _bytearray(synd[1:])
    for i in xrange(len(pos)):
        x = gf_pow(generator, erase_pos_reversed[i])
        for j in xrange(len(fsynd) - 1):
            fsynd[j] = gf_mul(fsynd[j], x) ^ fsynd[j + 1]
    return fsynd
def rs_correct_msg(msg_in, nsym, fcr=0, generator=2, erase_pos=None, only_erasures=False):
    '''Reed-Solomon main decoding function'''
    global field_charac
    if len(msg_in) > field_charac:
        raise ValueError("Message is too long (%i when max is %i)" % (len(msg_in), field_charac))
    msg_out = _bytearray(msg_in)
    if erase_pos is None:
        erase_pos = _bytearray()
    else:
        if isinstance(erase_pos, list):
            erase_pos = _bytearray(erase_pos)
        for e_pos in erase_pos:
            msg_out[e_pos] = 0
    if len(erase_pos) > nsym: raise ReedSolomonError("Too many erasures to correct")
    synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)
    if max(synd) == 0:
        return msg_out[:-nsym], msg_out[-nsym:], erase_pos
    if only_erasures:
        err_pos = _bytearray()
    else:
        fsynd = rs_forney_syndromes(synd, erase_pos, len(msg_out), generator)
        err_loc = rs_find_error_locator(fsynd, nsym, erase_count=len(erase_pos))
        err_pos = rs_find_errors(inverted(err_loc), len(msg_out), generator)
        if err_pos is None:
            raise ReedSolomonError("Could not locate error")
    msg_out = rs_correct_errata(msg_out, synd, erase_pos + err_pos, fcr, generator)
    synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)
    if max(synd) > 0:
        raise ReedSolomonError("Could not correct message")
    return msg_out[:-nsym], msg_out[-nsym:], erase_pos + err_pos
def rs_correct_msg_nofsynd(msg_in, nsym, fcr=0, generator=2, erase_pos=None, only_erasures=False):
    '''Reed-Solomon main decoding function, without using the modified Forney syndromes'''
    global field_charac
    if len(msg_in) > field_charac:
        raise ValueError("Message is too long (%i when max is %i)" % (len(msg_in), field_charac))
    msg_out = _bytearray(msg_in)
    if erase_pos is None:
        erase_pos = _bytearray()
    else:
        if isinstance(erase_pos, list):
            erase_pos = _bytearray(erase_pos)
        for e_pos in erase_pos:
            msg_out[e_pos] = 0
    if len(erase_pos) > nsym: raise ReedSolomonError("Too many erasures to correct")
    synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)
    if max(synd) == 0:
        return msg_out[:-nsym], msg_out[-nsym:], []
    erase_loc = None
    erase_count = 0
    if erase_pos:
        erase_count = len(erase_pos)
        msg_out_len = len(msg_out)
        erase_pos_reversed = [msg_out_len-1-eras for eras in erase_pos]
        erase_loc = rs_find_errata_locator(erase_pos_reversed, generator=generator)
    if only_erasures:
        err_loc = inverted(erase_loc)
    else:
        err_loc = rs_find_error_locator(synd, nsym, erase_loc=erase_loc, erase_count=erase_count)
        err_loc = inverted(err_loc)
    err_pos = rs_find_errors(err_loc, len(msg_out), generator)
    if err_pos is None:
        raise ReedSolomonError("Could not locate error")
    msg_out = rs_correct_errata(msg_out, synd, err_pos, fcr=fcr, generator=generator)
    synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)
    if max(synd) > 0:
        raise ReedSolomonError("Could not correct message")
    return msg_out[:-nsym], msg_out[-nsym:], erase_pos + err_pos
def rs_check(msg, nsym, fcr=0, generator=2):
    '''Returns true if the message + ecc has no error of false otherwise (may not always catch a wrong decoding or a wrong message, particularly if there are too many errors -- above the Singleton bound --, but it usually does)'''
    return ( max(rs_calc_syndromes(msg, nsym, fcr, generator)) == 0 )
class RSCodec(object):
    '''
    A Reed Solomon encoder/decoder. After initializing the object, use ``encode`` to encode a
    (byte)string to include the RS correction code, and pass such an encoded (byte)string to
    ``decode`` to extract the original message (if the number of errors allows for correct decoding).
    The ``nsym`` argument is the length of the correction code, and it determines the number of
    error bytes (if I understand this correctly, half of ``nsym`` is correctable)
    '''
    '''
    Modifications by rotorgit 2/3/2015:
    Added support for US FAA ADSB UAT RS FEC, by allowing user to specify
    different primitive polynomial and non-zero first consecutive root (fcr).
    For UAT/ADSB use, set fcr=120 and prim=0x187 when instantiating
    the class; leaving them out will default for previous values (0 and
    0x11d)
    '''
    def __init__(self, nsym=10, nsize=255, fcr=0, prim=0x11d, generator=2, c_exp=8, single_gen=True):
        '''Initialize the Reed-Solomon codec. Note that different parameters change the internal values (the ecc symbols, look-up table values, etc) but not the output result (whether your message can be repaired or not, there is no influence of the parameters).
        nsym : number of ecc symbols (you can repair nsym/2 errors and nsym erasures.
        nsize : maximum length of each chunk. If higher than 255, will use a higher Galois Field, but the algorithm's complexity and computational cost will raise quadratically...
        single_gen : if you want to use the same RSCodec for different nsym parameters (but nsize the same), then set single_gen=False. This is only required for encoding with various number of ecc symbols, as for decoding this is always possible even if single_gen=True.
        '''
        if nsize > 255 and c_exp <= 8:
            c_exp = int(math.log(2 ** (math.floor(math.log(nsize) / math.log(2)) + 1), 2))
        if c_exp != 8 and prim == 0x11d:
            prim = find_prime_polys(generator=generator, c_exp=c_exp, fast_primes=True, single=True)[0]
            if nsize == 255:
                nsize = int(2**c_exp - 1)
        if nsym >= nsize:
            raise ValueError('ECC symbols must be strictly less than the total message length (nsym < nsize).')
        self.nsym = nsym
        self.nsize = nsize
        self.fcr = fcr
        self.prim = prim
        self.generator = generator
        self.c_exp = c_exp
        self.gf_log, self.gf_exp, self.field_charac = init_tables(prim, generator, c_exp)
        if single_gen:
            self.gen = {}
            self.gen[nsym] = rs_generator_poly(nsym, fcr=fcr, generator=generator)
        else:
            self.gen = rs_generator_poly_all(nsize, fcr=fcr, generator=generator)
    def chunk(self, data, chunk_size):
        '''Split a long message into chunks
        DEPRECATED: inlined alternate form so that we can preallocate arrays and hence get faster results with JIT compilers such as PyPy.'''
        for i in xrange(0, len(data), chunk_size):
            chunk = data[i:i+chunk_size]
            yield chunk
    def encode(self, data, nsym=None):
        '''Encode a message (ie, add the ecc symbols) using Reed-Solomon, whatever the length of the message because we use chunking
        Optionally, can set nsym to encode with a different number of error correction symbols, but RSCodec must be initialized with single_gen=False first.
        slice_assign=True allows to speed up the loop quite significantly in JIT compilers such as PyPy by preallocating the output bytearray and slice assigning into it, instead of constantly extending an empty bytearray, but this only works in Python 3, not Python 2, hence is disabled by default for retrocompatibility.
        '''
        global gf_log, gf_exp, field_charac
        gf_log, gf_exp, field_charac = self.gf_log, self.gf_exp, self.field_charac
        nsize, fcr, generator = self.nsize, self.fcr, self.generator
        if not nsym:
            nsym = self.nsym
        gen = self.gen[nsym]
        if isinstance(data, str):
            data = _bytearray(data)
        chunk_size = int(nsize - nsym)
        total_chunks = int(math.ceil(float(len(data)) / float(chunk_size)))
        enc = _bytearray(total_chunks * nsize)
        for i in xrange(0, total_chunks):
            enc[i*nsize:(i+1)*nsize] = rs_encode_msg(data[i*chunk_size:(i+1)*chunk_size], nsym, fcr=fcr, generator=generator, gen=gen)
        return enc
    def decode(self, data, nsym=None, erase_pos=None, only_erasures=False):
        '''Repair a message, whatever its size is, by using chunking. May return a wrong result if number of errors > nsym because then too many errors to be corrected.
        Note that it returns a couple of vars: the repaired messages, and the repaired messages+ecc (useful for checking).
        Usage: rmes, rmesecc = RSCodec.decode(data).
        Optionally: can specify nsym to decode messages of different parameters, erase_pos with a list of erasures positions to double the number of erasures that can be corrected compared to unlocalized errors, only_erasures boolean to specify if we should only look for erasures, which speeds up and doubles the total correction power.
        '''
        global gf_log, gf_exp, field_charac
        gf_log, gf_exp, field_charac = self.gf_log, self.gf_exp, self.field_charac
        if isinstance(data, str):
            data = _bytearray(data)
        if isinstance(erase_pos, list):
            erase_pos = _bytearray(erase_pos)
        if not nsym:
            nsym = self.nsym
        nsize = self.nsize
        fcr = self.fcr
        generator = self.generator
        chunk_size = nsize
        total_chunks = int(math.ceil(float(len(data)) / float(chunk_size)))
        nmes = int(nsize-nsym)
        dec = _bytearray(total_chunks * nmes)
        dec_full = _bytearray(total_chunks * nsize)
        errata_pos_all = _bytearray()
        for i in xrange(0, total_chunks):
            if erase_pos is not None:
                e_pos = [x for x in erase_pos if x < nsize]
                erase_pos = [x - nsize for x in erase_pos if x >= nsize]
            else:
                e_pos = _bytearray()
            rmes, recc, errata_pos = rs_correct_msg(data[i*chunk_size:(i+1)*chunk_size], nsym, fcr=fcr, generator=generator, erase_pos=e_pos, only_erasures=only_erasures)
            dec[i*nmes:(i+1)*nmes] = rmes
            dec_full[i*nsize:(i+1)*nsize] = rmes
            dec_full[i*nsize + nmes:(i+1)*nsize + nmes] = recc
            errata_pos_all.extend(errata_pos)
        return dec, dec_full, errata_pos_all
    def check(self, data, nsym=None):
        '''Check if a message+ecc stream is not corrupted (or fully repaired). Note: may return a wrong result if number of errors > nsym.'''
        if not nsym:
            nsym = self.nsym
        if isinstance(data, str):
            data = _bytearray(data)
        nsize = self.nsize
        fcr = self.fcr
        generator = self.generator
        chunk_size = nsize
        total_chunks = int(math.ceil(float(len(data)) / float(chunk_size)))
        check = [False] * total_chunks
        for i in xrange(0, total_chunks):
            check[i] = rs_check(data[i*chunk_size:(i+1)*chunk_size], nsym, fcr=fcr, generator=generator)
        return check
    def maxerrata(self, nsym=None, errors=None, erasures=None, verbose=False):
        '''Return the Singleton Bound for the current codec, which is the max number of errata (errors and erasures) that the codec can decode/correct.
        Beyond the Singleton Bound (too many errors/erasures), the algorithm will try to raise an exception, but it may also not detect any problem with the message and return 0 errors.
        Hence why you should use checksums if your goal is to detect errors (as opposed to correcting them), as checksums have no bounds on the number of errors, the only limitation being the probability of collisions.
        By default, return a tuple wth the maximum number of errors (2nd output) OR erasures (2nd output) that can be corrected.
        If errors or erasures (not both) is specified as argument, computes the remaining **simultaneous** correction capacity (eg, if errors specified, compute the number of erasures that can be simultaneously corrected).
        Set verbose to True to get print a report.'''
        if not nsym:
            nsym = self.nsym
        maxerrors = int(nsym/2)
        maxerasures = nsym
        if erasures is not None and erasures >= 0:
            if erasures > maxerasures:
                raise ReedSolomonError("Specified number of errors or erasures exceeding the Singleton Bound!")
            maxerrors = int((nsym-erasures)/2)
            if verbose:
                print('This codec can correct up to %i errors and %i erasures simultaneously' % (maxerrors, erasures))
            return maxerrors, erasures
        if errors is not None and errors >= 0:
            if errors > maxerrors:
                raise ReedSolomonError("Specified number of errors or erasures exceeding the Singleton Bound!")
            maxerasures = int(nsym-(errors*2))
            if verbose:
                print('This codec can correct up to %i errors and %i erasures simultaneously' % (errors, maxerasures))
            return errors, maxerasures
        if verbose:
            print('This codec can correct up to %i errors and %i erasures independently' % (maxerrors, maxerasures))
        return maxerrors, maxerasures
"""
rsc = RSCodec(10)
test = rsc.encode([1,2,3,4])
print(test)
test = rsc.encode(bytearray([1,2,3,4]))
print(test)
test = rsc.encode(b'hello world')
print(test)
test = rsc.decode(b'hello world\xed%T\xc4\xfd\xfd\x89\xf3\xa8\xaa')[0]
print(test)
test = rsc.decode(b'heXlo worXd\xed%T\xc4\xfdX\x89\xf3\xa8\xaa')[0]
print(test)
test = rsc.decode(b'hXXXo worXd\xed%T\xc4\xfdX\x89\xf3\xa8\xaa')[0]
print(test)
test = rsc.decode(b'hXXXo worXd\xed%T\xc4\xfdXX\xf3\xa8\xaa')[0]
print(test)
"""

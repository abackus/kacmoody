def bernoulli(n):
    """Returns the nth Bernoulli number, with negative sign convention.
    Note that negative sign convention doesn't matter for Eisenstein numbers
    since it only affects B_1."""
    s = 0
    for u in range(n+1):
        for v in range(u+1):
            s = s + (-1)^v * binomial(u, v) * v^n / (u+1)
    return s

def eisenstein_coeff(k, ell):
    """Returns the terms 1, ..., ell of the normalized Eisenstein series E_2k = G_2k/2zeta(2k), for ell >= 0."""
    terms = [1]
    factor = (4 * k)/bernoulli(2 * k)
    while len(terms) < ell:
        terms.append(-factor * sigma(len(terms), 2*k - 1))
    return terms

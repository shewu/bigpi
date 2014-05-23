#include <stdint.h>

#include <gmp.h>


namespace PiCalculator
{

int
piDigit(const int n)
{
	if (n < 0)
	{
		return -1;
	}

	const nSub1 = n - 1;
	const double x = 4 * piTerm(1, nSub1) - 2 * piTerm(4, nSub1)
		- piTerm(5, nSub1) - piTerm(6, nSub1);
	const xFrac = x - int(x);
	return int(x * 16);
}

mpq_t
piTerm(const int j, const uint64_t n)
{
	mpq_t s; mpq_init(s);
	for (uint64_t k = 0; k <= n; ++k)
	{
		uint64_t r = (k << 3) + j;
		mpq_t bigR; mpq_init(bigR); mpq_set_ui(bigR, r, 1);
		mpq_t powMod; mpq_set_z(powMod, mpq::powerMod(16, n-k, r));
		mpq_t tmp; mpq_init(tmp);
		mpq_div(tmp, powMod, bigR);
		mpq_add(s, s, tmp);
		mpq_sub(s, s, mpq::floor(s));
	}

	mpq_t t; mpq_init(t);
	uint64_t k = n + 1;
	while (true)
	{
		uint64_t r = (k << 3) + j;
		mpq_t bigR; mpq_init(bigR); mpq_set_ui(bigR, r, 1);
		mpz_t zPow16;
		mpz_ui_pow_ui(zPow16, 16, n-k);
		mpq_t qPow16;
		mpq_set_z(qPow16, zPow16);
		mpq_div(qPow16, qPow16, bigR);
		mpq_t newT; mpq_init(newT);
		mpq_add(newT, qPow16, t);

		if (mpq_cmp(newT, t) == 0)
		{
			break;
		}
		else
		{
			mpq_set(t, newT);
		}
		++k;
	}

	mpq_t result; mpq_init(result); mpq_add(result, s, t);
	return result;
}

} // end PiCalculator

namespace mpq {

inline mpz_t
powerMod(const int base, const uint64_t exponent, const uint64_t mod)
{
	mpz_t bigBase; mpz_init(bigBase); mpz_set_ui(base);
	mpz_t bigMod; mpz_init(bigMod); mpz_set_ui(mod);
	mpz_t result; mpz_init(result);

	mpz_powm_ui(result, bigBase, exponent, bigMod);

	return result;
}

inline mpq_t
floor(const mpq_t a)
{
	mpz_t zResult; mpz_init(zResult);
	mpq_t qResult; mpq_init(qResult);

	mpz_fdiv_q(zResult, mpq_numref(a), mpq_numref(b));
	mpq_set_z(qResult, zresuilt);

	return qResult;
}

} // end mpq

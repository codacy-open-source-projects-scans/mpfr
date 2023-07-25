/* mpfr_strtofr -- set a floating-point number from a string

Copyright 2004-2023 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
https://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

/* FIXME: MPFR_SADD_OVERFLOW is sometimes called with
     MPFR_EXP_MIN+2, MPFR_EXP_MAX-2
   which is not supported by the current implementation (if an argument
   is >= 0 and the other one is < 0, a simple addition is done, so that
   it is not guaranteed that the result is between these two values).
   It is not clear whether one can find failing cases in practice, as
   previous MPFR_SADD_OVERFLOW calls exclude some underflow/overflow
   cases. One would need to analyze the full code and do a case-by-case
   study.
   Moreover, mpfr_check_range() is called with an exponent field that
   may be far from the MPFR_EMIN_MIN and MPFR_EMAX_MAX limits (which is
   discouraged), and if the exponent field happens to be equal to one of
   the special values, erratic behavior will occur. It happens that this
   is probably not possible with the current code, due to the
     MPFR_SADD_OVERFLOW (exp, exp, ysize_bits, [...])
   where ysize_bits >= GMP_NUMB_BITS > 4, but this could silently be
   broken with future code.
*/

#include <ctype.h>  /* For isspace */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

#define MPFR_MAX_BASE 62

struct parsed_string {
  int            negative; /* non-zero iff the number is negative */
  int            base;     /* base of the string */
  unsigned char *mantissa; /* raw significand (without any point) */
  unsigned char *mant;     /* stripped significand (without starting and
                              ending zeroes). This points inside the area
                              allocated for the mantissa field. */
  size_t         prec;     /* length of mant (zero for +/-0) */
  size_t         alloc;    /* allocation size of mantissa */
  mpfr_exp_t     exp_base; /* number of digits before the point, + exponent
                              except in case of binary exponent (exp_bin) */
  mpfr_exp_t     exp_bin;  /* binary exponent of the pxxx format for
                              base = 2 or 16 */
};

/* This table has been generated by the following program.
   For 2 <= b <= MPFR_MAX_BASE,
   RedInvLog2Table[b-2][0] / RedInvLog2Table[b-2][1]
   is an upper approximation to log(2)/log(b), no larger than 1.
   Note: these numbers must fit on 16 bits, thus unsigned int is OK.
*/
static const unsigned int RedInvLog2Table[MPFR_MAX_BASE-1][2] = {
  {1, 1},
  {53, 84},
  {1, 2},
  {4004, 9297},
  {53, 137},
  {2393, 6718},
  {1, 3},
  {665, 2108},
  {4004, 13301},
  {949, 3283},
  {53, 190},
  {5231, 19357},
  {2393, 9111},
  {247, 965},
  {1, 4},
  {4036, 16497},
  {665, 2773},
  {5187, 22034},
  {4004, 17305},
  {51, 224},
  {949, 4232},
  {3077, 13919},
  {53, 243},
  {73, 339},
  {5231, 24588},
  {665, 3162},
  {2393, 11504},
  {4943, 24013},
  {247, 1212},
  {3515, 17414},
  {1, 5},
  {4415, 22271},
  {4036, 20533},
  {263, 1349},
  {665, 3438},
  {1079, 5621},
  {5187, 27221},
  {2288, 12093},
  {4004, 21309},
  {179, 959},
  {51, 275},
  {495, 2686},
  {949, 5181},
  {3621, 19886},
  {3077, 16996},
  {229, 1272},
  {53, 296},
  {109, 612},
  {73, 412},
  {1505, 8537},
  {5231, 29819},
  {283, 1621},
  {665, 3827},
  {32, 185},
  {2393, 13897},
  {1879, 10960},
  {4943, 28956},
  {409, 2406},
  {247, 1459},
  {231, 1370},
  {3515, 20929} };
#if 0
#define N 8
int main ()
{
  unsigned long tab[N];
  int i, n, base;
  mpfr_t x, y;
  mpq_t q1, q2;
  int overflow = 0, base_overflow;

  mpfr_init2 (x, 200);
  mpfr_init2 (y, 200);
  mpq_init (q1);
  mpq_init (q2);

  for (base = 2 ; base < 63 ; base ++)
    {
      mpfr_set_ui (x, base, MPFR_RNDN);
      mpfr_log2 (x, x, MPFR_RNDN);
      mpfr_ui_div (x, 1, x, MPFR_RNDN);
      printf ("Base: %d x=%e ", base, mpfr_get_d1 (x));
      for (i = 0 ; i < N ; i++)
        {
          mpfr_floor (y, x);
          tab[i] = mpfr_get_ui (y, MPFR_RNDN);
          mpfr_sub (x, x, y, MPFR_RNDN);
          mpfr_ui_div (x, 1, x, MPFR_RNDN);
        }
      for (i = N-1 ; i >= 0 ; i--)
        if (tab[i] != 0)
          break;
      mpq_set_ui (q1, tab[i], 1);
      for (i = i-1 ; i >= 0 ; i--)
          {
            mpq_inv (q1, q1);
            mpq_set_ui (q2, tab[i], 1);
            mpq_add (q1, q1, q2);
          }
      printf("Approx: ", base);
      mpq_out_str (stdout, 10, q1);
      printf (" = %e\n", mpq_get_d (q1) );
      fprintf (stderr, "{");
      mpz_out_str (stderr, 10, mpq_numref (q1));
      fprintf (stderr, "UL, ");
      mpz_out_str (stderr, 10, mpq_denref (q1));
      fprintf (stderr, "UL},\n");
      if (mpz_cmp_ui (mpq_numref (q1), 1<<16-1) >= 0
          || mpz_cmp_ui (mpq_denref (q1), 1<<16-1) >= 0)
        overflow = 1, base_overflow = base;
    }

  mpq_clear (q2);
  mpq_clear (q1);
  mpfr_clear (y);
  mpfr_clear (x);
  if (overflow )
    printf ("OVERFLOW for base =%d!\n", base_overflow);
}
#endif


/* Compatible with any locale, but one still assumes that 'a', 'b', 'c',
   ..., 'z', and 'A', 'B', 'C', ..., 'Z' are consecutive values (like
   in any ASCII-based character set). */
static int
digit_value_in_base (int c, int base)
{
  int digit;

  MPFR_ASSERTD (base > 0 && base <= MPFR_MAX_BASE);

  if (c >= '0' && c <= '9')
    digit = c - '0';
  else if (c >= 'a' && c <= 'z')
    digit = (base >= 37) ? c - 'a' + 36 : c - 'a' + 10;
  else if (c >= 'A' && c <= 'Z')
    digit = c - 'A' + 10;
  else
    return -1;

  return MPFR_LIKELY (digit < base) ? digit : -1;
}

/* Compatible with any locale, but one still assumes that 'a', 'b', 'c',
   ..., 'z', and 'A', 'B', 'C', ..., 'Z' are consecutive values (like
   in any ASCII-based character set). */
/* TODO: support EBCDIC. */
static int
fast_casecmp (const char *s1, const char *s2)
{
  unsigned char c1, c2;

  do
    {
      c2 = *(const unsigned char *) s2++;
      if (c2 == '\0')
        return 0;
      c1 = *(const unsigned char *) s1++;
      if (c1 >= 'A' && c1 <= 'Z')
        c1 = c1 - 'A' + 'a';
    }
  while (c1 == c2);
  return 1;
}

/* Parse a string and fill pstr.
   Return the advanced ptr too.
   It returns:
      -1 if invalid string,
      0 if special string (like nan),
      1 if the string is OK.
      2 if overflows
   So it doesn't return the ternary value
   BUT if it returns 0 (NAN or INF), the ternary value is also '0'
   (ie NAN and INF are exact) */
static int
parse_string (mpfr_ptr x, struct parsed_string *pstr,
              const char **string, int base)
{
  const char *str = *string;
  unsigned char *mant;
  int point;
  int res = -1;  /* Invalid input return value */
  const char *prefix_str;
  int decimal_point;

  decimal_point = (unsigned char) MPFR_DECIMAL_POINT;

  /* Init variable */
  pstr->mantissa = NULL;

  /* Optional leading whitespace */
  /* For non-"C" locales, the ISO C standard allows isspace(0) to
     return true. So we need to stop explicitly on '\0'. */
  while (*str != '\0' && isspace ((unsigned char) *str))
    str++;

  /* An optional sign `+' or `-' */
  pstr->negative = (*str == '-');
  if (*str == '-' || *str == '+')
    str++;

  /* Can be case-insensitive NAN */
  if (fast_casecmp (str, "@nan@") == 0)
    {
      str += 5;
      goto set_nan;
    }
  if (base <= 16 && fast_casecmp (str, "nan") == 0)
    {
      str += 3;
    set_nan:
      /* Check for "(dummychars)" */
      if (*str == '(')
        {
          const char *s;
          for (s = str+1 ; *s != ')' ; s++)
            if (!(*s >= 'A' && *s <= 'Z')
                && !(*s >= 'a' && *s <= 'z')
                && !(*s >= '0' && *s <= '9')
                && *s != '_')
              break;
          if (*s == ')')
            str = s+1;
        }
      *string = str;
      MPFR_SET_NAN(x);
      /* MPFR_RET_NAN not used as the return value isn't a ternary value */
      __gmpfr_flags |= MPFR_FLAGS_NAN;
      return 0;
    }

  /* Can be case-insensitive INF */
  if (fast_casecmp (str, "@inf@") == 0)
    {
      str += 5;
      goto set_inf;
    }
  if (base <= 16 && fast_casecmp (str, "infinity") == 0)
    {
      str += 8;
      goto set_inf;
    }
  if (base <= 16 && fast_casecmp (str, "inf") == 0)
    {
      str += 3;
    set_inf:
      *string = str;
      MPFR_SET_INF (x);
      (pstr->negative) ? MPFR_SET_NEG (x) : MPFR_SET_POS (x);
      return 0;
    }

  /* If base=0 or 16, it may include '0x' prefix */
  prefix_str = NULL;
  if ((base == 0 || base == 16) && str[0]=='0'
      && (str[1]=='x' || str[1] == 'X'))
    {
      prefix_str = str;
      base = 16;
      str += 2;
    }
  /* If base=0 or 2, it may include '0b' prefix */
  if ((base == 0 || base == 2) && str[0]=='0'
      && (str[1]=='b' || str[1] == 'B'))
    {
      prefix_str = str;
      base = 2;
      str += 2;
    }
  /* Else if base=0, we assume decimal base */
  if (base == 0)
    base = 10;
  pstr->base = base;

  /* Alloc mantissa */
  pstr->alloc = (size_t) strlen (str) + 1;
  pstr->mantissa = (unsigned char*) mpfr_allocate_func (pstr->alloc);

  /* Read mantissa digits */
 parse_begin:
  mant = pstr->mantissa;
  point = 0;
  pstr->exp_base = 0;
  pstr->exp_bin  = 0;

  for (;;) /* Loop until an invalid character is read */
    {
      int c = (unsigned char) *str++;
      /* The cast to unsigned char is needed because of digit_value_in_base;
         decimal_point uses this convention too. */
      if (c == '.' || c == decimal_point)
        {
          if (MPFR_UNLIKELY(point)) /* Second '.': stop parsing */
            break;
          point = 1;
          continue;
        }
      c = digit_value_in_base (c, base);
      if (c == -1)
        break;
      MPFR_ASSERTN (c >= 0); /* c is representable in an unsigned char */
      *mant++ = (unsigned char) c;
      if (!point)
        pstr->exp_base ++;
    }
  str--; /* The last read character was invalid */

  /* Update the # of char in the mantissa */
  pstr->prec = mant - pstr->mantissa;
  /* Check if there are no characters in the mantissa (Invalid argument) */
  if (pstr->prec == 0)
    {
      /* Check if there was a prefix (in such a case, we have to read
         again the mantissa without skipping the prefix)
         The allocated mantissa is still big enough since we will
         read only 0, and we alloc one more char than needed.
         FIXME: Not really friendly. Maybe cleaner code? */
      if (prefix_str != NULL)
        {
          str = prefix_str;
          prefix_str = NULL;
          goto parse_begin;
        }
      goto end;
    }

  /* Valid entry */
  res = 1;
  MPFR_ASSERTD (pstr->exp_base >= 0);

  /* FIXME: In the code below (both cases), if the exponent from the
     string is large, it will be replaced by MPFR_EXP_MIN or MPFR_EXP_MAX,
     i.e. it will have a different value. This may not change the result
     in most cases, but there is no guarantee on very long strings when
     mpfr_exp_t is a 32-bit type, as the exponent could be brought back
     to the current exponent range. */

  /* an optional exponent (e or E, p or P, @) */
  if ( (*str == '@' || (base <= 10 && (*str == 'e' || *str == 'E')))
       && (!isspace((unsigned char) str[1])) )
    {
      char *endptr;
      /* the exponent digits are kept in ASCII */
      mpfr_exp_t sum;
      long read_exp = strtol (str + 1, &endptr, 10);
      if (endptr != str+1)
        str = endptr;
      sum =
        read_exp < MPFR_EXP_MIN ? (str = endptr, MPFR_EXP_MIN) :
        read_exp > MPFR_EXP_MAX ? (str = endptr, MPFR_EXP_MAX) :
        (mpfr_exp_t) read_exp;
      MPFR_SADD_OVERFLOW (sum, sum, pstr->exp_base,
                          mpfr_exp_t, mpfr_uexp_t,
                          MPFR_EXP_MIN, MPFR_EXP_MAX,
                          res = 2, res = 3);
      /* Since exp_base was positive, read_exp + exp_base can't
         do a negative overflow. */
      MPFR_ASSERTD (res != 3);
      pstr->exp_base = sum;
    }
  else if ((base == 2 || base == 16)
           && (*str == 'p' || *str == 'P')
           && (!isspace((unsigned char) str[1])))
    {
      char *endptr;
      long read_exp = strtol (str + 1, &endptr, 10);
      if (endptr != str+1)
        str = endptr;
      pstr->exp_bin =
        read_exp < MPFR_EXP_MIN ? (str = endptr, MPFR_EXP_MIN) :
        read_exp > MPFR_EXP_MAX ? (str = endptr, MPFR_EXP_MAX) :
        (mpfr_exp_t) read_exp;
    }

  /* Remove 0's at the beginning and end of mantissa[0..prec-1] */
  mant = pstr->mantissa;
  for ( ; (pstr->prec > 0) && (*mant == 0) ; mant++, pstr->prec--)
    if (MPFR_LIKELY (pstr->exp_base != MPFR_EXP_MIN))
      pstr->exp_base--;
  for ( ; (pstr->prec > 0) && (mant[pstr->prec - 1] == 0); pstr->prec--);
  pstr->mant = mant;

  /* Check if x = 0 */
  if (pstr->prec == 0)
    {
      MPFR_SET_ZERO (x);
      if (pstr->negative)
        MPFR_SET_NEG(x);
      else
        MPFR_SET_POS(x);
      res = 0;
    }

  *string = str;
 end:
  if (pstr->mantissa != NULL && res != 1)
    mpfr_free_func (pstr->mantissa, pstr->alloc);
  return res;
}

/* Transform a parsed string to a mpfr_t according to the rounding mode
   and the precision of x.
   Returns the ternary value. */
static int
parsed_string_to_mpfr (mpfr_ptr x, struct parsed_string *pstr, mpfr_rnd_t rnd)
{
  mpfr_prec_t precx, prec, ysize_bits, pstr_size;
  mpfr_exp_t exp;
  mp_limb_t *result;
  int count, exact;
  mp_size_t ysize, real_ysize, diff_ysize;
  int res, err;
  const int extra_limbs = GMP_NUMB_BITS >= 12 ? 1 : 2; /* see below */
  MPFR_ZIV_DECL (loop);
  MPFR_TMP_DECL (marker);

  MPFR_LOG_FUNC
    (("rnd=%d", rnd),
     ("", 0));

  /* initialize the working precision */
  precx = MPFR_GET_PREC (x);
  prec = precx + MPFR_INT_CEIL_LOG2 (precx);

  /* Compute the value y of the leading characters as long as rounding is not
     possible.
     Note: We have some integer overflow checking using MPFR_EXP_MIN and
     MPFR_EXP_MAX in this loop. Thanks to the large margin between these
     extremal values of the mpfr_exp_t type and the valid minimum/maximum
     exponents, such integer overflows would correspond to real underflow
     or overflow on the result (possibly except in huge precisions, which
     are disregarded here; anyway, in practice, such issues could occur
     only with 32-bit precision and exponent types). Such checks could be
     extended to real early underflow/overflow checking, in order to avoid
     useless computations in such cases; in such a case, be careful that
     the approximation errors need to be taken into account. */
  MPFR_TMP_MARK(marker);
  MPFR_ZIV_INIT (loop, prec);
  for (;;)
    {
      mp_limb_t *y0, *y;

      /* y will be regarded as a number with precision prec. */
      ysize = MPFR_PREC2LIMBS (prec);
      /* prec bits corresponds to ysize limbs */
      ysize_bits = (mpfr_prec_t) ysize * GMP_NUMB_BITS;
      MPFR_ASSERTD (ysize_bits >= prec);
      /* and to ysize_bits >= prec > precx bits. */
      /* We need to allocate one more limb as specified by mpn_set_str
         (a limb may be written in rp[rn]). Note that the manual of GMP
         up to 5.1.3 was incorrect on this point.
         See the following discussion:
         https://gmplib.org/list-archives/gmp-bugs/2013-December/003267.html */
      y0 = MPFR_TMP_LIMBS_ALLOC (2 * ysize + extra_limbs + 1);
      y = y0 + ysize; /* y has (ysize + extra_limbs + 1) allocated limbs */

      /* pstr_size is the number of bytes we want to read from pstr->mant
         to fill at least ysize full limbs with mpn_set_str.
         We must have base^(pstr_size-1) >= (2^(GMP_NUMB_BITS))^ysize
         (in the worst case, the first digit is one and all others are zero).
         i.e., pstr_size >= 1 + ysize*GMP_NUMB_BITS/log2(base)
          Since ysize ~ prec/GMP_NUMB_BITS and prec < Umax/2 =>
          ysize*GMP_NUMB_BITS can not overflow.
         We compute pstr_size = 1 + ceil(ysize_bits * Num / Den)
          where 1/log2(base) <= Num/Den <= 1
         It is not exactly ceil(1/log2(base)) but could be one more (base 2).
         Quite ugly since it tries to avoid overflow:
         let Num = RedInvLog2Table[pstr->base-2][0]
         and Den = RedInvLog2Table[pstr->base-2][1],
         and ysize_bits = a*Den+b,
         then ysize_bits * Num/Den = a*Num + (b * Num)/Den,
         thus ceil(ysize_bits * Num/Den) = a*Num + floor(b * Num + Den - 1)/Den

         Note: denoting m = pstr_size and n = ysize_bits, assuming we have
         m = 1 + ceil(n/log2(b)), i.e., b^(m-1) >= 2^n > b^(m-2), then
         b^(m-1)/2^n < b, and since we consider m characters of the input,
         the corresponding part is less than b^m < b^2*2^n.
         This implies that if b^2 < 2^GMP_NUMB_BITS, which for b <= 62 holds
         for GMP_NUMB_BITS >= 12, we have real_ysize <= ysize+1 below
         (this also implies that for GMP_NUMB_BITS >= 13, the number of bits
         of y[real_ysize-1] below is less than GMP_NUMB_BITS, thus
         count < GMP_NUMB_BITS).
         Warning: for GMP_NUMB_BITS=8, we can have real_ysize = ysize + 2!
         Hence the allocation above for ysize + extra_limbs limbs.
      */
      {
        unsigned int Num = RedInvLog2Table[pstr->base-2][0];
        unsigned int Den = RedInvLog2Table[pstr->base-2][1];
        MPFR_ASSERTD (Num <= Den && Den <= 65535); /* thus no overflow */
        pstr_size = (ysize_bits / Den) * Num
          + ((unsigned long) (ysize_bits % Den) * Num + Den - 1) / Den
          + 1;
        MPFR_ASSERTD (pstr_size <= 1 + ysize_bits);
      }

      /* Since pstr_size corresponds to at least ysize_bits bits,
         and ysize_bits >= prec, the weight of the neglected part of
         pstr->mant (if any) is < ulp(y) < ulp(x). */

      /* If the number of wanted bytes is more than what is available
         in pstr->mant, i.e. pstr->prec, reduce it to pstr->prec. */
      if (pstr_size > pstr->prec)
        pstr_size = pstr->prec;

      MPFR_ASSERTD (pstr_size >= 0);

      /* Convert str (potentially truncated to pstr_size) into binary.
         Note that pstr->mant is big endian, thus no offset is needed. */
      real_ysize = mpn_set_str (y, pstr->mant, pstr_size, pstr->base);

      /* See above for the explanation of the following assertion. */
      MPFR_ASSERTD (real_ysize <= ysize + extra_limbs);

      /* The Boolean "exact" will attempt to track exactness of the result:
         If it is true, then this means that the result is exact, allowing
         termination, even though the rounding test may not succeed.
         Conversely, if the result is exact, then "exact" will not
         necessarily be true at the end of the Ziv loop, but we will need
         to make sure that at some point, "exact" will be true in order to
         guarantee termination. FIXME: check that. */
      /* First, consider the part of the input string that has been ignored.
         Note that the trailing zeros have been removed in parse_string, so
         that if something has been ignored, it must be non-zero. */
      exact = pstr_size == pstr->prec;

      /* Normalize y and set the initial value of its exponent exp, which
         is 0 when y is not shifted.
         Since pstr->mant was normalized, mpn_set_str guarantees that
         the most significant limb is non-zero. */
      MPFR_ASSERTD (y[real_ysize - 1] != 0); /* mpn_set_str guarantees this */
      count_leading_zeros (count, y[real_ysize - 1]);
      diff_ysize = ysize - real_ysize;
      MPFR_LOG_MSG (("diff_ysize = %ld\n", (long) diff_ysize));
      if (diff_ysize >= 0)
        {
          /* We have enough limbs to store {y, real_ysize} exactly
             in {y, ysize}, so that we can do a left shift, without
             losing any information ("exact" will not change). */
          if (count != 0)
            mpn_lshift (y + diff_ysize, y, real_ysize, count);
          if (diff_ysize > 0)
            {
              if (count == 0)
                mpn_copyd (y + diff_ysize, y, real_ysize);
              MPN_ZERO (y, diff_ysize);
            }
          /* exp = negation of the total shift count, avoiding overflows. */
          exp = - ((mpfr_exp_t) diff_ysize * GMP_NUMB_BITS + count);
        }
      else
        {
          /* Shift {y, real_ysize} for (GMP_NUMB_BITS - count) bits to the
             right, and put the ysize most significant limbs into {y, ysize}.
             We have either real_ysize = ysize + 1 or real_ysize = ysize + 2
             (only possible with extra_limbs == 2). */
          MPFR_ASSERTD (diff_ysize == -1 ||
                        (extra_limbs == 2 && diff_ysize == -2));
          if (count != 0)
            {
              /* Before doing the shift, consider the limb that will entirely
                 be lost if real_ysize = ysize + 2. */
              exact = exact && (diff_ysize == -1 || y[0] == MPFR_LIMB_ZERO);
              /* mpn_rshift allows overlap, provided destination <= source */
              /* FIXME: The bits lost due to mpn_rshift are not taken
                 into account in the error analysis below! */
              if (mpn_rshift (y, y - (diff_ysize + 1), real_ysize,
                              GMP_NUMB_BITS - count) != MPFR_LIMB_ZERO)
                exact = 0; /* some non-zero bits have been shifted out */
            }
          else
            {
              /* the case real_ysize = ysize + 2 with count = 0 cannot happen
                 even with GMP_NUMB_BITS = 8 since 62^2 < 256^2/2 */
              MPFR_ASSERTD (diff_ysize == -1);
              exact = exact && y[0] == MPFR_LIMB_ZERO;
              /* copy {y+real_ysize-ysize, ysize} to {y, ysize} */
              mpn_copyi (y, y + 1, real_ysize - 1);
            }
          /* exp = shift count */
          /* TODO: add some explanations about what exp means exactly. */
          exp = GMP_NUMB_BITS * (- (mpfr_exp_t) diff_ysize) - count;
        }

      /* compute base^(exp_base - pstr_size) on n limbs */
      if (IS_POW2 (pstr->base))
        {
          /* Base: 2, 4, 8, 16, 32 */
          int pow2;
          mpfr_exp_t tmp;

          MPFR_LOG_MSG (("case 1 (base = power of 2)\n", 0));

          count_leading_zeros (pow2, (mp_limb_t) pstr->base);
          pow2 = GMP_NUMB_BITS - pow2 - 1; /* base = 2^pow2 */
          MPFR_ASSERTD (0 < pow2 && pow2 <= 5);
          /* exp += pow2 * (pstr->exp_base - pstr_size) + pstr->exp_bin
             with overflow checking
             and check that we can add/subtract 2 to exp without overflow */
          MPFR_SADD_OVERFLOW (tmp, pstr->exp_base, -pstr_size,
                              mpfr_exp_t, mpfr_uexp_t,
                              MPFR_EXP_MIN, MPFR_EXP_MAX,
                              goto overflow, goto underflow);
          if (tmp > 0 && MPFR_EXP_MAX / pow2 <= tmp)
            goto overflow;
          else if (tmp < 0 && MPFR_EXP_MIN / pow2 >= tmp)
            goto underflow;
          tmp *= pow2;
          MPFR_SADD_OVERFLOW (tmp, tmp, pstr->exp_bin,
                              mpfr_exp_t, mpfr_uexp_t,
                              MPFR_EXP_MIN, MPFR_EXP_MAX,
                              goto overflow, goto underflow);
          MPFR_SADD_OVERFLOW (exp, exp, tmp,
                              mpfr_exp_t, mpfr_uexp_t,
                              MPFR_EXP_MIN+2, MPFR_EXP_MAX-2,
                              goto overflow, goto underflow);

          result = y;
          err = 0;
        }
      /* case non-power-of-two-base, and pstr->exp_base > pstr_size */
      else if (pstr->exp_base > pstr_size)
        {
          mp_limb_t *z;
          mpfr_exp_t exp_z;

          MPFR_LOG_MSG (("case 2 (exp_base > pstr_size)\n", 0));

          result = MPFR_TMP_LIMBS_ALLOC (2 * ysize + 1);

          /* z = base^(exp_base-pstr_size) using space allocated at y-ysize */
          z = y0;
          /* NOTE: exp_base-pstr_size can't overflow since pstr_size > 0 */
          err = mpfr_mpn_exp (z, &exp_z, pstr->base,
                              pstr->exp_base - pstr_size, ysize);
          if (err == -2)
            goto overflow;
          exact = exact && (err == -1);

          /* If exact is non zero, then z equals exactly the value of the
             pstr_size most significant digits from pstr->mant, i.e., the
             only difference can come from the neglected pstr->prec-pstr_size
             least significant digits of pstr->mant.
             If exact is zero, then z is rounded toward zero with respect
             to that value. */

          /* multiply(y = 0.mant[0]...mant[pr-1])_base by base^(exp-g):
             since both y and z are rounded toward zero, so is "result" */
          mpn_mul_n (result, y, z, ysize);

          /* compute the error on the product */
          if (err == -1)
            err = 0;
          err ++;

          /* compute the exponent of y */
          /* exp += exp_z + ysize_bits with overflow checking
             and check that we can add/subtract 2 to exp without overflow */
          MPFR_SADD_OVERFLOW (exp_z, exp_z, ysize_bits,
                              mpfr_exp_t, mpfr_uexp_t,
                              MPFR_EXP_MIN, MPFR_EXP_MAX,
                              goto overflow, goto underflow);
          MPFR_SADD_OVERFLOW (exp, exp, exp_z,
                              mpfr_exp_t, mpfr_uexp_t,
                              MPFR_EXP_MIN+2, MPFR_EXP_MAX-2,
                              goto overflow, goto underflow);

          /* normalize result */
          if (MPFR_LIMB_MSB (result[2 * ysize - 1]) == 0)
            {
              mp_limb_t *r = result + ysize - 1;
              mpn_lshift (r, r, ysize + 1, 1);
              /* Overflow checking not needed */
              exp --;
            }

          /* if the low ysize limbs of {result, 2*ysize} are all zero,
             then the result is still "exact" (if it was before) */
          exact = exact && (mpn_scan1 (result, 0) >= ysize_bits);
          result += ysize;
        }
      /* case exp_base < pstr_size */
      else if (pstr->exp_base < pstr_size)
        {
          mp_limb_t *z;
          mpfr_exp_t exp_z;

          MPFR_LOG_MSG (("case 3 (exp_base < pstr_size)\n", 0));

          result = MPFR_TMP_LIMBS_ALLOC (3 * ysize + 1);

          /* y0 = y * K^ysize */
          MPN_ZERO (y0, ysize);

          /* pstr_size - pstr->exp_base can overflow */
          exp_z = pstr->exp_base == MPFR_EXP_MIN ?
            MPFR_EXP_MAX : -pstr->exp_base;  /* avoid integer overflow */
          MPFR_SADD_OVERFLOW (exp_z, pstr_size, exp_z,
                              mpfr_exp_t, mpfr_uexp_t,
                              MPFR_EXP_MIN, MPFR_EXP_MAX,
                              goto underflow, goto overflow);

          /* (z, exp_z) = base^(pstr_size - exp_base) */
          z = result + 2*ysize + 1;
          err = mpfr_mpn_exp (z, &exp_z, pstr->base, exp_z, ysize);

          /* Now {z, ysize} * 2^(exp_z_out - ysize_bits) is an approximation
             to base^exp_z_in (denoted b^e below), rounded toward zero, with:
             * if err = -1, the result is exact;
             * if err = -2, an overflow occurred in the computation of exp_z;
             * otherwise the error is bounded by 2^err ulps.
             Thus the exact value of b^e is between z and z + 2^err, where
             z is {z, ysize} properly scaled by a power of 2. Then the error
             will be:
               y/b^e - trunc(y/z) = eps1 + eps2
             with
               eps1 = y/b^e - y/z <= 0
               eps2 = y/z - trunc(y/z) >= 0
             thus the errors will (partly) compensate, giving a bound
             max(|eps1|,|eps2|).
             In addition, there is a 3rd error eps3 since y might be the
             conversion of only a part of the character string, and/or y
             might be truncated by the mpn_rshift call above:
               eps3 = exact_y/b^e - y/b^e >= 0.
          */
          if (err == -2)
            goto underflow; /* FIXME: Sure? */
          else if (err == -1)
            err = 0; /* see the note below */
          else
            exact = 0;

          /* exp -= exp_z + ysize_bits with overflow checking
             and check that we can add/subtract 2 to exp without overflow */
          MPFR_SADD_OVERFLOW (exp_z, exp_z, ysize_bits,
                              mpfr_exp_t, mpfr_uexp_t,
                              MPFR_EXP_MIN, MPFR_EXP_MAX,
                              goto underflow, goto overflow);
          MPFR_SADD_OVERFLOW (exp, exp, -exp_z,
                              mpfr_exp_t, mpfr_uexp_t,
                              MPFR_EXP_MIN+2, MPFR_EXP_MAX-2,
                              goto overflow, goto underflow);

          /* Compute the integer division y/z rounded toward zero.
             The quotient will be put at result + ysize (size: ysize + 1),
             and the remainder at result (size: ysize).
             Both the dividend {y, 2*ysize} and the divisor {z, ysize} are
             normalized, i.e., the most significant bit of their most
             significant limb is 1. */
          MPFR_ASSERTD (MPFR_LIMB_MSB (y0[2 * ysize - 1]) != 0);
          MPFR_ASSERTD (MPFR_LIMB_MSB (z[ysize - 1]) != 0);
          mpn_tdiv_qr (result + ysize, result, (mp_size_t) 0, y0,
                       2 * ysize, z, ysize);

          /* The truncation error of the mpn_tdiv_qr call (eps2 above) is at
             most 1 ulp. Idem for the error eps3, which has the same sign,
             thus eps2 + eps3 <= 2 ulps.
             FIXME: For eps3, this is not obvious and should be explained.
             For the error eps1 coming from the approximation to b^e,
             we have (still up to a power-of-2 normalization):
             y/z - y/b^e = y * (b^e-z) / (z * b^e) <= y * 2^err / (z * b^e).
             We have to convert that error in terms of ulp(trunc(y/z)).
             We first have ulp(trunc(y/z)) = ulp(y/z).

             FIXME: There must be some discussion about the exponents,
                    because up to a power of 2, 1/2 <= |y/z| < 1 and
                    1 <= |y/z| < 2 are equivalent and give no information.
                    Moreover 1/2 <= b^e < 1 has not been explained and may
                    hide mistakes since one may have 1/2 <= z < 1 < b^e.

             Since both y and z are normalized, the quotient
             {result+ysize, ysize+1} has exactly ysize limbs, plus maybe one
             bit (this corresponds to the MPFR_ASSERTD below):
             * if the quotient has exactly ysize limbs, then 1/2 <= |y/z| < 1
               (up to a power of 2) and since 1/2 <= b^e < 1, the error is at
               most 2^(err+1) ulps;
             * if the quotient has one extra bit, then 1 <= |y/z| < 2
               (up to a power of 2) and since 1/2 <= b^e < 1, the error is at
               most 2^(err+2) ulps; but since we will shift the result right
               below by one bit, the final error will be at most 2^(err+1) ulps
               too.

             Thus the error is:
             * at most 2^(err+1) ulps for eps1
             * at most 2 ulps for eps2 + eps3, which is of opposite sign
             and we can bound the error by 2^(err+1) ulps in all cases.

             Note: If eps1 was 0, the error would be bounded by 2 ulps,
             thus replacing err = -1 by err = 0 above was the right thing
             to do, since 2^(0+1) = 2.
          */
          MPFR_ASSERTD (result[2 * ysize] <= 1);

          err += 1; /* see above for the explanation of the +1 term */

          /* if the remainder of the division is zero, then the result is
             still "exact" if it was before */
          exact = exact && (mpn_popcount (result, ysize) == 0);

          /* normalize result */
          if (result[2 * ysize] == MPFR_LIMB_ONE)
            {
              mp_limb_t *r = result + ysize;

              exact = exact && ((*r & MPFR_LIMB_ONE) == 0);
              mpn_rshift (r, r, ysize + 1, 1);
              /* Overflow Checking not needed */
              exp ++;
            }
          result += ysize;
        }
      /* case exp_base = pstr_size: no multiplication or division needed */
      else
        {
          MPFR_LOG_MSG (("case 4 (exp_base = pstr_size)\n", 0));

          /* base^(exp-pr) = 1             nothing to compute */
          result = y;
          err = 0;
        }

      MPFR_LOG_MSG (("exact = %d, err = %d, precx = %Pd\n",
                     exact, err, precx));

      /* at this point, result is an approximation rounded toward zero
         of the pstr_size most significant digits of pstr->mant, with
         equality in case exact is non-zero. */

      /* test if rounding is possible, and if so exit the loop.
         Note: we also need to be able to determine the correct ternary value,
         thus we use the precx + (rnd == MPFR_RNDN) trick.
         For example if result = xxx...xxx111...111 and rnd = RNDN,
         then we know the correct rounding is xxx...xx(x+1), but we cannot know
         the correct ternary value. */
      if (exact || mpfr_round_p (result, ysize, ysize_bits - err - 1,
                                 precx + (rnd == MPFR_RNDN)))
        break;

      /* update the prec for next loop */
      MPFR_ZIV_NEXT (loop, prec);
    } /* loop */
  MPFR_ZIV_FREE (loop);

  /* round y */
  if (mpfr_round_raw (MPFR_MANT (x), result, ysize_bits,
                      pstr->negative, precx, rnd, &res))
    {
      /* overflow when rounding y */
      MPFR_MANT (x)[MPFR_LIMB_SIZE (x) - 1] = MPFR_LIMB_HIGHBIT;
      /* Overflow Checking not needed */
      exp ++;
    }

  /* Note: if exact <> 0, then the approximation {result, ysize} is exact,
     thus no double-rounding can occur:
     (a) either the ternary value res is non-zero, and it is the correct
         ternary value that we should return
     (b) or the ternary value res is zero, and we should return 0. */

  /* Set sign of x before exp since check_range needs a valid sign */
  (pstr->negative) ? MPFR_SET_NEG (x) : MPFR_SET_POS (x);

  /* DO NOT USE MPFR_SET_EXP. The exp may be out of range! */
  MPFR_SADD_OVERFLOW (exp, exp, ysize_bits,
                      mpfr_exp_t, mpfr_uexp_t,
                      MPFR_EXP_MIN, MPFR_EXP_MAX,
                      goto overflow, goto underflow);
  MPFR_EXP (x) = exp;
  res = mpfr_check_range (x, res, rnd);
  goto end;

 underflow:
  /* This is called when there is a huge overflow
     (Real expo < MPFR_EXP_MIN << __gmpfr_emin */
  if (rnd == MPFR_RNDN)
    rnd = MPFR_RNDZ;
  res = mpfr_underflow (x, rnd, (pstr->negative) ? -1 : 1);
  goto end;

 overflow:
  res = mpfr_overflow (x, rnd, (pstr->negative) ? -1 : 1);

 end:
  MPFR_TMP_FREE (marker);
  return res;
}

static void
free_parsed_string (struct parsed_string *pstr)
{
  mpfr_free_func (pstr->mantissa, pstr->alloc);
}

int
mpfr_strtofr (mpfr_ptr x, const char *string, char **end, int base,
              mpfr_rnd_t rnd)
{
  int res;
  struct parsed_string pstr;

  /* For base <= 36, parsing is case-insensitive. */
  MPFR_ASSERTN (base == 0 || (base >= 2 && base <= 62));

  /* If an error occurred, it must return 0. */
  MPFR_SET_ZERO (x);
  MPFR_SET_POS (x);

  MPFR_STAT_STATIC_ASSERT (MPFR_MAX_BASE >= 62);
  res = parse_string (x, &pstr, &string, base);
  /* If res == 0, then it was exact (NAN or INF),
     so it is also the ternary value */
  if (MPFR_UNLIKELY (res == -1))  /* invalid data */
    res = 0;  /* x is set to 0, which is exact, thus ternary value is 0 */
  else if (res == 1)
    {
      res = parsed_string_to_mpfr (x, &pstr, rnd);
      free_parsed_string (&pstr);
    }
  else if (res == 2)
    res = mpfr_overflow (x, rnd, (pstr.negative) ? -1 : 1);
  MPFR_ASSERTD (res != 3);
#if 0
  else if (res == 3)
    {
      /* This is called when there is a huge overflow
         (Real expo < MPFR_EXP_MIN << __gmpfr_emin */
      if (rnd == MPFR_RNDN)
        rnd = MPFR_RNDZ;
      res = mpfr_underflow (x, rnd, (pstr.negative) ? -1 : 1);
    }
#endif

  if (end != NULL)
    *end = (char *) string;
  return res;
}

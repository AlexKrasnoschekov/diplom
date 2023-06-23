/* asg6r.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int asg6r_c(real *a, real *b, real *x, real *t, integer *is,
        integer *n, integer *l)
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    extern /* Subroutine */ int asp0r_c(real *, real *, real *, real *,
            integer *, integer *, integer *, integer *);

    /* Parameter adjustments */
    --is;
    --t;
    --x;
    --b;
    a_dim1 = *n;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    asp0r_c(&a[a_offset], &b[1], &x[1], &t[1], &is[1], n, n, l);
    return 0;
} /* asg6r_c */

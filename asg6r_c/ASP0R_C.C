/* asp0r.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static real c_b14 = 1.f;

/* Subroutine */ int asp0r_c(real *a, real *b, real *x, real *t, integer *s,
        integer *n, integer *m, integer *p)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal), r_sign(real *, real *);

    /* Local variables */
    static integer c__, d__;
    static real e, f, g;
    static integer h__, i__, j, k;
    static real l;
    static integer o;
    static real r__, cf, bj, bn, cm, bx;


#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]

    /* Parameter adjustments */
    --t;
    --b;
    --s;
    --x;
    a_dim1 = *n;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    if (*p == 1) {
        goto L1;
    }
    goto L18;
L1:
    i__1 = *m;
    for (d__ = 1; d__ <= i__1; ++d__) {
        if (d__ != *n) {
            goto L2;
        }
        t[d__] = a_ref(d__, d__);
        a_ref(d__, d__) = 0.f;
        s[d__] = d__;
        goto L18;
L2:
        c__ = d__;
        e = 0.f;
        i__2 = *m;
        for (j = d__; j <= i__2; ++j) {
            l = 0.f;
            i__3 = *n;
            for (i__ = d__; i__ <= i__3; ++i__) {
                bx = (r__1 = a_ref(i__, j), abs(r__1));
                if (l < bx) {
                    l = bx;
                }
/* L3: */
            }
            bj = 1.f / l;
            f = 0.f;
            i__3 = *n;
            for (i__ = d__; i__ <= i__3; ++i__) {
                bx = a_ref(i__, j) * bj;
                f += bx * bx;
/* L4: */
            }
            l = (real)sqrt(f) * l;
            if (e >= l) {
                goto L5;
            }
            e = l;
            c__ = j;
L5:
            ;
        }
        if (c__ == d__) {
            goto L7;
        }
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
            cf = a_ref(i__, c__);
            a_ref(i__, c__) = a_ref(i__, d__);
            a_ref(i__, d__) = cf;
/* L6: */
        }
L7:
        s[d__] = c__;
        bn = 1.f / e;
        i__2 = *n;
        for (i__ = d__; i__ <= i__2; ++i__) {
            t[i__] = a_ref(i__, d__) * bn;
/* L8: */
        }
        if (t[d__] != 0.f) {
            goto L9;
        }
        t[d__] = 1.f;
        g = -e;
        goto L11;
L9:
        g = (real)r_sign(&c_b14, &t[d__]);
        bn = (r__1 = t[d__], abs(r__1)) + 1.f;
        t[d__] = g * bn;
        bn = 1.f / (real)sqrt(bn);
        i__2 = *n;
        for (i__ = d__; i__ <= i__2; ++i__) {
            t[i__] *= bn;
/* L10: */
        }
        g = -e * g;
L11:
        i__2 = *n;
        for (i__ = d__; i__ <= i__2; ++i__) {
            a_ref(i__, d__) = t[i__];
/* L12: */
        }
        if (d__ == *m) {
            goto L16;
        }
        o = d__ + 1;
        i__2 = *m;
        for (j = o; j <= i__2; ++j) {
            cm = 0.f;
            i__3 = *n;
            for (i__ = d__; i__ <= i__3; ++i__) {
                cm += a_ref(i__, j) * t[i__];
/* L13: */
            }
            i__3 = *n;
            for (i__ = d__; i__ <= i__3; ++i__) {
                a_ref(i__, j) = a_ref(i__, j) - cm * t[i__];
/* L14: */
            }
/* L15: */
        }
L16:
        t[d__] = g;
/* L17: */
    }
L18:
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
        if (j == *n) {
            goto L22;
        }
        r__ = 0.f;
        i__2 = *n;
        for (i__ = j; i__ <= i__2; ++i__) {
            r__ += b[i__] * a_ref(i__, j);
/* L19: */
        }
        i__2 = *n;
        for (i__ = j; i__ <= i__2; ++i__) {
            b[i__] -= r__ * a_ref(i__, j);
/* L20: */
        }
/* L21: */
    }
L22:
    i__ = *m;
    x[i__] = b[i__] / t[i__];
    if (*m == 1) {
        goto L27;
    }
L23:
    h__ = i__;
    --i__;
    r__ = 0.f;
    i__1 = *m;
    for (j = h__; j <= i__1; ++j) {
        r__ += a_ref(i__, j) * x[j];
/* L24: */
    }
    x[i__] = (b[i__] - r__) / t[i__];
    if (i__ > 1) {
        goto L23;
    }
    j = *m;
L25:
    --j;
    k = s[j];
    if (k == j) {
        goto L26;
    }
    r__ = x[j];
    x[j] = x[k];
    x[k] = r__;
L26:
    if (j > 1) {
        goto L25;
    }
L27:
    return 0;
} /* asp0r_c */

#undef a_ref

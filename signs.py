"""File: signs.py"""

import sympy

from sympy import *


#the variables we need
bmm1,bm1,bmm2,bm2,k1,k2,k,n,mup,mup1,mup2,e1,e2,fp1,fp2 = symbols('bbmm1,bm1,bmm2,bm2,k1,k2,k,n,mup,mup1,mup2,e1,e2,fp1,fp2')
	

#the sign for n-valent negative punctures of type 1. If some of the variables are undetermined, they can be passed as a symbol equal to the symbol in the expression
def p2typ1(bmm1_val, bm1_val, k1_val, mup_val, mup1_val):
    expr = bmm1 + bm1 +k1+(n + mup)*(bm1 + 1 + mup1)
    return expr.subs([(bmm1,bmm1_val), (bm1, bm1_val), (k1,k1_val), (mup, mup_val), (mup1, mup1_val)])

#the sign for n-valent negative punctures of type 2. If some of the variables are undetermined, they can be passed as a symbol equal to the symbol in the expression
def p2typ2(bmm1_val, bm1_val, k1_val, e1_val, mup_val):
    expr = bmm1 + bm1 +k1+ e1 + bm1*(n + mup + 1)
    return expr.subs([(bmm1,bmm1_val), (bm1, bm1_val), (k1,k1_val), (e1, e1_val),(mup, mup_val)])

def Y0(bmm1_val, bm1_val, bmm2_val, bm2_val, e1_val, e2_val, n_val, fp2_val, k2_val, mup2_val, k_val, k1_val, fp1_val, mup1_val):

    expr = bmm1 + bmm2 + bm1 + bm2 + e1*(bm2 + e2 + 1) + n*fp2 + k2 + bm1*(mup2 + fp2 + n + 1) + k + k1 + (n+fp1+ mup1)*(1 + mup2 + bm2)

    return expr.subs([(bmm1,bmm1_val), (bm1, bm1_val), (bmm2,bmm2_val), (bm2, bm2_val), (e1, e1_val), (e2, e2_val), (n, n_val), (k1,k1_val), (k2,k2_val), (k,k_val), (fp1, fp1_val), (fp2, fp2_val), (mup1, mup1_val), (mup2, mup2_val)])


#do everything modulo 2 instead

def p2typ1_mod2(bmm1_val, bm1_val, k1_val, mup_val, mup1_val):
    expr = bmm1 + bm1 +k1+(n + mup)*(bm1 + 1 + mup1)
    expr_eval = expr.subs([(bmm1,bmm1_val), (bm1, bm1_val), (k1,k1_val), (mup, mup_val), (mup1, mup1_val)])

    return expr_eval.as_poly(domain = 'GF(2)')


def p2typ2_mod2(bmm1_val, bm1_val, k1_val, e1_val, mup_val):
    expr = bmm1 + bm1 +k1+ e1 + bm1*(n + mup + 1)
    expr_eval = expr.subs([(bmm1,bmm1_val), (bm1, bm1_val), (k1,k1_val), (e1, e1_val),(mup, mup_val)])
    
    return expr_eval.as_poly(domain = 'GF(2)')


def Y0_mod2(bmm1_val, bm1_val, bmm2_val, bm2_val, e1_val, e2_val, n_val, fp2_val, k2_val, mup2_val, k_val, k1_val, fp1_val, mup1_val):

    expr = bmm1 + bmm2 + bm1 + bm2 + e1*(bm2 + e2 + 1) + n*fp2 + k2 + bm1*(mup2 + fp2 + n + 1) + k + k1 + (n+fp1+ mup1)*(1 + mup2 + bm2)


    expr_eval = expr.subs([(bmm1,bmm1_val), (bm1, bm1_val), (bmm2,bmm2_val), (bm2, bm2_val), (e1, e1_val), (e2, e2_val), (n, n_val), (k1,k1_val), (k2,k2_val), (k,k_val), (fp1, fp1_val), (fp2, fp2_val), (mup1, mup1_val), (mup2, mup2_val)])

 
    return expr_eval.as_poly(domain = 'GF(2)')

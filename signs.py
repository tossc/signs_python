"""File: signs.py"""

import sympy

from sympy import *


#the variables we need
bmm1,bm1,bmm2,bm2,k1,k2,k,n,mup,mup1,mup2,e1,e2,fp1,fp2,n_switch,mua,wua,bmm = symbols('bbmm1,bm1,bmm2,bm2,k1,k2,k,n,mup,mup1,mup2,e1,e2,fp1,fp2,n_switch,mua,wua,bmm')
	

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

    if expr_eval.is_constant():
        return expr_eval%2
    else:
        return trunc(expr_eval,2)

def p2typ2_mod2(bmm1_val, bm1_val, k1_val, e1_val, mup_val):
    expr = bmm1 + bm1 +k1+ e1 + bm1*(n + mup + 1)
    expr_eval = expr.subs([(bmm1,bmm1_val), (bm1, bm1_val), (k1,k1_val), (e1, e1_val),(mup, mup_val)])
    
    if expr_eval.is_constant():
        return expr_eval%2
    else:
        return trunc(expr_eval,2)

def Y0_mod2(bmm1_val, bm1_val, bmm2_val, bm2_val, e1_val, e2_val, n_val, fp2_val, k2_val, mup2_val, k_val, k1_val, fp1_val, mup1_val):



    expr = bmm1 + bmm2 + bm1 + bm2 + e1*(bm2 + e2 + 1) + n*fp2 + k2 + bm1*(mup2 + fp2 + n + 1) + k + k1 + (n+fp1+ mup1)*(1 + mup2 + bm2)


    expr_eval = expr.subs([(bmm1,bmm1_val), (bm1, bm1_val), (bmm2,bmm2_val), (bm2, bm2_val), (e1, e1_val), (e2, e2_val), (n, n_val), (k1,k1_val), (k2,k2_val), (k,k_val), (fp1, fp1_val), (fp2, fp2_val), (mup1, mup1_val), (mup2, mup2_val)])

 
    if expr_eval.is_constant():
        return expr_eval%2
    else:
        return trunc(expr_eval,2)

def Y1_mod2(bmm1_val, bm1_val, bmm2_val, bm2_val, e1_val, e2_val, n_val, fp2_val, k2_val, mup2_val, k_val, k1_val, fp1_val, mup1_val):

    expr = bmm1 + bmm2 + bm1 + bm2 + e1*(bm2 + e2 + 1) + n*fp2 + k2 + bm1*(mup2 + fp2 + n + 1) + k + k1 + (n+fp1+ mup1)*(1 + mup2 + bm2) + fp1 + fp2 + 1


    expr_eval = expr.subs([(bmm1,bmm1_val), (bm1, bm1_val), (bmm2,bmm2_val), (bm2, bm2_val), (e1, e1_val), (e2, e2_val), (n, n_val), (k1,k1_val), (k2,k2_val), (k,k_val), (fp1, fp1_val), (fp2, fp2_val), (mup1, mup1_val), (mup2, mup2_val)])

 
    if expr_eval.is_constant():
        return expr_eval%2
    else:
        return trunc(expr_eval,2)

#the sign of a switch
def switch_mod2(k_val, k1_val, fp1_val):

    expr = 1+ k + k1 + fp1 

    expr_eval = expr.subs([(k1,k1_val), (k,k_val), (fp1, fp1_val)])

 
    if expr_eval.is_constant():
        return expr_eval%2
    else:
        return trunc(expr_eval,2)


#for the last gluing

#case 1, when we do not have any conformal variations
def pos_noconf_mod2(n_switch_val, n_val, wua_val, mua_val):

    expr = (n_switch -1)*n_switch/2*(n+1) + 1 + n_switch*(n*wua + wua) + (1 + wua)*(n + mua)

    expr_eval = expr.subs([(n_switch,n_switch_val), (n,n_val), (wua,wua_val), (mua,mua_val)])

    if expr_eval.is_constant():
        return expr_eval%2
    else:
        return trunc(expr_eval,2)

#case 2, when we have a 1-valent positive puncture and conformal variations
def pos_1val_mod2(n_switch_val, n_val, wua_val, mua_val, bmm_val, bm1_val, k1_val):

    expr = (n_switch -1)*n_switch/2*(n+1) + bmm + 1 + n_switch*(n*wua + wua +1) +  wua*(n + mua + bm1 + 1) + (n + mua)*(bm1 + 1) + n + k1

    expr_eval = expr.subs([(n_switch,n_switch_val), (n,n_val), (wua,wua_val), (mua,mua_val), (bmm,bmm_val), (bm1,bm1_val), (k1,k1_val)])

    if expr_eval.is_constant():
        return expr_eval%2
    else:
        return trunc(expr_eval,2)


#case 3, when we have a 2-valent positive puncture 
def pos_2val_mod2(n_val, wua_val, mua_val, bm1_val, k1_val, mup1_val, bm2_val, e1_val, e2_val, k2_val, bmm1_val, bmm2_val):

    expr = mup1*(bm1 + bm2 + mua) + bm1*(n+1) + bm2*(n + mua) + e1*(bm2 + e2) + k1 + k2+ bmm1 + bmm2+ bm2

    expr_eval = expr.subs([(n,n_val), (wua,wua_val), (mua,mua_val),  (bm1,bm1_val), (k1,k1_val), (mup1, mup1_val), (bm2,bm2_val), (e1, e1_val), (e2,e2_val), (k2,k2_val), (bmm1,bmm1_val), (bmm2, bmm2_val)])

    if expr_eval.is_constant():
        return expr_eval%2
    else:
        return trunc(expr_eval,2)

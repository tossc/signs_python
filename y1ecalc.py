import sympy, signs

'''calculates the sign difference between trees with one positive puncture a and one negative puncture b, and the two different trees we get after a particular Legendrian isotopy, and these trees have in addition one end and one Y1-vertex.'''

from sympy import *

#unknown variables
n, wua_dim, wsb_dim, a, b = symbols('n,wua_dim,wsb_dim,a,b')

#sign for the tree with only two punctures
def straight():
    step_1 = signs.pos_noconf_mod2(0,n, wua_dim, a)

    return step_1

#sign for the tree with the end as b1, b as b2
def tree_II():
    step_1 = signs.Y1_mod2(0,0,0,0,1,0,n,n-wua_dim-1,n-wua_dim-1,b,n-wua_dim,n,n,1)
    step_2 = signs.pos_1val_mod2(0,n,wua_dim,a,1,1,n-wua_dim)

    res = (step_1 + step_2).subs(n**2,n)
    
    if res.is_constant():
        return res%2
    else:
        return trunc(res,2)



#sign for the tree with end as b2, b as b1
def tree_III():
    step_1 = signs.Y1_mod2(0,0,0,0,0,1,n,n,n,1,n-wua_dim,n-wua_dim-1,n-wua_dim-1,b)
    step_2 = signs.pos_1val_mod2(0,n,wua_dim,a,1,1,n-wua_dim)

    res = (step_1 + step_2).subs(n**2,n)
    
    if res.is_constant():
        return res%2
    else:
        return trunc(res,2)


#now we calculate the difference if we also include the sigmas and the intersection and end signs. We return expressions depending on the parity of the Maslov indicies of the punctures

s_Y1 = symbols('s_Y1')



def tot_diff_01():
    ''' difference between the unperturbed tree and the tree with the end as first puncture'''
    nu_end = 1
    intersect_diff = n*(wua_dim + 1) + 1
    sigma_diff = s_Y1(1,b)

    tot = straight() + tree_II() + nu_end + intersect_diff + sigma_diff

    tot1 = tot.subs([(b,a), (n**2,n)])

    if tot1.is_constant():
        return tot1%2
    else:
        return trunc(tot1,2)


def tot_diff_02():
    ''' difference between the unperturbed tree and the tree with the end as second puncture'''
    nu_end = 0
    intersect_diff =  1
    sigma_diff = s_Y1(b,1)

    tot = straight() + tree_III() + nu_end + intersect_diff + sigma_diff

    tot1 = tot.subs([(b,a), (n**2,n)])

    if tot1.is_constant():
        return tot1%2
    else:
        return trunc(tot1,2)

#Now we see what different relation we get by varying the value of the parity of mu(a)

def sigma_diff_from_01():
    for s in (0,1):
        p = tot_diff_01().subs(a,s)
        print(s,trunc(p,2))
        

def sigma_diff_from_02():
    for s in (0,1):
        p = tot_diff_02().subs(a,s)
        print(s,trunc(p,2))

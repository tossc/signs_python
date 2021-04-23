
import sympy, signs

'''calculates the sign differences between a straight flow line with positive puncture a and negative puncture b, and the trees we get by a Legendrian isotopy as described in Tobias trees-paper, and also in my long paper about signs of trees. Now with the assumption that we get different signs from the switch depending on if we traverses the switch in the upper sheet or lower, if the negative puncture is of odd Maslov index. In case the negative puncture is of even Maslov index it should not make a difference, but we put that in as a parameter too.'''

from sympy import *

#unknown variables
n, wua_dim, wsb_dim, a, b, k_dim = symbols('n,wua_dim,wsb_dim,a,b,k_dim')

#sign for the tree with only two punctures
def straight():
    step_1 = signs.pos_noconf_mod2(0,n, wua_dim, a)

    return step_1


#sign for the tree labeled Gamma_A in my long paper, and Gamma_1 in the calculations from 31/7 2019. This tree has one Y0-vertex, one end, one switch and one negative puncture b. Consider first the case when it is the lower section that was irotoped, so that the first vertex we encounter following the boundary of the standard domain is b

def gamma_A1():
    
    # first glue b-piece to the switch-piece, using the sign formula for switches
    step_1 = signs.switch_mod2(wsb_dim-1, wsb_dim,wsb_dim)

    #then glue this to the end-piece at the Y_0-vertex.     
    step_2 = signs.Y0_mod2(0,0,0,0,0,1,n,n,n,1,wsb_dim-1,wsb_dim-1,wsb_dim-1,b+1)

    # and finally glue the positive a-piece to the thing
    step_3 = signs.pos_1val_mod2(0,n,wua_dim,a,1,1,wsb_dim-1)

    tot = (step_1 + step_2 + step_3).subs([(b,a), (wsb_dim, n + 1 + wua_dim), (n**2,n)])

    if tot.is_constant():
        return tot%2
    else:
        return trunc(tot,2)

#sign for the tree similar to Gamma_A in my long paper, labeled by Gamma_2 in the calculations from 31/7 2019. This tree has one Y0-vertex, one end, one switch and one negative puncture b. It is now  the upper section that is isotoped, so that the first vertex we encounter following the boundary of the standard domain is the end 

def gamma_A2():
    
    # first glue b-piece to the switch-piece, using the sign formula for switches
    step_1 = signs.switch_mod2(wsb_dim-1, wsb_dim,wsb_dim)

    #then glue this to the end-piece at the Y_0-vertex.     
    step_2 = signs.Y0_mod2(0,0,0,0,1,0,n,wsb_dim-1,wsb_dim-1,b+1,wsb_dim-1,n,n,1)

    # and finally glue the positive a-piece to the thing
    step_3 = signs.pos_1val_mod2(0,n,wua_dim,a,1,1,wsb_dim-1)

    tot = (step_1 + step_2 + step_3).subs([(b,a), (wsb_dim, n + 1 + wua_dim), (n**2,n)])

    if tot.is_constant():
        return tot%2
    else:
        return trunc(tot,2)


#sign for the tree labeled Gamma_B in my long paper, and Gamma_3 in the calculations from 1/8 2019. This tree has one Y0-vertex, one Y1-vertex, two ends, one switch and one negative puncture b. Consider first the case when it is the lower section that is isotoped, so that the first vertex we encounter following the boundary of the standard domain is b

def gamma_B3():
    
    # first glue b-piece to the switch-piece, using the sign formula for switches
    step_1 = signs.switch_mod2(wsb_dim-1, wsb_dim,wsb_dim)

    #print(step_1) 

    #then glue this to the first end-piece at the Y_0-vertex.     
    step_2 = signs.Y0_mod2(0,0,0,0,0,1,n,n,n,1,wsb_dim-1,wsb_dim-1,wsb_dim-1,b+1)

    #print(step_2)

    #then glue this to the secont end-piece at the Y_1-vertex.
    step_3 = signs.Y1_mod2(1,1,0,0,1,1,n,n,n,1,k_dim, wsb_dim - 1, wsb_dim, b+2)

    #print(step_3)

    # and finally glue the positive a-piece to the thing
    step_4 = signs.pos_1val_mod2(0,n,wua_dim,a,2,2,k_dim)

    #print(step_4)

    tot = (step_1 + step_2 + step_3 + step_4).subs([(b,a), (wsb_dim, n + 1 + wua_dim), (n**2,n)])

    if tot.is_constant():
        return tot%2
    else:
        return trunc(tot,2)


#sign for the tree labeled Gamma_B in my long paper, and Gamma_4 in the calculations from 1/8 2019. This tree has one Y0-vertex, one Y1-vertex, two ends, one switch and one negative puncture b. We consider now the case when it is the upper section that is isotoped, so that the first vertex we encounter following the boundary of the standard domain is e2.

def gamma_B4():
    
    # first glue b-piece to the switch-piece, using the sign formula for switches
    step_1 = signs.switch_mod2(wsb_dim-1, wsb_dim,wsb_dim)

    #print(step_1) 

    #then glue this to the first end-piece at the Y_0-vertex.     
    step_2 = signs.Y0_mod2(0,0,0,0,1,0,n,wsb_dim-1,wsb_dim-1,b+1,wsb_dim-1,n,n,1)

    #print(step_2)

    #then glue this to the second end-piece at the Y_1-vertex.
    step_3 = signs.Y1_mod2(0,0,1,1,1,1,n, wsb_dim, wsb_dim-1, b, k_dim,n,n,1)

    #print(step_3)

    # and finally glue the positive a-piece to the thing
    step_4 = signs.pos_1val_mod2(0,n,wua_dim,a,1,2,k_dim)

    #print(step_4)

    tot = (step_1 + step_2 + step_3 + step_4).subs([(b,a), (wsb_dim, n + 1 + wua_dim), (n**2,n)])

    if tot.is_constant():
        return tot%2
    else:
        return trunc(tot,2)



#now we calculate the difference if we also include the sigmas and the intersection and end signs. We return expressions depending on the parity of the Maslov indicies of the punctures

s_Y1,s_Y0, s_switch, l ,u = symbols('s_Y1,s_Y0, s_switch,l,u')



def tot_diff_A1():
    ''' difference between the unperturbed tree and the tree with an Y0, end, switch and b-puncture, where the lower sheet is isotoped'''
    nu_end = 0
    intersect_diff = n + 1
    sigma_diff = s_Y0(b+1,1) + s_switch(b,l)

    tot = straight() + gamma_A1() + nu_end + intersect_diff + sigma_diff

    tot1 = tot.subs([(b,a), (n**2,n)])

    if tot1.is_constant():
        return tot1%2
    else:
        return trunc(tot1,2)


def tot_diff_A2():
    ''' difference between the unperturbed tree and the tree with an Y0, end, switch and b-puncture, where the upper sheet is isotoped. Might be a mistake in intersect_diff by one. Or maybe this tree does not exist'''
    nu_end = 1
    intersect_diff = n + n*wua_dim + 1
    sigma_diff = s_Y0(1,b+1) + s_switch(b,u)

    tot = straight() + gamma_A2() + nu_end + intersect_diff + sigma_diff

    tot1 = tot.subs([(b,a), (n**2,n)])

    if tot1.is_constant():
        return tot1%2
    else:
        return trunc(tot1,2)


def tot_diff_B3():
    ''' difference between the unperturbed tree and the tree with an Y0, Y1, two end, one switch and b-buncture, where the lower sheet is isotoped, meaning that the punctures are coming in the order bee. There seems to be no difference between the cases when the tree is on the upper half of the annulus we get from the cusps of the perturbed sheet or the lower'''
    nu_end = 0
    intersect_diff = n+1
    sigma_diff = s_Y0(b+1,1) + s_switch(b,l) + s_Y1(b,1)

    tot = straight() + gamma_B3() + nu_end + intersect_diff + sigma_diff

    tot1 = tot.subs([(b,a), (n**2,n)])

    if tot1.is_constant():
        return tot1%2
    else:
        return trunc(tot1,2)


def tot_diff_B4():
    ''' difference between the unperturbed tree and the tree with an Y0, Y1, two end, one switch and b-buncture, where the upper sheet is isotoped'''
    nu_end = 0
    intersect_diff = 1
    sigma_diff = s_Y0(1,b+1) + s_switch(b,u) + s_Y1(1,b)

    tot = straight() + gamma_B4() + nu_end + intersect_diff + sigma_diff

    tot1 = tot.subs([(b,a), (n**2,n)])

    if tot1.is_constant():
        return tot1%2
    else:
        return trunc(tot1,2)

#Now we see what different relation we get by varying the value of the parity of mu(a)

def sigma_diff_from_A1():
    for s in (0,1):
        p = tot_diff_A1().subs(a,s)
        print(s,trunc(p,2))
        
def sigma_diff_from_A2():
    for s in (0,1):
        p = tot_diff_A2().subs(a,s)
        print(s,trunc(p,2))

def sigma_diff_from_B3():
    for s in (0,1):
        p = tot_diff_B3().subs(a,s)
        print(s,trunc(p,2))


def sigma_diff_from_B4():
    for s in (0,1):
        p = tot_diff_B4().subs(a,s)
        print(s,trunc(p,2))
        


import sympy, signs


'''calculates the sign differences between trees with Y0-vertices under legendrian movs, as described in the appendix of the short article'''

from sympy import *


#unknown variables
n,b1,b2,b3,wsb1,wsb2,wsb3, k_dim, wua,a = symbols('n,b1,b2,b3,wsb1,wsb2,wsb3,k_dim,wua,a')

#signs for the first tree
def case_I1():

    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,wsb3,wsb3,b3,wsb3+wsb2+n,wsb2,wsb2,b2)

    #print('step 1 \n', step_1, '\n')

    step_2 = signs.Y0_mod2(0,0,1,1,0,0,n,wsb3+wsb2+n+1,wsb3+wsb2+n,b2+b3,k_dim,wsb1,wsb1,b1)

    #print('step 2 \n', step_2, '\n')

    finish = signs.pos_1val_mod2(0,n,wua,a,1,2,k_dim)

    #print('finish \n', finish, '\n')

    #res=(step_1 + step_2 + finish).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,2)])

    return trunc(step_1 + step_2 + finish,2)

#signs for the second tree
def case_I2():

    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,wsb2,wsb2,b2,wsb2+wsb1+n,wsb1,wsb1,b1)

    #print('step 1 \n', step_1, '\n')
    
    step_2 = signs.Y0_mod2(1,1,0,0,0,0,n,wsb3,wsb3,b3,k_dim,wsb1+wsb2+n,wsb1+wsb2+n+1,b1+b2)
 
    #print('step 2 \n', step_2, '\n')

    finish = signs.pos_1val_mod2(0,n,wua,a,2,2,k_dim)

    #print('finish \n', finish, '\n')
    
    return trunc(step_1 + step_2 + finish,2)


def difference():
   # diff = (case_I1() + case_I2()).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,n)])

    diff = case_I1() + case_I2()

    if diff.is_constant():
        return diff%2
    else:
        return trunc(diff,2)

s_Y0, ub1, lb1, ub2, lb2, ub3, lb3 = symbols('s_Y0, ub0, lb1, ub2, lb2, ub3, lb3' )

def tot_diff():
    intersect_diff = (n+1)*(1 + wsb3)
    #sigma_diff = s_Y0(b2,b3) + s_Y0(b1,b2+b3) + s_Y0(b1,b2) + s_Y0(b1+b2,b3)
    sigma_diff = s_Y0(b2,b3,lb3,ub2) + s_Y0(b1,b2+b3,lb3,ub1) + s_Y0(b1,b2,lb2,ub1) + s_Y0(b1+b2,b3,lb3,ub1)
    tot = (difference() + intersect_diff + sigma_diff).expand().subs(n**2,n)

    #tot = difference() + intersect_diff + sigma_diff

    #tot = difference() + n + 1 + n*wsb3 + wsb3 + s_Y0(b2,b3) + s_Y0(b1,b2+b3) + s_Y0(b1,b2) + s_Y0(b1+b2,b3)
 
    #if tot.is_constant:
    #    return tot%2
    #else:
    return trunc(tot.subs(n**2,n),2)

#Now we see what different relations we get for the sigmas when we vary the values for b1,b2,b3

def sigma_diff():
    p = tot_diff()
    for s1 in (0,1):
        p1 = p.subs(ub1,s1)
        for s2 in (0,1):
            p2 = p1.subs([(lb1,s2),(ub2,s2), (b1,s1 + s2)])
            for s3 in (0,1):
                p3 = p2.subs([(lb2,s3), (ub3,s3), (b2,s2 + s3)])
                for s4 in (0,1):
                    p4 = p3.subs([(lb3,s4), (b3,s3 + s4)])
                    #p1 = p.subs([(b1,s1), (b2,s2), (b3,s3)])
                    p5 = trunc(p4,2)
                #p2 = p1.subs(n,0)
                #return s1,s2,s3,p1
                    print(s1, s2,s3,s4,p5)

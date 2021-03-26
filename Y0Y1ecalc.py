
import sympy, signs


'''calculates the sign differences between trees with 2 Y0-vertices and the same configuration but where we have introduced an additional end and Y1-vertex via a Legendrian move'''

from sympy import *


#unknown variables
n,b1,b2,b3,wsb1,wsb2,wsb3, k_dim, wua,a,vk_dim,k1_dim = symbols('n,b1,b2,b3,wsb1,wsb2,wsb3,k_dim,wua,a,vk_dim,k1_dim')

#signs for the first tree, without end and Y1
def case_I1():

    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,wsb3,wsb3,b3,wsb3+wsb2+n,wsb2,wsb2,b2)

    #print('step 1 \n', step_1, '\n')

    step_2 = signs.Y0_mod2(0,0,1,1,0,0,n,wsb3+wsb2+n+1,wsb3+wsb2+n,b2+b3,k_dim,wsb1,wsb1,b1)

    #print('step 2 \n', step_2, '\n')

    finish = signs.pos_1val_mod2(0,n,wua,a,1,2,k_dim)

    #print('finish \n', finish, '\n')

    #res=(step_1 + step_2 + finish).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,2)])

    return trunc(step_1 + step_2 + finish,2)


#signs for the first tree, with end and Y1 inserted just before the positive puncture,the end coming as the first puncture along the boundary of the corresponding disk
def case_I1_pert1a():

    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,wsb3,wsb3,b3,wsb3+wsb2+n,wsb2,wsb2,b2)

    #print('step 1 \n', step_1, '\n')

    step_2 = signs.Y0_mod2(0,0,1,1,0,0,n,wsb3+wsb2+n+1,wsb3+wsb2+n,b2+b3,k_dim,wsb1,wsb1,b1)

    #print('step 2 \n', step_2, '\n')

    step_3 = signs.Y1_mod2(0,0,1,2,1,0,n,vk_dim,k_dim,a,k1_dim,n,n,1)

    finish = signs.pos_1val_mod2(0,n,wua,a,1,3,k1_dim)

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


#signs for the second tree, with end and Y1 just before the positive puncture, the end being the first puncture
def case_I2_pert1a():

    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,wsb2,wsb2,b2,wsb2+wsb1+n,wsb1,wsb1,b1)

    #print('step 1 \n', step_1, '\n')
    
    step_2 = signs.Y0_mod2(1,1,0,0,0,0,n,wsb3,wsb3,b3,k_dim,wsb1+wsb2+n,wsb1+wsb2+n+1,b1+b2)
 
    step_3 = signs.Y0_mod2(0,0,2,2,1,0,n,vk_dim,k_dim,a,k1_dim,n,n,1)

    #print('step 2 \n', step_2, '\n')

    finish = signs.pos_1val_mod2(0,n,wua,a,1,3,k1_dim)

    #print('finish \n', finish, '\n')
    
    return trunc(step_1 + step_2 + finish,2)



def difference():
   # diff = (case_I1() + case_I2()).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,n)])

    diff = case_I1() + case_I2()

    if diff.is_constant():
        return diff%2
    else:
        return trunc(diff,2)



def difference_pert1a():
   # diff = (case_I1() + case_I2()).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,n)])

    diff = case_I1_pert1a() + case_I2_pert1a()

    if diff.is_constant():
        return diff%2
    else:
        return trunc(diff,2)



s_Y0 = symbols('s_Y0')

def tot_diff():
    #intersect_diff = (n+1)*(1 + wsb3)
    #sigma_diff = s_Y0(b2,b3) + s_Y0(b1,b2+b3) + s_Y0(b1,b2) + s_Y0(b1+b2,b3)
    #tot = (difference() + intersect_diff + sigma_diff).expand().subs(n**2,n)

    #tot = difference() + intersect_diff + sigma_diff

    tot = difference() + n + 1 + n*wsb3 + wsb3 + s_Y0(b2,b3) + s_Y0(b1,b2+b3) + s_Y0(b1,b2) + s_Y0(b1+b2,b3)
 
    #if tot.is_constant:
    #    return tot%2
    #else:
    return trunc(tot.subs(n**2,n),2)

#Now we see what different relations we get for the sigmas when we vary the values for b1,b2,b3

def sigma_diff():
    for s1 in (0,1):
        for s2 in (0,1):
            for s3 in (0,1):
                p = tot_diff()
                p1 = p.subs([(b1,s1), (b2,s2), (b3,s3)])
                p2 = trunc(p1,2)
                #p2 = p1.subs(n,0)
                #return s1,s2,s3,p1
                print(s1, s2,s3,p2)

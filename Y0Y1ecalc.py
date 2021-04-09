
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
#first and second step is similar to the unperturbed tree
    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,wsb3,wsb3,b3,wsb3+wsb2+n,wsb2,wsb2,b2)

    #print('step 1 \n', step_1, '\n')

    step_2 = signs.Y0_mod2(0,0,1,1,0,0,n,wsb3+wsb2+n+1,wsb3+wsb2+n,b2+b3,k_dim,wsb1,wsb1,b1)

    #print('step 2 \n', step_2, '\n')
#this is the new gluing
    step_3 = signs.Y1_mod2(0,0,1,2,1,0,n,n+wua+1,k_dim,a,k1_dim,n,n,1)

    finish = signs.pos_1val_mod2(0,n,wua,a,1,3,k1_dim)

    #print('finish \n', finish, '\n')

    #res=(step_1 + step_2 + finish).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,2)])

    return trunc(step_1 + step_2 + step_3 + finish,2)



#signs for the first tree, with end and Y1 inserted just before the positive puncture,the end coming as the last puncture along the boundary of the corresponding disk
def case_I1_pert1b():
#first and second step is similar to the unperturbed tree
    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,wsb3,wsb3,b3,wsb3+wsb2+n,wsb2,wsb2,b2)

    #print('step 1 \n', step_1, '\n')

    step_2 = signs.Y0_mod2(0,0,1,1,0,0,n,wsb3+wsb2+n+1,wsb3+wsb2+n,b2+b3,k_dim,wsb1,wsb1,b1)

#print('step 2 \n', step_2, '\n')

#this is the new gluing
    step_3 = signs.Y1_mod2(1,2,0,0,0,1,n,n,n,1,k1_dim,k_dim,n+wua+1,a)

    finish = signs.pos_1val_mod2(0,n,wua,a,3,3,k1_dim)

    #print('finish \n', finish, '\n')

    #res=(step_1 + step_2 + finish).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,2)])

    return trunc(step_1 + step_2 + step_3 + finish,2)



#signs for the second tree, first the unperturbed one
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
#first two gluings are the same as in the unperturbed case
    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,wsb2,wsb2,b2,wsb2+wsb1+n,wsb1,wsb1,b1)

    #print('step 1 \n', step_1, '\n')
    
    step_2 = signs.Y0_mod2(1,1,0,0,0,0,n,wsb3,wsb3,b3,k_dim,wsb1+wsb2+n,wsb1+wsb2+n+1,b1+b2)
 
    step_3 = signs.Y1_mod2(0,0,2,2,1,0,n,n+wua+1,k_dim,a,k1_dim,n,n,1)

    #print('step 2 \n', step_2, '\n')

    finish = signs.pos_1val_mod2(0,n,wua,a,1,3,k1_dim)

    #print('finish \n', finish, '\n')
    
    return trunc(step_1 + step_2 + step_3 + finish,2)


#signs for the second tree, with end and Y1 just before the positive puncture, the end being the last puncture
def case_I2_pert1b():
#first two gluings are the same as in the unperturbed case
    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,wsb2,wsb2,b2,wsb2+wsb1+n,wsb1,wsb1,b1)

    #print('step 1 \n', step_1, '\n')
    
    step_2 = signs.Y0_mod2(1,1,0,0,0,0,n,wsb3,wsb3,b3,k_dim,wsb1+wsb2+n,wsb1+wsb2+n+1,b1+b2)

#this is the new gluing
    step_3 = signs.Y1_mod2(2,2,0,0,0,1,n,n,n,1,k1_dim,k_dim,n+wua+1,a)

    finish = signs.pos_1val_mod2(0,n,wua,a,3,3,k1_dim)

    #print('finish \n', finish, '\n')

    #res=(step_1 + step_2 + finish).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,2)])

    return trunc(step_1 + step_2 + step_3 + finish,2)

def difference():
   # diff = (case_I1() + case_I2()).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,n)])

    diff = case_I1() + case_I2()

    if diff.is_constant():
        return diff%2
    else:
        return trunc(diff,2)



def difference_pert1a():
#this gives the difference between the two perturbed trees, where the perturbation is made just before the positive puncture, and the end is the first puncture. We want the same difference as we get from the difference between the unperturbed trees
    
    
    # diff = (case_I1() + case_I2()).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,n)])

    diff = case_I1_pert1a() + case_I2_pert1a()

    if diff.is_constant():
        return diff%2
    else:
        return trunc(diff,2)


def difference_pert1b():
#this gives the difference between the two perturbed trees, where the perturbation is made just before the positive puncture, and the end is the last puncture. We want the same difference as we get from the difference between the unperturbed trees
    
    
    # diff = (case_I1() + case_I2()).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,n)])

    diff = case_I1_pert1b() + case_I2_pert1b()

    if diff.is_constant():
        return diff%2
    else:
        return trunc(diff,2)

s_Y0,s_Y1 = symbols('s_Y0,s_Y1')

def rel_difference_pert1a():
#this gives the difference between the perturbed and the unperturbed trees, where the perturbation is made just before the positive puncture, and the end is the first puncture. We have taken into account the difference in intersection signs
    
    
#for case I1
    stab_diff = case_I1_pert1a() + case_I1()
    intersect_diff = n*wua + n 
    sigma_diff = s_Y1(1,a)
    end_diff = 1 

    diff = stab_diff + intersect_diff + sigma_diff + end_diff

    if diff.is_constant():
        print('For the first tree we get the difference',  diff%2, '\n')
    else:
        print('For the first tree we get the difference',  trunc(diff.subs(n**2,n),2), '\n')

#for case I2
    stab_diff = case_I2_pert1a() + case_I2()
    intersect_diff = n*wua + n 
    sigma_diff = s_Y1(1,a)
    end_diff = 1 

    diff = stab_diff + intersect_diff + sigma_diff + end_diff

    if diff.is_constant():
        print('For the second tree we get the difference',  diff%2, '\n')
    else:
        print('For the second tree we get the difference',  trunc(diff.subs(n**2,n),2), '\n')


def rel_difference_pert1b():
#this gives the difference between the perturbed and the unperturbed trees, where the perturbation is made just before the positive puncture, and the end is the last puncture. We have taken into account the difference in intersection signs
    
    
#for case I1
    stab_diff = case_I1_pert1b() + case_I1()
    intersect_diff = 0
    sigma_diff = s_Y1(a,1)
    end_diff = 0

    diff = stab_diff + intersect_diff + sigma_diff + end_diff

    if diff.is_constant():
        print('For the first tree we get the difference',  diff%2, '\n')
    else:
        print('For the first tree we get the difference',  trunc(diff.subs(n**2,n),2), '\n')

#for case I2
    stab_diff = case_I2_pert1b() + case_I2()
    intersect_diff = 0
    sigma_diff = s_Y1(a,1)
    end_diff = 0

    diff = stab_diff + intersect_diff + sigma_diff + end_diff

    if diff.is_constant():
        print('For the second tree we get the difference',  diff%2, '\n')
    else:
        print('For the second tree we get the difference',  trunc(diff.subs(n**2,n),2), '\n')

def tot_diff_1a():
    '''This function returns the total difference between the unperturbed tree and the perturbed tree where the end comes as first puncture. From rel_difference it follows that we get the same difference from tree I1 and I2, thus we can use the data from tree I1. ''' 
    
    stab_diff = case_I1_pert1a() + case_I1()
    intersect_diff = n*wua + n 
    sigma_diff = s_Y1(1,a)
    end_diff = 1 

    diff = stab_diff + intersect_diff + sigma_diff + end_diff

    if diff.is_constant():
        return diff%2
    else:
        return  trunc(diff.subs(n**2,n),2)


 
def tot_diff_1b():
    '''This function returns the total difference between the unperturbed tree and the perturbed tree where the end comes as last puncture. From rel_difference it follows that we get the same difference from tree I1 and I2, thus we can use the data from tree I1. ''' 
    
    stab_diff = case_I1_pert1b() + case_I1()
    intersect_diff = 0
    sigma_diff = s_Y1(a,1)
    end_diff = 0

    diff = stab_diff + intersect_diff + sigma_diff + end_diff

    if diff.is_constant():
        return diff%2
    else:
        return  trunc(diff.subs(n**2,n),2)

#Now we see what different relations we get for the sigmas when we vary the values for a

#First we write function which takes an expression containing variable a and susbstitute this with a value val

def a_subs(expr,val):
    expr_subs = expr.subs(a,val)
    if expr_subs.is_constant():
        return expr_subs%2
    else:
        return trunc(expr_subs,2)


def sigma_diff():
    for s in (0,1):
        print('|a| = ',s, 'perturbation 1a' , a_subs(tot_diff_1a(),s))
        print('|a| = ',s, 'perturbation 1b' , a_subs(tot_diff_1b(),s))

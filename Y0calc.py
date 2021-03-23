import sympy, signs


'''calculates the sign differences between trees with Y0-vertices under legendrian movs, as described in the appendix of the short article'''

from sympy import *


#unknown variables
n,b1,b2,b3,wsb1,wsb2,wsb3, k_dim, wua,a = symbols('n,b1,b2,b3,wsb1,wsb2,wsb3,k_dim,wua,a')

#signs for the first tree
def case_I1():

    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,wsb3,wsb3,b3,wsb3+wsb2+n,wsb2,wsb2,b2)

    print('step 1 \n', step_1, '\n')

    step_2 = signs.Y0_mod2(0,0,1,1,0,0,n,wsb3+wsb2+n+1,wsb3+wsb2+n,b2+b3,k_dim,wsb1,wsb1,b1)

    print('step 2 \n', step_2, '\n')

    finish = signs.pos_1val_mod2(0,n,wua,a,1,2,k_dim)

    print('finish \n', finish, '\n')

    #res=(step_1 + step_2 + finish).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,2)])

    return trunc(step_1 + step_2 + finish,2)

#signs for the second tree
def case_I2():

    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,wsb2,wsb2,b2,wsb2+wsb1+n,wsb1,wsb1,b1)

    print('step 1 \n', step_1, '\n')
    
    step_2 = signs.Y0_mod2(1,1,0,0,0,0,n,wsb3,wsb3,b3,k_dim,wsb1+wsb2+n,wsb1+wsb2+n+1,b1+b2)
 
    print('step 2 \n', step_2, '\n')

    finish = signs.pos_1val_mod2(0,n,wua,a,2,2,k_dim)

    print('finish \n', finish, '\n')
    
    return trunc(step_1 + step_2 + finish,2)


def difference():
    diff = (case_I1() + case_I2()).subs([(a,b1+b2+b3), (wua,3*n+1+wsb1+wsb2+wsb3),(n**2,n)])

    if diff.is_constant():
        return diff%2
    else:
        return trunc(diff,2)

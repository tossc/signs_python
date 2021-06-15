import sympy, signs

'''calculates the sign differences between trees with 2-valent negative vertices under legendrian moves, as shown in pictures from 1 march 2021. We have not included the 1 coming from different orientations of the space of conformal variatons.'''

from sympy import *

#unknown variables
n,b1,b2,b3 = symbols('n,b1,b2,b3')

#signs for the first pair of trees, labeled I in the pictures
def case_I1():
    
    step_1 = signs.p2typ2_mod2(0,0,n,0,b3)
    
    step_2 = signs.Y0_mod2(0,0,1,1,0,0,n,1,0,b2+b3,0,n-1,n-1,b1)

    return step_1 + step_2

def case_I2():
    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,n,n,b2,n-1,n-1,n-1,b1)

    step_2 = signs.p2typ2_mod2(1,1,n-1,0,b3)

    return step_1 + step_2

#sign difference for the first pair of trees. 
def result_1():
    #print( ' The difference for tree I is', (case_I1() - case_I2()))
    return case_I1() - case_I2()


#signs for the second pair of trees, labeled II in the pictures
def case_II1():
    
    step_1 = signs.p2typ2_mod2(0,0,n,0,b3)
    
    step_2 = signs.p2typ1_mod2(1,1,0,b1,b2+b3)

    return step_1 + step_2

def case_II2():
    step_1 = signs.p2typ1_mod2(0,0,n,b1,b2)

    step_2 = signs.p2typ2_mod2(1,1,0,0,b3)

    return step_1 + step_2

#sign difference for the second pair of trees. 
def result_2():
    #print( ' The difference for tree II is', (case_II1() - case_II2()))
    return case_II1() - case_II2()


#signs for the third pair of trees, labeled III in the pictures
def case_III1():
    
    step_1 = signs.Y0_mod2(0,0,0,0,0,0,n,n-1,n-1,b3,n-1,n,n,b2)
    
    step_2 = signs.p2typ1_mod2(1,1,n-1,b1,b2+b3)

    return step_1 + step_2

def case_III2():
    step_1 = signs.p2typ1_mod2(0,0,n,b1,b2)

    step_2 = signs.Y0_mod2(1,1,0,0,0,0,n,n-1,n-1,b3,0,0,1,b1+b2)

    return step_1 + step_2

#total sign difference for the third pair of trees, except that we haven't taken into account the difference in Vcon. 
def result_3():
    #print( ' The difference for tree III is', (case_III1() - case_III2()))
    return case_III1() - case_III2()


#signs for the fourth pair of trees, labeled IV in the pictures
def case_IV1():
    
    step_1 = signs.p2typ1_mod2(0,0,n,b2,b3) 
   
    step_2 = signs.Y0_mod2(0,0,1,1,0,0,n,1,0,b2+b3,0,n,n,b1)
    

    return step_1 + step_2

def case_IV2():
    step_1 = signs.p2typ2_mod2(0,0,n,0,b2)

    step_2 = signs.Y0_mod2(1,1,0,0,0,0,n,n,n,b3,0,0,1,b1+b2)

    return step_1 + step_2

#total sign difference for the fourth pair of trees, except that we haven't taken into account the difference in Vcon. 
def result_4():
    #print( ' The difference for tree IV is', (case_IV1() - case_IV2()))
    return case_IV1() - case_IV2()

#now we calculate the difference if we also include the sigmas and the intersection signs. We return expressions depending on the parity of the Maslov indicies of the punctures. The Vcon difference is added to the intersect_diff


s_p2neg, s_Y0, ub1, lb1, ub2, lb2, ub3, lb3 = symbols('s_p2neg, s_Y0, ub1,lb1,ub2,lb2,ub3,lb3')

def tot_diff_1():
    intersect_diff = n+1+1
    sigma_diff = s_p2neg(b3,b2,2,lb3,ub2) + s_p2neg(b3,b1+b2,2,lb3,ub1) + s_Y0(b1,b2+b3,lb3,ub1) + s_Y0(b1,b2,lb2,ub1)
    tot = result_1() + intersect_diff + sigma_diff   
    p = Poly(tot, domain = 'GF(2)')

    return p.subs(n**2,n)



def tot_diff_2():
    intersect_diff = 1
    sigma_diff = s_p2neg(b3,b2,2,lb3,ub2) + s_p2neg(b3,b1+b2,2,lb3,ub1) + s_p2neg(b1,b2+b3,1,lb3,ub1) + s_p2neg(b1,b2,1,lb2,ub1)
    tot = result_2().subs(n,1) + intersect_diff + sigma_diff   
    p = Poly(tot, domain = 'GF(2)')

    return p.subs(n**2,n)

def tot_diff_3():
    intersect_diff = 1
    sigma_diff = s_p2neg(b1,b2+b3,1,lb3,ub1) + s_p2neg(b1,b2,1,lb2,ub1) + s_Y0(b2,b3,lb3,ub2) + s_Y0(b1+b2,b3,lb3,ub1)
    tot = result_3() + intersect_diff + sigma_diff   
    p = Poly(tot, domain = 'GF(2)')

    return p.subs(n**2,n)



def tot_diff_4():
    intersect_diff = 0
    sigma_diff = s_p2neg(b2,b3,1,lb3,ub1) + s_p2neg(b2,b1,2,lb2,ub1) + s_Y0(b1,b2+b3,lb3,ub1) + s_Y0(b1+b2,b3,lb3,ub1)
    tot = result_4() + intersect_diff + sigma_diff   
    p = Poly(tot, domain = 'GF(2)')

    return p.subs(n**2,n)


#Now we see what different relations we get for the sigmas when we vary the values for b1,b2,b3

def sigma_diff_from_1():
    p = tot_diff_1()
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

    #for s1 in (0,1):
     #   for s2 in (0,1):
      #      for s3 in (0,1):
       #         p1 = p.subs([(b1,s1), (b2,s2), (b3,s3)])
        #        p2 = trunc(p1,2)
         #       #p2 = p1.subs(n,0)
          #      #return s1,s2,s3,p1
           #     print(s1, s2,s3,p2)



def sigma_diff_from_2():
    p = tot_diff_2()
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



def sigma_diff_from_3():
    p = tot_diff_3()
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



def sigma_diff_from_4():
    p = tot_diff_4()
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
                    print(s1, s2,s3,s4,p5)

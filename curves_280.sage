# In this file we verify the correctness of the table of Section 4.5 of the article ENDOMORPHISM ALGEBRAS OF GEOMETRICALLY SPLIT GENUS 2 JACOBIANS OVER Q
x = QQ['x'].gen()
R.<x> = PolynomialRing(QQ)
# These are the Igusa-Clebsch invariants of the four genus 2 curves whose endomorphism algebra is M_2(Q(sqrt{-280})
list1 = [(8, 930077308872000/392189741895961, 39766891191113246400000/7766839758954136428541, 12269384610004081158912000000/3046074880411510714060466826533022901),
         ( 8, 2014876866302125056/154113765252174121,   42534277501280985644218638336/1512523584750862050519464525,1187580791506718351402800188163111781727731712/29137588084833905384686565664211273122819690625),
         ( 8, 1918868399919702000/29496637472153037481 , -2868910682451370697198520000/160198618224437412060058379371  ,70785331219835183108019840011782156332000000/4725320565306079067682171301646060745203880204451 ),
          ( 8, 737413625842542000/229236961684206961 , 610455144068957687671512000/109755670067550586025891209,-601665387889698876386708259413868000000/25160056333899554523781819246319940426505849 )
]

K.<gK> = NumberField(x^2 + 280)

P.<r,s,t>=QQ[]
for igs in list1:
    print('Igusa--Clebsch invariants: ')
    print(igs)
    MC = Mestre_conic(igs)
    MCK = MC.base_extend(K)
    print('Does the Mestre conic have a point over Q?', magma(MC).HasRationalPoint() )
    print('Does the Mestre conic have a point over Q(sqrt(-280))? ', magma(MCK).HasRationalPoint())
    # we loop looking for (n,n)-isogenies, for 2<= n <= 11
    # these are the igusa invariants provided by Kumar, we extracted them from the files that can be downloaded from the source files in the arxiv version of the article https://arxiv.org/abs/1412.2849 
    for n in [3,6,8]:
        if n == 3:
            I2 = t^2*8*(s^2 + 4*s - 2)
            I4 = t^4*4*(s^4 - 4*s^3 + 6*s^2 - 48*r*s - 4*s + 192*r + 1)
            I6 = t^6*8*(s^6 + 2*s^5 - 21*s^4 - 40*r*s^3 + 44*s^3 + 144*r*s^2 - 41*s^2 + 792*r*s + 18*s - 288*r^2 - 320*r - 3)
            I10 = t^10*2^(14)*r^3
        elif n == 6:
            A1 = 6*(36*s-r^2)
            A = -3*(r+6)^2*(144*r^2*s^4+1728*r*s^4+5184*s^4+5472*r^3*s^3+34560*r^2*s^3+10368*r*s^3 +3672*r^4*s^2+19296*r^3*s^2+4320*r^2*s^2+600*r^5*s+2496*r^4*s +288*r^3*s+25*r^6+60*r^5+36*r^4)
            B1 = -2*(3888*r^2*s^3-233280*r*s^3+139968*s^3-1404*r^4*s^2-12960*r^3*s^2 -190512*r^2*s^2-139968*r*s^2+72*r^6*s+972*r^5*s+6156*r^4*s -7776*r^3*s-11664*r^2*s-r^8-18*r^7-135*r^6-324*r^5-972*r^4)/3
            B  = 18*(r+6)^4 *(192*r^2*s^6+2304*r*s^6+6912*s^6-16704*r^3*s^5-96768*r^2*s^5+20736*r*s^5 -59184*r^4*s^4-144576*r^3*s^4+19008*r^2*s^4-27744*r^5*s^3-54912*r^4*s^3+5760*r^3*s^3-5676*r^6*s^2-9360*r^5*s^2 +720*r^4*s^2-564*r^7*s-576*r^6*s-48*r^5*s-21*r^8+12*r^7-4*r^6)
            B2 = 6144*r^7*(r+6)^6*(6*s+r)^3
            I2, I4, I6, I10 = [t^2*(-24*B1/A1),t^4*( -12*A),t^6*( 96*(A/A1)*B1-36*B), t^10*(-4*A1*B2)]
        elif n == 8:
            A1 = -64*(s-2)^2*s^7*(s+2*r)^8*(s^2+4*r)^2*(s^2+2*s+4*r)*(s^3+2*s^2+12*r*s+8*r^2)
            A = -4*(49*s^10+264*r*s^9+224*s^9+520*r^2*s^8+1224*r*s^8+256*s^8+480*r^3*s^7 +2624*r^2*s^7+1536*r*s^7+208*r^4*s^6+2880*r^3*s^6+3280*r^2*s^6 +1696*r^4*s^5+2688*r^3*s^5+384*r^5*s^4-128*r^6*s^3-1152*r^5*s^3 -512*r^6*s^2+64*r^8)/3
            B1 =  64*s^7*(s+2*r)^8*(s^2+2*s+4*r)*(3*s^15-2*s^14+84*r*s^13+20*s^13+32*r^2*s^12+80*r*s^12+120*s^12 +1248*r^2*s^11+784*r*s^11+112*s^11+1024*r^3*s^10+1984*r^2*s^10 +3072*r*s^10-160*s^10+272*r^4*s^9+9600*r^3*s^9+14464*r^2*s^9 -448*r*s^9-192*s^9+9312*r^4*s^8+29952*r^3*s^8+16128*r^2*s^8 -4352*r*s^8+384*s^8+3136*r^5*s^7+58304*r^4*s^7+80896*r^3*s^7 -28928*r^2*s^7+5376*r*s^7+128*r^6*s^6+51200*r^5*s^6+157824*r^4*s^6 -55296*r^3*s^6+28160*r^2*s^6+17664*r^6*s^5+209152*r^5*s^5 +5120*r^4*s^5+70656*r^3*s^5+1024*r^7*s^4+153600*r^6*s^4  +119808*r^5*s^4+116736*r^4*s^4+48128*r^7*s^3+189440*r^6*s^3 +159744*r^5*s^3+2048*r^8*s^2+139264*r^7*s^2+141312*r^6*s^2  +40960*r^8*s+61440*r^7*s+8192*r^8) /3
            B  = -8*(27*s^16+216*r*s^15+578*s^15+648*r^2*s^14+4032*r*s^14+4272*s^14 +864*r^3*s^13+9408*r^2*s^13+32328*r*s^13+10752*s^13+432*r^4*s^12 +5184*r^3*s^12+96048*r^2*s^12+91008*r*s^12+8624*s^12-10848*r^4*s^11 +132192*r^3*s^11+323136*r^2*s^11+80640*r*s^11-16128*r^5*s^10 +56328*r^4*s^10+628992*r^3*s^10+316416*r^2*s^10-5888*r^6*s^9 -55872*r^5*s^9+739200*r^4*s^9+684288*r^3*s^9-58176*r^6*s^8 +558144*r^5*s^8+915456*r^4*s^8-2304*r^7*s^7+324096*r^6*s^7  +829440*r^5*s^7+7296*r^8*s^6+216576*r^7*s^6+552064*r^6*s^6 +140544*r^8*s^5+267264*r^7*s^5+46080*r^9*s^4+86016*r^8*s^4 +3072*r^10*s^3+27648*r^9*s^3+12288*r^10*s^2-1024*r^12)/27
            B2 =  1/4
            I2, I4, I6, I10 = [t^2*(-24*B1/A1), t^4*(-12*A),t^6*( 96*(A/A1)*B1-36*B), t^10*(-4*A1*B2)]
        # these are the invariants of our curve, they are rational numbers
        I_2, I_4, I_6, I_10 = igs
        # we check if there are some values of r, s that give the above igusa invariants, up to multiplication by t in the weighted projective space
        I = P.ideal([I2.numerator() - I_2*I2.denominator() , I4.numerator() - I_4*I4.denominator(), I6.numerator() - I_6*I6.denominator(), I10.numerator() - I_10*I10.denominator()])
        # first we check if there are r and s in Qbar that give the right Igusa invariants
        if int(magma(I).VarietySizeOverAlgebraicClosure()) > 0:
            print('Found (',n,', ',n, ') isogeny!')
            # now we compute them over a concrete number field
            Qbarm = magma.AlgebraicClosure()
            pts_over_qbar = magma(I).Variety(Qbarm)
            pol = R(Qbarm.AbsolutePolynomial().Coefficients().sage())
            H.<g> = NumberField(pol)
            pts = magma(I).Variety(H).sage()
            rr_t = pts[0][0]
            ss_t = pts[0][1]     
            tt_t = pts[0][2]
            rr = H(rr_t.list())
            ss = H(ss_t.list())
            print(H)
            print('Values of r and s that give these Igusa-Clebsch invariants')
            print('r = ',rr)
            print('s = ',ss)
            assert rr.minpoly() == rr_t.minpoly()
            assert ss.minpoly() == ss_t.minpoly()
            # these are the formulas for the sum and product of the j-invariants of the quotient curves. Again, they are extracted from Kumar's files
            if n == 3:
                j1_plus_j2 = (2*ss^9-17*ss^8+64*ss^7-324*rr*ss^6-140*ss^6+1350*rr*ss^5+196*ss^5-2097*rr*ss^4-182*ss^4  +17496*rr^2*ss^3+1368*rr*ss^3+112*ss^3-23328*rr^2*ss^2-162*rr*ss^2-44*ss^2 +9720*rr^2*ss-198*rr*ss+10*ss-314928*rr^3-432*rr^2+63*rr-1)/rr^2
                zsq = (4*ss^6-20*ss^5+41*ss^4-432*rr*ss^3-44*ss^3-216*rr*ss^2+26*ss^2+576*rr*ss-8*ss +11664*rr^2-184*rr+1)
                j1_times_j2 = (ss^4-4*ss^3+6*ss^2+432*rr*ss-4*ss-288*rr+1)^3/rr^3
            if n == 6:
                j1_plus_j2 =-54*(1119744*(rr+6)^5*ss^9 -31104*rr*(rr+6)^4*(rr^2+162*rr-1944)*ss^8+ 31104*rr^2*(rr+6)^3*(4*rr^3+51*rr^2-3888*rr+43740)*ss^7-7776*rr^3*(rr+6)^2*(5*rr^4-726*rr^3+9540*rr^2+36936*rr-2099520)*ss^6-2592*rr^4*(rr+6)*(31*rr^5-777*rr^4-756*rr^3+208008*rr^2-5458752*rr-42515280)*ss^5-216*rr^5*(301*rr^6+1002*rr^5-32472*rr^4-1228608*rr^3-59731344*rr^2-629226144*rr-1836660096)*ss^4-216*rr^6*(50*rr^6-8601*rr^5-256410*rr^4-3320352*rr^3-32122656*rr^2-195570288*rr-459165024)*ss^3-54*rr^8*(5*rr+54)*(143*rr^4+4248*rr^3+61776*rr^2+419904*rr+944784)*ss^2+ 6*rr^11*(5*rr+54)^2*(11*rr+162)*ss+ rr^12*(5*rr+54)^3 )/(rr^11*(6*ss+rr)^4)
                j1_times_j2 =  729*( 144*(rr+6)^2*ss^4 -288*rr*(rr+6)*(rr+114)*ss^3 + 72*rr^2*(11*rr^2+108*rr+540)*ss^2+ 24*rr^3*(rr+18)*(5*rr+54)*ss + rr^4*(5*rr+54)^2 )^3/(rr^12*(6*ss+rr)^6)
                zsq = (36*ss-rr^2)*(1296*rr^2*ss^3+15552*rr*ss^3+46656*ss^3-36*rr^4*ss^2-6048*rr^3*ss^2-11664*rr^2*ss^2             +139968*rr*ss^2+132*rr^5*ss+2484*rr^4*ss+31104*rr^3*ss+104976*rr^2*ss+7*rr^6             +108*rr^5+972*rr^4)
            if n == 8:
                j1_plus_j2 = 32768*(4*ss^27+60*rr*ss^26+8*ss^26+387*rr^2*ss^25+184*rr*ss^25-48*ss^25+1391*rr^3*ss^24       +2174*rr^2*ss^24-720*rr*ss^24-96*ss^24+3000*rr^4*ss^23+15134*rr^3*ss^23             -6740*rr^2*ss^23-2208*rr*ss^23+192*ss^23+3828*rr^5*ss^22+65780*rr^4*ss^22             -42852*rr^3*ss^22-24200*rr^2*ss^22+2880*rr*ss^22+384*ss^22+2480*rr^6*ss^21             +185624*rr^5*ss^21-165856*rr^4*ss^21-184968*rr^3*ss^21+33040*rr^2*ss^21             +8832*rr*ss^21-256*ss^21+144*rr^7*ss^20+346595*rr^6*ss^20             -337376*rr^5*ss^20-1125168*rr^4*ss^20+293712*rr^3*ss^20+84640*rr^2*ss^20             -3840*rr*ss^20-512*ss^20-768*rr^8*ss^19+427358*rr^7*ss^19             -106364*rr^6*ss^19-5278400*rr^5*ss^19+1987456*rr^4*ss^19             +542880*rr^3*ss^19-49088*rr^2*ss^19-11776*rr*ss^19-320*rr^9*ss^18             +339128*rr^8*ss^18+1215760*rr^7*ss^18-18053964*rr^6*ss^18             +9957952*rr^5*ss^18+2316736*rr^4*ss^18-578240*rr^3*ss^18-90496*rr^2*ss^18             +163440*rr^9*ss^17+3311916*rr^8*ss^17-43967336*rr^7*ss^17             +35594976*rr^6*ss^17+5399424*rr^5*ss^17-3898880*rr^4*ss^17             -180608*rr^3*ss^17+42736*rr^10*ss^16+4418268*rr^9*ss^16             -76012872*rr^8*ss^16+89070336*rr^7*ss^16+1828560*rr^6*ss^16             -14510080*rr^5*ss^16+1347328*rr^4*ss^16+4704*rr^11*ss^15             +3466560*rr^10*ss^15-93622456*rr^9*ss^15+155248096*rr^8*ss^15             -26157408*rr^7*ss^15-29090752*rr^6*ss^15+10639360*rr^5*ss^15             +1617952*rr^11*ss^14-82592128*rr^10*ss^14+188526048*rr^9*ss^14             -78815424*rr^8*ss^14-21061888*rr^7*ss^14+35969216*rr^6*ss^14             +418880*rr^12*ss^13-52556992*rr^11*ss^13+160847232*rr^10*ss^13             -111294400*rr^9*ss^13+38234560*rr^8*ss^13+72901760*rr^7*ss^13             +47040*rr^13*ss^12-24236156*rr^12*ss^12+100604800*rr^11*ss^12             -79703296*rr^10*ss^12+121336512*rr^9*ss^12+96260480*rr^8*ss^12             -7944552*rr^13*ss^11+52570656*rr^12*ss^11-10556672*rr^11*ss^11             +155437568*rr^10*ss^11+86613632*rr^9*ss^11-1687504*rr^14*ss^10             +27433728*rr^13*ss^10+34603104*rr^12*ss^10+119947264*rr^11*ss^10             +58505216*rr^10*ss^10-174048*rr^15*ss^9+13906688*rr^14*ss^9             +38752192*rr^13*ss^9+66470016*rr^12*ss^9+37677056*rr^11*ss^9             +5551104*rr^15*ss^8+24100608*rr^14*ss^8+33684992*rr^13*ss^8             +26533952*rr^12*ss^8+1425920*rr^16*ss^7+9608704*rr^15*ss^7             +16784896*rr^14*ss^7+15668096*rr^13*ss^7+172032*rr^17*ss^6             +1678848*rr^16*ss^6+5877760*rr^15*ss^6+5395712*rr^14*ss^6             -504832*rr^17*ss^5+145408*rr^16*ss^5-31232*rr^15*ss^5-359936*rr^18*ss^4             -1171456*rr^17*ss^4-1111040*rr^16*ss^4-64512*rr^19*ss^3-645120*rr^18*ss^3             -579584*rr^17*ss^3-122880*rr^19*ss^2-100352*rr^18*ss^2+16384*rr^20*ss             +20480*rr^19*ss+8192*rr^21+8192*rr^20) /(ss^5*(ss+2*rr)^7*(ss^2+2*ss+4*rr)^8)
                j1_times_j2 =  67108864*(4*ss^10+24*rr*ss^9-16*ss^9+55*rr^2*ss^8-96*rr*ss^8+16*ss^8+60*rr^3*ss^7                -196*rr^2*ss^7+96*rr*ss^7+28*rr^4*ss^6-120*rr^3*ss^6+220*rr^2*ss^6                +136*rr^4*ss^5+288*rr^3*ss^5+264*rr^5*ss^4+240*rr^4*ss^4+112*rr^6*ss^3                +48*rr^5*ss^3-32*rr^6*ss^2+4*rr^8)^3 /(ss^8*(ss+2*rr)^8*(ss^2+2*ss+4*rr)^8)
                zsq = -(ss^3+2*ss^2+12*rr*ss+8*rr^2)*(ss+2*rr+2)*(27*ss^4+54*rr*ss^3-52*ss^3-48*rr*ss^2+44*ss^2+96*rr^2*ss+72*rr*ss-16*rr^3-16*rr^2)
            KH.<xH> = PolynomialRing(H)
            min_pol_j = xH^2 - xH*H(j1_plus_j2) + H(j1_times_j2)
            # we compute j1 and j2
            if min_pol_j.is_irreducible():
                HH.<gHH> = NumberField(min_pol_j)
                j1, j2 = [a[0] for a in min_pol_j.roots(HH)]
            elif len(min_pol_j.roots())==1: # in this case the two j-invariants are the same
                j1 = min_pol_j.roots()[0][0]
                j2 = j1
            else:
                j1, j2 = [a[0] for a in min_pol_j.roots()]
            E1 = EllipticCurve_from_j(j1)
            E2 = EllipticCurve_from_j(j2)
            print('Curve E1 has CM?', E1.has_cm())
            if E1.has_cm():
                print('CM discriminant of E1 = ', E1.cm_discriminant())
            print('Curve E2 has CM?', E2.has_cm())
            if E2.has_cm():
                print('CM discriminant of E2 = ', E2.cm_discriminant())
            break
    print('_____________________________________')


# ###########################################################
# # Output
# ###########################################################
# sage: %runfile curves_280.sage
# Igusa--Clebsch invariants: 
# (8, 930077308872000/392189741895961, 39766891191113246400000/7766839758954136428541, 12269384610004081158912000000/3046074880411510714060466826533022901)
# Does the Mestre conic have a point over Q? false
# Does the Mestre conic have a point over Q(sqrt(-280))?  true
# Found ( 3 ,  3 ) isogeny!
# Number Field in g with defining polynomial x^4 - 103682/59411343*x^2 + 1/3529707677063649
# Values of r and s that give these Igusa-Clebsch invariants
# r =  -2377704960008242957125/32*g^2 + 4149463607322291675/32
# s =  -22754544369/16*g^2 + 39710559/16
# Curve E1 has CM? True
# CM discriminant of E1 =  -280
# Curve E2 has CM? True
# CM discriminant of E2 =  -280
# _____________________________________
# Igusa--Clebsch invariants: 
# (8, 2014876866302125056/154113765252174121, 42534277501280985644218638336/1512523584750862050519464525, 1187580791506718351402800188163111781727731712/29137588084833905384686565664211273122819690625)
# Does the Mestre conic have a point over Q? false
# Does the Mestre conic have a point over Q(sqrt(-280))?  true
# Found ( 6 ,  6 ) isogeny!
# Number Field in g with defining polynomial x^4 + 60561630225789/78514652200*x^2 + 141495739192361398214601/24658202440347859360000
# Values of r and s that give these Igusa-Clebsch invariants
# r =  157029304400/77864953147443*g^2 - 1/207
# s =  -245358288125/856514484621873*g^2 + 7/36432
# Curve E1 has CM? True
# CM discriminant of E1 =  -280
# Curve E2 has CM? True
# CM discriminant of E2 =  -280
# _____________________________________
# Igusa--Clebsch invariants: 
# (8, 1918868399919702000/29496637472153037481, -2868910682451370697198520000/160198618224437412060058379371, 70785331219835183108019840011782156332000000/4725320565306079067682171301646060745203880204451)
# Does the Mestre conic have a point over Q? false
# Does the Mestre conic have a point over Q(sqrt(-280))?  true
# Found ( 6 ,  6 ) isogeny!
# Number Field in g with defining polynomial x^4 - 2523401259407875/84464166906432*x^2 + 245652324986738538567015625/28536781964790410612171882496
# Values of r and s that give these Igusa-Clebsch invariants
# r =  126696250359648/1802429471005625*g^2 - 3/460
# s =  -37539629736192/3965344836212375*g^2 + 7/2277
# Curve E1 has CM? True
# CM discriminant of E1 =  -280
# Curve E2 has CM? True
# CM discriminant of E2 =  -280
# _____________________________________
# Igusa--Clebsch invariants: 
# (8, 737413625842542000/229236961684206961, 610455144068957687671512000/109755670067550586025891209, -601665387889698876386708259413868000000/25160056333899554523781819246319940426505849)
# Does the Mestre conic have a point over Q? false
# Does the Mestre conic have a point over Q(sqrt(-280))?  true
# Found ( 8 ,  8 ) isogeny!
# Number Field in g with defining polynomial x^4 - 607281867005/70600011300864*x^2 + 28663778387478025/179437017444556482342977273856
# Values of r and s that give these Igusa-Clebsch invariants
# r =  60998409763946496/2116393255955*g^2 - 1654653000/6612799
# s =  -5083200813662208/2116393255955*g^2 + 151113348/6612799
# Curve E1 has CM? True
# CM discriminant of E1 =  -280
# Curve E2 has CM? True
# CM discriminant of E2 =  -280
# _____________________________________

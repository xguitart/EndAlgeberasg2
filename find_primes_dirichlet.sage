# This code generates Table 2 of the article ENDOMORPHISM ALGEBRAS OF GEOMETRICALLY SPLIT GENUS 2 JACOBIANS OVER Q
ds = [3,4,7,8,11,19,43,67,163]
results = [] # format [d, dp, p, a]
for i in range(len(ds)):
    dp = ds[i]
    for j in range(i+1, len(ds)):
        d = ds[j]
        if d == 3 or d == 4 or d == 8:
            continue
        if d == 7 and dp == 3:
            continue
        if d == 11 and dp == 3:
            continue
        if d == 7 and dp == 4:
            continue
        K.<rd> = QuadraticField(-d)
        found = 0
        for b in range(1,10000,2):
            pi = 1/2 + b/2*rd
            p = ZZ(pi.norm())
            if p == 2 or p == 3:
                continue
            if p.is_prime() and legendre_symbol(-dp,p) == -1:
                if found == 0:
                    first_p = p
                    first_a = pi.trace()
                found = found + 1
                if found == 2:
                    results.append((d,dp,first_p, first_a, p,pi.trace()))
                    break
        assert found == 2
                
# falten (8,3), (8,7), (7,3) i (11,3)

d = 8
dp = 3
K.<rd> = QuadraticField(-d)
found = 0
for b in range(1,1000,1):
    pi = 3 + b*rd
    p = ZZ(pi.norm())
    if p == 2 or p == 3:
        continue
    if p.is_prime() and legendre_symbol(-dp,p) == -1:
        if found == 0:
            first_p = p
            first_a = pi.trace()
        found = found + 1
        if found == 2:
            results.append((d,dp,first_p, first_a, p,pi.trace()))
            break
assert found == 2


d = 8
dp = 7
K.<rd> = QuadraticField(-d)
found = 0
for b in range(1,1000,1):
    pi = 1 + b*rd
    p = ZZ(pi.norm())
    if p == 2 or p == 3:
        continue
    if p.is_prime() and legendre_symbol(-dp,p) == -1:
        if found == 0:
            first_p = p
            first_a = pi.trace()
        found = found + 1
        if found == 2:
            results.append((d,dp,first_p, first_a, p,pi.trace()))
            break
assert found == 2


d = 7
dp = 3
K.<rd> = QuadraticField(-d)
found = 0
for b in range(1,1000,1):
    pi = 1 + b*rd
    p = ZZ(pi.norm())
    if p == 2 or p == 3:
        continue
    if p.is_prime() and legendre_symbol(-dp,p) == -1:
        if found == 0:
            first_p = p
            first_a = pi.trace()
        found = found + 1
        if found == 2:
            results.append((d,dp,first_p, first_a, p,pi.trace()))
            break
assert found == 2

d = 11
dp = 3
K.<rd> = QuadraticField(-d)
found = 0
for b in range(1,1000,2):
    pi = 3/2 + b/2*rd
    p = ZZ(pi.norm())
    if p == 2 or p == 3:
        continue
    if p.is_prime() and legendre_symbol(-dp,p) == -1:
        if found == 0:
            first_p = p
            first_a = pi.trace()
        found = found + 1
        if found == 2:
            results.append((d,dp,first_p, first_a, p,pi.trace()))
            break
assert found == 2

results.sort(key=lambda x: (x[0],x[1]))
for i in range(17):
    res = results[i]
    if i<16:
        res2 = results[i+17]
    if i == 16:
        res2 = [' ',' ',' ',' ',' ',' ']
    print('$(',res[0],',',res[1],')$',' & ',res[2],' & ',res[4],' & ',res[5], ' & &  ', '$(',res2[0],',',res2[1],')$',' & ',res2[2],' & ',res2[4],' & ',res2[5], '\\\\')
    print('\hline ')


# check that the results are correct
for res in results:
    d, dp, p1, a1, p2, a2 = res
    K.<rd> = QuadraticField(-d)
    F.<rdp> = QuadraticField(-dp)
    assert len(K.ideal(p1).factor()) == 2 and len(K.ideal(p2).factor()) == 2
    assert len(F.ideal(p1).factor()) == 1 and len(F.ideal(p2).factor()) == 1
    assert K.elements_of_norm(p1)[0].trace().abs() == a1
    assert K.elements_of_norm(p2)[0].trace().abs() == a2

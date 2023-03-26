import sys, re, math, copy, json, random
sys.path.append("./src")
from constants import *
from pathlib import Path
from sym import Sym
from operator import itemgetter

def erf(x):
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911
    sign = 1
    if x < 0:
        sign = -1
    x = abs(x)
    t = 1 / (1 + (p * x))
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-x * x)
    return sign * y

def cliffsDelta(ns1,ns2):
    if len(ns1) > 128: ns1 = samples(ns1,128)
    if len(ns2) > 128: ns2 = samples(ns2,128)
    n,gt,lt = 0,0,0
    for x in ns1:
        for y in ns2:
            n += 1
            if x > y:
                gt += 1
            if x < y:
                lt += 1
    return abs(lt - gt)/n > options['cliff']

def delta(i, other):
    e, y, z = 1E-32, i, other
    return abs(y.mu - z.mu) / ((e + y.sd ** 2 / y.n + z.sd ** 2 / z.n) ** .5)


def bootstrap(y0, z0, NUM):
    x, y, z, yhat, zhat = NUM(), NUM(), NUM(), [], []
    for y1 in y0:
        x.add(y1)
        y.add(y1)
    for z1 in z0:
        x.add(z1)
        z.add(z1)
        
    xmu, ymu, zmu = x.mu, y.mu, z.mu
    for y1 in y0: yhat.append(y1 - ymu + xmu)
    for z1 in z0: zhat.append(z1 - zmu + xmu)
    tobs = delta(y, z)
    n = 0

    for _ in range(1, options['bootstrap'] + 1):
        i = NUM()
        other = NUM()
        for y in samples(yhat).values():
            i.add(y)
        for z in samples(zhat).values():
            other.add(z)
        if delta(i, other) > tobs:
            n += 1
    return n / options['bootstrap'] >= options['conf']

def RX(t,s): 
    t = sorted(t)
    return {'name' : s or "", 'rank':0, 'n':len(t), 'show':"", 'has':t}

def coerce(s):
    if s == 'true':
        return True
    elif s == 'false':
        return False
    elif s.isdigit():
        return int(s)
    elif '.' in s and s.replace('.', '').isdigit():
        return float(s)
    else:
        return s

def eg(key, str, fun):
    egs[key] = fun
    global help
    help = help + '  -g '+ key + '\t' + str + '\n'

def rint(lo,hi, mSeed = None):
    return math.floor(0.5 + rand(lo,hi, mSeed))


def rand(lo, hi, mSeed = None):
    lo, hi = lo or 0, hi or 1
    global Seed
    Seed = 1 if mSeed else (16807 * Seed) % 2147483647
    return lo + (hi-lo) * Seed / 2147483647

def rnd(n, nPlaces = 3):
    mult = 10**nPlaces
    return math.floor(n * mult + 0.5) / mult

def csv(file, fun):
    t = []
    with open(file, 'r', encoding='utf-8') as file:
        for _, line in enumerate(file):
            row = list(map(coerce, line.strip().split(',')))
            t.append(row)
            fun(row)

def kap(t, fun):
    u = {}
    for v in t:
        k = t.index(v)
        v, k = fun(k,v)
        u[k or len(u)] = v
    return u

def dict_kap(t, fun):
    u = {}
    for k,v in t.items():
        v, k = fun(k,v)
        u[k or len(u)] = v
    return u

def cosine(a,b,c):
    den = 1 if c == 0 else 2*c
    x1 = (a**2 + c**2 - b**2) / den
    x2 = max(0, min(1, x1))
    y  = abs((a**2 - x2**2))**.5
    if isinstance(y, complex):
        print('a', a)
        print('x1', x1)
        print('x2', x2)
    return x2, y

def any(t):
    return t[rint(0, len(t) - 1)]

def many(t, n):
    arr = []
    for index in range(1, n + 1):
        arr.append(any(t))
    return arr

def show(node, what, cols, nPlaces, lvl = 0):
  if node:
    print('|..' * lvl, end = '')
    if not node.get('left'):
        print(node['data'].rows[-1].cells[-1])
    else:
        print(int(rnd(100*node['c'], 0)))
    show(node.get('left'), what,cols, nPlaces, lvl+1)
    show(node.get('right'), what,cols,nPlaces, lvl+1)

def deepcopy(t):
    return copy.deepcopy(t)

def oo(t):
    d = t.__dict__
    d['a'] = t.__class__.__name__
    d['id'] = id(t)
    d = dict(sorted(d.items()))
    print(d)

def merge(rx1, rx2):
    rx3 = RX([], rx1['name'])
    rx3['has'] = rx1['has'] + rx2['has']
    rx3['has'] = sorted(rx3['has'])
    rx3['n'] = len(rx3['has'])
    return rx3

def range_fun(at,txt,lo,hi=None):
    return {'at':at,'txt':txt,'lo':lo,'hi':lo or hi or lo,'y':Sym()}

def extend(range,n,s):
    range['lo'] = min(n, range['lo'])
    range['hi'] = max(n, range['hi'])
    range['y'].add(s)

def value(has,nB = None, nR = None, sGoal = None):
    sGoal,nB,nR = sGoal or True, nB or 1, nR or 1
    b,r = 0,0
    for x,n in has.items():
        if x==sGoal:
            b = b + n
        else:
            r = r + n
    b,r = b/(nB+1/float("inf")), r/(nR+1/float("inf"))
    return b**2/(b+r)


def merge2(col1,col2):
  new = merge(col1,col2)
  if new.div() <= (col1.div()*col1.n + col2.div()*col2.n)/new.n:
    return new
  
def div(t):
  t= t['has'] if t['has'] else t
  return (t[ len(t)*9//10 ] - t[ len(t)*1//10 ]) / 2.56

def mid(t):
    t = t['has'] if t['has'] else t
    n = (len(t) - 1) // 2
    return (t[n] + t[n + 1]) / 2 if len(t) % 2 == 0 else t[n + 1]

def sortRXS(rxs):
    for i, x in enumerate(rxs):
        for j, y in enumerate(rxs):
            if mid(x) < mid(y):
                rxs[j], rxs[i] = rxs[i], rxs[j]
    return rxs

def mergeAny(ranges0):
    def noGaps(t):
        for j in range(1,len(t)):
            t[j]['lo'] = t[j-1]['hi']
        t[0]['lo']  = float("-inf")
        t[len(t)-1]['hi'] =  float("inf")
        return t

    ranges1,j = [],0
    while j <= len(ranges0)-1:
        left = ranges0[j]
        right = None if j == len(ranges0)-1 else ranges0[j+1]
        if right:
            y = merge2(left['y'], right['y'])
            if y:
                j = j+1
                left['hi'], left['y'] = right['hi'], y
        ranges1.append(left)
        j = j+1
    return noGaps(ranges0) if len(ranges0)==len(ranges1) else mergeAny(ranges1)

def samples(t, n = None):
    u = {}
    for i in range(1, (n or len(t)) + 1):
        u[i] = t[random.randint(0, len(t) - 1)]
    return u

def gaussian(mu, sd):
    mu, sd = mu or 0, sd or 1
    sq, pi, log, cos, r = math.sqrt, math.pi, math.log, math.cos, random.random
    return mu + sd * sq(-2 * log(r())) * cos(2 * pi * r())
    
def scottKnot(rxs, NUM):
    def merges(i, j):
        out = RX([], rxs[i]['name'])
        for k in range(i, j + 1):
            out = merge(out, rxs[j])
        return out

    def same(lo, cut, hi):
        l = merges(lo, cut)
        r = merges(cut + 1, hi)
        return cliffsDelta(l['has'], r['has']) and bootstrap(l['has'], r['has'], NUM)

    def recurse(lo, hi, rank):
        b4 = merges(lo, hi)
        best = 0
        cut = None
        for j in range(lo, hi + 1):
            if j < hi:
                l = merges(lo, j)
                r = merges(j + 1, hi)
                now = (l['n'] * (mid(l) - mid(b4)) ** 2 + r['n'] * (mid(r) - mid(b4)) ** 2) / (l['n'] + r['n'])
                if now > best:
                    if abs(mid(l) - mid(r)) >= cohen:
                        cut, best = j, now
        if cut != None and not same(lo, cut, hi):
            rank = recurse(lo, cut, rank) + 1
            rank = recurse(cut + 1, hi, rank)
        else:
            for i in range(lo, hi + 1):
                rxs[i]['rank'] = rank
        return rank
    rxs = sortRXS(rxs)
    cohen = div(merges(0, len(rxs) - 1)) * options['cohen']
    recurse(0, len(rxs) - 1, 1)
    return rxs

def tiles(rxs):
    huge = float('inf')
    lo, hi = huge, float('-inf')
    for rx in rxs:
        lo, hi = min(lo, rx['has'][0]), max(hi, rx['has'][len(rx['has']) - 1])
    for rx in rxs:
        t, u = rx['has'], []

        def of(x, most):
            return int(max(0, min(most, x)))

        def at(x):
            return t[of(len(t) * x // 1, len(t))]

        def pos(x):
            return math.floor(of(options['width'] * (x - lo) / (hi - lo + 1E-32) // 1, options['width']))

        for i in range(0, options['width'] + 1):
            u.append(" ")
        a, b, c, d, e = at(.1), at(.3), at(.5), at(.7), at(.9)
        A, B, C, D, E = pos(a), pos(b), pos(c), pos(d), pos(e)
        for i in range(A, B + 1):
            u[i] = "-"
        for i in range(D, E + 1):
            u[i] = "-"
        u[options['width'] // 2] = "|"
        u[C] = "*"
        x = []
        for i in [a, b, c, d, e]:
            x.append(options['Fmt'].format(i))
        rx['show'] = ''.join(u) + str(x)
    return rxs
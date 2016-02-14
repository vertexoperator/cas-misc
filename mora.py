# -*- coding:utf-8 -*-
try:
  import cPickle as pickle
except:
  import pickle


from numbers import *
from itertools import takewhile
from fractions import Fraction
import copy
import math
from functools import reduce



class MonomialExpression:
    def __init__(self):
        self.varlist = []
    def __eq__(lhs,rhs):
        if isinstance(rhs,lhs.__class__) and lhs.varlist==rhs.varlist:
            return True
        else:
            return False
    def __ne__(lhs,rhs):
        return not(lhs==rhs)
    def __mul__(lhs,rhs):
        ret = MonomialExpression()
        ret.varlist = lhs.varlist + rhs.varlist
        return ret
    def __hash__(self):
        return hash(pickle.dumps(self.varlist))
    def __str__(self):
        res = []
        idx = 0
        for _ in range(len(self.varlist)):
            if idx>len(self.varlist)-1:break
            varname = self.varlist[idx]
            tmp = list(takewhile(lambda x:x==varname , self.varlist[idx:]))
            if len(tmp)==1:
                res.append(varname)
            else:
                res.append( "{0}^{1}".format(varname , len(tmp)) )
            idx+=len(tmp)
        return "*".join(res)
    def __pow__(self,n):
        assert(isinstance(n,int) and n>=0)
        ret = MonomialExpression()
        ret.varlist = self.varlist*n
        return ret
    def deg(self):
        return len(self.varlist)


class Expression(object):
    def __init__(self , val=None):
        if isinstance(val,Number):
            m = MonomialExpression()
            self.coeffs = {m:val}
        elif isinstance(val,MonomialExpression):
            self.coeffs = {val:1}
        else:
            self.coeffs = {}
    def simplify(self):
        keys = [m for m in self.coeffs if self.coeffs[m]==0]
        for m in keys:
            del self.coeffs[m]
    def __add__(lhs,rhs):
        ret = Expression()
        for (m,c) in lhs.coeffs.items():
            ret.coeffs[m] = ret.coeffs.get(m,0)+c
        if isinstance(rhs,Expression):
            for (m,c) in rhs.coeffs.items():
                ret.coeffs[m] = ret.coeffs.get(m,0)+c
        elif isinstance(rhs,Number):
            m = MonomialExpression()
            ret.coeffs[m] = ret.coeffs.get(m,0)+rhs
        ret.simplify()
        return ret
    def __neg__(self):
        ret = Expression()
        for (m,c) in self.coeffs.items():
            ret.coeffs[m] = ret.coeffs.get(m,0)-c
        return ret
    def __sub__(lhs,rhs):
        ret = Expression()
        for (m,c) in lhs.coeffs.items():
            ret.coeffs[m] = ret.coeffs.get(m,0)+c
        if isinstance(rhs,Expression):
            for (m,c) in rhs.coeffs.items():
                ret.coeffs[m] = ret.coeffs.get(m,0)-c
        elif isinstance(rhs,Number):
            m = MonomialExpression()
            ret.coeffs[m] = ret.coeffs.get(m,0)-rhs
        ret.simplify()
        return ret
    def __mul__(lhs,rhs):
        ret = Expression()
        if isinstance(rhs,Expression):
            for (m1,c1) in lhs.coeffs.items():
               for (m2,c2) in rhs.coeffs.items():
                   c3 = c1*c2
                   m3 = m1*m2
                   ret.coeffs[m3] = ret.coeffs.get(m3,0)+c3
        elif isinstance(rhs,Number):
            for (m1,c1) in lhs.coeffs.items():
                ret.coeffs[m1] = c1*rhs
        ret.simplify()
        return ret
    def __div__(lhs,rhs):
        if isinstance(rhs,Number):
           return lhs*Fraction(1,rhs)
        else:
           assert(False),("invalid division {0}/{1}".format(str(lhs),str(rhs)))
    def __rmul__(rhs,lhs):
        if isinstance(lhs,Expression):
            return (lhs*rhs)
        elif isinstance(lhs,Number):
            return (rhs*lhs)
    def __radd__(rhs,lhs):
        if isinstance(lhs,Expression):
            return (lhs+rhs)
        elif isinstance(lhs,Number):
            return (rhs+lhs)
    def __rsub__(rhs,lhs):
        if isinstance(lhs,Expression):
            return (lhs-rhs)
        elif isinstance(lhs,Number):
            return (-rhs+lhs)
    def __eq__(lhs,rhs):
        if isinstance(rhs,Expression):
            return (lhs-rhs==0)
        elif isinstance(rhs,Number) and rhs!=0:
            if len(lhs.coeffs)!=1:return False
            m,c = list(lhs.coeffs.items())[0]
            return (m.deg()==0 and c==rhs)
        elif isinstance(rhs,Number) and rhs==0:
            if len(lhs.coeffs)>1:
                return False
            elif len(lhs.coeffs)==1:
                m,c = list(lhs.coeffs.items())[0]
                return (m.deg()==0 and c==0)
            else:
                assert(len(lhs.coeffs)==0)
                return True
        else:
            return False
    def __pow__(self , n):
        return reduce(lambda x,y:x*y,[self]*n,1)
    def __nonzero__(self):
        return not (self==0)
    def __ne__(lhs , rhs):
        return not (lhs==rhs)
    def __str__(self):
        cs = [(m,c) for (m,c) in self.coeffs.items()]
        if len(cs)==0:return "0"
        cs.sort(key=lambda x:x[0].deg() , reverse=True)
        res = []
        for idx,(m,c) in enumerate(cs):
            if idx==0 and m.deg()==0:
                res.append( str(c) )
            elif idx==0 and c==1:
                res.append( str(m) )
            elif idx==0 and c==-1:
                res.append( "-{0}".format(str(m)) )
            elif idx==0:
                res.append( "{0} {1}".format(str(c) , str(m)) )
            elif m.deg()==0:
                mc = str(c)
                if mc[0]!="-" and mc[0]!="+":
                    res.append("+")
                    res.append(mc)
                else:
                    res.append("{0} {1}".format(mc[0],mc[1:]))
            elif c==1:
                res.append("+")
                res.append(str(m))
            elif c==-1:
                res.append("-")
                res.append(str(m))
            else:
                mc = "{0} {1}".format(str(c) , str(m))
                if mc[0]!="-" and mc[0]!="+":
                    res.append("+")
                    res.append(mc)
                else:
                    res.append("{0} {1}".format(mc[0],mc[1:]))
        return " ".join(res)
    def __repr__(self):
        return str(self)



class Symbol(Expression):
    def __init__(self , varname):
#       assert(isinstance(varname,basestring) and len(varname)>0)
       super(Symbol,self).__init__()
       m = MonomialExpression()
       m.varlist = [varname]
       self.coeffs[m] = 1




"""
>>> x,y,z=Symbol("x"),Symbol("y"),Symbol("z")
>>> [m2i(p.coeffs.keys()[0],[x,y]) for p in [Expression(1),x,y,x*x,x*y,y*x,y*y,x*x*x,x*x*y,x*y*x,x*y*y,y*x*x,y*x*y,y*y*x,y*y*y,x*x*x*x]]
[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

>>> [m2i(p.coeffs.keys()[0],[x,y,z]) for p in [Expression(1),x,y,z,x*x,x*y,x*z,y*x,y*y,y*z,z*x,z*y,z*z,x*x*x,x*x*y,x*x*z]]
[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

"""
def m2i(e , vlist):
    assert(len(vlist)>0)
    assert(isinstance(e,MonomialExpression)),("type error : {0}".format(repr(e)))
    if len(vlist)==1:
       varname = vlist[0].coeffs.keys()[0].varlist[0]
       assert(all([s==varname for s in e.varlist]))
       return len(e.varlist)
    else:
       ret = 0
       N = len(vlist)
       for vname in e.varlist:
          c = vlist.index(Symbol(vname))+1
          ret = ret*N+c
       return ret




def ilog(n , base):
   ret = 0
   n0 = 1
   while True:
      if n0>n:return (ret-1)
      ret+=1
      n0*=base


def ideg(n , Nvar):
   ret = 0
   n0 = 1
   while True:
      if n0>n*(Nvar-1)+1:return (ret-1)
      ret+=1
      n0*=Nvar


def im_decode(n , Nvar):
   ret = []
   deg = ideg(n , Nvar)
   base = ((Nvar**deg)-1)//(Nvar-1)
   rest = n-base
   for _ in range(deg):
      rem,mod = divmod(rest , Nvar)
      ret.append(int(mod))
      rest = (rest-mod)//Nvar
   ret.reverse()
   return ret



def im_encode(ls , Nvar):
    ret = 0
    for val in ls:
       c = val+1
       ret = ret*Nvar+c
    return ret



"""
>>> x,y,z=Symbol("x"),Symbol("y"),Symbol("z")
>>> [m2i(i2m(n,[x,y,z]),[x,y,z]) for n in range(25)]
[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]

"""
def i2m(n , vlist):
    assert(len(vlist)>0)
    assert(n>=0)
    assert(all([isinstance(v,Symbol) for v in vlist]))
    if len(vlist)==1:
       m = MonomialExpression()
       m.varlist = [vlist[0]]*n
       return m
    if n==0:return MonomialExpression()
    N = len(vlist)
    deg = ideg(n , N)
    base = ((N**deg)-1)//(N-1)
    assert(n>=base),(n,base)
    step = ((N**(deg-1))-1)//(N-1)
    rest = n-base
    ret = []
    for _ in range(deg):
        rem,mod = divmod(rest , N)
        ret.append(vlist[int(mod)])
        rest = (rest-mod)//N
    ret.reverse()
    return list(reduce(lambda x,y:x*y , ret , 1).coeffs.keys())[0]
 



def im_mul(n , m , Nvar):
   if Nvar==1:return (n+m)
   deg1 = ideg(n,Nvar)
   deg2 = ideg(m,Nvar)
   p1,p2 = Nvar**deg1 , Nvar**deg2
   base = (p1*p2-1)//(Nvar-1)
   base1 = (p1-1)//(Nvar-1)
   base2 = (p2-1)//(Nvar-1)
   return base + p2*(n-base1) + (m-base2)





def memoize(f):
   cache = {}
   def __fun__(*x):
     if not x in cache:cache[x] = f(*x)
     return cache[x]
   return __fun__


def im_lpop(n , Nvar):
   if n<=Nvar:return 0
   if Nvar==1:return (n-1)
   deg = ideg(n,Nvar)
   p = Nvar**(deg-1)
   base0 = (p-1)//(Nvar-1)
   lsym = (n-base0)//(p)
   return n - lsym*p




def im_rpop(n , Nvar):
   if n<=Nvar:return 0
   if Nvar==1:return (n-1)
   deg = ideg(n,Nvar)
   base = (Nvar**deg-1)//(Nvar-1)
   return (n-base)//Nvar + base//Nvar  #!!-- not equal to (n//Nvar)




"""
非可換多項式の"分散表現"
N変数(非可換)単項式を、(単項式順序を保つように)単一の自然数にencodeして表現する

"""
class DExpression(object):
    __slots__ = ["Nvar" , "coeffs" , "_simplified"]
    def __init__(self , Nvar):
        self.Nvar = Nvar
        self.coeffs = {0:0}
        self._simplified = False
    def simplify(self):
        if self._simplified:return self
        keys = [k for (k,c) in self.coeffs.items() if c==0 and k!=0]
        for k in keys:del self.coeffs[k]
        self._simplified = True
        return self
    def __add__(lhs,rhs):
        ret = DExpression(lhs.Nvar)
        ret.coeffs = copy.deepcopy(lhs.coeffs)
        if isinstance(rhs,DExpression):
            assert(lhs.Nvar==rhs.Nvar)
            for (m,c) in rhs.coeffs.items():
                ret.coeffs[m] = ret.coeffs.get(m,0)+c
        elif isinstance(rhs,Number):
            ret.coeffs[0] += rhs
            ret._simplified = lhs._simplified
        return ret.simplify()
    def __neg__(self):
        ret = DExpression(self.Nvar)
        for (m,c) in self.coeffs.items():
            ret.coeffs[m] = -c
        return ret.simplify()
    def __sub__(lhs,rhs):
        ret = DExpression(lhs.Nvar)
        ret.coeffs = copy.deepcopy(lhs.coeffs)
        if isinstance(rhs,DExpression):
            assert(lhs.Nvar==rhs.Nvar)
            for (m,c) in rhs.coeffs.items():
                ret.coeffs[m] = ret.coeffs.get(m,0)-c
        elif isinstance(rhs,Number):
            ret.coeffs[0] -= rhs
            ret._simplified = lhs._simplified
        return ret.simplify()
    def __isub__(self , rhs):
        for (m,c) in rhs.coeffs.items():
            c1 = self.coeffs.get(m,0) - c
            if m==0 or c1!=0:self.coeffs[m] = c1
            elif m!=0 and (m in self.coeffs):del self.coeffs[m]
        self._simplified = False
        return self
    def __mul__(lhs,rhs):
        ret = DExpression(lhs.Nvar)
        if isinstance(rhs,DExpression):
            assert(lhs.Nvar==rhs.Nvar)
            for (m1,c1) in lhs.coeffs.items():
               for (m2,c2) in rhs.coeffs.items():
                   c3 = c1*c2
                   m3 = im_mul(m1 , m2 , ret.Nvar)
                   ret.coeffs[m3] = ret.coeffs.get(m3,0)+c3
        elif isinstance(rhs,Number):
            for m1 in lhs.coeffs.keys():
                ret.coeffs[m1] = lhs.coeffs[m1]*rhs
        return ret.simplify()
    def __div__(lhs,rhs):
        if isinstance(rhs,Number):
           return lhs*Fraction(1,rhs)
        else:
           assert(False),("invalid division {0}/{1}".format(str(lhs),str(rhs)))
    def __rmul__(rhs,lhs):
        if isinstance(lhs,DExpression):
            return (lhs*rhs)
        elif isinstance(lhs,Number):
            return (rhs*lhs)
    def __radd__(rhs,lhs):
        if isinstance(lhs,DExpression):
            return (lhs+rhs)
        elif isinstance(lhs,Number):
            return (rhs+lhs)
    def __rsub__(rhs,lhs):
        if isinstance(lhs,DExpression):
            return (lhs-rhs)
        elif isinstance(lhs,Number):
            return (-rhs+lhs)
    def __eq__(lhs,rhs):
        if isinstance(rhs,DExpression):
            lhs.simplify()
            rhs.simplify()
            if len(lhs.coeffs)!=len(rhs.coeffs):return False
            for (m,c) in lhs.coeffs.items():
                if not m in rhs.coeffs or rhs.coeffs[m]!=c:return False
            return True
        elif isinstance(rhs,Number):
            return (list(lhs.coeffs.keys())==[0] and lhs.coeffs[0]==rhs)
        else:
            return False
    def __pow__(self,n):
        assert(isinstance(n,int) and n>=0)
        return reduce(lambda x,y:x*y,[self]*n,1)
    def __nonzero__(self):
        return not (self==0)
    def __ne__(lhs , rhs):
        return not (lhs==rhs)



def p2dp(p , varlist):
    ret = DExpression(len(varlist))
    for (m,c) in p.coeffs.items():
        k = m2i(m , varlist)
        ret.coeffs[k] = c
    return ret


def dp2p(p , varlist):
    ret = Expression()
    for (k,c) in p.coeffs.items():
       m = i2m(k , varlist)
       if c!=0:ret.coeffs[m] = c
    return ret





def im_lpops(n , Nvar , k):
   ret = []
   if k==0:return [n]
   if n<=Nvar:return [n,0]
   if Nvar==1:return range(n,max(n-k , 0)-1,-1)[::-1]
   deg = ilog(n*(Nvar-1)+1 , Nvar)
   assert(deg>=k),"im_lpops error"
   r = n
   p = Nvar**(deg-1)
   base0 = (p-1)//(Nvar-1)
   ret.append(r)
   for _ in range(k):
      lsym = (r-base0)//(p)
      r -= lsym*p
      ret.append(r)
      p = p//Nvar
      base0 = (p-1)//(Nvar-1)
   return ret


def im_div(m0 , m1 , Nvar):
    if m0<m1:return None
    if m0==m1:return (0,0)
    deg0 = ilog(m0*(Nvar-1)+1 , Nvar)
    deg1 = ilog(m1*(Nvar-1)+1 , Nvar)
    if deg0==deg1 and m0!=m1:return None
    assert(deg0>deg1)
    rm_list = im_lpops(m0 , Nvar , deg0)
    lm,degL = m0,deg0
    degR0 = deg0
    for _ in range(degL+1):
       if deg1+degL<=deg0:
           k = deg1+degL+degR0 - deg0
           rm = rm_list[k]
           if m0==im_mul(lm , im_mul(m1 ,rm , Nvar) , Nvar):
              return (lm,rm)
       lm = im_rpop(lm , Nvar)
       degL-=1
    return None



def im_div_all(m0 , m1 , Nvar):
    if m0<m1:return []
    if m0==m1:return [(0,0)]
    ret = []
    deg0 = ilog(m0*(Nvar-1)+1 , Nvar)
    deg1 = ilog(m1*(Nvar-1)+1 , Nvar)
    if deg0==deg1 and m0!=m1:return ret
    lm = im_rpop(m0 , Nvar)
    rm0 = im_lpop(m0 , Nvar)
    degL = deg0 - 1
    degR0 = deg0 - 1
    while degL>=0:
       if deg1+degL+degR0<deg0:return None
       if deg1+degL<=deg0 and deg1+degL+degR0>=deg0:
           rm = rm0
           degR = degR0
           while degR>=0:
             if deg1+degL+degR<deg0:break
             if deg1+degL+degR==deg0:
                if m0==im_mul(lm , im_mul(m1 ,rm , Nvar) , Nvar):
                    ret.append( (lm,rm) )
             if rm==0:break
             rm = im_lpop(rm , Nvar)
             degR-=1
       if lm==0:break
       lm = im_rpop(lm , Nvar)
       degL-=1
    return ret



def im_div_with_deg(m0 , m1 , Nvar):
    if m0<m1:return None
    if m0==m1:return (0,0,0,0)
    deg0 = ilog(m0*(Nvar-1)+1 , Nvar)
    deg1 = ilog(m1*(Nvar-1)+1 , Nvar)
    if deg0==deg1 and m0!=m1:return None
    assert(deg0>deg1)
    rm_list = im_lpops(m0 , Nvar , deg0)
    lm,degL = m0,deg0
    degR0 = deg0
    for _ in range(degL+1):
       if deg1+degL<=deg0:
           k = deg1+degL+degR0 - deg0
           rm = rm_list[k]
           if m0==im_mul(lm , im_mul(m1 ,rm , Nvar) , Nvar):
              return (lm,rm,degL,deg0-deg1-degL)
       lm = im_rpop(lm , Nvar)
       degL-=1
    return None


import time
def dp_nf(p , _G):
    def tip(p):
       return max(p.coeffs.keys())
    if len(p.coeffs)==1:return p
    h = copy.deepcopy(p)
    h.simplify()
    r = DExpression(p.Nvar)
    Nvar = p.Nvar
    G = [p.simplify() for p in _G if p!=0]
    HMG = [(p*Fraction(1,p.coeffs[tip(p)]) , tip(p)) for p in G]
    while True:
       m = tip(h)
       if m==0:
           r.coeffs[m] = h.coeffs[m]
           break
       c0 = h.coeffs[m]
       if c0==0:
           del h.coeffs[m]
           continue
       for (p1,m1) in HMG:
#           assert(p1.coeffs[m1]==1),(p1.coeffs[m1],m1)
           if m<m1:continue
           if m==m1:
               h -= p1*c0
               break
           rem = im_div_with_deg(m , m1 , Nvar)
           if rem!=None:
              um,vm,udeg,vdeg = rem
              for (m2,c2) in p1.coeffs.items():
                  deg2 = ideg(m2 , Nvar)
                  pu,pc,pv = Nvar**udeg,Nvar**deg2,Nvar**vdeg
                  m3 = pv*(m2-(pc-1)//(Nvar-1))+(vm-(pv-1)//(Nvar-1))
                  m3 = (pu*pc*pv-1)//(Nvar-1) + (pc*pv)*(um-(pu-1)//(Nvar-1)) + m3
                  #assert(m3==im_mul(um , im_mul(m2, vm,Nvar) , Nvar)),(m2,um,vm,m3)
                  c4 = h.coeffs.get(m3,0) - c0*c2
                  if c4!=0 or m3==0:
                     h.coeffs[m3] = c4
                  elif m3 in h.coeffs:
                     del h.coeffs[m3]
              #assert(im_mul(um , im_mul(m1,vm,Nvar),Nvar)==m)
              #assert(tip(h)<m),(tip(h),m)
              break
       else:
           del h.coeffs[m]
           r.coeffs[m] = c0
    return r





def dp_obs(p0 , p1):
   def tip(p):
      return max(p.coeffs.keys())
   ret = []
   assert(p0.Nvar==p1.Nvar)
   Nvar = p0.Nvar
   m0,m1 = tip(p0),tip(p1)
   if True:
       lm = im_rpop(m0 , Nvar)
       rm0 = im_lpop(m1 , Nvar)
       while lm>0:
           rm = rm0
           while rm>0:
              if im_mul(m0 , rm , Nvar)==im_mul(lm , m1 , Nvar):
                  assert(rm!=m1 and m0!=lm)
                  ret.append( (0 , p0 , rm , lm , p1 , 0) )
              rm = im_lpop(rm , Nvar)
           lm = im_rpop(lm , Nvar)
           rm0 = im_lpop(rm0 , Nvar)
   if m0!=m1:
       lm = im_rpop(m1 , Nvar)
       rm0 = im_lpop(m0 , Nvar)
       while lm>0:
           rm = rm0
           while rm>0:
              if im_mul(m1 , rm , Nvar)==im_mul(lm , m0 , Nvar):
                  assert(m1!=lm and rm!=m0) 
                  ret.append( (lm , p0 , 0 , 0 , p1 , rm) )
              rm = im_lpop(rm , Nvar)
           lm = im_rpop(lm , Nvar)
           rm0 = im_lpop(rm0 , Nvar)
   if m0>m1:  #-- central overlap
       for (lm,rm) in im_div_all(m0 , m1 , Nvar):
          ret.append( (0 , p0 , 0 , lm , p1 , rm) )
   if m0<m1:  #-- central overlap
       for (lm,rm) in im_div_all(m1 , m0 , Nvar):
          ret.append( (lm , p0 , rm , 0 , p1 , 0) )
   return ret




import time
def dp_mora(_G):
    NoSugar = False
    def tdeg(p,Nvar):
        return max([ideg(m,Nvar) for m in p.coeffs.keys()])
    def tip(p):
        return max(p.coeffs.keys())
    G = [q*Fraction(1,q.coeffs[tip(q)]) for q in _G if q!=0]
    B = []
    #-- inter-reduction
    while True:
       ok = False
       for n,p in enumerate(G):
          p2 = dp_nf(p , [q for q in G if q!=p and q!=0])
          if p2==0:
              G[n] *= 0
          elif p2!=p:
              G[n] = p2*Fraction(1 , p2.coeffs[tip(p2)])
              break
          elif n==len(G)-1:
              ok = True
              break
       if ok:break
    G = [p for p in G if p!=0]
    if len(G)==0:return []
    Nvar = G[0].Nvar
    G.sort(key=lambda x:tdeg(x,Nvar))
    masks = [True]*len(G)
    sugars = [tdeg(p,Nvar) for p in G]
    print("{0} interreduced bases".format(len(G)))
    for i,p in enumerate(G):
        for j,q in enumerate(G):
            if i<=j:
               for c in dp_obs(p,q):
                   cm = im_mul(c[0] , im_mul(tip(c[1]) , c[2] , Nvar) , Nvar)
                   if NoSugar:
                       s_obs = 0
                   else:
                       s_obs = max(tdeg(p,Nvar)+ideg(c[0],Nvar)+ideg(c[2],Nvar),
                                   tdeg(q,Nvar)+ideg(c[3],Nvar)+ideg(c[5],Nvar))
                   B.append( (c,cm,s_obs) )
    B.sort(key=lambda x:(x[2],x[1]) , reverse=True)
    Nobs = 0
    while len(B)>0:
        (lm0,p0,rm0,lm1,p1,rm1),_,s_h = B.pop()
        Nobs+=1
        assert(len([x for x in [lm0,rm0,lm1,rm1] if x==0])>=2)
        u0,v0 = DExpression(Nvar),DExpression(Nvar)
        u1,v1 = DExpression(Nvar),DExpression(Nvar)
        u0.coeffs[lm0] = 1
        v0.coeffs[rm0] = 1
        u1.coeffs[lm1] = 1
        v1.coeffs[rm1] = 1
        t0 = time.time()
        p2 = dp_nf(u0*p0*v0 - u1*p1*v1 , G)
        t1 = time.time()
        if t1-t0>5.0:
           print("1 obstruction removed (time:{0:.3f}) (rp={1})".format(t1-t0 , len(B)))
        if p2==0:continue
        p2 = p2*Fraction(1, p2.coeffs[tip(p2)])
        G.append(p2)
        masks.append( True )
        sugars.append( s_h )
        m2 = tip(p2)
        RED = []
        for (c,cm,s) in B:
            cm = im_mul(c[0],im_mul(tip(c[1]) , c[2] ,Nvar),Nvar)
            rem = im_div(cm , m2 , Nvar)
            if rem!=None and rem[0]!=0 and rem[1]!=0:RED.append( (c,cm,s) )
        for c in RED:
            B.remove( c )
        for i,(tf,p) in enumerate(zip(masks,G)):
            if not tf:continue
            for c in dp_obs(p,p2):
                cm = im_mul(c[0] , im_mul(tip(c[1]) , c[2] , Nvar) , Nvar)
                if NoSugar:
                   s_ps = 0
                else:
                   s_ps = max(sugars[i]+ideg(c[0],Nvar)+ideg(c[2],Nvar),
                              s_h+ideg(c[3],Nvar)+ideg(c[5],Nvar))
                B.append( (c,cm,s_ps) )
        for n,p in enumerate(G[:-1]):
            rem = im_div(tip(p) , tip(p2) , Nvar)
            if rem!=None:masks[n] = False
        B.sort(key=lambda x:(x[2],x[1]) , reverse=True)
        print ("{0} obstructions removed (new basis deg={1}) nb={2} nab={3} rp={4} sugar={5}".format(Nobs,ideg(tip(p2),Nvar),sum(masks),len(G),len(B),s_h))
    G = [p for (tf,p) in zip(masks,G) if tf]
    print("start reducing (total obstructions={0} , number of bases={1})\n".format(Nobs,len(G)))
    RG = []
    for n,p in enumerate(G):
        p = dp_nf(p,RG+G[n+1:])
        if p!=0:RG.append( p*Fraction(1,p.coeffs[tip(p)]) )
    assert(sum(masks)==len(RG)),"No of minimal bases!=No of reduced bases?"
    return RG




def mora(G , varlist):
    G1 = [p2dp(p,varlist) for p in G]
    G2 = dp_mora(G1)
    return [dp2p(p,varlist) for p in G2]


def p_nf(p , G , vars):
    r = dp_nf(p2dp(p , vars) , [p2dp(q,vars) for q in G])
    return dp2p(r,vars)



"""
S : Groebner bases for two-sided ideal
T : "S-compatible" Groebner bases for left ideal
"""
def gs_nf(p , S , T , vars):
    def tip(p):
        im = max(p2dp(p , vars).coeffs.keys())
        return i2m(im,vars)
    if p==0:return p
    p0 = p_nf(p , S , vars)
    if p0==0:return p0
    G = [p for p in T if p!=0]
    for p in G:p.simplify()
    r = Expression()
    HMG = [tip(p) for p in G]
    while True:
       if p0==0:break
       m = tip(p0)
       for (p1,m1) in zip(G,HMG):
           if m.deg()<m1.deg():continue
           if m.varlist[-m1.deg():]==m1.varlist:
               u = MonomialExpression()
               u.varlist = m.varlist[:-m1.deg()]
               c0,c1 = p0.coeffs[m],p1.coeffs[m1]
               p0 = (c1*p0 - c0*Expression(u)*p1)*Fraction(1,c1)
               p0 = p_nf(p0 , S , vars)
               break
       else:
           r += Expression(m)*p0.coeffs[m]
           del p0.coeffs[m]
    return r


def right_justified_completion(S , T , vars):
    def tip(p):
        im = max(p2dp(p , vars).coeffs.keys())
        return i2m(im,vars)
    G = []
    for p in T:
        if p==0:continue
        p.simplify()
        m = tip(p)
        G.append( (p*Fraction(1,p.coeffs[m]),m) )
    while True:
        G1 = []
        G0 = [x[0] for x in G]
        for (p0,m0) in G:
            for (p1,m1) in G:
                if m0.deg()<m1.deg():continue
                if m0.varlist[-m1.deg():]==m1.varlist:
                    u = MonomialExpression()
                    u.varlist = m0.varlist[:-m1.deg()]
                    h = p0 - Expression(u)*p1*Fraction(p0.coeffs[tip(p0)] , p1.coeffs[tip(p1)])
                    p2 = gs_nf(h , S , G0 , vars)
                    if p2!=0:
                        m2 = tip(p2)
                        p2 = p2*Fraction(1,p2.coeffs[m2])
                        G1.append( (p2,m2) )
        if len(G1)==0:
            RG = []
            G0 = [p for (p,_) in G]
            for n,p in enumerate(G0):
                p = gs_nf(p, S , RG+G0[n+1:],vars)
                if p!=0:RG.append( p*Fraction(1,p.coeffs[tip(p)]) )
            return RG
        else:
           G.extend(G1)


"""
Groebner-Shirshov pairs for (_S,_T)

"""
def gs(_S , _T , vars):
    def tip(p):
        im = max(p2dp(p , vars).coeffs.keys())
        return i2m(im,vars)
    S = mora(_S , vars)
    T = right_justified_completion(S , _T , vars)
    while True:
        B = []
        T1 = []
        for p in S:
            for q in T:
                B.extend( dp_obs(p2dp(p,vars) , p2dp(q,vars)) )
        for (ilm0,dp0,irm0,ilm1,dq1,irm1) in  B:
            """
               obstruction pairs for $(p,q) \in (S,T)$ should be the following forms:
               (1) p - bq
               (2) pa - bq
               (3) apb - q
            """
            if irm1!=0:continue
            lm0,rm0 = i2m(ilm0,vars),i2m(irm0,vars)
            lm1,rm1 = i2m(ilm1,vars),i2m(irm1,vars)
            h_l = Expression(lm0)*dp2p(dp0,vars)*Expression(rm0)
            h_r = Expression(lm1)*dp2p(dq1,vars)*Expression(rm1)
            assert(tip(h_l)==tip(h_r))
            c0,c1 = h_l.coeffs[tip(h_l)] , h_r.coeffs[tip(h_r)]
            h = gs_nf(h_l-h_r*Fraction(c0,c1) , S , T , vars)
            if h!=0:T1.append( h*Fraction(1 , h.coeffs[tip(h)]) )
        if len(T1)==0:
            break
        else:
            T2 = []
            for n,p in enumerate(T1):
                h = gs_nf(p , S , T+T2+T1[n+1:] , vars)
                if h!=0:
                    T2.append( h*Fraction(1 , h.coeffs[tip(h)]) )
            if len(T2)==0:break
            T.extend( T2 )
            T = right_justified_completion(S , T , vars)
    return (S,T)



def eqset(X , Y):
   if len(X)!=len(Y):return False
   for x in X:
      if not any([x==y for y in Y]):return False
   for y in Y:
      if not any([x==y for x in X]):return False
   return True



def homogenize(G , vars , t):
    assert(isinstance(t,Symbol))
    tname = list(t.coeffs.keys())[0].varlist[0]
    ret = []
    for p in G:
        d = max([len(m.varlist) for m in p.coeffs.keys()])
        q = Expression()
        for (m,c) in p.coeffs.items():
            m2 = MonomialExpression()
            m2.varlist = m.varlist + ([tname]*(d-len(m.varlist)))
            q.coeffs[m2] = c
        ret.append( q )
    for v in vars:
        ret.append( t*v-v*t )
    return ret



def test_mora():
   x,y,z,w = Symbol("x"),Symbol("y"),Symbol("z"),Symbol("w")
   assert( eqset( mora([x*y-x,y*x-y] , [x,y]) , [x*y-x,y*x-y,y*y-y,x*x-x]) )
   assert( eqset( mora([x*y*x-x,y*x-x,x*x-y] , [x,y]) , [x*x-x,-x+y]) )
   a,b,c = Symbol("a"),Symbol("b"),Symbol("c")
   assert( eqset( mora([a*b*c+c,b*a*b+a],[c,b,a]) , [a*b*c+c,b*a*b+a,a*c-b*c,-b*a*a+a*a*b]) )
   #-- mora's paper example2.7
   assert( eqset( mora([a*a*a-a , a*b*b*b-a , b*b*b-a] , [b,a]) , [b**3-a,a**2-a,a*b-b*a]) )
   #--
   assert( eqset( mora([a*b*a-b , b*a*b-b], [a,b]) , [b*b-a*b,b*a-a*b,a*a*b-b]) )
   #----
   assert( eqset( mora([a*a*b-1 , a*b*b-1],[a,b]) , [-a+b, a*a*a-1]) )
   #--
   assert( eqset( mora([x*x-y*y],[y,x]) , [x*x-y*y,-y*y*x+x*y*y]) )
   #--
   GB = [x*x-(x*y+7*y*y)*Fraction(1,3) , -y*x*y+x*y*x , 
         -Fraction(3,22)*y*x*y+Fraction(7,22)*(y**3)-Fraction(21,22)*y*y*x+x*y*y,
         -Fraction(8526,7601)*(y**4)+ y*y*x*y , Fraction(19307,22803)*(y**4)+(y**3)*x , y**5]
   RG = mora([-x*y-7*y*y+3*x*x ,x*y*x-y*x*y],[y,x])
   assert( eqset( RG , GB) ),RG
   #--
   GB = [a*b-c , b*c-a , c*a-b , -a*a+c*c , -a*a+b*b , a*a*a-c*b , b*a*a-a*c , a*a*c-b*a ,
         -a*c*b + c*b*a , -a*c*b + b*a*c]
   assert( eqset(mora([a*b-c,b*c-a,c*a-b],[a,b,c]) , GB) )
   #-- example from B.J. Keller's "Algorithms and Orders for Finding Noncommutative Groebner Bases"
   GB = [x*y*z*y-w , y*z*y*w , x*y*w*y - z , y*w*y*z , z*z , w*w , w*z*y*w , z*w*y*z]
   assert( eqset(mora([x*y*z*y-w , y*z*y*w , x*y*w*y-z,y*w*y*z],[x,y,z,w]) , GB) )


def test_gs():
    #-- sl_3 and its irrep L(1,0)
    #-- Note: The basis of this irrep can be {1,f1,f2*f1}
    f1,f2 = Symbol("f1"),Symbol("f2")
    S0 = [f2*f2*f1-2*f2*f1*f2+f1*f2*f2 , f2*f1*f1-2*f1*f2*f1+f1*f1*f2]
    T0 = [f1*f1,f2]   #-- vanishing of singular vectors
    S,T = gs(S0,T0,[f1,f2])
    assert( eqset(S , S0) )
    assert( eqset(T , [f1**2 , f2 , f1*f2*f1]) )
    #-- sl_3 and its irrep L(1,1)
    f1,f2 = Symbol("f1"),Symbol("f2")
    S0 = [f2*f2*f1-2*f2*f1*f2+f1*f2*f2 , f2*f1*f1-2*f1*f2*f1+f1*f1*f2]
    T0 = [f1*f1,f2*f2]   #-- vanishing of singular vectors
    S,T = gs(S0,T0,[f1,f2])
    assert( eqset(S , S0) )
    assert( eqset(T , [f1**2, f2**2, (f1**3)*f2, -Fraction(1,2)*(f1**2)*f2 + f1*f2*f1, f2*f1*f2*f1*f2, f1**2*f2*f1*f2]) )



#-- see "A presentation for the Virasoro and super-Virasoro algebras"
def virasoro():
   L3 = Symbol("L_3")
   K2 = Symbol("L_{-2}")
   L1 = (L3*K2-K2*L3)*Fraction(1,5)
   K1 = (L1*K2-K2*L1)*Fraction(1,3)
   L2 = (L3*K1-K1*L3)*Fraction(1,4)
   L0 = (L1*K1-K1*L1)*Fraction(1,2)
   L5 = (L3*L2-L2*L3)
   c = ((L2*K2-K2*L2)-4*L0)*Fraction(1,6)   #-- Virasoro central charge
   c1 = (L3*L0 -L0*L3)-3*L3
   c2 = (L0*K2-K2*L0)-2*K2
   c4 = (L2*L1-L1*L2)-L3
   c6 = (L5*L2-L2*L5) - L5
   I = [c1,c2,c4,c6]
   G = mora(I , [L3,K2])
   return(G)


def bershadsky():
    H,E,F,C=Symbol("H"),Symbol("E"),Symbol("F"),Symbol("C")
    c1 = H*E-E*H-E
    c2 = H*F-F*H+F
    c3 = E*F-F*E-H*H-C
    c4 = C*E-E*C
    c5 = C*F-F*C
    c6 = C*H-H*C
    I = [c1,c2,c3,c4,c5,c6]
    G = mora(I , [C,E,F,H])
    return(G)


def braid3_11():
   x,y,z = Symbol("x"),Symbol("y"),Symbol("z")
   generators = [y*x*y-z*y*z , x*y*x-z*x*y , z*x*z-y*z*x , x**4+y**3+z**3+x*y*z]
   return mora(generators , [x,y,z])



def braid62():
   x,y,z = Symbol("x"),Symbol("y"),Symbol("z")
   generators = [y*x*y-z*y*z , x*y*x-z*x*y , z*x*z-y*z*x , x**3-2*y**3+3*z**3-4*x*y*z+5*x*z**2-6*x*y**2+7*x*x*z-8*x*x*y]
   return mora(generators , [x,y,z])


def lp1_10():
   x,y,z = Symbol("x"),Symbol("y"),Symbol("z")
   generators = [z**4+y*x*y*x-x*y*y*x-3*z*y*x*z,x**3+y*x*y-x*y*x,z*y*x-x*y*z+z*x*z]
   return mora(generators , [z,y,x])


def lv2_15():
   x,y,z = Symbol("x"),Symbol("y"),Symbol("z")
   I = [x*y+y*z , x*x+x*x*y-y*x-y*y]
   return mora(I , [x,y,z])



if __name__=="__main__":
   test_mora()
   test_gs()


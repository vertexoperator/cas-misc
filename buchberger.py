# -*- coding:utf-8 -*-
try:
  import cPickle as pickle
except:
  import pickle

from numbers import *
from fractions import Fraction
import copy


class Monomial(object):
   def __init__(self):
      self.degs = {}
   def __hash__(self):
      return hash(pickle.dumps(self.degs))
   def __eq__(self , other):
      if len(self.degs)!=len(self.degs):
         return False
      else:
         for k,v in self.degs.items():
             if other.degs.get(k,0)!=v:
                return False
         return True
   def __mul__(self , other):
      ret = Monomial()
      vars = set(self.degs.keys()) | set(other.degs.keys())
      for var in vars:
         d = self.degs.get(var , 0) + other.degs.get(var , 0)
         if d!=0:
            ret.degs[var] = d
      return ret
   def __str__(self):
      if len(self.degs)==0:
         return str(1)
      else:
         terms = []
         for k,v in self.degs.items():
             if v==0:
                pass
             elif v==1:
                terms.append( str(k) )
             else:
                terms.append( "%s^%d" % (k,v) )
         return "*".join( terms )



class Polynomial(object):
   def __init__(self):
       self.coeffs = {}
   def __eq__(self , other):
       if not isinstance(other , Polynomial):
           if len(self.coeffs)>1:
                return False
           elif len(self.coeffs)==1:
                k,v = self.coeffs.items()[0]
                return (len(k.degs)==0 and v==other)
           else:
                return (other==0)
       else:
           for k,v in self.coeffs.items():
               if other.coeffs.get(k , 0)!=v:
                  return False
           for k,v in other.coeffs.items():
               if self.coeffs.get(k , 0)!=v:
                  return False
           return True
   def __ne__(self,other):
       return not(self==other)
   def __add__(self , other):
       if isinstance(other , Polynomial):
          ms = set(self.coeffs.keys()) | set(other.coeffs.keys())
          ret = Polynomial()
          for m in ms:
              c = self.coeffs.get(m , 0) + other.coeffs.get(m , 0)
              if c!=0:ret.coeffs[m] = c
          return ret
       else:
          ret = Polynomial()
          cm = Monomial()
          for m,c in self.coeffs.items():
              ret.coeffs[m] = c
          if other!=0:ret.coeffs[cm] = ret.coeffs.get(cm , 0) + other
          return ret
   def __radd__(self , other):
       return self.__add__(other)
   def __sub__(self , other):
       if isinstance(other , Polynomial):
          ms = set(self.coeffs.keys()) | set(other.coeffs.keys())
          ret = Polynomial()
          for m in ms:
              c = self.coeffs.get(m , 0) - other.coeffs.get(m , 0)
              if c!=0:ret.coeffs[m] = c
          return ret
       else:
          ret = Polynomial()
          cm = Monomial()
          for m,c in self.coeffs.items():
              ret.coeffs[m] = c
          if other!=0:ret.coeffs[cm] = ret.coeffs.get(cm , 0) - other
          return ret
   def __neg__(self):
       ret = Polynomial()
       for k,v in self.coeffs.items():
           ret.coeffs[k] = -v
       return ret
   def __rsub__(self , other):
       tmp = self.__sub__(other)
       return tmp.__neg__()
   def __mul__(self , other):
       if isinstance(other , Polynomial):
          ret = Polynomial()
          for m1,c1 in self.coeffs.items():
             for m2,c2 in other.coeffs.items():
                 m = m1*m2
                 c = c1*c2
                 if c!=0:ret.coeffs[m] = ret.coeffs.get(m , 0) + c
          return ret
       else:
          ret = Polynomial()
          for m,c in self.coeffs.items():
             c0 = c*other
             if c0!=0:ret.coeffs[m] = ret.coeffs.get(m , 0) + c0
          return ret
   def __rmul__(self , other):
       return self.__mul__(other)
   def __pow__(self , n):
       assert(isinstance(n,int))
       assert(n>=0)
       ret = 1
       for _ in range(n):
           ret *= self
       return ret
   """
   def __call__(self,**args):   #-- 代入操作
       vars = [str(v) for v in args.keys()]
       ucvars = [("v_%d" % n) for n in xrange(len(vars))]
       terms = []
       for (m,c) in self.coeffs.items():
           if len(m.degs)==0:
              terms.append( str(c) )
           else:
              mterms = []
              for (varname , degree) in m.degs.items():
                  if varname in vars:
                      mterms.append('%s**%d' % (ucvars[vars.index(varname)] , degree))
                  else:
                      mterms.append('Symbol("%s")**%d' % (str(varname) , degree))
              terms.append('%s*%s' % (str(c) , '*'.join(mterms)))
       if len(terms)==0:
           pyexpr = ("lambda %s:%s" % (",".join(ucvars) , "0"))
       else:
           pyexpr = ("lambda %s:%s" % (",".join(ucvars) , "+".join(terms)))
       return apply(eval(pyexpr) , args.values())
   """
   def __str__(self):
       terms = []
       for idx,(m,c) in enumerate(self.coeffs.items()):
           if c==0:
               pass
           elif idx==0 and len(m.degs.keys())==0:
               terms.append( str(c) )
           elif idx==0 and c==1:
               terms.append( str(m) )
           elif idx==0:
               terms.append( "%s*%s" % (str(c),str(m)) )
           elif len(m.degs.keys())==0:
               if c<0:
                  terms.append( "- " + str(-c) )
               else:
                  terms.append( "+ " + str(c) )
           elif c==1:
               terms.append( "+ " + str(m) )
           elif c==-1:
               terms.append( "- %s" % str(m) )
           elif c<0:
               terms.append( "- %s*%s" % (str(-c),str(m)) )
           else:
               terms.append( "+ %s*%s" % (str(c),str(m)) )
       if len(terms)!=0:
          return " ".join(terms)
       else:
          return "0"
   def __repr__(self):
       return str(self)


class Variable(Polynomial):
   def __init__(self , s):
      super(Variable , self).__init__()
      m = Monomial()
      m.degs[s] = 1
      self.coeffs[m] = 1
   @property
   def name(self):
      for m in self.coeffs.keys():
          for (s,d) in m.degs.items():
              if d==1:return s



def iadd(v1,v2):
    return tuple([(c1+c2) for (c1,c2) in zip(v1,v2)])


def isub(v1,v2):
    return tuple([(c1-c2) for (c1,c2) in zip(v1,v2)])


def idiv(v1,v2):
    ret = []
    for (c1,c2) in zip(v1,v2):
        c = c1-c2
        if c<0:return None
        ret.append(c)
    return tuple(ret)


class DPolynomial:
   __slots__ = ["Nvar","Nweight","coeffs","weights","weightMat","_normalized"]
   def __init__(self , Nvar , Nweight , weightMat):
       self.Nvar = Nvar
       self.Nweight = Nweight
       self.coeffs = {}
       self.weights = {}
       self.weightMat = weightMat
       self._normalized = False
   def normalize(self):
       if self._normalized:return self
       zero_keys = []
       for m,c in self.coeffs.items():
           if c==0:
               zero_keys.append( m )
           elif not m in self.weights:
               w = [0]*(self.Nweight)
               for i in range(self.Nweight):
                   for j in range(self.Nvar):
                       w[i] = w[i] + self.weightMat[i][j]*m[j]
               self.weights[m] = tuple(w)
       for m in self.weights.keys():
           if not m in self.coeffs:zero_keys.append( m )
       for m in zero_keys:
           if m in self.coeffs:del self.coeffs[m]
           if m in self.weights:del self.weights[m]
       self._normalized = True
       return self
   @property
   def tip(self):    #-- leading monomial
       self.normalize()
       if len(self.coeffs)==0:
          return tuple([0]*self.Nweight)
       else:
          m,w = max(self.weights.items() ,key=lambda x:x[1])
          return m
   def __add__(lhs , rhs):
       ret = copy.deepcopy(lhs)
       ret._normalized = False
       if isinstance(rhs,Number):
           m = tuple([0]*lhs.Nvar)
           ret.coeffs[m] = ret.coeffs.get(m,0)+rhs
           ret.weights[m] = tuple([0]*lhs.Nweight)
       elif isinstance(rhs,DPolynomial):
           for (m,c) in rhs.coeffs.items():
               ret.coeffs[m] = ret.coeffs.get(m,0)+c
               ret.weights[m] = rhs.weights[m]
       return ret.normalize()
   def __neg__(self):
       ret = copy.deepcopy(self)
       for (m,c) in ret.coeffs.items():
           ret.coeffs[m] = -c
       return ret.normalize()
   def __mul__(lhs,rhs):
       if isinstance(rhs,Number):
           if rhs==0:return DPolynomial(lhs.Nvar , lhs.Nweight , lhs.weightMat)
           ret = copy.deepcopy(lhs)
           ret._normalized = False
           for m,c in ret.coeffs.items():
               ret.coeffs[m] = c*rhs
       else:
           ret = DPolynomial(lhs.Nvar , lhs.Nweight , lhs.weightMat)
           for (m1,c1) in lhs.coeffs.items():
              if c1==0:continue
              for (m2,c2) in rhs.coeffs.items():
                if c2==0:continue
                m3 = iadd(m1,m2)
                ret.coeffs[m3] = ret.coeffs.get(m3,0) + (c1*c2)
                ret.weights[m3] = iadd(lhs.weights[m1] , rhs.weights[m2])
       return ret.normalize()
   def __sub__(lhs,rhs):
       ret = copy.deepcopy(lhs)
       ret._normalized = False
       if isinstance(rhs,Number):
           m = tuple([0]*lhs.Nvar)
           ret.coeffs[m] = ret.coeffs.get(m,0)-rhs
           ret.weights[m] = tuple([0]*lhs.Nweight)
       elif isinstance(rhs,DPolynomial):
           for (m,c) in rhs.coeffs.items():
               ret.coeffs[m] = ret.coeffs.get(m,0)-c
               ret.weights[m] = rhs.weights[m]
       return ret.normalize()
   def __eq__(lhs , rhs):
       lhs.normalize()
       if isinstance(rhs,DPolynomial):
           rhs.normalize()
           if len(lhs.coeffs)!=len(rhs.coeffs):return False
           if lhs.Nvar!=rhs.Nvar:return False
           for m,c in lhs.coeffs.items():
               if c!=rhs.coeffs.get(m,0):return False
           return True
       elif isinstance(rhs,Number):
           if len(lhs.coeffs)>1:
               return False
           if len(lhs.coeffs)==0:return (rhs==0)
           m,c = list(lhs.coeffs.items())[0]
           if c==rhs and all([x==0 for x in m]):return True
           return False
   def __ne__(lhs,rhs):
       return not (lhs==rhs)
   def __rmul__(lhs,rhs):
       return rhs*lhs
   def __radd__(lhs,rhs):
       return rhs+lhs
   def __rsub__(lhs,rhs):
       return (-rhs)+lhs


def dp2p(p , vars):
    assert(all([isinstance(s,Variable) for s in vars]))
    assert(p.Nvar==len(vars))
    ret = Polynomial()
    for deglist,c in p.coeffs.items():
        m = Monomial()
        for idx,deg in enumerate(deglist):
            if deg!=0:
                m.degs[vars[idx].name] = deg
        ret.coeffs[m] = c
    return ret


def p2dp(p , vars , weightMat):
    assert(all([isinstance(s,Variable) for s in vars]))
    Nweight = len(weightMat)
    ret = DPolynomial(len(vars) , Nweight , copy.deepcopy(weightMat))
    for m,c in p.coeffs.items():
        assert(all([(v in [s.name for s in vars]) for v in m.degs.keys()])),("unknown variable:{0} {1}".format(v , [str(c) for c in vars]))
        deglist = tuple([m.degs.get(v.name,0) for v in vars])
        ret.coeffs[deglist] = c
    ret.normalize()
    return ret


#-- graded reverse lexcographic order
def grevlex(Nvar):
    weightMat = [[0]*Nvar for _ in range(Nvar+1)]
    for i in range(Nvar+1):
        for j in range(Nvar):
           if i==0:weightMat[i][j] = 1
           elif i+j==Nvar:weightMat[i][j] = -1
    return weightMat


#-- pure lexicographic order
def purelex(Nvar):
    weightMat = [[0]*Nvar for _ in range(Nvar)]
    for i in range(Nvar):
        for j in range(Nvar):
           if i==j:weightMat[i][j] = 1
    return weightMat


#-- total degree lexicographic order
def deglex(Nvar):
    weightMat = [[0]*Nvar for _ in range(Nvar+1)]
    for i in range(Nvar+1):
        for j in range(Nvar):
           if i==0:weightMat[i][j] = 1
           elif i==j:weightMat[i][j] = 1
    return weightMat


def dp_nf(p , G):
    Nweight = p.Nweight
    h = copy.deepcopy(p)
    r = DPolynomial(h.Nvar , h.Nweight , h.weightMat)
    HMG = [(p , p.tip) for p in G if p!=0]
    if h==0:return 0
    while True:
        m = h.tip
        if h==0:break
        c_h = h.coeffs[m]
        if sum(m)==0:
            r.coeffs[m] = c_h
            break
        if c_h==0:
            del h.coeffs[m]
            continue
        w_h = h.weights[m]
        for (p1,m1) in HMG:
            rem = idiv(m,m1)
            if rem!=None:
                c1 = p1.coeffs[m1]
                for (m2,c2) in p1.coeffs.items():
                    m3 = iadd(rem , m2)
                    if c2!=0 and not m3 in h.weights:
                        h.weights[m3] = isub(iadd(p1.weights[m2] , w_h) , p1.weights[m1])
                    h.coeffs[m3] = h.coeffs.get(m3,0) - Fraction(c2*c_h,c1)
                    if h.coeffs[m3]==0:
                        del h.coeffs[m3]
                        del h.weights[m3]
                break
        else:
            del h.coeffs[m]
            r.weights[m] = w_h
            del h.weights[m]
            r.coeffs[m] = c_h
    return r


def p_nf(p , G , vars , order):
    G = [p2dp(q,vars,order) for q in G]
    r = p2dp(p,vars,order)
    r = dp_nf(r,G)
    return dp2p(r,vars)


import time
def dp_buchberger(_G):
    NoSugar = False
    def tdeg(p):
        return max( [sum(m) for m in p.coeffs.keys()] )
    nf_time = 0.0
    Nobs,ZR,NZR = 0,0,0
    G = [q*Fraction(1,q.coeffs[q.tip]) for q in _G if q!=0]
    #-- inter-reduce
    while True:
       ok = False
       for n,p in enumerate(G):
          p2 = dp_nf(p , [q for q in G if q!=p and q!=0])
          if p2==0:
              G[n] *= 0
          elif p2!=p:
              G[n] = p2*Fraction(1 , p2.coeffs[p2.tip])
              break
          elif n==len(G)-1:
              ok = True
              break
       if ok:break
    G = [p for p in G if p!=0]
    sugars = [tdeg(p) for p in G]
    masks = [True]*len(G)
    if len(G)==0:return []
    #-- find first obstructions
    B = []
    for i,p in enumerate(G):
       for j,q in enumerate(G):
           if i<j:
              mp , mq = p.tip , q.tip
              #-- Buchberger's product criterion to avoid unnecessary reduction
              if all([min(i1,i2)==0 for (i1,i2) in zip(mp,mq)]):continue
              mpq = tuple([max(i1,i2) for (i1,i2) in zip(mp,mq)])
              if NoSugar:
                  s_pq = 0
              else:
                  s_pq = max(sugars[i]+sum(mpq)-sum(mp) , sugars[j]+sum(mpq)-sum(mq))
              B.append( (i , j , idiv(mpq,mp) , idiv(mpq,mq) , mpq , s_pq) )
    B.sort(key=lambda x:(x[5],x[4]) , reverse=True)
    while len(B)>0:
        i,j,um,vm,_,s_h = B.pop()
        Nobs+=1
        p,q = G[i],G[j]
        tp = DPolynomial(p.Nvar , p.Nweight , p.weightMat)
        tq = DPolynomial(q.Nvar , q.Nweight , q.weightMat)
        tp.coeffs[um] = 1
        tq.coeffs[vm] = 1
        tp._normalized = False
        tq._normalized = False
        tp.normalize()
        tq.normalize()
        t0 = time.time()
        h = dp_nf(tp*p - tq*q , [t for (tf,t) in zip(masks,G) if tf])
        t1 = time.time()
        nf_time += (t1-t0)
        if h==0:
            ZR+=1
            continue
        NZR+=1
        h = h*Fraction(1 , h.coeffs[h.tip])
        #-- useless pair elimination of  Gebauer-Moeller
        RED = []
        htip = h.tip
        for obs in B:
            rem = idiv(obs[4] , htip)
            if rem==None:continue
            i,j,_,_,mfg,_ = obs
            ftip,gtip = G[i].tip,G[j].tip
            mfh = tuple([max(i1,i2) for (i1,i2) in zip(ftip,htip)])
            mgh = tuple([max(i1,i2) for (i1,i2) in zip(gtip,htip)])
            if mfh!=mfg and mgh!=mfg:
                 RED.append( obs )
        for obs in RED:
            B.remove(obs)
        #-- new obstructions
        for i,p in enumerate(G):
            if not masks[i]:continue
            mp , ms = p.tip , h.tip
            mps = tuple([max(i1,i2) for (i1,i2) in zip(mp,ms)])
            if iadd(mp,ms)==mps:continue
            if any([obs[0]<i and obs[1]==len(G) and obs[4]==mps for obs in B]):
                continue
            if any([tf and idiv(mps , g.tip)!=None for (tf,g) in zip(masks,G)]):
                if NoSugar:
                   s_ps = 0
                else:
                   s_ps = max(sugars[i]+sum(mps)-sum(mp) , s_h+sum(mps)-sum(ms))
                B.append( (i , len(G) , idiv(mps,mp) , idiv(mps,ms) , mps , s_ps) )
        for n,p in enumerate(G):
            rem = idiv(p.tip , h.tip)
            if rem!=None:masks[n] = False
        G.append( h )
        masks.append( True )
        sugars.append( s_h )
        B.sort(key=lambda x:(x[5],x[4]) , reverse=True)
        print("{0} obstruction killed nb={1} nab={2} rp={3}".format(Nobs ,sum(masks), len(G) ,len(B)))
    #-- find reduced basis
    G = [p for (tf,p) in zip(masks,G) if tf] 
    RG = []
    for n,p in enumerate(G):
        p = dp_nf(p,RG+G[n+1:])
        if p!=0:RG.append( p*Fraction(1,p.coeffs[p.tip]) )
    assert(len(G)==len(RG)),"G should be minimal bases"
    print("total obstructions={0} , NF time={1:.3f} ZR={2} NZR={3}\n".format(Nobs,nf_time,ZR,NZR))
    return RG




def groebner(_G , vars , order):
    G = [p2dp(q,vars,order) for q in _G]
    G = dp_buchberger(G)
    return [dp2p(q,vars) for q in G]


def eqset(X , Y):
   if len(X)!=len(Y):return False
   for x in X:
      if not any([x==y for y in Y]):return False
   for y in Y:
      if not any([x==y for x in X]):return False
   return True


if __name__=="__main__":
    x,y,z = Variable("x"),Variable("y"),Variable("z")
    GB = groebner([x*y+2 , x*x*x+x],[x,y],grevlex(2))
    assert( eqset(GB , [x - Fraction(1,2)*y , y*y+4]) ),"test-1 failed"
    assert( eqset(groebner([-x**3+y, x**2*y-z], [x,y,z], grevlex(3)),[x**3 - y, y*(x**2) - z, y**2 - x*z]) ),"test-2 failed"
    assert( eqset(groebner([x*x*x-2*x*y , x*x*y-2*y*y+x],[x,y],grevlex(2)),[x**2, y**2 - Fraction(1,2)*x, y*x]) ),"test-3 failed"
    assert( eqset(groebner([x*x*x-2*x*y , x*x*y-2*y*y+x],[x,y],purelex(2)),[-2*y**2 + x, y**3]) ),"test-4 failed"
    #-- test-5
    x,y,z,w = Variable("x"),Variable("y"),Variable("z"),Variable("w")
    GB5 = [y**2*x*z - w, y**2*x*w - z, z**2, w**2, z*w]
    assert( eqset(groebner([x*y*y*z-w , y*y*z*w,x*y*y*w-z],[x,y,z,w],grevlex(4)),GB5) )
    #-- katsura-2 benchmark
    u0,u1,u2 = Variable("u0"),Variable("u1"),Variable("u2")
    I = [u0+2*u2+2*u1-1,2*u1*u0+2*u1*u2-u1,u0**2-u0+2*u2**2+2*u1**2]
    GB = [2*u2 - 1 + 2*u1 + u0, Fraction(1,5)*u2 - Fraction(3,5)*u2**2 - Fraction(1,5)*u1 + u1**2, 
          Fraction(-2,5)*u2 + Fraction(6,5)*u2**2 - Fraction(1,10)*u1 + u1*u2, 
          Fraction(1,70)*u2 - Fraction(79,210)*u2**2 + Fraction(1,30)*u1 + u2**3]
    assert( eqset(groebner(I , [u0,u1,u2] , grevlex(3)) , GB) ),"katsura-2 failed"
    #--katsura-3
    u0,u1,u2,u3 = Variable("u0"),Variable("u1"),Variable("u2"),Variable("u3")
    I = [u0+2*u3+2*u2+2*u1-1,2*u2*u0+2*u1*u3-u2+u1**2,2*u1*u0+2*u2*u3+2*u1*u2-u1,u0**2-u0+2*u3**2+2*u2**2+2*u1**2]
    GB = groebner(I,[u0,u1,u2,u3],grevlex(4))
    assert(len(GB)==7),"katsura-3 failed"
    #-- homogenized katsura-4
    h,u0,u1,u2,u3,u4 = Variable("h"),Variable("u0"),Variable("u1"),Variable("u2"),Variable("u3"),Variable("u4")
    I = [u0+2*u3+2*u2+2*u1+2*u4-h,2*u3*u0-u3*h+2*u1*u2+2*u4*u1,
         2*u2*u0+2*u1*u3+(2*u4-h)*u2+u1**2,
         2*u1*u0+(2*u2+2*u4)*u3+2*u1*u2-u1*h,
         u0**2-u0*h+2*u3**2+2*u2**2+2*u1**2+2*u4**2]
    GB = groebner(I,[h,u0,u1,u2,u3,u4],grevlex(6))
    assert(len(GB)==13),"h-katsura-4 failed"


OVERVIEW
--------
commutative and noncommutative Groebner basis


Environment
-----------

- CPython 2.7.5
- pypy3  2.4.0 (Python3.2.5 compatible)


Benchmarks for Buchberger algorithm
-----------------------------------
- The base field is Q
- The timings measured on Core i5-2430M (2.40GHz) with 8GB RAM
- Examples are taken from http://invo.jinr.ru/ginv/benchmark.html


|               | buchberger.py  | Risa/Asir(gr) | Risa/Asir(nd_gr)|
|:--------------|:--------------:|:-------------:|:---------------:|
|cyclic-6       | 5.647sec       | 0.39sec       |  0.624sec       |
|katsura-6      | 4.361sec       | 1.06sec       |  0.905sec       |
|katsura-7      | 32.401sec      | 7.441sec      |  10.72sec       |
|eco-8          | 14.074sec      | 2.262sec      |  1.654sec       |
|hcyclic-7      | 665.529sec     | 22.6sec       |  9.578sec       |


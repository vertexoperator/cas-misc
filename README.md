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


|               | buchberger.py  | Risa/Asir(gr) | Risa/Asir(nd_gr)|Singular      |
|:--------------|:--------------:|:-------------:|:---------------:|:------------:|
|cyclic-6       | 4.368sec       | 0.39sec       |  0.624sec       | 0.130sec     |
|katsura-6      | 2.829sec       | 1.06sec       |  0.905sec       |              |
|katsura-7      | 24.174sec      | 7.441sec      |  10.72sec       | 1.430sec     |
|katsura-8      |                | 70.11sec      |  66.47sec       | 11.830sec    |
|eco-8          | 10.268sec      | 2.262sec      |  1.654sec       |              |
|eco-9          | 144.912sec     | 36.44sec      |  13.74sec       | 1.780sec     |
|hcyclic-7      | 606.613sec     | 22.6sec       |  9.578sec       | 3.430sec     |


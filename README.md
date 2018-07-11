# BCJR Algorithm in GF(2^n)
## Usage
```
python3 trellis.py --plot_sections 0 5 --input=test_data/file1.csv --b=1
python3 trellis.py --plot_sections 0 7 --input=test_data/file2.csv --b=1
python3 trellis.py --plot_sections 0 5 --input=test_data/file3.csv --b=2
```
## Introduction
This Python code can draw a trellis from any linear block code in GF(2^n)
- Input: csv
- Output: trellis drawn on the screen
A (5,3) linear block code over GF(2)
![file1](/img/file1.png)

# BCH bound
## Usage
```
python3 bound.py --n=5 --q=2 --verbose
python3 bound.py --n=5 --q=4 --verbose
python3 bound.py --n=15 --q=2 --verbose
python3 bound.py --n=15 --q=4 --verbose
```
## Flowchart
Take n=5, q=4 for example
1. Find m = 2, so extension field is GF(16)
2. Find cyclotonics group and its generator
  - \[0\]: x+1
  - \[3,12\]: x^2 + (a^2+a+1)x + 1
  - \[6,9\]: x^2 + (a^2+a)x + 1
3. Find mapping between GF(16) to GF(4)
  - a^0  = 1 = beta^0 = 1 on GF(4)
  - a^5  = a^2+a = beta^1 = b on GF(4)
  - a^10 = a^2+a+1 = beta^2 = b+1 on GF(4)
4. Find irreducible basis generator on GF(4)
  - x+1 -> x+1
  - x^2 + (a^2+a+1)x + 1 -> x^2 + (b+1)x + 1
  - x^2 + (a^2+a)x + 1 -> x^2 + bx + 1
5. Calculate BCH bound of all generator and its minimum weight. There are total 2^3-2 generators, excluding 0 and x^5-1 for every combination of generator basis
  - x^2 + bx + 1
    - Roots: a^6, a^9
    - BCH bound = 2
    - extension BCH bound = 3
    - Tzeng's bound = 3
    - Minimum weight = 3
  - x^2 + (b+1)x + 1
    - Roots: a^3, a^12
    - BCH bound = 2
    - extension BCH bound = 3
    - Tzeng's bound = 3
    - Minimum weight = 3
  - x^4 + x^3 + x^2 + x + 1
    - Roots: a^6, a^9, a^3, a^12
    - BCH bound = 2
    - extension BCH bound = 5
    - Tzeng's bound = 5
    - Minimum weight = 5
  - x + 1
    - Roots: a^6, a^9, a^3, a^12
    - BCH bound = 2
    - extension BCH bound = 2
    - Tzeng's bound = 2
    - Minimum weight = 2
  - x^3 + (b+1)x^2 + (b+1)x + 1
    - Roots: 1, a^6, a^9
    - BCH bound = 2
    - extension BCH bound = 4
    - Tzeng's bound = 4
    - Minimum weight = 4
  - x^3 + bx^2 + bx + 1
    - Roots: 1, a^3, a^12
    - BCH bound = 2
    - extension BCH bound = 4
    - Tzeng's bound = 4
    - Minimum weight = 4
  
# (Narrow-sense) BCH Code decoding
## Usage
There are several kinds of usage in this BCH decoding code
1. Use default r(x) and d0
  - r(x) =  1, 1, 1, 1, 1, 0, 0, 1, 1
  - d0 = 5
  ```
  python3 bch.py --verbose
  ```
2. Give r(x) and d0 from command line
  ```
  python3 bch.py --rx=111110011 --d0=5 --verbose
  ```

3. Give c(x) and use default e(x)
  - e(x) = x^3 + x^2
  ```
  python3 bch.py --cx=111010001 --verbose
  ```

4. Give c(x) and e(x)
  ```
  python3 bch.py --cx=111010001 --ex=1100 --verbose
  ```

## Flowchart
1. Define r(x).
2. Find error number v, which is max v such that det(M_v)=0 and M_v is a vxv error spectrum matrix with value M_v\[i\]\[j\] = r(alpha^(i+j)).
3. Calculate syndrone polynomial s(x), where s_i = r(alpha^i)
4. Calculate locator polynomial from Berlekamp-Massey Algorithm
5. Find roots of locator polynomial, which will be error locator
6. Calculate evaluation polynomial w(x) = sigma(x)\*s(x) mod x^v
7. Calculate error value Yi from Forney's algorithm
8. Recover e(x) from error locator pair (X,Y), where e(x) = sum(Yi\*x^(log(Xi)))

## Example
1. Define r(x). r(x) = x^8 + x^7 + x^6 + x^5 + x^4 + x + 1. From the result, we know that c(x) = x^8 + x^7 + x^6 + x^4 + 1, e(x) = x^5 + 1
2. Find error number v. There are two non-zero term in e(x) => v=2
3. Calculate syndrone polynomial s(x), where s_i = r(alpha^i)
  - s0 = r(a^0) = a^3 + 1
  - s1 = r(a^1) = a   + 1
  - s0 = r(a^2) = a^2
4. Calculate locator polynomial from Berlekamp-Massey Algorithm
  - From Berlekamp-Massey algorithm, we found that locator polynomial is (a^3+a^2)x^2 + (a^2)x + 1
5. Find roots of locator polynomial, which will be error locator
  - Roots #0 = a^2 + a + 1 = a^10 => X = (a^10)^(-1) = a^5
  - Roots #1 = a^3     + 1 = a^14 => X = (a^14)^(-1) = a^1
6. Calculate evaluation polynomial w(x) = sigma(x)\*s(x) mod x^v
  - Evaluation polynomial w(x) = a^2
7. Calculate error value Yi from Forney's algorithm
  - Y_0 = w(a^10)/sigma'(a^10) = (a^2)/(a^2) = 1
  - Y_1 = w(a^14)/sigma'(a^14) = (a^2)/(a^2) = 1
8. Recover e(x) from error locator pair (X,Y), where e(x) = sum(Yi\*x^(log(Xi)))
  - E_0 = (1,a^10) => e0(x) = x^5
  - E_1 = (1,a^14) => e1(x) = x
  - e(x) = e0(x) + e1(x) = x^5 + x
9. c(x) = r(x) - e(x) = x^8 + x^7 + x^6 + x^4 + 1

# Convolution code decoder
## Usage
```
python3 conv_code.py --gen=conv_csv/g.csv --out_seq=conv_csv/output.csv 
```
```
python3 conv_code.py --gen=conv_csv/g.csv --out_seq=conv_csv/output.csv --q=4
```
```
python3 conv_code.py --gen=conv_csv/g2.csv --out_seq=conv_csv/output.csv --k=2
```

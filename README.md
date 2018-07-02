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
python3 bound.py --n=5 --q=2 --verbose=1
python3 bound.py --n=5 --q=4
python3 bound.py --n=15 --q=2
python3 bound.py --n=15 --q=4
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
  

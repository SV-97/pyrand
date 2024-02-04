# pyrand
Pure rust implementation of (parts of) python's random module with compatible PRNG behaviour: seeding with equivalent values will yield identical PRNG output. So running 
```python
import random
rng = random.Random("Pizza")
rng.choices(range(20), k=5)
```
will yield `[3, 3, 15, 2, 16]` and running the equivalent
```rust
use pyrand::{PyMt19937, PySeedable, RandomChoiceIterator};
let rng = &mut PyMt19937::py_seed("Pizza");
assert_eq!(
    (0..20).choose(rng).take(5).collect::<Vec<_>>(),
    vec![3, 3, 15, 2, 16]
);
```
will also yield `[3, 3, 15, 2, 16]`.

## This seems pretty nieche. Why does this have to be a thing?

Assume you've written some python code using (pseudo-)random numbers like for example
```python
import random
random.seed(input("Enter seed: "))
# run some computations using random numbers
...
result = random.random()

print(f"Your result is: {result}")
```
and your users somehow end up depending on that result in some way. Then your app has to consistently return the same exact values if the user provides the same input. This binds you very closely to python and its current `random` implementation.

I actually ran into the unfortunate case of having some Python code where this was the case. I wanted to rewrite that code in Rust and kind off couldn't without also pulling in Python to handle the PRNG part (or statically linking parts of CPython or whatever) - which would've killed the whole rewrite. So long story short: I reimplemented the central functionality of python's random number generation code in rust, making sure to retain constistent output between the two implementations.

If you actually need this kind of compatibility for some of the functions I haven't implemented yet feel free to implement them and open a pull request - or file an issue on the GitHub repo.

## `rand` support

There are (optional) implementations for the `rand` traits but they might reimplement some "natively" available generators on top of the basic internal generator - and do so differently to how they're implemented in python. This means that if you use the `rand` interface to the generator you have to verify that you actually get matching output for the particular case you have.

### Implementation details

The implementation is based off [`_randommodule.c`](https://github.com/python/cpython/blob/3.12/Modules/_randommodule.c) and [`random.py`](https://github.com/python/cpython/blob/3.12/Lib/random.py) and I tried to retain most function and variablenames etc. from there. The basic algorithm is a [32-bit Mersenne Twister (MT19937)](https://en.wikipedia.org/wiki/Mersenne_Twister).
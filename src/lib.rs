#![feature(array_chunks)]

//! Implements a 32-bit mersenne twister specifically aiming to be able to
//! reproduce the results of python's `random` module.
//! Seeding this generator with the equivalent to a given python value will
//! produce identical values.
//!
//! This implementation is based on Python 3.13 and should be downwards compatible
//! at least down to Python 3.4.
//!
//! Function- and variable names as well as large parts of the implementation are
//! strongly oriented on the CPython source. I also took some comments from there.
//!
//! Note that while a `rand` implementation is given you shouldn't expect most methods
//! from `rand::Rand` to return the same values as their python counterparts.
//!
//! # Examples
//!
//! Choosing random values from a collection like
//! ```
//! use pyrand::{PyMt19937, PySeedable, RandomChoiceIterator};
//! let rng = &mut PyMt19937::py_seed("Pizza");
//! assert_eq!(
//!     (0..20).choose(rng).take(5).collect::<Vec<_>>(),
//!     vec![3, 3, 15, 2, 16],
//! );
//! ```
//! will return the same result as the equivalent python code
//! ```python
//! import random
//! rng = random.Random("Pizza")
//! assert rng.choices(range(20), k=5) == [3, 3, 15, 2, 16]
//! ```

use num::{BigInt, BigUint, One, Signed, Zero};
use sha2::{Digest, Sha512};

const MATRIX_A: u32 = 0x9908b0df;
const UPPER_MASK: u32 = 0x80000000;
const LOWER_MASK: u32 = 0x7fffffff;

const M_DEFAULT: usize = 397;

pub struct MTState<const N: usize> {
    mt: [u32; N], // state array
    mti: usize,   // index into state array
}
const N_PYTHON: usize = 624;
pub type PyMt19937 = MTState<N_PYTHON>;

impl<const N: usize> MTState<N> {
    pub fn init_genrand(seed: u32) -> Self {
        let mut mt = [0; N];
        mt[0] = seed;
        for mti in 1..N {
            mt[mti] = 1812433253_u32
                .wrapping_mul(mt[mti - 1] ^ (mt[mti - 1] >> 30))
                .wrapping_add(mti as u32);
        }
        MTState { mt, mti: N + 1 }
    }

    pub fn init_by_array(key: &[u32]) -> Self {
        let mut state = Self::init_genrand(19650218);
        let mt = &mut state.mt;

        let mut i = 1;
        let mut j = 0;
        for _ in (1..=key.len().max(N)).rev() {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1664525)))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= N {
                mt[0] = mt[N - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }

        for _ in (1..N).rev() {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1566083941)))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= N {
                mt[0] = mt[N - 1];
                i = 1;
            }
        }
        mt[0] = 0x80000000; // MSB is 1; assuring non-zero initial array
        state
    }

    pub fn genrand_uint32<const M: usize>(&mut self) -> u32 {
        let mt = &mut self.mt;
        let mag01 = [0, MATRIX_A];

        let mut y;
        if self.mti >= N {
            // generate N words at one time
            for kk in 0..(N - M) {
                y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[(y & 0x1) as usize];
            }

            for kk in N - M..N - 1 {
                y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] =
                    mt[kk.wrapping_add(M.wrapping_sub(N))] ^ (y >> 1) ^ mag01[(y & 0x1) as usize];
            }
            y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[(y & 0x1) as usize];
            self.mti = 0;
        }

        y = mt[self.mti];
        self.mti += 1;

        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;

        y
    }

    pub fn genrand_uint32_default(&mut self) -> u32 {
        self.genrand_uint32::<M_DEFAULT>()
    }

    pub fn genrand_res53<const M: usize>(&mut self) -> f32 {
        let a = self.genrand_uint32::<M>() >> 5;
        let b = self.genrand_uint32::<M>() >> 6;
        (a as f32 * 67108864.0 + b as f32) * (1.0 / 9007199254740992.0)
    }

    pub fn genrand_res53_default(&mut self) -> f32 {
        self.genrand_res53::<M_DEFAULT>()
    }

    pub fn getrandbits(&mut self, num_bits: usize) -> BigUint {
        if num_bits == 0 {
            Zero::zero() // Python code returns the integer 0 here
        } else if num_bits <= 32 {
            // Fast path
            BigUint::from(self.genrand_uint32_default() >> (32 - num_bits))
        } else {
            // python is using 32-bit words here to match the MT - we find
            // out how many of those we need
            let words = (num_bits - 1) / 32 + 1;
            let mut v = Vec::with_capacity(words);
            let mut remaining_bits = num_bits;
            if cfg!(target_endian = "big") {
                // to implement this check _random_Random_getrandbits_impl
                // in _randommodule.c
                unimplemented!(
                    "Random bit generation for big endian platforms isn't implemented yet"
                )
            } else {
                for _ in 0..words {
                    let r = self.genrand_uint32_default();
                    v.push(if remaining_bits < 32 {
                        r >> (32 - remaining_bits) /* Drop least significant bits */
                    } else {
                        r
                    });
                    remaining_bits = remaining_bits.saturating_sub(32);
                }
            }
            BigUint::new(v)
        }
    }

    /// Generates a random integer below a given (exclusive) upper bound
    fn randbelow_with_getrandbits(&mut self, excl_upper_bound: &BigUint) -> BigUint {
        let num_bits = excl_upper_bound.bits();
        let mut r = self.getrandbits(num_bits as usize); // 0 <= r < 2**num_bits
        while &r >= excl_upper_bound {
            r = self.getrandbits(num_bits as usize);
        }
        r
    }

    /// randrange version with only upper bound specified; lower is implicitly 0
    pub fn randrange_from_zero(&mut self, excl_upper_bound: &BigUint) -> Option<BigUint> {
        if excl_upper_bound.is_zero() {
            None
        } else {
            Some(self.randbelow_with_getrandbits(excl_upper_bound))
        }
    }

    /// randrange version with both lower and upper bounds specified
    pub fn randrange(
        &mut self,
        incl_lower: &BigInt,
        excl_upper: &BigInt,
        step: &BigInt,
    ) -> Option<BigInt> {
        // width is strictly positive
        let width = {
            let w = excl_upper.checked_sub(incl_lower)?;
            if w.is_positive() {
                Some(BigUint::try_from(w).unwrap())
            } else {
                // note that we in particular return none if w is zero
                None
            }
        }?;
        if step.is_one() {
            let r = BigInt::from(self.randbelow_with_getrandbits(&width));
            Some(incl_lower.clone() + r)
        } else if step.is_zero() {
            None
        } else {
            let n: BigInt = if step.is_positive() {
                (BigInt::from(width) + step - BigInt::one()) / step
            } else {
                (BigInt::from(width) + step + BigInt::one()) / step
            };
            if n.is_positive() {
                let r =
                    BigInt::from(self.randbelow_with_getrandbits(&BigUint::try_from(n).unwrap()));
                Some(incl_lower + step * r)
            } else {
                None
            }
        }
    }
}

impl PyMt19937 {
    /// Allows making randomly selecting elements from a collection in a manner compatible with
    /// python's `random.choices`. Weighted choices currently not implemented - should be a
    /// rather easy addition though (it's basically just a cumsum and binary search).
    pub fn choices<T>(&mut self, population: impl IntoIterator<Item = T>) -> Choices<T> {
        Choices {
            rng: self,
            population: population.into_iter().collect(),
        }
    }
}

pub struct Choices<'a, T> {
    rng: &'a mut PyMt19937,
    population: Vec<T>,
}

impl<'a, T: Clone> Iterator for Choices<'a, T> {
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        let n = self.population.len() as f32;
        let i = (self.rng.genrand_res53_default() * n).floor() as usize;
        Some(self.population[i].clone())
    }
}

/// Iterator extension trait to provide fluent interface for python choices
pub trait RandomChoiceIterator: Iterator {
    fn choose<'a>(&mut self, rng: &'a mut PyMt19937) -> Choices<'a, Self::Item>
    where
        Self: Sized,
    {
        rng.choices(self)
    }
}

impl<T: Iterator> RandomChoiceIterator for T {}

/// Trait for python-compatibly seedable Rngs
pub trait PySeedable<Seed> {
    /// Seed an Rng in a way that produces results compatible with python
    fn py_seed(seed: Seed) -> Self;
}

impl PySeedable<Vec<u8>> for PyMt19937 {
    /// Python `bytes` and `bytearrays` are special cased as arguments to random.seed.
    /// We implement the current (version 2) behaviour of this special casing.
    /// To do so we turn the bytes into "large ints" and use those for seeding
    fn py_seed(mut bytes: Vec<u8>) -> Self {
        let mut hasher = Sha512::new();
        hasher.update(&bytes);
        let digest = hasher.finalize();
        bytes.extend_from_slice(&digest);
        // This flip is necessary to account for how python stores its ints / does
        // the bytearray conversion.
        // Note that the python code really makes an int from the bytes and later
        // bytes from the int. We omit this conversion and merely do this flip instead.
        if cfg!(target_endian = "big") {
            // the python implementation flips the array on big endian platforms.
            // We just don't flip it in that case.
        } else {
            bytes.reverse();
        }

        // right pad with enough bytes to make the total length divisible by 4
        bytes.extend(vec![0; 4 - bytes.len() % 4]);
        // convert padded bytes into 32-bit values by grouping in blocks of 4
        let key: Vec<_> = bytes
            .array_chunks()
            .map(|b| u32::from_le_bytes(*b))
            .collect();
        let stripped = strip_suffix_iter(&key, &[0]).unwrap_or(&key);
        PyMt19937::init_by_array(stripped)
    }
}

impl PySeedable<&[u8]> for PyMt19937 {
    /// Python `bytes` and `bytearrays` are special cased as arguments to random.seed.
    /// We implement the current (version 2) behaviour of this special casing.
    /// To do so we turn the bytes into "large ints" and use those for seeding
    fn py_seed(bytes: &[u8]) -> Self {
        PyMt19937::py_seed(bytes.to_vec())
    }
}

impl PySeedable<String> for PyMt19937 {
    /// Python strings are special cased as arguments to random.seed.
    /// We implement the current (version 2) behaviour of this special casing.
    /// Also see [init_from_py_bytes].
    fn py_seed(string: String) -> Self {
        PyMt19937::py_seed(string.into_bytes())
    }
}

impl PySeedable<&str> for PyMt19937 {
    /// Python strings are special cased as arguments to random.seed.
    /// We implement the current (version 2) behaviour of this special casing.
    /// Also see [init_from_py_bytes].
    fn py_seed(string: &str) -> Self {
        PyMt19937::py_seed(string.as_bytes().to_owned())
    }
}

impl PySeedable<u8> for PyMt19937 {
    fn py_seed(seed: u8) -> Self {
        PyMt19937::init_by_array(&[seed as u32])
    }
}

impl PySeedable<u16> for PyMt19937 {
    fn py_seed(seed: u16) -> Self {
        PyMt19937::init_by_array(&[seed as u32])
    }
}

impl PySeedable<u32> for PyMt19937 {
    fn py_seed(seed: u32) -> Self {
        PyMt19937::init_by_array(&[seed])
    }
}

/// Strips the given suffix as often as possible.
fn strip_suffix_iter<'a, T>(s: &'a [T], suffix: &[T]) -> Option<&'a [T]>
where
    T: PartialEq,
{
    if suffix.is_empty() {
        Some(s)
    } else if let Some(mut stripped) = s.strip_suffix(suffix) {
        while let Some(more_stripped) = stripped.strip_suffix(suffix) {
            stripped = more_stripped;
        }
        Some(stripped)
    } else {
        None
    }
}

impl PySeedable<&[u32]> for PyMt19937 {
    fn py_seed(key: &[u32]) -> Self {
        // python trims off blocks that only contain zeroes
        let stripped = strip_suffix_iter(key, &[0]).unwrap_or(key);
        PyMt19937::init_by_array(stripped)
    }
}

impl PySeedable<u64> for PyMt19937 {
    fn py_seed(seed: u64) -> Self {
        let bytes = seed.to_le_bytes();
        let key = [
            u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]),
            u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]),
        ];
        PyMt19937::py_seed(key.as_slice())
    }
}

impl PySeedable<u128> for PyMt19937 {
    fn py_seed(seed: u128) -> Self {
        let bytes = seed.to_le_bytes();
        dbg!(bytes);
        let key = [
            u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]),
            u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]),
            u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]),
            u32::from_le_bytes([bytes[12], bytes[13], bytes[14], bytes[15]]),
        ];
        PyMt19937::py_seed(key.as_slice())
    }
}

#[cfg(feature = "rand")]
mod rand_impl {
    use std::mem::MaybeUninit;

    use super::{PyMt19937, PySeedable};
    use rand::{RngCore, SeedableRng};

    impl RngCore for PyMt19937 {
        fn next_u32(&mut self) -> u32 {
            self.genrand_uint32_default()
        }

        fn next_u64(&mut self) -> u64 {
            let mut bytes = [0; 64 / 8];
            bytes[..4].copy_from_slice(&self.genrand_uint32_default().to_le_bytes());
            bytes[4..].copy_from_slice(&self.genrand_uint32_default().to_le_bytes());
            u64::from_le_bytes(bytes)
        }

        fn fill_bytes(&mut self, dest: &mut [u8]) {
            if !dest.is_empty() {
                let mut bytes = MaybeUninit::uninit();
                for (i, b) in dest.iter_mut().enumerate() {
                    let j = i % 4;
                    if j == 0 {
                        bytes = MaybeUninit::new(self.genrand_uint32_default().to_le_bytes());
                    }
                    *b = unsafe { bytes.assume_init()[j] };
                }
            }
        }

        fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
            self.fill_bytes(dest);
            Ok(())
        }
    }

    impl SeedableRng for PyMt19937 {
        type Seed = Vec<u8>;
        fn from_seed(seed: Self::Seed) -> Self {
            Self::py_seed(seed)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod getrandbits {
        use super::*;

        #[test]
        fn seed_bytes() {
            // values to compare against are generated using
            // `random.seed(bytes([1,2,3])); [random.getrandbits(10) for _ in range(10)]`
            // using python 3.13
            let correct = [933, 396, 47, 660, 910, 1022, 111, 767, 433, 95];
            let mut twister = PyMt19937::py_seed(vec![1, 2, 3]);
            assert_eq!(
                &(0..10)
                    .map(|_| twister.getrandbits(10).try_into().unwrap())
                    .collect::<Vec<u32>>(),
                &correct
            )
        }

        #[test]
        fn seed_str() {
            // values to compare against are generated using
            // `random.seed("Hello World 123"); [random.getrandbits(10) for _ in range(10)]`
            // using python 3.13
            let correct = [847, 216, 778, 72, 909, 144, 201, 835, 719, 840];
            let mut twister = PyMt19937::py_seed("Hello World 123");
            assert_eq!(
                &(0..10)
                    .map(|_| twister.getrandbits(10).try_into().unwrap())
                    .collect::<Vec<u32>>(),
                &correct
            )
        }

        #[test]
        fn seed_u8() {
            // values to compare against are generated using
            // `random.seed(123); [random.getrandbits(10) for _ in range(10)]`
            // using python 3.13
            let correct = [53, 274, 89, 787, 417, 272, 110, 858, 922, 895];
            let mut twister = PyMt19937::py_seed(123_u8);
            assert_eq!(
                &(0..10)
                    .map(|_| twister.getrandbits(10).try_into().unwrap())
                    .collect::<Vec<u32>>(),
                &correct
            )
        }

        #[test]
        fn seed_u16() {
            // values to compare against are generated using
            // `random.seed(312); [random.getrandbits(10) for _ in range(10)]`
            // using python 3.13
            let correct = [665, 148, 672, 391, 33, 755, 151, 334, 994, 237];
            let mut twister = PyMt19937::py_seed(312_u16);
            assert_eq!(
                &(0..10)
                    .map(|_| twister.getrandbits(10).try_into().unwrap())
                    .collect::<Vec<u32>>(),
                &correct
            )
        }

        #[test]
        fn seed_u32() {
            // values to compare against are generated using
            // `random.seed(65_536); [random.getrandbits(10) for _ in range(10)]`
            // using python 3.13
            let correct = [426, 758, 954, 800, 651, 78, 434, 100, 755, 58];
            let mut twister = PyMt19937::py_seed(65_536_u32);
            assert_eq!(
                &(0..10)
                    .map(|_| twister.getrandbits(10).try_into().unwrap())
                    .collect::<Vec<u32>>(),
                &correct
            )
        }

        #[test]
        fn seed_u64() {
            // values to compare against are generated using
            // `random.seed((2<<31) + 51_231); [random.getrandbits(10) for _ in range(10)]`
            // using python 3.13
            let correct = [992, 608, 345, 718, 115, 503, 696, 707, 387, 763];
            let mut twister = PyMt19937::py_seed(2_u64.pow(32) + 51_231);
            assert_eq!(
                &(0..10)
                    .map(|_| twister.getrandbits(10).try_into().unwrap())
                    .collect::<Vec<u32>>(),
                &correct
            )
        }

        #[test]
        fn seed_u128() {
            // values to compare against are generated using
            // `random.seed((2<<63) + 51_239_123); [random.getrandbits(10) for _ in range(10)]`
            // using python 3.13
            let correct = [502, 568, 107, 39, 744, 290, 438, 109, 142, 290];
            let mut twister = PyMt19937::py_seed(2_u128.pow(64) + 51_239_123);
            assert_eq!(
                &(0..10)
                    .map(|_| twister.getrandbits(10).try_into().unwrap())
                    .collect::<Vec<u32>>(),
                &correct
            )
        }

        #[test]
        fn small_value_large_type() {
            // values to compare against are generated using
            // `random.seed(123); [random.getrandbits(10) for _ in range(10)]`
            // using python 3.13
            let correct = [53, 274, 89, 787, 417, 272, 110, 858, 922, 895];
            let mut twister = PyMt19937::py_seed(123_u128);
            assert_eq!(
                &(0..10)
                    .map(|_| twister.getrandbits(10).try_into().unwrap())
                    .collect::<Vec<u32>>(),
                &correct
            )
        }
    }

    #[test]
    fn choices() {
        // values to compare against are generated using
        // `random.seed("Magst du Pizza? Yope"); random.choices(range(20), k=10)`
        // using python 3.13
        let correct = [6, 2, 6, 8, 11, 7, 17, 7, 8, 5];
        let mut twister = PyMt19937::py_seed("Magst du Pizza? Yope");
        assert_eq!(
            twister.choices(0..20).take(10).collect::<Vec<_>>(),
            &correct
        )
    }

    #[test]
    fn randbelow() {
        // `random.seed("Pizza"); [random.randrange(100) for i in range(10)]`
        let correct = [25, 52, 24, 70, 14, 72, 70, 56, 43, 32];
        let mut twister = PyMt19937::py_seed("Pizza");
        assert_eq!(
            &(0..10)
                .map(|_| twister
                    .randbelow_with_getrandbits(&BigUint::from(100_u32))
                    .try_into()
                    .unwrap())
                .collect::<Vec<u32>>(),
            &correct
        );
        // `[random.randrange(10_000) for i in range(10)]`
        let correct = [2306, 7958, 1402, 2301, 2970, 7787, 5438, 6340, 7300, 140];
        assert_eq!(
            &(0..10)
                .map(|_| twister
                    .randbelow_with_getrandbits(&BigUint::from(10_000_u32))
                    .try_into()
                    .unwrap())
                .collect::<Vec<u32>>(),
            &correct
        );
    }

    #[test]
    fn rangrange_from_zero() {
        // `random.seed("Pizza"); [random.randrange(43_057) for i in range(10)]`
        let correct = [
            13098, 26831, 12744, 36225, 7605, 37058, 35879, 28724, 22476, 16773,
        ];
        let mut twister = PyMt19937::py_seed("Pizza");
        assert_eq!(
            &(0..10)
                .map(|_| twister
                    .randrange_from_zero(&BigUint::from(43_057_u32))
                    .unwrap()
                    .try_into()
                    .unwrap())
                .collect::<Vec<u32>>(),
            &correct
        );
    }

    #[test]
    fn rangrange_unit_step() {
        // `random.seed("Pizza"); [random.randrange(42, 43_057) for i in range(10)]`
        let correct = [
            13140, 26873, 12786, 36267, 7647, 37100, 35921, 28766, 22518, 16815,
        ];
        let mut twister = PyMt19937::py_seed("Pizza");
        assert_eq!(
            &(0..10)
                .map(|_| twister
                    .randrange(
                        &BigInt::from(42_u32),
                        &BigInt::from(43_057_u32),
                        &BigInt::one()
                    )
                    .unwrap()
                    .try_into()
                    .unwrap())
                .collect::<Vec<u32>>(),
            &correct
        );
    }

    #[test]
    fn rangrange_nonunit_step() {
        // `random.seed("Pizza"); [random.randrange(42, 43_057, 12) for i in range(10)]`
        let correct = [
            9858, 20154, 9594, 38886, 27210, 5742, 27834, 40998, 26946, 21582,
        ];
        let mut twister = PyMt19937::py_seed("Pizza");
        assert_eq!(
            &(0..10)
                .map(|_| twister
                    .randrange(
                        &BigInt::from(42_u32),
                        &BigInt::from(43_057_u32),
                        &BigInt::from(12_u32)
                    )
                    .unwrap()
                    .try_into()
                    .unwrap())
                .collect::<Vec<u32>>(),
            &correct
        );
    }

    #[test]
    fn rangrange_large() {
        // `random.seed("Pizza"); [random.randrange(42, 10**100) for i in range(10)]`
        let correct = [
            BigInt::parse_bytes(b"9579866485220884275812008134464035965241537083189421689381792593824445882708815065843427244115679591", 10).unwrap(),
            BigInt::parse_bytes(b"3888203122170713226666112762701632817564811813426175309426977673877400812229713204130602887840425185", 10).unwrap(),
            BigInt::parse_bytes(b"5671731618993807604597264481562670303958256589845472646965893604588043197507096527112390358858623451", 10).unwrap(),
            BigInt::parse_bytes(b"4194032766356154954455712807769223619918079113218561184406604612474027197687254690833120987693660483", 10).unwrap(),
            BigInt::parse_bytes(b"7476178950561835834867063441551291285800190891810969659540612285203169282520874331082648528181368436", 10).unwrap(),
            BigInt::parse_bytes(b"1861103480107524638919591893137299509305415151477751224752864145416695733183243089594346927767194340", 10).unwrap(),
            BigInt::parse_bytes(b"3922203142767393832538673954495653361515200423210709219982766035391585977477681283927291826719657524", 10).unwrap(),
            BigInt::parse_bytes(b"6128897628109692866058641368274476540862767831381945527933770759683809548646997692265219366175075688", 10).unwrap(),
            BigInt::parse_bytes(b"9584709203694536284250111641038241722576534138774799060713437581067222978206646281627216084226872580", 10).unwrap(),
            BigInt::parse_bytes(b"3033019269308996974263424226402724479542287417888143915688162688154376377154955676213661065686678850", 10).unwrap()
        ];
        let mut twister = PyMt19937::py_seed("Pizza");
        assert_eq!(
            &(0..10)
                .map(|_| twister
                    .randrange(
                        &BigInt::from(42_u32),
                        &BigInt::from(10_u32).pow(100),
                        &BigInt::one(),
                    )
                    .unwrap())
                .collect::<Vec<BigInt>>(),
            &correct
        );
    }

    #[cfg(feature = "rand")]
    mod randtest {
        use super::*;
        use rand::{Rng, SeedableRng};

        #[test]
        fn seed_bytes() {
            // values to compare against are generated using
            // `random.seed(bytes([1,2,3])); [random.getrandbits(32) for _ in range(10)]`
            // using python 3.13
            let correct = [
                3914323484, 1663236869, 198133927, 2769096544, 3820271845, 4288717126, 469685982,
                3218529643, 1817328981, 401606798,
            ];
            let mut twister = PyMt19937::from_seed(vec![1, 2, 3]);
            assert_eq!(&twister.gen::<[u32; 10]>(), &correct)
        }
    }
}

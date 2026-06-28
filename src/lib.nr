mod scalar_field;
mod test;
pub mod ed25519;

pub use crate::scalar_field::ScalarField;
use bignum::{BigNum, bignum::evaluate_quadratic_expression};
use std::ops::{Add, Neg, Sub};

/// A point on a Twisted Edwards curve whose coordinates live in an arbitrary
/// prime field represented by the `BigNum` type `B`.
///
/// `B` is the base field of the curve (the field the x/y coordinates live in).
/// `Params` supplies the curve constants `a`, `d` and the generator (all of type `B`).
pub struct Curve<B, Params> {
    pub x: B,
    pub y: B,
}

// ####################################################################################################################
// ####################################################################################################################
// ### T R A I T S
// ####################################################################################################################
// ####################################################################################################################

/// Parametrises a Twisted Edwards curve defined over the base field `B`.
pub trait TECurveParameterTrait<B>
where
    B: BigNum,
{
    fn a() -> B; // twisted edwards curve parameter a
    fn d() -> B; // twisted edwards curve parameter d
    fn gen() -> (B, B); // generator point x/y coordinates
}

/// Defines methods that a valid Curve implementation must satisfy
pub trait CurveTrait<B, Params>: Add + Sub + Eq + Neg + Default {
    fn new(x: B, y: B) -> Self;
    fn zero() -> Self;
    fn one() -> Self;
    fn dbl(self) -> Self;
    fn mul<let NScalarSlices: u32>(self, x: ScalarField<NScalarSlices>) -> Self;
    fn msm<let N: u32, let NScalarSlices: u32>(
        points: [Self; N],
        scalars: [ScalarField<NScalarSlices>; N],
    ) -> Self;

    fn eq(self, x: Self) -> bool {
        self == x
    }
    fn is_zero(self) -> bool {
        self == Self::zero()
    }

    fn is_on_curve(self) -> bool;
    fn assert_is_on_curve(self);
    fn assert_equal(self, other: Self);
}

// ####################################################################################################################
// ####################################################################################################################
// ### C O N S T R A I N E D    F U N C T I O N S
// ####################################################################################################################
// ####################################################################################################################
impl<B, Params> std::default::Default for Curve<B, Params>
where
    B: BigNum,
    Params: TECurveParameterTrait<B>,
{
    /// Returns point at infinity
    fn default() -> Self {
        Curve::zero()
    }
}

impl<B, Params> std::ops::Add for Curve<B, Params>
where
    B: BigNum,
    Params: TECurveParameterTrait<B>,
{
    /// Compute `self + other`
    fn add(self, other: Self) -> Self {
        Curve::add_internal(self, other, Params::a(), Params::d())
    }
}

impl<B, Params> std::ops::Neg for Curve<B, Params>
where
    B: BigNum,
{
    /// Negate a point
    fn neg(self) -> Self {
        Curve { x: self.x.neg(), y: self.y }
    }
}

impl<B, Params> std::ops::Sub for Curve<B, Params>
where
    B: BigNum,
    Params: TECurveParameterTrait<B>,
{
    /// Compute `self - other`
    fn sub(self, other: Self) -> Self {
        Curve::add_internal(self, other.neg(), Params::a(), Params::d())
    }
}

impl<B, Params> std::cmp::Eq for Curve<B, Params>
where
    B: BigNum,
{
    /// Compute `self == other`
    fn eq(self, other: Self) -> bool {
        (self.x == other.x) & (self.y == other.y)
    }
}

impl<B, Params> std::convert::From<(B, B)> for Curve<B, Params> {
    /// Construct from tuple of base field elements
    ///
    /// Use this method instead of `new` if you know x/y is on the curve
    fn from((x, y): (B, B)) -> Self {
        Curve { x, y }
    }
}

impl<B, Params> CurveTrait<B, Params> for Curve<B, Params>
where
    B: BigNum,
    Params: TECurveParameterTrait<B>,
{

    /// Construct a new point
    ///
    /// If you know the x/y coords form a valid point DO NOT USE THIS METHOD
    /// This method calls `assert_is_on_curve`.
    /// Instead, directly construct via Curve{x, y} or use from((x, y))
    fn new(x: B, y: B) -> Self {
        let result = Curve { x, y };
        result.assert_is_on_curve();
        result
    }

    /// Return the Identity element (point at infinity), which is `(0, 1)`
    fn zero() -> Self {
        Curve { x: B::zero(), y: B::one() }
    }

    /// Return the Generator of the group
    fn one() -> Self {
        let (x, y) = Params::gen();
        Curve { x, y }
    }

    /// Validate a point is on the curve
    ///
    /// cheaper than `is_on_curve` (assert is cheaper than returning a bool)
    ///
    /// Checks `a*x^2 + y^2 == 1 + d*x^2*y^2`
    fn assert_is_on_curve(self) {
        let a = Params::a();
        let d = Params::d();
        let t0 = self.x.sqr();
        let t1 = self.y.sqr();
        let t2 = a * t0 + t1;
        let t3 = B::one() + d * t0 * t1;
        assert(t2 == t3);
    }

    /// Constrain two points to equal each other
    fn assert_equal(self, other: Self) {
        assert(self.x == other.x);
        assert(self.y == other.y);
    }

    /// Return a bool that describes whether the point is on the curve
    ///
    /// If you don't need to handle the failure case, it is cheaper to call `assert_is_on_curve`
    fn is_on_curve(self) -> bool {
        let a = Params::a();
        let d = Params::d();
        let t0 = self.x.sqr();
        let t1 = self.y.sqr();
        let t2 = a * t0 + t1;
        let t3 = B::one() + d * t0 * t1;
        (t2 == t3)
    }

    /// Compute `self + self`
    fn dbl(self) -> Self {
        Curve::dbl_internal(self, Params::a(), Params::d())
    }

    /// Compute `self * scalar`
    ///
    /// Uses the Straus method via lookup tables.
    /// Assumes backend has an efficient implementation of a memory table abstraction
    /// i.e. `let x = table[y]` is efficient even if `y` is not known at compile time
    fn mul<let NScalarSlices: u32>(self: Self, scalar: ScalarField<NScalarSlices>) -> Self {
        // define a, d params locally to make code more readable (shouldn't affect performance)
        let a = Params::a();
        let d = Params::d();

        // Construct tables of precomputed point coordinates.
        let (table_x, table_y): ([B; 16], [B; 16]) = self.compute_straus_point_table(a, d);

        // Initialize the accumulator with the point that maps to the first (most significant) scalar slice
        let idx = scalar.base4_slices[0] as u32;
        let mut accumulator: Self = Curve { x: table_x[idx], y: table_y[idx] };

        // Execute a double-and-add subroutine
        // 1. Compute `accumulator = accumulator * 16`
        // 2. Extract 4-bits from the scalar multiplier and
        //    use them to retrieve the corresponding point from our point table
        // Note: this is similar to the "double and add" scalar multiplication method, except we use base16 instead of base2!
        for i in 1..NScalarSlices {
            accumulator = accumulator.dbl_internal(a, d);
            accumulator = accumulator.dbl_internal(a, d);
            accumulator = accumulator.dbl_internal(a, d);
            accumulator = accumulator.dbl_internal(a, d);
            let idx = scalar.base4_slices[i] as u32;
            let x = table_x[idx];
            let y = table_y[idx];
            accumulator = accumulator.add_internal(Curve { x, y }, a, d);
        }

        if (scalar.skew) {
            accumulator = accumulator - self;
        }
        accumulator
    }

    /// compute `points[0] * scalar[0] + ... + points[N-1] * scalar[N-1]`
    ///
    /// Is cheaper than `mul` when processing >1 point due to reduced number of point doublings
    /// uses the Straus MSM method via lookup tables.
    fn msm<let N: u32, let NScalarSlices: u32>(
        points: [Self; N],
        scalars: [ScalarField<NScalarSlices>; N],
    ) -> Self {
        let a = Params::a();
        let d = Params::d();

        // Generalized version of `mul` for multiple points.
        let mut point_tables: [([B; 16], [B; 16]); N] =
            [([B::zero(); 16], [B::zero(); 16]); N];
        for j in 0..N {
            point_tables[j] = points[j].compute_straus_point_table(a, d);
        }

        let idx = scalars[0].base4_slices[0] as u32;
        let mut accumulator: Self = Curve { x: point_tables[0].0[idx], y: point_tables[0].1[idx] };
        for j in 1..N {
            let idx = scalars[j].base4_slices[0] as u32;
            let P = Curve { x: point_tables[j].0[idx], y: point_tables[j].1[idx] };
            accumulator = accumulator.add_internal(P, a, d);
        }
        for i in 1..NScalarSlices {
            accumulator = accumulator.dbl_internal(a, d);
            accumulator = accumulator.dbl_internal(a, d);
            accumulator = accumulator.dbl_internal(a, d);
            accumulator = accumulator.dbl_internal(a, d);
            for j in 0..N {
                let idx = scalars[j].base4_slices[i] as u32;
                let x = point_tables[j].0[idx];
                let y = point_tables[j].1[idx];
                accumulator = accumulator.add_internal(Curve { x, y }, a, d);
            }
        }

        for j in 0..N {
            if (scalars[j].skew == true) {
                accumulator = accumulator - points[j];
            }
        }
        accumulator
    }
}

// ####################################################################################################################
// ####################################################################################################################
// ### H E L P E R    F U N C T I O N S
// ####################################################################################################################
// ####################################################################################################################
impl<B, Params> Curve<B, Params>
where
    B: BigNum,
{

    /// add two points together using the unified (complete) twisted Edwards addition law:
    ///
    ///     x3 = (x1*y2 + y1*x2) / (1 + d*x1*x2*y1*y2)
    ///     y3 = (y1*y2 - a*x1*x2) / (1 - d*x1*x2*y1*y2)
    ///
    /// This method takes `a` and `d` as arguments (rather than reading them from `Params`)
    /// because of a Noir limitation where `Params` cannot be accessed by an internal function
    /// called from another internal function.
    fn add_internal(self, other: Self, a: B, d: B) -> Self {
        let x1 = self.x;
        let y1 = self.y;
        let x2 = other.x;
        let y2 = other.y;

        // Compute every witness unconstrained (no constraints emitted here), then
        // pin them down below with 6 quadratic-expression reductions instead of the
        // ~9 a naive sequence of `*` / `/` would cost. The witness computation lives
        // in its own `unconstrained fn` so that this is a single Brillig call whose
        // six outputs are each constrained directly below (satisfying the soundness
        // lint), rather than a chain of intermediate unconstrained values.
        //
        // Safety: all six returned witnesses are constrained by the
        // `evaluate_quadratic_expression` calls that follow.
        let (x1x2, y1y2, e, t, x3, y3) =
            unsafe { __add_witnesses(x1, y1, x2, y2, a, d) };

        // 1. x1x2 == x1 * x2
        evaluate_quadratic_expression([[x1]], [[false]], [[x2]], [[false]], [x1x2], [true]);
        // 2. y1y2 == y1 * y2
        evaluate_quadratic_expression([[y1]], [[false]], [[y2]], [[false]], [y1y2], [true]);
        // 3. e == x1x2 * y1y2
        evaluate_quadratic_expression([[x1x2]], [[false]], [[y1y2]], [[false]], [e], [true]);
        // 4. t == d * e
        evaluate_quadratic_expression([[d]], [[false]], [[e]], [[false]], [t], [true]);
        // 5. x3*(1 + t) - x1*y2 - y1*x2 == 0   (x3*t and x3*1 split into product + linear)
        evaluate_quadratic_expression(
            [[x3], [x1], [y1]],
            [[false], [true], [true]],
            [[t], [y2], [x2]],
            [[false], [false], [false]],
            [x3],
            [false],
        );
        // 6. y3*(1 - t) - y1*y2 + a*x1*x2 == 0   (= y3 - y3*t - y1y2 + a*x1x2)
        evaluate_quadratic_expression(
            [[y3], [a]],
            [[true], [false]],
            [[t], [x1x2]],
            [[false], [false]],
            [y3, y1y2],
            [false, true],
        );

        Self { x: x3, y: y3 }
    }

    /// add a point to itself
    ///
    /// Uses the dedicated twisted Edwards doubling formula:
    ///
    ///     x3 = (2*x*y)         / (a*x^2 + y^2)
    ///     y3 = (y^2 - a*x^2)   / (2 - a*x^2 - y^2)
    ///
    /// The curve equation `a*x^2 + y^2 == 1 + d*x^2*y^2` lets us compute the denominators
    /// without the `d*x^2*y^2` product, so this is one reduction cheaper than reusing
    /// `add_internal`. `d` is unused but kept in the signature for a uniform call site.
    fn dbl_internal(self, a: B, _d: B) -> Self {
        let x = self.x;
        let y = self.y;

        // Safety: all five witnesses are constrained by the
        // `evaluate_quadratic_expression` calls below. The doubling formula is complete
        // for the supported curves, so the denominators are never zero.
        let (u, v, au, x3, y3) = unsafe { __dbl_witnesses(x, y, a) };

        // 1. u == x * x
        evaluate_quadratic_expression([[x]], [[false]], [[x]], [[false]], [u], [true]);
        // 2. v == y * y
        evaluate_quadratic_expression([[y]], [[false]], [[y]], [[false]], [v], [true]);
        // 3. au == a * u
        evaluate_quadratic_expression([[a]], [[false]], [[u]], [[false]], [au], [true]);
        // 4. x3*(au + v) - 2*x*y == 0
        evaluate_quadratic_expression(
            [[x3], [x3], [x], [x]],
            [[false], [false], [false], [false]],
            [[au], [v], [y], [y]],
            [[false], [false], [true], [true]],
            [],
            [],
        );
        // 5. y3*(2 - au - v) - (v - au) == 0   (= 2*y3 - y3*au - y3*v + au - v)
        evaluate_quadratic_expression(
            [[y3], [y3]],
            [[false], [false]],
            [[au], [v]],
            [[true], [true]],
            [y3, y3, au, v],
            [false, false, false, true],
        );

        Self { x: x3, y: y3 }
    }

    /// Compute a 4-bit lookup table of point multiples for the Straus windowed scalar multiplication algorithm.
    ///
    /// Table contains the signed odd multiples used by the wNAF-style scalar mul algorithm to
    /// minimize the total number of required point additions.
    ///
    /// We use `([B; 16], [B; 16])` (separate x and y tables) rather than `[Curve; 16]` so that a
    /// runtime index can be used directly as the argument to a memory lookup for each coordinate.
    fn compute_straus_point_table(self, a: B, d: B) -> ([B; 16], [B; 16]) {
        let mut table_x: [B; 16] = [B::zero(); 16];
        let mut table_y: [B; 16] = [B::zero(); 16];
        table_x[8] = self.x;
        table_y[8] = self.y;
        let D = self.dbl_internal(a, d);
        for i in 1..8 {
            let Q = Self { x: table_x[7 + i], y: table_y[7 + i] };
            let V = D.add_internal(Q, a, d);
            table_x[8 + i] = V.x;
            table_y[8 + i] = V.y;
        }
        for i in 0..8 {
            table_x[i] = table_x[15 - i].neg();
            table_y[i] = table_y[15 - i];
        }
        (table_x, table_y)
    }
}

/// Unconstrained witness computation for `add_internal`.
///
/// Returns `(x1x2, y1y2, e, t, x3, y3)` where `e = x1*x2*y1*y2`, `t = d*e` and
/// `(x3, y3)` is the twisted Edwards sum. All outputs are constrained by the caller.
unconstrained fn __add_witnesses<B>(x1: B, y1: B, x2: B, y2: B, a: B, d: B) -> (B, B, B, B, B, B)
where
    B: BigNum,
{
    let one = B::one();
    let x1x2 = x1.__mul(x2);
    let y1y2 = y1.__mul(y2);
    let e = x1x2.__mul(y1y2); // x1*x2*y1*y2
    let t = d.__mul(e); // d*x1*x2*y1*y2
    let num_x = x1.__mul(y2).__add(y1.__mul(x2)); // x1*y2 + y1*x2
    let num_y = y1y2.__sub(a.__mul(x1x2)); // y1*y2 - a*x1*x2
    let x3 = num_x.__div(one.__add(t)); // / (1 + t)
    let y3 = num_y.__div(one.__sub(t)); // / (1 - t)
    (x1x2, y1y2, e, t, x3, y3)
}

/// Unconstrained witness computation for `dbl_internal`.
///
/// Returns `(u, v, au, x3, y3)` where `u = x^2`, `v = y^2`, `au = a*x^2` and
/// `(x3, y3)` is the doubled point. All outputs are constrained by the caller.
unconstrained fn __dbl_witnesses<B>(x: B, y: B, a: B) -> (B, B, B, B, B)
where
    B: BigNum,
{
    let u = x.__mul(x); // x^2
    let v = y.__mul(y); // y^2
    let au = a.__mul(u); // a*x^2
    let xy = x.__mul(y);
    let num_x = xy.__add(xy); // 2*x*y
    let num_y = v.__sub(au); // y^2 - a*x^2
    let one = B::one();
    let two = one.__add(one);
    let x3 = num_x.__div(au.__add(v)); // / (a*x^2 + y^2)      == 1 + t
    let y3 = num_y.__div(two.__sub(au).__sub(v)); // / (2 - a*x^2 - y^2)  == 1 - t
    (u, v, au, x3, y3)
}

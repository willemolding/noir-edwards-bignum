use crate::Curve;
use crate::CurveTrait;
use crate::ed25519::Ed25519;
use crate::scalar_field::ScalarField;
use bignum::BigNum;
use std::ops::Neg;

/// Construct an ed25519 point directly from the 120-bit limb representation of its coordinates.
///
/// Unlike BabyJubJub, the ed25519 base field (2^255 - 19) does not fit in a native `Field`,
/// so coordinates have to be supplied as `BigNum` limbs rather than `Field` literals.
fn point(x: [u128; 3], y: [u128; 3]) -> Ed25519 {
    Curve { x: BigNum::from_limbs(x), y: BigNum::from_limbs(y) }
}

fn assert_eq_point(result: Ed25519, expected: Ed25519) {
    assert(result.x == expected.x);
    assert(result.y == expected.y);
}

#[test]
fn test_generator_on_curve() {
    Ed25519::one().assert_is_on_curve();
    assert(Ed25519::one().is_on_curve());
}

#[test]
fn test_add_sub() {
    let b = Ed25519::one();
    // 3B computed two ways
    let expected = b + (b + b);
    // 4B - B
    let result = b.dbl().dbl() - b;
    assert_eq_point(result, expected);
}

#[test]
fn test_scalar_mul_small() {
    // 9 * B, checked against a manual 8B + B and an independently computed vector
    let scalar_f: ScalarField<64> = ScalarField::from(9);
    let result = Ed25519::one().mul(scalar_f);

    let b = Ed25519::one();
    let mut manual = b.dbl().dbl().dbl();
    manual = manual + b;
    assert_eq_point(result, manual);

    let expected = point(
        [0xd80e5c2e4e0294d3e023065185715c, 0xc970c80071651bf336e06f9422b886, 0x357c],
        [0x789013236e4430c74e44845522f1c0, 0x23c2c2dd0df4b2befce956f2d2fd1f, 0x7f3d],
    );
    assert_eq_point(result, expected);
}

#[test]
fn test_scalar_mul_vector() {
    let scalar_f: ScalarField<64> = ScalarField::from(1234567890123456789);
    let result = Ed25519::one().mul(scalar_f);
    let expected = point(
        [0x2429f2f7a0427ce897e4729187bd2f, 0xd16f1a3e8e66055591a596e09d9acc, 0x5e3],
        [0x570534fb765a50f0f60b0c636a4047, 0xcbeda1207673018275d54d0b743d3c, 0x2eae],
    );
    assert_eq_point(result, expected);
}

#[test]
fn test_scalar_mul_large_vector() {
    // A scalar just below 2^252 (still < native field modulus).
    let scalar = 6453482971500802406899721436289478618597120925962473168869639197649048738061;
    let scalar_f: ScalarField<64> = ScalarField::from(scalar);
    let result = Ed25519::one().mul(scalar_f);
    let expected = point(
        [0xfac89db53a54f786b16ee3aaef44ee, 0x5f918207f7c0b743ac985bab703ef4, 0x3d05],
        [0xb2ccc3337196af2efd33cded083827, 0xf4dd8147ec646d713c4576687b184c, 0x775a],
    );
    assert_eq_point(result, expected);
}

#[test]
fn test_msm() {
    let b = Ed25519::one();
    let p2 = point(
        [0xd8e4515d7a45f83c5a14e2843ce0e, 0x384c9f5a046c3d043b7d1833e7ac08, 0x36ab],
        [0x390f51643851560e5f46ae6af8a3c9, 0xcdf3092329c21da25ee8c9a21f5697, 0x2260],
    ); // 2B
    let p3 = point(
        [0xe809f3bbf1ce955855981af50e4107, 0x8685fcf4bd4ee9e3ee194b1d810783, 0x1456],
        [0x882b46f7abf23d12c4c4b59f4062b8, 0x63e32b47d52f87ce6468dd36ad41f0, 0x31c5],
    ); // 7B

    let points = [b, p2, p3];
    let scalars = [
        ScalarField::<64>::from(1111111111111111111111111111111111111111),
        ScalarField::<64>::from(2222222222222222222222222222222222222222),
        ScalarField::<64>::from(3333333333333333333333333333333333333333),
    ];
    let result = Curve::msm(points, scalars);

    let expected = point(
        [0xf2b36a7282d3e08863181f6d408af5, 0xb8cc72c03f895a25592f2a5751cff0, 0x90a],
        [0x544f51dcd89cf2889f8979e752cf5, 0xd54d3ed804e374ed5d14809409d843, 0x30a6],
    );
    assert_eq_point(result, expected);
}

#[test]
fn test_full_range_scalar_needs_64_slices() {
    // ed25519's group order ell is *just above* 2^252, so a reduced scalar can exceed the
    // 2^252 range that `ScalarField<63>` covers. `ScalarField<64>` is required.
    //
    // (ell - 1) * B == -B, which we use as the expected value without an external vector.
    let ell_minus_1 = 7237005577332262213973186563042994240857116359379907606001950938285454250988;
    let scalar_f: ScalarField<64> = ScalarField::from(ell_minus_1);
    let result = Ed25519::one().mul(scalar_f);
    assert_eq_point(result, Ed25519::one().neg());
}

#[test]
fn test_scalar_field_64_boundaries() {
    let values: [Field; 6] = [
        0,
        1,
        0xffffffffffffffffffffffffffffffff,
        0x100000000000000000000000000000000,
        0x1000000000000000000000000000000000000000000000000000000000000000,
        -1,
    ];
    for x in values {
        let scalar_f: ScalarField<64> = ScalarField::from(x);
        assert(ScalarField::<64>::into(scalar_f) == x);
    }
}

#[test]
fn test_scalar_mul_64_zero() {
    let scalar_f: ScalarField<64> = ScalarField::from(0);
    assert(Ed25519::one().mul(scalar_f).is_zero());
}

#[test]
// Check the honest decomposition at group level: compare `mul` against a
// bitwise double-and-add reference at the boundary scalars. Rejection of
// aliased decompositions is covered by the `assert_wnaf_decomposition_binds`
// tests in `scalar_field.nr`.
fn test_full_scalar_mul_matches_bitwise() {
    let point = Ed25519::one();
    let values: [Field; 5] =
        [0, 1, 0xffffffffffffffffffffffffffffffff, 0x100000000000000000000000000000000, -1];
    for x in values {
        let scalar: ScalarField<64> = ScalarField::from(x);
        let mut expected = Ed25519::zero();
        for bit in x.to_be_bits::<254>() {
            expected = expected.dbl();
            expected = if bit { expected + point } else { expected };
        }
        assert_eq_point(point.mul(scalar), expected);
    }
}

#[test(should_fail_with = "Point is not on the curve")]
fn test_new_rejects_point_not_on_curve() {
    // The generator with y + 1.
    let _ = Ed25519::new(
        BigNum::from_limbs([
            0x2cc7609525a7b2c9562d608f25d51a,
            0x36d3cd6e53fec0a4e231fdd6dc5c69,
            0x2169,
        ]),
        BigNum::from_limbs([
            0x666666666666666666666666666659,
            0x666666666666666666666666666666,
            0x6666,
        ]),
    );
}

#[test(should_fail_with = "Point is not on the curve")]
fn test_assert_is_on_curve_rejects_point_not_on_curve() {
    point([0, 0, 0], [0, 0, 0]).assert_is_on_curve();
}

#[test]
fn test_assert_equal() {
    let b = Ed25519::one();
    b.dbl().assert_equal(b + b);
}

#[test(should_fail_with = "Points are not equal")]
fn test_assert_equal_rejects_different_x() {
    // `b` and `-b` share the y coordinate and differ in x.
    let b = Ed25519::one();
    b.assert_equal(b.neg());
}

#[test(should_fail_with = "Points are not equal")]
fn test_assert_equal_rejects_different_y() {
    // `(0, -1)` and `(0, 1)` share the x coordinate and differ in y.
    let minus_one = [0xffffffffffffffffffffffffffffec, 0xffffffffffffffffffffffffffffff, 0x7fff];
    point([0, 0, 0], minus_one).assert_equal(Ed25519::zero());
}
